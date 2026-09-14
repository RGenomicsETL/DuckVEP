#!/usr/bin/env Rscript
# Record-minimized replay: one complete source VCF record, unchanged alleles and
# model, through the existing three-engine comparison. This is not allele or
# transcript-model minimization. Each worker checkpoints one compact Parquet.
duckvep_replay_lanes <- c("duckvep_vep_csq", "fastvep_vep_csq", "native_tab17")
duckvep_replay_columns <- c("lane", "record_index", "alt_index", "Feature", "field", "actual", "expected")
duckvep_replay_context <- c("run_hash", "task_index", "profile", "original_record_index", "input_id",
  "header", "record", "input_sha256", "model_sha256", "status", "passed", "summary", "outputs", "error",
  "retained_directory")

duckvep_replay_indices <- function(index, limit) {
  text <- as.character(index)
  values <- suppressWarnings(as.numeric(text))
  stopifnot(!anyNA(text), all(grepl("^[1-9][0-9]*$", text)), !anyNA(values),
    all(values <= limit), all(values <= .Machine$integer.max))
  as.integer(values)
}

duckvep_replay_canonical <- function(x) {
  x <- x[, duckvep_replay_columns, drop = FALSE]
  x[] <- lapply(x, as.character)
  x <- x[do.call(order, c(unname(x), list(na.last = TRUE))), , drop = FALSE]
  rownames(x) <- NULL
  x
}

duckvep_replay_record <- function(lines, source, index) {
  records <- lines[!startsWith(lines, "#")]
  index <- duckvep_replay_indices(index, length(records))
  stopifnot(length(index) == 1L, nrow(source) == 1L,
    as.character(source$record_index) == as.character(index),
    as.character(source$alt_index) == "1")
  fields <- strsplit(records[[index]], "\t", fixed = TRUE)[[1L]]
  stopifnot(length(fields) >= 8L, identical(fields[[3L]], source$Uploaded_variation),
    identical(fields[c(1, 2, 4, 5)], unname(as.character(source[1L, c("CHROM", "POS", "REF", "ALT")]))))
  list(header = lines[startsWith(lines, "#")], record = records[[index]])
}

duckvep_replay_summary <- function(summary, profile, records) {
  # This replay corpus has exactly one annotation key per physical record in
  # each lane. Refuse other cardinalities instead of assuming singleton counts.
  checks <- c("missing_keys", "extra_keys", "actual_duplicate_keys", "expected_duplicate_keys", "invalid_keys",
    "source_duplicate_alleles", "source_invalid_alleles", "actual_missing_source_alleles",
    "expected_missing_source_alleles", "actual_unknown_alleles", "expected_unknown_alleles")
  counts <- c("input_records", "input_alleles", "actual_rows", "expected_rows", "union_keys",
    "compared_keys", "source_rows", "source_alleles")
  stopifnot(nrow(summary) == 3L,
    setequal(summary$comparison, duckvep_replay_lanes), all(summary$case == profile),
    all(unlist(summary[counts]) == records),
    all(summary$compared_fields == ifelse(summary$comparison == "native_tab17", 16L, 31L)),
    all(unlist(summary[checks]) == 0))
}

duckvep_replay_failures <- function(rows, source, keys, summary) {
  source_indices <- duckvep_replay_indices(source$record_index, nrow(source))
  stopifnot(setequal(source_indices, seq_len(nrow(source))), !anyDuplicated(source_indices),
    all(source$alt_index == 1L), !anyDuplicated(source$Uploaded_variation),
    nrow(keys) == nrow(source), !anyDuplicated(keys$Uploaded_variation),
    setequal(keys$Uploaded_variation, source$Uploaded_variation),
    all(rows$lane %in% duckvep_replay_lanes),
    !anyDuplicated(rows[c("lane", "record_index", "alt_index", "Feature", "field")]))
  at <- match(duckvep_replay_indices(rows$record_index, nrow(source)), source_indices)
  feature <- keys$Feature[match(source$Uploaded_variation[at], keys$Uploaded_variation)]
  stopifnot(!anyNA(at), all(rows$alt_index == source$alt_index[at]),
    identical(as.character(rows$Feature), as.character(feature)))
  for (lane in duckvep_replay_lanes) {
    stopifnot(sum(rows$lane == lane) == summary$field_failures[summary$comparison == lane])
  }
}

duckvep_replay_inventory <- function(con, evidence) {
  parquet <- function(path) DBI::dbGetQuery(con, paste("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, path), ")"))
  cases <- utils::read.csv(file.path(evidence, "cases.csv"))
  summaries <- utils::read.csv(file.path(evidence, "summary.csv"))
  stopifnot(nrow(cases) > 0L, !anyDuplicated(cases$case), all(cases$prepared),
    nrow(summaries) == 3L * nrow(cases),
    nrow(utils::read.csv(file.path(evidence, "errors.csv"))) == 0L)
  tasks <- list()
  for (profile in cases$case) {
    directory <- file.path(evidence, profile)
    summary <- summaries[summaries$case == profile, ]
    duckvep_replay_summary(summary, profile, cases$input_records[cases$case == profile])
    source_rows <- parquet(file.path(directory, "source_keys.parquet"))
    failures <- do.call(rbind, lapply(duckvep_replay_lanes, function(lane) {
      x <- parquet(file.path(directory, paste0("comparison_", lane), "field_failures.parquet"))
      cbind(lane = rep(lane, nrow(x)), x)
    }))
    stopifnot(nrow(source_rows) == cases$input_records[cases$case == profile])
    keys <- utils::read.delim(file.path(directory, "duckvep_native_tab17.tsv.gz"),
      colClasses = "character", quote = "", na.strings = NULL, check.names = FALSE)
    duckvep_replay_failures(failures, source_rows, keys, summary)
    lines <- readLines(file.path(directory, "input.vcf"))
    model <- duckvep_evidence_sha256(file.path(directory, paste0(profile, ".gff3")))
    for (index in sort(unique(duckvep_replay_indices(failures$record_index, nrow(source_rows))))) {
      row <- source_rows[source_rows$record_index == index, , drop = FALSE]
      input <- duckvep_replay_record(lines, row, index)
      tasks[[length(tasks) + 1L]] <- list(profile = profile, source = row,
        input = c(input$header, input$record), failures = failures[failures$record_index == index, ], model = model)
    }
  }
  stopifnot(sum(vapply(tasks, function(t) nrow(t$failures), 0L)) == sum(summaries$field_failures))
  list(tasks = tasks, cases = cases, summary = summaries)
}

duckvep_replay_validate <- function(con, directory, profile, source, expected, input, model, shared, identity, status) {
  if (status != 1L) stop("singleton comparison must exit 1 with its retained disagreements")
  parquet <- function(path) DBI::dbGetQuery(con, paste("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, path), ")"))
  child <- file.path(directory, profile)
  errors <- utils::read.csv(file.path(directory, "errors.csv"))
  summary <- utils::read.csv(file.path(directory, "summary.csv"))
  stopifnot(nrow(errors) == 0L)
  duckvep_replay_summary(summary, profile, 1L)
  receipt <- utils::read.delim(file.path(directory, "receipt.tsv"), colClasses = "character")
  values <- stats::setNames(receipt$value, receipt$field)
  required <- c(identity, mode = "replay", requested_cases = "1", completed_cases = "1",
    input_records = "1", input_alleles = "1", completed_comparisons = "3",
    shared_inputs_unchanged = "TRUE", errors = "0")
  stopifnot(!anyDuplicated(receipt$field), identical(unname(values[names(required)]), unname(required)))
  observed_source <- parquet(file.path(child, "source_keys.parquet"))
  singleton_source <- source
  singleton_source$record_index <- 1
  stopifnot(identical(lapply(observed_source, as.character), lapply(singleton_source, as.character)),
    identical(readLines(file.path(child, "input.vcf")), input),
    identical(duckvep_evidence_sha256(file.path(child, paste0(profile, ".gff3"))), model))
  before <- utils::read.delim(file.path(directory, "inputs_before.tsv"), colClasses = "character")
  after <- utils::read.delim(file.path(directory, "inputs_after.tsv"), colClasses = "character")
  at <- match(shared$id, before$id)
  stopifnot(!anyDuplicated(before$id), !anyNA(at),
    identical(shared$sha256, before$sha256[at]), identical(before, after))
  observed <- do.call(rbind, lapply(duckvep_replay_lanes, function(lane) {
    rows <- parquet(file.path(child, paste0("comparison_", lane), "field_failures.parquet"))
    stopifnot(all(rows$record_index == 1L), all(rows$alt_index == 1L),
      nrow(rows) == summary$field_failures[summary$comparison == lane])
    rows$record_index <- rep(source$record_index, nrow(rows))
    cbind(lane = rep(lane, nrow(rows)), rows)
  }))
  passed <- identical(duckvep_replay_canonical(expected), duckvep_replay_canonical(observed))
  files <- c("duckvep_native_tab17.tsv", "duckvep_vep_csq.tsv", "fastvep_native_tab17.tsv",
    "fastvep_vep_csq.vcf", "vep.vcf")
  hashes <- stats::setNames(vapply(file.path(child, files), duckvep_evidence_sha256, ""), files)
  hashes <- c(hashes, receipt.tsv = duckvep_evidence_sha256(file.path(directory, "receipt.tsv")))
  list(passed = passed, failures = observed, summary = jsonlite::toJSON(summary, dataframe = "rows"),
    outputs = jsonlite::toJSON(as.list(hashes), auto_unbox = TRUE),
    error = if (passed) "" else "retained failure multiset changed")
}

duckvep_replay_resume <- function(retained, tasks, assigned, run_hash) {
  stopifnot(identical(names(retained), c(duckvep_replay_context, duckvep_replay_columns)),
    is.logical(retained$passed), !anyNA(retained$passed), all(retained$run_hash == run_hash))
  indices <- duckvep_replay_indices(retained$task_index, length(tasks))
  stopifnot(all(indices %in% assigned))
  for (index in unique(indices)) {
    rows <- retained[indices == index, , drop = FALSE]
    context <- unique(rows[duckvep_replay_context])
    task <- tasks[[index]]
    input_bytes <- charToRaw(paste0(paste(task$input, collapse = "\n"), "\n"))
    stopifnot(nrow(context) == 1L, context$profile == task$profile,
      as.character(context$original_record_index) == as.character(task$source$record_index),
      context$input_id == task$source$Uploaded_variation,
      context$header == paste(head(task$input, -1L), collapse = "\n"),
      context$record == tail(task$input, 1L), context$model_sha256 == task$model,
      context$input_sha256 == digest::digest(input_bytes, algo = "sha256", serialize = FALSE))
    if (!context$passed) {
      stopifnot(nzchar(context$error), dir.exists(context$retained_directory))
      next
    }
    stopifnot(context$status == 1L, context$error == "", context$retained_directory == "",
      identical(duckvep_replay_canonical(rows), duckvep_replay_canonical(task$failures)))
    summary <- jsonlite::fromJSON(context$summary)
    duckvep_replay_summary(summary, task$profile, 1L)
    for (lane in duckvep_replay_lanes) {
      stopifnot(sum(rows$lane == lane) == summary$field_failures[summary$comparison == lane])
    }
    hashes <- jsonlite::fromJSON(context$outputs)
    stopifnot(length(hashes) == 6L, all(lengths(hashes) == 1L), !anyDuplicated(names(hashes)),
      setequal(names(hashes), c("duckvep_native_tab17.tsv", "duckvep_vep_csq.tsv",
        "fastvep_native_tab17.tsv", "fastvep_vep_csq.vcf", "vep.vcf", "receipt.tsv")),
      all(grepl("^[0-9a-f]{64}$", unlist(hashes))))
  }
  unique(indices)
}

main <- function() {
  option <- optparse::make_option
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    option("--evidence", help = "published original field campaign to replay"),
    option("--output", default = ""), option("--extension-receipt", dest = "extension_receipt"),
    option("--fastvep-build-receipt", dest = "fastvep_build_receipt"), option("--fastvep", default = ""),
    option("--vep-prefix", dest = "vep_prefix", default = ""), option("--jobs", type = "integer", default = 1L),
    option("--cpus", default = "", help = "comma-separated Linux CPU IDs, one per worker"),
    option("--shard", type = "integer", default = 1L), option("--shards", type = "integer", default = 1L),
    option("--resume", action = "store_true", default = FALSE)
  )))
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"))
  source(file.path(root, "scripts/duckvep_evidence.R"))
  source(file.path(root, "r/duckhtsbench/R/registry.R"))
  source(file.path(root, "r/duckhtsbench/R/stage.R"))
  source(file.path(root, "r/duckhtsbench/R/duckvep.R"))
  source(file.path(root, "r/duckhtsbench/R/fastvep.R"))
  stopifnot(!is.null(opt$evidence), nzchar(opt$output), !is.null(opt$extension_receipt), !is.null(opt$fastvep_build_receipt),
    nzchar(opt$fastvep), nzchar(opt$vep_prefix), !is.na(opt$jobs), opt$jobs >= 1L, opt$jobs <= 64L,
    !is.na(opt$shard), !is.na(opt$shards), opt$shard >= 1L, opt$shard <= opt$shards,
    grepl("^[0-9]+(,[0-9]+)*$", opt$cpus), nzchar(Sys.which("taskset")))
  cpus <- strsplit(opt$cpus, ",", fixed = TRUE)[[1L]]
  stopifnot(length(cpus) == opt$jobs, !anyDuplicated(as.integer(cpus)), !anyNA(as.integer(cpus)))
  paths <- normalizePath(c(opt$evidence, opt$extension_receipt, opt$fastvep_build_receipt,
    opt$fastvep, opt$vep_prefix), mustWork = TRUE)
  evidence <- paths[[1L]]
  extension <- file.path(root, "build/release/duckhts.duckdb_extension")
  revision <- duckvep_evidence_revision(root)
  self <- file.path(root, "benchmarks/benchmark_duckvep_fastvep_replay.R")
  duckvep_evidence_command("git", c("-C", root, "ls-files", "--error-unmatch", self), "replay driver must be committed")
  duckvep_evidence_assert_checkout(root, revision, allowed_outputs = opt$output)
  ext <- duckvep_evidence_read_extension_receipt(paths[[2L]], root, extension, revision)
  # Cold staging precedes fork; workers only validate and reuse these files.
  duckhts_bench_stage_repository_fixtures(root, "duckvep-projection")
  pins <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b",
    fastvep = "18177c26a0d1d2419fe43c3e8f6d4a0b5c4a3eb6")
  build <- duckhts_bench_read_fastvep_build(paths[[3L]], pins[["fastvep"]], paths[[4L]])
  identity <- c(source_revision = revision, build_binding = ext$binding,
    extension_sha256 = duckvep_evidence_sha256(extension), fastvep_sha256 = duckvep_evidence_sha256(paths[[4L]]),
    fastvep_binding = build[["binding"]], pins)
  manifest <- utils::read.delim(file.path(evidence, "artifacts.tsv"), colClasses = "character")
  stopifnot(!anyDuplicated(manifest$path), !any(grepl("(^/|(^|/)\\.\\.(/|$))", manifest$path)))
  inputs <- c(file.path(evidence, manifest$path), file.path(evidence, "artifacts.tsv"), self,
    extension, paths[2:4], file.path(dirname(paths[[3L]]),
      c(build[["log"]], build[["source_tree"]], build[["source_commit_object"]])))
  hashes <- vapply(inputs, duckvep_evidence_sha256, "")
  stopifnot(identical(unname(hashes[seq_len(nrow(manifest))]), manifest$sha256))
  con <- DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  inventory <- duckvep_replay_inventory(con, evidence)
  DBI::dbDisconnect(con, shutdown = TRUE) # No DuckDB handle crosses fork().
  tasks <- inventory$tasks
  shared <- utils::read.delim(file.path(evidence, "inputs_before.tsv"), colClasses = "character")
  selected <- which((seq_along(tasks) - 1L) %% opt$shards == opt$shard - 1L)
  stopifnot(length(selected) > 0L)
  binding <- c(revision = revision, identity, driver_sha256 = duckvep_evidence_sha256(self),
    original_manifest_sha256 = duckvep_evidence_sha256(file.path(evidence, "artifacts.tsv")),
    input_identity = digest::digest(hashes, algo = "sha256"), vep_prefix = paths[[5L]],
    jobs = opt$jobs, cpus = opt$cpus, shard = opt$shard, shards = opt$shards,
    original_records = sum(inventory$cases$input_records), original_failure_cells = sum(inventory$summary$field_failures),
    required_singletons = length(tasks), assigned_singletons = length(selected))
  output <- normalizePath(opt$output, mustWork = FALSE)
  receipt <- file.path(output, "run.tsv")
  if (opt$resume) {
    previous <- utils::read.delim(receipt, colClasses = "character")
    stopifnot(identical(stats::setNames(previous$value, previous$field), binding))
  } else {
    stopifnot(!file.exists(output), dir.create(output, recursive = TRUE))
    utils::write.table(data.frame(field = names(binding), value = unname(binding)), receipt,
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  build_sources <- c(paths[[2L]], paths[[3L]],
    file.path(dirname(paths[[3L]]), c(build[["log"]], build[["source_tree"]], build[["source_commit_object"]])))
  build_targets <- file.path(output,
    c("extension_build.tsv", "fastvep_build.tsv", build[["log"]], build[["source_tree"]],
      build[["source_commit_object"]]))
  stopifnot(endsWith(build[["log"]], ".log"), !anyDuplicated(basename(build_targets)),
    !any(basename(build_targets) %in% c("run.tsv", "completion.csv")))
  if (!opt$resume) stopifnot(all(file.copy(build_sources, build_targets)))
  stopifnot(identical(unname(vapply(build_sources, duckvep_evidence_sha256, "")),
    unname(vapply(build_targets, duckvep_evidence_sha256, ""))))
  run_hash <- duckvep_evidence_sha256(receipt)
  worker <- function(worker_index) {
    con <- DBI::dbConnect(duckdb::duckdb(config = list(threads = "1", memory_limit = "256MB")))
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    checkpoint <- file.path(output, sprintf("worker-%02d.parquet", worker_index))
    retained <- if (file.exists(checkpoint)) {
      stopifnot(identical(readLines(paste0(checkpoint, ".sha256")), duckvep_evidence_sha256(checkpoint)))
      DBI::dbGetQuery(con, paste("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, checkpoint), ")"))
    } else NULL
    assigned <- selected[(seq_along(selected) - 1L) %% opt$jobs == worker_index - 1L]
    completed <- if (is.null(retained)) integer() else duckvep_replay_resume(retained, tasks, assigned, run_hash)
    for (task_index in setdiff(assigned, completed)) {
      task <- tasks[[task_index]]
      scratch <- tempfile(sprintf("replay-%05d-", task_index), tmpdir = output)
      stopifnot(dir.create(scratch))
      input <- file.path(scratch, "singleton.vcf")
      writeLines(task$input, input)
      child <- file.path(scratch, "result")
      args <- c("-c", cpus[[worker_index]], Sys.which("Rscript"),
        file.path(root, "benchmarks/benchmark_duckvep_fastvep_field_conformance.R"),
        "--case", task$profile, "--replay-input", input, "--output", child,
        "--extension", extension, "--extension-receipt", paths[[2L]], "--fastvep", paths[[4L]],
        "--fastvep-build-receipt", paths[[3L]], "--vep-prefix", paths[[5L]])
      status <- suppressWarnings(system2("taskset", shQuote(args),
        stdout = file.path(scratch, "process.log"), stderr = file.path(scratch, "process.log")))
      result <- tryCatch(duckvep_replay_validate(con, child, task$profile, task$source,
        task$failures, task$input, task$model, shared, identity, status), error = function(e) {
          list(passed = FALSE, failures = NULL, summary = "", outputs = "", error = conditionMessage(e))
        })
      rows <- result$failures
      if (is.null(rows) || !nrow(rows)) rows <- as.data.frame(stats::setNames(
        rep(list(NA_character_), length(duckvep_replay_columns)), duckvep_replay_columns))
      rows <- duckvep_replay_canonical(rows)
      context <- data.frame(run_hash, task_index, profile = task$profile,
        original_record_index = task$source$record_index, input_id = task$source$Uploaded_variation,
        header = paste(head(task$input, -1L), collapse = "\n"), record = tail(task$input, 1L),
        input_sha256 = duckvep_evidence_sha256(input), model_sha256 = task$model,
        status, passed = result$passed, summary = result$summary, outputs = result$outputs,
        error = result$error, retained_directory = if (result$passed) "" else scratch)
      retained <- rbind(retained, cbind(context[rep(1L, nrow(rows)), ], rows))
      DBI::dbWriteTable(con, "checkpoint", retained, overwrite = TRUE)
      temporary <- tempfile("checkpoint-", tmpdir = output, fileext = ".parquet")
      DBI::dbExecute(con, paste("COPY checkpoint TO", DBI::dbQuoteString(con, temporary), "(FORMAT PARQUET)"))
      stopifnot(file.rename(temporary, checkpoint))
      writeLines(duckvep_evidence_sha256(checkpoint), paste0(checkpoint, ".sha256"))
      if (result$passed) unlink(scratch, recursive = TRUE)
      cat(sprintf("worker %d: singleton %d/%d %s\n", worker_index, task_index, length(tasks),
        if (result$passed) "record-minimized" else "FAILED; artifacts retained"))
    }
    if (is.null(retained)) return(data.frame(task_index = integer(), passed = logical()))
    unique(retained[c("task_index", "passed")])
  }
  results <- parallel::mclapply(seq_len(opt$jobs), worker, mc.cores = opt$jobs)
  stopifnot(!any(vapply(results, inherits, logical(1L), "try-error")))
  results <- do.call(rbind, results)
  stopifnot(!anyDuplicated(results$task_index), setequal(results$task_index, selected),
    identical(hashes, vapply(inputs, duckvep_evidence_sha256, "")))
  duckvep_evidence_assert_checkout(root, revision, allowed_outputs = output)
  completion <- data.frame(shard = opt$shard, shards = opt$shards,
    required_singletons = length(tasks), assigned_singletons = length(selected),
    completed_singletons = nrow(results), passed_singletons = sum(results$passed),
    assigned_records_minimized = all(results$passed),
    all_records_minimized = length(selected) == length(tasks) && all(results$passed),
    allele_or_model_minimized = FALSE)
  utils::write.csv(completion, file.path(output, "completion.csv"), row.names = FALSE)
  print(completion)
  if (!all(results$passed)) quit(status = 1L)
}

if (sys.nframe() == 0L) main()
