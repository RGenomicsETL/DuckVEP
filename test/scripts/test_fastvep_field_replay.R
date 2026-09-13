#!/usr/bin/env Rscript
# Network-free controls for exact record identity and failure-multiset replay.
local({
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_replay.R"))
  source(file.path(root, "scripts/duckvep_evidence.R"))
  directory <- tempfile("fastvep-replay-test-")
  stopifnot(dir.create(directory))
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(config = list(threads = "1")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  rejects <- function(expression) {
    result <- tryCatch(force(expression), error = function(e) FALSE)
    stopifnot(identical(result, FALSE) || (is.list(result) && identical(result$passed, FALSE)))
  }
  parquet <- function(x, path) {
    if (file.exists(path)) unlink(path)
    DBI::dbWriteTable(con, "fixture", x, overwrite = TRUE)
    DBI::dbExecute(con, paste("COPY fixture TO", DBI::dbQuoteString(con, path), "(FORMAT PARQUET)"))
  }
  table <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  source_row <- data.frame(record_index = 7L, alt_index = 1L, Uploaded_variation = "case_7",
    CHROM = "chr1", POS = 10L, REF = "A", ALT = "G")
  header <- c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  record <- "chr1\t10\tcase_7\tA\tG\t.\tPASS\t."
  lines <- c(header, rep("chr1\t9\tother\tC\tT\t.\tPASS\t.", 6L), record)
  stopifnot(identical(duckvep_replay_record(lines, source_row, 7L), list(header = header, record = record)))
  text_source <- source_row
  text_source$record_index <- "7"
  text_source$alt_index <- "1"
  stopifnot(identical(duckvep_replay_record(lines, text_source, "7"), list(header = header, record = record)),
    identical(sort(duckvep_replay_indices(c("11", "2", "10"), 11L)), c(2L, 10L, 11L)))
  for (index in c("0", "-1", "1.5", "1e0", "NA", NA, "2147483648")) {
    rejects(duckvep_replay_indices(index, .Machine$integer.max))
  }
  rejects(duckvep_replay_record(lines, source_row, 6L))
  for (column in c("Uploaded_variation", "CHROM", "POS", "REF", "ALT", "alt_index")) {
    wrong <- source_row
    wrong[[column]] <- "changed"
    rejects(duckvep_replay_record(lines, wrong, 7L))
  }
  profile <- "forward"
  child <- file.path(directory, profile)
  stopifnot(dir.create(child))
  writeLines(c(header, record), file.path(child, "input.vcf"))
  model_path <- file.path(child, "forward.gff3")
  writeLines("##gff-version 3", model_path)
  model <- duckvep_evidence_sha256(model_path)
  singleton <- source_row
  singleton$record_index <- 1L
  parquet(singleton, file.path(child, "source_keys.parquet"))
  errors <- data.frame(case = character(), stage = character(), message = character())
  utils::write.csv(errors, file.path(directory, "errors.csv"), row.names = FALSE)
  summary <- data.frame(case = profile, comparison = duckvep_replay_lanes,
    input_records = 1L, input_alleles = 1L, actual_rows = 1L, expected_rows = 1L, union_keys = 1L,
    missing_keys = 0L, extra_keys = 0L, actual_duplicate_keys = 0L, expected_duplicate_keys = 0L,
    invalid_keys = 0L, compared_keys = 1L, source_rows = 1L, source_alleles = 1L,
    source_duplicate_alleles = 0L, source_invalid_alleles = 0L, actual_missing_source_alleles = 0L,
    expected_missing_source_alleles = 0L, actual_unknown_alleles = 0L, expected_unknown_alleles = 0L,
    compared_fields = c(31L, 31L, 16L), field_failures = c(1L, 2L, 1L), passed = FALSE)
  utils::write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  shared <- data.frame(id = "reference", path = model_path, sha256 = model)
  table(shared, file.path(directory, "inputs_before.tsv"))
  table(shared, file.path(directory, "inputs_after.tsv"))
  identity <- c(source_revision = strrep("a", 40L), build_binding = "htslib_distclean_make_release",
    extension_sha256 = strrep("b", 64L), fastvep_sha256 = strrep("c", 64L),
    fastvep_binding = "cargo_verified_commit_tree_release_locked_offline")
  receipt <- c(identity, mode = "replay", requested_cases = "1", completed_cases = "1",
    input_records = "1", input_alleles = "1", completed_comparisons = "3", shared_inputs_unchanged = "TRUE", errors = "0")
  table(data.frame(field = names(receipt), value = unname(receipt)), file.path(directory, "receipt.tsv"))
  raw_outputs <- c("duckvep_native_tab17.tsv", "duckvep_vep_csq.tsv", "fastvep_native_tab17.tsv",
    "fastvep_vep_csq.vcf", "vep.vcf")
  for (name in raw_outputs) writeLines(name, file.path(child, name))
  expected <- data.frame(lane = duckvep_replay_lanes[c(1L, 2L, 2L, 3L)], record_index = 7L,
    alt_index = 1L, Feature = "TX1", field = c("HGVSc", "HGVSc", "HGVSp", "Consequence"),
    actual = c("a", "", NA, "-"), expected = c("b", "-", "", "missense_variant"))
  original_sources <- source_row[rep(1L, 7L), ]
  original_sources$record_index <- seq_len(7L)
  original_sources$Uploaded_variation <- paste0("case_", seq_len(7L))
  original_keys <- data.frame(Uploaded_variation = original_sources$Uploaded_variation, Feature = "TX1")
  duckvep_replay_failures(expected, original_sources, original_keys, summary)
  rejects(duckvep_replay_failures(expected[-1L, ], original_sources, original_keys, summary))
  for (column in c("record_index", "alt_index", "Feature")) {
    wrong <- expected
    wrong[[column]][[1L]] <- if (column == "Feature") "unknown_feature" else 8L
    rejects(duckvep_replay_failures(wrong, original_sources, original_keys, summary))
  }
  for (lane in duckvep_replay_lanes) {
    stopifnot(dir.create(file.path(child, paste0("comparison_", lane))))
    rows <- expected[expected$lane == lane, -1L]
    rows$record_index <- 1L
    parquet(rows, file.path(child, paste0("comparison_", lane), "field_failures.parquet"))
  }
  validate <- function(status = 1L) duckvep_replay_validate(con, directory, profile, source_row,
    expected, c(header, record), model, shared, identity, status)
  result <- validate()
  stopifnot(result$passed, identical(duckvep_replay_canonical(result$failures), duckvep_replay_canonical(expected)),
    identical(duckvep_replay_canonical(expected[4:1, ]), duckvep_replay_canonical(expected)))
  rejects(validate(0L))
  rejects(validate(2L))
  # Matched row counts do not excuse changed, dropped or duplicated field cells.
  path <- file.path(child, "comparison_fastvep_vep_csq", "field_failures.parquet")
  rows <- expected[expected$lane == "fastvep_vep_csq", -1L]
  rows$record_index <- 1L
  changed <- rows
  changed$actual[[1L]] <- "changed"
  changed_key <- rows
  changed_key$field[[1L]] <- "Consequence"
  for (mutation in list(changed, changed_key, rows[1L, ], rbind(rows, rows[1L, ]))) {
    parquet(mutation, path)
    rejects(validate())
  }
  parquet(rows, path)
  for (column in c("compared_keys", "actual_rows", "compared_fields", "missing_keys", "actual_unknown_alleles")) {
    wrong <- summary
    wrong[[column]][[1L]] <- wrong[[column]][[1L]] + 1L
    utils::write.csv(wrong, file.path(directory, "summary.csv"), row.names = FALSE)
    rejects(validate())
  }
  utils::write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  utils::write.csv(data.frame(case = profile, stage = "VEP", message = "injected child failure"),
    file.path(directory, "errors.csv"), row.names = FALSE)
  rejects(validate())
  utils::write.csv(errors, file.path(directory, "errors.csv"), row.names = FALSE)
  wrong <- singleton
  wrong$record_index <- 7L
  parquet(wrong, file.path(child, "source_keys.parquet"))
  rejects(validate())
  parquet(singleton, file.path(child, "source_keys.parquet"))
  wrong <- shared
  wrong$sha256 <- strrep("d", 64L)
  table(wrong, file.path(directory, "inputs_before.tsv"))
  table(wrong, file.path(directory, "inputs_after.tsv"))
  rejects(validate())
  table(shared, file.path(directory, "inputs_before.tsv"))
  table(shared, file.path(directory, "inputs_after.tsv"))
  wrong_receipt <- receipt
  wrong_receipt[["extension_sha256"]] <- strrep("d", 64L)
  table(data.frame(field = names(wrong_receipt), value = unname(wrong_receipt)), file.path(directory, "receipt.tsv"))
  rejects(validate())
  table(data.frame(field = names(receipt), value = unname(receipt)), file.path(directory, "receipt.tsv"))
  writeLines(c(header, sub("case_7", "case_8", record)), file.path(child, "input.vcf"))
  rejects(validate())
  writeLines(c(header, record), file.path(child, "input.vcf"))
  stopifnot(validate()$passed)
  task <- list(profile = profile, source = source_row, input = c(header, record), failures = expected, model = model)
  run_hash <- strrep("e", 64L)
  context <- data.frame(run_hash, task_index = 1L, profile, original_record_index = 7L, input_id = "case_7",
    header = paste(header, collapse = "\n"), record, input_sha256 = duckvep_evidence_sha256(file.path(child, "input.vcf")),
    model_sha256 = model, status = 1L, passed = TRUE, summary = result$summary, outputs = result$outputs,
    error = "", retained_directory = "")
  retained <- cbind(context[rep(1L, nrow(expected)), ], expected)
  checkpoint <- file.path(directory, "checkpoint.parquet")
  parquet(retained, checkpoint)
  retained <- DBI::dbGetQuery(con, paste("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, checkpoint), ")"))
  resume <- function(rows) duckvep_replay_resume(rows, list(task), 1L, run_hash)
  stopifnot(identical(resume(retained), 1L))
  # A rewritten checkpoint and matching checksum cannot certify lost cells or
  # a different record. Resume rechecks the content against the original pack.
  changed <- retained
  changed$actual[[1L]] <- "changed"
  for (mutation in list(changed, retained[-1L, ], rbind(retained, retained[1L, ]))) rejects(resume(mutation))
  for (column in c("run_hash", "profile", "original_record_index", "input_id", "header", "record",
    "model_sha256", "input_sha256", "summary", "outputs")) {
    changed <- retained
    changed[[column]] <- "changed"
    rejects(resume(changed))
  }
  changed <- retained
  changed$record[[1L]] <- "conflicting context"
  rejects(resume(changed))
  rejects(resume(retained[, -1L]))
  rejects(duckvep_replay_resume(retained, list(task), integer(), run_hash))
  changed <- retained
  changed$passed <- FALSE
  changed$error <- "injected child failure"
  changed$retained_directory <- directory
  stopifnot(identical(resume(changed), 1L), !any(changed$passed))
  changed$retained_directory <- file.path(directory, "missing-failure-evidence")
  rejects(resume(changed))
  rmd <- readLines(file.path(root, "benchmarks/benchmark_duckvep_fastvep.Rmd"))
  begin <- which(rmd == "```{r complete-field-singletons}")
  stopifnot(length(begin) == 1L)
  end <- which(seq_along(rmd) > begin & rmd == "```")[[1L]]
  gate <- parse(text = rmd[seq.int(begin + 1L, end - 1L)])
  definition <- Filter(function(x) is.call(x) && identical(x[[1L]], quote(`<-`)) &&
    identical(x[[2L]], quote(validate_replay_bundle)), gate)
  stopifnot(length(definition) == 1L)
  eval(definition[[1L]])
  evidence <- file.path(directory, "original-pack")
  bundle <- file.path(directory, "replay-pack")
  stopifnot(dir.create(evidence), dir.create(bundle), dir.create(file.path(evidence, profile)))
  original <- file.path(evidence, profile)
  parquet(original_sources, file.path(original, "source_keys.parquet"))
  writeLines(c(header, paste0("chr1\t10\tcase_", 1:7, "\tA\tG\t.\tPASS\t.")), file.path(original, "input.vcf"))
  stopifnot(file.copy(model_path, file.path(original, "forward.gff3")))
  compressed <- gzfile(file.path(original, "duckvep_native_tab17.tsv.gz"), "wt")
  table(original_keys, compressed)
  close(compressed)
  original_summary <- summary
  counts <- c("input_records", "input_alleles", "actual_rows", "expected_rows", "union_keys", "compared_keys", "source_rows", "source_alleles")
  original_summary[counts] <- 7L
  utils::write.csv(original_summary, file.path(evidence, "summary.csv"), row.names = FALSE)
  utils::write.csv(data.frame(case = profile, input_records = 7L, input_alleles = 7L, prepared = TRUE),
    file.path(evidence, "cases.csv"), row.names = FALSE)
  utils::write.csv(errors, file.path(evidence, "errors.csv"), row.names = FALSE)
  for (lane in duckvep_replay_lanes) {
    stopifnot(dir.create(file.path(original, paste0("comparison_", lane))))
    parquet(expected[expected$lane == lane, -1L], file.path(original, paste0("comparison_", lane), "field_failures.parquet"))
  }
  inventory <- duckvep_replay_inventory(con, evidence)
  stopifnot(length(inventory$tasks) == 1L, sum(inventory$cases$input_records) == 7L,
    sum(inventory$summary$field_failures) == 4L,
    identical(duckvep_replay_canonical(inventory$tasks[[1L]]$failures), duckvep_replay_canonical(expected)))
  wrong_original <- original_summary
  wrong_original$field_failures[[1L]] <- 0L
  utils::write.csv(wrong_original, file.path(evidence, "summary.csv"), row.names = FALSE)
  rejects(duckvep_replay_inventory(con, evidence))
  utils::write.csv(original_summary, file.path(evidence, "summary.csv"), row.names = FALSE)
  artifacts <- list.files(evidence, recursive = TRUE)
  table(data.frame(path = artifacts, sha256 = unname(vapply(file.path(evidence, artifacts), duckvep_evidence_sha256, ""))),
    file.path(evidence, "artifacts.tsv"))
  write_fields <- function(x, file) table(data.frame(field = names(x), value = unname(x)), file.path(bundle, file))
  pins <- c(vep = "57ea5c52340acc1f156267f810ad162e26597082", variation = "2fb834b987ede3824e200197a838ce11e91aeb4b",
    fastvep = "18177c26a0d1d2419fe43c3e8f6d4a0b5c4a3eb6")
  writeLines("synthetic completed build", file.path(bundle, "build.log"))
  source_proof <- file.path(root, "benchmarks/data/duckvep_fastvep/fastvep_verified_commit_tree")
  stopifnot(all(file.copy(file.path(source_proof, c("source-tree.txt", "source-commit.bin")), bundle)))
  write_fields(c(binding = "cargo_verified_commit_tree_release_locked_offline", source_commit = pins[["fastvep"]],
    cargo_lock_sha256 = strrep("f", 64L), toolchain = "1.98.1", rustc = "synthetic", cargo = "synthetic",
    rustflags = "-C target-cpu=native", command = "cargo build --release --locked --offline",
    executable_sha256 = identity[["fastvep_sha256"]], log = "build.log",
    log_sha256 = duckvep_evidence_sha256(file.path(bundle, "build.log")), exit_status = "0",
    source_tree = "source-tree.txt", source_tree_sha256 = duckvep_evidence_sha256(file.path(bundle, "source-tree.txt")),
    source_commit_object = "source-commit.bin"),
    "fastvep_build.tsv")
  write_fields(c(source_revision = identity[["source_revision"]], path = "/synthetic/build/release/duckhts.duckdb_extension",
    binding = identity[["build_binding"]], sha256 = identity[["extension_sha256"]]), "extension_build.tsv")
  run <- c(revision = identity[["source_revision"]], identity, pins,
    driver_sha256 = strrep("f", 64L), input_identity = strrep("f", 64L),
    original_manifest_sha256 = duckvep_evidence_sha256(file.path(evidence, "artifacts.tsv")),
    jobs = "1", shard = "1", shards = "1", original_records = "7", original_failure_cells = "4",
    required_singletons = "1", assigned_singletons = "1")
  write_fields(run, "run.tsv")
  retained$run_hash <- duckvep_evidence_sha256(file.path(bundle, "run.tsv"))
  seal_checkpoint <- function(rows) {
    path <- file.path(bundle, "worker-01.parquet")
    parquet(rows, path)
    writeLines(duckvep_evidence_sha256(path), paste0(path, ".sha256"))
  }
  seal_checkpoint(retained)
  completion <- data.frame(shard = 1L, shards = 1L, required_singletons = 1L, assigned_singletons = 1L,
    completed_singletons = 1L, passed_singletons = 1L, assigned_records_minimized = TRUE,
    all_records_minimized = TRUE, allele_or_model_minimized = FALSE)
  utils::write.csv(completion, file.path(bundle, "completion.csv"), row.names = FALSE)
  report <- function() validate_replay_bundle(bundle, evidence, c(records = 7L, tasks = 1L, cells = 4L))
  stopifnot(report()$retained_failure_cells == 4L)
  for (mutation in list(retained[-1L, ], rbind(retained, retained[1L, ]), transform(retained, actual = "changed"))) {
    seal_checkpoint(mutation)
    rejects(report())
  }
  seal_checkpoint(retained)
  wrong_completion <- completion
  wrong_completion$shards <- 2L
  utils::write.csv(wrong_completion, file.path(bundle, "completion.csv"), row.names = FALSE)
  rejects(report())
  utils::write.csv(completion, file.path(bundle, "completion.csv"), row.names = FALSE)
  for (field in c("original_manifest_sha256", "extension_sha256", "fastvep_sha256")) {
    wrong_run <- run
    wrong_run[[field]] <- strrep("0", 64L)
    write_fields(wrong_run, "run.tsv")
    resealed <- retained
    resealed$run_hash <- duckvep_evidence_sha256(file.path(bundle, "run.tsv"))
    seal_checkpoint(resealed)
    rejects(report())
  }
  write_fields(run, "run.tsv")
  seal_checkpoint(retained)
  stopifnot(report()$changed_or_missing_cells == 0L)
  proof <- file.path(dirname(bundle), "fastvep_verified_commit_tree")
  stopifnot(dir.create(proof))
  stopifnot(file.copy(file.path(bundle, "fastvep_build.tsv"), file.path(proof, "build.tsv")),
    all(file.copy(file.path(bundle, c("build.log", "source-tree.txt", "source-commit.bin")), proof)))
  verified <- read.delim(file.path(bundle, "fastvep_build.tsv"), colClasses = "character")
  for (binding in c("cargo_fresh_release_locked_offline", "cargo_verified_tree_release_locked_offline")) {
    removed <- c("source_commit_object", if (binding == "cargo_fresh_release_locked_offline")
      c("source_tree", "source_tree_sha256"))
    measured <- verified[!verified$field %in% removed, ]
    measured$value[measured$field == "binding"] <- binding
    write_fields(setNames(measured$value, measured$field), "fastvep_build.tsv")
    run[["fastvep_binding"]] <- binding
    write_fields(run, "run.tsv")
    retained$run_hash <- duckvep_evidence_sha256(file.path(bundle, "run.tsv"))
    seal_checkpoint(retained)
    recorded <- duckvep_evidence_sha256(file.path(bundle, "fastvep_build.tsv"))
    stopifnot(report()$changed_or_missing_cells == 0L,
      identical(duckvep_evidence_sha256(file.path(bundle, "fastvep_build.tsv")), recorded))
  }
  unlink(file.path(proof, "build.tsv"))
  rejects(suppressWarnings(report()))
  cat("FastVEP singleton replay and report identity, denominator and injected failure controls: OK\n")
})
