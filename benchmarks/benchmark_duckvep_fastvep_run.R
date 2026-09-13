#!/usr/bin/env Rscript

# Paired fresh-process observations. All engines include startup, model loading
# and real-file output; CSQ includes flattening the upstream VCF transport.
duckvep_fastvep_run_job <- function(job, root) {
  run <- function(executable, args, env = character()) {
    status <- system2(executable, shQuote(args), env = env)
    if (status != 0L) stop("benchmark command failed: ", executable, " (", status, ")")
  }
  if (job$engine == "duckvep") {
    args <- c(file.path(root, "benchmarks/benchmark_duckvep_fastvep_worker.R"),
      "--extension", job$extension, "--model", job$model, "--input", job$input,
      "--output", job$output, "--output-contract", job$contract,
      "--threads", job$threads, "--distance", job$distance, "--memory-limit", job$memory_limit)
    if (job$contract != "operational17") args <- c(args, "--gff3", job$gff3)
    if (job$contract == "vep_csq") args <- c(args, "--fasta", job$fasta)
    run("Rscript", args)
  } else {
    output <- if (job$contract == "vep_csq") paste0(job$output, ".vcf") else job$output
    args <- c("annotate", "--input", job$input, "--output", output,
      "--output-format", if (job$contract == "vep_csq") "vcf" else "tab",
      "--transcript-cache", job$cache, "--distance", job$distance, "--no-progress")
    if (job$contract == "vep_csq") args <- c(args, "--hgvs", "--fasta", job$fasta)
    run(job$fastvep, args, paste0("RAYON_NUM_THREADS=", job$threads))
    if (job$contract == "vep_csq") {
      run("Rscript", c(file.path(root, "benchmarks/benchmark_duckvep_fastvep_extract.R"),
        "--input", output, "--output", job$output, "--threads", job$threads,
        "--memory-limit", job$memory_limit))
      unlink(output)
    }
  }
}

main <- function() {
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  parser <- optparse::OptionParser(option_list = list(
    optparse::make_option("--job", default = ""),
    optparse::make_option("--output", default = ""),
    optparse::make_option("--work-dir", dest = "work_dir", default = dirname(tempdir())),
    optparse::make_option("--input-id", dest = "input_id", default = "variantkey_giab_hg002_v421"),
    optparse::make_option("--extension", default = "build/release/duckhts.duckdb_extension"),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = ""),
    optparse::make_option("--fastvep", default = ".sync/fastVEP/target/release/fastvep"),
    optparse::make_option("--checkout", default = ".sync/fastVEP"),
    optparse::make_option("--affinity-one", dest = "affinity_one", default = "2"),
    optparse::make_option("--affinity-four", dest = "affinity_four", default = "2,4,6,8"),
    optparse::make_option("--repetitions", type = "integer", default = 3L),
    optparse::make_option("--diagnostic", action = "store_true", default = FALSE)))
  opt <- optparse::parse_args(parser)
  if (nzchar(opt$job)) return(duckvep_fastvep_run_job(readRDS(opt$job), root))
  if (!nzchar(opt$output) || file.exists(opt$output) || opt$repetitions < 1L) {
    stop("a new --output directory and positive repetitions are required")
  }
  for (count in c(1L, 4L)) {
    affinity <- if (count == 1L) opt$affinity_one else opt$affinity_four
    cpus <- strsplit(affinity, ",", fixed = TRUE)[[1L]]
    if (length(cpus) != count || anyDuplicated(cpus) || any(!grepl("^[0-9]+$", cpus))) {
      stop("affinity must name exactly ", count, " distinct CPU IDs")
    }
  }
  source(file.path(root, "scripts/duckvep_evidence.R"))
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"))
  for (name in c("registry", "stage", "duckvep", "fastvep")) {
    source(file.path(root, "r/duckhtsbench/R", paste0(name, ".R")))
  }
  Sys.setenv(DUCKHTSBENCH_REGISTRY = file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"))
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  if (!opt$diagnostic) {
    if (!nzchar(opt$extension_receipt)) stop("published observations require --extension-receipt")
    duckvep_evidence_assert_checkout(root, revision)
    duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)
    if (opt$repetitions < 3L) stop("published paired observations require at least three repetitions")
  }
  staged <- duckhts_bench_stage_fastvep(root, opt$checkout, opt$fastvep)
  inputs <- c(input = opt$input_id, model = "duckvep_ensembl116_model",
    fasta = "ensembl116_grch38_fasta_fa", gff3 = "ensembl116_grch38_gff3")
  paths <- vapply(inputs, duckhts_bench_artifact_path, character(1L))
  if (!all(file.exists(paths))) stop("all registered benchmark inputs must be staged first")
  for (id in inputs) duckhts_bench_validate_identity(id, duckhts_bench_artifact_path(id))
  registry <- duckhts_bench_registry()
  identity <- duckhts_bench_identity_fields(
    registry$supplier_identity[registry$id == inputs[["model"]]])
  duckhts_bench_validate_duckvep_ensembl116_model(paths[["model"]], extension,
    identity[["source_manifest_sha256"]])
  if (!dir.create(opt$output, recursive = TRUE)) stop("could not create a new output directory")
  work <- tempfile("duckvep-fastvep-", tmpdir = normalizePath(opt$work_dir, mustWork = TRUE))
  dir.create(work)
  message("Working files: ", work)
  # Keep a failed job and its log for diagnosis; remove successful output files
  # only after a full-row receipt has been written.
  source_files <- file.path(root, "benchmarks", paste0("benchmark_duckvep_fastvep_",
    c("run", "worker", "fields", "extract", "receipt"), ".R"))
  bound_files <- c(paths, fasta_index = paste0(paths[["fasta"]], ".fai"),
    extension = extension, fastvep = normalizePath(opt$fastvep),
    cache = staged[["cache"]], cache_receipt = staged[["receipt"]],
    stats::setNames(source_files, basename(source_files)))
  hashes <- vapply(bound_files, duckvep_evidence_sha256, character(1L))
  utils::write.csv(data.frame(artifact = names(bound_files), path = unname(bound_files), sha256 = hashes),
    file.path(opt$output, "inputs.csv"), row.names = FALSE)
  if (!file.copy(staged[["receipt"]], file.path(opt$output, "fastvep_cache_receipt.tsv"))) {
    stop("could not retain cache provenance")
  }
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
  counts <- DBI::dbGetQuery(con, paste0("SELECT count(*)::VARCHAR records,
    sum(len(ALT))::VARCHAR alt_alleles,
    sum(len(list_filter(ALT, a -> regexp_full_match(REF, '[ACGTNacgtn]+')
      AND regexp_full_match(a, '[ACGTNacgtn]+') AND upper(REF) <> upper(a))))::VARCHAR
      eligible_literal_alleles FROM read_bcf(", DBI::dbQuoteString(con, paths[["input"]]),
    ", scan_mode := 'sequential', decompression_threads := 0)"))
  DBI::dbDisconnect(con, shutdown = TRUE)
  if (is.na(counts$records) || as.numeric(counts$records) == 0 ||
      is.na(counts$eligible_literal_alleles) || as.numeric(counts$eligible_literal_alleles) == 0) {
    stop("benchmark input must contain records and eligible literal ALT alleles")
  }
  metadata <- c(source_revision = revision,
    binding = if (opt$diagnostic) "diagnostic_unbound" else "source_bound",
    fastvep_source_revision = duckvep_evidence_command("git", c("-C", opt$checkout,
      "rev-parse", "HEAD"), "cannot identify FastVEP source"),
    fastvep_version = duckvep_evidence_command(opt$fastvep, "--version", "cannot identify FastVEP executable"),
    input_id = opt$input_id, input_records = counts$records, input_alt_alleles = counts$alt_alleles,
    eligible_literal_alleles = counts$eligible_literal_alleles,
    r_version = R.version.string, duckdb_version = as.character(utils::packageVersion("duckdb")),
    distance = "5000", memory_limit = "4GB", supplementary_providers = "none",
    affinity_one = opt$affinity_one, affinity_four = opt$affinity_four)
  utils::write.csv(data.frame(field = names(metadata), value = unname(metadata)),
    file.path(opt$output, "metadata.csv"), row.names = FALSE)
  configurations <- data.frame(engine = c(rep("duckvep", 3), rep("fastvep", 2)),
    contract = c("operational17", "native_tab17", "vep_csq", "native_tab17", "vep_csq"))
  completed <- character()
  for (threads in c(1L, 4L)) for (run in seq_len(opt$repetitions)) {
    # Alternate tool order across repetitions to expose drift in a shared host.
    order <- if (run %% 2L) seq_len(nrow(configurations)) else rev(seq_len(nrow(configurations)))
    for (i in order) {
      job <- c(as.list(configurations[i, ]), as.list(paths), list(extension = extension,
        fastvep = normalizePath(opt$fastvep), cache = staged[["cache"]], threads = threads,
        distance = 5000L, memory_limit = "4GB", output = file.path(work, "output.tsv")))
      label <- paste(job$engine, job$contract, threads, run, sep = "_")
      message("Running ", label)
      job_path <- file.path(work, "job.rds")
      saveRDS(job, job_path)
      timing <- file.path(opt$output, paste0(label, ".time"))
      log <- file.path(opt$output, paste0(label, ".log"))
      affinity <- if (threads == 1L) opt$affinity_one else opt$affinity_four
      args <- c("-v", "-o", timing, "taskset", "-c", affinity, "Rscript",
        file.path(root, "benchmarks/benchmark_duckvep_fastvep_run.R"), "--job", job_path)
      status <- system2("/usr/bin/time", shQuote(args), stdout = log, stderr = log)
      if (status != 0L) stop("benchmark failed; files retained at ", work, "; log: ", log)
      skip <- duckvep_fastvep_tab_header(job$output, duckvep_fastvep_fields(job$contract))
      receipt <- file.path(opt$output, paste0(label, ".csv"))
      args <- c(file.path(root, "benchmarks/benchmark_duckvep_fastvep_receipt.R"),
        "--input", job$output, "--tool", job$engine, "--output-contract", job$contract,
        "--threads", threads, "--run", run, "--timing-file", timing, "--output", receipt)
      args <- c(args, "--skip-lines", skip)
      if (system2("Rscript", shQuote(args)) != 0L) stop("receipt failed; output retained at ", work)
      result <- utils::read.csv(receipt, colClasses = "character")
      if (nrow(result) != 1L || is.na(result$row_count) || as.numeric(result$row_count) == 0) {
        stop("empty or invalid output receipt; output retained at ", work)
      }
      completed <- c(completed, receipt)
      unlink(job$output)
    }
  }
  if (!identical(hashes, vapply(bound_files, duckvep_evidence_sha256, character(1L)))) {
    stop("benchmark source, executable or input changed during measurement")
  }
  if (!opt$diagnostic) duckvep_evidence_assert_checkout(root, revision)
  duckhts_bench_fastvep_source(opt$checkout, metadata[["fastvep_source_revision"]])
  stopifnot(length(completed) == nrow(configurations) * 2L * opt$repetitions)
  observations <- do.call(rbind, lapply(completed, utils::read.csv, colClasses = "character"))
  groups <- split(observations, paste(observations$tool, observations$output_contract))
  fingerprints <- c("row_count", "xor_hash", "low32_sum", "high32_sum")
  for (name in names(groups)) {
    if (nrow(unique(groups[[name]][fingerprints])) != 1L) {
      stop("output multiset changes across repeats or thread counts: ", name)
    }
  }
  artifacts <- list.files(opt$output, full.names = TRUE)
  utils::write.csv(data.frame(file = basename(artifacts),
    sha256 = vapply(artifacts, duckvep_evidence_sha256, character(1L))),
    file.path(opt$output, "artifacts.csv"), row.names = FALSE)
  utils::write.csv(data.frame(source_revision = revision, completed = TRUE,
    binding = metadata[["binding"]], observations = length(completed), repetitions = opt$repetitions,
    artifacts_sha256 = duckvep_evidence_sha256(file.path(opt$output, "artifacts.csv"))),
    file.path(opt$output, "completion.csv"), row.names = FALSE)
  unlink(work, recursive = TRUE)
  message("Paired observations complete: ", opt$output)
}

if (sys.nframe() == 0L) main()
