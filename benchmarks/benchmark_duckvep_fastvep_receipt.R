#!/usr/bin/env Rscript

# Produce one output-file receipt for the DuckVEP/FastVEP comparison. The
# fingerprint is independent of row order but remains sensitive to every
# projected column. Callers retain one row per timed or diagnostic run.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})

op <- OptionParser()
op <- add_option(op, "--input", default = "")
op <- add_option(op, "--tool", default = "")
op <- add_option(
  op,
  "--output-contract",
  dest = "output_contract",
  default = ""
)
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(op, "--run", type = "integer", default = 1L)
op <- add_option(op, "--timing-file", dest = "timing_file", default = "")
op <- add_option(
  op,
  "--skip-lines",
  dest = "skip_lines",
  type = "integer",
  default = 0L
)
op <- add_option(op, "--output", default = "")
op <- add_option(op, "--source-map", dest = "source_map", default = "")
op <- add_option(op, "--source-sha256", dest = "source_sha256", default = "")
op <- add_option(op, "--coverage-output", dest = "coverage_output", default = "")
op <- add_option(op, "--memory-limit", dest = "memory_limit", default = "4GB")
op <- add_option(op, "--max-spill", dest = "max_spill", default = "8GB")
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (!nzchar(opt$input) || !file.exists(opt$input)) {
  die("--input must name an existing output file")
}
if (!nzchar(opt$tool) || !nzchar(opt$output_contract) || !nzchar(opt$output)) {
  die("--tool, --output-contract, and --output are required")
}
if (opt$threads < 1L || opt$run < 1L || opt$skip_lines < 0L) {
  die("--threads and --run must be positive; --skip-lines must be non-negative")
}
coverage_options <- nzchar(c(opt$source_map, opt$source_sha256, opt$coverage_output))
if (any(coverage_options) && !all(coverage_options)) {
  die("--source-map, --source-sha256 and --coverage-output are required together")
}
check_coverage <- all(coverage_options)
if (check_coverage && (!file.exists(opt$source_map) || !grepl("^[0-9a-f]{64}$", opt$source_sha256))) {
  die("source map must exist and source SHA256 must be a lowercase digest")
}

input <- normalizePath(opt$input)
drv <- duckdb(dbdir = ":memory:")
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
input_sql <- as.character(dbQuoteString(con, input))
dbExecute(con, paste("SET threads =", opt$threads))
dbExecute(con, paste("SET memory_limit =", dbQuoteString(con, opt$memory_limit)))
dbExecute(con, paste("SET max_temp_directory_size =", dbQuoteString(con, opt$max_spill)))
dbExecute(con, paste("SET temp_directory =", dbQuoteString(con, file.path(tempdir(), "duckdb"))))
if (check_coverage) {
  script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")])
  source(file.path(dirname(script), "benchmark_duckvep_fastvep_fields.R"))
  relation <- duckvep_fastvep_tab_relation(con, input,
    duckvep_fastvep_transport_fields(opt$output_contract))
} else {
  relation <- glue("read_csv({input_sql}, delim = '\\t', header = true,
    skip = {opt$skip_lines}, quote = '', escape = '', all_varchar = true)")
}
dbExecute(con, paste("CREATE TEMP VIEW fastvep_receipt_output AS SELECT * FROM", relation))

fingerprint <- dbGetQuery(
  con,
  glue(
    "WITH rows AS (
       SELECT hash(row(*COLUMNS(*)))::UBIGINT AS h
       FROM fastvep_receipt_output
     ), receipt AS (
       SELECT
         count(*)::UBIGINT AS row_count,
         bit_xor(h)::UBIGINT AS xor_hash,
         sum((h & 4294967295::UBIGINT)::HUGEINT)::HUGEINT AS low32_sum,
         sum((h >> 32)::HUGEINT)::HUGEINT AS high32_sum
       FROM rows
     )
     SELECT
       row_count::VARCHAR AS row_count,
       xor_hash::VARCHAR AS xor_hash,
       low32_sum::VARCHAR AS low32_sum,
       high32_sum::VARCHAR AS high32_sum
     FROM receipt"
  )
)

file_sha256 <- function(path) {
  result <- system2("sha256sum", shQuote(path), stdout = TRUE)
  value <- strsplit(result, " +")[[1L]][1L]
  if (!is.null(attr(result, "status")) || !grepl("^[0-9a-f]{64}$", value)) {
    die("cannot hash artifact: {path}")
  }
  value
}
sha256 <- file_sha256(input)
lines <- strsplit(system2("wc", c("-l", shQuote(input)), stdout = TRUE), " +")[[1L]]
lines <- lines[nzchar(lines)][[1L]]
timing_file <- if (nzchar(opt$timing_file)) basename(opt$timing_file) else ""

receipt <- data.frame(
  tool = opt$tool,
  output_contract = opt$output_contract,
  threads = opt$threads,
  run = opt$run,
  timing_file = timing_file,
  row_count = fingerprint$row_count,
  bytes = as.character(file.info(input)$size),
  lines = lines,
  sha256 = sha256,
  multiset_checked = TRUE,
  fingerprint_scope = "full_row",
  xor_hash = fingerprint$xor_hash,
  low32_sum = fingerprint$low32_sum,
  high32_sum = fingerprint$high32_sum,
  stringsAsFactors = FALSE
)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(receipt, opt$output, row.names = FALSE, quote = TRUE)
if (check_coverage) {
  dir.create(dirname(opt$coverage_output), recursive = TRUE, showWarnings = FALSE)
  failures <- paste0(sub("\\.csv$", "", opt$coverage_output), ".failures.parquet")
  coverage <- duckvep_fastvep_source_coverage(con, "fastvep_receipt_output",
    opt$source_map, opt$output_contract, failures)
  if (as.character(coverage$output_rows) != fingerprint$row_count) {
    die("source coverage changed the final output row denominator")
  }
  coverage <- cbind(data.frame(tool = opt$tool, output_contract = opt$output_contract,
    threads = opt$threads, run = opt$run, scope = "final_output",
    input_sha256 = opt$source_sha256, output_sha256 = sha256,
    source_map_sha256 = file_sha256(opt$source_map),
    failures_file = basename(failures), failures_sha256 = file_sha256(failures)), coverage)
  utils::write.csv(coverage, opt$coverage_output, row.names = FALSE, quote = TRUE)
  if (!coverage$passed) die("final output failed physical source-ALT coverage; see {opt$coverage_output}")
}
