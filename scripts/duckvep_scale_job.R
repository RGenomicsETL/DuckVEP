#!/usr/bin/env Rscript
# One scale job: load the model under the native ceilings, then annotate a panel through the
# public builders and stream each result to Parquet. Normally started by duckvep_scale_run.sh
# inside a cgroup; it can also be run directly.
#
#   Rscript scripts/duckvep_scale_job.R --panel PANEL.parquet --out JOBDIR [options]
#
# Options (defaults are the agreed ceilings):
#   --modes compact,complete17   --threads 6          --model PATH        --extension PATH
#   --memory-limit 8GB           --max-temp 32GiB     --native-budget-mib 4096
#   --workers 6                  --scratch-mib 128    --emit-mib 256
#   --limit-rows N               --keep-output        --retry-after-capacity   --job-id ID
#   --regulation (default) | --no-regulation   resident RegulatoryFeature and MotifFeature intervals
#   --temp-dir DIR (DuckDB spill; default JOBDIR/spill)
#
# Writes JOBDIR/result.tsv (one row per mode) and never leaves a partial Parquet behind. The
# outcome is ok, capacity_error (an explicit native or DuckDB capacity error) or failed.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })

parse_args <- function(argv) {
  flags <- c("keep-output", "retry-after-capacity", "regulation", "no-regulation")
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a)
    key <- substring(a, 3L)
    if (key %in% flags) { out[[key]] <- TRUE; i <- i + 1L; next }
    if (i == length(argv)) stop("missing value for ", a)
    out[[key]] <- argv[[i + 1L]]
    i <- i + 2L
  }
  out
}
opt <- parse_args(commandArgs(TRUE))
get_opt <- function(key, default = NULL) if (is.null(opt[[key]])) default else opt[[key]]
script_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "scripts"
})
repo <- normalizePath(file.path(script_dir, ".."))
panel <- get_opt("panel"); job_dir <- get_opt("out")
if (is.null(panel) || is.null(job_dir)) stop("--panel and --out are required")
modes <- strsplit(get_opt("modes", "compact,complete17"), ",", fixed = TRUE)[[1L]]
if (!all(modes %in% c("compact", "complete17"))) stop("--modes: compact and/or complete17")
job_id <- get_opt("job-id", "1")
threads <- as.integer(get_opt("threads", "6"))
model <- get_opt("model", Sys.getenv("DUCKVEP_SCALE_MODEL",
  "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb"))
extension <- get_opt("extension", file.path(repo, "build/release/extension/duckvep/duckvep.duckdb_extension"))
memory_limit <- get_opt("memory-limit", "8GB")
max_temp <- get_opt("max-temp", "32GiB")
budget_mib <- as.numeric(get_opt("native-budget-mib", "4096"))
workers <- as.integer(get_opt("workers", "6"))
scratch_mib <- as.numeric(get_opt("scratch-mib", "128"))
emit_mib <- as.numeric(get_opt("emit-mib", "256"))
limit_rows <- get_opt("limit-rows")
keep_output <- isTRUE(opt[["keep-output"]])
retry <- isTRUE(opt[["retry-after-capacity"]])
regulation <- !isTRUE(opt[["no-regulation"]])
for (f in c(panel, model, extension)) if (!file.exists(f)) stop("missing file: ", f)
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
spill <- get_opt("temp-dir", file.path(job_dir, "spill"))
dir.create(spill, showWarnings = FALSE, recursive = TRUE)
result_path <- file.path(job_dir, "result.tsv")
unlink(result_path)

canonical <- file.path(repo, "benchmarks/data/scale_contracts/canonical-metadata.parquet")
source(file.path(repo, "benchmarks/duckvep_field_projection.R"), local = TRUE)

con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
q <- function(x) as.character(dbQuoteString(con, x))
run <- function(sql) invisible(dbExecute(con, sql))
get <- function(sql) dbGetQuery(con, sql)
status_kib <- function(key) {
  line <- grep(paste0("^", key, ":"), readLines("/proc/self/status", warn = FALSE), value = TRUE)
  as.numeric(sub(paste0("^", key, ":[[:space:]]*([0-9]+).*"), "\\1", line))
}
classify <- function(msg) {
  if (grepl("capacity error", msg, ignore.case = TRUE) ||
      grepl("Out of Memory Error|max_temp_directory_size|failed to offload|Could not allocate|not enough memory",
        msg, ignore.case = TRUE)) "capacity_error" else "failed"
}
budget_table <- function() {
  b <- get("SELECT owner, current_bytes, high_water_bytes FROM duckvep_native_budget()")
  b[order(b$owner), ]
}
budget_columns <- function(prefix, b, column) {
  v <- as.list(as.numeric(b[[column]])); names(v) <- paste0(prefix, "_", b$owner)
  v
}
fmt <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) "" else
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE, digits = 15) else as.character(x)
emit_row <- function(fields) {
  line <- vapply(fields, fmt, character(1L))
  if (!file.exists(result_path)) writeLines(paste(names(line), collapse = "\t"), result_path)
  cat(paste(gsub("[\t\n]", " ", line), collapse = "\t"), "\n", file = result_path, append = TRUE, sep = "")
}

job_start <- proc.time()[["elapsed"]]
regulation_features <- NA_real_
base <- list(job = job_id, panel = panel, threads = threads, memory_limit = memory_limit,
  max_temp = max_temp, native_budget_mib = budget_mib, workers = workers,
  scratch_mib = scratch_mib, emit_mib = emit_mib)
row_of <- function(mode, outcome, reason, extra = list()) {
  c(base, list(mode = mode, outcome = outcome, reason = reason,
    rss_peak_kib = status_kib("VmHWM"), job_elapsed_s = round(proc.time()[["elapsed"]] - job_start, 3)),
    extra)
}

load_model <- function() {
  interval <- if (regulation) paste0(",\n interval_feature_query := 'SELECT regulation_feature_index, seq_region, feature_start, feature_end, feature_kind FROM duckvep_bench_model.duckvep_regulation_features ORDER BY seq_region, feature_start, regulation_feature_index'") else ""
  run(paste0("SELECT * FROM duckvep_model_load('grch38',
 'SELECT seq_region, sequence_length FROM duckvep_bench_model.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM duckvep_bench_model.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM duckvep_bench_model.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM duckvep_bench_model.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM duckvep_bench_model.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 transcript_coverage_complete := TRUE", interval, ")"))
}
set_ceilings <- function(budget_bytes) {
  run(paste0("SELECT duckvep_native_budget_set(", format(budget_bytes, scientific = FALSE), ")"))
  run(paste0("SELECT duckvep_worker_limits_set(", workers, ", ",
    format(scratch_mib * 1048576, scientific = FALSE), ", ",
    format(emit_mib * 1048576, scientific = FALSE), ", 0)"))
}
prepare_events <- function(need_complete) {
  input <- paste0("read_parquet(", q(panel), ")")
  if (!is.null(limit_rows)) input <- paste0("(SELECT * FROM ", input, " LIMIT ", as.integer(limit_rows), ")")
  run(paste0("CREATE OR REPLACE TEMP TABLE events AS SELECT ",
    "row_number() OVER (ORDER BY seq_region, position, reference, alternate)::UBIGINT AS event_index, ",
    "seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, ",
    "NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, ",
    "NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM ", input,
    " ORDER BY seq_region, position, reference, alternate"))
  run("CREATE OR REPLACE TEMP VIEW ordered_events AS SELECT * FROM events ORDER BY seq_region, position, event_index")
  if (need_complete) {
    run("CREATE OR REPLACE TEMP TABLE field_events AS WITH anchored AS (
 SELECT e.*, r.name AS chrom, length(reference) != length(alternate) AND left(reference,1)=left(alternate,1) AS strip_anchor
 FROM events e JOIN duckvep_bench_model.duckvep_sequence_regions r USING(seq_region)
 ) SELECT * EXCLUDE(strip_anchor), event_index AS record_index, 1::BIGINT AS alt_index,
 NULL::VARCHAR AS variant_id, [alternate] AS alternates, reference AS uploaded_reference,
 CASE WHEN strip_anchor THEN coalesce(nullif(substr(reference,2),''),'-') ELSE reference END AS native_reference,
 [CASE WHEN strip_anchor THEN coalesce(nullif(substr(alternate,2),''),'-') ELSE alternate END] AS native_alternates,
 chrom || ':' || (position + strip_anchor::UBIGINT)::VARCHAR || CASE WHEN position + strip_anchor::UBIGINT = position + length(reference) - 1 THEN ''
 ELSE '-' || (position + length(reference) - 1)::VARCHAR END AS native_location FROM anchored")
    run("CREATE OR REPLACE TEMP VIEW field_ordered_events AS SELECT * FROM field_events ORDER BY seq_region, position, record_index, alt_index")
    run(paste0("CREATE OR REPLACE TEMP TABLE field_metadata AS SELECT transcript_index, NULL::VARCHAR AS symbol, ",
      "canonical, NULL::VARCHAR AS tsl, NULL::VARCHAR AS appris, NULL::VARCHAR AS ccds FROM read_parquet(",
      q(canonical), ")"))
  }
  get("SELECT count(*)::DOUBLE AS n FROM events")$n
}
mode_query <- function(mode) {
  if (mode == "compact") return("SELECT * FROM query(duckvep_annotate_sql('ordered_events', 'grch38'))")
  complete <- duckvep_field_field_query(con, "native_tab17", include_identity = TRUE, model_name = "grch38")
  sub("SELECT p.record_index", "SELECT p.event_index, p.record_index", complete, fixed = TRUE)
}
run_mode <- function(mode, panel_rows, sequence) {
  final <- file.path(job_dir, paste0("output-", mode, ".parquet"))
  partial <- paste0(final, ".partial")
  unlink(c(final, partial))
  seconds <- NA_real_
  on.exit(unlink(partial), add = TRUE)
  run("SELECT duckvep_native_budget_reset_high_water()")
  outcome <- tryCatch({
    seconds <- unname(system.time(run(paste0("COPY (", mode_query(mode), ") TO ", q(partial),
      " (FORMAT PARQUET, COMPRESSION ZSTD)")))[["elapsed"]])
    NULL
  }, error = function(e) conditionMessage(e))
  b <- budget_table()
  extra <- c(list(panel_rows = panel_rows, sequence = sequence), budget_columns("native_hw", b, "high_water_bytes"),
    budget_columns("native_cur", b, "current_bytes"))
  if (!is.null(outcome)) {
    unlink(partial)
    emit_row(row_of(mode, classify(outcome), substr(gsub("[\r\n]+", " ", outcome), 1L, 400L),
      c(extra, list(published = FALSE))))
    return(invisible(FALSE))
  }
  sum_row <- get(paste0("SELECT count(*)::VARCHAR AS n, sum(hash(t)::HUGEINT)::VARCHAR AS s, ",
    "bit_xor(hash(t))::VARCHAR AS x FROM read_parquet(", q(partial), ") t"))
  bytes <- file.size(partial)
  file.rename(partial, final)
  if (!keep_output) unlink(final)
  emit_row(row_of(mode, "ok", "", c(extra, list(mode_s = round(seconds, 3),
    alleles_per_s = round(panel_rows / seconds, 1), out_rows = sum_row$n, out_bytes = bytes,
    hash_sum = sum_row$s, hash_xor = sum_row$x, published = TRUE))))
  invisible(TRUE)
}

fail_all <- function(outcome, reason, extra = list()) {
  for (mode in modes) emit_row(row_of(mode, outcome, substr(gsub("[\r\n]+", " ", reason), 1L, 400L),
    c(list(published = FALSE), extra)))
}

rc <- tryCatch({
  run(paste0("LOAD ", q(normalizePath(extension))))
  run(paste0("SET threads=", threads))
  run(paste0("SET memory_limit=", q(memory_limit)))
  run(paste0("SET temp_directory=", q(spill)))
  run(paste0("SET max_temp_directory_size=", q(max_temp)))
  run("SET preserve_insertion_order=false")
  run("SET enable_progress_bar=false")
  run("INSTALL json")
  run("LOAD json")
  run(paste0("ATTACH ", q(normalizePath(model)), " AS duckvep_bench_model (READ_ONLY)"))
  if (regulation) regulation_features <- get("SELECT count(*)::DOUBLE AS n FROM duckvep_bench_model.duckvep_regulation_features")$n
  base$regulation <- regulation; base$regulation_features <- if (regulation) regulation_features else 0
  set_ceilings(budget_mib * 1048576)
  baseline <- budget_table()
  load_s <- NA_real_
  loaded <- tryCatch({
    load_s <- unname(system.time(load_model())[["elapsed"]])
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(loaded)) {
    after <- budget_table()
    extra <- c(list(load_error = "yes", native_total_baseline = baseline$current_bytes[baseline$owner == "total"],
      native_total_after_error = after$current_bytes[after$owner == "total"]))
    if (retry) {
      # The same connection must stay usable: restore the ceiling, reload and annotate.
      set_ceilings(4096 * 1048576)
      again <- tryCatch({
        load_model()
        n <- prepare_events(FALSE)
        get("SELECT count(*)::DOUBLE AS n FROM query(duckvep_annotate_sql('ordered_events', 'grch38'))")$n > 0 && n > 0
      }, error = function(e) FALSE)
      extra$retry_ok <- isTRUE(again)
    }
    fail_all(classify(loaded), loaded, extra)
    1L
  } else {
    hw <- budget_table()
    panel_rows <- prepare_events("complete17" %in% modes)
    ok_all <- TRUE
    for (mode in modes) {
      ok <- run_mode(mode, panel_rows, sequence = which(modes == mode))
      ok_all <- ok_all && ok
    }
    if (ok_all) 0L else 1L
  }
}, error = function(e) {
  fail_all(classify(conditionMessage(e)), conditionMessage(e))
  1L
})
try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
# Load time and load-phase native high-water are reported once per job in a side file.
if (exists("load_s") && exists("hw")) {
  writeLines(c(paste(c("load_s", names(budget_columns("native_load_hw", hw, "high_water_bytes"))), collapse = "\t"),
    paste(c(fmt(round(load_s, 3)), vapply(budget_columns("native_load_hw", hw, "high_water_bytes"), fmt, "")), collapse = "\t")),
    file.path(job_dir, "load.tsv"))
}
quit(save = "no", status = rc)
