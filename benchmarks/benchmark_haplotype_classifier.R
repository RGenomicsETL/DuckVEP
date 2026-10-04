#!/usr/bin/env Rscript
# Synthetic classifier workload for duckvep_haplotypes() on immutable extension builds (#2).
#
#   Rscript benchmarks/benchmark_haplotype_classifier.R [--copies N] [--runs R] [--threads T] EXT [EXT ...]
#
# Each EXT is a built duckvep.duckdb_extension. The script copies it to a private read-only file (content hash
# recorded), so rebuilding the working tree cannot change what is measured, and runs every extension x
# workload in a fresh R/DuckDB process R times (median reported, default 5 runs). Workloads replicate the
# committed base-R fixtures, so nothing external is needed:
#   frame     test/data/haplotype/frame_*        180 transcripts x N copies, 324 phased calls x N (slice 4 workload)
#   startstop test/data/haplotype/startstop_*    348 transcripts x N copies, 522 phased calls x N (slice 5 workload)
# Every copy has its own contigs (seq_region), so copies never interact; the model is the full transcript and exon
# tables, loaded once per process and timed separately from the query. The measured query materializes the complete
# duckvep_haplotypes() output (all columns, carriers, provenance) into a table. The checksum covers row count, status
# and reason counts and every SO set and IMPACT, so two builds can be compared for identical output or for the
# expected differences. The model addresses seq_region with 16 bits, so N x transcripts must stay below 65536
# (the default 180 copies is 32,400 and 62,640 transcripts). One thread by default (the contract gates single-core throughput).
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) return(default)
  value <- args[i + 1L]; args <<- args[-c(i, i + 1L)]; value
}
worker <- opt("worker", NA_character_)
copies <- as.integer(opt("copies", 180L))
runs <- as.integer(opt("runs", 5L))
threads <- as.integer(opt("threads", 1L))
here <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
root <- normalizePath(file.path(dirname(here), ".."))
data <- file.path(root, "test", "data", "haplotype")
workloads <- list(frame = list(prefix = "frame", transcripts = 180L),
                  startstop = list(prefix = "startstop", transcripts = 348L))

run_worker <- function(extension, workload, copies, threads) {
  suppressPackageStartupMessages({ library(DBI); library(duckdb) })
  spec <- workloads[[workload]]
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true", threads = as.character(threads))))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, sprintf("LOAD '%s'", extension))
  path <- function(suffix) file.path(data, paste0(spec$prefix, suffix))
  n <- spec$transcripts
  stopifnot(n * copies < 65536L)
  dbExecute(con, sprintf(paste("CREATE TABLE base_tx AS SELECT * FROM read_csv('%s', delim='\t', header=true)"), path("_transcripts.tsv")))
  dbExecute(con, sprintf(paste("CREATE TABLE base_exons AS SELECT * FROM read_csv('%s', delim='\t', header=true)"), path("_exons.tsv")))
  dbExecute(con, sprintf(paste("CREATE TABLE base_vcf AS SELECT string_split(line, chr(9)) f FROM",
    "(SELECT unnest(string_split(content, chr(10))) line FROM read_text('%s')) WHERE line<>'' AND NOT starts_with(line,'#')"), path(".vcf")))
  dbExecute(con, sprintf(paste("CREATE TABLE copies AS SELECT range AS copy FROM range(%d)"), copies))
  model_time <- system.time(dbGetQuery(con, sprintf(paste(
    "SELECT loaded FROM duckvep_model_load('bench',",
    "'SELECT i::UINTEGER seq_region FROM range(%d) t(i)',",
    "'SELECT (c.copy*%d + t.seq_region)::UINTEGER transcript_index,(c.copy*%d + t.seq_region)::UINTEGER seq_region,",
    "t.transcript_start::UBIGINT transcript_start,t.transcript_end::UBIGINT transcript_end,t.strand::TINYINT strand,",
    "(c.copy*%d + t.seq_region)::UINTEGER gene_index,3::UBIGINT transcript_flags,t.transcript_start::UBIGINT cds_start,",
    "t.transcript_end::UBIGINT cds_end,t.cds_sequence::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''''::BLOB pre_cds_sequence,''''::BLOB post_cds_sequence FROM base_tx t, copies c ORDER BY 1',",
    "'SELECT (c.copy*%d + e.seq_region)::UINTEGER transcript_index,e.exon_start::UBIGINT exon_start,e.exon_end::UBIGINT exon_end,",
    "e.exon_cdna_start::UBIGINT exon_cdna_start,e.exon_cdna_end::UBIGINT exon_cdna_end,e.phase::TINYINT phase,",
    "e.end_phase::TINYINT end_phase FROM base_exons e, copies c ORDER BY 1,e.exon_cdna_start')"),
    n * copies, n, n, n, n)))[["elapsed"]]
  dbExecute(con, sprintf(paste(
    "CREATE TABLE calls AS SELECT (c.copy*1000000 + row_number() OVER (PARTITION BY c.copy ORDER BY v.rowid))::BIGINT AS event_index,",
    "(c.copy*%d + t.seq_region)::INT AS seq_region,v.f[2]::BIGINT AS position,v.f[4] AS reference,v.f[5] AS alternate,1 AS alt_index,",
    "(c.copy*%d + t.seq_region)::INT AS transcript_index,0 AS sample_index,",
    "(CASE v.f[10] WHEN '1|0' THEN [1,0] WHEN '0|1' THEN [0,1] ELSE [1,1] END)::INTEGER[] AS alleles,",
    "[false,true]::BOOLEAN[] AS phase_before,NULL::BIGINT AS phase_set",
    "FROM (SELECT *,rowid FROM base_vcf) v JOIN base_tx t ON t.\"case\"=v.f[1], copies c"), n, n))
  input_records <- dbGetQuery(con, "SELECT count(*) n FROM calls")$n
  query_time <- system.time(dbExecute(con, paste(
    "CREATE TABLE result AS SELECT * FROM duckvep_haplotypes('SELECT * FROM calls ORDER BY seq_region, position, event_index', 'bench')")))[["elapsed"]]
  summary <- dbGetQuery(con, paste(
    "SELECT count(*) AS n_rows, sum(len(carriers)) AS n_carriers, count(*) FILTER (WHERE prediction_status='predicted') AS n_predicted,",
    "count(*) FILTER (WHERE prediction_status='eligible_classifier_pending') AS n_pending,",
    "sum(hash(coalesce(array_to_string(list_sort(haplotype_consequences), ','), '') || '/' || coalesce(haplotype_impact, '') ||",
    "'/' || prediction_status || '/' || prediction_reason)::HUGEINT) AS checksum FROM result"))
  cat(sprintf("RESULT\t%s\t%d\t%d\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%s\n", workload, copies, input_records, model_time, query_time,
    summary$n_rows, summary$n_carriers, summary$n_predicted, summary$n_pending, format(summary$checksum, scientific = FALSE)))
}

if (!is.na(worker)) {
  run_worker(args[1L], worker, copies, threads)
  quit(save = "no")
}
options(width = 220)
stopifnot(length(args) >= 1L, all(file.exists(args)))
private <- tempfile("duckvep-bench-")
dir.create(private)
copies_of <- vapply(args, function(path) {
  # DuckDB derives the entry point from the file name, so each copy keeps it in its own directory.
  dir.create(file.path(private, sprintf("ext%02d", match(path, args))))
  target <- file.path(private, sprintf("ext%02d", match(path, args)), "duckvep.duckdb_extension")
  stopifnot(file.copy(path, target))
  Sys.chmod(target, "0444")
  target
}, "")
hashes <- vapply(copies_of, function(path) substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L), "")
cat(sprintf("copies=%d runs=%d threads=%d\n", copies, runs, threads))
results <- list()
for (i in seq_along(copies_of)) for (workload in names(workloads)) {
  rows <- lapply(seq_len(runs), function(r) {
    out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(here), "--worker", workload,
      "--copies", copies, "--threads", threads, shQuote(copies_of[[i]])), stdout = TRUE, stderr = TRUE)
    line <- grep("^RESULT\t", out, value = TRUE)
    if (length(line) != 1L) stop("worker failed:\n", paste(out, collapse = "\n"))
    strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
  })
  m <- do.call(rbind, rows)
  stopifnot(length(unique(m[, 10L])) == 1L)
  results[[length(results) + 1L]] <- data.frame(extension = args[[i]], sha256 = substr(hashes[[i]], 1L, 12L),
    workload = workload, input_records = as.integer(m[1L, 3L]), model_load_s = median(as.numeric(m[, 4L])),
    query_s = median(as.numeric(m[, 5L])), min_query_s = min(as.numeric(m[, 5L])),
    records_per_s = round(as.integer(m[1L, 3L]) / median(as.numeric(m[, 5L]))),
    output_rows = as.integer(m[1L, 6L]), carriers = as.integer(m[1L, 7L]), predicted = as.integer(m[1L, 8L]),
    pending = as.integer(m[1L, 9L]), checksum = m[1L, 10L])
  print(results[[length(results)]], row.names = FALSE)
}
unlink(private, recursive = TRUE, force = TRUE)
