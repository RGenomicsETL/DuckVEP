#!/usr/bin/env Rscript
# One pinned-core measurement of GRCh38 model load and a 1M-allele annotation:
# wall seconds, process peak RSS and, where the extension has the native budget,
# the native high-water mark by owner.
#
#   taskset -c 2 Rscript benchmarks/scale_budget.R EXTENSION [gnomad1m|bench]
#
# EXTENSION must be an immutable copy named duckvep.duckdb_extension.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(TRUE)
if (!length(args) %in% 1:2) stop("usage: Rscript benchmarks/scale_budget.R EXTENSION [gnomad1m|bench]")
extension <- normalizePath(args[[1L]], mustWork = TRUE)
corpus <- if (length(args) == 2L) args[[2L]] else "gnomad1m"
stopifnot(corpus %in% c("gnomad1m", "bench"))
model <- Sys.getenv("DUCKVEP_SCALE_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
panel <- "/root/duckvep/data/gnomad-v4.1/panels/genomes-v1-77d2cdff65171780/panel-1000000.parquet"
metadata <- "benchmarks/data/scale_contracts/canonical-metadata.parquet"
source("benchmarks/duckvep_field_projection.R", local = TRUE)
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
q <- function(x) as.character(dbQuoteString(con, x))
run <- function(sql) invisible(dbExecute(con, sql))
get <- function(sql) dbGetQuery(con, sql)
elapsed <- function(code) unname(system.time(force(code))[["elapsed"]])
status <- function(key) {
  line <- grep(paste0("^", key, ":"), readLines("/proc/self/status", warn = FALSE), value = TRUE)
  as.numeric(sub(paste0("^", key, ":[[:space:]]*([0-9]+).*"), "\\1", line))
}
emit <- function(key, value) cat(sprintf("%s\t%s\n", key, format(value, scientific = FALSE, trim = TRUE)))
run(paste0("LOAD ", q(extension)))
run("SET threads=1")
run("SET memory_limit='8GB'")
run(paste0("SET temp_directory=", q(file.path(tempdir(), "scale-budget-spill"))))
run("SET preserve_insertion_order=false")
run(paste0("ATTACH ", q(model), " AS duckvep_bench_model (READ_ONLY)"))
input <- if (corpus == "bench") "duckvep_bench_model.bench_variants" else paste0("read_parquet(", q(panel), ")")
run(paste0("CREATE TEMP TABLE events AS SELECT row_number() OVER (ORDER BY seq_region, position, reference, alternate)::UBIGINT AS event_index, ",
  "seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, ",
  "NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, ",
  "NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, ",
  "NULL::UBIGINT AS mate_position FROM ", input, " ORDER BY seq_region, position, reference, alternate"))
run("CREATE TEMP VIEW ordered_events AS SELECT * FROM events ORDER BY seq_region, position, event_index")
run("CREATE TEMP TABLE field_events AS WITH anchored AS (
 SELECT e.*, r.name AS chrom, length(reference) != length(alternate) AND
 left(reference,1)=left(alternate,1) AS strip_anchor
 FROM events e JOIN duckvep_bench_model.duckvep_sequence_regions r USING(seq_region)
 ) SELECT * EXCLUDE(strip_anchor), event_index AS record_index, 1::BIGINT AS alt_index,
 NULL::VARCHAR AS variant_id, [alternate] AS alternates, reference AS uploaded_reference,
 CASE WHEN strip_anchor THEN coalesce(nullif(substr(reference,2),''),'-') ELSE reference END AS native_reference,
 [CASE WHEN strip_anchor THEN coalesce(nullif(substr(alternate,2),''),'-') ELSE alternate END] AS native_alternates,
 chrom || ':' || (position + strip_anchor::UBIGINT)::VARCHAR || CASE WHEN
 position + strip_anchor::UBIGINT = position + length(reference) - 1 THEN '' ELSE
 '-' || (position + length(reference) - 1)::VARCHAR END AS native_location FROM anchored")
run("CREATE TEMP VIEW field_ordered_events AS SELECT * FROM field_events ORDER BY seq_region, position, record_index, alt_index")
run(paste0("CREATE TEMP TABLE field_metadata AS SELECT transcript_index, NULL::VARCHAR AS symbol, ",
  "canonical, NULL::VARCHAR AS tsl, NULL::VARCHAR AS appris, NULL::VARCHAR AS ccds FROM read_parquet(", q(metadata), ")"))
budget <- function(label) {
  b <- tryCatch(get("SELECT owner, current_bytes, high_water_bytes FROM duckvep_native_budget()"), error = function(e) NULL)
  if (is.null(b)) return(invisible())
  for (i in seq_len(nrow(b))) {
    emit(paste0(label, "_current_mib_", b$owner[i]), round(b$current_bytes[i] / 1048576, 1))
    emit(paste0(label, "_highwater_mib_", b$owner[i]), round(b$high_water_bytes[i] / 1048576, 1))
  }
}
emit("rss_before_load_kib", status("VmRSS"))
load_s <- elapsed(run("SELECT * FROM duckvep_model_load('grch38',
 'SELECT seq_region, sequence_length FROM duckvep_bench_model.duckvep_sequence_regions',
 'SELECT * FROM duckvep_bench_model.duckvep_transcripts',
 'SELECT * FROM duckvep_bench_model.duckvep_exons',
 mature_mirna_query := 'SELECT * FROM duckvep_bench_model.duckvep_mature_mirna',
 peptide_edit_query := 'SELECT * FROM duckvep_bench_model.duckvep_peptide_edits',
 transcript_coverage_complete := TRUE)"))
emit("load_s", round(load_s, 3))
emit("peak_rss_after_load_kib", status("VmHWM"))
invisible(budget("load"))
tryCatch(run("SELECT duckvep_native_budget_reset_high_water()"), error = function(e) NULL)
n_in <- get("SELECT count(*) AS n FROM events")$n
emit("input_rows", n_in)
compact_s <- elapsed(n <- get("SELECT count(*) AS n FROM query(duckvep_annotate_sql('ordered_events', 'grch38'))")$n)
emit("compact_rows", n); emit("compact_s", round(compact_s, 3))
complete_s <- elapsed(n <- get("SELECT count(*) AS n FROM query(duckvep_annotate_projected_sql('field_ordered_events', 'grch38'))")$n)
emit("complete17_rows", n); emit("complete17_s", round(complete_s, 3))
emit("peak_rss_after_annotation_kib", status("VmHWM"))
invisible(budget("annotate"))
