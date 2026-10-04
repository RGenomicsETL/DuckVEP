#!/usr/bin/env Rscript
# One-core five-replicate wall times for the fixed scale workloads.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(TRUE)
if (!length(args) %in% 2:3 || !args[[2L]] %in% c("bench", "gnomad1m") ||
    (length(args) == 3L && args[[3L]] != "--fingerprint-only")) {
  stop("usage: Rscript benchmarks/scale_fusion_timing.R <source checkout> <bench|gnomad1m> [--fingerprint-only]")
}
source_root <- normalizePath(args[[1L]])
corpus <- args[[2L]]
model <- "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb"
metadata <- file.path(source_root, "benchmarks/data/scale_contracts/canonical-metadata.parquet")
panel <- "/root/duckvep/data/gnomad-v4.1/panels/genomes-v1-77d2cdff65171780/panel-1000000.parquet"
extension <- file.path(source_root, "build/release/extension/duckvep/duckvep.duckdb_extension")
source(file.path(source_root, "benchmarks/duckvep_field_projection.R"), local = TRUE)
version <- if (identical(source_root, normalizePath(getwd()))) "after" else "main"
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
q <- function(x) as.character(dbQuoteString(con, x))
run <- function(sql) invisible(dbExecute(con, sql))
get <- function(sql) dbGetQuery(con, sql)
elapsed <- function(code) unname(system.time(force(code))[["elapsed"]])
run(paste0("LOAD ", q(extension)))
run("SET threads=1")
run("SET memory_limit='8GB'")
run("SET preserve_insertion_order=false")
run(paste0("SET temp_directory=", q(file.path(tempdir(), "scale-timing-spill"))))
run(paste0("ATTACH ", q(model), " AS duckvep_bench_model (READ_ONLY)"))
input <- if (corpus == "bench") "duckvep_bench_model.bench_variants" else
  paste0("read_parquet(", q(panel), ")")
run(paste0("CREATE TEMP TABLE events AS SELECT row_number() OVER (ORDER BY seq_region, position, reference, alternate)::UBIGINT AS event_index, ",
  "seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, ",
  "NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, ",
  "NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, ",
  "NULL::UBIGINT AS mate_position FROM ", input,
  " ORDER BY seq_region, position, reference, alternate"))
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
load_sql <- "SELECT * FROM duckvep_model_load('grch38',
 'SELECT seq_region, sequence_length FROM duckvep_bench_model.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM duckvep_bench_model.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM duckvep_bench_model.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM duckvep_bench_model.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM duckvep_bench_model.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 transcript_coverage_complete := TRUE)"
compact <- "SELECT event_index, transcript_index, gene_index, consequence_mask, region_mask,
 impact_code, status_code, reason_code, cdna_position, cds_position, protein_position,
 reference_amino_acid_code, alternate_amino_acid_code, nmd_prediction_code,
 nmd_escape_reasons, regulation_feature_index, overlap_object_code
 FROM query(duckvep_annotate_sql('ordered_events', 'grch38'))"
compact_count <- "SELECT count(*) AS n FROM query(duckvep_annotate_sql('ordered_events', 'grch38'))"
complete <- if (version == "main")
  duckvep_field_field_query(con, "native_tab17", include_identity = TRUE) else
  duckvep_field_field_query(con, "native_tab17", include_identity = TRUE, model_name = "grch38")
if (version == "main") {
  complete <- sub("'field_comparison'", "'grch38'", complete, fixed = TRUE)
  complete_count <- "SELECT count(*) AS n FROM (
   SELECT unnest(_duckvep_annotate_small_projected('grch38', seq_region, position,
     reference, alternate, 5000, 5000)) FROM field_ordered_events)"
} else {
  complete_count <- "SELECT count(*) AS n FROM query(duckvep_annotate_projected_sql(
    'field_ordered_events', 'grch38'))"
}
fingerprint_only <- length(args) == 3L
if (fingerprint_only) run(load_sql)
results <- list()
for (i in seq_len(if (fingerprint_only) 0L else 5L)) {
  if (i > 1L) run("SELECT duckvep_model_drop('grch38')")
  load_time <- elapsed(run(load_sql))
  for (contract in c("compact", "complete17")) {
    annotation_sql <- if (contract == "compact") compact_count else complete_count
    output_sql <- if (contract == "compact") compact else complete
    annotation_time <- elapsed(n <- get(annotation_sql)$n)
    path <- tempfile(paste0("scale-", contract, "-"), fileext = ".parquet")
    total_time <- elapsed(run(paste0("COPY (", output_sql, ") TO ", q(path),
      " (FORMAT PARQUET, COMPRESSION UNCOMPRESSED)")))
    bytes <- file.info(path)$size
    unlink(path)
    results[[length(results) + 1L]] <- data.frame(version = version, corpus = corpus,
      contract = contract, repetition = i, input_rows = get("SELECT count(*) AS n FROM events")$n,
      output_rows = n, load_s = load_time, annotation_s = annotation_time,
      formatting_sink_s = total_time - annotation_time, total_s = total_time,
      sink_bytes = bytes)
  }
  print(results[[length(results)]], row.names = FALSE)
}
out <- do.call(rbind, results)
for (contract in c("compact", "complete17")) {
  query <- if (contract == "compact") compact else complete
  fingerprint <- get(paste0("SELECT count(*) AS rows, sum(hash(a)::HUGEINT)::VARCHAR AS hash_sum, ",
    "bit_xor(hash(a))::VARCHAR AS hash_xor FROM (", query, ") a"))
  message(version, " ", corpus, " ", contract, " fingerprint: ",
    paste(unlist(fingerprint), collapse = " "))
}
if (!fingerprint_only) write.table(out, file = paste0("/tmp/scale-fusion-", version, "-", corpus, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
