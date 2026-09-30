#!/usr/bin/env Rscript

# End-to-end DuckVEP worker for the FastVEP comparison. The outer benchmark
# process pins this worker and records GNU time. This worker deliberately
# includes extension/model loading, sequential VCF decoding, explicit sorting,
# annotation, text projection, and COPY to a real local file.

suppressMessages({
  library(DBI)
  library(duckdb)
  library(glue)
  library(optparse)
})

root <- tryCatch(
  system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE),
  error = function(e) "."
)

op <- OptionParser()
op <- add_option(
  op,
  "--extension",
  default = file.path(root, "build", "release", "duckhts.duckdb_extension")
)
op <- add_option(
  op,
  "--duckhts-extension",
  dest = "duckhts_extension",
  default = "",
  help = "DuckHTS extension, loaded in a separate in-process instance, that decodes the VCF and GFF3 when --extension is the standalone duckvep extension"
)
op <- add_option(op, "--model", default = "")
op <- add_option(op, "--input", default = "")
op <- add_option(op, "--output", default = "")
op <- add_option(op, "--output-contract", dest = "output_contract", default = "operational17")
op <- add_option(op, "--fasta", default = "")
op <- add_option(op, "--gff3", default = "")
op <- add_option(op, "--include-identity", dest = "include_identity", action = "store_true", default = FALSE)
op <- add_option(op, "--threads", type = "integer", default = 1L)
op <- add_option(op, "--distance", type = "integer", default = 5000L)
op <- add_option(op, "--memory-limit", dest = "memory_limit", default = "4GB")
op <- add_option(op, "--max-spill", dest = "max_spill", default = "8GB")
op <- add_option(
  op,
  "--profile-json",
  dest = "profile_json",
  default = "",
  help = "optional DuckDB JSON profile path for the measured COPY query"
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
root <- normalizePath(root[[1L]], mustWork = TRUE)
source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"), local = TRUE)
invisible(duckvep_fastvep_fields(opt$output_contract))
if (opt$output_contract == "operational17" &&
    (nzchar(opt$fasta) || nzchar(opt$gff3) || opt$include_identity)) {
  die("--fasta, --gff3 and --include-identity require a field-complete output contract")
}
if (opt$output_contract == "vep_csq" && (!nzchar(opt$fasta) || !nzchar(opt$gff3))) {
  die("vep_csq requires --fasta and --gff3 for HGVS and source-labelled symbols")
}

required <- c(opt$extension, opt$model, opt$input,
  if (nzchar(opt$duckhts_extension)) opt$duckhts_extension,
  c(opt$fasta, opt$gff3)[nzchar(c(opt$fasta, opt$gff3))])
missing <- required[!nzchar(required) | !file.exists(required)]
if (length(missing) != 0L) {
  die("missing input(s):\n{paste(missing, collapse = '\n')}")
}
if (!nzchar(opt$output)) {
  die("--output is required")
}
if (opt$threads < 1L || opt$threads > 1024L) {
  die("--threads must be from 1 through 1024")
}
if (opt$distance < 0L || opt$distance > 2^32 - 1) {
  die("--distance must fit an unsigned 32-bit integer")
}

extension <- normalizePath(opt$extension)
model <- normalizePath(opt$model)
input <- normalizePath(opt$input)
output <- normalizePath(opt$output, mustWork = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

drv <- duckdb(
  dbdir = ":memory:",
  config = list(allow_unsigned_extensions = "true")
)
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql_q <- function(x) as.character(dbQuoteString(con, x))

invisible(dbExecute(con, glue("LOAD {sql_q(extension)}")))
# DuckHTS embeds an older DuckVEP whose function names collide with the
# standalone extension, so it never shares a database instance with it. When
# --duckhts-extension is given, a second in-process instance decodes the VCF and
# the GFF3 into Parquet (inside this timed process, on the same pinned CPUs) and
# read_bcf/read_gff are bound to those files for the field builders.
if (nzchar(opt$duckhts_extension)) {
  stage_dir <- file.path(tempdir(), "duckhts-stage")
  dir.create(stage_dir, recursive = TRUE)
  staged_vcf <- file.path(stage_dir, "source.parquet")
  staged_gff <- file.path(stage_dir, "gff.parquet")
  reader_drv <- duckdb(dbdir = ":memory:", config = list(allow_unsigned_extensions = "true"))
  reader <- dbConnect(reader_drv)
  reader_q <- function(x) as.character(dbQuoteString(reader, x))
  invisible(dbExecute(reader, glue("LOAD {reader_q(normalizePath(opt$duckhts_extension))}")))
  invisible(dbExecute(reader, glue("PRAGMA threads={opt$threads}")))
  invisible(dbExecute(reader, glue("SET memory_limit = {reader_q(opt$memory_limit)}")))
  invisible(dbExecute(reader, "SET preserve_insertion_order = true"))
  # Records stay in file order (the unindexed VCF has no contig-parallel scan);
  # with more than one thread, htslib decompression uses worker threads.
  decode_threads <- if (opt$threads > 1L) opt$threads else 0L
  invisible(dbExecute(reader, glue("COPY (SELECT CHROM, POS, ID, REF, ALT
    FROM read_bcf({reader_q(input)}, scan_mode := 'sequential',
      decompression_threads := {decode_threads}))
    TO {reader_q(staged_vcf)} (FORMAT PARQUET)")))
  if (nzchar(opt$gff3)) {
    # The field builders read only these attribute keys from gene and transcript
    # lines. read_gff's attributes_map builds a map for all 10.7M lines before
    # any filter applies, so read the raw attribute string of the 0.7M wanted
    # lines and split out just these keys (URL-decoded values; checked equal to
    # attributes_map on every wanted line).
    gff_keys <- "'ID', 'Name', 'tag', 'transcript_support_level', 'ccdsid'"
    invisible(dbExecute(reader, glue("COPY (SELECT feature, map_from_entries(list_transform(
        list_filter(string_split(attributes, ';'), x -> split_part(x, '=', 1) IN ({gff_keys})),
        x -> {{'k': split_part(x, '=', 1), 'v': url_decode(substr(x, strpos(x, '=') + 1))}}))
        AS attributes_map
      FROM read_gff({reader_q(normalizePath(opt$gff3))}, scan_mode := 'sequential')
      WHERE feature NOT IN ('exon', 'CDS', 'chromosome', 'biological_region',
        'five_prime_UTR', 'three_prime_UTR')) TO {reader_q(staged_gff)} (FORMAT PARQUET)")))
  }
  dbDisconnect(reader, shutdown = TRUE)
  invisible(dbExecute(con, glue("CREATE TEMP MACRO read_bcf(path, scan_mode := 'sequential',
    decompression_threads := 0) AS TABLE SELECT * FROM read_parquet({sql_q(staged_vcf)})")))
  if (nzchar(opt$gff3)) {
    invisible(dbExecute(con, glue("CREATE TEMP MACRO read_gff(path, attributes_map := TRUE,
      scan_mode := 'sequential') AS TABLE SELECT * FROM read_parquet({sql_q(staged_gff)})")))
  }
}
invisible(dbExecute(con, glue("PRAGMA threads={opt$threads}")))
invisible(dbExecute(con, glue("SET memory_limit = {sql_q(opt$memory_limit)}")))
invisible(dbExecute(con, glue("SET max_temp_directory_size = {sql_q(opt$max_spill)}")))
invisible(dbExecute(con, glue("SET temp_directory = {sql_q(file.path(tempdir(), 'duckdb'))}")))
invisible(dbExecute(con, "SET preserve_insertion_order = false"))
if (nzchar(opt$profile_json)) {
  profile_json <- normalizePath(opt$profile_json, mustWork = FALSE)
  dir.create(dirname(profile_json), recursive = TRUE, showWarnings = FALSE)
  invisible(dbExecute(con, "PRAGMA enable_profiling = 'json'"))
  invisible(dbExecute(
    con,
    glue("PRAGMA profiling_output = {sql_q(profile_json)}")
  ))
}
invisible(dbExecute(
  con,
  glue("ATTACH {sql_q(model)} AS duckvep_bench_model (READ_ONLY)")
))
relations_source <- file.path(root, "r/duckhtsbench/R/duckvep_relations.R")
if (file.exists(relations_source)) {
  source(relations_source, local = TRUE)
  model_relations <- duckhts_bench_duckvep_relations(con, "duckvep_bench_model")
} else {
  # Standalone repository: use the flat duckvep_* relations when the model
  # carries them, otherwise project the nested model_regions/model_transcripts
  # exactly as the DuckHTS benchmark helper does.
  flat <- dbGetQuery(con, "SELECT table_name FROM information_schema.tables
    WHERE table_catalog = 'duckvep_bench_model' AND table_schema = 'main'
      AND starts_with(table_name, 'duckvep_')")$table_name
  relation_names <- c("duckvep_sequence_regions", "duckvep_transcripts",
    "duckvep_exons", "duckvep_mature_mirna", "duckvep_peptide_edits")
  if (all(relation_names %in% flat)) {
    model_relations <- as.list(setNames(
      paste0("duckvep_bench_model.main.", relation_names), relation_names))
  } else {
    regions <- "duckvep_bench_model.main.model_regions"
    transcripts <- "duckvep_bench_model.main.model_transcripts"
    projections <- c(
      duckvep_sequence_regions = paste("SELECT seq_region, sequence_length, seq_region_name AS name FROM", regions),
      duckvep_transcripts = paste("SELECT * FROM", transcripts),
      duckvep_exons = paste("SELECT transcript_index, exon.* FROM", transcripts,
        "CROSS JOIN UNNEST(exons) AS u(exon)"),
      duckvep_mature_mirna = paste("SELECT transcript_index, region.* FROM", transcripts,
        "CROSS JOIN UNNEST(mature_mirna_regions) AS u(region)"),
      duckvep_peptide_edits = paste("SELECT transcript_index, edit.* FROM", transcripts,
        "CROSS JOIN UNNEST(peptide_edits) AS u(edit)"))
    model_relations <- as.list(setNames(
      paste0("(", projections, ") AS duckvep_prepared"), names(projections)))
  }
}
invisible(dbExecute(
  con,
  glue("CREATE TEMP TABLE duckvep_bench_regions AS
   SELECT seq_region, name
   FROM {model_relations[['duckvep_sequence_regions']]}")
))
invisible(dbExecute(
  con,
  "CREATE TEMP TABLE duckvep_bench_transcript_labels AS
   SELECT transcript_index, gene_stable_id, transcript_stable_id, strand
   FROM duckvep_bench_model.model_transcripts"
))

region_query <- paste(
  if (nzchar(opt$fasta)) "SELECT seq_region, sequence_length, name AS seq_region_name"
  else "SELECT seq_region, sequence_length",
  "FROM", model_relations[["duckvep_sequence_regions"]],
  "ORDER BY seq_region"
)
transcript_query <- paste(
  "SELECT transcript_index, seq_region, transcript_start, transcript_end,",
  "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
  "codon_table, pre_cds_sequence, post_cds_sequence",
  "FROM", model_relations[["duckvep_transcripts"]],
  "ORDER BY seq_region, transcript_start, transcript_index"
)
exon_query <- paste(
  "SELECT transcript_index, exon_start, exon_end, exon_cdna_start,",
  "exon_cdna_end, phase, end_phase",
  "FROM", model_relations[["duckvep_exons"]],
  "ORDER BY transcript_index, exon_cdna_start"
)
mature_mirna_query <- paste(
  "SELECT transcript_index, mature_mirna_start, mature_mirna_end",
  "FROM", model_relations[["duckvep_mature_mirna"]],
  "ORDER BY transcript_index, mature_mirna_start"
)
peptide_edit_query <- paste(
  "SELECT transcript_index, protein_position, alternate_amino_acid",
  "FROM", model_relations[["duckvep_peptide_edits"]],
  "ORDER BY transcript_index, protein_position"
)

loaded <- dbGetQuery(
  con,
  glue(
    "SELECT loaded FROM duckvep_model_load(
       'fastvep_comparison',
       {sql_q(region_query)}, {sql_q(transcript_query)}, {sql_q(exon_query)},
       mature_mirna_query := {sql_q(mature_mirna_query)},
       peptide_edit_query := {sql_q(peptide_edit_query)},
       {if (nzchar(opt$fasta)) paste0('reference_fasta := ', sql_q(normalizePath(opt$fasta)), ',') else ''}
       transcript_coverage_complete := TRUE)"
  )
)$loaded
if (length(loaded) != 1L || !isTRUE(loaded[[1L]])) {
  die("DuckVEP model load failed")
}
if (opt$output_contract == "operational17") {
  invisible(dbExecute(con, "DETACH duckvep_bench_model"))
} else {
  duckvep_fastvep_prepare_fields(con, input, opt$output_contract, opt$distance, opt$gff3)
  query <- duckvep_fastvep_field_query(
    con, opt$output_contract, opt$include_identity, opt$distance
  )
  invisible(dbExecute(con, glue("COPY ({query}) TO {sql_q(output)}
    (FORMAT CSV, DELIMITER E'\\t', HEADER TRUE, QUOTE '', ESCAPE '')")))
  quit(status = 0L)
}

# FastVEP's native tab writer has a fixed 17-column contract. DuckVEP emits the
# same columns here so real-file costs are comparable. Fields not owned by the
# current rich consequence row remain '-' rather than being fabricated. This
# projection intentionally does not load regulatory/motif intervals because
# FastVEP's transcript cache does not contain them.
copy_query <- glue(
  "COPY (
     WITH source AS (
       SELECT
         row_number() OVER ()::UBIGINT AS record_index,
         CHROM AS chrom, POS::UBIGINT AS position, ID AS variant_id,
         REF AS reference, ALT AS alternates
       FROM read_bcf(
         {sql_q(input)}, scan_mode := 'sequential', decompression_threads := 0
       )
     ),
     alleles AS (
       SELECT
         s.record_index, a.alt_index, s.chrom, s.position, s.variant_id,
         s.reference, a.alternate
       FROM source s
       CROSS JOIN UNNEST(s.alternates) WITH ORDINALITY AS a(alternate, alt_index)
       WHERE regexp_full_match(s.reference, '[ACGTNacgtn]+')
         AND regexp_full_match(a.alternate, '[ACGTNacgtn]+')
         AND upper(s.reference) <> upper(a.alternate)
     ),
     prepared AS (
       SELECT
         a.record_index, a.alt_index, r.seq_region, a.chrom, a.position,
         a.variant_id, upper(a.reference) AS reference,
         upper(a.alternate) AS alternate
       FROM alleles a
       JOIN duckvep_bench_regions r
         ON r.name = regexp_replace(a.chrom, '^chr', '')
     ),
     annotated AS (
       SELECT
         v.*,
         unnest(_duckvep_annotate_small_rich(
           'fastvep_comparison', v.seq_region, v.position,
           v.reference, v.alternate, {opt$distance}
         )) AS annotation
       FROM (
         SELECT *
         FROM prepared
         ORDER BY seq_region, position, record_index, alt_index
       ) v
     )
     SELECT
       coalesce(
         a.variant_id,
         concat(a.chrom, ':', a.position, ':', a.reference, ':', a.alternate)
       ) AS Uploaded_variation,
       concat(a.chrom, ':', a.position) AS Location,
       a.alternate AS Allele,
       coalesce(t.gene_stable_id, '-') AS Gene,
       coalesce(t.transcript_stable_id, '-') AS Feature,
       CASE WHEN a.annotation.transcript_index IS NULL THEN '-' ELSE 'Transcript' END
         AS Feature_type,
       replace(a.annotation.consequence, '&', ',') AS Consequence,
       coalesce(a.annotation.cdna_position::VARCHAR, '') AS cDNA_position,
       coalesce(a.annotation.cds_position::VARCHAR, '') AS CDS_position,
       coalesce(a.annotation.protein_position::VARCHAR, '') AS Protein_position,
       CASE
         WHEN a.annotation.reference_amino_acid IS NULL
           AND a.annotation.alternate_amino_acid IS NULL THEN '-'
         ELSE concat(
           coalesce(a.annotation.reference_amino_acid, ''), '/',
           coalesce(a.annotation.alternate_amino_acid, '')
         )
       END AS Amino_acids,
       '-' AS Codons,
       '-' AS Existing_variation,
       a.annotation.impact AS IMPACT,
       '-' AS DISTANCE,
       coalesce(t.strand::VARCHAR, '-') AS STRAND,
       '-' AS FLAGS
     FROM annotated a
     LEFT JOIN duckvep_bench_transcript_labels t
       ON t.transcript_index = a.annotation.transcript_index
   ) TO {sql_q(output)} (
     FORMAT CSV, DELIMITER E'\\t', HEADER TRUE, QUOTE '', ESCAPE ''
   )"
)

invisible(dbExecute(con, copy_query))
