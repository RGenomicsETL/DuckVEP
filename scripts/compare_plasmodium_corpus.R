#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: compare_plasmodium_corpus.R MODEL_DB VCF PROVENANCE_TSV VEP_JSON OUTPUT_DIR", call. = FALSE)
}
model <- normalizePath(args[[1L]], mustWork = TRUE)
vcf <- normalizePath(args[[2L]], mustWork = TRUE)
provenance <- normalizePath(args[[3L]], mustWork = TRUE)
oracle <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
con <- dbConnect(duckdb(shared_home = FALSE,
                        config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
q <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, paste0("LOAD ", q(normalizePath(
  "build/release/extension/duckvep/duckvep.duckdb_extension", mustWork = TRUE))))
dbExecute(con, paste0("ATTACH ", q(model), " AS model (READ_ONLY)"))
dbExecute(con, "SET threads=4")
dbExecute(con, "CREATE TABLE duckvep_sequence_regions AS
  SELECT seq_region, sequence_length, seq_region_name AS name FROM model.model_regions")
dbExecute(con, "CREATE TABLE duckvep_transcripts AS
  SELECT transcript_index, seq_region, transcript_start, transcript_end,
         strand, gene_index, transcript_flags, cds_start, cds_end,
         cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence
  FROM model.model_transcripts")
dbExecute(con, "CREATE TABLE duckvep_exons AS
  SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start,
         e.exon_cdna_end, e.phase, e.end_phase
  FROM model.model_transcripts, unnest(exons) AS x(e)")
dbExecute(con, "CREATE TABLE duckvep_transcript_names AS
  SELECT transcript_index, transcript_stable_id AS transcript_id,
         transcript_biotype AS biotype, codon_table, strand, seq_region_name
  FROM model.model_transcripts")
loaded <- dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load(
  'plasmodium',
  ", q("SELECT seq_region, sequence_length FROM duckvep_sequence_regions ORDER BY seq_region"), ",
  ", q("SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index"), ",
  ", q("SELECT transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase FROM duckvep_exons ORDER BY transcript_index, exon_cdna_start"), ",
  transcript_coverage_complete := true)"))$loaded
stopifnot(identical(loaded, TRUE))

v <- read.delim(vcf, comment.char = "#", header = FALSE,
                col.names = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"),
                stringsAsFactors = FALSE)
p <- read.delim(provenance, stringsAsFactors = FALSE)
stopifnot(nrow(v) == nrow(p), setequal(v$id, p$id), !anyDuplicated(v$id))
dbWriteTable(con, "variants", v, temporary = TRUE)
dbWriteTable(con, "provenance", p, temporary = TRUE)
dbExecute(con, "CREATE TEMP VIEW events AS SELECT v.*, r.seq_region
  FROM variants v JOIN duckvep_sequence_regions r ON v.chrom = r.name")
stopifnot(dbGetQuery(con, "SELECT count(*) n FROM events")$n == nrow(v))
dbExecute(con, "CREATE TEMP TABLE duck_pairs AS
  WITH annotation AS (
    SELECT v.id, unnest(_duckvep_annotate_small_rich(
      'plasmodium', v.seq_region, v.pos::UBIGINT, v.ref, v.alt, 5000::UBIGINT)) AS a
    FROM (SELECT * FROM events ORDER BY seq_region, pos, id) v
  ), terms AS (
    SELECT a.id, n.transcript_id AS tx, unnest(string_split(a.a.consequence, '&')) term,
           a.a.status AS status, a.a.reason AS reason
    FROM annotation a JOIN duckvep_transcript_names n
      ON n.transcript_index = a.a.transcript_index
  )
  SELECT id, tx, list_aggregate(list_sort(list_distinct(list(term))), 'string_agg', '&')
    AS consequence, string_agg(DISTINCT status, '&') AS status,
    string_agg(DISTINCT coalesce(reason, ''), '&') AS reason
  FROM terms GROUP BY id, tx")
records <- lapply(readLines(oracle, warn = FALSE),
                  jsonlite::fromJSON, simplifyVector = FALSE)
stopifnot(length(records) == nrow(v), setequal(v$id,
          vapply(records, `[[`, character(1L), "id")))
oracle_rows <- lapply(records, function(record) {
  lapply(record$transcript_consequences, function(tc) {
    data.frame(id = record$id, tx = tc$transcript_id,
               consequence = paste(sort(unique(unlist(tc$consequence_terms))),
                                   collapse = "&"))
  })
})
oracle_rows <- do.call(rbind, unlist(oracle_rows, recursive = FALSE))
dbWriteTable(con, "oracle_rows", oracle_rows, temporary = TRUE)
dbExecute(con, "CREATE TEMP TABLE oracle_pairs AS
  WITH terms AS (SELECT id, tx, unnest(string_split(consequence, '&')) term
                FROM oracle_rows)
  SELECT id, tx, list_aggregate(list_sort(list_distinct(list(term))),
                                'string_agg', '&') AS consequence
  FROM terms GROUP BY id, tx")
dbExecute(con, "CREATE TEMP TABLE differential AS
  SELECT coalesce(o.id, d.id) AS id, coalesce(o.tx, d.tx) AS tx,
         o.consequence AS oracle_terms, d.consequence AS duckvep_terms,
         d.status, d.reason,
         CASE WHEN o.id IS NULL THEN 'duckvep_only'
              WHEN d.id IS NULL THEN 'oracle_only'
              WHEN o.consequence != d.consequence THEN 'terms_differ'
              ELSE 'exact' END AS verdict
  FROM oracle_pairs o FULL OUTER JOIN duck_pairs d USING (id, tx)")
dbExecute(con, "CREATE TEMP VIEW annotated_diff AS
  SELECT x.*, v.chrom, v.pos, v.ref, v.alt,
         CASE WHEN length(v.ref) = length(v.alt) THEN
           CASE WHEN length(v.ref) = 1 THEN 'SNV' ELSE 'MNV' END
           WHEN length(v.ref) < length(v.alt) THEN 'insertion' ELSE 'deletion' END AS allele_shape,
         t.codon_table, t.biotype, t.strand,
         CASE WHEN v.chrom LIKE '%MIT%' THEN 'mitochondrion'
              WHEN v.chrom LIKE '%API%' THEN 'apicoplast' ELSE 'nuclear' END AS contig_class
  FROM differential x JOIN variants v ON x.id = v.id
  LEFT JOIN duckvep_transcript_names t ON x.tx = t.transcript_id")
counts <- dbGetQuery(con, "SELECT verdict, count(*) n FROM differential GROUP BY verdict ORDER BY verdict")
write.csv(counts, file.path(output, "verdicts.csv"), row.names = FALSE)
write.csv(dbGetQuery(con, "SELECT status, reason, count(*) n FROM duck_pairs
                        GROUP BY ALL ORDER BY 1, 2"),
          file.path(output, "duckvep_status.csv"), row.names = FALSE)
for (key in c("codon_table", "contig_class", "biotype", "strand", "allele_shape")) {
  stmt <- paste0("SELECT ", key, ", verdict, count(*) n FROM annotated_diff GROUP BY ALL ORDER BY 1, 2")
  write.csv(dbGetQuery(con, stmt), file.path(output, paste0(key, ".csv")), row.names = FALSE)
}
dbExecute(con, "CREATE TEMP VIEW term_diff AS
  SELECT *, unnest(string_split(coalesce(oracle_terms, ''), '&')) AS term,
         'oracle' AS side FROM annotated_diff
  UNION ALL
  SELECT *, unnest(string_split(coalesce(duckvep_terms, ''), '&')) AS term,
         'duckvep' AS side FROM annotated_diff")
write.csv(dbGetQuery(con, "SELECT side, term, verdict, count(*) n FROM term_diff
                        WHERE term <> '' GROUP BY ALL ORDER BY 1, 2, 3"),
          file.path(output, "so_term.csv"), row.names = FALSE)
write.csv(dbGetQuery(con, "SELECT * FROM annotated_diff WHERE verdict <> 'exact'
                        ORDER BY id, tx"), file.path(output, "disagreements.csv"),
          row.names = FALSE, na = "")
write.csv(dbGetQuery(con, "SELECT p.*, x.tx AS annotated_tx, x.oracle_terms,
                         x.duckvep_terms, x.verdict
                         FROM provenance p JOIN differential x USING (id)
                         WHERE p.source = 'codon_witness' ORDER BY p.table, p.id, x.tx"),
          file.path(output, "codon_witness_pairs.csv"), row.names = FALSE)
witness <- dbGetQuery(con, "SELECT p.table, x.oracle_terms, x.duckvep_terms,
                              x.verdict FROM provenance p
                        JOIN differential x ON x.id = p.id AND x.tx = p.tx
                        WHERE p.source = 'codon_witness'")
stopifnot(nrow(witness) == sum(p$source == 'codon_witness'),
          all(witness$verdict == 'exact'),
          all(witness$oracle_terms[witness$table == 1L] == 'stop_gained'),
          all(witness$oracle_terms[witness$table == 4L] == 'synonymous_variant'),
          all(witness$oracle_terms[witness$table == 11L] == 'start_lost'),
          all(c(1L, 4L, 11L) %in% witness$table))
print(counts)
cat("variants:", nrow(v), "oracle pairs:",
    dbGetQuery(con, "SELECT count(*) n FROM oracle_pairs")$n,
    "DuckVEP pairs:", dbGetQuery(con, "SELECT count(*) n FROM duck_pairs")$n, "\n")
