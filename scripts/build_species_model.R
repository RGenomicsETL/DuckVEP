#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(digest)
})
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(7L, 8L)) {
  stop("usage: build_species_model.R CORE_DB REFERENCE_PARQUET FASTA MODEL_DB ASSEMBLY SOURCE_VERSION EXTENSION [CORE_INPUT_SQL]", call. = FALSE)
}
stopifnot(all(file.exists(args[c(1L, 2L, 3L, 7L)])), !file.exists(args[[4L]]))
if (length(args) == 8L) stopifnot(file.exists(args[[8L]]))
source_sha <- digest(file = args[[1L]], algo = "sha256")
reference_sha <- digest(file = args[[3L]], algo = "sha256")
con <- dbConnect(duckdb(shared_home = FALSE,
                        config = list(allow_unsigned_extensions = "true")),
                 dbdir = args[[4L]])
on.exit(dbDisconnect(con, shutdown = TRUE))
q <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, 'SET threads=2')
dbExecute(con, paste0('LOAD ', q(normalizePath(args[[7L]]))))
dbExecute(con, paste0('ATTACH ', q(normalizePath(args[[1L]])),
  ' AS core_snapshot (READ_ONLY)'))
filter <- paste0("is_current=1; stable_id present; biotype!=artifact; !readthrough_tra; ",
  if (args[[5L]] %in% c("BDGP6.54", "TAIR10")) "toplevel FASTA regions" else "primary-assembly FASTA regions")
if (length(args) == 8L) {
  dbExecute(con, paste(readLines(args[[8L]], warn = FALSE), collapse = '\n'))
  if (args[[5L]] == 'GRCm39') {
    filter <- paste0(filter, '; quarantine miRNA 35-56 on ENSMUST00000175144 (53 nt cDNA)')
  }
} else {
  dbExecute(con, 'CREATE SCHEMA core_input')
  for (name in c('attrib_type', 'coord_system', 'seq_region', 'seq_region_attrib',
                 'gene', 'transcript', 'transcript_attrib', 'translation',
                 'translation_attrib', 'exon', 'exon_transcript')) {
    dbExecute(con, paste0('CREATE VIEW core_input.', name,
                          ' AS FROM core_snapshot.source_core.', name))
  }
}
dbExecute(con, paste0('CREATE TABLE reference_chunks AS FROM read_parquet(',
  q(normalizePath(args[[2L]])), ')'))
dbExecute(con, paste0("CREATE TABLE model_regions AS FROM query(duckvep_ensembl_regions_sql(",
  q('core_input'), ", 'reference_chunks', ", q(args[[5L]]), "))"))
dbExecute(con, paste0("CREATE TABLE model_transcripts AS FROM query(duckvep_ensembl_transcripts_sql(",
  q('core_input'), ", 'reference_chunks', ", q(args[[5L]]), "))"))
dbExecute(con, paste0("CREATE TABLE model_receipt AS FROM query(duckvep_model_receipt_sql(",
  "'model_regions', 'model_transcripts', 'ensembl_core', ", q(args[[6L]]), ", ",
  q(args[[5L]]), ", ", q(source_sha), ", ", q(reference_sha), ", ",
  q(filter), "))"))
dbExecute(con, 'DROP TABLE reference_chunks')
dbExecute(con, 'DROP SCHEMA core_input CASCADE')
dbExecute(con, 'CHECKPOINT')
print(dbGetQuery(con, 'SELECT * FROM model_receipt'))
print(dbGetQuery(con, 'SELECT codon_table, count(*) n FROM model_transcripts GROUP BY ALL ORDER BY 1'))
