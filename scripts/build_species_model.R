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
core_schema <- 'core_snapshot.source_core'
filter <- "is_current=1; stable_id present; biotype!=artifact; !readthrough_tra; primary-assembly FASTA regions"
if (length(args) == 8L) {
  dbExecute(con, paste(readLines(args[[8L]], warn = FALSE), collapse = '\n'))
  core_schema <- 'core_input'
  filter <- paste0(filter, '; quarantine miRNA 35-56 on ENSMUST00000175144 (53 nt cDNA)')
}
dbExecute(con, paste0('CREATE TABLE reference_chunks AS FROM read_parquet(',
  q(normalizePath(args[[2L]])), ')'))
dbExecute(con, paste0("CREATE TABLE model_regions AS FROM duckvep_ensembl_regions(",
  q(core_schema), ", 'reference_chunks', ", q(args[[5L]]), ")"))
dbExecute(con, paste0("CREATE TABLE model_transcripts AS FROM duckvep_ensembl_transcripts(",
  q(core_schema), ", 'reference_chunks', ", q(args[[5L]]), ")"))
dbExecute(con, paste0("CREATE TABLE model_receipt AS FROM duckvep_model_receipt(",
  "'model_regions', 'model_transcripts', 'ensembl_core', ", q(args[[6L]]), ", ",
  q(args[[5L]]), ", ", q(source_sha), ", ", q(reference_sha), ", ",
  q(filter), ")"))
dbExecute(con, 'DROP TABLE reference_chunks')
if (length(args) == 8L) dbExecute(con, 'DROP SCHEMA core_input CASCADE')
dbExecute(con, 'CHECKPOINT')
print(dbGetQuery(con, 'SELECT * FROM model_receipt'))
print(dbGetQuery(con, 'SELECT codon_table, count(*) n FROM model_transcripts GROUP BY ALL ORDER BY 1'))
