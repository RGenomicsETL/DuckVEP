#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: check_r_v2_haplotypes.R <v1-extension> <v2-extension>")
library(DBI)
library(Rduckvep)

check_host <- function(extension, v2) {
  con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
    config = list(allow_unsigned_extensions = "true")))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste("LOAD", dbQuoteString(con, normalizePath(extension))))
  stopifnot(rduckvep_load_model(con, "haplotypes",
    "SELECT 0::UINTEGER seq_region",
    "SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region,
      100::UBIGINT transcript_start, 111::UBIGINT transcript_end, 1::TINYINT strand,
      0::UINTEGER gene_index, 3::UBIGINT transcript_flags,
      100::UBIGINT cds_start, 111::UBIGINT cds_end, 'ATGGCTGCTTAA'::BLOB cds_sequence,
      1::UTINYINT codon_table, ''::BLOB pre_cds_sequence, ''::BLOB post_cds_sequence",
    "SELECT 0::UINTEGER transcript_index, 100::UBIGINT exon_start, 111::UBIGINT exon_end,
      1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase"))
  create <- if (v2) "CREATE TEMP TABLE" else "CREATE TABLE"
  dbExecute(con, paste(create, "calls AS SELECT 1::UBIGINT event_index, 0::UINTEGER seq_region,
    104::UBIGINT AS position, 'C' AS reference, 'A' AS alternate,
    1::UINTEGER alt_index, 0::UINTEGER transcript_index, 0::UBIGINT sample_index,
    [0,1]::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before, NULL::BIGINT phase_set"))
  query <- "SELECT * FROM calls"
  rows <- rduckvep_haplotypes(con, query, "haplotypes", hgvs = TRUE)
  stopifnot(nrow(rows) > 0L)
  rduckvep_haplotypes(con, query, "haplotypes", hgvs = TRUE, table_name = 'output"table')
  stopifnot(identical(dbReadTable(con, 'output"table'), rows))
  budget <- "SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner"
  before <- dbGetQuery(con, budget)
  failure <- tryCatch(rduckvep_haplotypes(con, "SELECT * FROM missing_calls", "haplotypes",
    table_name = 'output"table', overwrite = TRUE), error = conditionMessage)
  stopifnot(is.character(failure), grepl("missing_calls", failure),
    identical(dbReadTable(con, 'output"table'), rows),
    identical(dbGetQuery(con, budget), before))
  failure <- tryCatch(getFromNamespace(".duckvep_create_table", "Rduckvep")(
    con, 'output"table', "SELECT error('output execution failed')", TRUE), error = conditionMessage)
  stopifnot(is.character(failure), grepl("output execution failed", failure),
    identical(dbReadTable(con, 'output"table'), rows))
  failure <- tryCatch(rduckvep_haplotypes(con, query, "haplotypes", table_name = 'output"table'),
    error = conditionMessage)
  stopifnot(is.character(failure), grepl("already exists", failure))
  rduckvep_haplotypes(con, query, "haplotypes", hgvs = TRUE,
    table_name = 'output"table', overwrite = TRUE)
  stopifnot(identical(dbReadTable(con, 'output"table'), rows),
    identical(dbGetQuery(con, budget), before),
    dbGetQuery(con, "SELECT duckvep_model_drop('haplotypes')")[[1L]])
  rows
}

v1 <- check_host(args[[1L]], FALSE)
v2 <- check_host(args[[2L]], TRUE)
stopifnot(identical(v1, v2))
cat("R haplotypes: equal v1/v2 rows, caller-local v2 input, atomic replacement, and cleanup passed\n")
