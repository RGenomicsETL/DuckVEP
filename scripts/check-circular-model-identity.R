#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !file.exists(args[1L]) || !file.exists(args[2L])) {
  stop("usage: check-circular-model-identity.R MODEL_COPY EXTENSION_COPY")
}

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

con <- dbConnect(
  duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true"))
)
on.exit(dbDisconnect(con, shutdown = TRUE))
quote_string <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, paste("LOAD", quote_string(normalizePath(args[2L]))))
dbExecute(con, paste(
  "ATTACH", quote_string(normalizePath(args[1L])), "AS stored (READ_ONLY)"
))

old <- dbGetQuery(con, "SELECT * FROM stored.model_receipt")
if (nrow(old) != 1L) {
  stop("stored model must contain one receipt")
}
dbExecute(con, paste(
  "CREATE TEMP VIEW regions_with_topology AS SELECT *,",
  "false AS circular FROM stored.model_regions"
))
dbExecute(con, "CREATE TEMP VIEW transcripts AS FROM stored.model_transcripts")
regulation <- dbGetQuery(con, paste(
  "SELECT count(*) > 0 AS present FROM information_schema.tables",
  "WHERE table_catalog = 'stored' AND table_name = 'duckvep_regulation_features'"
))$present
if (regulation) {
  dbExecute(con, paste(
    "CREATE TEMP VIEW regulation_features AS",
    "FROM stored.duckvep_regulation_features"
  ))
}
parameters <- c(
  "regions_with_topology", "transcripts", old$source_name,
  old$source_version, old$assembly, old$source_manifest_sha256,
  old$reference_sha256, old$transcript_filter
)
sql <- paste0(
  "SELECT model_sha256 FROM query(duckvep_model_receipt_sql(",
  paste(vapply(parameters, quote_string, character(1L)), collapse = ", "),
  ifelse(regulation, ", {regulation_features_table: 'regulation_features'}", ""),
  "))"
)
observed <- dbGetQuery(con, sql)$model_sha256
if (!identical(observed, old$model_sha256)) {
  stop(sprintf("stored model fingerprint changed: %s != %s", observed, old$model_sha256))
}
message("Stored model fingerprint unchanged: ", observed)
