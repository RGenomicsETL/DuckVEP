library(DBI)
library(duckdb)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript benchmarks/scale_model_load.R EXTENSION MODEL")
extension <- normalizePath(args[[1L]], mustWork = TRUE)
model <- normalizePath(args[[2L]], mustWork = TRUE)
con <- dbConnect(duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
quoted <- function(value) as.character(dbQuoteString(con, value))
run <- function(sql) invisible(dbExecute(con, sql))
run(paste0("LOAD ", quoted(extension)))
run("SET threads=1")
run("SET memory_limit='8GB'")
run(paste0("ATTACH ", quoted(model), " AS source_model (READ_ONLY)"))
start <- proc.time()[["elapsed"]]
loaded <- dbGetQuery(con, "SELECT loaded FROM duckvep_model_load('grch38',
 'SELECT seq_region, sequence_length FROM source_model.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM source_model.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM source_model.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM source_model.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM source_model.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 transcript_coverage_complete := TRUE)")
elapsed <- proc.time()[["elapsed"]] - start
stopifnot(identical(loaded$loaded, TRUE))
status <- readLines("/proc/self/status", warn = FALSE)
peak <- grep("^VmHWM:", status, value = TRUE)
cat(sprintf("load_seconds\t%.3f\npeak_rss_kib\t%s\n", elapsed,
  sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", peak)))
dbDisconnect(con, shutdown = TRUE)
