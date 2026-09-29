library(DBI)
library(duckdb)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:3) stop("usage: Rscript benchmarks/scale_model_load.R EXTENSION MODEL [ordered|unordered]")
mode <- if (length(args) == 3L) args[[3L]] else "ordered"
stopifnot(mode %in% c("ordered", "unordered"))
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
q <- function(table, keys) {
  sql <- paste("SELECT * FROM source_model.", table, sep = "")
  if (mode == "ordered") sql <- paste(sql, "ORDER BY", keys)
  quoted(sql)
}
regions <- quoted(paste("SELECT seq_region, sequence_length FROM source_model.duckvep_sequence_regions",
  if (mode == "ordered") "ORDER BY seq_region" else ""))
loaded <- dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('grch38', ", regions, ", ",
  q("duckvep_transcripts", "seq_region, transcript_start, transcript_index"), ", ",
  q("duckvep_exons", "transcript_index, exon_cdna_start"), ", ",
  "mature_mirna_query := ", q("duckvep_mature_mirna", "transcript_index, mature_mirna_start"), ", ",
  "peptide_edit_query := ", q("duckvep_peptide_edits", "transcript_index, protein_position"), ", ",
  "transcript_coverage_complete := TRUE)"))
elapsed <- proc.time()[["elapsed"]] - start
stopifnot(identical(loaded$loaded, TRUE))
status <- readLines("/proc/self/status", warn = FALSE)
peak <- grep("^VmHWM:", status, value = TRUE)
cat(sprintf("load_seconds\t%.3f\npeak_rss_kib\t%s\n", elapsed,
  sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", peak)))
dbDisconnect(con, shutdown = TRUE)
