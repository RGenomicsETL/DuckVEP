library(DBI)
library(duckdb)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
blocks <- strsplit(paste(readLines("test/sql/duckvep_compatible_overlap.test"), collapse = "\n"), "\n\n", fixed = TRUE)[[1L]]
setup <- blocks[startsWith(blocks, "statement ok\n")][1:2]
setup <- substring(setup, nchar("statement ok\n") + 1L)
probe <- function(extension) {
  con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste("LOAD", dbQuoteString(con, normalizePath(extension))))
  for (sql in setup) {
    if (startsWith(sql, "SELECT")) dbGetQuery(con, sql) else dbExecute(con, sql)
  }
  queries <- c(
    "SELECT event_index,seq_region,position,reference,[alternate] alternates,transcript_index,sample_index,'0|1' gt FROM overlap_calls",
    "SELECT event_index,seq_region,position,reference,[alternate] alternates,transcript_index,sample_index,'0/1' gt FROM overlap_calls")
  lapply(queries, function(calls) dbGetQuery(con, paste0("SELECT * FROM duckvep_haplotypes(",
    dbQuoteString(con, calls), ",'compatible_overlap',input_mode:='source_records',phase_policy:='vep_compat',hgvs:=true) ORDER BY ALL")))
}
old <- probe(args[1L])
new <- probe(args[2L])
stopifnot(identical(old, new))
cat("Legacy raw-source vep_compat: identical complete typed rows for phased and unphased overlapping inputs on both strands\n")
