args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
extension <- normalizePath(args[[1]], mustWork = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
fixture <- file.path(dirname(script), "../data/projection/parallel.sql")
local({
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE,
    config = list(allow_unsigned_extensions = "true", threads = "1")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
  DBI::dbExecute(con, paste(readLines(fixture), collapse = "\n"))
  sql <- DBI::dbGetQuery(con, paste(
    "SELECT duckvep_transcript_projection_sql('projection_events',",
    "'projection_annotations', 'projection_transcripts') AS sql"))$sql[[1]]
  for (threads in c(1L, 2L, 4L)) {
    DBI::dbExecute(con, paste("SET threads =", threads))
    if (threads > 1L) DBI::dbExecute(con, "ANALYZE")
    stopifnot(identical(DBI::dbGetQuery(con,
      "SELECT current_setting('disabled_optimizers') AS value")$value, ""))
    DBI::dbExecute(con, paste("CREATE TABLE result", threads, " AS ", sql, sep = ""))
    stopifnot(DBI::dbGetQuery(con, paste0("SELECT count(*) AS n FROM result", threads))$n == 692444)
  }
  DBI::dbExecute(con, "SET disabled_optimizers = 'build_side_probe_side'")
  DBI::dbExecute(con, paste("CREATE TABLE reference AS", sql))
  for (threads in c(1L, 2L, 4L)) {
    difference <- DBI::dbGetQuery(con, sprintf(paste(
      "WITH differences AS ((FROM result%d EXCEPT ALL FROM reference)",
      "UNION ALL (FROM reference EXCEPT ALL FROM result%d))",
      "SELECT count(*) AS n FROM differences"), threads, threads))$n
    stopifnot(difference == 0)
  }
  cat("Projection passed: 692444 rows at 1, 2 and 4 threads; complete multisets identical.\n")
})
