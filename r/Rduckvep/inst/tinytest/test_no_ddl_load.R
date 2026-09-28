library(DBI)
library(Rduckvep)

# Windows keeps the database file locked until the driver is finalized, so
# collect it before the same file is reopened.
close_db <- function(con) {
  dbDisconnect(con, shutdown = TRUE)
  invisible(gc())
}

local({
  path <- tempfile(fileext = ".duckdb")
  on.exit(unlink(c(path, paste0(path, ".wal"))))
  writer <- dbConnect(duckdb::duckdb(dbdir = path, shared_home = FALSE))
  dbExecute(writer, "CREATE TABLE caller_data AS SELECT 42 AS value")
  close_db(writer)

  reader <- dbConnect(duckdb::duckdb(dbdir = path, read_only = TRUE,
                                     shared_home = FALSE, allow_extensions = TRUE,
                                     config = list(allow_unsigned_extensions = "true")))
  rduckvep_load(reader)
  expect_equal(dbGetQuery(reader, "SELECT value FROM caller_data")$value, 42)
  expect_true(grepl("caller_data", rduckvep_annotate_sql(reader, "caller_data", "model"),
                    fixed = TRUE))
  close_db(reader)

  writer <- dbConnect(duckdb::duckdb(dbdir = path, shared_home = FALSE,
                                     allow_extensions = TRUE,
                                     config = list(allow_unsigned_extensions = "true")))
  dbExecute(writer, "CREATE MACRO duckvep_annotate(x) AS x || ':caller'")
  rduckvep_load(writer)
  expect_equal(dbGetQuery(writer, "SELECT duckvep_annotate('kept') AS value")$value,
               "kept:caller")
  dbExecute(writer, "CHECKPOINT")
  close_db(writer)

  plain <- dbConnect(duckdb::duckdb(dbdir = path, shared_home = FALSE))
  on.exit(close_db(plain), add = TRUE)
  expect_equal(dbGetQuery(plain, "SELECT duckvep_annotate('kept') AS value")$value,
               "kept:caller")
  expect_equal(dbGetQuery(plain, paste(
    "SELECT count(*) AS n FROM duckdb_functions()",
    "WHERE database_name = 'main' AND function_name LIKE 'duckvep_%'",
    "AND function_name <> 'duckvep_annotate'"))$n, 0)
  expect_equal(dbGetQuery(plain, paste(
    "SELECT count(*) AS n FROM information_schema.tables",
    "WHERE table_schema = 'main' AND table_name LIKE 'duckvep_%'"))$n, 0)
})
