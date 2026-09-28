library(DBI)
library(Rduckvep)

# Windows keeps a database file locked until its instance shuts down, and
# closing the last connection does not shut it down while the driver object
# is still referenced. Keep the driver and shut it down explicitly before the
# same file is reopened.
open_db <- function(path, ...) {
  drv <- duckdb::duckdb(dbdir = path, ...)
  list(drv = drv, con = dbConnect(drv))
}
close_db <- function(db) {
  dbDisconnect(db$con)
  duckdb::duckdb_shutdown(db$drv)
}

local({
  path <- tempfile(fileext = ".duckdb")
  on.exit(unlink(c(path, paste0(path, ".wal"))))
  ext_config <- list(allow_unsigned_extensions = "true")

  writer <- open_db(path, shared_home = FALSE)
  dbExecute(writer$con, "CREATE TABLE caller_data AS SELECT 42 AS value")
  close_db(writer)

  reader <- open_db(path, read_only = TRUE, shared_home = FALSE,
                    allow_extensions = TRUE, config = ext_config)
  rduckvep_load(reader$con)
  expect_equal(dbGetQuery(reader$con, "SELECT value FROM caller_data")$value, 42)
  expect_true(grepl("caller_data",
                    rduckvep_annotate_sql(reader$con, "caller_data", "model"),
                    fixed = TRUE))
  close_db(reader)

  writer <- open_db(path, shared_home = FALSE, allow_extensions = TRUE,
                    config = ext_config)
  dbExecute(writer$con, "CREATE MACRO duckvep_annotate(x) AS x || ':caller'")
  rduckvep_load(writer$con)
  expect_equal(dbGetQuery(writer$con, "SELECT duckvep_annotate('kept') AS value")$value,
               "kept:caller")
  dbExecute(writer$con, "CHECKPOINT")
  close_db(writer)

  plain <- open_db(path, shared_home = FALSE)
  on.exit(close_db(plain), add = TRUE)
  expect_equal(dbGetQuery(plain$con, "SELECT duckvep_annotate('kept') AS value")$value,
               "kept:caller")
  expect_equal(dbGetQuery(plain$con, paste(
    "SELECT count(*) AS n FROM duckdb_functions()",
    "WHERE database_name = 'main' AND function_name LIKE 'duckvep_%'",
    "AND function_name <> 'duckvep_annotate'"))$n, 0)
  expect_equal(dbGetQuery(plain$con, paste(
    "SELECT count(*) AS n FROM information_schema.tables",
    "WHERE table_schema = 'main' AND table_name LIKE 'duckvep_%'"))$n, 0)
})
