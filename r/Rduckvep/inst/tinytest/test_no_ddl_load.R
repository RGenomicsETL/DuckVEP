library(DBI)
library(Rduckvep)

# Each reopen of the database file runs in its own R process. Windows keeps a
# DuckDB file locked by the process that opened it, and process exit is the
# only release that is reliable on every platform. It is also the scenario
# these gates protect: another process opening the file read-only.
rscript <- file.path(R.home("bin"),
                     if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_phase <- function(path, code) {
  script <- tempfile(fileext = ".R")
  result <- tempfile(fileext = ".rds")
  on.exit(unlink(c(script, result)))
  writeLines(c(
    "library(DBI)",
    "library(Rduckvep)",
    sprintf("path <- %s", deparse(normalizePath(path, mustWork = FALSE))),
    "ext_config <- list(allow_unsigned_extensions = 'true')",
    code,
    sprintf("saveRDS(result, %s)", deparse(result))
  ), script)
  # system2(env =) is not honoured on Windows; the child inherits R_LIBS.
  old_libs <- Sys.getenv("R_LIBS", unset = NA)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  on.exit(if (is.na(old_libs)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old_libs),
          add = TRUE)
  output <- system2(rscript, c("--vanilla", shQuote(script)), stdout = TRUE,
                    stderr = TRUE)
  if (!file.exists(result)) stop(paste(output, collapse = "\n"), call. = FALSE)
  readRDS(result)
}

local({
  path <- tempfile(fileext = ".duckdb")
  on.exit(unlink(c(path, paste0(path, ".wal"))))

  run_phase(path, c(
    "con <- dbConnect(duckdb::duckdb(dbdir = path))",
    "dbExecute(con, 'CREATE TABLE caller_data AS SELECT 42 AS value')",
    "dbDisconnect(con, shutdown = TRUE)",
    "result <- TRUE"
  ))

  read_only <- run_phase(path, c(
    "con <- dbConnect(duckdb::duckdb(dbdir = path, read_only = TRUE,",
    "  allow_extensions = TRUE, config = ext_config))",
    "rduckvep_load(con)",
    "result <- list(",
    "  value = dbGetQuery(con, 'SELECT value FROM caller_data')$value,",
    "  sql = rduckvep_annotate_sql(con, 'caller_data', 'model'))",
    "dbDisconnect(con, shutdown = TRUE)"
  ))
  expect_equal(read_only$value, 42)
  expect_true(grepl("caller_data", read_only$sql, fixed = TRUE))

  kept <- run_phase(path, c(
    "con <- dbConnect(duckdb::duckdb(dbdir = path, allow_extensions = TRUE,",
    "  config = ext_config))",
    "dbExecute(con, \"CREATE MACRO duckvep_annotate(x) AS x || ':caller'\")",
    "rduckvep_load(con)",
    "result <- dbGetQuery(con, \"SELECT duckvep_annotate('kept') AS value\")$value",
    "dbExecute(con, 'CHECKPOINT')",
    "dbDisconnect(con, shutdown = TRUE)"
  ))
  expect_equal(kept, "kept:caller")

  plain <- dbConnect(duckdb::duckdb(dbdir = path))
  on.exit(dbDisconnect(plain, shutdown = TRUE), add = TRUE)
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
