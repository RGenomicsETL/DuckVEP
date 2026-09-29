library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)

# Every wrapper reports a missing relation as an R error and leaves the connection usable.
missing <- list(
  function() rduckvep_prepare_sv_geometry(con, "no_such_table"),
  function() rduckvep_prepare_breakend_pairs(con, "no_such_table"),
  function() rduckvep_prepare_breakend_fusion(con, "no_such_table", "also_missing"),
  function() rduckvep_prepare_structural_hgvs(con, "no_such_table", "also_missing"),
  function() rduckvep_prepare_expansionhunter(con, "no_such_table", "also_missing"))
for (call in missing) {
  expect_error(call(), "no_such_table")
  expect_identical(DBI::dbGetQuery(con, "SELECT 1 AS one")$one, 1L)
}
# NA and vector relation names are refused before any SQL is built.
expect_error(rduckvep_prepare_sv_geometry(con, NA_character_))
expect_error(rduckvep_prepare_sv_geometry(con, c("a", "b")))
expect_error(rduckvep_prepare_sv_geometry(con, "x", bogus = 1), "unknown option")

# NULL rows and extreme values become explicit non-ok rows, not statement errors.
dbExecute(con, "CREATE TEMP TABLE hostile (event_index BIGINT, pos BIGINT, ref VARCHAR,
  alt VARCHAR, info VARCHAR)")
dbExecute(con, "INSERT INTO hostile VALUES (0, NULL, NULL, NULL, NULL),
  (1, 9223372036854775807, 'N', '<DEL>', 'END=9223372036854775807'),
  (2, -9223372036854775808, '', '<DEL>', 'END=5'),
  (3, 100, 'N', '<DEL>', 'END=200;CIPOS=-9223372036854775808,9223372036854775807'),
  (4, 100, 'N', '<DEL>', 'END=200;CIPOS=-2,3')")
rows <- rduckvep_prepare_sv_geometry(con, "hostile")
expect_identical(rows$status, c(rep("unsupported_geometry", 4L), "ok"))
expect_identical(DBI::dbGetQuery(con, "SELECT 1 AS one")$one, 1L)
dbDisconnect(con, shutdown = TRUE)
