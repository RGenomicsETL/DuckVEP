if (requireNamespace("duckdb", quietly = TRUE)) local({
  con <- rduckvep_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  expect_true(is.data.frame(DBI::dbGetQuery(con, "SELECT duckvep_phase_call([0,1], [true,false]) AS phase")))
  expect_identical(DBI::dbGetQuery(con,
    "SELECT duckvep_repeat_alleles([{unit:'A',count:1}], [{unit:'C',count:2}], TRUE).alternate AS allele")$allele,
    "CC")
}
