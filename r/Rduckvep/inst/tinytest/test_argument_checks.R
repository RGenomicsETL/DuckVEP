library(tinytest)
library(DBI)
# Argument checks run before any DuckVEP SQL, so the extension is not loaded here.
con <- dbConnect(duckdb::duckdb(shared_home = FALSE))

# A string argument that is not one nonempty string is refused by name.
expect_error(rduckvep_haplotypes(con, NA_character_, "m"), "calls_query must be one nonempty string")
expect_error(rduckvep_haplotypes(con, "SELECT 1", ""), "model_name must be one nonempty string")
expect_error(rduckvep_coding_transcripts(con, c("a", "b"), "m"), "events_query must be one nonempty string")
expect_error(rduckvep_coding_calls(con, "m", 1L), "path must be one nonempty string")
expect_error(rduckvep_load_model(con, "m", "SELECT 1", NULL, "SELECT 1"),
             "transcripts_query must be one nonempty string")
expect_error(rduckvep_save_model(con, "m", ""), "path must be one nonempty string")
expect_error(rduckvep_restore_model(con, NA_character_, "p"), "name must be one nonempty string")

# The output table arguments are checked first.
expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", overwrite = NA), "overwrite must be TRUE or FALSE")
expect_error(rduckvep_coding_calls(con, "m", "p", overwrite = "yes"), "overwrite must be TRUE or FALSE")
expect_error(rduckvep_coding_transcripts(con, "SELECT 1", "m", table_name = ""),
             "table_name must be one nonempty string")
dbExecute(con, "CREATE TABLE present AS SELECT 1 AS one")
expect_error(rduckvep_coding_calls(con, "m", "p", table_name = "present"), "Table already exists: present")

# Haplotype options and workspace capacities.
expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", hgvs = NA), "hgvs must be TRUE or FALSE")
expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", input_mode = "source_records"),
             "source_records requires phase_policy='vep_compat'")
expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", "strict", 5), "capacities must have unique SQL parameter names")
expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", max_ploidy = 2, max_ploidy = 3),
             "capacities must have unique SQL parameter names")
expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", `max ploidy` = 2),
             "capacities must have unique SQL parameter names")
for (bad in list(0, -1, 1.5, Inf, NA_real_, 2^54, "2", c(1, 2))) {
  expect_error(rduckvep_haplotypes(con, "SELECT 1", "m", max_ploidy = bad),
               "max_ploidy must be one positive, exactly representable integer")
}
expect_identical(Rduckvep:::.duckvep_capacity_literals(list(max_ploidy = 2, workspace_limit = 2^40)),
                 list(max_ploidy = "2", workspace_limit = "1099511627776"))
expect_identical(Rduckvep:::.duckvep_capacity_literals(list()), list())

dbDisconnect(con, shutdown = TRUE)
