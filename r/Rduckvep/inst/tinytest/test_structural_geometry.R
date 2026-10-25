library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
# read_bcf-shaped input: a list ALT and typed INFO columns; END as text exercises
# the whole-number rule.
dbExecute(con, "CREATE TEMP TABLE sv_input AS SELECT row_number() OVER () AS event_index, * FROM (VALUES
  (13546123, 'G', ['<INS>'], '13546123', NULL::INTEGER[], NULL::INTEGER[], NULL::INTEGER[], NULL),
  (13546123, 'G', ['<INS>'], '13546123', NULL, [-2, 0], [0, 2], 'ATG'),
  (13546123, 'G', ['GATG'], NULL, NULL, NULL, NULL, NULL),
  (20, 'A', ['<DEL>'], '22', NULL, [3, -3], NULL, NULL),
  (20, 'AC', ['<DEL>'], '22', NULL, NULL, NULL, NULL),
  (20, 'A', ['<DEL>'], '20.5', NULL, NULL, NULL, NULL),
  (20, 'A', ['<DEL>'], '19', NULL, NULL, NULL, NULL),
  (2147483647, 'A', ['AT'], NULL, NULL, NULL, NULL, NULL),
  (20, 'A', ['<DEL>'], '40', [-5], NULL, NULL, NULL)
) t(\"POS\", \"REF\", \"ALT\", \"INFO_END\", \"INFO_SVLEN\", \"INFO_CIPOS\", \"INFO_CIEND\", \"INFO_SEQ\")")
rows <- rduckvep_prepare_sv_geometry(con, "sv_input")
expect_identical(rows$status[1:3], rep("ok", 3))
expect_identical(rows$status[4:8], rep("unsupported_geometry", 5))
expect_identical(rows$mode[2], "structural")
expect_equal(rows$nominal_start[1], rows$nominal_start[2])
expect_equal(rows$nominal_end[1], rows$nominal_end[2])
expect_equal(rows$outer_start[2], 13546122)
expect_equal(rows$outer_end[2], 13546125)
expect_true(is.na(rows$inserted_sequence[2]))
expect_identical(rows$source_sequence[2], "ATG")
expect_identical(rows$inserted_sequence[3], "ATG")
expect_identical(rows$mode[3], "literal_insertion")
# SVLEN takes precedence over END, as VEP reads them.
expect_identical(rows$status[9], "ok")
expect_equal(c(rows$nominal_start[9], rows$nominal_end[9]), c(21, 25))
dbDisconnect(con, shutdown = TRUE)
