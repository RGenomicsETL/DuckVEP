library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
dbExecute(con, "CREATE TEMP TABLE sv_input AS SELECT row_number() OVER () AS event_index, * FROM (VALUES
  (13546123, 'G', '<INS>', 'SVTYPE=INS;END=13546123'),
  (13546123, 'G', '<INS>', 'SVTYPE=INS;END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG'),
  (13546123, 'G', 'GATG', '.'),
  (20, 'A', '<DEL>', 'END=22;CIPOS=3,-3'),
  (20, 'AC', '<DEL>', 'END=22'),
  (20, 'A', '<DEL>', 'END=20.5'),
  (20, 'A', '<DEL>', 'END=19'),
  (2147483647, 'A', 'AT', '.')
) t(pos, ref, alt, info)")
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
dbDisconnect(con, shutdown = TRUE)
