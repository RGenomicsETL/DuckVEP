library(tinytest)

nominal <- rduckvep_prepare_sv_geometry(13546123, "G", "<INS>",
  "SVTYPE=INS;END=13546123")
confidence <- rduckvep_prepare_sv_geometry(13546123, "G", "<INS>",
  "SVTYPE=INS;END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG")
literal <- rduckvep_prepare_sv_geometry(13546123, "G", "GATG", ".")
expect_identical(confidence$status, "ok")
expect_identical(confidence$mode, "structural")
expect_equal(confidence$nominal_start, nominal$nominal_start)
expect_equal(confidence$nominal_end, nominal$nominal_end)
expect_equal(confidence$outer_start, 13546122)
expect_equal(confidence$outer_end, 13546125)
expect_true(is.na(confidence$inserted_sequence))
expect_identical(confidence$source_sequence, "ATG")
expect_identical(literal$inserted_sequence, "ATG")
expect_identical(literal$mode, "literal_insertion")
expect_identical(rduckvep_prepare_sv_geometry(20, "A", "<DEL>",
  "END=22;CIPOS=3,-3")$status, "unsupported_geometry")
expect_identical(rduckvep_prepare_sv_geometry(20, "AC", "<DEL>",
  "END=22")$status, "unsupported_geometry")
expect_identical(rduckvep_prepare_sv_geometry(20, "A", "<DEL>",
  "END=20.5")$status, "unsupported_geometry")
expect_identical(rduckvep_prepare_sv_geometry(20, "A", "<DEL>",
  "END=19")$status, "unsupported_geometry")
expect_identical(rduckvep_prepare_sv_geometry(2147483647, "A", "AT",
  ".")$status, "unsupported_geometry")
