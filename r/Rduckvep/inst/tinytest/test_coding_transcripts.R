local({
  con <- rduckvep_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  model <- paste0("SELECT loaded FROM duckvep_model_load('tinytest_discovery', ",
    "'SELECT 0::UINTEGER seq_region', ",
    "'SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,111::UBIGINT cds_end,''ATGGCTGCTTAA''::BLOB cds_sequence,1::UTINYINT codon_table,''''::BLOB pre_cds_sequence,''''::BLOB post_cds_sequence', ",
    "'SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase')")
  expect_true(DBI::dbGetQuery(con, model)$loaded)
  events <- paste("SELECT * FROM (VALUES",
    "(1::UBIGINT, 0::UINTEGER, 104::UBIGINT, 'C', 'A'),",   # substitution in the CDS
    "(2::UBIGINT, 0::UINTEGER, 99::UBIGINT, 'TCC', 'T'),",  # deletes CDS bases 100 and 101
    "(3::UBIGINT, 0::UINTEGER, 111::UBIGINT, 'AGG', 'A'),", # anchor is the last CDS base, deleted bases are beyond it
    "(4::UBIGINT, 0::UINTEGER, 300::UBIGINT, 'C', 'A'),",   # intergenic
    "(5::UBIGINT, 0::UINTEGER, 104::UBIGINT, 'C', '<DEL>')",# symbolic
    ") AS v(event_index, seq_region, \"position\", reference, alternate)")
  pairs <- rduckvep_coding_transcripts(con, events, "tinytest_discovery")
  expect_true(is.data.frame(pairs))
  expect_equal(sort(pairs$event_index), c(1, 2))
  expect_true(all(pairs$transcript_index == 0))
  expect_true(rduckvep_coding_transcripts(con, events, "tinytest_discovery", table_name = "tinytest_pairs"))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM tinytest_pairs")$n, 2)
  expect_error(rduckvep_coding_transcripts(con, events, "tinytest_discovery", table_name = "tinytest_pairs"), "already exists")
  expect_error(rduckvep_coding_transcripts(con, events, "no_such_model"), "unknown model")
  expect_error(rduckvep_coding_transcripts(con, "", "tinytest_discovery"), "events_query")
  expect_true(DBI::dbGetQuery(con, "SELECT duckvep_model_drop('tinytest_discovery') AS dropped")$dropped)
})
