local({
  con <- rduckvep_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  model <- paste0("SELECT loaded FROM duckvep_model_load('tinytest', ",
    "'SELECT 0::UINTEGER seq_region', ",
    "'SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,111::UBIGINT cds_end,''ATGGCTGCTTAA''::BLOB cds_sequence,1::UTINYINT codon_table,''''::BLOB pre_cds_sequence,''''::BLOB post_cds_sequence', ",
    "'SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase')")
  expect_true(DBI::dbGetQuery(con, model)$loaded)
  calls <- paste("SELECT 1::UBIGINT event_index,0::UINTEGER seq_region,",
    "104::UBIGINT AS \"position\",'C' AS reference,'A' AS alternate,",
    "1::UINTEGER alt_index,0::UINTEGER transcript_index,",
    "0::UBIGINT sample_index,[0,1]::INTEGER[] alleles,",
    "[false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set")
  haplotypes <- rduckvep_haplotypes(con, calls, "tinytest")
  expect_true(is.data.frame(haplotypes))
  expect_true("carriers" %in% names(haplotypes))
  expect_true(nrow(haplotypes) > 0L)
  expect_true(DBI::dbGetQuery(con, "SELECT duckvep_model_drop('tinytest') AS dropped")$dropped)
})
