local({
  con <- rduckvep_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  model <- paste0("SELECT loaded FROM duckvep_model_load('tinytest_calls', ",
    "'SELECT 0::UINTEGER seq_region, 1000::UBIGINT sequence_length, ''chrT'' seq_region_name', ",
    "'SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,111::UBIGINT cds_end,''ATGGCTGCTTAA''::BLOB cds_sequence,1::UTINYINT codon_table,''''::BLOB pre_cds_sequence,''''::BLOB post_cds_sequence', ",
    "'SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase')")
  expect_true(DBI::dbGetQuery(con, model)$loaded)
  vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(vcf), add = TRUE)
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set\">",
    "##contig=<ID=chrT,length=1000>",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", sep = "\t"),
    paste("chrT", 50, ".", "C", "A", ".", ".", ".", "GT:PS", "0|1:1", sep = "\t"),      # outside the CDS
    paste("chrT", 104, ".", "C", "A,G", ".", ".", ".", "GT:PS", "1|2:77", sep = "\t"),  # two ALT alleles
    paste("chrT", 106, ".", "C", "T", ".", ".", ".", "GT:PS", "0/1:PATMAT", sep = "\t"),
    paste("other", 104, ".", "C", "A", ".", ".", ".", "GT:PS", "0|1:1", sep = "\t")),  # unknown contig
    vcf)
  calls <- rduckvep_coding_calls(con, "tinytest_calls", vcf)
  expect_true(is.data.frame(calls))
  expect_equal(nrow(calls), 3L)
  expect_equal(sort(calls$event_index), c(128, 129, 192))
  expect_equal(calls$phase_set[calls$event_index == 128], 77)
  expect_true(is.na(calls$phase_set[calls$event_index == 192]))
  expect_equal(unlist(calls$alleles[calls$event_index == 129]), c(1L, 2L))
  expect_equal(unlist(calls$phase_before[calls$event_index == 192]), c(FALSE, FALSE))
  expect_true(rduckvep_coding_calls(con, "tinytest_calls", vcf, table_name = "tinytest_coding_calls"))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM tinytest_coding_calls")$n, 3)
  expect_error(rduckvep_coding_calls(con, "tinytest_calls", vcf, table_name = "tinytest_coding_calls"), "already exists")
  haplotypes <- rduckvep_haplotypes(con, "SELECT * FROM tinytest_coding_calls", "tinytest_calls")
  expect_true(nrow(haplotypes) > 0L)
  expect_error(rduckvep_coding_calls(con, "tinytest_calls", tempfile(fileext = ".vcf")), "cannot open")
  expect_error(rduckvep_coding_calls(con, "no_such_model", vcf), "require a loaded model")
  expect_error(rduckvep_coding_calls(con, "tinytest_calls", ""), "path")
  expect_true(DBI::dbGetQuery(con, "SELECT duckvep_model_drop('tinytest_calls') AS dropped")$dropped)
})
