library(tinytest)
library(DBI)

local({
  extension_path <- Sys.getenv("DUCKVEP_EXTENSION", unset = "")
  con <- if (nzchar(extension_path)) rduckvep_connect(extension_path = extension_path) else rduckvep_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  regions <- "SELECT i::UINTEGER seq_region FROM range(4) t(i)"
  transcripts <- paste(
    "SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,100::UBIGINT transcript_start,",
    "123::UBIGINT transcript_end,CASE WHEN i=1 THEN -1 ELSE 1 END::TINYINT strand,",
    "i::UINTEGER gene_index,0::UBIGINT transcript_flags,100::UBIGINT cds_start,123::UBIGINT cds_end,",
    "(CASE WHEN i=1 THEN 'TTG' WHEN i=3 THEN 'NTG' ELSE 'CTG' END ||",
    "'GCTGAACTGAAACGCAGCTAA')::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence,NULL::BLOB cdna_sequence",
    "FROM range(4) t(i)")
  exons <- paste(
    "SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,123::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,24::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM range(4) t(i)")
  load <- paste0("SELECT loaded FROM duckvep_model_load('annotated_initiation',",
    paste(vapply(c(regions, transcripts, exons), dbQuoteString, con = con, FUN.VALUE = ""),
      collapse = ","), ")")
  expect_true(dbGetQuery(con, load)$loaded)
  dbExecute(con, paste(
    "CREATE TABLE annotated_calls AS SELECT * FROM (VALUES",
    "(1::BIGINT,0::INTEGER,100::BIGINT,'C','T',1::INTEGER,0::INTEGER),",
    "(2,1,123,'A','G',1,1),(3,0,103,'G','A',1,0),",
    "(4,1,120,'C','T',1,1),(5,0,101,'T','C',1,0),",
    "(6,2,103,'G','A',1,2),(7,3,103,'G','A',1,3))",
    "t(event_index,seq_region,position,reference,alternate,alt_index,transcript_index),",
    "LATERAL (SELECT event_index::INTEGER sample_index,[1,0]::INTEGER[] alleles,",
    "[false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set) gt"))
  got <- rduckvep_haplotypes(con, "SELECT * FROM annotated_calls", "annotated_initiation")
  carriers <- do.call(rbind, got$carrier_predictions)
  carriers <- carriers[order(carriers$sample_index), ]
  expect_identical(carriers$prediction_status,
    c(rep("predicted", 6), "unsupported_context"))
  expect_identical(carriers$prediction_reason,
    c(rep("supported_domain", 6), "noncanonical_start"))
  expect_true("TTGGCTGAACTGAAACGCAGCTAA" %in% got$cds[got$transcript_index == 0])
  expect_true("CTGGCTGAACTGAAACGCAGCTAA" %in% got$cds[got$transcript_index == 1])
  expect_true("LAELKRS*" %in% got$protein[got$transcript_index == 0])
  expect_true("LAELKRS*" %in% got$protein[got$transcript_index == 1])
  snapshot <- tempfile(fileext = ".dvsnap")
  on.exit(unlink(snapshot), add = TRUE)
  expect_true(rduckvep_save_model(con, "annotated_initiation", snapshot))
  expect_true(rduckvep_restore_model(con, "annotated_restored", snapshot))
  restored <- rduckvep_haplotypes(con, "SELECT * FROM annotated_calls", "annotated_restored")
  expect_identical(restored, got)
})
