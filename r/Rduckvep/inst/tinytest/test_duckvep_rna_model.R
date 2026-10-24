library(tinytest)
library(DBI)

local({
  artifact <- Sys.getenv("DUCKVEP_EXTENSION", "")
  con <- rduckvep_connect(extension_path = if (nzchar(artifact)) artifact else NULL)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  regions <- "SELECT 0::UINTEGER seq_region, 12::UBIGINT sequence_length"
  transcripts <- paste(
    "SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "1::UBIGINT transcript_start,12::UBIGINT transcript_end,-1::TINYINT strand,",
    "0::UINTEGER gene_index,259::UBIGINT transcript_flags,",
    "4::UBIGINT cds_start,12::UBIGINT cds_end,NULL::BLOB cds_sequence,",
    "NULL::UTINYINT codon_table,NULL::BLOB pre_cds_sequence,",
    "NULL::BLOB post_cds_sequence,'ATGCCCTAAAGG'::BLOB cdna_sequence,",
    "1::UBIGINT rna_cds_cdna_start,9::UBIGINT rna_cds_cdna_end,",
    "2::UTINYINT rna_codon_table")
  exons <- paste("SELECT 0::UINTEGER transcript_index,1::UBIGINT exon_start,",
    "12::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase")
  programme <- paste("SELECT 0::UINTEGER transcript_index,4::UINTEGER cdna_start,",
    "5::UINTEGER cdna_end,'CC'::VARCHAR expected_substrate,",
    "'TT'::VARCHAR alternate_sequence,'substitution'::VARCHAR operation_kind,",
    "'Ensembl116:_rna_edit'::VARCHAR programme_name,",
    "'transcript_attrib:SeqEdit'::VARCHAR evidence_name,",
    "0::UTINYINT boundary_overrides")
  alleles <- paste("SELECT 41::BIGINT AS event_index,0::UINTEGER AS seq_region,6::BIGINT AS position,",
    "'A'::VARCHAR AS reference,'G'::VARCHAR AS alternate,1::INTEGER AS alt_index,",
    "0::UINTEGER AS transcript_index,1::INTEGER AS sample_index,",
    "[1,1]::INTEGER[] AS alleles,[false,true]::BOOLEAN[] AS phase_before,NULL::BIGINT AS phase_set")
  expect_true(rduckvep_load_model(con, "rna_wrapper", regions, transcripts, exons,
    rna_edit_query = programme))
  fingerprint <- dbGetQuery(con,
    "SELECT _duckvep_model_fingerprint('rna_wrapper')::VARCHAR AS fingerprint")
  result <- rduckvep_haplotypes(con, alleles, "rna_wrapper")
  expect_equal(nrow(result), 1L)
  # The programme is part of the transcript, as VEP reads it: MF* to MFQ* loses the stop.
  expect_identical(result$prediction_status, "predicted")
  expect_identical(result$prediction_reason, "supported_domain")
  expect_identical(result$prediction_reference_protein, "MF*")
  expect_identical(result$prediction_protein, "MFQ*")
  expect_true(all(is.na(result$cds)))
  expect_true(all(is.na(result$protein)))
  expect_identical(result$haplotype_consequences[[1L]], "stop_lost")
  expect_identical(result$haplotype_impact, "HIGH")
  expect_identical(result$nmd_prediction, "unknown")
  expect_identical(dbGetQuery(con,
    "SELECT _duckvep_model_fingerprint('rna_wrapper')::VARCHAR AS fingerprint"), fingerprint)
  snapshot <- tempfile(fileext = ".dvsnap")
  on.exit(unlink(snapshot), add = TRUE)
  expect_true(rduckvep_save_model(con, "rna_wrapper", snapshot))
  expect_true(rduckvep_restore_model(con, "rna_wrapper_copy", snapshot))
  expect_identical(rduckvep_haplotypes(con, alleles, "rna_wrapper_copy"), result)
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('rna_wrapper')")[[1L]])
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('rna_wrapper_copy')")[[1L]])
})
