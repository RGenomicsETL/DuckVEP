local({
  extension <- Sys.getenv("DUCKVEP_EXTENSION")
  con <- rduckvep_connect(extension_path = if (nzchar(extension)) extension else NULL)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,117::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "117::UBIGINT cds_end,'ATGGCTGCTGAAGGTTAA'::BLOB cds_sequence,1::UTINYINT codon_table")
  ex <- paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,",
    "117::UBIGINT exon_end,1::UBIGINT exon_cdna_start,18::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase")
  expect_true(rduckvep_load_model(con, "copy_proof", "SELECT 0::UINTEGER seq_region", tx, ex))
  calls <- function(gt, sample = 0L) {
    gt <- as.character(DBI::dbQuoteString(con, gt))
    sprintf(paste("SELECT event_index::UBIGINT event_index,0::UINTEGER seq_region,",
      "position::UBIGINT AS position,'G' AS reference,alternates,0::UINTEGER transcript_index,",
      "%d::UINTEGER sample_index,gt FROM (VALUES (1,103,['A'],%s),",
      "(2,112,['C'],%s)) t(event_index,position,alternates,gt)"), sample, gt[1L], gt[2L])
  }
  replay <- function(sql) rduckvep_haplotypes(con, sql, "copy_proof",
    phase_policy = "vep_compat", input_mode = "source_records", hgvs = TRUE)
  complete <- calls(c("1/1", "1|1"))
  alone <- replay(complete)
  expect_identical(alone$prediction_status, "predicted")
  expect_identical(alone$prediction_reason, "supported_domain")
  expect_identical(alone$cds, "ATGACTGCTGAACGTTAA")
  expect_identical(alone$haplotype_consequences[[1L]], "missense_variant")
  expect_equal(nrow(alone$carriers[[1L]]), 2L)

  for (gt in list(c("0/1", "0/1"), c("1", "1"), c("1/1/0", "1/1"), c("1/1", "1/."))) {
    mixed <- replay(paste(complete, "UNION ALL", calls(gt, 1L)))
    selected <- vapply(mixed$carriers, function(x) any(x$sample_index == 0L), logical(1))
    observed <- mixed[selected, , drop = FALSE]
    rownames(observed) <- NULL
    expect_identical(observed, alone)
    for (prediction in mixed$carrier_predictions[!selected]) {
      expect_true(all(prediction$prediction_status == "unsupported_context"))
      expect_true(all(prediction$prediction_reason == "non_strict_phase_policy"))
    }
  }
  for (ref in c("0", "0/0/0", "0/0")) {
    for (gt in list(c(ref, "1/1"), c("1/1", ref))) {
      value <- replay(calls(gt))
      expected <- if (ref == "0/0") "predicted" else "unsupported_context"
      expect_true(all(value$prediction_status == expected))
      expect_true(all(vapply(value$carrier_predictions,
        function(x) all(x$prediction_status == expected), logical(1))))
    }
  }
  expect_equal(nrow(replay(calls(c("0/0", "0/0")))), 0L)
  expect_error(replay(calls(c("2/2", "1/1"))), "raw GT status 3")
  leading <- replay(calls(c("|1|1", "1/1")))
  expect_true(all(leading$prediction_status == "unsupported_context"))
  expect_true(all(leading$prediction_reason == "non_strict_phase_policy"))

  DBI::dbExecute(con, paste("CREATE TABLE copy_input AS", complete))
  before <- DBI::dbGetQuery(con, "SELECT * FROM copy_input ORDER BY event_index")
  expect_identical(replay("SELECT * FROM copy_input"), alone)
  expect_identical(DBI::dbGetQuery(con, "SELECT * FROM copy_input ORDER BY event_index"), before)
})
