library(tinytest)
library(DBI)

local({
  extension <- Sys.getenv("DUCKVEP_EXTENSION")
  con <- rduckvep_connect(extension_path = if (nzchar(extension)) extension else NULL)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  reference <- "ATGGCTGCTGAAGGTTAA"
  quote <- function(x) as.character(dbQuoteString(con, x))
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,117::UBIGINT transcript_end,(1-2*i)::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "117::UBIGINT cds_end,'ATGGCTGCTGAAGGTTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2) t(i)")
  ex <- paste("SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,117::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,18::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM range(2) t(i)")
  expect_true(rduckvep_load_model(con, "single_raw",
    "SELECT i::UINTEGER seq_region FROM range(2) t(i)", tx, ex))
  columns <- c("transcript_index", "cds", "protein", "sequence_flags", "evidence_flags",
    "projection_status", "sequence_status", "edit_count", "carrier_count", "carriers",
    "contributors", "coding_blocks", "cds_differences", "protein_differences", "stop_in_displaced_frame",
    "hgvsc", "hgvsc_status", "hgvsp", "hgvsp_status", "prediction_policy", "prediction_status",
    "prediction_reason", "contributor_provenance", "normalized_edits", "carrier_predictions",
    "haplotype_consequences", "haplotype_impact", "nmd_rule", "nmd_prediction", "nmd_stop_position",
    "nmd_junction_position", "nmd_contributors", "nmd_exceptions", "prediction_reference_protein",
    "prediction_protein", "nominal_length_diff")
  calls <- function(gt, strand, multiallelic = FALSE, reverse = FALSE, starts = c(4L, 13L, 7L), swap_alts = FALSE, unused_alt = FALSE) {
    # Each column is one original record; each row is one sample.
    gt <- as.matrix(gt)
    rows <- character()
    for (sample in seq_len(nrow(gt))) for (j in seq_len(ncol(gt))) {
      position <- if (strand == 1L) 99L + starts[j] else 118L - starts[j]
      ref <- if (strand == 1L) "G" else "C"
      alt <- if (j == 1L) "A" else if (j == 2L) "C" else "T"
      alts <- if (j == 1L && multiallelic) if (swap_alts) c("A", "T") else c("T", "A") else alt
      if (j == 1L && multiallelic && unused_alt) alts <- c(alts, "C")
      if (strand == -1L) alts <- chartr("ACGT", "TGCA", alts)
      rows <- c(rows, sprintf("(%d,%d,%d,%s,[%s],%d,%d,%s,%s)", j,
        as.integer(strand == -1L), position, quote(ref), paste(quote(alts), collapse = ","),
        as.integer(strand == -1L), sample - 1L, quote(gt[sample, j]), quote(paste(sample, j))))
    }
    if (reverse) rows <- rev(rows)
    paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternates,transcript_index,sample_index,gt,call_tag)")
  }
  run <- function(query) {
    dbExecute(con, paste("CREATE OR REPLACE TABLE raw_public_input AS", query))
    original <- dbGetQuery(con, "SELECT * FROM raw_public_input ORDER BY sample_index,event_index")
    result <- rduckvep_haplotypes(con, "SELECT * FROM raw_public_input", "single_raw",
      phase_policy = "vep_compat", input_mode = "source_records", hgvs = TRUE)
    expect_identical(names(result), columns)
    expect_identical(dbGetQuery(con,
      "SELECT * FROM raw_public_input ORDER BY sample_index,event_index"), original)
    result
  }
  copies <- function(result) {
    do.call(rbind, lapply(seq_len(nrow(result)), function(i) {
      cbind(result$carriers[[i]], cds = result$cds[i],
        status = result$carrier_predictions[[i]]$prediction_status,
        reason = result$carrier_predictions[[i]]$prediction_reason)
    }))
  }
  expected <- function(gt, multiallelic = FALSE, starts = c(4L, 13L, 7L), swap_alts = FALSE) {
    # Apply original allele slots to the model, independently of emitted leaves.
    dna <- rep(reference, 2L)
    for (j in seq_along(gt)) {
      alleles <- as.integer(strsplit(gt[j], "[/|]")[[1L]])
      alts <- if (j == 1L && multiallelic) if (swap_alts) c("A", "T") else c("T", "A") else
        if (j == 1L) "A" else if (j == 2L) "C" else "T"
      for (lane in 1:2) if (alleles[lane] != 0L) {
        substr(dna[lane], starts[j], starts[j]) <- alts[alleles[lane]]
      }
    }
    dna
  }
  check <- function(result, dna, sample = 0L) {
    lanes <- copies(result)
    lanes <- lanes[lanes$sample_index == sample, ]
    expect_equal(nrow(lanes), 2L)
    expect_equal(sort(lanes$haplotype_lane), 1:2)
    expect_equal(sort(lanes$cds), sort(dna))
    expect_equal(lanes$cds[order(lanes$haplotype_lane)], dna)
    expect_true(all(lanes$ploidy == 2L & is.na(lanes$phase_set)))
    expect_true(all(lanes$status == "predicted" & lanes$reason == "supported_domain"))
  }
  expect_identical(expected(c("0/1", "1/1")),
    c("ATGGCTGCTGAACGTTAA", "ATGACTGCTGAACGTTAA"))
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) {
    for (spelling in c("0/1", "1/0", "0|1", "1|0")) {
      gt <- c(spelling, "1/1", "0|0")
      query <- calls(matrix(gt, nrow = 1L), strand, reverse = reverse)
      result <- run(query)
      check(result, expected(gt))
      expect_true(all(result$nominal_length_diff == 0))
      # HomoREF is background, not an ALT contributor.
      expect_true(all(vapply(result$contributors, function(x) all(x$event_index != 3), TRUE)))
    }
    for (spelling in c("0/2", "2/0")) {
      gt <- c(spelling, "1|1")
      result <- run(calls(matrix(gt, nrow = 1L), strand, TRUE, reverse))
      check(result, expected(gt, TRUE))
      expect_true(all(vapply(result$contributors, function(x) {
        all(x$alt_index[x$event_index == 1] == 2)
      }, TRUE)))
    }
    # The heterozygote may follow the homoALT in transcript order as well.
    gt <- c("0/1", "1/1")
    starts <- c(13L, 4L)
    check(run(calls(matrix(gt, nrow = 1L), strand, reverse = reverse, starts = starts)),
      expected(gt, starts = starts))
  }
  bad <- c("0/1", "1", "1/1/1", "./1", "|0/1", "/1|0", "0", "0/0/0", "./.")
  for (strand in c(1L, -1L)) for (before in c(FALSE, TRUE)) for (spelling in bad) {
    invalid <- if (before) c(spelling, "0/1") else c("0/1", spelling)
    gt <- rbind(c("0/1", "1/1"), invalid, c("0/0", "1/1"), c("1/1", "1/1"))
    query <- calls(gt, strand)
    result <- run(query)
    for (sample in c(0L, 2L, 3L)) check(result, expected(gt[sample + 1L, ]), sample)
    lanes <- copies(result)
    refused <- lanes[lanes$sample_index == 1L, ]
    expect_true(nrow(refused) > 0L)
    expect_true(all(refused$status == "unsupported_context"))
    expect_true(all(refused$reason == "non_strict_phase_policy"))
  }
  for (strand in c(1L, -1L)) {
    # Original ALT ordinals, not sorted output paths, determine real slots.
    for (spelling in c("1/2", "2/1", "1|2", "2|1"))
      for (swap in c(FALSE, TRUE)) for (unused in c(FALSE, TRUE)) {
      remapped <- if (swap) chartr("12", "21", spelling) else spelling
      gt <- c(remapped, "1/1", "0/0")
      result <- run(calls(matrix(gt, nrow = 1L), strand, TRUE,
        swap_alts = swap, unused_alt = unused))
      dna <- expected(gt, TRUE, swap_alts = swap)
      check(result, dna)
      residue <- vapply(dna, function(x)
        switch(substr(x, 4L, 6L), TCT = "Ser", ACT = "Thr"), "")
      protein <- paste0("M", substr(residue, 1L, 1L), "AER*")
      expect_equal(sort(result$protein), sort(unname(protein)))
      expect_true(all(result$hgvsc_status == "unproven_cis" & is.na(result$hgvsc)))
      expect_equal(sort(result$hgvsp), sort(paste0("p.[(Ala2", residue, ";Gly5Arg)]")))
      expect_true(all(result$hgvsp_status == "ok"))
      expect_true(all(result$nominal_length_diff == 0))
      lanes <- copies(result)
      expect_equal(lanes$cds[order(lanes$haplotype_lane)],
        expected(c(spelling, "1/1", "0/0"), TRUE))
      source_slots <- as.integer(strsplit(gt[1L], "[/|]")[[1L]])
      for (i in seq_len(nrow(result))) {
        source <- result$contributors[[i]]
        lane <- result$carriers[[i]]$haplotype_lane
        expect_equal(source$alt_index[source$event_index == 1L], source_slots[lane])
        expect_equal(sort(source$event_index), c(1, 2))
      }
    }
    for (spelling in c("1/2", "2/1")) {
      gt <- c(spelling, "1/1")
      starts <- c(13L, 4L)
      check(run(calls(matrix(gt, nrow = 1L), strand, TRUE, reverse = TRUE,
        starts = starts)), expected(gt, TRUE, starts = starts))
    }
    for (homo in c("1/1", "2/2", "3/3")) {
      gt <- c(homo, "1/1")
      result <- run(calls(matrix(gt, nrow = 1L), strand, TRUE, unused_alt = TRUE))
      dna <- if (homo == "3/3") expected(c("1/1", "1/1"), TRUE) else expected(gt, TRUE)
      if (homo == "3/3") substr(dna, 4L, 4L) <- "C"
      check(result, dna)
    }
    # Pipes alone carry neither an original PS identity nor a cross-record slot proof.
    phased <- run(calls(matrix(c("0|1", "0|1"), nrow = 1L), strand))
    phased_lanes <- copies(phased)
    expect_true(all(phased_lanes$status == "unsupported_context"))
    expect_true(all(phased_lanes$reason == "non_strict_phase_policy"))
    # One sample's refusal cannot censor a different sample's original copies.
    result <- run(calls(rbind(c("1/2", "1/1"), c("1/2", "0/1")), strand, TRUE))
    check(result, expected(c("1/2", "1/1"), TRUE))
    refused <- copies(result)
    refused <- refused[refused$sample_index == 1L, ]
    expect_true(all(refused$status == "unsupported_context"))
    # An ALT-only leaf does not acquire a second ALT carrier for the absent REF lane.
    gt <- "0/1"
    result <- run(calls(matrix(gt, nrow = 1L), strand))
    lanes <- copies(result)
    alt <- lanes[lanes$haplotype_lane == 2L, ]
    expect_equal(nrow(alt), 1L)
    expect_equal(alt$cds, expected(gt)[expected(gt) != reference])
    expect_true(all(alt$status == "predicted"))
    ref <- lanes[lanes$haplotype_lane == 1L, ]
    expect_true(all(ref$cds == reference))
    expect_true(all(ref$status == "unsupported_context"))
    expect_equal(sum(lanes$cds == alt$cds), 1L)
  }
  own <- function(result, sample) {
    rows <- list()
    for (i in seq_len(nrow(result))) {
      carriers <- result$carriers[[i]]
      for (j in which(carriers$sample_index == sample)) {
        one <- result[i, , drop = FALSE]
        one$carrier_count[] <- 1
        one$carriers[[1L]] <- carriers[j, , drop = FALSE]
        one$carrier_predictions[[1L]] <- one$carrier_predictions[[1L]][j, , drop = FALSE]
        rownames(one) <- rownames(one$carriers[[1L]]) <- rownames(one$carrier_predictions[[1L]]) <- NULL
        rows[[as.character(carriers$haplotype_lane[j])]] <- one
      }
    }
    rows[order(names(rows))]
  }
  for (strand in c(1L, -1L)) {
    query <- calls(matrix(c("1/2", "1/1"), nrow = 1L), strand, TRUE)
    alone <- run(query)
    mixed <- run(calls(rbind(c("1/2", "1/1"), c("1/2", "0/1"), c("1/1", "1/1")), strand, TRUE))
    expect_identical(own(mixed, 0L), own(alone, 0L))
    expect_error(run(paste0("SELECT * REPLACE(1 AS event_index) FROM (", query, ")")),
      "inconsistent source record identity")
    expect_error(run(paste0("SELECT * FROM (", query, ") WHERE event_index=1 UNION ALL ",
      "SELECT * REPLACE('2/1' AS gt) FROM (", query, ") WHERE event_index=1")),
      "inconsistent source record identity or source GT")
    expect_identical(run(query), alone)
  }
  query <- calls(rbind(c("1/2", "1/1"), c("2/1", "1/1"), c("0/1", "0/1")), 1L, TRUE)
  baseline <- run(query)
  limited <- paste0("SELECT * FROM duckvep_haplotypes(", quote(query), ", 'single_raw',",
    "phase_policy='vep_compat',input_mode='source_records',hgvs=true) LIMIT 1")
  expect_equal(nrow(dbGetQuery(con, limited)), 1L)
  expect_equal(run(query), baseline)
  expect_equal(run(query), baseline)
})
