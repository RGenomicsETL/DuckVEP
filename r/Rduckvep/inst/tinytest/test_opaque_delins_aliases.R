library(tinytest)
library(DBI)

local({
  extension <- Sys.getenv("DUCKVEP_EXTENSION")
  con <- rduckvep_connect(extension_path = if (nzchar(extension)) extension else NULL)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  reference <- "ATGGCTGCTGAAGGTTAA"
  quote <- function(x) as.character(dbQuoteString(con, x))
  rc <- function(x) {
    bases <- strsplit(x, "", fixed = TRUE)[[1L]]
    paste(rev(chartr("ACGTacgt", "TGCAtgca", bases)), collapse = "")
  }
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,117::UBIGINT transcript_end,(1-2*i)::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "117::UBIGINT cds_end,'ATGGCTGCTGAAGGTTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2) t(i)")
  ex <- paste("SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,117::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,18::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM range(2) t(i)")
  expect_true(rduckvep_load_model(con, "opaque_aliases",
    "SELECT i::UINTEGER seq_region FROM range(2) t(i)", tx, ex))
  calls <- function(sources, strand, reverse = FALSE) {
    position <- if (strand == 1L) 99L + sources$start else 118L - sources$start - nchar(sources$ref) + 1L
    ref <- if (strand == 1L) sources$ref else vapply(sources$ref, rc, "")
    alt <- if (strand == 1L) sources$alt else vapply(sources$alt, rc, "")
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      sources$event, as.integer(strand == -1L), position, quote(ref), quote(alt), as.integer(strand == -1L))
    if (reverse) rows <- rev(rows)
    paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
  }
  biological <- c("cds", "protein", "sequence_flags", "projection_status", "sequence_status",
    "prediction_status", "prediction_reason", "haplotype_consequences", "haplotype_impact",
    "nmd_prediction", "nmd_rule", "nmd_stop_position", "nmd_junction_position", "nmd_exceptions",
    "prediction_reference_protein", "prediction_protein", "nominal_length_diff")
  witnesses <- data.frame(start = c(4L, 4L, 4L, 4L, 1L, 16L),
    ref = c("GCTGCT", "GCT", "GCT", "GCT", "ATG", "TAA"),
    alt = c("AAA", "TTTCCC", "TAAGAA", "AA", "GCTCCC", "GGGCCC"),
    protein = c("MKEG*", "MFPAEG*", "M*", "MKLKV", "APAAEG*", "MAAEGGP"))
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) {
    sources <- data.frame(event = 1:2, start = c(3L, 4L),
      ref = c("G", "GCT"), alt = c("GCCC", "TTTCCC"))
    result <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(result$prediction_status, "predicted")
    expect_identical(result$cds, "ATGCCCTTTCCCGCTGAAGGTTAA")
    expect_equal(result$nominal_length_diff, 6)
    expect_equal(sort(unique(unlist(result$coding_blocks[[1L]]$event_indices))), 1:2)
    expect_equal(result$contributors[[1L]]$position[order(result$contributors[[1L]]$event_index)],
      if (strand == 1L) c(102, 103) else c(115, 112))
  }
  for (strand in c(1L, -1L)) for (i in seq_len(nrow(witnesses))) {
    w <- witnesses[i, ]
    source <- data.frame(event = 1L, start = w$start, ref = w$ref, alt = w$alt)
    control <- rduckvep_haplotypes(con, calls(source, strand), "opaque_aliases", hgvs = TRUE)
    expected <- paste0(substr(reference, 1L, w$start - 1L), w$alt,
      substr(reference, w$start + nchar(w$ref), nchar(reference)))
    expect_equal(nrow(control), 1L)
    expect_identical(control$cds, expected)
    expect_identical(control$protein, w$protein)
    for (count in c(2L, 3L)) for (padding in c(FALSE, TRUE)) for (reverse in c(FALSE, TRUE)) {
      sources <- source[rep(1L, count), ]
      sources$event <- seq_len(count)
      if (padding && w$start > 1L) {
        anchor <- substr(reference, w$start - 1L, w$start - 1L)
        sources$start[2L] <- sources$start[2L] - 1L
        sources$ref[2L] <- paste0(anchor, sources$ref[2L])
        sources$alt[2L] <- paste0(anchor, sources$alt[2L])
      }
      sources$ref[count] <- tolower(sources$ref[count])
      sources$alt[count] <- tolower(sources$alt[count])
      individual <- lapply(seq_len(count), function(j)
        rduckvep_haplotypes(con, calls(sources[j, ], strand), "opaque_aliases", hgvs = TRUE))
      result <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
      expect_equal(nrow(result), 1L)
      expect_identical(result[biological], control[biological])
      expect_identical(result$prediction_status, "predicted")
      expect_identical(result$hgvsc_status, "ok")
      expect_identical(result$hgvsp_status, "ok")
      original <- result$contributors[[1L]]
      original <- original[order(original$event_index), ]
      expect_equal(original$event_index, sources$event)
      expect_identical(original$reference, if (strand == 1L) sources$ref else unname(vapply(sources$ref, rc, "")))
      expect_identical(original$alternate, if (strand == 1L) sources$alt else unname(vapply(sources$alt, rc, "")))
      expect_equal(result$contributor_provenance[[1L]]$role, rep("applied", count))
      expect_equal(sort(unique(unlist(result$coding_blocks[[1L]]$event_indices))), seq_len(count))
      # A source's standalone operands and nomenclature are unaffected by joint replay.
      for (j in seq_len(count)) {
        single <- rduckvep_haplotypes(con, calls(sources[j, ], strand), "opaque_aliases", hgvs = TRUE)
        expect_identical(single, individual[[j]])
        expect_identical(single[biological], control[biological])
      }
    }
  }

  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) {
    # The opaque payload introduces an early stop; a disjoint later source stays post-stop.
    sources <- data.frame(event = 1:3, start = c(4L, 4L, 13L),
      ref = c("GCT", "GCT", "G"), alt = c("TAAGAA", "TAAGAA", "C"))
    result <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(result$cds, "ATGTAAGAAGCTGAACGTTAA")
    expect_identical(result$protein, "M*")
    roles <- result$contributor_provenance[[1L]]
    expect_identical(roles$role[order(roles$event_index)], c("applied", "applied", "post_stop"))
    expect_identical(result$hgvsp, "p.(Ala2Ter)")
    expect_false(3 %in% result$nmd_contributors[[1L]])
    sources$start[3L] <- 18L
    sources$ref[3L] <- "A"
    sources$alt[3L] <- "AC"
    control <- rduckvep_haplotypes(con, calls(sources[c(1L, 3L), ], strand), "opaque_aliases", hgvs = TRUE)
    result <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(result$cds, "ATGTAAGAAGCTGAAGGTTAAC")
    expect_identical(result[c(biological, "hgvsc", "hgvsp")], control[c(biological, "hgvsc", "hgvsp")])
    roles <- result$contributor_provenance[[1L]]
    expect_identical(roles$role[order(roles$event_index)], c("applied", "applied", "post_stop"))
    expect_equal(sort(result$nmd_contributors[[1L]]), 1:2)
    terminal <- data.frame(event = 1:2, start = 18L, ref = "A", alt = "AC")
    control <- rduckvep_haplotypes(con, calls(terminal[1L, ], strand), "opaque_aliases", hgvs = TRUE)
    result <- rduckvep_haplotypes(con, calls(terminal, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(result$cds, paste0(reference, "C"))
    expect_identical(result[biological], control[biological])
    expect_identical(control$hgvsc_status, "not_compound")
    expect_identical(control$hgvsp_status, "not_applicable")
    expect_identical(result$hgvsc, "c.18_*1insC")
    expect_identical(result$hgvsp, "p.(=)")
    expect_identical(result$contributor_provenance[[1L]]$role, rep("post_stop", 2L))
    terminal$alt[2L] <- "AG"
    refused <- rduckvep_haplotypes(con, calls(terminal, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(refused$prediction_status, "unsupported_overlap")
    terminal$ref[] <- "T"
    terminal$alt[] <- "TC"
    refused <- rduckvep_haplotypes(con, calls(terminal, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(refused$prediction_reason, "reference_mismatch")
    # No base ancestry is assigned to a distinct mutation inside the replacement.
    sources <- data.frame(event = 1:3, start = c(4L, 4L, 5L),
      ref = c("GCT", "GCT", "C"), alt = c("TTTCCC", "TTTCCC", "T"))
    refused <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(refused$prediction_status, "unsupported_overlap")
    expect_true(is.na(refused$cds))
    # An unchanged but incorrect padding base does not establish immutable-reference identity.
    sources <- data.frame(event = 1:2, start = c(4L, 3L),
      ref = c("GCT", "AGCT"), alt = c("TTTCCC", "ATTTCCC"))
    refused <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(refused$prediction_status, "unsupported_context")
    expect_identical(refused$prediction_reason, "reference_mismatch")
    sources <- rbind(sources, data.frame(event = 3L, start = 18L, ref = "A", alt = "AC"))
    refused <- rduckvep_haplotypes(con, calls(sources, strand, reverse), "opaque_aliases", hgvs = TRUE)
    expect_identical(refused$prediction_status, "unsupported_context")
    expect_identical(refused$prediction_reason, "reference_mismatch")
  }

  long_reference <- paste0("ATGAAAGCTGAA", strrep("CTG", 25L), "TAA")
  tx <- paste0("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,189::UBIGINT transcript_end,(1-2*i)::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "189::UBIGINT cds_end,", quote(long_reference), "::BLOB cds_sequence,",
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2) t(i)")
  ex <- paste("SELECT transcript_index::UINTEGER transcript_index,exon_start::UBIGINT exon_start,",
    "exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,",
    "exon_cdna_end::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM",
    "(VALUES (0,100,159,1,60),(0,160,189,61,90),(1,130,189,1,60),(1,100,129,61,90))",
    "v(transcript_index,exon_start,exon_end,exon_cdna_start,exon_cdna_end)")
  expect_true(rduckvep_load_model(con, "opaque_geometry",
    "SELECT i::UINTEGER seq_region FROM range(2) t(i)", tx, ex))
  for (strand in c(1L, -1L)) for (terminal in c(FALSE, TRUE)) for (count in 1:2) {
    ids <- c(seq_len(count), 3L, if (terminal) 4L)
    positions <- c(rep(if (strand == 1L) 103L else 181L, count), if (strand == 1L) 109L else 180L,
      if (terminal) if (strand == 1L) 189L else 100L)
    refs <- c(rep("AAAGCT", count), "G", if (terminal) "A")
    alts <- c(rep("TAA", count), "A", if (terminal) "AC")
    if (strand == -1L) {
      refs <- vapply(refs, rc, "")
      alts <- vapply(alts, rc, "")
    }
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      ids, as.integer(strand == -1L), positions, quote(refs), quote(alts), as.integer(strand == -1L))
    input <- paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    result <- rduckvep_haplotypes(con, input, "opaque_geometry", hgvs = TRUE)
    expect_identical(result$cds, paste0("ATGTAAAAA", strrep("CTG", 25L), "TAA", if (terminal) "C" else ""))
    expect_identical(result$protein, "M*")
    expect_identical(result$prediction_status, "predicted")
    expect_identical(result$nmd_prediction, "trigger")
    expect_equal(result$nmd_stop_position, 6)
    expect_equal(result$nmd_junction_position, 57)
    expect_equal(result$nmd_junction_position - result$nmd_stop_position, 51)
    expect_equal(result$nominal_length_diff, -3 + terminal)
    expect_equal(sort(result$nmd_contributors[[1L]]), seq_len(count))
    roles <- result$contributor_provenance[[1L]]
    expect_identical(roles$role[order(roles$event_index)], c(rep("applied", count), rep("post_stop", 1L + terminal)))
    if (count == 1L) control <- result else {
      expect_identical(result[c(biological, "hgvsc", "hgvsp")], control[c(biological, "hgvsc", "hgvsp")])
      expect_equal(sort(unique(unlist(result$coding_blocks[[1L]]$event_indices))), ids)
    }
  }
  extended_reference <- paste0("ATGAAA", strrep("CTG", 167L), "TAA")
  tx <- paste0("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,609::UBIGINT transcript_end,(1-2*i)::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "609::UBIGINT cds_end,", quote(extended_reference), "::BLOB cds_sequence,",
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2) t(i)")
  ex <- paste("SELECT transcript_index::UINTEGER transcript_index,exon_start::UBIGINT exon_start,",
    "exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,",
    "exon_cdna_end::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM",
    "(VALUES (0,100,159,1,60),(0,160,609,61,510),(1,550,609,1,60),(1,100,549,61,510))",
    "v(transcript_index,exon_start,exon_end,exon_cdna_start,exon_cdna_end)")
  expect_true(rduckvep_load_model(con, "opaque_long_exon",
    "SELECT i::UINTEGER seq_region FROM range(2) t(i)", tx, ex))
  for (strand in c(1L, -1L)) for (count in 1:2) {
    ids <- seq_len(count + 1L)
    starts <- c(rep(58L, count), 61L)
    refs <- c(rep("CTG", count), "CTG")
    alts <- c(rep("CTGAAA", count), "TAA")
    positions <- if (strand == 1L) starts + 99L else 611L - starts - nchar(refs)
    if (strand == -1L) {
      refs <- vapply(refs, rc, "")
      alts <- vapply(alts, rc, "")
    }
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      ids, as.integer(strand == -1L), positions, quote(refs), quote(alts), as.integer(strand == -1L))
    input <- paste0("SELECT * FROM (VALUES ", paste(rev(rows), collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    result <- rduckvep_haplotypes(con, input, "opaque_long_exon", hgvs = TRUE)
    expect_identical(result$cds, paste0(substr(extended_reference, 1L, 60L), "AAA", "TAA",
      substr(extended_reference, 64L, 510L)))
    expect_identical(result$prediction_status, "predicted")
    expect_identical(result$nmd_prediction, "escape")
    expect_equal(result$nmd_stop_position, 66)
    expect_equal(result$nmd_junction_position, 63)
    expect_true(grepl("long_exon", result$nmd_exceptions[[1L]], fixed = TRUE))
    if (count == 1L) long_control <- result else
      expect_identical(result[c(biological, "hgvsc", "hgvsp")], long_control[c(biological, "hgvsc", "hgvsp")])
  }
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) {
    # Equal spliced insertion sites at a junction have different genomic anchors.
    left <- substr(long_reference, 60L, 60L)
    right <- substr(long_reference, 61L, 61L)
    sources <- data.frame(event = 1:2, start = c(60L, 61L),
      ref = c(left, right), alt = c(paste0(left, "CCC"), paste0("CCC", right)))
    positions <- if (strand == 1L) 99L + sources$start else 190L - sources$start
    refs <- if (strand == 1L) sources$ref else vapply(sources$ref, rc, "")
    alts <- if (strand == 1L) sources$alt else vapply(sources$alt, rc, "")
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      sources$event, as.integer(strand == -1L), positions, quote(refs), quote(alts),
      as.integer(strand == -1L))
    input <- function(r) paste0("SELECT * FROM (VALUES ", paste(rows[r], collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    single <- rduckvep_haplotypes(con, input(1L), "opaque_geometry", hgvs = TRUE)
    other <- rduckvep_haplotypes(con, input(2L), "opaque_geometry", hgvs = TRUE)
    expect_identical(single$cds, other$cds)
    expect_identical(single$cds, paste0(substr(long_reference, 1L, 60L), "CCC",
      substr(long_reference, 61L, 90L)))
    pair <- rduckvep_haplotypes(con, input(if (reverse) 2:1 else 1:2), "opaque_geometry", hgvs = TRUE)
    if (strand == 1L) {
      expect_identical(pair$nmd_prediction, "not_applicable")
      expect_identical(pair$prediction_status, "predicted")
      expect_identical(pair$cds, paste0(substr(long_reference, 1L, 60L), "CCCCCC",
        substr(long_reference, 61L, 90L)))
      expect_equal(pair$nominal_length_diff, 6)
      expect_equal(sort(unique(unlist(pair$coding_blocks[[1L]]$event_indices))), 1:2)
    } else {
      expect_identical(pair$prediction_status, "unsupported_overlap")
      expect_identical(pair$prediction_reason, "ambiguous_same_gap_insertions")
      expect_true(is.na(pair$cds))
    }
    duplicate_rows <- rows
    triple_rows <- c(rows, sub("^\\(1,", "(3,", rows[1L]))
    triple <- rduckvep_haplotypes(con, paste0("SELECT * FROM (VALUES ",
      paste(triple_rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)"), "opaque_geometry", hgvs = TRUE)
    if (strand == 1L) {
      expect_identical(triple$cds, pair$cds)
      expect_equal(sort(unique(unlist(triple$coding_blocks[[1L]]$event_indices))), 1:3)
    } else expect_identical(triple$prediction_status, "unsupported_overlap")
    duplicate_rows[2L] <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      2L, as.integer(strand == -1L), positions[1L], quote(refs[1L]), quote(alts[1L]),
      as.integer(strand == -1L))
    duplicates <- rduckvep_haplotypes(con, paste0("SELECT * FROM (VALUES ",
      paste(duplicate_rows[if (reverse) 2:1 else 1:2], collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)"), "opaque_geometry", hgvs = TRUE)
    expect_identical(duplicates$cds, single$cds)
    expect_identical(duplicates$prediction_status, "predicted")
    expect_equal(sort(unique(unlist(duplicates$coding_blocks[[1L]]$event_indices))), 1:2)
  }
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) {
    # The c.60/61 junction is adjacent in the edited CDS but belongs to two exon-local edits.
    ids <- 1:4
    positions <- if (strand == 1L) c(103L, 157L, 157L, 160L) else c(181L, 130L, 130L, 129L)
    refs <- c("AAAGCT", "CTG", "CTG", "C")
    alts <- c("TAA", "AAAAAA", "AAAAAA", "T")
    if (strand == -1L) {
      refs <- vapply(refs, rc, "")
      alts <- vapply(alts, rc, "")
    }
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      ids, as.integer(strand == -1L), positions, quote(refs), quote(alts), as.integer(strand == -1L))
    input <- function(r) paste0("SELECT * FROM (VALUES ", paste(rows[r], collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    baseline <- rduckvep_haplotypes(con, input(c(1L, 2L, 4L)), "opaque_geometry", hgvs = TRUE)
    result <- rduckvep_haplotypes(con, input(if (reverse) rev(ids) else ids), "opaque_geometry", hgvs = TRUE)
    expected <- paste0("ATGTAA", substr(long_reference, 10L, 57L), "AAAAAA", "T",
      substr(long_reference, 62L, 90L))
    expect_identical(result$cds, expected)
    expect_identical(result$cds, baseline$cds)
    expect_identical(result$protein, "M*")
    expect_identical(result$prediction_status, "predicted")
    expect_identical(result$nmd_prediction, "trigger")
    expect_equal(result$nmd_stop_position, 6)
    expect_equal(result$nmd_junction_position, 60)
    expect_equal(result$nominal_length_diff, 0)
    expect_identical(result[c(biological, "hgvsc", "hgvsp")], baseline[c(biological, "hgvsc", "hgvsp")])
    expect_equal(sort(unique(unlist(result$coding_blocks[[1L]]$event_indices))), ids)
    expect_identical(result$contributor_provenance[[1L]]$role[order(result$contributor_provenance[[1L]]$event_index)],
      c("applied", "post_stop", "post_stop", "post_stop"))
    bad_ref <- if (strand == 1L) "G" else rc("G")
    bad_rows <- rows
    bad_rows[4L] <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      4L, as.integer(strand == -1L), positions[4L], quote(bad_ref), quote(alts[4L]),
      as.integer(strand == -1L))
    bad_input <- paste0("SELECT * FROM (VALUES ", paste(bad_rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    mismatched <- rduckvep_haplotypes(con, bad_input, "opaque_geometry", hgvs = TRUE)
    expect_identical(mismatched$prediction_reason, "reference_mismatch")
    expect_true(is.na(mismatched$cds))
  }
})
