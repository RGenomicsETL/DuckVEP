library(tinytest)
library(DBI)

local({
  extension <- Sys.getenv("DUCKVEP_EXTENSION")
  con <- rduckvep_connect(extension_path = if (nzchar(extension)) extension else NULL)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  reference <- "ATGGCTGAATTTTAA"
  # A reference-position assignment oracle validates each source against the immutable reference.
  union_cds <- function(position, ref, alt) {
    bases <- strsplit(reference, "", fixed = TRUE)[[1L]]
    assigned <- rep(NA_character_, length(bases))
    for (i in seq_along(position)) {
      at <- position[i] + seq_len(nchar(ref[i])) - 1L
      stopifnot(identical(paste(bases[at], collapse = ""), ref[i]))
      value <- strsplit(alt[i], "", fixed = TRUE)[[1L]]
      changed <- bases[at] != value
      at <- at[changed]; value <- value[changed]
      if (any(!is.na(assigned[at]) & assigned[at] != value)) return(NA_character_)
      assigned[at] <- value
    }
    bases[!is.na(assigned)] <- assigned[!is.na(assigned)]
    paste(bases, collapse = "")
  }
  complement <- c(A = "T", C = "G", G = "C", T = "A")
  rc <- function(x) paste(rev(unname(complement[strsplit(x, "", fixed = TRUE)[[1L]]])), collapse = "")
  quote <- function(x) as.character(dbQuoteString(con, x))
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,114::UBIGINT transcript_end,(1-2*i)::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "114::UBIGINT cds_end,'ATGGCTGAATTTTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2) t(i)")
  exon <- paste("SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,114::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,15::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM range(2) t(i)")
  expect_equal(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('compatible_overlap',",
    "'SELECT i::UINTEGER seq_region FROM range(2) t(i)',", quote(tx), ",", quote(exon), ")"))$loaded, TRUE)
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) for (scenario in c("target", "implicit_ps", "redundant", "shared_unchanged", "multi_block", "multi_redundant", "split_same_block", "padded_adjacent", "padded_disjoint", "neutral_inside_overlap", "conflict", "indel", "unphased", "partial", "reuse")) {
    position <- c(4L, 5L); ref <- c("GC", "CT"); alt <- c("AT", "TT")
    if (scenario == "redundant") { position[2L] <- 4L; ref[2L] <- "GC"; alt[2L] <- "AT" }
    if (scenario == "shared_unchanged") { position <- c(4L, 4L); ref <- c("GCT", "GCT"); alt <- c("ATT", "ATT") }
    if (scenario %in% c("multi_block", "multi_redundant")) {
      ref[1L] <- "GCTGAATT"; alt[1L] <- "ATTGAACT"
      if (scenario == "multi_redundant") {
        position <- c(position, position[1L]); ref <- c(ref, ref[1L]); alt <- c(alt, alt[1L])
      }
    }
    if (scenario == "split_same_block") { position <- c(4L, 4L); ref <- c("GCT", "G"); alt <- c("ACA", "A") }
    if (scenario == "padded_adjacent") alt[2L] <- "CA"
    if (scenario == "padded_disjoint") alt[1L] <- "AC"
    if (scenario == "neutral_inside_overlap") { ref[1L] <- "GCT"; alt <- c("ACA", "TA") }
    if (scenario == "conflict") alt[2L] <- "GT"
    if (scenario == "indel") alt[2L] <- "T"
    genomic_position <- if (strand == 1L) 99L + position else 116L - position - nchar(ref)
    genomic_ref <- if (strand == 1L) ref else vapply(ref, rc, "")
    genomic_alt <- if (strand == 1L) alt else vapply(alt, rc, "")
    phase <- if (scenario == "unphased") "[false,false]" else "[false,true]"
    alleles <- if (scenario == "partial") "[NULL,1]::INTEGER[]" else "[0,1]"
    phase_set <- if (scenario == "implicit_ps") "NULL::BIGINT" else "77::BIGINT"
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,%s,%s,%s)",
      seq_along(position), as.integer(strand == -1L), genomic_position, quote(genomic_ref), quote(genomic_alt),
      as.integer(strand == -1L), alleles, phase, phase_set)
    if (reverse) rows <- rev(rows)
    calls <- paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    result <- rduckvep_haplotypes(con, calls, "compatible_overlap", phase_policy = "strict", hgvs = TRUE)
    original <- dbGetQuery(con, calls)
    expect_equal(original$alleles, rep(list(if (scenario == "partial") c(NA_integer_, 1L) else c(0L, 1L)), length(position)))
    expect_equal(original$phase_before, rep(list(if (scenario == "unphased") c(FALSE, FALSE) else c(FALSE, TRUE)), length(position)))
    expect_equal(original$phase_set, rep(if (scenario == "implicit_ps") NA_real_ else 77, length(position)))
    expect_equal(dbGetQuery(con, calls), original)
    for (i in seq_len(nrow(result))) {
      carriers <- result$carriers[[i]]
      expect_equal(carriers$sample_index, rep(0, nrow(carriers)))
      expect_equal(carriers$phase_set, rep(if (scenario %in% c("unphased", "implicit_ps")) NA_real_ else 77, nrow(carriers)))
      expect_equal(carriers$ploidy, rep(2, nrow(carriers)))
      expect_true(all(carriers$haplotype_lane %in% 1:2))
      if (!scenario %in% c("partial", "unphased")) expect_equal(carriers$haplotype_lane, 2)
      sources <- result$contributors[[i]]
      sources <- sources[order(sources$event_index), ]
      expect_equal(sources$event_index, seq_along(position))
      expect_equal(sources$position, genomic_position)
      expect_identical(sources$reference, unname(genomic_ref))
      expect_identical(sources$alternate, unname(genomic_alt))
      if (scenario != "partial") {
        expected_evidence <- if (scenario == "unphased") 4 else 1
        expect_equal(sources$evidence_flags, rep(expected_evidence, length(position)))
        expect_equal(result$evidence_flags[i], expected_evidence)
      } else expect_true(all(bitwAnd(sources$evidence_flags, result$evidence_flags[i]) == sources$evidence_flags))
    }
    if (scenario == "target" && !reverse) {
      duplicate <- paste0(calls, " UNION ALL ", calls)
      expect_error(rduckvep_haplotypes(con, duplicate, "compatible_overlap", phase_policy = "strict"))
    }
    if (scenario == "partial") {
      known <- result$prediction_status == "predicted"
      expect_equal(sum(known), 1L)
      expect_true(all(result$prediction_status[!known] == "incomplete_input"))
      expect_true(all(result$prediction_reason[!known] == "missing_call"))
      expect_true(all(is.na(result$cds[!known])))
      result <- result[known, , drop = FALSE]
    }
    admitted <- scenario %in% c("target", "implicit_ps", "redundant", "shared_unchanged", "multi_block", "multi_redundant", "split_same_block", "padded_adjacent", "padded_disjoint", "neutral_inside_overlap", "partial", "reuse")
    if (admitted) {
      expect_equal(nrow(result), 1L)
      expect_identical(result$cds, union_cds(position, ref, alt))
      multi <- scenario %in% c("multi_block", "multi_redundant")
      split <- scenario == "split_same_block"
      triplet <- scenario %in% c("padded_adjacent", "neutral_inside_overlap")
      expect_identical(result$protein, if (multi) "MIEL*" else if (split) "MTEF*" else "MIEF*")
      expect_identical(result$hgvsc, if (multi) "c.[4_5delinsAT;10T>C]" else if (split) "c.4_6delinsACA" else if (triplet) "c.4_6delinsATA" else "c.4_5delinsAT")
      expect_identical(result$hgvsp, if (multi) "p.[(Ala2Ile;Phe4Leu)]" else if (split) "p.(Ala2Thr)" else "p.(Ala2Ile)")
      blocks <- result$coding_blocks[[1L]]
      expect_equal(nrow(blocks), if (multi) 2 else 1)
      expect_equal(sort(blocks$event_indices[[1L]]), seq_along(position))
      if (multi) {
        expect_equal(sort(blocks$event_indices[[2L]]), if (scenario == "multi_redundant") c(1, 3) else 1)
        edits <- result$normalized_edits[[1L]]
        expect_equal(edits$block_index, ifelse(edits$cds_start < 7, 0, 1))
      }
      expect_true(all(vapply(blocks$event_indices, function(x) !anyDuplicated(x), TRUE)))
    } else {
      expect_true(all(is.na(result$cds)))
      expect_true(all(vapply(result$haplotype_consequences, function(x) !length(x) || all(is.na(x)), TRUE)))
      expect_true(all(result$nmd_prediction == "unknown"))
      expect_true(all(is.na(result$hgvsc)))
      expect_true(all(is.na(result$hgvsp)))
      if (scenario == "conflict") expect_identical(result$prediction_status, "edit_conflict")
      expect_true(all(vapply(result$contributors, nrow, 0L) == 2L))
    }
  }
  # Unchanged context between a source's edit islands is not block membership.
  for (strand in c(1L, -1L)) for (lower in c(FALSE, TRUE)) for (reverse in c(FALSE, TRUE)) {
    position <- c(4L, 5L, 7L)
    ref <- c("GCTGAATT", "CT", "G"); alt <- c("ATTGAACT", "TT", "A")
    genomic_position <- if (strand == 1L) 99L + position else 116L - position - nchar(ref)
    genomic_ref <- if (strand == 1L) ref else vapply(ref, rc, "")
    genomic_alt <- if (strand == 1L) alt else vapply(alt, rc, "")
    if (lower) { genomic_ref <- tolower(genomic_ref); genomic_alt <- tolower(genomic_alt) }
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      seq_along(position), as.integer(strand == -1L), genomic_position, quote(genomic_ref), quote(genomic_alt),
      as.integer(strand == -1L))
    if (reverse) rows <- rev(rows)
    calls <- paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,sample_index,alleles,phase_before,phase_set)")
    result <- rduckvep_haplotypes(con, calls, "compatible_overlap", hgvs = TRUE)
    expect_identical(result$cds, union_cds(position, ref, alt))
    expect_identical(result$protein, "MIKL*")
    expect_equal(lapply(result$coding_blocks[[1L]]$event_indices, sort), list(c(1, 2), 3, 1))
    edits <- result$normalized_edits[[1L]]
    expect_equal(edits$block_index[edits$event_index == 1], c(0, 2))
    expect_equal(edits$block_index[edits$event_index == 2], 0)
    expect_equal(edits$block_index[edits$event_index == 3], 1)
  }
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) for (mixed in c(FALSE, TRUE)) {
    position <- c(4L, 4L, 8L)
    ref <- c(if (mixed) "GCTGA" else "GCT", "GCT", "A")
    alt <- c(if (mixed) "TAAGC" else "TAA", "TAA", "C")
    genomic_position <- if (strand == 1L) 99L + position else 116L - position - nchar(ref)
    genomic_ref <- if (strand == 1L) ref else vapply(ref, rc, "")
    genomic_alt <- if (strand == 1L) alt else vapply(alt, rc, "")
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)",
      seq_along(position), as.integer(strand == -1L), genomic_position, quote(genomic_ref), quote(genomic_alt),
      as.integer(strand == -1L))
    if (reverse) rows <- rev(rows)
    calls <- paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,sample_index,alleles,phase_before,phase_set)")
    result <- rduckvep_haplotypes(con, calls, "compatible_overlap", hgvs = TRUE)
    expect_identical(result$prediction_status, "predicted")
    expect_identical(result$cds, union_cds(position, ref, alt))
    expect_identical(result$protein, "M*")
    provenance <- result$contributor_provenance[[1L]]
    expect_identical(provenance$role[order(provenance$event_index)], c("applied", "applied", "post_stop"))
    expect_identical(result$nmd_prediction, "escape")
    expect_equal(sort(result$nmd_contributors[[1L]]), c(1, 2))
    edits <- result$normalized_edits[[1L]]
    expect_equal(edits$block_index[edits$event_index == 3], 1)
    expect_equal(lapply(result$coding_blocks[[1L]]$event_indices, sort),
      list(c(1, 2), if (mixed) c(1, 3) else 3))
  }
  lower_calls <- paste("SELECT * FROM (VALUES",
    "(1,0,103,'gc','at',1,0,0,[0,1],[false,true],77::BIGINT),",
    "(2,0,104,'ct','tt',1,0,0,[0,1],[false,true],77::BIGINT))",
    "v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
    "sample_index,alleles,phase_before,phase_set)")
  baseline <- rduckvep_haplotypes(con, lower_calls, "compatible_overlap", phase_policy = "strict", hgvs = TRUE)
  with_missing <- paste0(lower_calls, " UNION ALL SELECT * REPLACE (1::UBIGINT AS sample_index,",
    "[NULL,NULL]::INTEGER[] AS alleles) FROM (", lower_calls, ") missing_sample")
  combined <- rduckvep_haplotypes(con, with_missing, "compatible_overlap", phase_policy = "strict", hgvs = TRUE)
  known <- combined[vapply(combined$carriers, function(x) any(x$sample_index == 0), logical(1L)), , drop = FALSE]
  row.names(known) <- NULL
  expect_identical(known, baseline)
  with_split <- paste0(lower_calls, " UNION ALL SELECT * REPLACE (1::UBIGINT AS sample_index) FROM (",
    lower_calls, ") split_sample UNION ALL SELECT 3,0,107,'A','C',1,0,1,[0,1],[false,true],78::BIGINT")
  for (direction in c("ASC", "DESC")) {
    calls <- paste0("SELECT * FROM (", with_split, ") s ORDER BY event_index ", direction, ",sample_index ", direction)
    combined <- rduckvep_haplotypes(con, calls, "compatible_overlap", phase_policy = "strict", hgvs = TRUE)
    complete <- vapply(combined$carriers, function(x) any(x$sample_index == 0), logical(1L))
    known <- combined[complete, , drop = FALSE]
    row.names(known) <- NULL
    expect_identical(known, baseline)
    expect_true(all(combined$prediction_status[!complete] == "incomplete_input"))
    expect_true(all(combined$prediction_reason[!complete] == "unresolved_cross_ps_phase"))
  }
  reverse_calls <- paste("SELECT * FROM (VALUES",
    "(2,1,109,'ag','aa',1,1,0,[0,1],[false,true],77::BIGINT),",
    "(1,1,110,'gc','at',1,1,0,[0,1],[false,true],77::BIGINT))",
    "v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
    "sample_index,alleles,phase_before,phase_set)")
  for (name in c("overlap_lower", "overlap_mito")) {
    sequence <- if (name == "overlap_lower") "atggctgaattttaa" else "ATGGCTGAATTTAGG"
    model_tx <- paste0("SELECT * REPLACE (", quote(sequence), "::BLOB AS cds_sequence,",
      if (name == "overlap_mito") 2L else 1L, "::UTINYINT AS codon_table) FROM (", tx, ") t")
    expect_equal(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load(", quote(name), ",",
      "'SELECT i::UINTEGER seq_region FROM range(2) t(i)',", quote(model_tx), ",", quote(exon), ")"))$loaded, TRUE)
    for (strand in c(1L, -1L)) {
      calls <- if (strand == 1L) lower_calls else reverse_calls
      collapsed_calls <- paste0("SELECT * FROM (", calls, ") c WHERE event_index = 1")
      result <- rduckvep_haplotypes(con, calls, name, phase_policy = "strict", hgvs = TRUE)
      collapsed <- rduckvep_haplotypes(con, collapsed_calls, name, phase_policy = "strict", hgvs = TRUE)
      expect_identical(result$cds, collapsed$cds)
      expect_identical(result$protein, collapsed$protein)
      expect_identical(result$protein_differences, collapsed$protein_differences)
      expect_identical(result$hgvsc, "c.4_5delinsAT")
      expect_identical(result$hgvsp, "p.(Ala2Ile)")
      expect_identical(result$contributors[[1L]]$reference, if (strand == 1L) c("gc", "ct") else c("ag", "gc"))
      expect_identical(result$contributors[[1L]]$alternate, if (strand == 1L) c("at", "tt") else c("aa", "at"))
    }
  }
})
