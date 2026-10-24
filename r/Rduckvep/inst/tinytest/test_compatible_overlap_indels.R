library(tinytest)
library(DBI)

local({
  extension <- Sys.getenv("DUCKVEP_EXTENSION")
  con <- rduckvep_connect(extension_path = if (nzchar(extension)) extension else NULL)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  reference <- "ATGGCTGAATAA"
  bases <- strsplit(reference, "", fixed = TRUE)[[1L]]
  # Reference-indexed oracle over minimal CDS islands: a base is kept, substituted or absent, and an
  # interbase site carries at most one inserted string. Identical assertions coincide.
  union_cds <- function(start, ref, alt) {
    value <- bases; absent <- logical(length(bases)); inserted <- character(length(bases) + 1L)
    for (i in seq_along(start)) {
      at <- start[i] + seq_len(nchar(ref[i])) - 1L
      stopifnot(identical(paste(bases[at], collapse = ""), ref[i]))
      if (!nchar(ref[i])) {
        stopifnot(inserted[start[i]] %in% c("", alt[i])); inserted[start[i]] <- alt[i]
      } else if (!nchar(alt[i])) absent[at] <- TRUE
      else value[at] <- strsplit(alt[i], "", fixed = TRUE)[[1L]]
    }
    paste0(paste0(inserted[seq_along(bases)], ifelse(absent, "", value), collapse = ""), inserted[length(bases) + 1L])
  }
  complement <- c(A = "T", C = "G", G = "C", T = "A")
  rc <- function(x) paste(rev(unname(complement[strsplit(x, "", fixed = TRUE)[[1L]]])), collapse = "")
  quote <- function(x) as.character(dbQuoteString(con, x))
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,111::UBIGINT transcript_end,(1-2*i)::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "111::UBIGINT cds_end,'ATGGCTGAATAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2) t(i)")
  exon <- paste("SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM range(2) t(i)")
  expect_equal(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('overlap_indels',",
    "'SELECT i::UINTEGER seq_region FROM range(2) t(i)',", quote(tx), ",", quote(exon), ")"))$loaded, TRUE)
  # Anchored CDS-orientation records, their minimal islands, and Mutalyzer 3.1.1 joint names.
  witnesses <- list(
    duplicate_deletion = list(pos = c(3L, 3L), ref = c("GGCT", "GGCT"), alt = c("G", "G"),
      start = c(4L, 4L), island_ref = c("GCT", "GCT"), island_alt = c("", ""),
      hgvsc = "c.5_7del", hgvsp = "p.(Ala2del)"),
    nested_deletion = list(pos = c(3L, 4L), ref = c("GGCTGAA", "GCT"), alt = c("G", "G"),
      start = c(4L, 5L), island_ref = c("GCTGAA", "CT"), island_alt = c("", ""),
      hgvsc = "c.4_9del", hgvsp = "p.(Ala2_Glu3del)"),
    overlapping_deletion = list(pos = c(4L, 5L), ref = c("GCT", "CTG"), alt = c("G", "C"),
      start = c(5L, 6L), island_ref = c("CT", "TG"), island_alt = c("", ""),
      hgvsc = "c.5_7del", hgvsp = "p.(Ala2del)"),
    duplicate_insertion = list(pos = c(6L, 6L), ref = c("T", "T"), alt = c("TCCC", "TCCC"),
      start = c(7L, 7L), island_ref = c("", ""), island_alt = c("CCC", "CCC"),
      hgvsc = "c.6_7insCCC", hgvsp = "p.(Ala2_Glu3insPro)"),
    insertion_before_substitution = list(pos = c(6L, 7L), ref = c("T", "G"), alt = c("TCCC", "A"),
      start = c(7L, 7L), island_ref = c("", "G"), island_alt = c("CCC", "A"),
      hgvsc = "c.7delinsCCCA", hgvsp = "p.(Glu3delinsProLys)"),
    # Records apply at their own coordinates, also inside a repeat where an equivalent placement exists.
    nested_repeat_deletion = list(pos = c(3L, 4L), ref = c("GGCT", "GCT"), alt = c("G", "G"),
      start = c(4L, 5L), island_ref = c("GCT", "CT"), island_alt = c("", ""),
      hgvsc = "c.5_7del", hgvsp = "p.(Ala2del)"),
    insertion_before_substitution_in_repeat = list(pos = c(6L, 7L), ref = c("T", "G"), alt = c("TGGG", "A"),
      start = c(7L, 7L), island_ref = c("", "G"), island_alt = c("GGG", "A"),
      hgvsc = "c.7_8insGGA", hgvsp = "p.(Glu3delinsGlyLys)"),
    removed_substituted_base = list(pos = c(4L, 5L), ref = c("GCT", "C"), alt = c("G", "A"),
      status = "edit_conflict"),
    # The deletion record removes c.8 and the substitution changes c.8; c.9del would be an equivalent
    # deletion, but it is not what the record says.
    deleted_base_substituted_in_repeat = list(pos = c(7L, 8L), ref = c("GA", "A"), alt = c("G", "T"),
      status = "edit_conflict"))
  for (strand in c(1L, -1L)) for (reverse in c(FALSE, TRUE)) for (name in names(witnesses)) {
    w <- witnesses[[name]]
    genomic_position <- if (strand == 1L) 99L + w$pos else 112L - w$pos - nchar(w$ref) + 1L
    genomic_ref <- if (strand == 1L) w$ref else vapply(w$ref, rc, "")
    genomic_alt <- if (strand == 1L) w$alt else vapply(w$alt, rc, "")
    rows <- sprintf("(%d,%d,%d,%s,%s,1,%d,0,[0,1],[false,true],77::BIGINT)", seq_along(w$pos),
      as.integer(strand == -1L), genomic_position, quote(genomic_ref), quote(genomic_alt), as.integer(strand == -1L))
    if (reverse) rows <- rev(rows)
    calls <- paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") v(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,alleles,phase_before,phase_set)")
    result <- rduckvep_haplotypes(con, calls, "overlap_indels", phase_policy = "strict", hgvs = TRUE)
    expect_equal(nrow(result), 1L, info = name)
    sources <- result$contributors[[1L]]
    sources <- sources[order(sources$event_index), ]
    expect_equal(sources$position, genomic_position, info = name)
    expect_identical(sources$reference, unname(genomic_ref), info = name)
    expect_identical(sources$alternate, unname(genomic_alt), info = name)
    roles <- result$contributor_provenance[[1L]]$role
    if (is.null(w$status)) {
      expect_identical(result$prediction_status, "predicted", info = name)
      expect_identical(result$cds, union_cds(w$start, w$island_ref, w$island_alt), info = name)
      expect_identical(result$hgvsc, w$hgvsc, info = name)
      expect_identical(result$hgvsp, w$hgvsp, info = name)
      expect_identical(roles, c("applied", "applied"), info = name)
      blocks <- result$coding_blocks[[1L]]
      expect_equal(nrow(blocks), 1L, info = name)
      expect_equal(sort(blocks$event_indices[[1L]]), c(1, 2), info = name)
      edits <- result$normalized_edits[[1L]]
      edits <- edits[order(edits$event_index), ]
      expect_equal(edits$cds_start, w$start, info = name)
      # Normalized islands keep the source's genomic orientation.
      expect_identical(edits$reference, if (strand == 1L) w$island_ref else unname(vapply(w$island_ref, rc, "")), info = name)
      expect_equal(edits$block_index, c(0, 0), info = name)
    } else {
      expect_identical(result$prediction_status, w$status, info = name)
      expect_true(is.na(result$cds) && is.na(result$hgvsc) && is.na(result$hgvsp), info = name)
      expect_identical(roles, c("unapplied", "unapplied"), info = name)
    }
  }
})
