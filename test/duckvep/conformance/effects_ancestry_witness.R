#!/usr/bin/env Rscript

main <- function() {
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION"), mustWork = TRUE)
  con <- Rduckvep::rduckvep_connect(extension_path = extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  run <- function(reference, events, strand = 1L) {
    n <- nchar(reference)
    model <- paste0("effects_", reference, "_", strand)
    if (strand == -1L) {
      events$pos <- n - events$pos - nchar(events$ref) + 2L
      revcomp <- function(x) chartr("ACGT", "TGCA",
        paste(rev(strsplit(x, "", fixed = TRUE)[[1L]]), collapse = ""))
      events$ref <- vapply(events$ref, revcomp, character(1L))
      events$alt <- vapply(events$alt, revcomp, character(1L))
    }
    tx <- sprintf(paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
      "100::UBIGINT transcript_start,%d::UBIGINT transcript_end,%d::TINYINT strand,",
      "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
      "%d::UBIGINT cds_end,%s::BLOB cds_sequence,1::UTINYINT codon_table"),
      99L + n, strand, 99L + n, q(reference))
    ex <- sprintf(paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,",
      "%d::UBIGINT exon_end,1::UBIGINT exon_cdna_start,%d::UBIGINT exon_cdna_end,",
      "0::TINYINT phase,0::TINYINT end_phase"), 99L + n, n)
    if (!model %in% loaded) {
      stopifnot(Rduckvep::rduckvep_load_model(con, model, "SELECT 0::UINTEGER seq_region", tx, ex))
      loaded <<- c(loaded, model)
    }
    rows <- vapply(seq_len(nrow(events)), function(i)
      sprintf("(%d,%d,%s,%s)", events$id[i], events$pos[i] + 99L,
        q(events$ref[i]), q(events$alt[i])), character(1L))
    sql <- paste0("SELECT event_index,0 seq_region,position,reference,alternate,",
      "1 alt_index,0 transcript_index,0 sample_index,[1,1] alleles,",
      "[true,true] phase_before,NULL::BIGINT phase_set FROM (VALUES ",
      paste(rows, collapse = ","), ") v(event_index,position,reference,alternate)")
    result <- Rduckvep::rduckvep_haplotypes(con, sql, model)
    stopifnot(nrow(result) == 1L, identical(result$prediction_status, "predicted"))
    result
  }
  loaded <- character()
  cases <- list(
    opaque = data.frame(id = 1L, pos = 4L, ref = "GCT", alt = "TTTCCC"),
    opaque_alias = data.frame(id = 1:2, pos = 4L, ref = "GCT", alt = "TTTCCC"),
    joined = data.frame(id = 1:3, pos = c(4L, 6L, 11L),
      ref = c("GCT", "TGAA", "TT"), alt = c("TTTCCC", "T", "AA")),
    joined_alias = data.frame(id = 1:4, pos = c(11L, 4L, 6L, 4L),
      ref = c("TT", "GCT", "TGAA", "GCT"), alt = c("AA", "TTTCCC", "T", "TTTCCC")),
    stop = data.frame(id = 1L, pos = 11L, ref = "TT", alt = "AA"),
    balanced_descriptor = data.frame(id = 1:2, pos = c(4L, 11L),
      ref = c("GCTGAA", "TT"), alt = c("TTTCCC", "AA")))
  observed <- lapply(cases, function(events) run("ATGGCTGAATTTTAA", events))
  effects <- function(x) sort(x$haplotype_consequences[[1L]])
  expect <- function(x, cds, protein, terms) {
    stopifnot(identical(x$prediction_status, "predicted"),
      identical(x$cds, cds), identical(x$protein, protein),
      identical(effects(x), sort(terms)))
  }
  expect(observed$opaque, "ATGTTTCCCGAATTTTAA", "MFPEF*", "protein_altering_variant")
  stopifnot(identical(effects(observed$opaque), effects(observed$opaque_alias)))
  expect(observed$joined, "ATGTTTCCCTAATAA", "MFP*", "stop_gained")
  stopifnot(identical(effects(observed$joined), effects(observed$joined_alias)))
  expect(observed$stop, "ATGGCTGAATAATAA", "MAE*", "stop_gained")
  expect(observed$balanced_descriptor, "ATGTTTCCCTAATAA", "MFP*",
    c("stop_gained", "missense_variant"))
  message("opaque, aliases, joined and balanced source assertions passed")
  extended <- data.frame(id = 1:4, pos = c(4L, 7L, 9L, 14L),
    ref = c("G", "GCT", "TGAA", "TT"), alt = c("A", "TTTCCC", "T", "AA"))
  result <- run("ATGGCTGCTGAATTTTAA", extended)
  expect(result, "ATGACTTTTCCCTAATAA", "MTFP*",
    c("stop_gained", "missense_variant"))
  shuffled <- extended[c(4L, 2L, 1L, 3L, 2L), ]
  shuffled$id <- seq_len(nrow(shuffled))
  stopifnot(identical(effects(result), effects(run("ATGGCTGCTGAATTTTAA", shuffled))))
  message("independent upstream substitution and disjoint order assertions passed")
  deletion_stop <- data.frame(id = 1:4, pos = c(3L, 3L, 7L, 10L),
    ref = c("GGCT", "GGCT", "GAA", "T"), alt = c("G", "G", "TAA", "C"))
  opaque_stop <- data.frame(id = 1:2, pos = c(4L, 10L),
    ref = c("GCTG", "T"), alt = c("T", "C"))
  # A deleted original codon precedes the reached stop at original bases 7..9.
  # The opaque assertion supplies the same molecule without that codon ancestry.
  ref_bases <- strsplit("ATGGCTGAATTTTAA", "", fixed = TRUE)[[1L]]
  retained <- c(1:3, 7:15)
  molecular <- ref_bases[retained]
  molecular[match(c(7L, 10L), retained)] <- c("T", "C")
  stopifnot(identical(paste(molecular, collapse = ""), "ATGTAACTTTAA"),
    identical(retained[4:6], 7:9), identical(paste(ref_bases[4:6], collapse = ""), "GCT"))
  for (strand in c(1L, -1L)) {
    for (events in list(deletion_stop, deletion_stop[c(4L, 2L, 3L, 1L), ],
      deletion_stop[c(1L, 3L, 4L), ])) {
      result <- run("ATGGCTGAATTTTAA", events, strand)
      expect(result, "ATGTAACTTTAA", "M*", c("stop_gained", "inframe_deletion"))
      stopifnot(identical(result$nmd_prediction, "escape"), result$nmd_stop_position == 6L)
    }
    result <- run("ATGGCTGAATTTTAA", opaque_stop, strand)
    expect(result, "ATGTAACTTTAA", "M*", "stop_gained")
    expect(run("ATGGCTGAATTTTAA", deletion_stop[1L, , drop = FALSE], strand),
      "ATGGAATTTTAA", "MEF*", "inframe_deletion")
    expect(run("ATGGCTGAATTTTAA", deletion_stop[3:4, ], strand),
      "ATGGCTTAACTTTAA", "MA*", "stop_gained")
    post_stop <- data.frame(id = 1:2, pos = c(4L, 9L),
      ref = c("GCT", "ATTT"), alt = c("TAA", "A"))
    expect(run("ATGGCTGAATTTTAA", post_stop, strand),
      "ATGTAAGAATAA", "M*", "stop_gained")
  }
  message("source-distinct deletion before immediate stop and opaque/post-stop controls passed on both strands")
  frame <- data.frame(id = 1:3, pos = c(3L, 12L, 18L),
    ref = c("G", "G", "AG"), alt = c("GT", "GCCC", "A"))
  control <- run("ATGGCAGAACCGTTTCAAGGATCCTAA", frame[2L, , drop = FALSE])
  message("standalone +3: ", control$protein, " / ", paste(effects(control), collapse = ","))
  for (case in list(frame, frame[c(1L, 3L), ],
    transform(frame, alt = c("GT", "GTAA", "A")))) {
    result <- run("ATGGCAGAACCGTTTCAAGGATCCTAA", case)
    stopifnot(identical(effects(result), "protein_altering_variant"))
    message("frame restored: ", result$cds, " / ", result$protein)
  }
  early <- run("ATGGCTGAACCGTTTCAAGGATCCTAA", frame)
  expect(early, "ATGTGCTGAACCGCCCTTTCAAGATCCTAA", "MC*",
    c("stop_gained", "frameshift_variant"))
  late <- rbind(frame, data.frame(id = 4:5, pos = 23:24,
    ref = c("C", "C"), alt = c("A", "A")))
  result <- run("ATGGCAGAACCGTTTCAAGGATCCTAA", late)
  expect(result, "ATGTGCAGAACCGCCCTTTCAAGATAATAA", "MCRTALSR*", "stop_gained")
  message("before/after frame-restoration stop assertions passed")
}
main()
