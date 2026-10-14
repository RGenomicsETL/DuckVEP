#!/usr/bin/env Rscript

main <- function() {
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION"), mustWork = TRUE)
  con <- Rduckvep::rduckvep_connect(extension_path = extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  reference <- "ATGGCTGCTGAAGGTTAA"
  translate <- function(dna) {
    code <- c(ATG = "M", GCT = "A", GCC = "A", ACT = "T", ACC = "T",
      GAA = "E", GGT = "G", CGT = "R", TAA = "*", TGA = "*", GTT = "V",
      CAA = "Q", CTG = "L", TCT = "S", AGT = "S", TGT = "C")
    bases <- strsplit(dna, "", fixed = TRUE)[[1L]]
    codons <- vapply(seq.int(1L, length(bases) - 2L, by = 3L),
      function(i) paste0(bases[i:(i + 2L)], collapse = ""), character(1L))
    amino <- unname(code[codons])
    stopifnot(!anyNA(amino))
    end <- match("*", amino)
    paste0(amino[seq_len(if (is.na(end)) length(amino) else end)], collapse = "")
  }
  stopifnot(identical(translate(reference), "MAAEG*"))
  for (strand in c(1L, -1L)) {
    model <- paste0("reached_", if (strand == 1L) "plus" else "minus")
    tx <- sprintf(paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
      "100::UBIGINT transcript_start,117::UBIGINT transcript_end,%d::TINYINT strand,",
      "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
      "117::UBIGINT cds_end,%s::BLOB cds_sequence,1::UTINYINT codon_table"),
      strand, q(reference))
    ex <- paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,",
      "117::UBIGINT exon_end,1::UBIGINT exon_cdna_start,18::UBIGINT exon_cdna_end,",
      "0::TINYINT phase,0::TINYINT end_phase")
    stopifnot(Rduckvep::rduckvep_load_model(con, model,
      "SELECT 0::UINTEGER seq_region", tx, ex))
    run <- function(changes, expected_dna, expected_terms, expected_roles = NULL) {
      complement <- c(A = "T", T = "A", G = "C", C = "G")
      rows <- vapply(seq_len(nrow(changes)), function(i) {
        position <- changes$position[i]
        ref <- changes$ref[i]
        alt <- changes$alt[i]
        if (strand == -1L) {
          position <- 217L - position
          ref <- unname(complement[ref])
          alt <- unname(complement[alt])
        }
        sprintf("(%d,%d,%s,%s)", changes$id[i], position, q(ref), q(alt))
      }, character(1L))
      sql <- paste0("SELECT event_index,0 seq_region,position,reference,alternate,",
        "1 alt_index,0 transcript_index,0 sample_index,[1,1] alleles,",
        "[true,true] phase_before,NULL::BIGINT phase_set FROM (VALUES ",
        paste(rows, collapse = ","), ") v(event_index,position,reference,alternate)")
      observed <- Rduckvep::rduckvep_haplotypes(con, sql, model)
      stopifnot(nrow(observed) == 1L,
        identical(observed$prediction_status, "predicted"),
        identical(observed$cds, expected_dna),
        identical(observed$protein, translate(expected_dna)),
        setequal(observed$haplotype_consequences[[1L]], expected_terms),
        length(observed$haplotype_consequences[[1L]]) == length(expected_terms))
      if (!is.null(expected_roles))
        stopifnot(identical(observed$contributor_provenance[[1L]]$role[
          order(as.integer(observed$contributor_provenance[[1L]]$event_index))], expected_roles))
      observed
    }
    substitutions <- function(position, ref, alt, id = seq_along(position))
      data.frame(id = id, position = position, ref = ref, alt = alt)
    dna <- function(position, alt) {
      bases <- strsplit(reference, "", fixed = TRUE)[[1L]]
      bases[position - 99L] <- alt
      paste0(bases, collapse = "")
    }
    three <- substitutions(c(103L, 109L, 112L), c("G", "G", "G"), c("A", "T", "C"))
    witness <- run(three, dna(c(103L, 109L, 112L), c("A", "T", "C")),
      c("stop_gained", "missense_variant"), c("applied", "applied", "post_stop"))
    stopifnot(identical(as.integer(witness$nmd_stop_position), 12L),
      identical(witness$nmd_prediction, "escape"))
    run(three[2:3, ], dna(c(109L, 112L), c("T", "C")),
      "stop_gained", c("applied", "post_stop"))
    run(substitutions(c(105L, 109L), c("T", "G"), c("C", "T")),
      dna(c(105L, 109L), c("C", "T")), c("stop_gained", "synonymous_variant"))
    run(substitutions(c(103L, 105L, 109L), c("G", "T", "G"), c("A", "C", "T")),
      dna(c(103L, 105L, 109L), c("A", "C", "T")),
      c("stop_gained", "missense_variant"))
    run(three[c(3L, 1L, 2L), ], dna(c(103L, 109L, 112L), c("A", "T", "C")),
      c("stop_gained", "missense_variant"))
    run(substitutions(c(103L, 115L), c("G", "T"), c("A", "C")),
      dna(c(103L, 115L), c("A", "C")), c("stop_lost", "missense_variant"))
    lost <- run(substitutions(100L, "A", "C"), dna(100L, "C"), "start_lost")
    stopifnot(identical(lost$nmd_prediction, "unknown"))

    # Two individually missense substitutions in one codon yield final Ser.
    serine <- "ATGTCTGCTGAAGGTTAA"
    ser_model <- paste0(model, "_serine")
    ser_tx <- sub(reference, serine, tx, fixed = TRUE)
    stopifnot(Rduckvep::rduckvep_load_model(con, ser_model,
      "SELECT 0::UINTEGER seq_region", ser_tx, ex))
    ser_changes <- substitutions(c(103L, 104L), c("T", "C"), c("A", "G"))
    ser_sql <- paste0("SELECT event_index,0 seq_region,position,reference,alternate,",
      "1 alt_index,0 transcript_index,0 sample_index,[1,1] alleles,",
      "[true,true] phase_before,NULL::BIGINT phase_set FROM (VALUES ",
      paste(vapply(seq_len(nrow(ser_changes)), function(i) {
        position <- ser_changes$position[i]
        ref <- ser_changes$ref[i]
        alt <- ser_changes$alt[i]
        if (strand == -1L) {
          position <- 217L - position
          ref <- unname(c(A = "T", T = "A", G = "C", C = "G")[ref])
          alt <- unname(c(A = "T", T = "A", G = "C", C = "G")[alt])
        }
        sprintf("(%d,%d,%s,%s)", i, position, q(ref), q(alt))
      }, character(1L)), collapse = ","), ") v(event_index,position,reference,alternate)")
    ser <- Rduckvep::rduckvep_haplotypes(con, ser_sql, ser_model)
    stopifnot(nrow(ser) == 1L, identical(ser$cds, "ATGAGTGCTGAAGGTTAA"),
      identical(translate(ser$cds), "MSAEG*"),
      identical(ser$haplotype_consequences[[1L]], "synonymous_variant"))
  }
}

main()
