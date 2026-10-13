#!/usr/bin/env Rscript

main <- function() {
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION", "build/release/duckvep.duckdb_extension"), mustWork = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste("LOAD", quote(extension)))

  # Independent literal codon oracle: reference-programme position 3 (TGA -> U)
  # persists only across an intact reference codon. The call itself supplies no
  # evidence that this recoding happens in vivo.
  reference <- "ATGGCTTGAGCATAA"
  reverse_complement <- function(x) {
    bases <- strsplit(x, "", fixed = TRUE)[[1L]]
    paste(rev(chartr("ACGT", "TGCA", bases)), collapse = "")
  }
  model <- function(strand, sequence = reference, position = 3L, suffix = "",
                    code = "_selenocysteine", amino = "U") {
    name <- paste0(if (strand == 1L) "lift_plus" else "lift_minus", suffix)
    n <- nchar(sequence)
    queries <- c(
      paste0("SELECT 0::UINTEGER seq_region, ", n, "::UBIGINT sequence_length, 'chr1'::VARCHAR seq_region_name"),
      paste0("SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region, 1::UBIGINT transcript_start, ",
        n, "::UBIGINT transcript_end, ", strand,
        "::TINYINT strand, 0::UINTEGER gene_index, 0::UBIGINT transcript_flags, ",
        "1::UBIGINT cds_start, ", n, "::UBIGINT cds_end, ", quote(sequence),
        "::BLOB cds_sequence, 1::UTINYINT codon_table"),
      paste0("SELECT 0::UINTEGER transcript_index, 1::UBIGINT exon_start, ", n,
        "::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, ", n,
        "::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase")
    )
    peptide <- paste0("SELECT 0::UINTEGER transcript_index, ", position,
      "::UINTEGER protein_position, ", quote(amino), "::VARCHAR alternate_amino_acid,",
      quote(code), "::VARCHAR edit_code")
    loaded <- DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load(", quote(name), ",",
      paste(quote(queries), collapse = ","), ",peptide_edit_query:=", quote(peptide), ")"))
    stopifnot(isTRUE(loaded$loaded))
    name
  }
  call <- function(position, ref, alt, strand, event, sequence = reference) {
    if (strand == -1L) {
      position <- nchar(sequence) - (position + nchar(ref) - 1L) + 1L
      ref <- reverse_complement(ref)
      alt <- reverse_complement(alt)
    }
    paste0("SELECT ", event, "::UBIGINT record_index, ", event,
      "::UBIGINT event_index, 0::UINTEGER seq_region, ", position,
      "::UBIGINT \"position\", ", quote(ref), "::VARCHAR reference, ",
      quote(alt), "::VARCHAR alternate, 1::UINTEGER alt_index,",
      "0::UINTEGER transcript_index, 0::UINTEGER sample_index,",
      "[0,1]::INTEGER[] alleles, [TRUE]::BOOLEAN[] phase_before, 1::BIGINT phase_set")
  }
  predict <- function(name, calls) {
    DBI::dbGetQuery(con, paste0(
      "SELECT prediction_status,prediction_reason,prediction_reference_protein,",
      "prediction_protein,protein,cds,nmd_prediction,",
      "haplotype_consequences IS NULL AS no_so FROM duckvep_haplotypes(",
      quote(paste(calls, collapse = " UNION ALL ")), ",", quote(name),
      ",phase_policy:='strict')"))
  }
  check <- function(observed, expected, raw, expected_cds,
                    reference_protein = "MAUA*") {
    stopifnot(nrow(observed) == 1L,
      identical(observed$prediction_status, "conditional_assumed_recoding"),
      identical(observed$prediction_reason, "assumed_recoding_programme"),
      identical(observed$prediction_reference_protein, reference_protein),
      identical(observed$prediction_protein, expected),
      identical(observed$protein, raw),
      identical(as.character(observed$cds), expected_cds),
      identical(observed$nmd_prediction, "unknown"), isTRUE(observed$no_so))
  }
  for (strand in c(1L, -1L)) {
    name <- model(strand)
    if (identical(Sys.getenv("DUCKVEP_PADDING_WITNESS"), "1")) {
      for (programme in list(c("_selenocysteine", "U"), c("_stop_codon_rt", "W"))) {
        context <- model(strand, suffix = paste0("_witness_", programme[[2L]]),
          code = programme[[1L]], amino = programme[[2L]])
        for (variant in list(c(10L, "G", "T"), c(7L, "TGAGCA", "TGATCA"), c(9L, "AGCA", "A"))) {
          result <- predict(context, call(as.integer(variant[[1L]]), variant[[2L]],
            variant[[3L]], strand, 30L))
          message(paste(strand, programme[[1L]], paste(variant, collapse = "/"),
            result$prediction_status, result$prediction_protein, sep = " | "))
        }
      }
      next
    }
    fingerprint <- DBI::dbGetQuery(con, paste0("SELECT _duckvep_model_fingerprint(",
      quote(name), ")"))[[1L]]
    # Insert a complete CCT codon before Sec; codon 3 becomes codon 4.
    check(predict(name, call(6L, "T", "TCCT", strand, 1L)),
      "MAPUA*", "MAP*", "ATGGCTCCTTGAGCATAA")
    check(predict(name, call(6L, "TTGA", "TCCTTGA", strand, 14L)),
      "MAPUA*", "MAP*", "ATGGCTCCTTGAGCATAA")
    # Remove codon 2, preserving the original Sec at the new codon 2.
    check(predict(name, call(3L, "GGCT", "G", strand, 2L)),
      "MUA*", "M*", "ATGTGAGCATAA")
    check(predict(name, call(9L, "A", "ACCT", strand, 11L)),
      "MAUPA*", "MA*", "ATGGCTTGACCTGCATAA")
    # Genomic left anchoring retains CDS 9 on plus and CDS 12 on minus.
    check(predict(name, call(9L, "AGCA", "A", strand, 12L)),
      if (strand == 1L) "MAU*" else "MA*", "MA*", "ATGGCTTGATAA")
    check(predict(name, call(4L, "GCTT", "T", strand, 13L)),
      "MUA*", "M*", "ATGTGAGCATAA")
    for (programme in list(c("_selenocysteine", "U", "MAUA*"),
                           c("_stop_codon_rt", "W", "MAWA*"))) {
      context <- model(strand, suffix = paste0("_context_", programme[[2L]]),
        code = programme[[1L]], amino = programme[[2L]])
      minimal <- predict(context, call(10L, "G", "T", strand, 21L))
      padded <- predict(context, call(7L, "TGAGCA", "TGATCA", strand, 22L))
      check(minimal, sub("A\\*$", "S*", programme[[3L]]), "MA*", "ATGGCTTGATCATAA",
        programme[[3L]])
      check(padded, minimal$prediction_protein, "MA*", "ATGGCTTGATCATAA",
        programme[[3L]])
    }
    # A stop before the intact original Sec prevents it from being reached.
    check(predict(name, call(4L, "GCT", "TAA", strand, 3L)),
      "M*", "M*", "ATGTAATGAGCATAA")
    # Remove the programme codon itself: an identical later TGA is not Sec.
    removed <- predict(name, call(6L, "TTGA", "T", strand, 4L))
    stopifnot(nrow(removed) == 1L, identical(removed$prediction_status,
      "conditional_assumed_recoding"), identical(removed$prediction_protein, "MAA*"),
      identical(as.character(removed$cds), "ATGGCTGCATAA"))
    # A changed TGA cannot inherit U by original or shifted ordinal.
    check(predict(name, call(8L, "G", "C", strand, 5L)),
      "MASA*", "MASA*", "ATGGCTTCAGCATAA")
    # Two independent source edits restore phase. The last Sec codon is
    # original-reference bases, not an equal-length aggregate descriptor.
    restored <- "ATGGCTGCCTGAGCATAA"
    restored_model <- model(strand, restored, 4L, "_restored")
    restored_calls <- lapply(list(call(3L, "G", "GC", strand, 6L, restored),
      call(8L, "CC", "C", strand, 7L, restored)),
      sub, pattern = "[TRUE]::BOOLEAN[]", replacement = "[TRUE,TRUE]::BOOLEAN[]",
      fixed = TRUE)
    check(predict(restored_model, restored_calls), "MRCUA*", "MRC*",
      "ATGCGCTGCTGAGCATAA", "MAAUA*")
    readthrough <- model(strand, suffix = "_rt", code = "_stop_codon_rt", amino = "W")
    check(predict(readthrough, call(6L, "T", "TCCT", strand, 8L)),
      "MAPWA*", "MAP*", "ATGGCTCCTTGAGCATAA", "MAWA*")
    initial <- "GTGGCTTAA"
    start_model <- model(strand, initial, 1L, "_start", "initial_met", "M")
    check(predict(start_model, call(3L, "G", "GCCT", strand, 9L, initial)),
      "MPA*", "VPA*", "GTGCCTGCTTAA", "MA*")
    check(predict(start_model, call(1L, "G", "A", strand, 10L, initial)),
      "MA*", "MA*", "ATGGCTTAA", "MA*")
    stopifnot(identical(fingerprint, DBI::dbGetQuery(con,
      paste0("SELECT _duckvep_model_fingerprint(", quote(name), ")"))[[1L]]))
  }
  message("conditional indel lift: OK")
}
main()
