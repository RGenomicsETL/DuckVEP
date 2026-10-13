#!/usr/bin/env Rscript

# Study of original-call copy identifiability and raw/decoded admission.
main <- function() {
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION"), mustWork = TRUE)
  con <- Rduckvep::rduckvep_connect(extension_path = extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  tx <- paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,117::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "117::UBIGINT cds_end,'ATGGCTGCTGAAGGTTAA'::BLOB cds_sequence,1::UTINYINT codon_table")
  ex <- paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,",
    "117::UBIGINT exon_end,1::UBIGINT exon_cdna_start,18::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase")
  stopifnot(Rduckvep::rduckvep_load_model(con, "completion_study",
    "SELECT 0::UINTEGER seq_region", tx, ex))
  genotypes <- list(homozygous = c("1/1", "1/1"),
    one_heterozygote = c("0/1", "1/1"), two_heterozygotes = c("0/1", "0/1"))
  reference <- strsplit("ATGGCTGCTGAAGGTTAA", "", fixed = TRUE)[[1L]]
  observed <- lapply(names(genotypes), function(name) {
    gt <- genotypes[[name]]
    alleles <- vapply(strsplit(gt, "/", fixed = TRUE), as.integer, integer(2L))
    heterozygous <- which(alleles[1L, ] != alleles[2L, ])
    flips <- if (length(heterozygous) > 0L) expand.grid(rep(list(0:1), length(heterozygous))) else data.frame(flip = 0L)
    pairs <- vapply(seq_len(nrow(flips)), function(i) {
      copies <- alleles
      for (j in seq_along(heterozygous)) {
        if (flips[i, j] == 1L) copies[, heterozygous[j]] <- rev(copies[, heterozygous[j]])
      }
      sequences <- apply(copies, 1L, function(copy) {
        dna <- reference
        dna[c(4L, 13L)[copy == 1L]] <- c("A", "C")[copy == 1L]
        paste0(dna, collapse = "")
      })
      paste(sort(sequences), collapse = "/")
    }, character(1L))
    raw_sql <- sprintf(paste("SELECT event_index,0 seq_region,position,reference,alternates,",
      "0 transcript_index,0 sample_index,gt FROM (VALUES",
      "(1,103,'G',['A'],%s),(2,112,'G',['C'],%s))",
      "v(event_index,position,reference,alternates,gt)"), quote(gt[1L]), quote(gt[2L]))
    decoded_sql <- sprintf(paste("SELECT event_index,0 seq_region,position,reference,alternate,",
      "1 alt_index,0 transcript_index,0 sample_index,alleles,",
      "[false,false] phase_before,NULL::BIGINT phase_set FROM (VALUES",
      "(1,103,'G','A',[%s]),(2,112,'G','C',[%s]))",
      "v(event_index,position,reference,alternate,alleles)"),
      paste(alleles[, 1L], collapse = ","), paste(alleles[, 2L], collapse = ","))
    raw <- Rduckvep::rduckvep_haplotypes(con, raw_sql, "completion_study",
      phase_policy = "vep_compat", input_mode = "source_records")
    decoded <- Rduckvep::rduckvep_haplotypes(con, decoded_sql, "completion_study")
    arrangements <- tryCatch(DBI::dbGetQuery(con,
      paste0("SELECT * FROM duckvep_haplotype_arrangements(", quote(decoded_sql),
        ", 'completion_study')")), error = identity)
    if (length(heterozygous) < ncol(alleles)) {
      stopifnot(inherits(arrangements, "error"),
        grepl("requires non-NULL diploid biallelic heterozygous alt_events rows",
          conditionMessage(arrangements), fixed = TRUE))
      count <- NA_integer_
      hypothetical <- NA
    } else {
      stopifnot(!inherits(arrangements, "error"))
      count <- length(unique(arrangements$hypothesis_id))
      hypothetical <- any(arrangements$prediction_semantics == "hypothetical_assignment")
      stopifnot(count == length(unique(pairs)))
    }
    stopifnot(nrow(raw) > 0L, nrow(decoded) > 0L)
    if (length(heterozygous) <= 1L) {
      stopifnot(all(raw$prediction_status == "predicted"),
        identical(paste(sort(rep(raw$cds, raw$carrier_count)), collapse = "/"), unique(pairs)))
    }
    data.frame(case = name, independent_pairs = length(unique(pairs)),
      native_arrangements = count,
      raw_status = paste(sort(unique(raw$prediction_status)), collapse = ";"),
      raw_reason = paste(sort(unique(raw$prediction_reason)), collapse = ";"),
      decoded_status = paste(sort(unique(decoded$prediction_status)), collapse = ";"),
      hypothetical = hypothetical)
  })
  result <- do.call(rbind, observed)
  stopifnot(identical(result$independent_pairs, c(1L, 1L, 2L)),
    identical(result$raw_status, c("predicted", "predicted", "unsupported_context")),
    identical(result$decoded_status, c("predicted", "predicted", "incomplete_input")),
    identical(result$hypothetical, c(NA, NA, TRUE)))
  write.csv(result, stdout(), row.names = FALSE)

  effects_sql <- paste("SELECT event_index,0 seq_region,position,'G' reference,alternate,",
    "1 alt_index,0 transcript_index,0 sample_index,[1,1] alleles,",
    "[true,true] phase_before,NULL::BIGINT phase_set FROM (VALUES (1,103,'A'),(2,109,'T'),(3,112,'C'))",
    "v(event_index,position,alternate)")
  effects <- Rduckvep::rduckvep_haplotypes(con, effects_sql, "completion_study")
  alternate <- reference
  alternate[c(4L, 10L, 13L)] <- c("A", "T", "C")
  # Standard-code codons before and including the first stop: ATG, ACT, GCT, TAA.
  stopifnot(nrow(effects) == 1L,
    identical(effects$cds, paste0(alternate, collapse = "")),
    identical(effects$protein, "MTA*"),
    identical(effects$haplotype_consequences[[1L]], "stop_gained"),
    identical(effects$contributor_provenance[[1L]]$role, c("applied", "applied", "post_stop")))
  message("Subeffect witness: MAAEG* -> MTA*; current terms=stop_gained; ",
    "roles=applied/applied/post_stop. Ala2Thr is upstream of termination; ",
    "the Gly5Arg DNA change is downstream and is not a translated subeffect.")

  delins_sql <- paste("SELECT event_index,0 seq_region,103 AS position,'GCTGCT' reference,",
    "'AAA' alternate,1 alt_index,0 transcript_index,0 sample_index,[1,1] alleles,",
    "[true,true] phase_before,NULL::BIGINT phase_set FROM range(1,%d) t(event_index)")
  single <- Rduckvep::rduckvep_haplotypes(con, sprintf(delins_sql, 2L), "completion_study")
  duplicate <- Rduckvep::rduckvep_haplotypes(con, sprintf(delins_sql, 3L), "completion_study")
  joint <- paste0(c(reference[1:3], rep("A", 3L), reference[10:18]), collapse = "")
  stopifnot(nrow(single) == 1L, nrow(duplicate) == 1L,
    identical(single$cds, joint), identical(single$protein, "MKEG*"),
    identical(single$prediction_status, "predicted"),
    identical(duplicate$prediction_status, "predicted"),
    identical(duplicate$cds, single$cds), identical(duplicate$protein, single$protein),
    identical(duplicate$nominal_length_diff, single$nominal_length_diff),
    identical(as.integer(duplicate$contributor_provenance[[1L]]$event_index), 1:2))
  message("Delins witness: c.4_9delinsAAA gives ", joint,
    " / MKEG*; duplicate identical assertions apply once and retain both sources.")
}

main()
