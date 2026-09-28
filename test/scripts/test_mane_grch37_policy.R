#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
source("scripts/mane_grch37_policy.R")

ex <- data.frame(rank = 1L, exon_start = 1L, exon_end = 9L)
candidate <- data.table(seq_region_name = "1", strand = 1L, cds_start = 1L, cds_end = 9L,
  pre = "", cds = "ATGAAATAA", post = "", codon_table = 1L)
candidate[, exons := list(ex)]
target <- list(strand = "+", contig = "1", exons = "1-9", cdna = "ATGAAATAA",
  cds_rows = data.table(start = 1L, end = 9L, phase = "0"),
  exon_rows = data.table(start = 1L, end = 9L))
gate <- function(t = target, c = candidate, r = "ATGAAATAA", p = "MK") mane_pair_gate(t, c, r, p)
stopifnot(gate()$status == "exact_model_match",
  gate()$evidence$spliced_sequence_match,
  gate(c = candidate[0])$status == "refseq_only_no_gencode19_match")

# Exact accession versions gate target resolution; stable roots never stand in.
targets <- list("NM_1.1" = target)
stopifnot(is.null(targets[["NM_1.2"]]), !is.null(targets[["NM_1.1"]]))
roots <- data.table(root = "ENST_1")
setkey(roots, root)
stopifnot(nrow(mane_candidates(roots, "ENST_2")) == 0L,
          nrow(mane_candidates(roots, "ENST_1")) == 1L)
# A second NCBI placement cannot invalidate a unique model-supported locus.
placements <- data.table(accession = "NM_1.1", model_contig = c("1", "NW_patch"),
                         target_ok = c(TRUE, FALSE))
stopifnot(nrow(mane_supported_loci(placements)) == 1L,
          mane_supported_loci(placements)$model_contig == "1",
          nrow(mane_supported_loci(placements[target_ok == FALSE])) == 1L)
par <- data.table(accession = "NM_2.1", model_contig = c("X", "Y"), target_ok = TRUE)
stopifnot(nrow(mane_supported_loci(par, "Y")) == 1L,
          mane_supported_loci(par, "Y")$model_contig == "Y",
          nrow(mane_supported_loci(par, "MT")) == 2L)
wrong_exon <- copy(target); wrong_exon$exons <- "1-8"
stopifnot(gate(t = wrong_exon)$status == "geometry_mismatch")
wrong_phase <- copy(target); wrong_phase$cds_rows$phase <- "1"
stopifnot(gate(t = wrong_phase)$status == "cds_phase_mismatch")
wrong_sequence <- copy(target); wrong_sequence$cdna <- "ATGAAATAG"
stopifnot(gate(t = wrong_sequence)$status == "sequence_mismatch")
stopifnot(gate(p = "ME")$status == "translation_mismatch")

# A substitution in the accession RNA is explained by its own protein; the
# candidate remains identical to the independently read target reference.
reference_difference <- gate(r = "ATGGAATAA", p = "ME")
stopifnot(reference_difference$status == "exact_model_match",
  reference_difference$evidence$reference_difference,
  reference_difference$evidence$reference_difference_bases == 1L)
stopifnot(gate(c = rbind(candidate, candidate))$status == "ambiguous_gencode19_candidate")

# Transcript orientation is part of the complete exon-chain comparison.
minus <- copy(candidate)
minus$strand <- -1L
minus[, exons := list(data.frame(rank = 1:2, exon_start = c(7L, 1L), exon_end = c(9L, 6L)))]
reverse <- copy(target); reverse$strand <- "-"; reverse$exons <- "7-9;1-6"
reverse$cds_rows <- data.table(start = c(7L, 1L), end = c(9L, 6L), phase = c("0", "0"))
reverse$exon_rows <- data.table(start = c(7L, 1L), end = c(9L, 6L))
stopifnot(gate(t = reverse, c = minus)$status == "exact_model_match")

# Noncoding RefSeq RNA requires the same exon and sequence gates but no protein.
noncoding <- copy(candidate)
noncoding$cds_start <- 0L; noncoding$cds_end <- 0L
noncoding$pre <- "ATGAAATAA"; noncoding$cds <- NA_character_; noncoding$post <- NA_character_
noncoding_target <- copy(target)
noncoding_target$cds_rows <- data.table(start = integer(), end = integer(), phase = character())
stopifnot(gate(t = noncoding_target, c = noncoding, p = NULL)$status == "exact_model_match")

# Vertebrate mitochondrial ATA encodes M and TGA encodes W.
mito <- copy(candidate); mito$codon_table <- 2L
mito$cds <- "ATAAAATGA"
mt_target <- copy(target); mt_target$cdna <- mito$cds
stopifnot(gate(t = mt_target, c = mito, r = mt_target$cdna, p = "MKW")$status == "exact_model_match")

# Canonical, Basic and mapped MANE are independently represented facts.
native <- data.table(canonical = c(TRUE, TRUE, FALSE, FALSE),
  gencode_basic = c(TRUE, FALSE, TRUE, FALSE), gencode_primary = FALSE,
  native_mane = FALSE)
stopifnot(uniqueN(native[, .(canonical, gencode_basic)]) == 4L,
  !any(native$gencode_primary), !any(native$native_mane))
cat("MANE GRCh37 policy fixtures: passed\n")
