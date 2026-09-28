# Pure, accession-scoped MANE validation. Only the caller resolves RefSeq accessions.
mane_candidates <- function(model, accession_root) model[list(accession_root), nomatch = 0L]

mane_supported_loci <- function(rows, candidate_contig = NULL) {
  eligible <- rows[target_ok == TRUE]
  if (nrow(eligible) == 0L) return(rows)
  if (!is.null(candidate_contig) && any(eligible$model_contig == candidate_contig))
    return(eligible[model_contig == candidate_contig])
  eligible
}

mane_model_cds_geometry <- function(ex, cds_start, cds_end) {
  if (is.na(cds_start) || is.na(cds_end) || cds_start == 0) return("")
  ex <- ex[order(ex$rank), ]
  start <- pmax(ex$exon_start, cds_start)
  end <- pmin(ex$exon_end, cds_end)
  keep <- start <= end
  start <- start[keep]; end <- end[keep]
  lens <- end - start + 1
  prior <- c(0, head(cumsum(lens), -1L))
  paste(paste(start, end, (3 - prior %% 3) %% 3, sep = "-"), collapse = ";")
}

mane_cds_from_rna <- function(target, ref_rna) {
  ex <- target$exon_rows
  cds <- target$cds_rows
  if (!nrow(cds)) return("")
  offsets <- c(0L, head(cumsum(ex$end - ex$start + 1L), -1L))
  starts <- vapply(seq_len(nrow(cds)), function(i) {
    j <- which(ex$start <= cds$start[i] & ex$end >= cds$end[i])[1]
    if (is.na(j)) return(NA_integer_)
    offsets[j] + if (target$strand == "+") cds$start[i] - ex$start[j] + 1L else ex$end[j] - cds$end[i] + 1L
  }, integer(1))
  if (anyNA(starts)) return(NA_character_)
  paste0(substring(ref_rna, starts, starts + cds$end - cds$start), collapse = "")
}

mane_protein_matches <- function(target, ref_rna, protein, codon_table) {
  if (is.null(protein)) return(nrow(target$cds_rows) == 0L)
  if (nrow(target$cds_rows) == 0L) return(FALSE)
  cds <- mane_cds_from_rna(target, ref_rna)
  if (is.na(cds) || nchar(cds) == 0L || nchar(cds) %% 3L != 0L) return(FALSE)
  if (is.na(codon_table)) return(NA)
  code <- if (codon_table == 2L) Biostrings::getGeneticCode("SGC1") else Biostrings::GENETIC_CODE
  translated <- as.character(Biostrings::translate(Biostrings::DNAString(cds),
    genetic.code = code, if.fuzzy.codon = "X"))
  identical(sub("\\*$", "", translated), sub("\\*$", "", protein))
}

mane_pair_gate <- function(target, candidates, ref_rna, protein) {
  evidence <- list(exon_chain_match = NA, utr_exon_chain_match = NA,
    cds_phase_match = NA, spliced_sequence_match = NA,
    translated_sequence_match = NA, reference_difference = NA,
    reference_difference_bases = NA_integer_, cds_reference_difference_bases = NA_integer_)
  finish <- function(status, selected = NULL) list(status = status, selected = selected, evidence = evidence)
  if (!nrow(candidates)) return(finish("refseq_only_no_gencode19_match"))
  strand <- if (target$strand == "+") "1" else "-1"
  locus <- candidates$seq_region_name == target$contig & as.character(candidates$strand) == strand
  candidates <- candidates[which(locus)]
  if (!nrow(candidates)) return(finish("geometry_mismatch"))
  exon <- vapply(candidates$exons, function(ex)
    identical(paste(paste(ex$exon_start, ex$exon_end, sep = "-"), collapse = ";"), target$exons),
    logical(1))
  evidence$exon_chain_match <- any(exon)
  evidence$utr_exon_chain_match <- any(exon)
  target_cds <- paste(paste(target$cds_rows$start, target$cds_rows$end,
                            target$cds_rows$phase, sep = "-"), collapse = ";")
  cds <- vapply(seq_len(nrow(candidates)), function(i)
    identical(mane_model_cds_geometry(candidates$exons[[i]], candidates$cds_start[i],
                                      candidates$cds_end[i]), target_cds), logical(1))
  evidence$cds_phase_match <- any(cds)
  reference <- vapply(seq_len(nrow(candidates)), function(i) {
    transcript <- paste0(ifelse(is.na(candidates$pre[i]), "", candidates$pre[i]),
                         ifelse(is.na(candidates$cds[i]), "", candidates$cds[i]),
                         ifelse(is.na(candidates$post[i]), "", candidates$post[i]))
    identical(transcript, target$cdna)
  }, logical(1))
  evidence$spliced_sequence_match <- any(reference)

  if (!is.null(ref_rna) && nchar(ref_rna) == nchar(target$cdna)) {
    difference <- sum(strsplit(ref_rna, "", fixed = TRUE)[[1]] !=
                      strsplit(target$cdna, "", fixed = TRUE)[[1]])
    evidence$reference_difference_bases <- difference
    evidence$reference_difference <- difference > 0L
    if (nrow(target$cds_rows)) {
      ref_cds <- mane_cds_from_rna(target, ref_rna)
      genome_cds <- mane_cds_from_rna(target, target$cdna)
      if (!is.na(ref_cds) && !is.na(genome_cds) && nchar(ref_cds) == nchar(genome_cds))
        evidence$cds_reference_difference_bases <- sum(
          strsplit(ref_cds, "", fixed = TRUE)[[1]] !=
          strsplit(genome_cds, "", fixed = TRUE)[[1]])
    }
  }
  # Strict matches retain RNA-to-protein validation; coding-only matches
  # validate the reference-genome CDS against the pinned RefSeq protein.
  strict <- which(exon & cds & reference)
  coding <- which(cds & !exon & nrow(target$cds_rows) > 0L)
  scope <- if (length(strict) > 0L) strict else if (length(coding) > 0L) coding else seq_len(nrow(candidates))
  sequence <- if (length(strict) > 0L) ref_rna else target$cdna
  if (!is.null(sequence) && (!is.null(protein) || !nrow(target$cds_rows)))
    evidence$translated_sequence_match <- any(vapply(scope, function(i)
      mane_protein_matches(target, sequence, protein, candidates$codon_table[i]), logical(1)))
  if (length(strict) > 0L) {
    if (is.na(evidence$reference_difference_bases)) return(finish("sequence_mismatch"))
    if (!isTRUE(evidence$translated_sequence_match)) return(finish("translation_mismatch"))
    if (length(strict) != 1L) return(finish("ambiguous_gencode19_candidate"))
    return(finish("exact_model_match", candidates[strict]))
  }
  if (length(coding) > 0L && sum(cds) != 1L) return(finish("ambiguous_gencode19_candidate"))
  if (length(coding) == 1L && isTRUE(evidence$translated_sequence_match) &&
      identical(evidence$cds_reference_difference_bases, 0L) && !is.null(protein)) {
    evidence$exon_chain_match <- FALSE
    evidence$utr_exon_chain_match <- FALSE
    evidence$spliced_sequence_match <- reference[coding]
    return(finish("cds_exact_utr_differs", candidates[coding]))
  }
  if (!any(exon)) return(finish("geometry_mismatch"))
  if (!any(exon & cds)) return(finish("cds_phase_mismatch"))
  if (!any(exon & cds & reference)) return(finish("sequence_mismatch"))
  finish("translation_mismatch")
}
