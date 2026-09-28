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
  code <- if (codon_table == 2L) Biostrings::getGeneticCode("SGC1") else Biostrings::GENETIC_CODE
  translated <- as.character(Biostrings::translate(Biostrings::DNAString(cds),
    genetic.code = code, if.fuzzy.codon = "X"))
  identical(sub("\\*$", "", translated), sub("\\*$", "", protein))
}

mane_pair_gate <- function(target, candidates, ref_rna, protein) {
  evidence <- list(exon_chain_match = FALSE, cds_phase_match = FALSE,
    spliced_sequence_match = FALSE, translated_sequence_match = FALSE,
    reference_difference = FALSE, reference_difference_bases = NA_integer_)
  reject <- function(status) list(status = status, selected = NULL, evidence = evidence)
  if (!nrow(candidates)) return(reject("refseq_only_no_gencode19_match"))
  strand <- if (target$strand == "+") "1" else "-1"
  geometry <- vapply(seq_len(nrow(candidates)), function(i) {
    x <- candidates[i]
    if (!identical(x$seq_region_name, target$contig) ||
        !identical(as.character(x$strand), strand)) return(FALSE)
    ex <- x$exons[[1]]
    identical(paste(paste(ex$exon_start, ex$exon_end, sep = "-"), collapse = ";"), target$exons)
  }, logical(1))
  evidence$exon_chain_match <- any(geometry)
  if (!any(geometry)) return(reject("geometry_mismatch"))
  candidates <- candidates[geometry]
  cds <- vapply(seq_len(nrow(candidates)), function(i) {
    x <- candidates[i]
    mane_model_cds_geometry(x$exons[[1]], x$cds_start, x$cds_end) ==
      paste(paste(target$cds_rows$start, target$cds_rows$end,
                  target$cds_rows$phase, sep = "-"), collapse = ";")
  }, logical(1))
  evidence$cds_phase_match <- any(cds)
  if (!any(cds)) return(reject("cds_phase_mismatch"))
  candidates <- candidates[cds]
  reference <- vapply(seq_len(nrow(candidates)), function(i) {
    x <- candidates[i]
    transcript <- paste0(ifelse(is.na(x$pre), "", x$pre),
                         ifelse(is.na(x$cds), "", x$cds),
                         ifelse(is.na(x$post), "", x$post))
    identical(transcript, target$cdna)
  }, logical(1))
  evidence$spliced_sequence_match <- any(reference)
  if (!any(reference)) return(reject("sequence_mismatch"))
  candidates <- candidates[reference]
  if (is.null(ref_rna) || nchar(ref_rna) != nchar(target$cdna)) return(reject("sequence_mismatch"))
  difference <- sum(strsplit(ref_rna, "", fixed = TRUE)[[1]] !=
                    strsplit(target$cdna, "", fixed = TRUE)[[1]])
  evidence$reference_difference_bases <- difference
  evidence$reference_difference <- difference > 0L
  evidence$translated_sequence_match <- mane_protein_matches(target, ref_rna, protein,
    candidates$codon_table[1])
  if (!evidence$translated_sequence_match) return(reject("translation_mismatch"))
  if (nrow(candidates) != 1L) return(reject("ambiguous_gencode19_candidate"))
  list(status = "exact_model_match", selected = candidates[1], evidence = evidence)
}
