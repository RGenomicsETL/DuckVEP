#!/usr/bin/env Rscript

reverse_complement <- function(sequence) {
  bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  paste0(rev(c(A = "T", C = "G", G = "C", T = "A")[bases]), collapse = "")
}

compose_substitutions <- function(reference, edits) {
  result <- strsplit(reference, "", fixed = TRUE)[[1L]]
  for (edit in edits) {
    position <- edit[["position"]]
    stopifnot(identical(result[[position]], edit[["reference"]]))
    result[[position]] <- edit[["alternate"]]
  }
  paste0(result, collapse = "")
}

classify_final <- function(reference, alternate, coding_start, coding_end) {
  changed <- which(strsplit(reference, "", fixed = TRUE)[[1L]] !=
                     strsplit(alternate, "", fixed = TRUE)[[1L]])
  if (length(changed) == 0L) return("no_change")
  if (is.na(coding_start)) return("non_coding_transcript_exon_variant")
  if (all(changed < coding_start)) return("5_prime_UTR_variant")
  if (all(changed > coding_end)) return("3_prime_UTR_variant")
  "unsupported_mixed_feature"
}

reference <- "AACCGGTTAA"
coding_start <- 4L
coding_end <- 7L
utr5 <- compose_substitutions(reference, list(
  list(position = 2L, reference = "A", alternate = "T")
))
utr3 <- compose_substitutions(reference, list(
  list(position = 9L, reference = "A", alternate = "C")
))
restored <- compose_substitutions(reference, list(
  list(position = 2L, reference = "A", alternate = "T"),
  list(position = 2L, reference = "T", alternate = "A")
))
reverse_genomic <- reverse_complement(reference)
reverse_genomic_alternate <- compose_substitutions(reverse_genomic, list(
  list(position = 9L, reference = "T", alternate = "A")
))
stopifnot(
  identical(utr5, "ATCCGGTTAA"),
  identical(utr3, "AACCGGTTCA"),
  identical(restored, reference),
  identical(reverse_complement(reverse_genomic_alternate), utr5),
  identical(classify_final(reference, utr5, coding_start, coding_end),
            "5_prime_UTR_variant"),
  identical(classify_final(reference, utr3, coding_start, coding_end),
            "3_prime_UTR_variant"),
  identical(classify_final(reference, "AATCGGTTAA", NA_integer_, NA_integer_),
            "non_coding_transcript_exon_variant"),
  identical(classify_final(reference, restored, coding_start, coding_end), "no_change"),
  identical(classify_final(reference, "ATCCGGTTCA", coding_start, coding_end),
            "unsupported_mixed_feature")
)
