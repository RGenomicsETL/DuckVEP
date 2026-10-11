#!/usr/bin/env Rscript

reverse_complement <- function(sequence) {
  bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  paste0(rev(c(A = "T", C = "G", G = "C", T = "A")[bases]), collapse = "")
}

# Independent string construction on immutable, one-based reference positions.
compose_rna <- function(reference, edits, capacity = Inf) {
  edits <- edits[order(vapply(edits, `[[`, integer(1L), "position"))]
  cursor <- 1L
  pieces <- character()
  delta <- 0L
  previous_site <- 0L
  for (edit in edits) {
    position <- edit$position
    nr <- nchar(edit$reference)
    stopifnot(position >= cursor, position > previous_site,
      position <= nchar(reference) + 1L,
      position + nr - 1L <= nchar(reference),
      identical(substr(reference, position, position + nr - 1L), edit$reference))
    if (position > cursor) pieces <- c(pieces, substr(reference, cursor, position - 1L))
    pieces <- c(pieces, edit$alternate)
    cursor <- position + nr
    previous_site <- position + max(1L, nr) - 1L
    delta <- delta + nchar(edit$alternate) - nr
  }
  if (cursor <= nchar(reference)) pieces <- c(pieces, substr(reference, cursor, nchar(reference)))
  sequence <- paste0(pieces, collapse = "")
  stopifnot(nchar(sequence) == nchar(reference) + delta, nchar(sequence) <= capacity)
  list(sequence = sequence, delta = delta)
}

rna_edit_cases <- function(base) list(
  insertion = list(list(position = 2L, reference = base, alternate = paste0(base, "G"))),
  deletion = list(list(position = 2L, reference = paste0(base, base), alternate = base)),
  delins = list(list(position = 2L, reference = paste0(base, base), alternate = "GTG")),
  compensating_change = list(list(position = 2L, reference = base, alternate = paste0(base, "G")),
    list(position = 5L, reference = paste0(base, base), alternate = base)),
  compensating_no_change = list(list(position = 2L, reference = base, alternate = paste0(base, base)),
    list(position = 5L, reference = paste0(base, base), alternate = base)))

reference <- "AAAAAAATGAAATAACCCCC"
for (name in names(rna_edit_cases("A"))) {
  edits <- rna_edit_cases("A")[[name]]
  result <- compose_rna(reference, edits)
  delta <- c(insertion = 1L, deletion = -1L, delins = 1L,
    compensating_change = 0L, compensating_no_change = 0L)[[name]]
  stopifnot(result$delta == delta, substr(result$sequence, 7L + delta, 15L + delta) == "ATGAAATAA")
  genomic <- reverse_complement(reference)
  reversed <- lapply(edits, function(edit) list(
    position = nchar(reference) - edit$position - nchar(edit$reference) + 2L,
    reference = reverse_complement(edit$reference), alternate = reverse_complement(edit$alternate)))
  stopifnot(reverse_complement(compose_rna(genomic, reversed)$sequence) == result$sequence)
  if (name == "compensating_no_change") stopifnot(result$sequence == reference)
  short <- try(compose_rna(reference, edits, nchar(result$sequence) - 1L), silent = TRUE)
  stopifnot(inherits(short, "try-error"))
}
for (bad in list(
  list(list(position = 2L, reference = "C", alternate = "CG")),
  list(list(position = 2L, reference = "AA", alternate = "C"),
    list(position = 3L, reference = "A", alternate = "T"))))
  stopifnot(inherits(try(compose_rna(reference, bad), silent = TRUE), "try-error"))
stopifnot(compose_rna(reference, list(list(position = 21L, reference = "", alternate = "T")))$sequence ==
  paste0(reference, "T"))
cat("Independent RNA string oracle: positions, deltas, sequence, both strands, compensation, conflict and capacity verified\n")
