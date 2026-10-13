# Independent SeqEdit and standard-code translation fixture for the native stream tests.
raw <- strsplit("AAAAAAATGAAAGGGTAA", "", fixed = TRUE)[[1L]]
seqedit <- function(sequence, start, end, alternate) {
  c(sequence[seq_len(start - 1L)], strsplit(alternate, "", fixed = TRUE)[[1L]],
    sequence[seq.int(end + 1L, length(sequence))])
}
translate <- function(sequence, start = 7L, end = length(sequence), table = 1L) {
  genetic_code <- c(
    TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S",
    TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W",
    CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P",
    CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R",
    ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T",
    AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R",
    GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A",
    GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G")
  if (table == 2L) genetic_code[c("AGA", "AGG", "ATA", "TGA")] <- c("*", "*", "M", "W")
  cds <- sequence[seq.int(start, end)]
  triplets <- vapply(seq.int(1L, length(cds) - 2L, 3L),
                     function(i) paste0(cds[i:(i + 2L)], collapse = ""), "")
  amino_acids <- unname(genetic_code[triplets])
  first_stop <- match("*", amino_acids)
  paste0(amino_acids[seq_len(if (is.na(first_stop)) length(amino_acids) else first_stop)],
         collapse = "")
}
complement <- function(bases) chartr("ACGT", "TGCA", bases)
rna_reference <- seqedit(raw, 10L, 10L, "G")
for (strand in c(1L, -1L)) {
  genomic <- if (strand == 1L) raw else rev(complement(raw))
  genomic_position <- if (strand == 1L) 13L else length(raw) - 12L
  genomic_ref <- if (strand == 1L) "G" else "C"
  genomic_alt <- if (strand == 1L) "A" else "T"
  stopifnot(identical(genomic[genomic_position], genomic_ref))
  genomic[genomic_position] <- genomic_alt
  variant <- if (strand == 1L) genomic else rev(complement(genomic))
  rna_variant <- seqedit(variant, 10L, 10L, "G")
  stopifnot(identical(paste0(rna_reference, collapse = ""), "AAAAAAATGGAAGGGTAA"),
            identical(paste0(rna_variant, collapse = ""), "AAAAAAATGGAAAGGTAA"),
            identical(translate(rna_reference), "MEG*"),
            identical(translate(rna_variant), "MER*"),
            identical(translate(rna_reference, table = 2L), "MEG*"),
            identical(translate(rna_variant, table = 2L), "ME*"))
}
ensembl_raw <- strsplit("ATGCCCTAAAGG", "", fixed = TRUE)[[1L]]
ensembl_variant <- ensembl_raw
ensembl_variant[7L] <- "C"
ensembl_reference_rna <- seqedit(ensembl_raw, 4L, 5L, "TT")
ensembl_variant_rna <- seqedit(ensembl_variant, 4L, 5L, "TT")
stopifnot(identical(paste0(ensembl_reference_rna[1L:9L], collapse = ""), "ATGTTCTAA"),
          identical(paste0(ensembl_variant_rna[1L:9L], collapse = ""), "ATGTTCCAA"),
          identical(translate(ensembl_reference_rna, 1L, 12L, 2L), "MF*"),
          identical(translate(ensembl_variant_rna, 1L, 12L, 2L), "MFQ*"))
cat("RNA/CDS/peptides: ATGGAAGGGTAA MEG* / ATGGAAAGGTAA MER*; Ensembl MF* / MFQ*\n")
