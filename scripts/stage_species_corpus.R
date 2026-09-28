#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(Rsamtools)
  library(GenomicRanges)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(5L, 6L) || !nzchar(args[[5L]])) {
  stop("usage: stage_species_corpus.R MODEL_DB FASTA OUTPUT_VCF SEED SPECIES [SAMPLE_PER_GROUP]", call. = FALSE)
}
species <- args[[5L]]
sample_size <- if (length(args) == 6L) as.integer(args[[6L]]) else 32L
stopifnot(!is.na(sample_size), sample_size > 0L)
model <- normalizePath(args[[1L]], mustWork = TRUE)
fasta <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
seed <- as.integer(args[[4L]])
stopifnot(!is.na(seed), seed >= 0L, file.exists(paste0(fasta, ".fai")),
          !file.exists(output))
con <- dbConnect(duckdb(shared_home = FALSE), dbdir = model, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE))
set.seed(seed)
regions <- dbGetQuery(con, "SELECT seq_region_name AS chrom, sequence_length FROM model_regions")
transcripts <- dbGetQuery(con, "
  SELECT transcript_stable_id AS tx, seq_region_name AS chrom,
         transcript_start AS start, transcript_end AS end, strand,
         transcript_biotype AS biotype, codon_table,
         cds_start, cds_end, CAST(cds_sequence AS VARCHAR) AS cds,
         exons[1].exon_start AS exon_start, exons[1].exon_end AS exon_end,
         length(exons) AS exon_count
  FROM model_transcripts ORDER BY transcript_stable_id")
tables <- sort(unique(transcripts$codon_table[!is.na(transcripts$codon_table)]))
stopifnot(length(tables) > 0L)

# One bounded sample per observed contig/biotype/strand, with independent
# positions for each allele shape. The seed also controls grouping ranks.
groups <- split(seq_len(nrow(transcripts)),
                interaction(transcripts$chrom, transcripts$biotype,
                            transcripts$strand, drop = TRUE))
chosen <- unlist(lapply(groups, function(i) i[sample.int(length(i), min(length(i), sample_size))]),
                 use.names = FALSE)
fa <- FaFile(fasta)
open(fa)
on.exit(close(fa), add = TRUE)
base_at <- function(chrom, start, end = start) {
  as.character(scanFa(fa, param = GRanges(chrom, IRanges(start, end))))[[1L]]
}
variants <- list()
add <- function(chrom, pos, reference, alternate, source, table, tx) {
  if (pos < 1L || !grepl("^[ACGT]+$", reference) ||
      !grepl("^[ACGT]+$", alternate) || identical(reference, alternate)) return()
  variants[[length(variants) + 1L]] <<- data.frame(
    chrom = chrom, pos = pos, ref = reference, alt = alternate,
    source = source, table = table, tx = tx)
}
other <- function(base) sample(setdiff(c("A", "C", "G", "T"), base), 1L)
for (i in chosen) {
  t <- transcripts[i, ]
  width <- t$end - t$start + 1L
  if (width < 10L) next
  pos <- t$start + 1L + sample.int(width - 5L, 1L)
  ref <- base_at(t$chrom, pos, pos + 2L)
  if (!grepl("^[ACGT]{3}$", ref)) next
  add(t$chrom, pos, substr(ref, 1L, 1L), other(substr(ref, 1L, 1L)),
      "SNV", t$codon_table, t$tx)
  add(t$chrom, pos, substr(ref, 1L, 2L),
      paste0(other(substr(ref, 1L, 1L)), other(substr(ref, 2L, 2L))),
      "MNV", t$codon_table, t$tx)
  add(t$chrom, pos, substr(ref, 1L, 1L),
      paste0(substr(ref, 1L, 1L), other(substr(ref, 1L, 1L))),
      "insertion", t$codon_table, t$tx)
  add(t$chrom, pos, ref, substr(ref, 1L, 1L),
      "deletion", t$codon_table, t$tx)
}

# Translation-table witnesses are mapped from single-exon CDS coordinates.
# TGG-to-TGA witnesses distinguish table 1 from tables 2 and 4;
# table 11 also admits ATG-to-GTG initiation witnesses.
for (table in tables) {
  candidates <- transcripts[transcripts$codon_table == table &
                              transcripts$exon_count == 1L &
                              !is.na(transcripts$cds) &
                              nchar(transcripts$cds) >= 12L, ]
  if (table == 11L) {
    candidates <- candidates[startsWith(candidates$cds, "ATG"), ]
  } else {
    candidates <- candidates[grepl("TGG", candidates$cds, fixed = TRUE), ]
  }
  stopifnot(nrow(candidates) > 0L)
  candidates <- candidates[order(candidates$tx), ]
  for (j in seq_len(min(nrow(candidates), if (table == 11L) 30L else 12L))) {
    t <- candidates[j, ]
    offset <- if (table == 11L) 0L else {
      codons <- substring(t$cds, seq.int(1L, nchar(t$cds) - 2L, by = 3L),
                          seq.int(3L, nchar(t$cds), by = 3L))
      hit <- which(codons == "TGG")
      if (!length(hit)) next
      (hit[[1L]] - 1L) * 3L + 2L
    }
    pos <- if (t$strand == 1L) t$cds_start + offset else t$cds_end - offset
    ref <- base_at(t$chrom, pos)
    expected <- if (table == 11L) "A" else "G"
    alt <- if (table == 11L) "G" else "A"
    if (t$strand == -1L) {
      expected <- as.character(complement(DNAString(expected)))
      alt <- as.character(complement(DNAString(alt)))
    }
    stopifnot(identical(ref, expected))
    add(t$chrom, pos, ref, alt, "codon_witness", table, t$tx)
  }
}
variants <- do.call(rbind, variants)
variants <- variants[order(variants$chrom, variants$pos, variants$ref,
                           variants$alt, variants$source), ]
variants <- variants[!duplicated(variants[c("chrom", "pos", "ref", "alt")]), ]
prefix <- paste0(toupper(substr(strsplit(species, "_", fixed = TRUE)[[1L]], 1L, 1L)),
                 collapse = "")
variants$id <- sprintf("%s%06d", prefix, seq_len(nrow(variants)))
missing_biotypes <- setdiff(unique(transcripts$biotype),
                            transcripts$biotype[match(unique(variants$tx), transcripts$tx)])
if (length(missing_biotypes)) stop("missing biotypes: ", paste(missing_biotypes, collapse = ", "))
stopifnot(all(c("SNV", "MNV", "insertion", "deletion", "codon_witness") %in% variants$source),
          all(tables %in% variants$table[variants$source == "codon_witness"]))
writeLines(c("##fileformat=VCFv4.2",
             paste0("##reference=", fasta),
             paste0("##duckvep_seed=", seed),
             paste0("##contig=<ID=", regions$chrom, ",length=", regions$sequence_length, ">"),
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), output)
write.table(data.frame(variants$chrom, variants$pos, variants$id,
                       variants$ref, variants$alt, ".", "PASS", "."),
            output, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE, append = TRUE)
write.table(variants, paste0(output, ".provenance.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(table(variants$source, variants$table, useNA = "ifany"))
cat("variants:", nrow(variants), "biotypes:",
    paste(sort(unique(transcripts$biotype)), collapse = ","), "\n")
