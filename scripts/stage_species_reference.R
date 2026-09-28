#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(Rsamtools)
  library(GenomicRanges)
  library(Biostrings)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: stage_species_reference.R FASTA DATABASE OUTPUT_PARQUET", call. = FALSE)
}
fasta <- normalizePath(args[[1L]], mustWork = TRUE)
stopifnot(file.exists(paste0(fasta, ".fai")), !file.exists(args[[3L]]))
con <- dbConnect(duckdb(shared_home = FALSE), dbdir = args[[2L]])
on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, 'CREATE TABLE reference_chunks (chrom VARCHAR, "start" BIGINT, "end" BIGINT, seq VARCHAR)')
fa <- FaFile(fasta)
open(fa)
on.exit(close(fa), add = TRUE)
idx <- scanFaIndex(fa)
rows <- vector("list", 32L)
used <- 0L
for (i in seq_along(idx)) {
  chrom <- as.character(seqnames(idx)[i])
  length <- width(idx)[i]
  for (begin in seq.int(1L, length, by = 1000000L)) {
    end <- min(begin + 999999L, length)
    seq <- as.character(scanFa(fa, param = GRanges(chrom, IRanges(begin, end))))[[1L]]
    stopifnot(nchar(seq) == end - begin + 1L)
    used <- used + 1L
    rows[[used]] <- data.frame(chrom = chrom, start = begin - 1L,
                               end = end, seq = seq)
    if (used == length(rows)) {
      dbAppendTable(con, "reference_chunks", do.call(rbind, rows))
      used <- 0L
    }
  }
}
if (used > 0L) {
  dbAppendTable(con, "reference_chunks", do.call(rbind, rows[seq_len(used)]))
}
dbExecute(con, paste0('COPY reference_chunks TO ',
  as.character(dbQuoteString(con, args[[3L]])), ' (FORMAT PARQUET)'))
print(dbGetQuery(con, 'SELECT count(*) chunks, count(DISTINCT chrom) contigs, sum("end" - "start") bases FROM reference_chunks'))
