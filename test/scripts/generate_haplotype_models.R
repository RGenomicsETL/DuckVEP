#!/usr/bin/env Rscript
# Derive the DuckVEP model tables (transcripts and exons, TSV) for a coding-v1 fixture straight from its
# FASTA and GFF3, the same files bcftools csq reads. Base R only. Usage: generate_haplotype_models.R prefix...
# where prefix is 'vertical', 'same_codon', 'frame' or 'startstop' under test/data/haplotype. Every transcript is a single-CDS
# protein_coding mRNA with a table-1 CDS spanning its exons (no UTR).
root <- 'test/data/haplotype/'
revcomp <- function(s) paste(rev(c(A = 'T', C = 'G', G = 'C', T = 'A')[strsplit(s, '')[[1L]]]), collapse = '')
for (prefix in commandArgs(TRUE)) {
    lines <- readLines(paste0(root, prefix, '.fa'))
    names_ <- sub('^>', '', lines[c(TRUE, FALSE)]); genome <- setNames(lines[c(FALSE, TRUE)], names_)
    features <- read.delim(paste0(root, prefix, '.gff3'), header = FALSE, comment.char = '#', stringsAsFactors = FALSE,
        col.names = c('seq', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'))
    transcripts <- list(); exons <- list()
    for (name in names_) {
        cds_rows <- features[features$seq == name & features$type == 'CDS', ]
        strand <- unique(cds_rows$strand); stopifnot(length(strand) == 1L, strand %in% c('+', '-'))
        cds_rows <- cds_rows[order(cds_rows$start, decreasing = strand == '-'), ]
        pieces <- substring(genome[[name]], cds_rows$start, cds_rows$end)
        if (strand == '-') pieces <- vapply(pieces, revcomp, character(1L))
        lengths <- nchar(pieces); ends <- cumsum(lengths); starts <- ends - lengths + 1L
        # GFF3 CDS phase = bases to skip to reach a codon start; Ensembl model phase = bases already consumed.
        stopifnot(all(cds_rows$phase == (3L - (starts - 1L) %% 3L) %% 3L))
        region <- match(name, names_) - 1L
        transcripts[[length(transcripts) + 1L]] <- data.frame(case = name, seq_region = region,
            strand = if (strand == '-') -1L else 1L, transcript_start = min(cds_rows$start),
            transcript_end = max(cds_rows$end), cds_sequence = paste(pieces, collapse = ''))
        exons[[length(exons) + 1L]] <- data.frame(seq_region = region, exon_start = cds_rows$start,
            exon_end = cds_rows$end, exon_cdna_start = starts, exon_cdna_end = ends,
            phase = (starts - 1L) %% 3L, end_phase = ends %% 3L)
    }
    write_checked <- function(frame, path) {
        text <- capture.output(write.table(frame, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
        if (file.exists(path) && !identical(readLines(path, warn = FALSE), text))
            stop('model differs from independent generator: ', path)
        writeLines(text, path)
    }
    write_checked(do.call(rbind, transcripts), paste0(root, prefix, '_transcripts.tsv'))
    write_checked(do.call(rbind, exons), paste0(root, prefix, '_exons.tsv'))
    cat(prefix, ':', length(transcripts), 'transcripts,', sum(vapply(exons, nrow, 1L)), 'exons\n')
}
