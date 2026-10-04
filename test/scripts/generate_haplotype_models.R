#!/usr/bin/env Rscript
# Derive the DuckVEP model tables (transcripts and exons, TSV) for a duckvep-coding fixture straight from its
# FASTA and GFF3, the same files bcftools csq reads. Base R only. Usage: generate_haplotype_models.R prefix...
# where prefix is 'vertical', 'same_codon', 'frame', 'startstop' or 'nmd' under test/data/haplotype. Every transcript is a
# single-CDS protein_coding mRNA with a table-1 CDS. Exons come from the GFF3 exon rows, so a UTR (exon rows wider than the
# CDS rows, or exons with no CDS) puts the CDS inside the spliced transcript: exon_cdna_* are mRNA (cDNA) coordinates with
# the 5' UTR included, the Ensembl phase is -1 where the exon boundary is not coding, and the transcripts table gains
# cds_start/cds_end (genomic CDS bounds) and pre_cds_sequence/post_cds_sequence (the transcript-orientation 5'/3' UTR
# bases). Fixtures without a UTR keep their original columns.
root <- 'test/data/haplotype/'
revcomp <- function(s) paste(rev(c(A = 'T', C = 'G', G = 'C', T = 'A')[strsplit(s, '')[[1L]]]), collapse = '')
for (prefix in commandArgs(TRUE)) {
    lines <- readLines(paste0(root, prefix, '.fa'))
    names_ <- sub('^>', '', lines[c(TRUE, FALSE)]); genome <- setNames(lines[c(FALSE, TRUE)], names_)
    features <- read.delim(paste0(root, prefix, '.gff3'), header = FALSE, comment.char = '#', stringsAsFactors = FALSE,
        col.names = c('seq', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'))
    transcripts <- list(); exons <- list(); any_utr <- FALSE
    for (name in names_) {
        cds_rows <- features[features$seq == name & features$type == 'CDS', ]
        exon_rows <- features[features$seq == name & features$type == 'exon', ]
        strand <- unique(cds_rows$strand); stopifnot(length(strand) == 1L, strand %in% c('+', '-'))
        minus <- strand == '-'
        cds_rows <- cds_rows[order(cds_rows$start, decreasing = minus), ]
        exon_rows <- exon_rows[order(exon_rows$start, decreasing = minus), ]
        pieces <- substring(genome[[name]], cds_rows$start, cds_rows$end)
        if (minus) pieces <- vapply(pieces, revcomp, character(1L))
        lengths <- nchar(pieces); ends <- cumsum(lengths); starts <- ends - lengths + 1L
        # GFF3 CDS phase = bases to skip to reach a codon start; Ensembl model phase = bases already consumed.
        stopifnot(all(cds_rows$phase == (3L - (starts - 1L) %% 3L) %% 3L))
        region <- match(name, names_) - 1L
        utr <- !(nrow(exon_rows) == nrow(cds_rows) && all(exon_rows$start == cds_rows$start) &&
                 all(exon_rows$end == cds_rows$end))
        if (!utr) {
            row <- data.frame(case = name, seq_region = region,
                strand = if (minus) -1L else 1L, transcript_start = min(cds_rows$start),
                transcript_end = max(cds_rows$end), cds_sequence = paste(pieces, collapse = ''),
                cds_start = min(cds_rows$start), cds_end = max(cds_rows$end),
                pre_cds_sequence = NA_character_, post_cds_sequence = NA_character_)
            exon_frame <- data.frame(seq_region = region, exon_start = cds_rows$start,
                exon_end = cds_rows$end, exon_cdna_start = starts, exon_cdna_end = ends,
                phase = (starts - 1L) %% 3L, end_phase = ends %% 3L)
        } else {
            elen <- exon_rows$end - exon_rows$start + 1L
            eend <- cumsum(elen); estart <- eend - elen + 1L
            # Coding bases of each exon (CDS rows lie inside one exon each) and those before it.
            coding <- vapply(seq_len(nrow(exon_rows)), function(k) sum(pmax(0L, pmin(exon_rows$end[k], cds_rows$end) -
                pmax(exon_rows$start[k], cds_rows$start) + 1L)), integer(1L))
            before <- cumsum(coding) - coding
            # Whether the exon's first and last base (transcript orientation) are coding.
            first_base <- if (minus) exon_rows$end else exon_rows$start
            last_base <- if (minus) exon_rows$start else exon_rows$end
            in_cds <- function(g) vapply(g, function(x) any(x >= cds_rows$start & x <= cds_rows$end), logical(1L))
            row <- data.frame(case = name, seq_region = region,
                strand = if (minus) -1L else 1L, transcript_start = min(exon_rows$start),
                transcript_end = max(exon_rows$end), cds_sequence = paste(pieces, collapse = ''),
                cds_start = min(cds_rows$start), cds_end = max(cds_rows$end),
                pre_cds_sequence = NA_character_, post_cds_sequence = NA_character_)
            epieces <- substring(genome[[name]], exon_rows$start, exon_rows$end)
            if (minus) epieces <- vapply(epieces, revcomp, character(1L))
            mrna <- paste(epieces, collapse = '')
            g <- if (minus) cds_rows$end[1L] else cds_rows$start[1L]
            k <- which(g >= exon_rows$start & g <= exon_rows$end)
            offset0 <- estart[k] - 1L + (if (minus) exon_rows$end[k] - g else g - exon_rows$start[k])
            cds_text <- paste(pieces, collapse = '')
            stopifnot(substr(mrna, offset0 + 1L, offset0 + nchar(cds_text)) == cds_text)
            # An empty flank is written as NA ('.'), so no row ends in an empty field.
            flank5 <- substr(mrna, 1L, offset0); flank3 <- substring(mrna, offset0 + nchar(cds_text) + 1L)
            row$pre_cds_sequence <- if (nzchar(flank5)) flank5 else NA_character_
            row$post_cds_sequence <- if (nzchar(flank3)) flank3 else NA_character_
            exon_frame <- data.frame(seq_region = region, exon_start = exon_rows$start,
                exon_end = exon_rows$end, exon_cdna_start = estart, exon_cdna_end = eend,
                phase = ifelse(in_cds(first_base), before %% 3L, -1L),
                end_phase = ifelse(in_cds(last_base), (before + coding) %% 3L, -1L))
        }
        any_utr <- any_utr || utr
        transcripts[[length(transcripts) + 1L]] <- row
        exons[[length(exons) + 1L]] <- exon_frame
    }
    write_checked <- function(frame, path) {
        text <- capture.output(write.table(frame, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
        if (file.exists(path) && !identical(readLines(path, warn = FALSE), text))
            stop('model differs from independent generator: ', path)
        writeLines(text, path)
    }
    transcripts <- do.call(rbind, transcripts)
    if (!any_utr) transcripts <- transcripts[, setdiff(names(transcripts), c('cds_start', 'cds_end', 'pre_cds_sequence', 'post_cds_sequence'))]
    write_checked(transcripts, paste0(root, prefix, '_transcripts.tsv'))
    write_checked(do.call(rbind, exons), paste0(root, prefix, '_exons.tsv'))
    cat(prefix, ':', nrow(transcripts), 'transcripts,', sum(vapply(exons, nrow, 1L)), 'exons\n')
}
