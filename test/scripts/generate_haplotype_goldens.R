#!/usr/bin/env Rscript
# Self-contained coding-transcript fixture. No DuckVEP functions or model data.
out <- 'test/data/haplotype'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base <- paste0('ATG', 'AGC', paste(rep('CAA', 37L), collapse = ''), 'TAA')
stopifnot(nchar(base) == 120L, substr(base, 49L, 51L) == 'CAA')
frame <- base
substr(frame, 24L, 26L) <- 'TAA' # overwritten below: reference remains coding
substr(frame, 22L, 27L) <- 'CATAAA' # displaced reading frame makes TAA at 24..26
stopifnot(substr(frame, 24L, 26L) == 'TAA', substr(frame, 22L, 24L) == 'CAT')
v <- function(pos, ref, alt) c(pos = pos, ref = ref, alt = alt)
cases <- list(
  cis = list(seq = base, j = NA_integer_, lane = 0L,
             edits = list(v(4, 'A', 'T'), v(5, 'G', 'C'))),
  trans = list(seq = base, j = NA_integer_, lane = c(0L, 1L),
               edits = list(v(4, 'A', 'T'), v(5, 'G', 'C'))),
  frame_open = list(seq = frame, j = NA_integer_, lane = 0L,
                    edits = list(v(12, 'A', 'AT'))),
  restored_before_stop = list(seq = frame, j = NA_integer_, lane = 0L,
                              edits = list(v(12, 'A', 'AT'), v(19, 'CA', 'C'))),
  restored_after_stop = list(seq = frame, j = NA_integer_, lane = 0L,
                             edits = list(v(12, 'A', 'AT'), v(34, 'CA', 'C'))),
  start_lost = list(seq = base, j = NA_integer_, lane = 0L, edits = list(v(1, 'A', 'C'))),
  stop_lost = list(seq = base, j = NA_integer_, lane = 0L, edits = list(v(118, 'T', 'C'))),
  stop_retained = list(seq = base, j = NA_integer_, lane = 0L, edits = list(v(120, 'A', 'G'))),
  nmd49 = list(seq = base, j = 100L, lane = 0L, edits = list(v(49, 'C', 'T'))),
  nmd50 = list(seq = base, j = 101L, lane = 0L, edits = list(v(49, 'C', 'T'))),
  nmd51 = list(seq = base, j = 102L, lane = 0L, edits = list(v(49, 'C', 'T'))),
  intronless = list(seq = base, j = NA_integer_, lane = 0L, edits = list(v(49, 'C', 'T'))),
  shifted_junction = list(seq = base, j = 100L, lane = 0L,
                          edits = list(v(49, 'C', 'T'), v(60, 'A', 'ACAA')))
)
# Standard code; alternatives are translated from the edited cDNA, not from
# independently called alleles.
codons <- c('TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
            'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
            'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
            'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG')
aa <- strsplit('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '')[[1L]]
code <- setNames(aa, codons)
translate <- function(s) {
    starts <- seq.int(1L, nchar(s) - 2L, 3L)
    unname(code[substring(s, starts, starts + 2L)])
}
edit_lane <- function(seq, edits, lanes, lane) {
    for (i in rev(seq_along(edits))) {
        if (!(lane %in% lanes[[i]])) next
        e <- edits[[i]]; p <- as.integer(e['pos'])
        stopifnot(substr(seq, p, p + nchar(e['ref']) - 1L) == e['ref'])
        seq <- paste0(substr(seq, 1L, p - 1L), e['alt'],
                      substr(seq, p + nchar(e['ref']), nchar(seq)))
    }
    seq
}
records <- character(); gff <- c('##gff-version 3'); fasta <- character()
gold <- list()
for (name in names(cases)) {
    item <- cases[[name]]; j <- item$j
    cDNA <- item$seq
    if (is.na(j)) {
        genome <- cDNA; exon <- matrix(c(1L, nchar(cDNA)), ncol = 2L)
    } else {
        exon <- matrix(c(1L, j, 200L, 199L + nchar(cDNA) - j), ncol = 2L, byrow = TRUE)
        genome <- paste0(substr(cDNA, 1L, j), strrep('G', 199L - j),
                         substr(cDNA, j + 1L, nchar(cDNA)))
    }
    fasta <- c(fasta, paste0('>', name), genome)
    end <- nchar(genome)
    gene <- paste0('gene:G_', name); tx <- paste0('transcript:T_', name)
    gff <- c(gff, paste(name, 'fixture', 'gene', 1, end, '.', '+', '.', paste0('ID=', gene, ';biotype=protein_coding'), sep = '\t'),
             paste(name, 'fixture', 'mRNA', 1, end, '.', '+', '.', paste0('ID=', tx, ';Parent=', gene, ';biotype=protein_coding'), sep = '\t'))
    for (k in seq_len(nrow(exon))) {
        phase <- if (k == 1L) 0L else (3L - j %% 3L) %% 3L
        for (feature in c('exon', 'CDS'))
            gff <- c(gff, paste(name, 'fixture', feature, exon[k, 1L], exon[k, 2L], '.', '+',
                                if (feature == 'CDS') phase else '.', paste0('Parent=', tx), sep = '\t'))
    }
    for (i in seq_along(item$edits)) {
        e <- item$edits[[i]]; p <- as.integer(e['pos'])
        stopifnot(substr(cDNA, p, p + nchar(e['ref']) - 1L) == e['ref'])
        genomic <- if (is.na(j) || p <= j) p else p - j + 199L
        lane <- if (length(item$lane) == 1L) item$lane else item$lane[i]
        records <- c(records, paste(name, genomic, '.', e['ref'], e['alt'], '.', 'PASS', '.', 'GT',
                                    if (lane == 0L) '1|0' else '0|1', sep = '\t'))
    }
    for (lane in 0:1) {
        assigned <- if (length(item$lane) == 1L) rep(item$lane, length(item$edits)) else as.list(item$lane)
        if (!any(vapply(assigned, function(x) lane %in% x, logical(1L)))) {
            gold[[length(gold)+1L]] <- data.frame(case = name, lane = lane, so_terms = '',
              nmd = 'not_applicable', first_stop_end = 120L, junction = j)
            next
        }
        edited <- edit_lane(cDNA, item$edits, assigned, lane)
        orig <- translate(cDNA); protein <- translate(edited)
        first <- match('*', protein); refstop <- match('*', orig)
        startlost <- substr(edited, 1L, 3L) != 'ATG'
        gained <- !is.na(first) && first < refstop
        lost <- is.na(first) && !startlost
        retained <- !gained && !lost && substr(edited, nchar(edited)-2L, nchar(edited)) !=
          substr(cDNA, nchar(cDNA)-2L, nchar(cDNA)) && !startlost
        # A displaced frame is pathogenic if a stop occurs before restoration
        # or its offset is nonzero at coding-sequence exhaustion.
        offsets <- cumsum(vapply(item$edits, function(e) nchar(e['alt']) - nchar(e['ref']), integer(1L)))
        edited_stops <- vapply(item$edits, function(e) as.integer(e['pos']) + nchar(e['alt']) - 1L, integer(1L)) +
          c(0L, head(offsets, -1L))
        frame_at_stop <- !is.na(first) && any(offsets %% 3L != 0L & edited_stops <= first * 3L &
                           c(edited_stops[-1L], nchar(edited)) >= first * 3L)
        frame_open <- any(offsets %% 3L != 0L) && (frame_at_stop || tail(offsets, 1L) %% 3L != 0L)
        terms <- if (startlost) 'start_lost' else c(if (frame_open) 'frameshift_variant',
          if (gained) 'stop_gained', if (lost) 'stop_lost')
        if (!length(terms) && retained) terms <- 'stop_retained_variant'
        if (!length(terms)) terms <- if (identical(protein, orig)) 'synonymous_variant' else
            if (any(offsets %% 3L != 0L) || any(offsets != 0L)) 'protein_altering_variant' else 'missense_variant'
        junction <- if (is.na(j)) NA_integer_ else j + sum(vapply(seq_along(item$edits), function(i)
          if (lane %in% assigned[[i]] && as.integer(item$edits[[i]]['pos']) <= j)
            nchar(item$edits[[i]]['alt']) - nchar(item$edits[[i]]['ref']) else 0L, integer(1L)))
        if (startlost || lost) nmd <- 'unknown' else if (!gained) nmd <- 'not_applicable' else if (is.na(j)) nmd <- 'escape' else
            nmd <- if (junction - 3L * first > 50L) 'trigger' else 'escape'
        gold[[length(gold)+1L]] <- data.frame(case = name, lane = lane, so_terms = paste(terms, collapse = ','),
                                              nmd = nmd, first_stop_end = if (is.na(first)) NA_integer_ else first*3L,
                                              junction = junction)
    }
}
expected <- c(cis = 'synonymous_variant', trans = 'missense_variant',
  frame_open = 'frameshift_variant,stop_gained', restored_before_stop = 'protein_altering_variant',
  restored_after_stop = 'frameshift_variant,stop_gained', start_lost = 'start_lost',
  stop_lost = 'stop_lost', stop_retained = 'stop_retained_variant', nmd49 = 'stop_gained',
  nmd50 = 'stop_gained', nmd51 = 'stop_gained', intronless = 'stop_gained',
  shifted_junction = 'stop_gained')
golden <- do.call(rbind, gold)
stopifnot(identical(unname(golden$so_terms[golden$lane == 0L]), unname(expected)),
          identical(golden$nmd[golden$case %in% c('nmd49', 'nmd50', 'nmd51') & golden$lane == 0L],
                    c('escape', 'escape', 'trigger')),
          golden$nmd[golden$case == 'shifted_junction' & golden$lane == 0L] == 'trigger',
          all(golden$so_terms[golden$lane == 1L & golden$case != 'trans'] == ''))
verified_write <- function(lines, path) {
    lines <- unname(lines)
    if (file.exists(path) && !identical(readLines(path, warn = FALSE), lines))
        stop('fixture differs from independent generator: ', path)
    writeLines(lines, path)
}
verified_write(fasta, file.path(out, 'vertical.fa'))
verified_write(gff, file.path(out, 'vertical.gff3'))
verified_write(c('##fileformat=VCFv4.2', vapply(names(cases), function(n)
  paste0('##contig=<ID=', n, ',length=', nchar(fasta[match(paste0('>', n), fasta)+1L]), '>'), character(1L)),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS', records), file.path(out, 'vertical.vcf'))
golden_text <- capture.output(write.table(golden, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
verified_write(golden_text, file.path(out, 'vertical_goldens.tsv'))
