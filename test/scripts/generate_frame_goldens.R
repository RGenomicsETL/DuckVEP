#!/usr/bin/env Rscript
# Frame opening/restoration and stop-gain (coding-v1 slice 4) fixtures and goldens. Base R only: no
# DuckVEP code, no csq. One 36-codon coding transcript is laid out on both strands and in one-, two- and
# three-exon genomic arrangements (introns inside codons 12 and 18, and at the codon 25 boundary). Every
# scenario is written in transcript (cDNA) coordinates, mapped to genome VCF records, and its expected
# consequence is computed by editing the cDNA, translating it with the standard code, and evaluating the
# reading-frame history base by base on the edited sequence (not from any DuckVEP structure).
#
# The transcript is designed so the +1 (and -2) displaced frame has stops at cDNA 27 and 102, while the
# -1 (and +2) displaced frame has no stop before the CDS runs out: the same nominal +1/-1 history can
# be rescued or terminated early depending on where the edits sit and what they insert.
out <- 'test/data/haplotype'
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cds <- 'ATGAGTCCCATCCGAAACTTCCTCCGTGATTCGGAATTACAAAAACATGTCTCGGGCCGTTGCGCCCCAAAAGCTATCCTCTCCCTACATGGCCCGTTGTGTAATTAA'
stopifnot(nchar(cds) == 108L)
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
revcomp <- function(s) {
    chars <- strsplit(s, '')[[1L]]
    paste(rev(c(A = 'T', C = 'G', G = 'C', T = 'A')[chars]), collapse = '')
}
sub_at <- function(s, p, n) substr(s, p, p + n - 1L)
ref_protein <- translate(cds)
stopifnot(sum(ref_protein == '*') == 1L, ref_protein[36L] == '*')

# ---- layouts: exon lengths in transcript order -------------------------------------------------
layouts <- list(
  plus_1exon = list(strand = '+', exons = 108L),
  plus_2exon = list(strand = '+', exons = c(52L, 56L)),        # codon 18 (52-54) straddles the intron
  plus_3exon = list(strand = '+', exons = c(35L, 40L, 33L)),   # codon 12 (34-36) straddles; codon 25 boundary
  minus_1exon = list(strand = '-', exons = 108L),
  minus_2exon = list(strand = '-', exons = c(52L, 56L)),
  minus_3exon = list(strand = '-', exons = c(35L, 40L, 33L)))
intron <- 'GTAAGTCTGATCGAGTTCGACTTGACGTTAAGCTTACCAG'
flank <- 'CTAGGATCCTAGCATG'

# ---- scenarios: edits are (pos, ref, alt) in transcript orientation, VCF-anchored for indels --------
v <- function(pos, ref, alt, lane = 0L) list(pos = pos, ref = ref, alt = alt, lane = lane)
b <- function(p) sub_at(cds, p, 1L)
ins <- function(p, bases, lane = 0L) v(p, b(p), paste0(b(p), bases), lane)   # insert after base p
del <- function(p, n, lane = 0L) v(p, sub_at(cds, p, n + 1L), b(p), lane)     # delete bases p+1 .. p+n
snv <- function(p, alt, lane = 0L) v(p, b(p), alt, lane)
# slice 4 = decided by the frame/stop-gain classifier; slice 5 = start/terminal effects (pending).
# expect_a/expect_b are the intended lane results (documentation checked against the policy below).
scenarios <- list(
  frame_open_early_stop = list(slice = 4L, edits = list(ins(10, 'G'))),
  frame_open_runs_off = list(slice = 4L, edits = list(del(50, 1))),
  frame_open_late_stop = list(slice = 4L, edits = list(ins(30, 'T'))),
  frame_open_plus2_runs_off = list(slice = 4L, edits = list(ins(40, 'AC'))),
  frame_open_minus2_late_stop = list(slice = 4L, edits = list(del(40, 2))),
  frame_open_before_terminal = list(slice = 4L, edits = list(ins(103, 'G'))),
  restored_rescued = list(slice = 4L, edits = list(ins(30, 'T'), del(60, 1))),
  restored_early_stop = list(slice = 4L, edits = list(ins(10, 'G'), del(60, 1))),
  restored_rescued_wide = list(slice = 4L, edits = list(ins(30, 'T'), del(90, 1))),
  restored_early_stop_wide = list(slice = 4L, edits = list(ins(10, 'G'), del(90, 1))),
  restored_near_introns = list(slice = 4L, edits = list(ins(34, 'T'), del(53, 1))),
  restored_minus_plus = list(slice = 4L, edits = list(del(20, 1), ins(70, 'C'))),
  restored_plus2_minus2 = list(slice = 4L, edits = list(ins(31, 'AC'), del(55, 2))),
  restored_minus2_plus2 = list(slice = 4L, edits = list(del(30, 2), ins(60, 'CG'))),
  restoring_insertion_creates_stop = list(slice = 4L, edits = list(del(30, 2), ins(60, 'GA'))),
  restored_identical_peptide_early = list(slice = 4L, edits = list(del(17, 1), ins(20, 'T'))),
  restored_identical_peptide_late = list(slice = 4L, edits = list(del(82, 1), ins(86, 'T'))),
  restored_two_pairs = list(slice = 4L, edits = list(ins(30, 'T'), del(45, 1), ins(60, 'A'), del(78, 1))),
  restored_three_insertions = list(slice = 4L, edits = list(ins(30, 'T'), ins(50, 'A'), ins(70, 'C'))),
  restored_ins1_ins2 = list(slice = 4L, edits = list(ins(30, 'T'), ins(60, 'AC'))),
  restored_del1_del2 = list(slice = 4L, edits = list(del(30, 1), del(60, 2))),
  stop_gained_snv = list(slice = 4L, edits = list(snv(40, 'T'))),
  stop_gained_cis_pair = list(slice = 4L, edits = list(snv(91, 'T'), snv(93, 'A'))),
  stop_gained_inframe_insertion = list(slice = 4L, edits = list(ins(60, 'TAA'))),
  stop_gained_post_stop_frame = list(slice = 4L, edits = list(snv(40, 'T'), ins(80, 'G'))),
  stop_gained_after_restoration = list(slice = 4L, edits = list(ins(30, 'T'), del(60, 1), snv(98, 'A'))),
  early_stop_then_missense = list(slice = 4L, edits = list(ins(10, 'G'), snv(64, 'T'))),
  trans_control_of_pair = list(slice = 4L, edits = list(ins(30, 'T', 0L), del(60, 1, 1L))),
  hom_alt_frame_open = list(slice = 4L, edits = list(ins(30, 'T', 0:1))),
  start_codon_insertion = list(slice = 5L, edits = list(ins(2, 'C'))))

# ---- independent policy ------------------------------------------------------------------------
# Islands: strip the VCF anchor (common prefix) so an edit is only its differing bases.
island <- function(e) {
    r <- e$ref; a <- e$alt; k <- 0L
    while (nchar(r) > 0L && nchar(a) > 0L && substr(r, 1L, 1L) == substr(a, 1L, 1L) && nchar(r) != nchar(a)) {
        r <- substring(r, 2L); a <- substring(a, 2L); k <- k + 1L
    }
    list(start = e$pos + k, ref_len = nchar(r), alt = a)
}
# Walk the ascending islands over the reference, building the edited sequence and, for every edited
# base, whether the reading frame is displaced there. A retained base is displaced when the cumulative
# length change is not a multiple of three; an inserted/replacement base is displaced when either the
# entering or the leaving frame is displaced.
walk <- function(edits, lane) {
    mine <- Filter(function(e) lane %in% e$lane, edits)
    isl <- lapply(mine, island)
    isl <- isl[order(vapply(isl, function(x) x$start, numeric(1L)))]
    seq_out <- character(); displaced <- logical(); at <- 1L; offset <- 0L
    retained <- function(upto) {
        if (upto >= at) {
            seq_out <<- c(seq_out, strsplit(substr(cds, at, upto), '')[[1L]])
            displaced <<- c(displaced, rep(offset %% 3L != 0L, upto - at + 1L))
        }
    }
    for (e in isl) {
        retained(e$start - 1L)
        d <- nchar(e$alt) - e$ref_len
        if (nchar(e$alt)) {
            seq_out <- c(seq_out, strsplit(e$alt, '')[[1L]])
            displaced <- c(displaced, rep(offset %% 3L != 0L || (offset + d) %% 3L != 0L, nchar(e$alt)))
        }
        at <- e$start + e$ref_len; offset <- offset + d
    }
    retained(nchar(cds))
    list(alt = paste(seq_out, collapse = ''), displaced = displaced, nominal = vapply(isl,
        function(x) (nchar(x$alt) - x$ref_len) %% 3L != 0L, logical(1L)), n = length(isl),
        kinds = vapply(isl, function(x) if (x$ref_len == 0L) 'ins' else if (!nchar(x$alt)) 'del' else
            if (x$ref_len == nchar(x$alt)) 'sub' else 'mixed', character(1L)))
}
impact_of <- function(terms) {
    if (!length(terms)) return('')
    levels <- c(stop_gained = 'HIGH', start_lost = 'HIGH', stop_lost = 'HIGH', frameshift_variant = 'HIGH',
                missense_variant = 'MODERATE', inframe_insertion = 'MODERATE', inframe_deletion = 'MODERATE',
                protein_altering_variant = 'MODERATE', synonymous_variant = 'LOW', stop_retained_variant = 'LOW')
    order <- c('HIGH', 'MODERATE', 'LOW')
    order[min(match(levels[terms], order))]
}
# Returns the reduced whole-protein term set, or NA when the lane is outside this slice (start/terminal
# codon effects, or a first stop overlapping the reference terminator).
classify <- function(edits, lane) {
    if (!any(vapply(edits, function(e) lane %in% e$lane, logical(1L)))) return(list(terms = character(), protein = ''))
    w <- walk(edits, lane); alt <- w$alt; L <- nchar(alt)
    protein <- translate(alt); first <- match('*', protein)
    visible <- if (is.na(first)) protein else protein[seq_len(first)]
    text <- paste(visible, collapse = '')
    if (substr(alt, 1L, 3L) != 'ATG') return(list(terms = 'start_lost', protein = text))
    terms <- if (is.na(first)) {
        if (L %% 3L != 0L) 'frameshift_variant' else NA_character_
    } else {
        last <- 3L * first
        stop_bases <- (last - 2L):last
        if (last == L && L %% 3L == 0L) { # the reference terminator, read in frame
            if (any(w$displaced[stop_bases])) 'frameshift_variant' else {
                a <- protein[-length(protein)]; r <- ref_protein[-length(ref_protein)]
                if (identical(a, r)) 'synonymous_variant' else if (any(w$nominal)) 'protein_altering_variant' else
                    if (all(w$kinds == 'sub')) 'missense_variant' else if (all(w$kinds == 'ins')) 'inframe_insertion' else
                    if (all(w$kinds == 'del')) 'inframe_deletion' else 'protein_altering_variant'
            }
        } else if (last <= L - 3L) c('stop_gained', if (any(w$displaced[stop_bases])) 'frameshift_variant')
        else NA_character_
    }
    list(terms = terms, protein = text)
}

# ---- genome layout and VCF mapping -------------------------------------------------------------
records <- list(); gff <- '##gff-version 3'; fasta <- character(); gold <- list()
contigs <- character(); lengths <- integer()
for (layout_name in names(layouts)) {
    layout <- layouts[[layout_name]]; minus <- layout$strand == '-'
    ends <- cumsum(layout$exons); starts <- c(1L, head(ends, -1L) + 1L)
    pieces <- character()
    for (k in seq_along(ends)) pieces <- c(pieces, substr(cds, starts[k], ends[k]))
    transcript_genome <- paste0(flank, paste(pieces, collapse = intron), flank)
    idx <- function(t) { k <- findInterval(t - 1L, ends) + 1L; nchar(flank) + t + (k - 1L) * nchar(intron) }
    n <- nchar(transcript_genome)
    genome <- if (minus) revcomp(transcript_genome) else transcript_genome
    G <- function(t) if (minus) n + 1L - idx(t) else idx(t)
    for (scenario_name in names(scenarios)) {
        scenario <- scenarios[[scenario_name]]; edits <- scenario$edits
        name <- paste(layout_name, scenario_name, sep = '__')
        contigs <- c(contigs, name); lengths <- c(lengths, n); fasta <- c(fasta, paste0('>', name), genome)
        lo <- min(G(1L), G(108L)); hi <- max(G(1L), G(108L))
        gene <- paste0('gene:G_', name); tx <- paste0('transcript:T_', name)
        gff <- c(gff,
          paste(name, 'fixture', 'gene', lo, hi, '.', layout$strand, '.', paste0('ID=', gene, ';biotype=protein_coding'), sep = '\t'),
          paste(name, 'fixture', 'mRNA', lo, hi, '.', layout$strand, '.', paste0('ID=', tx, ';Parent=', gene, ';biotype=protein_coding'), sep = '\t'))
        for (k in seq_along(ends)) {
            e1 <- G(starts[k]); e2 <- G(ends[k]); a <- min(e1, e2); bb <- max(e1, e2)
            phase <- (3L - (starts[k] - 1L) %% 3L) %% 3L
            for (feature in c('exon', 'CDS'))
                gff <- c(gff, paste(name, 'fixture', feature, a, bb, '.', layout$strand,
                                    if (feature == 'CDS') phase else '.', paste0('Parent=', tx), sep = '\t'))
        }
        for (e in edits) {
            r <- nchar(e$ref); al <- nchar(e$alt)
            anchored <- r != al && (startsWith(e$alt, e$ref) || startsWith(e$ref, e$alt))
            if (!minus) {
                pos <- G(e$pos); ref <- e$ref; alt <- e$alt
            } else if (!anchored) {
                pos <- G(e$pos + r - 1L); ref <- revcomp(e$ref); alt <- revcomp(e$alt)
            } else if (r > al) { # deletion of transcript bases pos+al .. pos+r-1
                del_first <- e$pos + al; del_last <- e$pos + r - 1L
                low <- G(del_last); high <- G(del_first)
                pos <- low - 1L; ref <- substr(genome, low - 1L, high); alt <- substr(genome, low - 1L, low - 1L)
            } else { # insertion after transcript base pos
                pos <- G(e$pos + 1L)
                ref <- substr(genome, pos, pos)
                alt <- paste0(ref, revcomp(substring(e$alt, r + 1L)))
            }
            gt <- if (length(e$lane) == 2L) '1|1' else if (e$lane == 0L) '1|0' else '0|1'
            records[[length(records) + 1L]] <- data.frame(chrom = name, pos = pos, ref = ref, alt = alt, gt = gt,
                order = length(records))
        }
        for (lane in 0:1) {
            res <- classify(edits, lane)
            stopifnot(!anyNA(res$terms))
            gold[[length(gold) + 1L]] <- data.frame(case = name, layout = layout_name, scenario = scenario_name,
                lane = lane, so_terms = paste(res$terms, collapse = ','), impact = impact_of(res$terms),
                protein = res$protein, classifier_slice = scenario$slice)
        }
    }
}
records <- do.call(rbind, records)
golden <- do.call(rbind, gold)
stopifnot(!anyDuplicated(contigs), nrow(golden) == 2L * length(contigs))

verified_write <- function(lines, path) {
    lines <- unname(lines)
    if (file.exists(path) && !identical(readLines(path, warn = FALSE), lines))
        stop('fixture differs from independent generator: ', path)
    writeLines(lines, path)
}
tsv <- function(frame) capture.output(write.table(frame, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
verified_write(fasta, file.path(out, 'frame.fa'))
offset <- 0L; fai <- character()
for (i in seq_along(contigs)) {
    offset <- offset + nchar(contigs[i]) + 2L
    fai <- c(fai, paste(contigs[i], lengths[i], offset, lengths[i], lengths[i] + 1L, sep = '\t'))
    offset <- offset + lengths[i] + 1L
}
verified_write(fai, file.path(out, 'frame.fa.fai'))
verified_write(gff, file.path(out, 'frame.gff3'))
vcf <- records[order(match(records$chrom, contigs), records$pos, records$order), ]
verified_write(c('##fileformat=VCFv4.2',
  vapply(seq_along(contigs), function(i) paste0('##contig=<ID=', contigs[i], ',length=', lengths[i], '>'), character(1L)),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS',
  paste(vcf$chrom, vcf$pos, '.', vcf$ref, vcf$alt, '.', 'PASS', '.', 'GT', vcf$gt, sep = '\t')),
  file.path(out, 'frame.vcf'))
verified_write(tsv(golden), file.path(out, 'frame_goldens.tsv'))
cat(length(contigs), 'contigs;', nrow(vcf), 'VCF records;', nrow(golden), 'lane goldens\n')
if (nzchar(Sys.getenv('FRAME_SUMMARY'))) {
    one <- golden[golden$layout == 'plus_1exon' & golden$so_terms != '' , c('scenario', 'lane', 'so_terms', 'protein')]
    print(one, right = FALSE, row.names = FALSE)
}
