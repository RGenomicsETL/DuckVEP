#!/usr/bin/env Rscript
# Start and stop (coding-v1 slice 5) fixtures and goldens. Base R only: no DuckVEP code, no csq. The
# 36-codon coding transcript of the frame fixture (its +1 and -2 displaced frames read stops at cDNA 27 and
# 102, its -1 and +2 frames read none before the CDS runs out) is laid out on both strands in one-, two- and
# three-exon arrangements (introns inside codons 12 and 18 and at the codon 25 boundary). Every scenario is
# written in transcript (cDNA) coordinates, mapped to genome VCF records, and its expected consequence is
# computed by editing the cDNA, translating it with the standard code, and evaluating the reading-frame
# history base by base on the edited sequence and where the reference terminator (cDNA 106-108) went.
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
    if (nchar(s) < 3L) return(character())
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
snv <- function(p, alt, lane = 0L) { stopifnot(b(p) != alt); v(p, b(p), alt, lane) }
# Every scenario is decided by the start/stop classifier or shares its rules (slice 5). expect_a/expect_b are
# documentation only: the policy below computes every golden.
scenarios <- list(
  # -- start codon
  start_snv_first = list(slice = 5L, edits = list(snv(1, 'C'))),
  start_snv_last = list(slice = 5L, edits = list(snv(3, 'A'))),
  start_mnv = list(slice = 5L, edits = list(v(1, 'ATG', 'CCC'))),
  start_deletion_frame = list(slice = 5L, edits = list(del(1, 2))),
  start_deletion_inframe = list(slice = 5L, edits = list(del(1, 3))),
  start_insertion_inside = list(slice = 5L, edits = list(ins(2, 'CC'))),
  start_insertion_after_start = list(slice = 5L, edits = list(ins(3, 'GAA'))),
  start_rescued_by_insertion_a = list(slice = 5L, edits = list(ins(1, 'TGA'))),
  start_rescued_by_insertion_b = list(slice = 5L, edits = list(ins(2, 'GAT'))),
  start_deletion_within_exon = list(slice = 5L, edits = list(del(1, 30))),
  start_lost_suppresses_missense = list(slice = 5L, edits = list(snv(1, 'C'), snv(50, 'A'))),
  start_lost_suppresses_frame = list(slice = 5L, edits = list(snv(2, 'C'), ins(30, 'T'))),
  start_lost_suppresses_stop_gain = list(slice = 5L, edits = list(snv(1, 'T'), snv(40, 'T'))),
  start_lost_suppresses_stop_lost = list(slice = 5L, edits = list(snv(1, 'C'), snv(106, 'C'))),
  start_lost_trans_control = list(slice = 5L, edits = list(snv(1, 'C', 0L), snv(50, 'A', 1L))),
  start_lost_hom_alt = list(slice = 5L, edits = list(snv(1, 'C', 0:1))),
  start_intact_frame_open_early_stop = list(slice = 5L, edits = list(ins(3, 'G'))),
  # -- stop lost and run-off frames
  stop_lost_snv_first = list(slice = 5L, edits = list(snv(106, 'C'))),
  stop_lost_snv_second = list(slice = 5L, edits = list(snv(107, 'C'))),
  stop_lost_snv_third = list(slice = 5L, edits = list(snv(108, 'C'))),
  stop_lost_mnv = list(slice = 5L, edits = list(v(106, 'TAA', 'TGG'))),
  stop_lost_terminator_deleted = list(slice = 5L, edits = list(del(105, 3))),
  stop_lost_last_codon_and_terminator_deleted = list(slice = 5L, edits = list(del(102, 6))),
  stop_lost_deletion_partial = list(slice = 5L, edits = list(del(105, 1))),
  stop_lost_insertion_inside = list(slice = 5L, edits = list(ins(106, 'CCC'))),
  stop_lost_insertion_before_plus1 = list(slice = 5L, edits = list(ins(105, 'G'))),
  stop_lost_insertion_before_plus2 = list(slice = 5L, edits = list(ins(105, 'GG'))),
  run_off_deletion = list(slice = 5L, edits = list(del(50, 1))),
  run_off_insertion_plus2 = list(slice = 5L, edits = list(ins(40, 'AC'))),
  run_off_inside_codon_35 = list(slice = 5L, edits = list(ins(103, 'G'))),
  run_off_trans_control = list(slice = 5L, edits = list(del(50, 1, 0L), snv(40, 'A', 1L))),
  stop_lost_hom_alt = list(slice = 5L, edits = list(snv(106, 'C', 0:1))),
  stop_lost_and_missense = list(slice = 5L, edits = list(snv(106, 'C'), snv(50, 'A'))),
  # -- stop retained
  stop_retained_tag = list(slice = 5L, edits = list(snv(108, 'G'))),
  stop_retained_tga = list(slice = 5L, edits = list(snv(107, 'G'))),
  stop_retained_plus_synonymous = list(slice = 5L, edits = list(snv(107, 'G'), snv(9, 'A'))),
  stop_retained_plus_missense = list(slice = 5L, edits = list(snv(107, 'G'), snv(50, 'A'))),
  stop_retained_plus_inframe_insertion = list(slice = 5L, edits = list(snv(108, 'G'), ins(60, 'GGG'))),
  stop_retained_trans_control = list(slice = 5L, edits = list(snv(107, 'G', 0L), snv(50, 'A', 1L))),
  stop_retained_hom_alt = list(slice = 5L, edits = list(snv(107, 'G', 0:1))),
  stop_retained_after_frame_pair = list(slice = 5L, edits = list(snv(108, 'G'), ins(30, 'T'), del(60, 1))),
  stop_retained_displaced_insertion = list(slice = 5L, edits = list(ins(106, 'G'))),
  # -- stops near the terminator and controls that keep it intact
  inframe_insertion_before_terminator = list(slice = 5L, edits = list(ins(105, 'GGG'))),
  stop_inserted_before_terminator = list(slice = 5L, edits = list(ins(105, 'TAA'))),
  stop_gained_in_codon_35 = list(slice = 5L, edits = list(snv(103, 'T'), snv(105, 'A'))),
  inframe_deletion_before_terminator = list(slice = 5L, edits = list(del(102, 3))),
  # -- rescued terminations
  rescued_termination_by_restoring_insertion = list(slice = 5L, edits = list(del(50, 1), ins(104, 'C'))),
  rescued_termination_restore_inside_terminator = list(slice = 5L, edits = list(del(50, 1), ins(106, 'C'))),
  rescued_termination_by_deletion_in_terminator = list(slice = 5L, edits = list(ins(50, 'A'), del(106, 1))),
  rescued_termination_by_terminator_deletion = list(slice = 5L, edits = list(ins(50, 'A'), del(105, 1))),
  rescued_termination_after_stop_lost = list(slice = 5L, edits = list(snv(106, 'C'), del(50, 1), ins(70, 'C'))),
  # -- edits after the first stop are contributors, not expressed effects
  post_stop_terminator_lost = list(slice = 5L, edits = list(snv(40, 'T'), snv(106, 'C'))),
  post_stop_terminator_retained = list(slice = 5L, edits = list(snv(40, 'T'), snv(107, 'G'))),
  post_stop_run_off = list(slice = 5L, edits = list(ins(10, 'G'), del(105, 1))),
  post_stop_two_edits = list(slice = 5L, edits = list(snv(40, 'T'), ins(70, 'C'), snv(108, 'C'))),
  post_stop_restoring_deletion = list(slice = 5L, edits = list(ins(10, 'G'), del(60, 1))),
  restoring_insertion_creates_stop = list(slice = 5L, edits = list(del(30, 2), ins(60, 'GA'))),
  post_stop_after_start_lost = list(slice = 5L, edits = list(snv(1, 'C'), snv(40, 'T'), snv(106, 'C'))))

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
    # in_window flags the edited bases that hold the reference terminator: retained bases of the last
    # three reference positions and every base of an island that reaches them (an insertion inside them
    # counts, one just before them does not).
    in_window <- logical(); term0 <- nchar(cds) - 3L
    retained <- function(upto) {
        if (upto >= at) {
            seq_out <<- c(seq_out, strsplit(substr(cds, at, upto), '')[[1L]])
            displaced <<- c(displaced, rep(offset %% 3L != 0L, upto - at + 1L))
            in_window <<- c(in_window, at:upto > term0)
        }
    }
    for (e in isl) {
        retained(e$start - 1L)
        d <- nchar(e$alt) - e$ref_len
        if (nchar(e$alt)) {
            seq_out <- c(seq_out, strsplit(e$alt, '')[[1L]])
            displaced <- c(displaced, rep(offset %% 3L != 0L || (offset + d) %% 3L != 0L, nchar(e$alt)))
            in_window <- c(in_window, rep(e$start - 1L + e$ref_len > term0 && e$start - 1L < nchar(cds), nchar(e$alt)))
        }
        at <- e$start + e$ref_len; offset <- offset + d
    }
    retained(nchar(cds))
    list(alt = paste(seq_out, collapse = ''), displaced = displaced,
        window = if (any(in_window)) which(in_window)[1L] else length(seq_out) + 1L, nominal = vapply(isl,
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
# Returns the reduced whole-protein term set of the edited cDNA (contract section 2), in this order:
# the edited CDS not beginning with ATG is start_lost alone; no stop is stop_lost (plus frameshift when
# the length is not a multiple of three, i.e. the frame is displaced when the CDS runs out; nothing is
# ever extended downstream); a first stop starting before the terminator's window that truncates the reference
# peptide is stop_gained (an inserted stop codon after the unchanged peptide is judged as the terminal codon; plus
# frameshift when one of its bases is displaced); a first stop at the window is a frameshift when
# displaced, else identical peptide is stop_retained (terminal codon changed or moved) or synonymous
# (unchanged), and a changed peptide is protein_altering (any nominal frame edit), missense, inframe
# insertion/deletion or protein_altering.
classify <- function(edits, lane) {
    if (!any(vapply(edits, function(e) lane %in% e$lane, logical(1L)))) return(list(terms = character(), protein = ''))
    w <- walk(edits, lane); alt <- w$alt; L <- nchar(alt)
    protein <- translate(alt); first <- match('*', protein)
    visible <- if (is.na(first)) protein else protein[seq_len(first)]
    text <- paste(visible, collapse = '')
    if (substr(alt, 1L, 3L) != 'ATG') return(list(terms = 'start_lost', protein = text))
    terms <- if (is.na(first)) {
        c('stop_lost', if (L %% 3L != 0L) 'frameshift_variant')
    } else {
        last <- 3L * first
        stop_bases <- (last - 2L):last
        a <- protein[seq_len(first - 1L)]; r <- ref_protein[-length(ref_protein)]
        if (last - 2L < w$window && !identical(a, r)) c('stop_gained', if (any(w$displaced[stop_bases])) 'frameshift_variant')
        else if (any(w$displaced[stop_bases])) 'frameshift_variant' else {
            if (identical(a, r)) {
                if (last - 2L != w$window || substr(alt, last - 2L, last) != substring(cds, nchar(cds) - 2L))
                    'stop_retained_variant' else 'synonymous_variant'
            } else if (any(w$nominal)) 'protein_altering_variant' else
                if (all(w$kinds == 'sub')) 'missense_variant' else if (all(w$kinds == 'ins')) 'inframe_insertion' else
                if (all(w$kinds == 'del')) 'inframe_deletion' else 'protein_altering_variant'
        }
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
        if (!is.null(scenario$only) && !(layout_name %in% scenario$only)) next
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
verified_write(fasta, file.path(out, 'startstop.fa'))
offset <- 0L; fai <- character()
for (i in seq_along(contigs)) {
    offset <- offset + nchar(contigs[i]) + 2L
    fai <- c(fai, paste(contigs[i], lengths[i], offset, lengths[i], lengths[i] + 1L, sep = '\t'))
    offset <- offset + lengths[i] + 1L
}
verified_write(fai, file.path(out, 'startstop.fa.fai'))
verified_write(gff, file.path(out, 'startstop.gff3'))
vcf <- records[order(match(records$chrom, contigs), records$pos, records$order), ]
verified_write(c('##fileformat=VCFv4.2',
  vapply(seq_along(contigs), function(i) paste0('##contig=<ID=', contigs[i], ',length=', lengths[i], '>'), character(1L)),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS',
  paste(vcf$chrom, vcf$pos, '.', vcf$ref, vcf$alt, '.', 'PASS', '.', 'GT', vcf$gt, sep = '\t')),
  file.path(out, 'startstop.vcf'))
verified_write(tsv(golden), file.path(out, 'startstop_goldens.tsv'))
cat(length(contigs), 'contigs;', nrow(vcf), 'VCF records;', nrow(golden), 'lane goldens\n')
if (nzchar(Sys.getenv('STARTSTOP_SUMMARY'))) {
    one <- golden[golden$layout == 'plus_1exon' & golden$so_terms != '' , c('scenario', 'lane', 'so_terms', 'protein')]
    print(one, right = FALSE, row.names = FALSE)
}
