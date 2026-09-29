#!/usr/bin/env Rscript
# Same-codon (coding-v1 slice 3) fixtures and goldens. Base R only: no DuckVEP code, no csq.
# One 24-codon coding transcript is laid out on both strands and in one-, two- and three-exon
# genomic arrangements (codons that straddle an intron included). Every scenario is written in
# transcript (cDNA) coordinates, mapped to genome VCF records, and its expected consequence is
# computed by editing the cDNA and translating it with the standard code.
out <- 'test/data/haplotype'
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cds <- paste(c('ATG', 'GCT', 'GAA', 'CTG', 'AAA', 'CGC', 'TCA', 'GGT', 'CCA', 'TTC', 'GAT', 'CAC',
               'TGG', 'ACC', 'GTA', 'AGC', 'ATT', 'CAG', 'GAG', 'TAC', 'TGT', 'CTC', 'GCC', 'TAA'),
             collapse = '')
stopifnot(nchar(cds) == 72L)

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

# ---- layouts: exon lengths in transcript order -------------------------------------------------
layouts <- list(
  plus_1exon = list(strand = '+', exons = 72L),
  plus_2exon = list(strand = '+', exons = c(46L, 26L)),  # codon 16 (46-48) straddles the intron
  plus_3exon = list(strand = '+', exons = c(31L, 31L, 10L)), # codons 11 (31-33) and 21 (61-63) straddle
  minus_1exon = list(strand = '-', exons = 72L),
  minus_2exon = list(strand = '-', exons = c(46L, 26L)),
  minus_3exon = list(strand = '-', exons = c(31L, 31L, 10L)))
intron <- 'GTAAGTCTGATCGAGTTCGACTTGACGTTAAGCTTACCAG'
flank <- 'CTAGGATCCTAGCATG'

# ---- scenarios: edits are (pos, ref, alt) in transcript orientation; lanes 0/1 --------------------
v <- function(pos, ref, alt, lane = 0L) list(pos = pos, ref = ref, alt = alt, lane = lane)
# slice 3 = decided by the same-codon classifier; slice 4 = a new stop before the terminator (decided by the
# frame/stop-gain classifier); slice 5 = start/terminal effects (pending until then).
scenarios <- list(
  cis_syn_pair_only_together = list(slice = 3L, edits = list(v(46, 'A', 'T'), v(47, 'G', 'C'))),
  trans_control_of_syn_pair = list(slice = 3L, edits = list(v(46, 'A', 'T', 0L), v(47, 'G', 'C', 1L))),
  cis_mis_pair = list(slice = 3L, edits = list(v(4, 'G', 'A'), v(5, 'C', 'A'))),
  cis_syn_pair_r = list(slice = 3L, edits = list(v(16, 'C', 'A'), v(18, 'C', 'A'))),
  cis_mnv_record = list(slice = 3L, edits = list(v(34, 'CA', 'GT'))),
  cis_three_snv_codon = list(slice = 3L, edits = list(v(37, 'T', 'C'), v(38, 'G', 'C'), v(39, 'G', 'C'))),
  hom_alt_snv = list(slice = 3L, edits = list(v(5, 'C', 'A', 0:1))),
  ins_at_codon_boundary = list(slice = 3L, edits = list(v(12, 'G', 'GTCT'))),
  ins_inside_codon_pure = list(slice = 3L, edits = list(v(4, 'G', 'GCTG'))),
  ins_inside_codon_altering = list(slice = 3L, edits = list(v(4, 'G', 'GAAA'))),
  del_whole_codon = list(slice = 3L, edits = list(v(12, 'GAAA', 'G'))),
  del_inside_codons_pure = list(slice = 3L, edits = list(v(52, 'CAGG', 'C'))),
  delins_inframe = list(slice = 3L, edits = list(v(7, 'GAACTG', 'GGG'))),
  far_syn_and_missense = list(slice = 3L, edits = list(v(9, 'A', 'G'), v(5, 'C', 'G'))),
  far_syn_and_syn = list(slice = 3L, edits = list(v(9, 'A', 'G'), v(27, 'A', 'G'))),
  far_missense_and_missense = list(slice = 3L, edits = list(v(5, 'C', 'G'), v(38, 'G', 'C'))),
  two_far_insertions = list(slice = 3L, edits = list(v(12, 'G', 'GTCT'), v(39, 'G', 'GACT'))),
  insertion_and_deletion_far = list(slice = 3L, edits = list(v(12, 'G', 'GTCT'), v(54, 'GGAG', 'G'))),
  two_far_deletions = list(slice = 3L, edits = list(v(12, 'GAAA', 'G'), v(54, 'GGAG', 'G'))),
  stop_created_by_pair = list(slice = 4L, edits = list(v(29, 'T', 'A'), v(30, 'C', 'A'))),
  start_lost_by_pair = list(slice = 5L, edits = list(v(1, 'A', 'C'), v(2, 'T', 'C'))),
  stop_lost_by_pair = list(slice = 5L, edits = list(v(70, 'T', 'C'), v(71, 'A', 'G'))),
  stop_retained_single = list(slice = 5L, edits = list(v(72, 'A', 'G'))))

# ---- independent policy for the frame-preserving scenarios above ------------------------------
# Identical peptide -> synonymous. Otherwise, by the normalized edit path: substitutions only ->
# missense; every edit a pure insertion -> inframe_insertion; every edit a pure deletion ->
# inframe_deletion; anything else that preserves the frame -> protein_altering_variant.
edit_lane <- function(seq, edits, lane) {
    for (e in edits[order(-vapply(edits, function(x) x$pos, numeric(1L)))]) {
        if (!(lane %in% e$lane)) next
        stopifnot(sub_at(seq, e$pos, nchar(e$ref)) == e$ref)
        seq <- paste0(substr(seq, 1L, e$pos - 1L), e$alt, substring(seq, e$pos + nchar(e$ref)))
    }
    seq
}
impact_of <- function(terms) {
    if (!length(terms)) return('')
    levels <- c(stop_gained = 'HIGH', start_lost = 'HIGH', stop_lost = 'HIGH', frameshift_variant = 'HIGH',
                missense_variant = 'MODERATE', inframe_insertion = 'MODERATE', inframe_deletion = 'MODERATE',
                protein_altering_variant = 'MODERATE', synonymous_variant = 'LOW', stop_retained_variant = 'LOW')
    order <- c('HIGH', 'MODERATE', 'LOW')
    order[min(match(levels[terms], order))]
}
ref_protein <- translate(cds)
classify <- function(edits, lane) {
    mine <- Filter(function(e) lane %in% e$lane, edits)
    if (!length(mine)) return(list(terms = character(), protein = ''))
    alt <- edit_lane(cds, edits, lane)
    stopifnot(all(vapply(mine, function(e) (nchar(e$alt) - nchar(e$ref)) %% 3L, numeric(1L)) == 0),
              nchar(alt) %% 3L == 0L)
    protein <- translate(alt)
    first <- match('*', protein)
    visible <- if (is.na(first)) protein else protein[seq_len(first)]
    ref_n <- length(ref_protein)
    terms <- if (substr(alt, 1L, 3L) != 'ATG') 'start_lost' else if (is.na(first)) 'stop_lost' else
        if (first < length(protein)) 'stop_gained' else
        if (sub_at(alt, nchar(alt) - 2L, 3L) != sub_at(cds, 70L, 3L)) 'stop_retained_variant' else character()
    if (!length(terms)) {
        a <- protein[-length(protein)]; r <- ref_protein[-ref_n]
        subs <- all(vapply(mine, function(e) nchar(e$ref) == nchar(e$alt), logical(1L)))
        ins <- all(vapply(mine, function(e) nchar(e$ref) < nchar(e$alt) && startsWith(e$alt, e$ref),
                          logical(1L)))
        del <- all(vapply(mine, function(e) nchar(e$alt) < nchar(e$ref) && startsWith(e$ref, e$alt),
                          logical(1L)))
        terms <- if (identical(a, r)) 'synonymous_variant' else if (subs) 'missense_variant' else
            if (ins) 'inframe_insertion' else if (del) 'inframe_deletion' else 'protein_altering_variant'
    }
    list(terms = terms, protein = paste(visible, collapse = ''))
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
        lo <- min(G(1L), G(72L)); hi <- max(G(1L), G(72L))
        gene <- paste0('gene:G_', name); tx <- paste0('transcript:T_', name)
        gff <- c(gff,
          paste(name, 'fixture', 'gene', lo, hi, '.', layout$strand, '.', paste0('ID=', gene, ';biotype=protein_coding'), sep = '\t'),
          paste(name, 'fixture', 'mRNA', lo, hi, '.', layout$strand, '.', paste0('ID=', tx, ';Parent=', gene, ';biotype=protein_coding'), sep = '\t'))
        for (k in seq_along(ends)) {
            e1 <- G(starts[k]); e2 <- G(ends[k]); a <- min(e1, e2); b <- max(e1, e2)
            phase <- (3L - (starts[k] - 1L) %% 3L) %% 3L
            for (feature in c('exon', 'CDS'))
                gff <- c(gff, paste(name, 'fixture', feature, a, b, '.', layout$strand,
                                    if (feature == 'CDS') phase else '.', paste0('Parent=', tx), sep = '\t'))
        }
        # VCF records (genome, forward strand, left-anchored for indels).
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
            gold[[length(gold) + 1L]] <- data.frame(case = name, layout = layout_name, scenario = scenario_name,
                lane = lane, so_terms = paste(res$terms, collapse = ','), impact = impact_of(res$terms),
                protein = res$protein, classifier_slice = scenario$slice)
        }
    }
}
records <- do.call(rbind, records)
golden <- do.call(rbind, gold)
stopifnot(!anyDuplicated(contigs), nrow(golden) == 2L * length(contigs))

# Every scenario keeps the reference in the other lane empty (the trans and hom_alt scenarios excepted).
stopifnot(all(golden$so_terms[golden$lane == 1L & golden$scenario == 'cis_mis_pair'] == ''))

verified_write <- function(lines, path) {
    lines <- unname(lines)
    if (file.exists(path) && !identical(readLines(path, warn = FALSE), lines))
        stop('fixture differs from independent generator: ', path)
    writeLines(lines, path)
}
tsv <- function(frame) capture.output(write.table(frame, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
verified_write(fasta, file.path(out, 'same_codon.fa'))
offset <- 0L; fai <- character()
for (i in seq_along(contigs)) {
    offset <- offset + nchar(contigs[i]) + 2L
    fai <- c(fai, paste(contigs[i], lengths[i], offset, lengths[i], lengths[i] + 1L, sep = '\t'))
    offset <- offset + lengths[i] + 1L
}
verified_write(fai, file.path(out, 'same_codon.fa.fai'))
verified_write(gff, file.path(out, 'same_codon.gff3'))
vcf <- records[order(match(records$chrom, contigs), records$pos, records$order), ]
verified_write(c('##fileformat=VCFv4.2',
  vapply(seq_along(contigs), function(i) paste0('##contig=<ID=', contigs[i], ',length=', lengths[i], '>'), character(1L)),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS',
  paste(vcf$chrom, vcf$pos, '.', vcf$ref, vcf$alt, '.', 'PASS', '.', 'GT', vcf$gt, sep = '\t')),
  file.path(out, 'same_codon.vcf'))
verified_write(tsv(golden), file.path(out, 'same_codon_goldens.tsv'))
cat(length(contigs), 'contigs;', nrow(vcf), 'VCF records;', nrow(golden), 'lane goldens\n')
