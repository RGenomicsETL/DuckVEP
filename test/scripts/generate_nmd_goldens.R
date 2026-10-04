#!/usr/bin/env Rscript
# NMD attribution (rule ejc50) fixtures and goldens. Base R only: no DuckVEP code and no
# csq (csq reports NMD_transcript only as a biotype marker and predicts no NMD). Each scenario edits the spliced
# mRNA (5' UTR + CDS + 3' UTR), translates the edited CDS with the standard code, follows every base's exon through
# the edits, and applies the rule to the edited transcript:
#   * a newly premature stop (the first stop of the edited CDS ends before the reference terminator and the peptide
#     before it is not the reference peptide) -> trigger when J - S > 50 else escape, intronless -> escape;
#     S = final nucleotide of the first stop codon, J = final nucleotide of the penultimate exon, both 1-based in
#     edited mRNA coordinates (5' UTR included), so an indel before or inside the penultimate exon moves J and an
#     indel before the stop moves S, whether or not that edit is the one that creates the stop;
#   * a lost start, or an edited CDS with no stop (a frame that runs off, stop_lost) -> unknown;
#   * any other known termination (synonymous, missense, inframe, stop_retained) -> not_applicable;
#   * a lane with no ALT carries no carrier: not_applicable (no DuckVEP row exists in strict decoded output).
# Attribution: the edits applied up to and including the stop codon (an edit island starting inside or before it),
# by VCF position; post-stop edits are listed apart and never attributed.
# Layouts: both strands, 1 to 3 (and 4) exons, each with and without a UTR (u1: 13-nt 5' UTR, 17-nt 3' UTR). An exon
# boundary of 0 or of the CDS length (u1 only) is a UTR-only exon. Edits never sit next to an exon boundary (checked),
# so an insertion cannot be attributed to either exon.
out <- 'test/data/haplotype'
dir.create(out, recursive = TRUE, showWarnings = FALSE)

caa <- paste0('ATG', 'AGC', paste(rep('CAA', 37L), collapse = ''), 'TAA')
mix <- 'ATGAGTCCCATCCGAAACTTCCTCCGTGATTCGGAATTACAAAAACATGTCTCGGGCCGTTGCGCCCCAAAAGCTATCCTCTCCCTACATGGCCCGTTGTGTAATTAA'
cdsets <- list(caa = caa, mix = mix)
stopifnot(nchar(caa) == 120L, nchar(mix) == 108L)
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
revcomp <- function(s) paste(rev(c(A = 'T', C = 'G', G = 'C', T = 'A')[strsplit(s, '')[[1L]]]), collapse = '')
sub_at <- function(s, p, n) substr(s, p, p + n - 1L)
utr_variants <- list(u0 = list(five = '', three = ''),
                     u1 = list(five = 'GGCGCTTCCGGAC', three = 'CCTGAGGATCCAGTTTA'))
intron <- 'GTAAGTCTGATCGAGTTCGACTTGACGTTAAGCTTACCAG'
flank <- 'CTAGGATCCTAGCATG'

# ---- scenarios -----------------------------------------------------------------------------------
# Edits are (pos, ref, alt) in CDS orientation and coordinates, VCF-anchored for indels. `exons` are the cumulative exon
# ends in CDS coordinates (the last exon runs to the end of the mRNA); `expect` documents the intended lane-0 result and
# distance D = J - S (NA when there is no junction or nothing premature) and is checked against the policy.
E <- function(cds) list(
    b = function(p) sub_at(cds, p, 1L),
    v = function(pos, ref, alt, lane = 0L) list(pos = pos, ref = ref, alt = alt, lane = lane),
    ins = function(p, bases, lane = 0L) list(pos = p, ref = sub_at(cds, p, 1L), alt = paste0(sub_at(cds, p, 1L), bases), lane = lane),
    del = function(p, n, lane = 0L) list(pos = p, ref = sub_at(cds, p, n + 1L), alt = sub_at(cds, p, 1L), lane = lane),
    snv = function(p, alt, lane = 0L) { stopifnot(sub_at(cds, p, 1L) != alt); list(pos = p, ref = sub_at(cds, p, 1L), alt = alt, lane = lane) })
sc <- function(cds, exons, edits, expect, utr = c('u0', 'u1')) list(cds = cds, exons = exons, edits = edits, expect = expect, utr = utr)
x <- function(nmd, d = NA_integer_) list(nmd = nmd, d = d)
scenarios <- list(
  # -- threshold: stop codon 17 (CDS 49-51, S = 51) against a penultimate exon ending 49, 50 and 51 bases later
  d49 = sc('caa', c(100L, 120L), function(e) list(e$snv(49, 'T')), x('escape', 49L)),
  d50 = sc('caa', c(101L, 120L), function(e) list(e$snv(49, 'T')), x('escape', 50L)),
  d51 = sc('caa', c(102L, 120L), function(e) list(e$snv(49, 'T')), x('trigger', 51L)),
  d52 = sc('caa', c(103L, 120L), function(e) list(e$snv(49, 'T')), x('trigger', 52L)),
  d0_stop_ends_exon = sc('caa', c(51L, 120L), function(e) list(e$snv(49, 'T')), x('escape', 0L)),
  # -- the same distances with the stop in the penultimate exon of three, and in the first exon of three
  pen_d50 = sc('caa', c(40L, 101L, 120L), function(e) list(e$snv(49, 'T')), x('escape', 50L)),
  pen_d51 = sc('caa', c(40L, 102L, 120L), function(e) list(e$snv(49, 'T')), x('trigger', 51L)),
  first_of_three_trigger = sc('caa', c(30L, 100L, 120L), function(e) list(e$snv(10, 'T')), x('trigger', 88L)),
  first_of_three_escape = sc('caa', c(60L, 90L, 120L), function(e) list(e$snv(49, 'T')), x('escape', 39L)),
  pen_exon_escape = sc('caa', c(30L, 100L, 120L), function(e) list(e$snv(88, 'T')), x('escape', 10L)),
  # -- stop in the last exon (J is before S) and in the last of four exons
  last_exon_stop = sc('caa', c(60L, 120L), function(e) list(e$snv(100, 'T')), x('escape', -42L)),
  last_of_three = sc('caa', c(30L, 60L, 120L), function(e) list(e$snv(100, 'T')), x('escape', -42L)),
  four_exons_escape = sc('caa', c(20L, 40L, 60L, 120L), function(e) list(e$snv(10, 'T')), x('escape', 48L)),
  four_exons_first_far = sc('caa', c(20L, 90L, 105L, 120L), function(e) list(e$snv(10, 'T')), x('trigger', 93L)),
  # -- intronless
  intronless = sc('caa', 120L, function(e) list(e$snv(49, 'T')), x('escape')),
  intronless_frame = sc('mix', 108L, function(e) list(e$ins(10, 'G')), x('escape')),
  intronless_early_stop = sc('caa', 120L, function(e) list(e$snv(10, 'T')), x('escape')),
  # -- an upstream indel moves S and J together: the distance is the unedited one (a S-only or J-only shift gets it wrong)
  up_ins3_d51 = sc('caa', c(102L, 120L), function(e) list(e$ins(30, 'GGG'), e$snv(49, 'T')), x('trigger', 51L)),
  up_ins3_d50 = sc('caa', c(101L, 120L), function(e) list(e$ins(30, 'GGG'), e$snv(49, 'T')), x('escape', 50L)),
  up_del3_d50 = sc('caa', c(101L, 120L), function(e) list(e$del(30, 3), e$snv(49, 'T')), x('escape', 50L)),
  up_del3_d51 = sc('caa', c(102L, 120L), function(e) list(e$del(30, 3), e$snv(49, 'T')), x('trigger', 51L)),
  up_ins6_three_exons = sc('caa', c(20L, 102L, 120L), function(e) list(e$ins(30, 'GGGGGG'), e$snv(49, 'T')), x('trigger', 51L)),
  up_indel_in_pen_exon = sc('caa', c(20L, 102L, 120L), function(e) list(e$ins(30, 'GGG'), e$snv(49, 'T')), x('trigger', 51L)),
  up_del_in_first_exon_of_three = sc('caa', c(20L, 101L, 120L), function(e) list(e$del(12, 3), e$snv(49, 'T')), x('escape', 50L)),
  # -- an indel after the stop inside the penultimate exon moves J only; one in the last exon moves neither
  post_ins_pen_shifts_j = sc('caa', c(60L, 100L, 120L), function(e) list(e$snv(49, 'T'), e$ins(80, 'GGG')), x('trigger', 52L)),
  post_del_pen_shifts_j = sc('caa', c(60L, 103L, 120L), function(e) list(e$snv(49, 'T'), e$del(80, 3)), x('escape', 49L)),
  post_ins_pen_two_exon = sc('caa', c(100L, 120L), function(e) list(e$snv(49, 'T'), e$ins(80, 'GGG')), x('trigger', 52L)),
  post_indel_last_exon_d49 = sc('caa', c(60L, 100L, 120L), function(e) list(e$snv(49, 'T'), e$ins(110, 'GGG')), x('escape', 49L)),
  post_indel_last_exon_d51 = sc('caa', c(60L, 102L, 120L), function(e) list(e$snv(49, 'T'), e$del(110, 3)), x('trigger', 51L)),
  post_ins_moves_j_no_flip = sc('caa', c(60L, 102L, 120L), function(e) list(e$snv(49, 'T'), e$ins(80, 'GGG')), x('trigger', 54L)),
  post_snv_after_stop = sc('caa', c(102L, 120L), function(e) list(e$snv(49, 'T'), e$snv(70, 'G')), x('trigger', 51L)),
  post_snv_last_exon = sc('caa', c(100L, 120L), function(e) list(e$snv(49, 'T'), e$snv(110, 'G')), x('escape', 49L)),
  # -- frameshifts with a premature stop (mix: the +1 frame reads a stop ending at edited cDNA 30; J follows the insertion)
  fs_d49 = sc('mix', c(78L, 108L), function(e) list(e$ins(10, 'G')), x('escape', 49L)),
  fs_d50 = sc('mix', c(79L, 108L), function(e) list(e$ins(10, 'G')), x('escape', 50L)),
  fs_d51 = sc('mix', c(80L, 108L), function(e) list(e$ins(10, 'G')), x('trigger', 51L)),
  fs_three_exons = sc('mix', c(20L, 80L, 108L), function(e) list(e$ins(10, 'G')), x('trigger', 51L)),
  fs_deletion_minus2 = sc('mix', c(80L, 108L), function(e) list(e$del(10, 2)), x('trigger', 51L)),
  fs_stop_by_insertion = sc('caa', c(101L, 120L), function(e) list(e$ins(30, 'TAAG')), x('trigger', NA_integer_)),
  fs_post_stop_edits = sc('mix', c(80L, 108L), function(e) list(e$ins(10, 'G'), e$ins(90, 'C'), e$snv(100, 'A')), x('trigger', 51L)),
  # -- a frame that runs off the CDS has no termination: unknown (also with a junction)
  run_off_deletion = sc('mix', c(60L, 108L), function(e) list(e$del(50, 1)), x('unknown')),
  run_off_plus2 = sc('mix', c(60L, 108L), function(e) list(e$ins(40, 'AC')), x('unknown')),
  run_off_intronless = sc('mix', 108L, function(e) list(e$del(50, 1)), x('unknown')),
  # -- lost initiation and lost termination are unknown; known termination is not_applicable
  start_lost = sc('caa', c(100L, 120L), function(e) list(e$snv(1, 'C')), x('unknown')),
  start_lost_with_stop = sc('caa', c(100L, 120L), function(e) list(e$snv(1, 'C'), e$snv(49, 'T')), x('unknown')),
  stop_lost = sc('caa', c(100L, 120L), function(e) list(e$snv(118, 'C')), x('unknown')),
  missense = sc('caa', c(100L, 120L), function(e) list(e$snv(49, 'G')), x('not_applicable')),
  synonymous = sc('caa', c(100L, 120L), function(e) list(e$snv(51, 'G')), x('not_applicable')),
  stop_retained = sc('caa', c(100L, 120L), function(e) list(e$snv(120, 'G')), x('not_applicable')),
  inframe_deletion = sc('caa', c(100L, 120L), function(e) list(e$del(49, 3)), x('not_applicable')),
  inframe_insertion = sc('caa', c(100L, 120L), function(e) list(e$ins(49, 'GGG')), x('not_applicable')),
  frame_restored = sc('mix', c(70L, 108L), function(e) list(e$ins(30, 'T'), e$del(60, 1)), x('not_applicable')),
  # -- the haplotype decides, never a single allele: a stop made only by a cis pair, and its trans controls
  cis_pair_stop = sc('caa', c(57L, 120L), function(e) list(e$snv(4, 'T'), e$snv(6, 'A')), x('trigger', 51L)),
  cis_pair_stop_d50 = sc('caa', c(56L, 120L), function(e) list(e$snv(4, 'T'), e$snv(6, 'A')), x('escape', 50L)),
  trans_pair_halves = sc('caa', c(57L, 120L), function(e) list(e$snv(4, 'T', 0L), e$snv(6, 'A', 1L)), x('not_applicable')),
  trans_stop_and_missense = sc('caa', c(102L, 120L), function(e) list(e$snv(49, 'T', 0L), e$snv(70, 'G', 1L)), x('trigger', 51L)),
  trans_missense_only = sc('caa', c(102L, 120L), function(e) list(e$snv(50, 'G', 0L), e$snv(20, 'T', 1L)), x('not_applicable')),
  hom_alt_stop = sc('caa', c(102L, 120L), function(e) list(e$snv(49, 'T', 0:1)), x('trigger', 51L)),
  frame_pair_stop_second = sc('caa', c(101L, 120L), function(e) list(e$ins(12, 'T'), e$del(20, 1), e$snv(49, 'T')), x('escape', 50L)),
  # -- UTR-only exons (u1 only): the junction between the last coding exon and a 3' UTR exon counts
  utr_last_exon = sc('caa', c(120L, 120L), function(e) list(e$snv(49, 'T')), x('trigger', 69L), utr = 'u1'),
  utr_first_exon = sc('caa', c(0L, 120L), function(e) list(e$snv(49, 'T')), x('escape', -51L), utr = 'u1'),
  utr_first_and_last = sc('caa', c(0L, 100L, 120L, 120L), function(e) list(e$snv(49, 'T')), x('trigger', 69L), utr = 'u1'))

# ---- independent policy --------------------------------------------------------------------------
# An island is the differing part of a VCF record: the anchor (common prefix of an indel) is stripped.
island <- function(e) {
    r <- e$ref; a <- e$alt; k <- 0L
    while (nchar(r) > 0L && nchar(a) > 0L && substr(r, 1L, 1L) == substr(a, 1L, 1L) && nchar(r) != nchar(a)) {
        r <- substring(r, 2L); a <- substring(a, 2L); k <- k + 1L
    }
    list(start = e$pos + k, ref_len = nchar(r), alt = a)
}
# Edit the mRNA. Every edited base keeps the exon of the base it came from; an inserted or replacement base takes the
# exon of the bases it stands beside, which the boundary check makes unambiguous.
edit_mrna <- function(mrna, label, edits, lane, u5) {
    mine <- Filter(function(e) lane %in% e$lane, edits)
    isl <- lapply(mine, island)
    for (i in seq_along(isl)) { isl[[i]]$start <- isl[[i]]$start + u5; isl[[i]]$id <- i }
    isl <- isl[order(vapply(isl, function(z) z$start, numeric(1L)))]
    seq_out <- character(); lab_out <- integer(); at <- 1L; shift <- 0L; starts_out <- integer(); ids <- integer()
    for (z in isl) {
        lo <- max(1L, z$start - 1L); hi <- min(nchar(mrna), z$start + z$ref_len)
        if (length(unique(label[lo:hi])) != 1L) stop(sprintf('edit at mRNA %d next to an exon boundary (%s)', z$start, paste(label[lo:hi], collapse = ',')))   # no edit next to an exon boundary
        if (z$start > at) {
            seq_out <- c(seq_out, strsplit(substr(mrna, at, z$start - 1L), '')[[1L]])
            lab_out <- c(lab_out, label[at:(z$start - 1L)])
        }
        starts_out <- c(starts_out, length(seq_out)); ids <- c(ids, z$id)
        if (nchar(z$alt)) {
            seq_out <- c(seq_out, strsplit(z$alt, '')[[1L]])
            lab_out <- c(lab_out, rep(label[min(z$start, nchar(mrna))], nchar(z$alt)))
        }
        at <- z$start + z$ref_len
    }
    if (at <= nchar(mrna)) {
        seq_out <- c(seq_out, strsplit(substr(mrna, at, nchar(mrna)), '')[[1L]])
        lab_out <- c(lab_out, label[at:nchar(mrna)])
    }
    list(seq = paste(seq_out, collapse = ''), label = lab_out, island_start0 = starts_out, island_id = ids, n = length(isl),
         change = vapply(isl, function(z) nchar(z$alt) - z$ref_len, integer(1L)),
         ref_start = vapply(isl, function(z) z$start, numeric(1L)),
         ref_end = vapply(isl, function(z) z$start + z$ref_len - 1L, numeric(1L)))
}
policy <- function(cds, u5, u3, mrna, label, edits, lane, n_exons, ends) {
    if (!any(vapply(edits, function(e) lane %in% e$lane, logical(1L))))
        return(list(carrier = FALSE, nmd = 'not_applicable', premature = FALSE, S = NA_integer_, J = NA_integer_,
                    D = NA_integer_, applied = integer(), post = integer(), attributed = integer()))
    w <- edit_mrna(mrna, label, edits, lane, u5)
    edited_cds <- substr(w$seq, u5 + 1L, nchar(w$seq) - u3)
    protein <- translate(edited_cds); first <- match('*', protein)
    ref_protein <- translate(cds); nref <- length(ref_protein)
    res <- list(carrier = TRUE, nmd = 'unknown', premature = FALSE, S = NA_integer_, J = NA_integer_, D = NA_integer_,
                applied = integer(), post = integer(), attributed = integer())
    if (substr(edited_cds, 1L, 3L) != 'ATG' || is.na(first)) return(res)
    # The reference terminator (last three CDS bases) sits after every edit that ends before it.
    before <- sum(w$change[w$ref_end <= u5 + nchar(cds) - 3L])
    term_start <- nchar(cds) - 2L + before
    res$premature <- (3L * first - 2L) < term_start && !identical(protein[seq_len(first - 1L)], ref_protein[-nref])
    if (!res$premature) { res$nmd <- 'not_applicable'; return(res) }
    S <- u5 + 3L * first
    J <- if (n_exons > 1L) { m <- which(w$label == n_exons - 1L); stopifnot(length(m) > 0L, length(which(w$label == n_exons)) > 0L); max(m) } else NA_integer_
    res$S <- S; res$J <- J
    res$D <- if (is.na(J)) NA_integer_ else as.integer(J - S)
    res$nmd <- if (n_exons == 1L) 'escape' else if (J - S > 50L) 'trigger' else 'escape'
    # island start in edited CDS coordinates, 0-based: post_stop when it is at or after the end of the stop codon
    post <- w$island_start0 - u5 >= 3L * first
    res$applied <- w$island_id[!post]; res$post <- w$island_id[post]
    # NMD attribution: the applied edits plus the post-stop edits that changed length at or before the penultimate
    # exon's last base (reference mRNA coordinates), i.e. exactly the edits that moved J.
    moved <- if (n_exons > 1L) w$change != 0L & w$ref_start <= ends[n_exons - 1L] else logical(length(post))
    res$attributed <- w$island_id[!post | (post & moved)]
    res
}

# ---- genome layout and VCF mapping ---------------------------------------------------------------
records <- list(); gff <- '##gff-version 3'; fasta <- character(); gold <- list()
contigs <- character(); lengths <- integer()
for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    cds <- cdsets[[scenario$cds]]; L <- nchar(cds)
    edits <- scenario$edits(E(cds))
    for (strand in c('+', '-')) for (utr_name in scenario$utr) {
        minus <- strand == '-'; five <- utr_variants[[utr_name]]$five; three <- utr_variants[[utr_name]]$three
        u5 <- nchar(five); u3 <- nchar(three)
        mrna <- paste0(five, cds, three)
        cuts <- scenario$exons                       # cumulative exon ends in CDS coordinates
        cuts <- cuts[-length(cuts)]
        if (utr_name == 'u0') stopifnot(!any(cuts <= 0L | cuts >= L))
        ends <- c(u5 + cuts, nchar(mrna)); starts <- c(1L, head(ends, -1L) + 1L)
        stopifnot(all(ends > starts - 1L), !is.unsorted(ends, strictly = TRUE))
        n_exons <- length(ends)
        label <- findInterval(seq_len(nchar(mrna)) - 1L, ends) + 1L
        pieces <- substring(mrna, starts, ends)
        transcript_genome <- paste0(flank, paste(pieces, collapse = intron), flank)
        idx <- function(t) { k <- findInterval(t - 1L, ends) + 1L; nchar(flank) + t + (k - 1L) * nchar(intron) }
        n <- nchar(transcript_genome)
        genome <- if (minus) revcomp(transcript_genome) else transcript_genome
        G <- function(t) if (minus) n + 1L - idx(t) else idx(t)
        name <- paste(scenario_name, if (minus) 'minus' else 'plus', utr_name, sep = '__')
        contigs <- c(contigs, name); lengths <- c(lengths, n); fasta <- c(fasta, paste0('>', name), genome)
        lo <- min(G(1L), G(nchar(mrna))); hi <- max(G(1L), G(nchar(mrna)))
        gene <- paste0('gene:G_', name); tx <- paste0('transcript:T_', name)
        gff <- c(gff,
          paste(name, 'fixture', 'gene', lo, hi, '.', strand, '.', paste0('ID=', gene, ';biotype=protein_coding'), sep = '\t'),
          paste(name, 'fixture', 'mRNA', lo, hi, '.', strand, '.', paste0('ID=', tx, ';Parent=', gene, ';biotype=protein_coding'), sep = '\t'))
        for (k in seq_along(ends)) {
            e1 <- G(starts[k]); e2 <- G(ends[k])
            gff <- c(gff, paste(name, 'fixture', 'exon', min(e1, e2), max(e1, e2), '.', strand, '.', paste0('Parent=', tx), sep = '\t'))
            c1 <- max(starts[k], u5 + 1L); c2 <- min(ends[k], u5 + L)
            if (c1 <= c2) {
                g1 <- G(c1); g2 <- G(c2)
                phase <- (3L - (c1 - u5 - 1L) %% 3L) %% 3L
                gff <- c(gff, paste(name, 'fixture', 'CDS', min(g1, g2), max(g1, g2), '.', strand, phase, paste0('Parent=', tx), sep = '\t'))
            }
        }
        edit_pos <- integer(length(edits))
        for (i in seq_along(edits)) {
            e <- edits[[i]]; p <- e$pos + u5
            r <- nchar(e$ref); al <- nchar(e$alt)
            anchored <- r != al && (startsWith(e$alt, e$ref) || startsWith(e$ref, e$alt))
            if (!minus) {
                pos <- G(p); ref <- e$ref; alt <- e$alt
            } else if (!anchored) {
                pos <- G(p + r - 1L); ref <- revcomp(e$ref); alt <- revcomp(e$alt)
            } else if (r > al) { # deletion of transcript bases pos+al .. pos+r-1
                del_first <- p + al; del_last <- p + r - 1L
                low <- G(del_last); high <- G(del_first)
                pos <- low - 1L; ref <- substr(genome, low - 1L, high); alt <- substr(genome, low - 1L, low - 1L)
            } else { # insertion after transcript base pos
                pos <- G(p + 1L)
                ref <- substr(genome, pos, pos)
                alt <- paste0(ref, revcomp(substring(e$alt, r + 1L)))
            }
            gt <- if (length(e$lane) == 2L) '1|1' else if (e$lane == 0L) '1|0' else '0|1'
            edit_pos[i] <- pos
            records[[length(records) + 1L]] <- data.frame(chrom = name, pos = pos, ref = ref, alt = alt, gt = gt,
                order = length(records))
        }
        for (lane in 0:1) {
            res <- policy(cds, u5, u3, mrna, label, edits, lane, n_exons, ends)
            if (lane == 0L && !is.null(scenario$expect)) {
                if (!(identical(res$nmd, scenario$expect$nmd) &&
                      (is.na(scenario$expect$d) || identical(res$D, scenario$expect$d) || (n_exons == 1L && is.na(res$D)))))
                    stop(sprintf('%s: policy %s D=%s S=%s J=%s, expected %s D=%s', name, res$nmd, res$D, res$S, res$J,
                                 scenario$expect$nmd, scenario$expect$d))
            }
            gold[[length(gold) + 1L]] <- data.frame(case = name, scenario = scenario_name,
                strand = strand, utr = utr_name, exons = n_exons, lane = lane,
                has_carrier = res$carrier, nmd = res$nmd, premature = res$premature,
                stop_end = res$S, junction = res$J, distance = res$D,
                applied_positions = if (length(res$applied)) paste(sort(edit_pos[res$applied]), collapse = ';') else '.',
                nmd_positions = if (length(res$attributed)) paste(sort(edit_pos[res$attributed]), collapse = ';') else '.',
                post_stop_positions = if (length(res$post)) paste(sort(edit_pos[res$post]), collapse = ';') else '.')
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
verified_write(fasta, file.path(out, 'nmd.fa'))
offset <- 0L; fai <- character()
for (i in seq_along(contigs)) {
    offset <- offset + nchar(contigs[i]) + 2L
    fai <- c(fai, paste(contigs[i], lengths[i], offset, lengths[i], lengths[i] + 1L, sep = '\t'))
    offset <- offset + lengths[i] + 1L
}
verified_write(fai, file.path(out, 'nmd.fa.fai'))
verified_write(gff, file.path(out, 'nmd.gff3'))
vcf <- records[order(match(records$chrom, contigs), records$pos, records$order), ]
verified_write(c('##fileformat=VCFv4.2',
  vapply(seq_along(contigs), function(i) paste0('##contig=<ID=', contigs[i], ',length=', lengths[i], '>'), character(1L)),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS',
  paste(vcf$chrom, vcf$pos, '.', vcf$ref, vcf$alt, '.', 'PASS', '.', 'GT', vcf$gt, sep = '\t')),
  file.path(out, 'nmd.vcf'))
verified_write(tsv(golden), file.path(out, 'nmd_goldens.tsv'))
cat(length(contigs), 'contigs;', nrow(vcf), 'VCF records;', nrow(golden), 'lane goldens\n')
