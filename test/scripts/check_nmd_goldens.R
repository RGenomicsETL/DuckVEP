#!/usr/bin/env Rscript
# NMD attribution (coding-v1 slice 6, rule ejc50-v1) contract check, offline, base R, no DuckVEP code and no csq.
#  1. regenerate the fixtures and goldens (fails if a committed file differs) and the model tables;
#  2. re-derive every lane golden a second, structurally different way from the committed FASTA, GFF3 and VCF alone:
#     apply the lane's VCF records to the genome base by base, keep each base's exon, CDS membership and reference-CDS origin,
#     splice and orient the edited transcript, translate the edited CDS, and apply the rule to the edited spliced coordinates
#     (the generator instead edits the mRNA in transcript orientation). Prediction, S, J, the distance and the attributed
#     and post-stop edit positions must equal nmd_goldens.tsv exactly;
#  3. check that the fixture covers what the contract's slice 6 gate names (J - S = 49/50/51 and the boundary 0 and 52; stop in the
#     last and in the penultimate exon; intronless; an upstream indel shifting S and J together; a post-stop indel in the
#     penultimate exon shifting J only and one in the last exon shifting neither; a frameshift with a premature stop and one that
#     runs off; both strands; 1 to 4 exons; with and without a UTR).
stopifnot(system2('Rscript', c('--vanilla', 'test/scripts/generate_nmd_goldens.R')) == 0L,
          system2('Rscript', c('--vanilla', 'test/scripts/generate_haplotype_models.R', 'nmd')) == 0L)
root <- 'test/data/haplotype/'
codons <- c('TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
            'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
            'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
            'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG')
code <- setNames(strsplit('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '')[[1L]], codons)
translate <- function(s) if (nchar(s) < 3L) character() else
    unname(code[substring(s, seq.int(1L, nchar(s) - 2L, 3L), seq.int(3L, nchar(s), 3L))])
comp <- c(A = 'T', C = 'G', G = 'C', T = 'A')
revcomp <- function(s) paste(rev(comp[strsplit(s, '')[[1L]]]), collapse = '')

lines <- readLines(paste0(root, 'nmd.fa'))
genome <- setNames(lines[c(FALSE, TRUE)], sub('^>', '', lines[c(TRUE, FALSE)]))
gff <- read.delim(paste0(root, 'nmd.gff3'), header = FALSE, comment.char = '#', stringsAsFactors = FALSE,
    col.names = c('seq', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'))
vcf <- read.delim(paste0(root, 'nmd.vcf'), header = FALSE, comment.char = '#', stringsAsFactors = FALSE,
    col.names = c('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'gt'))
gold <- read.delim(paste0(root, 'nmd_goldens.tsv'), stringsAsFactors = FALSE, na.strings = '.', colClasses = c(
    stop_end = 'integer', junction = 'integer', distance = 'integer', applied_positions = 'character', post_stop_positions = 'character'))
gold$applied_positions[is.na(gold$applied_positions)] <- '.'; gold$post_stop_positions[is.na(gold$post_stop_positions)] <- '.'

derive <- function(name, lane) {
    exon <- gff[gff$seq == name & gff$type == 'exon', ]; cds <- gff[gff$seq == name & gff$type == 'CDS', ]
    minus <- unique(exon$strand) == '-'
    g <- strsplit(genome[[name]], '')[[1L]]; n <- length(g)
    exon_no <- integer(n); in_cds <- logical(n)
    ex <- exon[order(exon$start, decreasing = minus), ]            # transcript order
    for (k in seq_len(nrow(ex))) exon_no[ex$start[k]:ex$end[k]] <- k
    for (k in seq_len(nrow(cds))) in_cds[cds$start[k]:cds$end[k]] <- TRUE
    # origin: 1-based index of the base in the reference CDS (transcript orientation), NA outside it
    origin <- rep(NA_integer_, n); walk <- which(in_cds); if (minus) walk <- rev(walk); origin[walk] <- seq_along(walk)
    ref_cds <- paste(if (minus) comp[g[walk]] else g[walk], collapse = '')
    rows <- vcf[vcf$chrom == name & (vcf$gt == '1|1' | (lane == 0L & vcf$gt == '1|0') | (lane == 1L & vcf$gt == '0|1')), ]
    if (!nrow(rows)) return(list(carrier = FALSE))
    # per-base state; an edit replaces genome bases [start, start + ref_len) by alt, in descending order so coordinates hold
    st <- list(base = g, exon = exon_no, cds = in_cds, origin = origin, tag = integer(n), after = integer(n))
    order_desc <- order(rows$pos, decreasing = TRUE)
    for (i in order_desc) {
        r <- rows[i, ]; ref <- r$ref; alt <- r$alt; k <- 0L
        while (nchar(ref) > 0L && nchar(alt) > 0L && substr(ref, 1L, 1L) == substr(alt, 1L, 1L) && nchar(ref) != nchar(alt)) {
            ref <- substring(ref, 2L); alt <- substring(alt, 2L); k <- k + 1L
        }
        s <- r$pos + k; rl <- nchar(ref)
        stopifnot(identical(paste(st$base[s:(s + max(rl, 1L) - 1L)][seq_len(rl)], collapse = ''), ref))
        anchor <- if (rl) s else s - 1L                            # base whose exon an inserted base takes
        lo <- max(1L, s - 1L); hi <- min(length(st$base), s + rl)
        if (nchar(ref) != nchar(alt) && length(unique(st$exon[lo:hi])) != 1L) stop(sprintf('%s: edit at genome %d (%s) is next to an exon boundary: %s', name, r$pos, r$ref, paste(st$exon[lo:hi], collapse = ',')))
        new <- if (nchar(alt)) strsplit(alt, '')[[1L]] else character()
        m <- length(new); before <- seq_len(s - 1L); after <- seq.int(s + rl, length(st$base))
        if (s + rl > length(st$base)) after <- integer()
        pick <- function(v, fill) c(v[before], rep(fill, m), v[after])
        e <- st$exon[anchor]; c_ <- st$cds[anchor]
        st <- list(base = c(st$base[before], new, st$base[after]), exon = pick(st$exon, e), cds = pick(st$cds, c_),
                   origin = pick(st$origin, NA_integer_), tag = pick(st$tag, i), after = c(st$after[before], rep(0L, m), st$after[after]))
        if (!m && rl) {                                            # a pure deletion tags both neighbours
            if (s > 1L) st$after[s - 1L] <- i; if (s <= length(st$base)) st$after[s] <- i
        }
    }
    keep <- which(st$exon > 0L); if (minus) keep <- rev(keep)
    tx_base <- if (minus) unname(comp[st$base[keep]]) else st$base[keep]
    tx_exon <- st$exon[keep]; tx_cds <- st$cds[keep]; tx_origin <- st$origin[keep]; tx_tag <- st$tag[keep]; tx_after <- st$after[keep]
    cds_at <- which(tx_cds); first_cds <- cds_at[1L]; stopifnot(identical(cds_at, seq.int(first_cds, length.out = length(cds_at))))
    u5 <- first_cds - 1L
    edited <- paste(tx_base[cds_at], collapse = '')
    protein <- translate(edited); first <- match('*', protein); ref_protein <- translate(ref_cds)
    n_exons <- nrow(ex)
    res <- list(carrier = TRUE, nmd = 'unknown', premature = FALSE, S = NA_integer_, J = NA_integer_, D = NA_integer_,
                applied = '.', post = '.')
    if (substr(edited, 1L, 3L) != 'ATG' || is.na(first)) return(res)
    term_bases <- which(!is.na(tx_origin[cds_at]) & tx_origin[cds_at] > nchar(ref_cds) - 3L)
    term_start <- if (length(term_bases)) min(term_bases) else nchar(edited) + 1L
    res$premature <- (3L * first - 2L) < term_start && !identical(protein[seq_len(first - 1L)], ref_protein[-length(ref_protein)])
    if (!res$premature) { res$nmd <- 'not_applicable'; return(res) }
    res$S <- u5 + 3L * first
    if (n_exons > 1L) res$J <- max(which(tx_exon == n_exons - 1L))
    res$D <- if (is.na(res$J)) NA_integer_ else res$J - res$S
    res$nmd <- if (n_exons == 1L) 'escape' else if (res$J - res$S > 50L) 'trigger' else 'escape'
    # An edit is post-stop when its altered block starts at or after the end of the stop codon in the edited CDS.
    start <- vapply(seq_len(nrow(rows)), function(i) {
        at <- which(tx_tag[cds_at] == i)
        if (length(at)) min(at) - 1L else max(which(tx_after[cds_at] == i))   # 0-based; a deletion: the base after it
    }, numeric(1L))
    res$applied <- if (any(start < 3L * first)) paste(sort(rows$pos[start < 3L * first]), collapse = ';') else '.'
    res$post <- if (any(start >= 3L * first)) paste(sort(rows$pos[start >= 3L * first]), collapse = ';') else '.'
    res
}

bad <- character(); counts <- c(carriers = 0L, reference = 0L)
for (i in seq_len(nrow(gold))) {
    r <- gold[i, ]; d <- derive(r$case, r$lane)
    if (!d$carrier) {
        counts[['reference']] <- counts[['reference']] + 1L
        if (r$has_carrier || r$nmd != 'not_applicable') bad <- c(bad, r$case)
        next
    }
    counts[['carriers']] <- counts[['carriers']] + 1L
    ok <- r$has_carrier && identical(d$nmd, r$nmd) && identical(d$premature, r$premature) &&
        identical(as.integer(d$S), as.integer(r$stop_end)) && identical(as.integer(d$J), as.integer(r$junction)) &&
        identical(as.integer(d$D), as.integer(r$distance)) &&
        (!d$premature || (identical(d$applied, r$applied_positions) && identical(d$post, r$post_stop_positions))) &&
        (d$premature || (r$applied_positions == '.' && r$post_stop_positions == '.'))
    if (!ok) bad <- c(bad, paste(r$case, r$lane))
}
if (length(bad)) stop('genome-derived NMD differs from the goldens: ', paste(head(bad, 10L), collapse = ', '))
stopifnot(nrow(gold) == 468L, counts[['carriers']] == 250L, counts[['reference']] == 218L)

# ---- coverage of the slice 6 gate --------------------------------------------------------------
car <- gold[gold$has_carrier, ]
dist <- function(s) sort(unique(car$distance[car$scenario %in% s]))
stopifnot(identical(dist(c('d49')), 49L), identical(dist('d50'), 50L), identical(dist('d51'), 51L), identical(dist('d52'), 52L),
          identical(dist('d0_stop_ends_exon'), 0L),
          all(car$nmd[car$scenario == 'd49'] == 'escape'), all(car$nmd[car$scenario == 'd50'] == 'escape'),
          all(car$nmd[car$scenario == 'd51'] == 'trigger'), all(car$nmd[car$scenario == 'd52'] == 'trigger'),
          all(car$distance[car$scenario %in% c('last_exon_stop', 'last_of_three')] < 0L),          # stop in the last exon
          all(car$nmd[car$scenario %in% c('pen_d50', 'pen_exon_escape')] == 'escape'), all(car$nmd[car$scenario == 'pen_d51'] == 'trigger'),
          all(car$exons[car$scenario == 'intronless'] == 1L), all(car$nmd[grepl('^intronless', car$scenario)] == 'escape'),
          all(is.na(car$junction[grepl('^intronless', car$scenario)])),
          all(car$distance[car$scenario %in% c('up_ins3_d51', 'up_del3_d51', 'up_ins3_d50', 'up_del3_d50')] %in% c(50L, 51L)),
          all(car$nmd[car$scenario %in% c('up_ins3_d51', 'up_del3_d51')] == 'trigger'),
          all(car$nmd[car$scenario %in% c('up_ins3_d50', 'up_del3_d50')] == 'escape'),
          all(car$distance[car$scenario == 'post_ins_pen_shifts_j'] == 52L), all(car$distance[car$scenario == 'post_del_pen_shifts_j'] == 49L),
          all(car$distance[car$scenario == 'post_indel_last_exon_d49'] == 49L), all(car$distance[car$scenario == 'post_indel_last_exon_d51'] == 51L),
          any(car$scenario == 'fs_d51' & car$nmd == 'trigger'), all(car$nmd[grepl('^run_off', car$scenario)] == 'unknown'),
          setequal(unique(car$strand), c('+', '-')), setequal(unique(car$exons), 1:4), setequal(unique(car$utr), c('u0', 'u1')),
          setequal(unique(car$nmd), c('unknown', 'not_applicable', 'escape', 'trigger')),
          all(car$nmd[car$premature] %in% c('escape', 'trigger')), !any(car$premature[car$nmd %in% c('unknown', 'not_applicable')]),
          all(car$applied_positions[car$premature] != '.'), any(car$post_stop_positions != '.'))
cat('nmd goldens: 250 carrier lanes and 218 reference lanes re-derived from the genome; the slice 6 gate is covered\n')
