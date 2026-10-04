#!/usr/bin/env Rscript
# Same-codon contract check, offline by default.
#  1. regenerate the base-R fixtures/goldens (fails if a committed fixture differs);
#  2. verify the pinned csq output against its receipt hashes, and re-run csq only when the
#     pinned bcftools is installed;
#  3. map the csq lane consequences through csq_so_map_v1.tsv and require exact agreement with the
#     goldens, apart from the divergences declared below and in test/data/haplotype/README.md;
#  4. check every csq amino-acid change against the base-R peptides.
# HAPLOTYPE_CSQ_REFRESH=1 rewrites same_codon_csq.tsv and its receipt (needs the pinned bcftools).
stopifnot(system2('Rscript', c('--vanilla', 'test/scripts/generate_same_codon_goldens.R')) == 0L)
root <- 'test/data/haplotype/'
sha <- function(path) substr(system2('sha256sum', path, stdout = TRUE), 1L, 64L)
csq <- Sys.getenv('HAPLOTYPE_CSQ', '/usr/local/bin/bcftools')
pinned_version <- '1.23.1-70-g6dbd8fef+htslib-1.22.1'
run_csq <- function() {
    vcf <- tempfile(fileext = '.vcf'); on.exit(unlink(vcf))
    stopifnot(system2(csq, c('csq', '-p', 's', '-n', '1024', '-f', paste0(root, 'same_codon.fa'),
        '-g', paste0(root, 'same_codon.gff3'), '-Ov', '-o', vcf, paste0(root, 'same_codon.vcf')),
        stderr = tempfile()) == 0L)
    system2(csq, c('query', '-f', shQuote('%CHROM\t%POS\t%INFO/BCSQ\t[%TBCSQ]\n'), vcf), stdout = TRUE)
}
have_csq <- file.exists(csq) && identical(system2(csq, '--version-only', stdout = TRUE), pinned_version)
if (nzchar(Sys.getenv('HAPLOTYPE_CSQ_REFRESH'))) {
    stopifnot(have_csq)
    writeLines(run_csq(), paste0(root, 'same_codon_csq.tsv'))
    key <- c(bcftools_commit = '6dbd8fef51e529755a4a81544075dc6a43ff3cd3',
             htslib_commit = 'c1f35d67dd5ff1e226d94abe7b850f28c60f5910',
             bcftools_version = pinned_version, htslib_version = '1.22.1', phase_mode = 's', ncsq = '1024',
             input_sha256 = sha(paste0(root, 'same_codon.vcf')),
             reference_sha256 = sha(paste0(root, 'same_codon.fa')),
             gff3_sha256 = sha(paste0(root, 'same_codon.gff3')),
             output_sha256 = sha(paste0(root, 'same_codon_csq.tsv')))
    writeLines(c('field\tvalue', paste(names(key), key, sep = '\t')), paste0(root, 'same_codon_csq.oracle.tsv'))
}
receipt <- read.delim(paste0(root, 'same_codon_csq.oracle.tsv'), check.names = FALSE)
pin <- setNames(receipt$value, receipt$field)
files <- c(input = 'same_codon.vcf', reference = 'same_codon.fa', gff3 = 'same_codon.gff3',
           output = 'same_codon_csq.tsv')
for (k in names(files)) stopifnot(identical(sha(paste0(root, files[[k]])), pin[[paste0(k, '_sha256')]]))
stopifnot(pin[['bcftools_version']] == pinned_version, pin[['phase_mode']] == 's', pin[['ncsq']] == '1024')
pinned <- readLines(paste0(root, 'same_codon_csq.tsv'))
if (have_csq) {
    stopifnot(identical(run_csq(), pinned))
} else message('pinned bcftools not found; checking the goldens against the committed csq output')

map <- read.delim(paste0(root, 'csq_so_map_v1.tsv'), check.names = FALSE, na.strings = character())
lookup <- setNames(map$so_term[map$kind == 'csq'], map$csq_term[map$kind == 'csq'])
gold <- read.delim(paste0(root, 'same_codon_goldens.tsv'), na.strings = character(), check.names = FALSE,
                   colClasses = c(lane = 'integer'))
transcripts <- read.delim(paste0(root, 'same_codon_transcripts.tsv'), stringsAsFactors = FALSE)
cds <- unique(transcripts$cds_sequence); stopifnot(length(cds) == 1L)
codons <- c('TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',
            'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',
            'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',
            'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG')
code <- setNames(strsplit('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '')[[1L]], codons)
ref_protein <- paste(code[substring(cds, seq.int(1L, nchar(cds) - 2L, 3L), seq.int(3L, nchar(cds), 3L))],
                     collapse = '')

# Per-case, per-lane csq consequence terms and amino-acid strings (TBCSQ column 4 = lane 0, 5 = lane 1).
rows <- strsplit(pinned, '\t', fixed = TRUE)
stopifnot(all(lengths(rows) == 5L))
lane_terms <- list(); lane_aa <- list()
for (row in rows) {
    for (lane in 0:1) {
        key <- paste(row[1L], lane, sep = '/')
        for (record in strsplit(row[4L + lane], ',', fixed = TRUE)[[1L]]) {
            if (record %in% c('.', '') || startsWith(record, '@') || startsWith(record, '*')) next
            fields <- strsplit(record, '|', fixed = TRUE)[[1L]]
            lane_terms[[key]] <- union(lane_terms[[key]], strsplit(fields[1L], '&', fixed = TRUE)[[1L]])
            if (length(fields) >= 6L) lane_aa[[key]] <- c(lane_aa[[key]], fields[6L])
        }
    }
}
unmapped <- setdiff(unlist(lane_terms), names(lookup))
stopifnot(!length(unmapped))
# The reduced set omits the lower-severity splice_region subeffect csq adds at exon edges, because
# splicing is fixed to the selected model (contract section 2). No other term is dropped.
mapped <- function(terms) sort(setdiff(unname(lookup[terms]), 'splice_region_variant'))
keys <- paste(gold$case, gold$lane, sep = '/')
csq_so <- vapply(keys, function(k) paste(mapped(lane_terms[[k]]), collapse = ','), character(1L))
policy_so <- vapply(gold$so_terms, function(x) paste(sort(strsplit(x, ',', fixed = TRUE)[[1L]]), collapse = ','),
                    character(1L))

# Declared divergences between the reduced whole-protein policy and csq, keyed by scenario. csq
# labels an inframe change by its length change alone and reports every codon's effect; the policy
# uses the normalized edit path and one whole-protein category.
declared <- c(
  delins_inframe = 'csq labels the 6>3 replacement inframe_deletion by length; the policy has a mixed replacement',
  far_syn_and_missense = 'csq reports the synonymous and the missense codon; the policy reports one category',
  insertion_and_deletion_far = 'csq reports both inframe terms; the policy has mixed frame-preserving edits')
csq_expected <- c(delins_inframe = 'inframe_deletion', far_syn_and_missense = 'missense_variant,synonymous_variant',
                  insertion_and_deletion_far = 'inframe_deletion,inframe_insertion')
differs <- csq_so != policy_so
divergent <- unique(gold$scenario[differs])
stopifnot(setequal(divergent, names(declared)))
for (name in names(declared)) {
    lanes <- gold$scenario == name & gold$lane == 0L
    stopifnot(all(differs[lanes]), all(csq_so[lanes] == csq_expected[[name]]),
              all(!differs[gold$scenario == name & gold$lane == 1L]))
}

# Every csq amino-acid change must describe the base-R peptides (reference at its position, alternate
# in the edited protein), so the SO comparison cannot pass on a mislabelled change. Declared csq
# quirk: for the second indel of two indels in one lane, csq prints an alternate residue that is not
# in the edited protein (e.g. 18QE>17R for the deletion of E, where the edited residue 17 is Q). The
# SO term is right; only its alternate-side string is unreliable, so that side is skipped for
# these two scenarios while their reference side and SO sets are still checked.
alt_quirk <- c('two_far_deletions', 'insertion_and_deletion_far')
checked <- 0L
for (i in seq_len(nrow(gold))) {
    for (aa in lane_aa[[keys[i]]]) {
        m <- regmatches(aa, regexec('^([0-9]+)([A-Z*]+)(?:>([0-9]+)([A-Z*]+))?$', aa))[[1L]]
        stopifnot(length(m) >= 3L)
        pos <- as.integer(m[2L])
        stopifnot(substr(ref_protein, pos, pos + nchar(m[3L]) - 1L) == m[3L])
        if (nzchar(m[4L]) && !(gold$scenario[i] %in% alt_quirk)) stopifnot(substr(gold$protein[i], as.integer(m[4L]),
            as.integer(m[4L]) + nchar(m[5L]) - 1L) == m[5L])
        checked <- checked + 1L
    }
}

expected <- data.frame(case = gold$case, lane = gold$lane, csq_so_terms = unname(csq_so),
    declared_divergence = ifelse(differs, declared[gold$scenario], '.'), stringsAsFactors = FALSE)
lines <- capture.output(write.table(expected, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
path <- paste0(root, 'same_codon_csq_expected.tsv')
if (nzchar(Sys.getenv('HAPLOTYPE_CSQ_REFRESH'))) writeLines(lines, path)
stopifnot(identical(readLines(path), lines))
lanes3 <- gold$classifier_slice == 3L
stopifnot(sum(lanes3) > 0L, all(csq_so[lanes3 & !differs] == policy_so[lanes3 & !differs]))
cat(sprintf(paste('%d lane goldens over %d layouts; %d csq lane strings mapped; %d agree exactly, %d in %d declared',
                  'divergences; %d csq amino-acid changes verified against base-R peptides\n'),
            nrow(gold), length(unique(gold$layout)), nrow(gold), sum(!differs), sum(differs), length(declared),
            checked))
