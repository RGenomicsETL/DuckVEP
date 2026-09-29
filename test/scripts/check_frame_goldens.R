#!/usr/bin/env Rscript
# Frame opening/restoration and stop-gain (coding-v1 slice 4) contract check, offline by default.
#  1. regenerate the base-R fixtures/goldens (fails if a committed fixture differs);
#  2. verify the pinned csq output against its receipt hashes, and re-run csq only when the pinned
#     bcftools is installed;
#  3. map the csq lane consequences through csq_so_map_v1.tsv and require exact agreement with the
#     goldens apart from the divergences declared below and in test/data/haplotype/README.md;
#  4. check every csq amino-acid change against the base-R peptides.
# HAPLOTYPE_CSQ_REFRESH=1 rewrites frame_csq.tsv, its receipt and frame_csq_expected.tsv (needs the pinned bcftools).
stopifnot(system2('Rscript', c('--vanilla', 'test/scripts/generate_frame_goldens.R')) == 0L)
root <- 'test/data/haplotype/'
sha <- function(path) substr(system2('sha256sum', path, stdout = TRUE), 1L, 64L)
csq <- Sys.getenv('HAPLOTYPE_CSQ', '/usr/local/bin/bcftools')
pinned_version <- '1.23.1-70-g6dbd8fef+htslib-1.22.1'
run_csq <- function() {
    vcf <- tempfile(fileext = '.vcf'); on.exit(unlink(vcf))
    stopifnot(system2(csq, c('csq', '-p', 's', '-n', '1024', '-f', paste0(root, 'frame.fa'),
        '-g', paste0(root, 'frame.gff3'), '-Ov', '-o', vcf, paste0(root, 'frame.vcf')),
        stderr = tempfile()) == 0L)
    system2(csq, c('query', '-f', shQuote('%CHROM\t%POS\t%INFO/BCSQ\t[%TBCSQ]\n'), vcf), stdout = TRUE)
}
have_csq <- file.exists(csq) && identical(system2(csq, '--version-only', stdout = TRUE), pinned_version)
if (nzchar(Sys.getenv('HAPLOTYPE_CSQ_REFRESH'))) {
    stopifnot(have_csq)
    writeLines(run_csq(), paste0(root, 'frame_csq.tsv'))
    key <- c(bcftools_commit = '6dbd8fef51e529755a4a81544075dc6a43ff3cd3',
             htslib_commit = 'c1f35d67dd5ff1e226d94abe7b850f28c60f5910',
             bcftools_version = pinned_version, htslib_version = '1.22.1', phase_mode = 's', ncsq = '1024',
             input_sha256 = sha(paste0(root, 'frame.vcf')),
             reference_sha256 = sha(paste0(root, 'frame.fa')),
             gff3_sha256 = sha(paste0(root, 'frame.gff3')),
             output_sha256 = sha(paste0(root, 'frame_csq.tsv')))
    writeLines(c('field\tvalue', paste(names(key), key, sep = '\t')), paste0(root, 'frame_csq.oracle.tsv'))
}
receipt <- read.delim(paste0(root, 'frame_csq.oracle.tsv'), check.names = FALSE)
pin <- setNames(receipt$value, receipt$field)
files <- c(input = 'frame.vcf', reference = 'frame.fa', gff3 = 'frame.gff3', output = 'frame_csq.tsv')
for (k in names(files)) stopifnot(identical(sha(paste0(root, files[[k]])), pin[[paste0(k, '_sha256')]]))
stopifnot(pin[['bcftools_version']] == pinned_version, pin[['phase_mode']] == 's', pin[['ncsq']] == '1024')
pinned <- readLines(paste0(root, 'frame_csq.tsv'))
if (have_csq) {
    stopifnot(identical(run_csq(), pinned))
} else message('pinned bcftools not found; checking the goldens against the committed csq output')

map <- read.delim(paste0(root, 'csq_so_map_v1.tsv'), check.names = FALSE, na.strings = character())
lookup <- setNames(map$so_term[map$kind == 'csq'], map$csq_term[map$kind == 'csq'])
gold <- read.delim(paste0(root, 'frame_goldens.tsv'), na.strings = character(), check.names = FALSE,
                   colClasses = c(lane = 'integer'))
transcripts <- read.delim(paste0(root, 'frame_transcripts.tsv'), stringsAsFactors = FALSE)
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
stopifnot(!length(setdiff(unlist(lane_terms), names(lookup))))
# The reduced set omits the lower-severity splice_region subeffect csq adds at exon edges, because
# splicing is fixed to the selected model (contract section 2). No other term is dropped.
mapped <- function(terms) sort(setdiff(unname(lookup[terms]), 'splice_region_variant'))
keys <- paste(gold$case, gold$lane, sep = '/')
csq_so <- vapply(keys, function(k) paste(mapped(lane_terms[[k]]), collapse = ','), character(1L))
policy_so <- vapply(gold$so_terms, function(x) paste(sort(strsplit(x, ',', fixed = TRUE)[[1L]]), collapse = ','),
                    character(1L))

# Declared divergences between the reduced whole-protein policy and csq, keyed by scenario and lane.
# 'both' means every ALT-carrying lane of the scenario. The value asserted is the mapped csq set, so a
# change in the pinned csq output cannot pass as the same divergence.
single_frameshift <- paste('csq reports a single frameshift record as frameshift only and does not add the',
                           'stop_gained of the early stop it creates; the policy reports both HIGH terms')
restored <- paste('csq adds the per-record lower-severity missense to inframe_altering when the frame is',
                  'restored; the policy reports the one whole-protein category protein_altering_variant')
run_off <- paste('csq reports a frame that runs off the CDS with no stop as frameshift only; the contract has',
                 'the reference termination abolished as well, so the policy reports frameshift_variant and stop_lost')
declared <- list(
  frame_open_early_stop = list(lane = 0L, csq = 'frameshift_variant', why = single_frameshift),
  frame_open_late_stop = list(lane = 0L, csq = 'frameshift_variant', why = single_frameshift),
  frame_open_minus2_late_stop = list(lane = 0L, csq = 'frameshift_variant', why = single_frameshift),
  early_stop_then_missense = list(lane = 0L, csq = 'frameshift_variant', why = single_frameshift),
  trans_control_of_pair = list(lane = 0L, csq = 'frameshift_variant', why = single_frameshift),
  hom_alt_frame_open = list(lane = 0:1, csq = 'frameshift_variant', why = single_frameshift),
  frame_open_runs_off = list(lane = 0L, csq = 'frameshift_variant', why = run_off),
  frame_open_plus2_runs_off = list(lane = 0L, csq = 'frameshift_variant', why = run_off),
  frame_open_before_terminal = list(lane = 0L, csq = 'frameshift_variant', why = run_off),
  `trans_control_of_pair:run_off` = list(lane = 1L, csq = 'frameshift_variant', why = run_off),
  restored_rescued = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_rescued_wide = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_near_introns = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_minus_plus = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_plus2_minus2 = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_minus2_plus2 = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_two_pairs = list(lane = 0L, csq = 'missense_variant,protein_altering_variant', why = restored),
  restored_identical_peptide_early = list(lane = 0L, csq = 'protein_altering_variant,synonymous_variant',
    why = paste('csq adds inframe_altering to an identical peptide after a restored frame; the contract maps',
                'an identical peptide to synonymous_variant')),
  restored_identical_peptide_late = list(lane = 0L, csq = 'protein_altering_variant,synonymous_variant',
    why = paste('csq adds inframe_altering to an identical peptide after a restored frame; the contract maps',
                'an identical peptide to synonymous_variant')),
  restored_three_insertions = list(lane = 0L, csq = 'inframe_insertion',
    why = paste('csq labels a net +3 of three frame-shifting insertions inframe_insertion by length;',
                'the policy has a restored frame after displacement, protein_altering_variant')),
  restored_ins1_ins2 = list(lane = 0L, csq = 'inframe_insertion',
    why = paste('csq labels a net +3 of a +1 and a +2 insertion inframe_insertion by length; the policy',
                'has a restored frame after displacement, protein_altering_variant')),
  restored_del1_del2 = list(lane = 0L, csq = 'inframe_deletion',
    why = paste('csq labels a net -3 of a -1 and a -2 deletion inframe_deletion by length; the policy',
                'has a restored frame after displacement, protein_altering_variant')),
  stop_gained_inframe_insertion = list(lane = 0L, csq = 'inframe_insertion,stop_gained',
    why = 'csq also reports the inframe_insertion carrying the stop; the reduced set reports the HIGH stop_gained only'),
  stop_gained_after_restoration = list(lane = 0L, csq = 'missense_variant,protein_altering_variant,stop_gained',
    why = paste('csq reports the rescued pair upstream of the new stop as well; the reduced set reports the',
                'HIGH stop_gained only, and no frameshift because the stop is outside the displaced frame')))
differs <- csq_so != policy_so
carrying <- policy_so != '' | csq_so != ''
divergent <- unique(gold$scenario[differs])
scenario_of <- function(name) sub(':.*$', '', name)
stopifnot(setequal(divergent, unique(scenario_of(names(declared)))))
declared_text <- rep('.', nrow(gold))
for (name in names(declared)) {
    d <- declared[[name]]
    scenario <- scenario_of(name)
    lanes <- gold$scenario == scenario & gold$lane %in% d$lane
    stopifnot(all(differs[lanes]), all(csq_so[lanes] == d$csq))
    declared_text[lanes] <- d$why
}
# Every divergent lane must be covered by a declaration (a scenario may carry several, one per lane group).
covered <- rep(FALSE, nrow(gold))
for (name in names(declared)) covered <- covered | (gold$scenario == scenario_of(name) & gold$lane %in% declared[[name]]$lane)
stopifnot(all(covered[differs]))

# Every csq amino-acid change must describe the base-R peptides (reference at its position, alternate in
# the edited protein), so the SO comparison cannot pass on a mislabelled change. Declared csq quirk: for
# a frameshifting deletion and for some compound records csq's alternate string differs from the edited
# protein at the residues next to the deletion (e.g. 20L for an edited residue 20R). The SO term is
# right; only the alternate side of the scenarios below is skipped, their reference side is still checked
# and each listed scenario must really show the mismatch, so the list cannot go stale.
alt_quirk <- c('restored_del1_del2', 'restored_minus_plus', 'restored_rescued_wide', 'trans_control_of_pair')
checked <- 0L; alt_checked <- 0L; quirk_seen <- character()
for (i in seq_len(nrow(gold))) {
    aas <- lane_aa[[keys[i]]]
    for (j in seq_along(aas)) {
        m <- regmatches(aas[j], regexec('^([0-9]+)([A-Z*]+)(?:>([0-9]+)([A-Z*]+))?$', aas[j]))[[1L]]
        stopifnot(length(m) >= 3L)
        pos <- as.integer(m[2L])
        stopifnot(substr(ref_protein, pos, pos + nchar(m[3L]) - 1L) == m[3L])
        checked <- checked + 1L
        if (nzchar(m[4L])) {
            same <- substr(gold$protein[i], as.integer(m[4L]), as.integer(m[4L]) + nchar(m[5L]) - 1L) == m[5L]
            if (gold$scenario[i] %in% alt_quirk) {
                if (!same) quirk_seen <- union(quirk_seen, gold$scenario[i])
            } else {
                stopifnot(same)
                alt_checked <- alt_checked + 1L
            }
        }
    }
}
stopifnot(setequal(quirk_seen, alt_quirk))

expected <- data.frame(case = gold$case, lane = gold$lane, csq_so_terms = unname(csq_so),
    declared_divergence = declared_text, stringsAsFactors = FALSE)
lines <- capture.output(write.table(expected, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
path <- paste0(root, 'frame_csq_expected.tsv')
if (nzchar(Sys.getenv('HAPLOTYPE_CSQ_REFRESH'))) writeLines(lines, path)
stopifnot(identical(readLines(path), lines))
decided <- gold$classifier_slice %in% 4:5 & gold$so_terms != ''
stopifnot(sum(decided) > 0L, all(csq_so[decided & !differs] == policy_so[decided & !differs]))
cat(sprintf(paste('%d lane goldens over %d layouts; %d ALT-carrying lanes mapped from csq; %d agree exactly,',
                  '%d in %d declared divergences; %d csq amino-acid changes checked (%d alternate sides)\n'),
            nrow(gold), length(unique(gold$layout)), sum(carrying), sum(carrying & !differs), sum(differs),
            length(declared), checked, alt_checked))
