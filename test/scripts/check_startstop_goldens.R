#!/usr/bin/env Rscript
# Start and stop (coding-v1 slice 5) contract check, offline by default.
#  1. regenerate the base-R fixtures/goldens (fails if a committed fixture differs);
#  2. verify the pinned csq output against its receipt hashes, and re-run csq only when the pinned
#     bcftools is installed;
#  3. map the csq lane consequences through csq_so_map_v1.tsv and require exact agreement with the
#     goldens apart from the divergences declared below and in test/data/haplotype/README.md;
#  4. check every csq amino-acid change against the base-R peptides.
# HAPLOTYPE_CSQ_REFRESH=1 rewrites startstop_csq.tsv, its receipt and startstop_csq_expected.tsv (needs the pinned bcftools).
stopifnot(system2('Rscript', c('--vanilla', 'test/scripts/generate_startstop_goldens.R')) == 0L)
root <- 'test/data/haplotype/'
sha <- function(path) substr(system2('sha256sum', path, stdout = TRUE), 1L, 64L)
csq <- Sys.getenv('HAPLOTYPE_CSQ', '/usr/local/bin/bcftools')
pinned_version <- '1.23.1-70-g6dbd8fef+htslib-1.22.1'
run_csq <- function() {
    vcf <- tempfile(fileext = '.vcf'); on.exit(unlink(vcf))
    stopifnot(system2(csq, c('csq', '-p', 's', '-n', '1024', '-f', paste0(root, 'startstop.fa'),
        '-g', paste0(root, 'startstop.gff3'), '-Ov', '-o', vcf, paste0(root, 'startstop.vcf')),
        stderr = tempfile()) == 0L)
    system2(csq, c('query', '-f', shQuote('%CHROM\t%POS\t%INFO/BCSQ\t[%TBCSQ]\n'), vcf), stdout = TRUE)
}
have_csq <- file.exists(csq) && identical(system2(csq, '--version-only', stdout = TRUE), pinned_version)
if (nzchar(Sys.getenv('HAPLOTYPE_CSQ_REFRESH'))) {
    stopifnot(have_csq)
    writeLines(run_csq(), paste0(root, 'startstop_csq.tsv'))
    key <- c(bcftools_commit = '6dbd8fef51e529755a4a81544075dc6a43ff3cd3',
             htslib_commit = 'c1f35d67dd5ff1e226d94abe7b850f28c60f5910',
             bcftools_version = pinned_version, htslib_version = '1.22.1', phase_mode = 's', ncsq = '1024',
             input_sha256 = sha(paste0(root, 'startstop.vcf')),
             reference_sha256 = sha(paste0(root, 'startstop.fa')),
             gff3_sha256 = sha(paste0(root, 'startstop.gff3')),
             output_sha256 = sha(paste0(root, 'startstop_csq.tsv')))
    writeLines(c('field\tvalue', paste(names(key), key, sep = '\t')), paste0(root, 'startstop_csq.oracle.tsv'))
}
receipt <- read.delim(paste0(root, 'startstop_csq.oracle.tsv'), check.names = FALSE)
pin <- setNames(receipt$value, receipt$field)
files <- c(input = 'startstop.vcf', reference = 'startstop.fa', gff3 = 'startstop.gff3', output = 'startstop_csq.tsv')
for (k in names(files)) stopifnot(identical(sha(paste0(root, files[[k]])), pin[[paste0(k, '_sha256')]]))
stopifnot(pin[['bcftools_version']] == pinned_version, pin[['phase_mode']] == 's', pin[['ncsq']] == '1024')
pinned <- readLines(paste0(root, 'startstop_csq.tsv'))
if (have_csq) {
    stopifnot(identical(run_csq(), pinned))
} else message('pinned bcftools not found; checking the goldens against the committed csq output')

map <- read.delim(paste0(root, 'csq_so_map_v1.tsv'), check.names = FALSE, na.strings = character())
lookup <- setNames(map$so_term[map$kind == 'csq'], map$csq_term[map$kind == 'csq'])
gold <- read.delim(paste0(root, 'startstop_goldens.tsv'), na.strings = character(), check.names = FALSE,
                   colClasses = c(lane = 'integer'))
transcripts <- read.delim(paste0(root, 'startstop_transcripts.tsv'), stringsAsFactors = FALSE)
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

# Declared divergences between the reduced whole-protein policy and csq, keyed by scenario and lane (a
# `name:tag` key gives a second declaration for the same scenario). `csq` is the asserted mapped csq set, one
# value or c(plus =, minus =) when csq depends on the strand representation; `layouts` restricts a declaration to
# the layouts matching that regex. Every lane whose sets differ must be covered, and every declared lane must
# differ, so a change in the pinned csq output cannot pass as the same divergence.
why <- c(
  suppressed = paste('start loss suppresses other biological predictions (contract section 2), so the policy',
                     'reports start_lost alone; csq also reports the other effect'),
  start_retained = paste('csq adds start_retained_variant (LOW, outside the reduced set) when the start codon is',
                         'kept; the policy reports the whole-protein category only'),
  start_deleted = paste('csq calls a deletion that removes bases 2-4 of the start codon start_retained_variant',
                        'although the edited CDS begins AGT; the policy classifies the edited sequence: start_lost'),
  start_anchor = paste('csq anchors an insertion after the last base of the start codon inside the codon: it',
                       'reports start_lost on the plus strand; the policy uses the normalized edit (an insertion',
                       'between codons 1 and 2 leaves the start codon intact)'),
  single_frameshift = paste('csq reports a single frameshift record as frameshift only and does not add the',
                            'stop_gained of the early stop it creates; the policy reports both HIGH terms'),
  run_off = paste('csq reports a frame that runs off the CDS with no stop as frameshift only; the contract has',
                  'the reference termination abolished as well, so the policy reports frameshift_variant and stop_lost'),
  run_off_plus = paste('csq reports the run-off frame as frameshift only on the plus strand and as',
                       'stop_lost&frameshift on the minus strand: its indel anchor differs by strand (the contract',
                       'frame SO is normalized-edit-path-sensitive); the policy reports both on every layout'),
  in_frame_edit = paste('csq also reports the in-frame insertion or deletion that carries the abolished',
                        'termination; the reduced set reports the HIGH stop_lost only'),
  lower = paste('csq also reports the lower-severity effect (missense, synonymous or in-frame) of another edit;',
                'the reduced set reports the one highest whole-protein category'),
  restored = paste('csq adds the per-record lower-severity missense to inframe_altering when the frame is',
                   'restored; the policy reports the one whole-protein category protein_altering_variant'),
  retained_added = paste('csq adds stop_retained_variant for the terminal codon next to an expressed effect; the',
                         'policy reports the expressed effect only (stop_retained is LOW and lower classes give way)'),
  post_stop = paste('csq reports an edit after the first stop as its own effect; the policy keeps it as a',
                    'post_stop contributor, not an expressed effect (contract section 2)'),
  inserted_stop = paste('a stop codon inserted next to the terminator leaves the peptide unchanged and is',
                        'sequence-equivalent to a stop inserted after it: the policy reports stop_retained_variant;',
                        'csq reports stop_gained on the plus strand and stop_retained on the minus strand'),
  anchor_terminal = paste('csq anchors an insertion before the terminator inside the terminator codon on the minus',
                          'strand and adds stop_retained_variant; the policy uses the normalized edit'),
  rescued_terminator_deletion = paste('the normalized deletion of one base of TT (cDNA 105-106) has several placements;',
                                      'csq reports a restored frame on the plus strand and misses the stop the +1 frame',
                                      'reads at cDNA 102, and adds stop_retained on the minus strand; the policy reports',
                                      'the stop gained in the displaced frame'))
d <- function(lane, csq, reason, layouts = NULL) list(lane = lane, csq = csq, why = why[[reason]], layouts = layouts)
declared <- list(
  start_deletion_within_exon = d(0L, 'inframe_deletion,start_lost', 'suppressed'),
  start_lost_suppresses_frame = d(0L, 'frameshift_variant,start_lost', 'suppressed'),
  start_lost_suppresses_missense = d(0L, 'missense_variant,start_lost', 'suppressed'),
  start_lost_suppresses_stop_gain = d(0L, 'start_lost,stop_gained', 'suppressed'),
  start_lost_suppresses_stop_lost = d(0L, 'start_lost,stop_lost', 'suppressed'),
  post_stop_after_start_lost = d(0L, 'start_lost,stop_gained,stop_lost', 'suppressed'),
  start_rescued_by_insertion_a = d(0L, 'inframe_insertion,start_retained_variant', 'start_retained'),
  start_rescued_by_insertion_b = d(0L, 'inframe_insertion,start_retained_variant', 'start_retained'),
  start_deletion_inframe = d(0L, 'start_retained_variant', 'start_deleted'),
  start_insertion_after_start = d(0L, 'inframe_insertion,start_lost', 'start_anchor', '^plus'),
  start_intact_frame_open_early_stop = d(0L, 'start_lost', 'start_anchor', '^plus'),
  `start_intact_frame_open_early_stop:minus` = d(0L, 'frameshift_variant', 'single_frameshift', '^minus'),
  run_off_deletion = d(0L, 'frameshift_variant', 'run_off'),
  run_off_insertion_plus2 = d(0L, 'frameshift_variant', 'run_off'),
  run_off_inside_codon_35 = d(0L, 'frameshift_variant', 'run_off'),
  run_off_trans_control = d(0L, 'frameshift_variant', 'run_off'),
  stop_lost_deletion_partial = d(0L, 'frameshift_variant', 'run_off_plus', '^plus'),
  stop_lost_insertion_before_plus1 = d(0L, 'frameshift_variant', 'run_off_plus', '^plus'),
  stop_lost_insertion_before_plus2 = d(0L, 'frameshift_variant', 'run_off_plus', '^plus'),
  stop_lost_terminator_deleted = d(0L, 'inframe_deletion,stop_lost', 'in_frame_edit'),
  stop_lost_last_codon_and_terminator_deleted = d(0L, 'inframe_deletion,stop_lost', 'in_frame_edit'),
  stop_lost_insertion_inside = d(0L, 'inframe_insertion,stop_lost', 'in_frame_edit'),
  stop_lost_and_missense = d(0L, 'missense_variant,stop_lost', 'lower'),
  stop_retained_plus_synonymous = d(0L, 'stop_retained_variant,synonymous_variant', 'lower'),
  stop_retained_plus_missense = d(0L, 'missense_variant,stop_retained_variant', 'retained_added'),
  stop_retained_plus_inframe_insertion = d(0L, 'inframe_insertion,stop_retained_variant', 'retained_added'),
  stop_retained_after_frame_pair = d(0L, 'missense_variant,protein_altering_variant,stop_retained_variant', 'retained_added'),
  stop_retained_displaced_insertion = d(0L, 'frameshift_variant,stop_retained_variant', 'retained_added'),
  stop_inserted_before_terminator = d(0L, c(plus = 'inframe_insertion,stop_gained',
    minus = 'inframe_insertion,stop_retained_variant'), 'inserted_stop'),
  inframe_insertion_before_terminator = d(0L, 'inframe_insertion,stop_retained_variant', 'anchor_terminal', '^minus'),
  rescued_termination_after_stop_lost = d(0L, 'missense_variant,protein_altering_variant,stop_lost', 'lower'),
  rescued_termination_restore_inside_terminator = d(0L, 'protein_altering_variant,stop_lost', 'lower'),
  rescued_termination_by_restoring_insertion = d(0L, 'missense_variant,protein_altering_variant', 'restored'),
  rescued_termination_by_deletion_in_terminator = d(0L, 'frameshift_variant,stop_gained,stop_retained_variant', 'retained_added'),
  rescued_termination_by_terminator_deletion = d(0L, c(plus = 'missense_variant,protein_altering_variant',
    minus = 'frameshift_variant,stop_gained,stop_retained_variant'), 'rescued_terminator_deletion'),
  post_stop_terminator_lost = d(0L, 'stop_gained,stop_lost', 'post_stop'),
  post_stop_terminator_retained = d(0L, 'stop_gained,stop_lost', 'post_stop'),
  post_stop_run_off = d(0L, 'frameshift_variant,stop_gained,stop_retained_variant', 'post_stop', '^minus'))
scenario_of <- function(name) sub(':.*$', '', name)
declared_csq <- function(dec, layout) if (is.null(names(dec$csq))) rep(dec$csq, length(layout)) else
    unname(dec$csq[sub('_.*$', '', layout)])
declared_rows <- function(name, dec) gold$scenario == scenario_of(name) & gold$lane %in% dec$lane &
    (if (is.null(dec$layouts)) TRUE else grepl(dec$layouts, gold$layout))
differs <- csq_so != policy_so
carrying <- policy_so != '' | csq_so != ''
divergent <- unique(gold$scenario[differs])
stopifnot(setequal(divergent, unique(scenario_of(names(declared)))))
declared_text <- rep('.', nrow(gold)); covered <- rep(FALSE, nrow(gold))
for (name in names(declared)) {
    dec <- declared[[name]]
    lanes <- declared_rows(name, dec)
    ok <- csq_so[lanes] == declared_csq(dec, gold$layout)[lanes]
    if (!all(ok)) stop('declaration ', name, ': ', paste(unique(paste(gold$layout, csq_so)[lanes][!ok]), collapse = '; '))
    stopifnot(any(lanes), all(differs[lanes]), !any(covered & lanes))
    declared_text[lanes] <- dec$why
    covered <- covered | lanes
}
stopifnot(all(covered[differs]), !any(covered & !differs))

# Every csq amino-acid change must describe the base-R peptides (reference at its position, alternate in
# the edited protein), so the SO comparison cannot pass on a mislabelled change. Declared csq quirk: for
# a frameshifting deletion and for some compound records csq's alternate string differs from the edited
# protein at the residues next to the deletion (e.g. 20L for an edited residue 20R). The SO term is
# right; only the alternate side of the scenarios below is skipped, their reference side is still checked
# and each listed scenario must really show the mismatch, so the list cannot go stale.
alt_quirk <- c('stop_lost_deletion_partial', 'inframe_deletion_before_terminator',
               'rescued_termination_by_terminator_deletion')
unknown_seen <- 0L
checked <- 0L; alt_checked <- 0L; quirk_seen <- character()
for (i in seq_len(nrow(gold))) {
    aas <- lane_aa[[keys[i]]]
    for (j in seq_along(aas)) {
        m <- regmatches(aas[j], regexec('^([0-9]+)([A-Z*]+)(?:>([0-9]+)([A-Z*?]+))?$', aas[j]))[[1L]]
        if (length(m) < 3L) stop('unparsed csq amino-acid change: ', aas[j], ' for ', keys[i])
        pos <- as.integer(m[2L])
        stopifnot(substr(ref_protein, pos, pos + nchar(m[3L]) - 1L) == m[3L])
        checked <- checked + 1L
        if (nzchar(m[4L]) && m[5L] == '?') {
            # csq leaves the alternate residue unknown where the stop was abolished and reads nothing past the
            # CDS: the policy protein has no residue there either (no downstream extension is invented).
            stopifnot(nchar(gold$protein[i]) < as.integer(m[4L]), gold$protein[i] != '')
            unknown_seen <- unknown_seen + 1L
        } else if (nzchar(m[4L])) {
            same <- substr(gold$protein[i], as.integer(m[4L]), as.integer(m[4L]) + nchar(m[5L]) - 1L) == m[5L]
            if (gold$scenario[i] %in% alt_quirk) {
                if (!same) quirk_seen <- union(quirk_seen, gold$scenario[i])
            } else {
                if (!same) stop('alt aa mismatch: ', gold$case[i], ' ', aas[j], ' vs ', gold$protein[i])
                alt_checked <- alt_checked + 1L
            }
        }
    }
}
stopifnot(setequal(quirk_seen, alt_quirk))

expected <- data.frame(case = gold$case, lane = gold$lane, csq_so_terms = unname(csq_so),
    declared_divergence = declared_text, stringsAsFactors = FALSE)
lines <- capture.output(write.table(expected, sep = '\t', row.names = FALSE, quote = FALSE, na = '.'))
path <- paste0(root, 'startstop_csq_expected.tsv')
if (nzchar(Sys.getenv('HAPLOTYPE_CSQ_REFRESH'))) writeLines(lines, path)
stopifnot(identical(readLines(path), lines))
decided <- gold$classifier_slice %in% 4:5 & gold$so_terms != ''
stopifnot(sum(decided) > 0L, all(csq_so[decided & !differs] == policy_so[decided & !differs]))
cat(sprintf(paste('%d lane goldens over %d layouts; %d ALT-carrying lanes mapped from csq; %d agree exactly,',
                  '%d in %d declared divergences; %d csq amino-acid changes checked (%d alternate sides, %d unknown-residue ?)\n'),
            nrow(gold), length(unique(gold$layout)), sum(carrying), sum(carrying & !differs), sum(differs),
            length(declared), checked, alt_checked, unknown_seen))
