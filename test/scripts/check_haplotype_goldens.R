#!/usr/bin/env Rscript
# Independent base-R policy check and pinned executable csq fixture check.
source_file <- 'test/scripts/generate_haplotype_goldens.R'
stopifnot(system2('Rscript', source_file) == 0L)
stopifnot(system2('Rscript', c('test/scripts/check_haplotype_csq_map.R',
                              'test/data/haplotype/csq_so_map_v1.tsv')) == 0L)
root <- 'test/data/haplotype/'
receipt <- read.delim(paste0(root, 'vertical_csq.oracle.tsv'), check.names = FALSE)
pin <- setNames(receipt$value, receipt$field)
files <- c(input = 'vertical.vcf', reference = 'vertical.fa',
           gff3 = 'vertical.gff3', output = 'vertical_csq.tsv')
for (key in names(files)) {
    actual <- substr(system2('sha256sum', paste0(root, files[[key]]), stdout = TRUE), 1L, 64L)
    stopifnot(identical(actual, pin[[paste0(key, '_sha256')]]))
}
csq <- '/usr/local/bin/bcftools'
stopifnot(system2(csq, '--version-only', stdout = TRUE) == pin[['bcftools_version']],
          pin[['phase_mode']] == 's', pin[['ncsq']] == '345')
vcf <- tempfile(fileext = '.vcf')
status <- system2(csq, c('csq', '-p', 's', '-n', '345', '-f', paste0(root, 'vertical.fa'),
                         '-g', paste0(root, 'vertical.gff3'), '-Ov', '-o', vcf,
                         paste0(root, 'vertical.vcf')), stderr = tempfile())
stopifnot(status == 0L)
observed <- system2(csq, c('query', '-f', shQuote('%CHROM\t%POS\t%INFO/BCSQ\t[%TBCSQ]\n'),
                          vcf), stdout = TRUE)
unlink(vcf)
stopifnot(identical(observed, readLines(paste0(root, 'vertical_csq.tsv'))))
map <- read.delim(paste0(root, 'csq_so_map_v1.tsv'), check.names = FALSE)
lookup <- setNames(map$so_term[map$kind == 'csq'], map$csq_term[map$kind == 'csq'])
expected_csq <- c(cis = 'synonymous', trans = 'missense', frame_open = 'frameshift',
  restored_before_stop = 'inframe_altering,missense',
  restored_after_stop = 'frameshift,stop_gained', start_lost = 'start_lost',
  stop_lost = 'stop_lost', stop_retained = 'stop_retained',
  nmd49 = 'stop_gained', nmd50 = 'stop_gained', nmd51 = 'stop_gained',
  intronless = 'stop_gained', shifted_junction = 'stop_gained')
rows <- strsplit(observed, '\t', fixed = TRUE)
terms <- setNames(vector('list', length(expected_csq)), names(expected_csq))
for (row in rows) {
    stopifnot(length(row) == 5L)
    if (!startsWith(row[4L], '@') && row[4L] != '.') {
        # FORMAT/TBCSQ contains the left (lane 0) consequence for this fixture.
        first <- strsplit(row[4L], '|', fixed = TRUE)[[1L]][1L]
        if (!startsWith(first, '*')) terms[[row[1L]]] <- union(terms[[row[1L]]],
          strsplit(first, '&', fixed = TRUE)[[1L]])
    }
}
for (name in names(expected_csq))
    stopifnot(identical(sort(terms[[name]]), sort(strsplit(expected_csq[[name]], ',', fixed = TRUE)[[1L]])))
lane1 <- vapply(rows, function(row) row[5L], character(1L))
stopifnot(sum(lane1 != '.') == 1L,
          startsWith(lane1[which(lane1 != '.')], 'missense|G_trans|T_trans|'))
stopifnot(all(unlist(terms) %in% names(lookup)),
          any(grepl('^shifted_junction\t60\t\\*inframe_insertion', observed)),
          any(grepl('^cis\t5\t@4\t@4\t\\.$', observed)),
          any(grepl('2S>2C', observed, fixed = TRUE)),
          any(grepl('2S>2T', observed, fixed = TRUE)))
gold <- read.delim(paste0(root, 'vertical_goldens.tsv'), na.strings = character(),
                   check.names = FALSE)
stopifnot(nrow(gold) == 26L,
  all(gold$so_terms[gold$lane == 1L & gold$case != 'trans'] == ''),
  identical(gold$nmd[gold$case %in% c('nmd49', 'nmd50', 'nmd51') & gold$lane == 0L],
            c('escape', 'escape', 'trigger')))
csq_so <- vapply(terms, function(x) paste(sort(unique(unname(lookup[x]))), collapse = ','), character(1L))
policy_so <- setNames(gold$so_terms[gold$lane == 0L], gold$case[gold$lane == 0L])
policy_so <- vapply(policy_so, function(x) paste(sort(strsplit(x, ',', fixed = TRUE)[[1L]]), collapse = ','), character(1L))
stopifnot(setequal(names(csq_so), names(policy_so)))
differences <- names(csq_so)[csq_so != policy_so]
stopifnot(setequal(differences, c('frame_open', 'restored_before_stop')))
cat('26 base-R lane goldens; 18 pinned csq record strings; 2 declared reduced-SO divergences\n')
