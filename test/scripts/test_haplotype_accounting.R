#!/usr/bin/env Rscript
scratch <- tempfile('csq-account-')
dir.create(scratch)
line <- function(pos, annotation, gt = '1|0', ps = 'A', ref = 'C', alt = 'T',
                 lane0 = annotation, lane1 = '.')
  paste('1', pos, ref, alt, annotation, gt, ps, lane0, lane1, sep = '\t')
coding <- function(tx) paste0('missense|G|', tx, '|protein_coding|+|2Q>2*|4C>T')
a <- c(line(1, coding('ENST123')), line(2, 'intron|G||lncRNA'),
       line(3, coding('ENST124')), line(4, coding('ENST123'), '0/1', '.'),
       line(5, '@1'), line(6, coding('ENST123'), ref = 'N'),
       line(7, '.', ps = 'A'), line(8, coding('ENST125'), ps = 'A'),
       line(9, coding('ENST125'), ps = 'B'), line(10, coding('ENST126')),
       line(10, coding('ENST126'), alt = 'G'))
b <- a
b[7L] <- line(7, coding('ENST123'))
a[4L] <- line(4, coding('ENST123'), '0/1', '.', lane0 = '.')
strict <- file.path(scratch, 'strict.tsv'); domain <- file.path(scratch, 'domain.tsv')
parity <- file.path(scratch, 'parity.tsv'); counts <- file.path(scratch, 'counts.tsv')
writeLines(a, strict); writeLines(b, domain)
writeLines(c('tx\tparity', 'ENST123\tidentical', 'ENST124\tgeometry_mismatch',
             'ENST125\tidentical', 'ENST126\tidentical'), parity)
script <- 'scripts/account_haplotype_csq.R'
run <- function() suppressWarnings(system2('Rscript', c(script, strict, domain, parity, counts),
                                             stdout = TRUE, stderr = TRUE))
stopifnot(is.null(attr(run(), 'status')))
observed <- setNames(read.delim(counts)$count, read.delim(counts)$category)
stopifnot(identical(unname(observed[c('compared', 'outside_csq_domain', 'model_disagreement',
                                      'phase_disagreement')]), c(1L, 5L, 1L, 4L)),
          length(readLines(gzfile(paste0(counts, '.records.tsv.gz')))) == length(a) + 1L)
writeLines(rev(b), domain)
stopifnot(!is.null(attr(run(), 'status')))
writeLines(b[-11L], domain)
stopifnot(!is.null(attr(run(), 'status')))
unlink(scratch, recursive = TRUE)
cat('accounting partitions all records; reordered and missing streams fail closed\n')
