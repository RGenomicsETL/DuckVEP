#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 1:2)
map <- read.delim(args[1], check.names = FALSE, na.strings = character())
stopifnot(!anyNA(map), !anyDuplicated(map$csq_term[map$kind %in% c('csq', 'reference')]),
          all(nzchar(map$csq_term[map$kind == 'csq'])),
          all(map$so_term[map$kind == 'vep_only'] == 'no_csq_equivalent'),
          all(map$so_term[map$kind == 'reference'] == 'unsupported_context'))
metadata <- readLines('src/kernel/src/duckvep_so_metadata.inc')
vep <- sub('.*\\{ "([^"]+)".*', '\\1', grep('= \\{ "', metadata, value = TRUE))
collapsed <- unlist(strsplit(map$collapsed_vep_terms[map$collapsed_vep_terms != '.'], ',', fixed = TRUE))
stopifnot(setequal(vep, collapsed), !anyDuplicated(collapsed),
          all(map$so_term[map$kind == 'csq'] %in% c(vep, 'unsupported_context')),
          sum(map$kind == 'csq') == 20L)
if (length(args) == 2L) {
    # bcftools query -f '%INFO/BCSQ\n' result: require complete coverage, not
    # a hand-selected subset of the observed vocabulary.
    input <- suppressWarnings(file(args[2], open = 'rt'))
    observed <- character()
    repeat {
        lines <- readLines(input, n = 10000L)
        if (!length(lines)) break
        fields <- unlist(strsplit(lines[lines != '.'], ',', fixed = TRUE))
        terms <- sub('\\|.*', '', fields)
        terms <- sub('^\\*', '', terms)
        observed <- union(observed, unlist(strsplit(terms[!grepl('^@|^$', terms)], '&', fixed = TRUE)))
    }
    close(input)
    missing <- setdiff(observed, map$csq_term[map$kind == 'csq'])
    if (length(missing)) stop('unmapped csq terms: ', paste(missing, collapse = ', '))
    cat('observed csq terms:', paste(sort(observed), collapse = ', '), '\n')
}
cat('csq→SO v1: 20 consequence terms; 41 VEP SO terms accounted for\n')
