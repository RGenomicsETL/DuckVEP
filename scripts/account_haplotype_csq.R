#!/usr/bin/env Rscript
# Inputs are positional bcftools query streams from the same VCF, with
# -p s and -p a respectively. The a stream detects phase-hidden coding
# overlap; only the s stream can authorize a comparison.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
strict <- args[1]
domain <- args[2]
parity <- read.delim(args[3], check.names = FALSE)
model <- setNames(parity$parity, parity$tx)
coding <- '(synonymous|missense|coding_sequence|inframe_deletion|inframe_insertion|inframe_altering|frameshift|stop_gained|stop_lost|start_lost|stop_retained|start_retained)'
# BCSQ's compound references (@position and *term) cannot be classified
# independently of their linked source records.
tx_for <- function(info) {
    if (!grepl('|protein_coding|', info, fixed = TRUE)) return(character())
    fields <- strsplit(info, ',', fixed = TRUE)[[1L]]
    fields <- fields[grepl('|protein_coding|', fields, fixed = TRUE) &
                     grepl(paste0('(^|&)', coding), sub('^\\*', '', fields)) &
                     !grepl('^@|^\\*', fields)]
    if (!length(fields)) return(character())
    parts <- strsplit(fields, '|', fixed = TRUE)
    ids <- unique(vapply(parts, function(x) if (length(x) >= 4L &&
      grepl('^ENST[0-9]+$', x[3L])) x[3L] else '', character(1L)))
    ids[nzchar(ids)]
}
# PS is relevant for heterozygous calls; homozygotes have unambiguous lanes.
het <- function(gt) {
    alleles <- strsplit(gt, '[/|]')[[1L]]
    length(alleles) == 2L && all(grepl('^[0-9]+$', alleles)) && alleles[1L] != alleles[2L]
}
phase <- new.env(hash = TRUE, parent = emptyenv())
last_end <- new.env(hash = TRUE, parent = emptyenv())
overlap <- new.env(hash = TRUE, parent = emptyenv())
read_pairs <- function(fun) {
    a <- file(strict, 'rt'); b <- file(domain, 'rt')
    on.exit({close(a); close(b)})
    repeat {
        sa <- readLines(a, n = 2000L); sb <- readLines(b, n = 2000L)
        if (length(sa) != length(sb)) stop('streams have different record counts')
        if (!length(sa)) break
        fun(sa, sb)
    }
}
read_pairs(function(sa, sb) {
    records <- strsplit(sb, '\t', fixed = TRUE)
    for (r in records) {
        if (length(r) < 7L || length(r) > 9L) stop('malformed domain record')
        r <- c(r, rep('', 9L - length(r)))
        lanes <- lapply(r[8:9], tx_for)
        ids <- unique(unlist(lanes, use.names = FALSE))
        if (!length(ids)) next
        right <- as.integer(r[2L]) + nchar(r[3L]) - 1L
        for (lane in seq_along(lanes)) for (id in lanes[[lane]]) {
            key <- paste(id, lane, sep = ':')
            previous <- last_end[[key]]
            if (!is.null(previous) && as.integer(r[2L]) <= previous) overlap[[id]] <- TRUE
            last_end[[key]] <- if (is.null(previous)) right else max(previous, right)
        }
        if (!het(r[6L]) || !grepl('|', r[6L], fixed = TRUE)) next
        ps <- if (r[7L] == '.') 'implicit' else r[7L]
        for (id in ids) {
            previous <- phase[[id]]
            if (is.null(previous)) phase[[id]] <- ps
            else if (previous != ps) phase[[id]] <- 'conflict'
        }
    }
})
counts <- setNames(integer(4), c('compared', 'outside_csq_domain', 'model_disagreement', 'phase_disagreement'))
reasons <- new.env(hash = TRUE, parent = emptyenv())
ledger <- gzfile(paste0(args[4], '.records.tsv.gz'), 'wt')
writeLines('record_index\tchrom\tpos\tref\talt\tcategory\treason', ledger)
ordinal <- 0L
read_pairs(function(sa, sb) {
    x <- strsplit(sa, '\t', fixed = TRUE)
    y <- strsplit(sb, '\t', fixed = TRUE)
    lines <- character(length(x))
    for (i in seq_along(x)) {
        s <- x[[i]]; d <- y[[i]]
        if (length(s) < 7L || length(s) > 9L || length(d) < 7L || length(d) > 9L ||
            !identical(s[1:4], d[1:4]) || !identical(s[6:7], d[6:7]))
            stop('strict/domain stream drift')
        s <- c(s, rep('', 9L - length(s)))
        d <- c(d, rep('', 9L - length(d)))
        ids <- tx_for(paste(d[8:9], collapse = ',')); gt <- d[6L]
        eligible <- grepl('^[ACGT]+$', d[3L]) && nchar(d[3L]) <= 50L &&
                    all(grepl('^[ACGT]+$', strsplit(d[4L], ',', fixed = TRUE)[[1L]])) &&
                    all(nchar(strsplit(d[4L], ',', fixed = TRUE)[[1L]]) <= 50L)
        if (!eligible || !length(ids) || any(vapply(ids, function(id) isTRUE(overlap[[id]]), logical(1L)))) {
            label <- 'outside_csq_domain'
            reason <- if (!eligible) 'unsupported_allele' else if (!length(ids))
              if (grepl('^@', d[5L])) 'compound_reference_only' else 'no_carried_coding_consequence'
              else 'overlapping_coding_edits'
        } else if (any(is.na(model[ids]) | model[ids] != 'identical')) {
            label <- 'model_disagreement'; reason <- 'exon_or_CDS_geometry'
        } else if (!grepl('^[0-9]+[|/][0-9]+$', gt) ||
                   (het(gt) && (!grepl('|', gt, fixed = TRUE) ||
                     any(vapply(ids, function(id) identical(phase[[id]], 'conflict'), logical(1L)))))) {
            label <- 'phase_disagreement'; reason <- 'ambiguous_lane_or_phase_set'
        } else if (!setequal(ids, tx_for(paste(s[8:9], collapse = ',')))) {
            label <- 'phase_disagreement'; reason <- 'strict_csq_coding_set_differs'
        } else {
            label <- 'compared'; reason <- 'coding_parity_and_phase'
        }
        ordinal <<- ordinal + 1L
        lines[i] <- paste(ordinal, paste(d[1:4], collapse = '\t'), label, reason, sep = '\t')
        counts[label] <<- counts[label] + 1L
        old <- reasons[[reason]]
        reasons[[reason]] <- if (is.null(old)) 1L else old + 1L
    }
    writeLines(lines, ledger)
})
close(ledger)
out <- data.frame(category = names(counts), count = unname(counts))
write.table(out, args[4], sep = '\t', row.names = FALSE, quote = FALSE)
cat('records:', sum(counts), '\n')
print(out, row.names = FALSE)
cat('reasons:\n')
print(sort(unlist(as.list(reasons)), decreasing = TRUE))
