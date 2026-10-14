#!/usr/bin/env Rscript
# Full-codon indels in the final edited molecule are expressed only before its first stop.
source('test/scripts/generate_nmd_goldens.R', local = TRUE)

expected <- c(up_ins3_d51 = 'inframe_insertion', up_ins3_d50 = 'inframe_insertion',
              up_del3_d50 = 'inframe_deletion', up_del3_d51 = 'inframe_deletion',
              up_ins6_three_exons = 'inframe_insertion',
              up_indel_in_pen_exon = 'inframe_insertion',
              up_del_in_first_exon_of_three = 'inframe_deletion')
for (name in names(expected)) {
    scenario <- scenarios[[name]]
    cds <- cdsets[[scenario$cds]]
    edits <- scenario$edits(E(cds))
    w <- edit_mrna(cds, rep(1L, nchar(cds)), edits, 0L, 0L)
    protein <- translate(w$seq)
    first <- match('*', protein)
    atom <- island(edits[[1L]])
    stopifnot(!is.na(first), sub_at(w$seq, 3L * first - 2L, 3L) == 'TAA',
              atom$start < 3L * first - 2L,
              edits[[2L]]$pos == 49L,
              identical(protein[first], '*'))
    if (expected[[name]] == 'inframe_insertion') {
        stopifnot(atom$ref_len == 0L, nchar(atom$alt) %% 3L == 0L,
                  atom$start %% 3L == 1L,
                  sub_at(w$seq, atom$start, nchar(atom$alt)) == atom$alt,
                  all(translate(atom$alt) == 'G'))
    } else {
        stopifnot(atom$ref_len == 3L, atom$alt == '',
                  atom$start %% 3L == 1L,
                  sub_at(cds, atom$start, 3L) == 'CAA')
    }
    lanes <- golden[golden$scenario == name & golden$has_carrier, ]
    stopifnot(nrow(lanes) == 4L, all(lanes$premature))
}
cat('28 mapped in-frame indel lanes precede the independently translated TAA stop\n')
