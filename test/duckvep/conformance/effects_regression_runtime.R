#!/usr/bin/env Rscript

main <- function() {
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION"), mustWork = TRUE)
  con <- Rduckvep::rduckvep_connect(extension_path = extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  run <- function(name, reference, flags, table, changes) {
    n <- nchar(reference)
    padded <- bitwAnd(flags, 16L) != 0L
    end <- 99L + n - as.integer(padded)
    tx <- sprintf(paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
      "100::UBIGINT transcript_start,%d::UBIGINT transcript_end,1::TINYINT strand,",
      "0::UINTEGER gene_index,%d::UBIGINT transcript_flags,100::UBIGINT cds_start,",
      "%d::UBIGINT cds_end,%s::BLOB cds_sequence,%d::UTINYINT codon_table"),
      end, flags, end, quote(reference), table)
    ex <- sprintf(paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,",
      "%d::UBIGINT exon_end,1::UBIGINT exon_cdna_start,%d::UBIGINT exon_cdna_end,",
      "%d::TINYINT phase,0::TINYINT end_phase"),
      end, n - as.integer(padded), as.integer(padded))
    loaded <- DBI::dbGetQuery(con, sprintf("SELECT loaded FROM duckvep_model_load(%s,%s,%s,%s)",
      quote(name), quote("SELECT 0::UINTEGER seq_region"), quote(tx), quote(ex)))
    stopifnot(identical(loaded$loaded, TRUE))
    rows <- vapply(seq_len(nrow(changes)), function(i)
      sprintf("(%d,%d,%s,%s)", changes$id[i], changes$position[i] + 99L,
        quote(changes$ref[i]), quote(changes$alt[i])), character(1L))
    sql <- paste0("SELECT event_index,0 seq_region,position,reference,alternate,",
      "1 alt_index,0 transcript_index,0 sample_index,[1,1] alleles,",
      "[true,true] phase_before,NULL::BIGINT phase_set FROM (VALUES ",
      paste(rows, collapse = ","), ") v(event_index,position,reference,alternate)")
    result <- Rduckvep::rduckvep_haplotypes(con, sql, name)
    stopifnot(nrow(result) == 1L, identical(result$prediction_status, "predicted"))
    result
  }
  change <- function(position, ref, alt)
    data.frame(id = seq_along(position), position = position, ref = ref, alt = alt)
  check <- function(result, cds, protein, terms, impact) {
    stopifnot(identical(result$cds, cds), identical(result$protein, protein),
      setequal(result$haplotype_consequences[[1L]], terms),
      length(result$haplotype_consequences[[1L]]) == length(terms),
      identical(result$haplotype_impact, impact))
  }
  unknown <- run("effects_unknown_ntg", "NTGGCTGAATAA", 19L, 1L,
    change(1L, "T", "C"))
  check(unknown, "NCGGCTGAATAA", "XAE*", "coding_sequence_variant", "MODIFIER")
  stopifnot(identical(unknown$contributor_provenance[[1L]]$role, "applied"),
    identical(unknown$contributors[[1L]]$reference, "T"),
    identical(unknown$contributors[[1L]]$alternate, "C"))
  no_change <- tryCatch(run("effects_lower_noop", "atggctgaataa", 3L, 1L,
    change(1L, "a", "A")), error = identity)
  stopifnot(inherits(no_change, "error"),
    grepl("invalid input or candidate", conditionMessage(no_change), fixed = TRUE))
  changed <- run("effects_lower_changed", "atggctgaataa", 3L, 1L,
    change(6L, "t", "c"))
  check(changed, "ATGGCCGAATAA", "MAE*", "synonymous_variant", "LOW")
  stopifnot(identical(changed$contributors[[1L]]$reference, "t"),
    identical(changed$contributors[[1L]]$alternate, "c"))
  mixed <- run("effects_mixed_start", "atggctgaataa", 3L, 2L,
    change(3L, "g", "a"))
  check(mixed, "ATAGCTGAATAA", "MAE*",
    c("synonymous_variant", "start_retained_variant"), "LOW")
  stopifnot(identical(mixed$contributors[[1L]]$reference, "g"),
    identical(mixed$contributors[[1L]]$alternate, "a"))
  missense <- run("effects_lower_missense", "atggctgaataa", 3L, 1L,
    change(5L, "c", "t"))
  check(missense, "ATGGTTGAATAA", "MVE*", "missense_variant", "MODERATE")
  start <- run("effects_mito_start", "ATGGCTGAATAA", 3L, 2L,
    change(3L, "G", "A"))
  check(start, "ATAGCTGAATAA", "MAE*",
    c("synonymous_variant", "start_retained_variant"), "LOW")
  stopifnot(identical(start$contributor_provenance[[1L]]$role, "applied"),
    identical(start$contributors[[1L]]$reference, "G"),
    identical(start$contributors[[1L]]$alternate, "A"))
  overlap <- run("effects_lower_overlap", "atggctgaattttaa", 3L, 1L,
    change(c(4L, 5L), c("gc", "ct"), c("at", "tt")))
  check(overlap, "ATGATTGAATTTTAA", "MIEF*", "missense_variant", "MODERATE")
  stopifnot(identical(overlap$contributors[[1L]]$reference, c("gc", "ct")),
    identical(overlap$contributors[[1L]]$alternate, c("at", "tt")),
    identical(overlap$contributor_provenance[[1L]]$role, c("applied", "applied")))
  reference <- "ATGAGTCCCATCCGAAACTTCCTCCGTGATTCGGAATTACAAAAACATGTCTCGGGCCGTTGCGCCCCAAAAGCTATCCTCTCCCTACATGGCCCGTTGTGTAATTAA"
  replay <- function(reference, edits) {
    for (i in rev(seq_len(nrow(edits)))) {
      position <- edits$position[i]
      stopifnot(identical(substr(reference, position,
        position + nchar(edits$ref[i]) - 1L), edits$ref[i]))
      reference <- paste0(substr(reference, 1L, position - 1L), edits$alt[i],
        substr(reference, position + nchar(edits$ref[i]), nchar(reference)))
    }
    reference
  }
  for (lane in c("early", "late")) {
    edits <- if (lane == "early") change(c(17L, 20L), c("AC", "T"), c("A", "TT"))
      else change(c(82L, 86L), c("TC", "T"), c("T", "TT"))
    observed <- run(paste0("effects_frame_", lane), reference, 3L, 1L, edits)
    check(observed, replay(reference, edits), "MSPIRNFLRDSELQKHVSGRCAPKAILSLHGPLCN*",
      "synonymous_variant", "LOW")
    stopifnot(!identical(observed$cds, reference),
      identical(observed$contributor_provenance[[1L]]$role, c("applied", "applied")),
      identical(observed$contributors[[1L]]$reference, edits$ref))
  }
  reach <- run("effects_stop_reach", "ATGGCTGCTGAAGGTTAA", 3L, 1L,
    change(c(4L, 10L, 13L), c("G", "G", "G"), c("A", "T", "C")))
  check(reach, "ATGACTGCTTAACGTTAA", "MTA*",
    c("stop_gained", "missense_variant"), "HIGH")
  stopifnot(identical(reach$contributor_provenance[[1L]]$role,
    c("applied", "applied", "post_stop")))
  print(data.frame(case = c("unknown", "lower_changed", "mito_start"),
    cds = c(unknown$cds, changed$cds, start$cds),
    protein = c(unknown$protein, changed$protein, start$protein)))
}
main()
