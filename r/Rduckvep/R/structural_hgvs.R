#' Prepare structural HGVS descriptions
#'
#' The events relation is shaped like DuckHTS `read_bcf()` output:
#' `event_index`, `CHROM`, `POS`, `REF`, `ALT` and optionally `INFO_END`,
#' `INFO_SVLEN`, `INFO_CIPOS`, `INFO_CIEND` and `INFO_IMPRECISE`; the reference
#' relation has `event_index` and `reference_sequence` from `POS` through the
#' end, which is `POS + |SVLEN|`, else `END`. Exact-span symbolic DEL, DUP, DUP:TANDEM and INV
#' alleles receive an unshifted genomic description and an equivalent literal
#' allele for transcript HGVS. Everything else is explicitly `unavailable` or
#' `unsupported`, with a stable reason. VEP 116 itself emits no HGVS for
#' symbolic structural alleles or breakends.
#'
#' @param con A connection with DuckVEP loaded.
#' @param events_table Name of the caller's VCF record relation.
#' @param reference_table Name of the reference sequence relation.
#' @param ... Named options for the native builder (`max_span`).
#' @return A data frame with one row per input event.
#' @export
rduckvep_prepare_structural_hgvs <- function(con, events_table, reference_table, ...) {
  sql <- .duckvep_builder_sql(con, "duckvep_prepare_structural_hgvs_sql",
                              list(events_table, reference_table), list(...))
  DBI::dbGetQuery(con, paste0("SELECT * FROM query(",
                             as.character(DBI::dbQuoteString(con, sql)), ")"))
}
