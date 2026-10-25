#' Prepare VCF structural geometry and insertion provenance
#'
#' The input relation is shaped like DuckHTS `read_bcf()` output, one row per
#' event: `event_index`, `POS`, `REF`, `ALT` (a list, with the 1-based
#' `alt_index` choosing the allele, or one allele) and optionally `INFO_END`,
#' `INFO_SVLEN`, `INFO_CIPOS`, `INFO_CIEND` and `INFO_SEQ`, in any letter case
#' and as scalars, lists or text. A symbolic allele spans `POS + 1` to
#' `POS + |SVLEN|`, else `END`, as VEP reads it. Symbolic alleles retain
#' breakpoint confidence bounds separately from the nominal coordinates.
#' Literal anchored insertions carry inserted bases; a symbolic `INFO/SEQ` is
#' source provenance, not a literal insertion.
#'
#' @param con A connection with DuckVEP loaded.
#' @param input_table Name of the caller's input relation.
#' @param ... Named options for the native builder.
#' @return A data frame with one typed preparation row per input event.
#' @export
rduckvep_prepare_sv_geometry <- function(con, input_table, ...) {
  sql <- .duckvep_builder_sql(con, "duckvep_prepare_sv_geometry_sql",
                              list(input_table), list(...))
  DBI::dbGetQuery(con, paste0("SELECT * FROM query(",
                             as.character(DBI::dbQuoteString(con, sql)), ")"))
}
