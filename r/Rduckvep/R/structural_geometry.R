#' Prepare VCF structural geometry and insertion provenance
#'
#' The input relation has `event_index`, `pos`, `ref`, `alt`, and `info`.
#' Symbolic alleles retain breakpoint confidence bounds separately from the
#' nominal coordinates. Literal anchored insertions carry inserted bases;
#' symbolic INFO/SEQ is source provenance, not a literal insertion.
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
