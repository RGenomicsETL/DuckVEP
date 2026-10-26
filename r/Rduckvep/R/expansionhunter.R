#' Prepare ExpansionHunter v5 repeat alleles
#'
#' The input relation has one row per event and sample, as DuckHTS
#' `read_bcf(..., tidy_format = TRUE)` gives them: `event_index`, `REF`, `ALT`,
#' the one-based `alt_index`, and optionally `INFO_END`, `INFO_REF`, `INFO_RL`,
#' `INFO_RU`, `FORMAT_GT`, `FORMAT_SO`, and either `FORMAT_REPCN` with
#' `FORMAT_REPCI` or `FORMAT_CN` with `FORMAT_CI`. The reference relation has `event_index`
#' and the literal `reference_sequence` from POS+1 through END. Each input
#' event retains its source, preparation reason and exactness decision.
#'
#' @param con A connection with DuckVEP loaded.
#' @param input_table Name of the caller's input relation.
#' @param reference_table Name of the literal reference relation.
#' @param ... Named options for the native builder.
#' @return A data frame with typed preparation rows.
#' @export
rduckvep_prepare_expansionhunter <- function(con, input_table, reference_table, ...) {
  sql <- .duckvep_builder_sql(con, "duckvep_prepare_expansionhunter_sql",
                              list(input_table, reference_table), list(...))
  DBI::dbGetQuery(con, paste0("SELECT * FROM query(",
                             as.character(DBI::dbQuoteString(con, sql)), ")"))
}
