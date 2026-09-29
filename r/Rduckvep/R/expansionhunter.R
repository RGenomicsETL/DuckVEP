#' Prepare ExpansionHunter v5 repeat alleles
#'
#' The input relation has `event_index`, `info`, `format`, `sample`, `ref`,
#' `alt`, and one-based `alt_index`. The reference relation has `event_index`
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
