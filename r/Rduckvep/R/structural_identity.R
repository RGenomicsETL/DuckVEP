#' Validate BND mate and event identity
#'
#' The input relation has `event_index`, `chrom`, `pos`, `id`, `ref`, `alt` and
#' `info`. The native builder returns one row per physical source record and
#' never merges a record with its mate. `status` is `reciprocal`, `unproven`,
#' `conflict`, `invalid` or `not_applicable`; `reason` is a stable identifier.
#' ALT syntax alone never proves a fusion, phase or inserted-only sequence:
#' `fusion_status` is always `not_asserted` and `phase_status`
#' `not_evaluated`.
#'
#' @param con A connection with DuckVEP loaded.
#' @param events_table Name of the caller's raw VCF record relation.
#' @param ... Named options for the native builder.
#' @return A data frame with one identity row per input record.
#' @export
rduckvep_prepare_breakend_pairs <- function(con, events_table, ...) {
  sql <- .duckvep_builder_sql(con, "duckvep_prepare_breakend_pairs_sql",
                              list(events_table), list(...))
  DBI::dbGetQuery(con, paste0("SELECT * FROM query(",
                             as.character(DBI::dbQuoteString(con, sql)), ")"))
}

#' Combine BND identity with endpoint genes
#'
#' `pairs_table` holds the result of [rduckvep_prepare_breakend_pairs()] and
#' `genes_table` has `event_index` and `gene_id`, one row per gene overlapped
#' by that physical endpoint. Each source record keeps its own row and gene
#' list. `candidate_partner_genes` records partner-gene evidence only: the
#' result never asserts a fusion, frame or phase.
#'
#' @inheritParams rduckvep_prepare_breakend_pairs
#' @param pairs_table Name of a relation with the identity builder columns.
#' @param genes_table Name of an `event_index`, `gene_id` relation.
#' @return A data frame with one row per physical source record.
#' @export
rduckvep_prepare_breakend_fusion <- function(con, pairs_table, genes_table, ...) {
  sql <- .duckvep_builder_sql(con, "duckvep_prepare_breakend_fusion_sql",
                              list(pairs_table, genes_table), list(...))
  DBI::dbGetQuery(con, paste0("SELECT * FROM query(",
                             as.character(DBI::dbQuoteString(con, sql)), ")"))
}
