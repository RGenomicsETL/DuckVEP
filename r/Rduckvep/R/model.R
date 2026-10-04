#' Load a resident transcript model, on either DuckDB C API host
#'
#' The v1 host runs the three (plus optional) relation queries on a private
#' connection with `duckvep_model_load()`. The v2 host cannot: its model sink is
#' a COPY format, so the same load is a sequence of `COPY` statements into
#' staging followed by `duckvep_model_publish()`, and those statements come from
#' `duckvep_model_load_sql()`. This function picks the host by what the loaded
#' extension provides and hides the difference; the result is the same resident
#' model. One visible difference remains: on v2 the queries run on `con` in its
#' own transaction, so temporary tables and uncommitted rows are visible to
#' them; on v1 the queries run on a separate connection and see only committed,
#' permanent relations.
#'
#' @param con A DBI connection with DuckVEP loaded.
#' @param name Model name (non-empty, not yet loaded).
#' @param regions_query,transcripts_query,exons_query SELECT statements that
#'   return the typed, ordered relations (see the function reference).
#' @param mature_mirna_query,peptide_edit_query,interval_feature_query Optional
#'   SELECT statements for the optional relations.
#' @param reference_fasta Optional path of an indexed reference FASTA.
#' @param transcript_coverage_complete Optional logical (`NA` is an error).
#' @return `TRUE`, invisibly, when the model is loaded.
#' @export
rduckvep_load_model <- function(con, name, regions_query, transcripts_query, exons_query,
                                mature_mirna_query = NULL, peptide_edit_query = NULL,
                                interval_feature_query = NULL, reference_fasta = NULL,
                                transcript_coverage_complete = NULL) {
  required <- list(name, regions_query, transcripts_query, exons_query)
  if (!all(vapply(required, function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x), NA))) {
    stop("name and the three queries must be non-empty strings", call. = FALSE)
  }
  optional <- list(mature_mirna_query = mature_mirna_query, peptide_edit_query = peptide_edit_query,
                   interval_feature_query = interval_feature_query, reference_fasta = reference_fasta,
                   transcript_coverage_complete = transcript_coverage_complete)
  optional <- optional[!vapply(optional, is.null, NA)]
  v2 <- nrow(DBI::dbGetQuery(con, "SELECT 1 FROM duckdb_functions()
    WHERE function_name = 'duckvep_model_publish' LIMIT 1")) > 0L
  if (v2) {
    args <- vapply(required, function(value) .duckvep_builder_literal(con, value), "")
    if (length(optional)) {
      fields <- paste0(as.character(DBI::dbQuoteIdentifier(con, names(optional))), " := ",
                       vapply(optional, function(value) .duckvep_builder_literal(con, value), ""))
      args <- c(args, paste0("struct_pack(", paste(fields, collapse = ", "), ")"))
    }
    statements <- DBI::dbGetQuery(con, paste0("SELECT unnest(duckvep_model_load_sql(",
                                              paste(args, collapse = ", "), ")) AS statement"))$statement
    for (statement in statements) DBI::dbExecute(con, statement)
  } else {
    args <- vapply(required, function(value) .duckvep_builder_literal(con, value), "")
    named <- vapply(names(optional), function(key) {
      paste0(key, " := ", .duckvep_builder_literal(con, optional[[key]]))
    }, "")
    DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load(",
                                paste(c(args, unname(named)), collapse = ", "), ")"))
  }
  invisible(TRUE)
}

.duckvep_model_path_call <- function(con, func, name, path) {
  if (!all(vapply(list(name, path), function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x), NA))) {
    stop("name and path must be non-empty strings", call. = FALSE)
  }
  DBI::dbGetQuery(con, paste0("SELECT ", func, "(", .duckvep_builder_literal(con, name), ", ",
                              .duckvep_builder_literal(con, path), ") AS done"))
  invisible(TRUE)
}

#' Save and restore a resident model as a snapshot file
#'
#' `rduckvep_save_model()` writes a loaded model's native arrays to one file with
#' `duckvep_model_save()`. `rduckvep_restore_model()` loads a model from that
#' file with `duckvep_model_restore()`: the file is mapped read-only instead of
#' copied, so restoring is several times faster than loading the relations, and
#' every R session that restores the same file shares one copy of it in memory.
#' A snapshot is validated when it is restored and is tied to the DuckVEP build
#' that wrote it; a build with a different model layout refuses it. It records
#' the path of the model's reference FASTA, not its bytes.
#'
#' @param con A DBI connection with DuckVEP loaded.
#' @param name Model name: a loaded model to save, or a new name to restore under.
#' @param path Path of the snapshot file.
#' @return `TRUE`, invisibly.
#' @export
rduckvep_save_model <- function(con, name, path) {
  .duckvep_model_path_call(con, "duckvep_model_save", name, path)
}

#' @rdname rduckvep_save_model
#' @export
rduckvep_restore_model <- function(con, name, path) {
  .duckvep_model_path_call(con, "duckvep_model_restore", name, path)
}
