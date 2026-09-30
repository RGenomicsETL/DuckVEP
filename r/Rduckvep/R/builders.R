#' Generate DuckVEP SQL on the caller's connection
#'
#' These functions invoke the extension's native SQL builders. The returned SQL
#' can be passed to `DBI::dbGetQuery()` or used in a `CREATE TABLE AS` statement.
#' Required relation names are passed as strings; the native builders quote them.
#' Optional named arguments are forwarded as a DuckDB STRUCT. Load the extension
#' on `con` with [rduckvep_load()] before calling a builder.
#'
#' @param con A DBI connection with DuckVEP loaded.
#' @param events_table Event relation name.
#' @param model_name Loaded model name.
#' @param annotations_table Annotation relation name.
#' @param transcripts_table Transcript relation name.
#' @param reference_table Reference sequence chunk relation (`chrom`, `start`, `end`, `seq`).
#' @param core_schema Ensembl core schema name.
#' @param reference_chunks_table Reference chunk relation name.
#' @param assembly Assembly name.
#' @param funcgen_schema Ensembl funcgen schema name.
#' @param regions_table Canonical region relation name.
#' @param source_name,source_version Source identity strings.
#' @param source_manifest_sha256,reference_sha256 Source and reference digests.
#' @param transcript_filter Transcript selection description.
#' @param ... Named builder options (for annotation: `hgvs`, `rich`,
#'   `upstream_distance`, `downstream_distance`; for Ensembl regions and
#'   transcripts: `species_id`; for receipts: `regulation_features_table`; for
#'   the loss-of-function relation: `gerp`, `ancestor`, `phylocsf`,
#'   `min_intron_size`, `gerp_end_trunc_cutoff`, `check_complete_cds`).
#'   Other builders accept no options. Each value must be a length-one string,
#'   logical, integer, double or `NA`.
#' @return A length-one SQL string.
#' @name rduckvep_builders
NULL

.duckvep_builder_literal <- function(con, value) {
  if (length(value) != 1L || !is.atomic(value) ||
      !is.null(dim(value)) ||
      !(is.character(value) || is.logical(value) || is.numeric(value))) {
    stop("Builder arguments must be scalar strings, logicals or numbers", call. = FALSE)
  }
  if (is.na(value)) return("NULL")
  if (is.character(value)) return(as.character(DBI::dbQuoteString(con, value)))
  if (is.logical(value)) return(if (value) "true" else "false")
  if (!is.finite(value)) stop("Builder numbers must be finite", call. = FALSE)
  format(value, scientific = FALSE, trim = TRUE)
}

.duckvep_builder_sql <- function(con, builder, required, options) {
  args <- vapply(required, function(value) .duckvep_builder_literal(con, value), "")
  if (length(options)) {
    keys <- names(options)
    if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
      stop("Builder options must have unique nonempty names", call. = FALSE)
    }
    fields <- paste0(as.character(DBI::dbQuoteIdentifier(con, keys)), " := ",
                     vapply(options, function(value) .duckvep_builder_literal(con, value), ""))
    args <- c(args, paste0("struct_pack(", paste(fields, collapse = ", "), ")"))
  }
  call <- paste0("SELECT ", builder, "(", paste(args, collapse = ", "), ") AS sql")
  DBI::dbGetQuery(con, call)$sql[[1L]]
}

#' @rdname rduckvep_builders
#' @export
rduckvep_annotate_sql <- function(con, events_table, model_name, ...) {
  .duckvep_builder_sql(con, "duckvep_annotate_sql", list(events_table, model_name), list(...))
}

#' @rdname rduckvep_builders
#' @export
rduckvep_transcript_projection_sql <- function(con, events_table, annotations_table,
                                                transcripts_table, ...) {
  .duckvep_builder_sql(con, "duckvep_transcript_projection_sql",
                       list(events_table, annotations_table, transcripts_table), list(...))
}

#' @rdname rduckvep_builders
#' @export
rduckvep_lof_sql <- function(con, annotations_table, transcripts_table, reference_table, ...) {
  .duckvep_builder_sql(con, "duckvep_lof_sql",
                       list(annotations_table, transcripts_table, reference_table), list(...))
}

#' @rdname rduckvep_builders
#' @export
rduckvep_ensembl_regions_sql <- function(con, core_schema, reference_chunks_table,
                                        assembly, ...) {
  .duckvep_builder_sql(con, "duckvep_ensembl_regions_sql",
                       list(core_schema, reference_chunks_table, assembly), list(...))
}

#' @rdname rduckvep_builders
#' @export
rduckvep_ensembl_transcripts_sql <- function(con, core_schema, reference_chunks_table,
                                           assembly, ...) {
  .duckvep_builder_sql(con, "duckvep_ensembl_transcripts_sql",
                       list(core_schema, reference_chunks_table, assembly), list(...))
}

#' @rdname rduckvep_builders
#' @export
rduckvep_ensembl_regulation_features_sql <- function(con, funcgen_schema, regions_table, ...) {
  .duckvep_builder_sql(con, "duckvep_ensembl_regulation_features_sql",
                       list(funcgen_schema, regions_table), list(...))
}

#' @rdname rduckvep_builders
#' @export
rduckvep_model_receipt_sql <- function(con, regions_table, transcripts_table,
                                      source_name, source_version, assembly,
                                      source_manifest_sha256, reference_sha256,
                                      transcript_filter, ...) {
  .duckvep_builder_sql(con, "duckvep_model_receipt_sql",
                       list(regions_table, transcripts_table, source_name, source_version,
                            assembly, source_manifest_sha256, reference_sha256,
                            transcript_filter), list(...))
}

#' Annotate events using a loaded DuckVEP model
#'
#' @inheritParams rduckvep_builders
#' @return A data frame of annotated events.
#' @export
rduckvep_annotate <- function(con, events_table, model_name, ...) {
  sql <- rduckvep_annotate_sql(con, events_table, model_name, ...)
  DBI::dbGetQuery(con, sql)
}
