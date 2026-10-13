#' Evaluate a supplied RNA product programme
#'
#' `reference` has `start1`, `dna`, and `strand`. `exons` and `blocks` are
#' data frames with `start1` and `end1`; blocks are in transcription order.
#' `sources` is a data frame with `source_id`, `start1`, `reference`, and
#' `alternate`. `programme` has `id`, `evidence_id`, `assumption_id`,
#' `expected_reference_rna`, `annotated_start1`, `alternative_start1`,
#' `alternative_trigger_id`, and `genetic_code`. Supply zero and an empty
#' trigger when there is no alternative initiation programme.
#'
#' Genomic positions and source IDs are unsigned 64-bit coordinates. Supply
#' decimal character values for integers outside R's exact double range. The
#' R result represents these unsigned 64-bit values as decimal character data.
#' The edit set is hypothetical, not evidence that calls were observed in cis.
#' The supplied RNA is validated against the reference and mapped blocks.
#' Peptides begin with methionine at the declared valid initiator and end at
#' the first reached stop; the returned RNA retains the complete mapped sequence.
#' `conditional_nmd` remains unknown without shared NMD eligibility rules.
#' This experimental versioned API supports disjoint, fixed-length substitutions;
#' it does not infer an RNA product or scan for alternative initiation.
#' @param conn A DBI connection with the DuckVEP extension loaded.
#' @param reference,exons,blocks,sources,programme Declared scenario inputs.
#' @return A one-row data frame with typed source and effect LIST columns.
#' @export
#' @examples
#' local({
#'   con <- rduckvep_connect()
#'   on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
#'   reference <- list(start1 = 1L, dna = "ATGAAATAA", strand = 1L)
#'   exons <- data.frame(start1 = 1L, end1 = 9L)
#'   sources <- data.frame(source_id = 0L, start1 = 4L,
#'                         reference = "A", alternate = "G")
#'   programme <- list(id = "annotated", evidence_id = "supplied",
#'     assumption_id = "joint_assignment", expected_reference_rna = reference$dna,
#'     annotated_start1 = 1L, alternative_start1 = 0L,
#'     alternative_trigger_id = "", genetic_code = 1L)
#'   answer <- rduckvep_product_scenario_v1(con, reference, exons, exons,
#'                                         sources, programme)
#'   answer[c("programme_id", "allele_evidence", "product_evidence", "conditional_nmd")]
#' })
rduckvep_product_scenario_v1 <- function(conn, reference, exons, blocks,
                                           sources, programme) {
  u64 <- function(x) {
    if (is.double(x) && !inherits(x, "integer64")) {
      if (any(!is.finite(x) | x != trunc(x) | x > 2^53 | x < 0)) {
        stop("Supply exact nonnegative genomic integers as decimal character values")
      }
      value <- format(x, scientific = FALSE, trim = TRUE, digits = 22L)
    } else {
      value <- as.character(x)
    }
    if (anyNA(value) || any(!grepl("^[0-9]+$", value))) {
      stop("Genomic integers must be unsigned decimal values")
    }
    value
  }
  intervals <- function(x) {
    stopifnot(is.data.frame(x), nrow(x) > 0L,
              all(c("start1", "end1") %in% names(x)))
    list(u64(x$start1), u64(x$end1))
  }
  stopifnot(is.data.frame(sources), nrow(sources) > 0L,
            all(c("source_id", "start1", "reference", "alternate") %in% names(sources)))
  e <- intervals(exons)
  b <- intervals(blocks)
  values <- list(u64(reference$start1), reference$dna, reference$strand,
                 e[[1L]], e[[2L]], b[[1L]], b[[2L]],
                 u64(sources$source_id), u64(sources$start1),
                 sources$reference, sources$alternate,
                 programme$id, programme$evidence_id, programme$assumption_id,
                 programme$expected_reference_rna, programme$annotated_start1,
                 programme$alternative_start1, programme$alternative_trigger_id,
                 programme$genetic_code)
  # A one-element outer list binds each entire vector as one SQL LIST value.
  values[4:11] <- lapply(values[4:11], list)
  sql <- paste0(
    "SELECT converted.* FROM (SELECT STRUCT_UPDATE(p, ",
    "sources := list_transform(p.sources, lambda s: STRUCT_PACK(",
    "source_id := s.source_id::VARCHAR, start1 := s.start1::VARCHAR, ",
    "reference := s.reference, alternate := s.alternate)), ",
    "source_effects := list_transform(p.source_effects, lambda s: STRUCT_PACK(",
    "source_id := s.source_id::VARCHAR, changed_dna := s.changed_dna, ",
    "changed_rna_bases := s.changed_rna_bases, ",
    "genomic_start1 := s.genomic_start1::VARCHAR, genomic_end1 := s.genomic_end1::VARCHAR, ",
    "rna_first1 := s.rna_first1, rna_last1 := s.rna_last1)), ",
    "reference_start1 := p.reference_start1::VARCHAR) AS converted ",
    "FROM (SELECT duckvep_product_scenario_v1(STRUCT_PACK(",
    "reference := STRUCT_PACK(start1 := ?::UBIGINT, dna := ?::VARCHAR, strand := ?::TINYINT), ",
    "exons := list_transform(list_zip(?::UBIGINT[], ?::UBIGINT[]), ",
    "lambda x: STRUCT_PACK(start1 := x[1], end1 := x[2])), ",
    "blocks := list_transform(list_zip(?::UBIGINT[], ?::UBIGINT[]), ",
    "lambda x: STRUCT_PACK(start1 := x[1], end1 := x[2])), ",
    "sources := list_transform(list_zip(?::UBIGINT[], ?::UBIGINT[], ",
    "?::VARCHAR[], ?::VARCHAR[]), lambda x: STRUCT_PACK(",
    "source_id := x[1], start1 := x[2], reference := x[3], alternate := x[4])), ",
    "programme := STRUCT_PACK(id := ?::VARCHAR, evidence_id := ?::VARCHAR, ",
    "assumption_id := ?::VARCHAR, expected_reference_rna := ?::VARCHAR, ",
    "annotated_start1 := ?::UINTEGER, alternative_start1 := ?::UINTEGER, ",
    "alternative_trigger_id := ?::VARCHAR, genetic_code := ?::UTINYINT)",
    ")) AS p) AS scenario) AS decoded")
  DBI::dbGetQuery(conn, sql, params = values)
}
