#!/usr/bin/env Rscript

main <- function(repo) {
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION",
    file.path(repo, "build/release/duckvep.duckdb_extension")), mustWork = TRUE)
  directory <- tempfile("phase-arrangements-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  driver <- duckdb::duckdb(dbdir = file.path(directory, "calls.duckdb"),
    config = list(allow_unsigned_extensions = "true"), shared_home = FALSE)
  connection <- DBI::dbConnect(driver)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE, after = FALSE)
  quote <- function(x) as.character(DBI::dbQuoteString(connection, x))
  DBI::dbExecute(connection, paste("LOAD", quote(extension)))
  transcript <- paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "1::UBIGINT transcript_start,15::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,1::UBIGINT cds_start,",
    "15::UBIGINT cds_end,'ATGAAACCCGGGTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence")
  exons <- paste0("SELECT 0::UINTEGER transcript_index,1::UBIGINT exon_start,",
    "15::UBIGINT exon_end,1::UBIGINT exon_cdna_start,15::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase")
  loaded <- DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(",
    quote("phase_runtime"), ",", quote("SELECT 0::UINTEGER seq_region"), ",",
    quote(transcript), ",", quote(exons), ")"))
  stopifnot(loaded$loaded[[1L]] == 1L)
  run <- function(rows, ...) {
    DBI::dbExecute(connection, "DROP TABLE IF EXISTS phase_calls")
    DBI::dbExecute(connection, paste0("CREATE TABLE phase_calls AS SELECT * FROM (VALUES ",
      paste(rows, collapse = ","), ") v(event_index,seq_region,position,reference,alternate,",
      "alt_index,transcript_index,sample_index,alleles,phase_before,phase_set)"))
    relation <- paste0("duckvep_haplotype_arrangements(",
      quote("SELECT * FROM phase_calls"), ",", quote("phase_runtime"),
      if (length(list(...)) > 0L) paste0(",", paste(..., collapse = ",")) else "", ")")
    DBI::dbGetQuery(connection, paste0("SELECT * FROM ", relation))
  }
  # Independent exhaustive states for these two substitutions: reference MKPG*,
  # singles M* and MI PG*, and cis MLPG*.
  unphased <- c(
    "(1,0,4,'A','T',1,0,7,[0,1],NULL::BOOLEAN[],NULL::BIGINT)",
    "(2,0,5,'A','T',1,0,7,[0,1],NULL::BOOLEAN[],NULL::BIGINT)")
  two <- run(unphased)
  stopifnot(nrow(two) == 8L, identical(sort(unique(two$protein)), sort(c("M*", "MIPG*", "MLPG*", "MKPG*"))),
    all(two$original_allele0 == 0L), all(two$original_allele1 == 1L),
    all(is.na(two$original_phase_set)), any(two$prediction_semantics == "reference_replay"),
    all(two$prediction_status %in% c("eligible_classifier_pending", "predicted")),
    all(two$prediction_semantics[two$contributes] == "hypothetical_assignment"))
  three <- run(c(unphased, "(3,0,8,'C','T',1,0,7,[0,1],NULL::BOOLEAN[],NULL::BIGINT)"))
  stopifnot(nrow(three) == 24L, length(unique(three$hypothesis_id)) == 4L)
  same_ps <- run(sub("NULL::BOOLEAN\\[\\],NULL::BIGINT", "[false,true],11::BIGINT", unphased))
  cross_ps <- run(c(
    "(1,0,4,'A','T',1,0,7,[0,1],[false,true],11::BIGINT)",
    "(2,0,5,'A','T',1,0,7,[0,1],[false,true],12::BIGINT)"))
  stopifnot(length(unique(same_ps$hypothesis_id)) == 1L,
    length(unique(cross_ps$hypothesis_id)) == 2L,
    identical(sort(unique(cross_ps$original_phase_set)), c(11, 12)))
  compensating <- run(c(
    "(4,0,4,'AA','A',1,0,7,[0,1],NULL::BOOLEAN[],NULL::BIGINT)",
    "(5,0,7,'C','CCC',1,0,7,[0,1],NULL::BOOLEAN[],NULL::BIGINT)"))
  stopifnot(nrow(compensating) == 8L, any(compensating$protein != "MKPG*"),
    any(compensating$prediction_semantics == "reference_replay"))
  exhausted <- try(run(unphased, "max_arrangements:=1", "max_replays:=4"), silent = TRUE)
  stopifnot(inherits(exhausted, "try-error"))
  repeated <- run(unphased)
  stopifnot(identical(two[, names(two)], repeated[, names(two)]))
  DBI::dbGetQuery(connection, "SELECT duckvep_model_drop('phase_runtime')")
  cat("Native phase arrangements: unphased cis/trans, three-site, phase-block, compensating-indel, reference, capacity and reuse cases verified\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: phase_arrangements_runtime.R [repository]")
main(normalizePath(if (length(args) > 0L) args[[1L]] else ".", mustWork = TRUE))
