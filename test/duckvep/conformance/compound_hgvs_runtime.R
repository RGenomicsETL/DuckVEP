#!/usr/bin/env Rscript

main <- function(repo) {
  fixtures <- file.path(repo, "test/duckvep/conformance/compound_hgvs_oracle_fixtures")
  extension <- normalizePath(Sys.getenv("DUCKVEP_EXTENSION",
    file.path(repo, "build/release/duckvep.duckdb_extension")), mustWork = TRUE)
  cases <- lapply(readLines(file.path(fixtures, "cases.jsonl")), jsonlite::fromJSON)
  cases <- Filter(function(case) !is.null(case$expected_native_hgvsc), cases)
  stopifnot(length(cases) >= 4L)
  fasta <- readLines(file.path(fixtures, "reference.fa"))
  headers <- which(startsWith(fasta, ">"))
  sequences <- vapply(seq_along(headers), function(i) {
    end <- if (i == length(headers)) length(fasta) else headers[[i + 1L]] - 1L
    paste(fasta[seq.int(headers[[i]] + 1L, end)], collapse = "")
  }, character(1L))
  names(sequences) <- substring(fasta[headers], 2L)
  gff <- read.delim(file.path(fixtures, "model.gff3"), header = FALSE, comment.char = "#",
    stringsAsFactors = FALSE, col.names = c("reference", "source", "type", "start", "end",
      "score", "strand", "phase", "attributes"))
  complement <- c(A = "T", C = "G", G = "C", T = "A")
  reverse_complement <- function(sequence) {
    paste(rev(unname(complement[strsplit(sequence, "", fixed = TRUE)[[1L]]])), collapse = "")
  }
  directory <- tempfile("compound-hgvs-runtime-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  driver <- duckdb::duckdb(dbdir = file.path(directory, "calls.duckdb"),
    config = list(allow_unsigned_extensions = "true"), shared_home = FALSE)
  connection <- DBI::dbConnect(driver)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE, after = FALSE)
  sql_quote <- function(value) as.character(DBI::dbQuoteString(connection, value))
  DBI::dbExecute(connection, paste("LOAD", sql_quote(extension)))
  for (case in cases) {
    strand <- if (case$strand == "+") 1L else -1L
    genomic <- sequences[[case$reference_id]]
    tx <- gff[gff$reference == case$reference_id & gff$type == "mRNA", ]
    cds <- gff[gff$reference == case$reference_id & gff$type == "CDS", ]
    exon <- gff[gff$reference == case$reference_id & gff$type == "exon", ]
    stopifnot(nrow(tx) == 1L, nrow(cds) == 1L, nrow(exon) == 1L,
      tx$strand == case$strand, cds$phase == "0")
    oriented <- function(begin, end) {
      if (begin > end) return("")
      sequence <- substr(genomic, begin, end)
      if (strand == -1L) sequence <- paste(rev(complement[
        strsplit(sequence, "", fixed = TRUE)[[1L]]]), collapse = "")
      sequence
    }
    reference <- oriented(cds$start, cds$end)
    left <- oriented(tx$start, cds$start - 1L)
    right <- oriented(cds$end + 1L, tx$end)
    pre <- if (strand == 1L) left else right
    post <- if (strand == 1L) right else left
    model <- paste0("runtime_", case$case_id)
    transcript <- sprintf(paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
      "%d::UBIGINT transcript_start,%d::UBIGINT transcript_end,%d::TINYINT strand,",
      "0::UINTEGER gene_index,3::UBIGINT transcript_flags,%d::UBIGINT cds_start,",
      "%d::UBIGINT cds_end,%s::BLOB cds_sequence,1::UTINYINT codon_table,",
      "%s::BLOB pre_cds_sequence,%s::BLOB post_cds_sequence"), tx$start, tx$end, strand,
      cds$start, cds$end, sql_quote(reference), sql_quote(pre), sql_quote(post))
    exons <- sprintf(paste0("SELECT 0::UINTEGER transcript_index,%d::UBIGINT exon_start,",
      "%d::UBIGINT exon_end,1::UBIGINT exon_cdna_start,%d::UBIGINT exon_cdna_end,",
      "0::TINYINT phase,0::TINYINT end_phase"), exon$start, exon$end, exon$end - exon$start + 1L)
    loaded <- DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(",
      sql_quote(model), ",", sql_quote("SELECT 0::UINTEGER seq_region"), ",",
      sql_quote(transcript), ",", sql_quote(exons), ")"))
    stopifnot(loaded$loaded[[1L]] == 1L)
    edits <- case$edits
    input_positions <- edits$position
    input_references <- edits$reference
    input_alternates <- edits$alternate
    for (index in seq_len(nrow(edits))) {
      if (!nzchar(input_alternates[[index]])) {
        stopifnot(input_positions[[index]] > 1L)
        anchor <- substr(reference, input_positions[[index]] - 1L, input_positions[[index]] - 1L)
        input_positions[[index]] <- input_positions[[index]] - 1L
        input_references[[index]] <- paste0(anchor, input_references[[index]])
        input_alternates[[index]] <- anchor
      }
    }
    edit_widths <- nchar(input_references)
    positions <- if (strand == 1L) cds$start - 1L + input_positions else
      cds$end - input_positions - edit_widths + 2L
    refs <- if (strand == 1L) input_references else vapply(
      input_references, reverse_complement, character(1L))
    alts <- if (strand == 1L) input_alternates else vapply(
      input_alternates, reverse_complement, character(1L))
    rows <- sprintf("(%d,0,%d,%s,%s,1,0,0,[1,1],[false,true],NULL::BIGINT)",
      seq_len(nrow(edits)), positions, sql_quote(refs), sql_quote(alts))
    DBI::dbExecute(connection, "DROP TABLE IF EXISTS runtime_calls")
    DBI::dbExecute(connection, paste0("CREATE TABLE runtime_calls AS SELECT * FROM (VALUES ",
      paste(rows, collapse = ","), ") v(event_index,seq_region,position,reference,alternate,",
      "alt_index,transcript_index,sample_index,alleles,phase_before,phase_set)"))
    relation <- paste0("duckvep_haplotypes(", sql_quote("SELECT * FROM runtime_calls"), ",",
      sql_quote(model), ",hgvs:=true)")
    schema <- DBI::dbGetQuery(connection, paste0("DESCRIBE SELECT * FROM ", relation))
    stopifnot(nrow(schema) == 36L, tail(schema$column_name, 1L) == "nominal_length_diff",
    all(c("prediction_reference_protein", "prediction_protein") %in% schema$column_name))
    result <- DBI::dbGetQuery(connection, paste0("SELECT cds,protein,hgvsc,hgvsc_status,hgvsp,",
      "list_sort(list_transform(contributors,lambda c: c.event_index))::VARCHAR source_ids ",
      "FROM ", relation))
    replay_cds <- function(sequence) {
      for (index in order(edits$position, decreasing = TRUE)) {
        position <- edits$position[[index]]
        ref <- edits$reference[[index]]
        alt <- edits$alternate[[index]]
        observed <- substr(sequence, position, position + nchar(ref) - 1L)
        if (!identical(observed, ref)) {
          stop(case$case_id, ": expected ", ref, " at c.", position, ", got ", observed,
            call. = FALSE)
        }
        sequence <- paste0(substr(sequence, 1L, position - 1L), alt,
          substr(sequence, position + nchar(ref), nchar(sequence)))
      }
      sequence
    }
    expected_cds <- replay_cds(reference)
    if (!is.null(case$expected_final_cds)) {
      stopifnot(identical(expected_cds, case$expected_final_cds))
    }
    expected_ids <- paste0("[", paste(seq_len(nrow(edits)), collapse = ", "), "]")
    check_names <- c("one result", "hgvsc status", "hgvsc", "CDS", "contributors")
    checks <- c(nrow(result) == 1L, result$hgvsc_status[[1L]] == "ok",
      identical(result$hgvsc[[1L]], case$expected_native_hgvsc),
      identical(result$cds[[1L]], expected_cds),
      identical(result$source_ids[[1L]], expected_ids))
    if (!is.null(case$expected_native_hgvsp)) {
      check_names <- c(check_names, "hgvsp")
      checks <- c(checks, identical(result$hgvsp[[1L]], case$expected_native_hgvsp))
    }
    if (!is.null(case$expected_protein_predicted)) {
      check_names <- c(check_names, "protein")
      checks <- c(checks, identical(result$protein[[1L]], case$expected_protein_predicted))
    }
    if (!all(checks)) {
      failed <- check_names[!checks]
      detail <- if ("hgvsc" %in% failed) paste0("; expected ", case$expected_native_hgvsc,
        ", got ", result$hgvsc[[1L]]) else ""
      stop(case$case_id, ": ", paste(failed, collapse = ", "), detail, call. = FALSE)
    }
    if (case$case_id == "plus_two_missense") {
      ambiguous <- DBI::dbGetQuery(connection, paste0("SELECT hgvsc,hgvsc_status FROM ",
        "duckvep_haplotypes(", sql_quote(paste0("SELECT event_index,seq_region,position,",
          "reference,[alternate] alternates,transcript_index,sample_index,'0/1' gt ",
          "FROM runtime_calls")), ",", sql_quote(model),
        ",input_mode:='source_records',phase_policy:='vep_compat',hgvs:=true) WHERE edit_count >= 2"))
      stopifnot(nrow(ambiguous) > 0L, all(is.na(ambiguous$hgvsc)),
        all(ambiguous$hgvsc_status == "unproven_cis"))
    }
    DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", sql_quote(model), ")"))
  }
  cat("Native compound HGVS:", length(cases),
    "standards-directed DNA/protein cases; both strands, CDS replay, contributors and ambiguity verified\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: compound_hgvs_runtime.R [repository]")
main(normalizePath(if (length(args) == 1L) args[[1L]] else ".", mustWork = TRUE))
