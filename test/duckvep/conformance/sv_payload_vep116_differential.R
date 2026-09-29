#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(jsonlite) })
source("r/Rduckvep/R/builders.R")
source("r/Rduckvep/R/structural_geometry.R")
vcf <- "test/duckvep/conformance/data/sv_payload_grch38.vcf"
fasta <- Sys.getenv("DUCKVEP_GRCH38_FASTA", "/root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
cache <- Sys.getenv("DUCKVEP_GRCH38_VEP_CACHE", "/root/.cache/duckhts/vep/cache-grch38-chr21")
model <- Sys.getenv("DUCKVEP_GRCH38_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
stopifnot(identical(tail(system2("samtools", c("faidx", fasta,
  "21:13546123-13546123"), stdout = TRUE), 1L), "G"))
lines <- grep("^21\\t", readLines(vcf), value = TRUE)
records <- lapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1L]])
directory <- tempfile("sv-vep116-")
dir.create(directory)
binary <- file.path(directory, "duckvep.duckdb_extension")
stopifnot(file.copy(Sys.getenv("DUCKVEP_EXTENSION_FILE", "build/release/duckvep.duckdb_extension"), binary))
con <- dbConnect(duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
q <- function(value) as.character(dbQuoteString(con, value))
dbExecute(con, paste("LOAD", q(binary)))
prepare <- function(rows) {
  dbWriteTable(con, "sv_prepare_input", data.frame(
    event_index = seq_along(rows) - 1L,
    pos = vapply(rows, function(x) as.numeric(x[2L]), 0),
    ref = vapply(rows, function(x) x[4L], ""),
    alt = vapply(rows, function(x) x[5L], ""),
    info = vapply(rows, function(x) x[8L], "")), temporary = TRUE, overwrite = TRUE)
  rduckvep_prepare_sv_geometry(con, "sv_prepare_input")
}
geometry <- prepare(records)
stopifnot(all(geometry$status == "ok"),
  identical(geometry$nominal_start, rep(13546124, 3)),
  identical(geometry$nominal_end, rep(13546123, 3)),
  identical(geometry$outer_start[2L], 13546122),
  identical(geometry$inner_start[2L], 13546124),
  identical(geometry$inner_end[2L], 13546123),
  identical(geometry$outer_end[2L], 13546125),
  identical(geometry$inserted_sequence[3L], "ATG"),
  identical(geometry$source_sequence[2L], "ATG"),
  is.na(geometry$inserted_sequence[2L]))
local({
  on.exit({dbDisconnect(con, shutdown = TRUE); unlink(directory, recursive = TRUE)}, add = TRUE)
  oracle <- file.path(directory, "oracle.json")
  rc <- system2("scripts/run_species_vep116_docker.sh", c("homo_sapiens", "GRCh38",
    "116", cache, fasta, vcf, oracle))
  stopifnot(identical(rc, 0L))
  vep <- lapply(readLines(oracle), fromJSON, simplifyVector = FALSE)
  stopifnot(length(vep) == 3L)
  dbExecute(con, paste("ATTACH", q(model), "AS m (READ_ONLY)"))
  region <- dbGetQuery(con, "SELECT seq_region FROM m.bench_regions WHERE chrom='21'")$seq_region
  stopifnot(length(region) == 1L)
  queries <- c(
    "SELECT seq_region, sequence_length FROM m.bench_regions ORDER BY seq_region",
    paste("SELECT transcript_index, seq_region, transcript_start, transcript_end,",
      "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
      "codon_table, pre_cds_sequence, post_cds_sequence FROM m.bench_transcripts",
      "ORDER BY transcript_index"),
    paste("SELECT transcript_index, exon_start, exon_end, exon_cdna_start,",
      "exon_cdna_end, phase, end_phase FROM m.bench_exons",
      "ORDER BY transcript_index, exon_cdna_start"))
  load <- paste0("SELECT loaded FROM duckvep_model_load('sv_payload',",
    paste(vapply(queries, q, ""), collapse = ","), ")")
  stopifnot(isTRUE(dbGetQuery(con, load)$loaded))
  dbExecute(con, paste0("CREATE TEMP TABLE sv_events AS ",
    paste(vapply(seq_along(records), function(i) {
      row <- records[[i]]
      paste0("SELECT ", i - 1L, "::UBIGINT AS event_index, ", region,
        "::UINTEGER AS seq_region, ", row[2L], "::UBIGINT AS position, ",
        q(row[4L]), "::VARCHAR AS reference, ", q(row[5L]),
        "::VARCHAR AS alternate, ",
        if (i == 3L) "NULL" else paste0(row[2L], "::UBIGINT"),
        " AS end_position, ",
        if (i == 3L) "NULL" else "'INS'", "::VARCHAR AS structural_type, ",
        "NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, ",
        "NULL::UBIGINT AS mate_position")
    }, ""), collapse = " UNION ALL ")))
  dbWriteTable(con, "sv_provenance", data.frame(event_index = 0:2, geometry),
    temporary = TRUE)
  duck <- dbGetQuery(con, paste(
    "SELECT a.event_index, t.transcript_stable_id AS feature, a.consequence,",
    "a.duckvep_status, p.outer_start, p.outer_end, p.inserted_sequence,",
    "p.source_sequence FROM query(duckvep_annotate_sql('sv_events', 'sv_payload',",
    "struct_pack(rich := true))) a LEFT JOIN m.model_transcripts t",
    "USING (transcript_index) JOIN sv_provenance p USING (event_index)",
    "ORDER BY event_index, feature"))
  stopifnot(all(duck$duckvep_status == "supported"),
    all(duck$outer_start[duck$event_index == 1] == 13546122),
    all(duck$outer_end[duck$event_index == 1] == 13546125),
    all(duck$source_sequence[duck$event_index == 1] == "ATG"),
    all(duck$inserted_sequence[duck$event_index == 2] == "ATG"))
  pairs <- function(x) {
    terms <- unlist(lapply(x$transcript_consequences, function(t) paste(
      t$transcript_id, sort(unlist(t$consequence_terms)), sep = ":")), use.names = FALSE)
    sort(terms)
  }
  duck_pairs <- function(index) {
    rows <- duck[duck$event_index == index, , drop = FALSE]
    sort(unlist(lapply(seq_len(nrow(rows)), function(j) paste(
      rows$feature[j], sort(strsplit(rows$consequence[j], "&", fixed = TRUE)[[1L]]),
      sep = ":")), use.names = FALSE))
  }
  for (i in seq_along(vep)) {
    stopifnot(vep[[i]]$start == geometry$nominal_start[i],
      vep[[i]]$end == geometry$nominal_end[i])
  }
  confidence_vcf <- "test/duckvep/conformance/data/structural_confidence_grch38.vcf"
  confidence_json <- file.path(directory, "confidence.json")
  stopifnot(identical(system2("scripts/run_species_vep116_docker.sh",
    c("homo_sapiens", "GRCh38", "116", cache, fasta, confidence_vcf,
      confidence_json)), 0L))
  confidence_vep <- lapply(readLines(confidence_json), fromJSON,
    simplifyVector = FALSE)
  confidence_rows <- grep("^21\\t", readLines(confidence_vcf), value = TRUE)
  stopifnot(length(confidence_vep) == 12L, length(confidence_rows) == 12L)
  for (i in seq.int(1L, 12L, by = 2L)) {
    nominal <- strsplit(confidence_rows[[i]], "\t", fixed = TRUE)[[1L]]
    uncertain <- strsplit(confidence_rows[[i + 1L]], "\t", fixed = TRUE)[[1L]]
    pair <- prepare(list(nominal, uncertain))
    a <- pair[1L, ]
    b <- pair[2L, ]
    stopifnot(identical(a$status, "ok"), identical(b$status, "ok"),
      is.na(a$outer_start), !is.na(b$outer_start), !is.na(b$outer_end),
      a$nominal_start == b$nominal_start,
      a$nominal_end == b$nominal_end,
      a$nominal_start == confidence_vep[[i]]$start,
      a$nominal_end == confidence_vep[[i]]$end,
      b$nominal_start == confidence_vep[[i + 1L]]$start,
      b$nominal_end == confidence_vep[[i + 1L]]$end,
      identical(pairs(confidence_vep[[i]]), pairs(confidence_vep[[i + 1L]])))
  }
  stopifnot(identical(pairs(vep[[1L]]), pairs(vep[[2L]])),
    identical(duck_pairs(0), duck_pairs(1)))
  target <- "ENST00000427446"
  term <- function(x) unlist(x$consequence_terms)
  vep_term <- function(i) term(Filter(function(x) identical(x$transcript_id, target),
    vep[[i]]$transcript_consequences)[[1L]])
  duck_term <- function(i) duck$consequence[duck$event_index == i & duck$feature == target]
  stopifnot(identical(sort(strsplit(duck_term(0), "&", fixed = TRUE)[[1L]]),
      sort(vep_term(1))),
    identical(sort(strsplit(duck_term(2), "&", fixed = TRUE)[[1L]]),
      sort(vep_term(3))),
    identical(vep[[3L]]$allele_string, "-/ATG"))
  cat("VEP 116: 6 typed nominal/CI SV pairs; DuckVEP: INS confidence and literal payload predicates\n")
})
