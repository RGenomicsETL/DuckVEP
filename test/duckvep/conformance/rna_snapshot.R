#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

run <- function() {
extension <- normalizePath(Sys.getenv("DUCKVEP_RNA_TEST_EXTENSION", "build/release/duckvep.duckdb_extension"), mustWork = TRUE)
fixture <- normalizePath("test/data/duckvep/rna/raw-cdna-v2.dvsnap", mustWork = TRUE)
con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE))
on.exit(dbDisconnect(con, shutdown = TRUE))
q <- function(value) as.character(dbQuoteString(con, value))
sql <- function(query) dbGetQuery(con, query)
dbExecute(con, paste("LOAD", q(extension)))
origin <- sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")
tmp <- tempfile("duckvep-rna-wire-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

stopifnot(sql(paste0("SELECT duckvep_model_restore('v2', ", q(fixture), ")"))[[1L]])
fp <- sql("SELECT CAST(_duckvep_model_fingerprint('v2') AS VARCHAR)")[[1L]]
stopifnot(identical(as.character(fp), "7835186884573121370"))
wire <- file.path(tmp, "v3.dvsnap")
stopifnot(sql(paste0("SELECT duckvep_model_save('v2', ", q(wire), ")"))[[1L]])
stopifnot(sql("SELECT duckvep_model_drop('v2')")[[1L]])
base <- sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")
bytes <- readBin(wire, "raw", n = file.info(wire)$size)
stopifnot(identical(bytes[9:16], c(as.raw(3L), raw(7L))))

for (mode in c("truncated", "corrupt")) {
  path <- file.path(tmp, paste0(mode, ".dvsnap"))
  payload <- bytes
  if (mode == "truncated") {
    payload <- payload[-length(payload)]
  } else {
    at <- 41L
    payload[at] <- as.raw(bitwXor(as.integer(payload[at]), 1L))
  }
  writeBin(payload, path)
  outcome <- tryCatch(sql(paste0("SELECT duckvep_model_restore('bad', ", q(path), ")"))[[1L]],
                      error = conditionMessage)
  stopifnot(is.character(outcome) || identical(outcome, FALSE),
            !sql("SELECT duckvep_model_drop('bad')")[[1L]])
  stopifnot(identical(base, sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")))
}
stopifnot(sql(paste0("SELECT duckvep_model_restore('v3', ", q(wire), ")"))[[1L]])
stopifnot(identical(as.character(sql("SELECT CAST(_duckvep_model_fingerprint('v3') AS VARCHAR)")[[1L]]), as.character(fp)))
stopifnot(sql("SELECT duckvep_model_drop('v3')")[[1L]])
legacy_fixture <- Sys.getenv("DUCKVEP_RNA_V1_FIXTURE")
if (nzchar(legacy_fixture)) {
  legacy_fixture <- normalizePath(legacy_fixture, mustWork = TRUE)
  stopifnot(sql(paste0("SELECT duckvep_model_restore('v1', ", q(legacy_fixture), ")"))[[1L]])
  legacy_fp <- as.character(sql("SELECT CAST(_duckvep_model_fingerprint('v1') AS VARCHAR)")[[1L]])
  legacy_v3 <- file.path(tmp, "legacy-v3.dvsnap")
  stopifnot(sql(paste0("SELECT duckvep_model_save('v1', ", q(legacy_v3), ")"))[[1L]])
  stopifnot(identical(readBin(legacy_fixture, "raw", n = 56L)[49:56],
                      readBin(legacy_v3, "raw", n = 56L)[49:56]))
  stopifnot(sql("SELECT duckvep_model_drop('v1')")[[1L]])
  stopifnot(sql(paste0("SELECT duckvep_model_restore('legacy-v3', ", q(legacy_v3), ")"))[[1L]])
  stopifnot(identical(legacy_fp, as.character(sql("SELECT CAST(_duckvep_model_fingerprint('legacy-v3') AS VARCHAR)")[[1L]])))
  stopifnot(sql("SELECT duckvep_model_drop('legacy-v3')")[[1L]])
}

dbExecute(con, "CREATE TABLE authority_regions AS SELECT 0::UINTEGER seq_region, 12::UBIGINT sequence_length")
dbExecute(con, "CREATE TABLE authority_transcripts AS SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region, 1::UBIGINT transcript_start, 12::UBIGINT transcript_end, -1::TINYINT strand, 0::UINTEGER gene_index, 0::UBIGINT transcript_flags, 4::UBIGINT cds_start, 12::UBIGINT cds_end, 'ATGCCCTAA'::BLOB cds_sequence, 2::UTINYINT codon_table, ''::BLOB pre_cds_sequence, 'AGG'::BLOB post_cds_sequence, 'ATGCCCTAAAGG'::BLOB cdna_sequence, 1::UBIGINT rna_cds_cdna_start, 9::UBIGINT rna_cds_cdna_end, 2::UTINYINT rna_codon_table")
dbExecute(con, "CREATE TABLE authority_exons AS SELECT 0::UINTEGER transcript_index, 1::UBIGINT exon_start, 12::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase")
load <- "SELECT loaded FROM duckvep_model_load('authority', 'SELECT * FROM authority_regions', 'SELECT * FROM authority_transcripts', 'SELECT * FROM authority_exons')"
stopifnot(sql(load)[[1L]])
authority_wire <- file.path(tmp, "authority.dvsnap")
stopifnot(sql(paste0("SELECT duckvep_model_save('authority', ", q(authority_wire), ")"))[[1L]])
authority_fp <- as.character(sql("SELECT CAST(_duckvep_model_fingerprint('authority') AS VARCHAR)")[[1L]])
stopifnot(sql("SELECT duckvep_model_drop('authority')")[[1L]])
stopifnot(sql(paste0("SELECT duckvep_model_restore('authority', ", q(authority_wire), ")"))[[1L]])
stopifnot(identical(authority_fp, as.character(sql("SELECT CAST(_duckvep_model_fingerprint('authority') AS VARCHAR)")[[1L]])))
stopifnot(sql("SELECT duckvep_model_drop('authority')")[[1L]])
base <- sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")

root <- normalizePath(".")
checker <- file.path(tmp, "rehash")
source <- file.path(root, "test/duckvep/conformance/rna_snapshot_rehash.c")
includes <- file.path(root, c("src", "src/include", "src/kernel/include", "third_party/cgranges", "third_party/htslib", "duckdb_capi"))
stopifnot(system2(Sys.getenv("CC", "cc"), c("-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections", paste0("-I", shQuote(includes)), shQuote(source), "-o", shQuote(checker))) == 0L)
number <- function(raw, at) sum(as.numeric(raw[at + 0:7]) * 256^(0:7))
for (mode in c("start", "end", "code", "presence", "length", "offset", "short")) {
  payload <- readBin(authority_wire, "raw", n = file.info(authority_wire)$size)
  offsets <- 609L + 8L * 51:53
  if (mode %in% c("start", "end", "code")) {
    pos <- number(payload, offsets[[match(mode, c("start", "end", "code"))]]) + 1L
    payload[[pos]] <- as.raw(if (mode == "code") 99L else 0L)
  } else if (mode == "presence") {
    payload[[159L]] <- as.raw(bitwAnd(as.integer(payload[[159L]]), bitwNot(8L)))
  } else if (mode == "length") {
    payload[[1057L + 8L * 52L]] <- as.raw(0L)
  } else if (mode == "offset") {
    payload[609L + 8L * 51L + 0:7] <- as.raw(rep(255L, 8L))
  } else {
    payload <- payload[-length(payload)]
  }
  path <- file.path(tmp, paste0("authority-", mode, ".dvsnap"))
  writeBin(payload, path)
  if (mode != "short" && mode != "offset") stopifnot(system2(checker, shQuote(path)) == 0L)
  failure <- tryCatch(sql(paste0("SELECT duckvep_model_restore('bad', ", q(path), ")"))[[1L]], error = conditionMessage)
  if (!is.character(failure) && !identical(failure, FALSE)) stop(mode, " was accepted")
  stopifnot(!sql("SELECT duckvep_model_drop('bad')")[[1L]],
            identical(base, sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")))
}
stopifnot(sql(paste0("SELECT duckvep_model_restore('authority', ", q(authority_wire), ")"))[[1L]],
          sql("SELECT duckvep_model_drop('authority')")[[1L]])
stopifnot(identical(origin, sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")))
cat("RNA authority snapshots, migration, rejection and ownership: OK\n")
}
run()
