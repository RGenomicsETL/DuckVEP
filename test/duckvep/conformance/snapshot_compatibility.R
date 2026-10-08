#!/usr/bin/env Rscript
# Published v1 and combined v2 snapshot wire and ownership conformance.
suppressPackageStartupMessages(library(DBI))
args <- commandArgs(trailingOnly = TRUE)
legacy_mode <- length(args) == 2L && args[[1L]] == "write-v1"

run <- function(directory, extension, legacy) {
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(dbQuoteString(con, x))
  sql <- function(x) dbGetQuery(con, x)
  exec <- function(x) invisible(dbExecute(con, x))
  exec(paste("LOAD", q(normalizePath(extension, mustWork = TRUE))))
  budget <- function() sql("SELECT owner, current_bytes FROM duckvep_native_budget() ORDER BY owner")
  fasta <- file.path(directory, "reference.fa")
  if (legacy) {
    sequence <- paste0("ATGTGAGAATAA", "AAAAAAA", "CCCAAATTT", "AA")
    writeLines(c(">chr1", sequence), fasta, useBytes = TRUE)
    writeLines(paste("chr1", 30, 6, 30, 31, sep = "\t"), paste0(fasta, ".fai"))
  }
  exec("CREATE TABLE regions AS SELECT 1::UINTEGER seq_region, 30::UBIGINT sequence_length, 'chr1'::VARCHAR seq_region_name")
  exec("CREATE TABLE transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, 1::UBIGINT transcript_start, 12::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 0::UBIGINT transcript_flags, 1::UBIGINT cds_start, 12::UBIGINT cds_end, 'ATGTGAGAATAA'::BLOB cds_sequence, 1::UTINYINT codon_table")
  exec("CREATE TABLE exons AS SELECT 0::UINTEGER transcript_index, 1::UBIGINT exon_start, 12::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase")
  exec("CREATE TABLE edits AS SELECT 0::UINTEGER transcript_index, 2::UINTEGER protein_position, 'U'::VARCHAR alternate_amino_acid")
  exec("CREATE TABLE events AS SELECT i::UBIGINT event_index, 1::UINTEGER seq_region, p::UBIGINT AS position, r::VARCHAR reference, a::VARCHAR alternate, NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change, NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position FROM (VALUES (0,4,'T','C'),(1,7,'G','A'),(2,21,'C','T')) v(i,p,r,a)")
  load <- function(name, typed = FALSE, full = FALSE, curated = FALSE) {
    tx <- "SELECT * FROM transcripts"
    ex <- "SELECT * FROM exons"
    if (full) {
      tx <- paste0("SELECT *, ''::BLOB pre_cds_sequence, ''::BLOB post_cds_sequence, cds_sequence cdna_sequence FROM transcripts UNION ALL SELECT 1::UINTEGER, 1::UINTEGER, 20::UBIGINT, 28::UBIGINT, 1::TINYINT, 1::UINTEGER, 0::UBIGINT, NULL::UBIGINT, NULL::UBIGINT, NULL::BLOB, NULL::UTINYINT, NULL::BLOB, NULL::BLOB, 'CCCAAATTT'::BLOB")
      ex <- paste0(ex, " UNION ALL SELECT 1::UINTEGER,20::UBIGINT,28::UBIGINT,1::UBIGINT,9::UBIGINT,0::TINYINT,0::TINYINT")
    }
    edit <- if (typed) "SELECT *, '_selenocysteine'::VARCHAR edit_code FROM edits" else "SELECT * FROM edits"
    sql(paste0("SELECT * FROM duckvep_model_load(", q(name), ", ", q("SELECT * FROM regions"), ", ", q(tx), ", ", q(ex), ", reference_fasta := ", q(fasta), if (curated) paste0(", peptide_edit_query := ", q(edit)), ")"))
  }
  annotate <- function(model) {
    # Columns shared by the published writer and current reader, without display aliases.
    sql(paste0("SELECT event_index, transcript_index, consequence_mask, region_mask, impact_code, status_code, reason_code, cdna_position, cds_position, protein_position, reference_amino_acid_code, alternate_amino_acid_code, transcript_hgvs, protein_hgvs, transcript_hgvs_status, transcript_hgvs_reason, protein_hgvs_status, protein_hgvs_reason FROM query(duckvep_annotate_sql('events', ", q(model), ", struct_pack(hgvs := true))) ORDER BY event_index, transcript_index NULLS FIRST"))
  }
  save <- function(model, path) stopifnot(sql(paste0("SELECT duckvep_model_save(", q(model), ",", q(path), ") ok"))$ok)
  restore <- function(model, path) stopifnot(sql(paste0("SELECT duckvep_model_restore(", q(model), ",", q(path), ") ok"))$ok)
  drop <- function(model) exec(paste0("SELECT duckvep_model_drop(", q(model), ")"))
  bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)
  fingerprint <- function(path) bytes(path)[49:56]
  if (legacy) {
    for (name in c("plain", "curated")) {
      print(load(name, curated = name == "curated"))
      path <- file.path(directory, paste0(name, ".v1"))
      save(name, path)
      stopifnot(identical(bytes(path)[9:16], c(as.raw(1), raw(7))))
      saveRDS(annotate(name), file.path(directory, paste0(name, ".rds")))
      cat(name, "v1 fingerprint bytes:", paste(format(fingerprint(path)), collapse = ""), "\n")
    }
    return(invisible(NULL))
  }
  # Re-checksum structurally edited fixtures with the implementation's wire checksum.
  # This lets required-array and slice checks run beyond checksum rejection.
  root <- normalizePath(".")
  native <- file.path(directory, "checksum.c")
  binary <- file.path(directory, "checksum")
  writeLines(c(
    paste0('#include "', root, '/src/core/duckvep_core_snapshot.c"'),
    'int main(int argc, char **argv) {',
    '  FILE *file; unsigned char *data; long size; uint64_t hash;',
    '  snap_v1_header_t *v1; const uint64_t *offset, *bytes; size_t sections, i;',
    '  if (argc != 2 || (file = fopen(argv[1], "rb+")) == NULL) return 1;',
    '  if (fseek(file, 0, SEEK_END) || (size = ftell(file)) < 0 || fseek(file, 0, SEEK_SET)) return 2;',
    '  data = malloc((size_t)size); if (!data || fread(data, 1, (size_t)size, file) != (size_t)size) return 3;',
    '  v1 = (snap_v1_header_t *)data;',
    '  if (v1->version == 1) { offset = v1->offset; bytes = v1->bytes; sections = SNAP_V1_SECTIONS; }',
    '  else { snap_header_t *v2 = (snap_header_t *)data; offset = v2->offset; bytes = v2->bytes; sections = SNAP_SECTIONS; }',
    '  v1->checksum = 0; hash = snap_hash(v1->version, data, (size_t)v1->header_bytes);',
    '  for (i = 0; i < sections; i++) hash = snap_hash(hash, data + offset[i], (size_t)bytes[i]);',
    '  v1->checksum = hash;',
    '  if (fseek(file, 0, SEEK_SET) || fwrite(data, 1, (size_t)size, file) != (size_t)size) return 4;',
    '  free(data); return fclose(file) != 0;',
    '}'
  ), native)
  includes <- file.path(root, c("src", "src/include", "src/kernel/include", "third_party/cgranges", "third_party/htslib", "duckdb_capi"))
  stopifnot(system2(Sys.getenv("CC", "cc"), c("-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections", paste0("-I", shQuote(includes)), shQuote(native), "-o", shQuote(binary))) == 0L)
  u64 <- function(x) as.raw(floor(x / 256^(0:7)) %% 256)
  number <- function(raw, at) sum(as.numeric(raw[at + 0:7]) * 256^(0:7))
  for (name in c("plain", "curated")) {
    path <- file.path(directory, paste0(name, ".v1"))
    for (i in seq_len(4L)) {
      restore(name, path)
      stopifnot(identical(annotate(name), readRDS(file.path(directory, paste0(name, ".rds")))))
      copy <- file.path(directory, paste0(name, ".v2"))
      save(name, copy)
      stopifnot(identical(fingerprint(path), fingerprint(copy)))
      wire <- bytes(copy)
      stopifnot(all(wire[105:112] == as.raw(0)), bitwAnd(as.integer(wire[[129L]]), 4L) == 0L,
                bitwAnd(as.integer(wire[[138L]]), 128L) == 0L,
                bitwAnd(as.integer(wire[[139L]]), 1L) == 0L,
                bitwAnd(as.integer(wire[[141L]]), 20L) == 0L)
      restore("legacy_v2", copy)
      stopifnot(identical(annotate(name), annotate("legacy_v2")))
      drop("legacy_v2")
      drop(name)
    }
    restore(name, path)
  }
  print(load("combined", typed = TRUE, full = TRUE, curated = TRUE))
  combined <- file.path(directory, "combined.v2")
  save("combined", combined)
  expected <- annotate("combined")
  restore("combined_copy", combined)
  stopifnot(identical(expected, annotate("combined_copy")))
  save("combined_copy", file.path(directory, "combined_copy.v2"))
  stopifnot(identical(bytes(combined), bytes(file.path(directory, "combined_copy.v2"))))
  wire <- bytes(combined)
  code_offset <- number(wire, 497L + 8L * 34L)
  stopifnot(bitwAnd(as.integer(wire[[129L]]), 4L) == 4L,
            bitwAnd(as.integer(wire[[141L]]), 20L) == 20L,
            wire[[code_offset + 1L]] == as.raw(2))
  cat("combined fingerprint bytes:", paste(format(fingerprint(combined)), collapse = ""), "\n")
  stopifnot(expected$transcript_hgvs[expected$event_index == 2 & expected$transcript_index == 1] == "n.2C>T",
            expected$reference_amino_acid_code[expected$event_index == 0 & expected$transcript_index == 0] == 85L)
  if ("conditional" %in% args) {
    exec("CREATE TABLE calls AS SELECT *, 0::UINTEGER transcript_index, 1::UINTEGER alt_index, 0::UINTEGER sample_index, 1::BIGINT phase_set, [1,0]::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before FROM events WHERE event_index = 1")
    conditional <- function(model) sql(paste0("SELECT * FROM duckvep_haplotypes('SELECT * FROM calls', ", q(model), ", hgvs := true)"))
    legacy_prediction <- conditional("curated")
    stopifnot(any(grepl("unsupported", legacy_prediction$prediction_status)),
              any(grepl("curated", legacy_prediction$prediction_reason)))
    stopifnot(identical(conditional("combined"), conditional("combined_copy")))
    cat("PASS: legacy conditional curation unsupported; combined conditional operands preserved\n")
  }
  # Raw wire corruption must fail before a model is installed; the same name remains reusable.
  for (path in c(file.path(directory, "plain.v1"), combined)) {
    raw <- bytes(path)
    v1 <- identical(raw[[9L]], as.raw(1))
    counts <- if (v1) 8L else 9L
    sections <- if (v1) 40L else 44L
    width_start <- 8L * (9L + counts) + 1L
    offset_start <- width_start + 8L * sections
    length_start <- offset_start + 8L * sections
    # Remove one nonempty original array, repack offsets, and authenticate the result.
    # The exon phase and end-phase arrays remain required in both wire versions.
    for (slot in c(0L, if (v1) 24L else 26L, if (v1) 25L else 27L, if (!v1) 15L)) {
      omitted <- raw
      present_at <- 8L * (8L + counts) + 1L + slot %/% 8L
      omitted[[present_at]] <- as.raw(bitwAnd(as.integer(omitted[[present_at]]), bitwXor(255L, 2L^(slot %% 8L))))
      streams <- vector("list", sections)
      at <- number(raw, 25L)
      at <- ceiling(at / 64) * 64
      for (section in seq_len(sections) - 1L) {
        n <- number(raw, length_start + 8L * section)
        start <- number(raw, offset_start + 8L * section)
        streams[[section + 1L]] <- if (section == slot || n == 0) raw() else raw[start + seq_len(n)]
        n <- length(streams[[section + 1L]])
        omitted[offset_start + 8L * section + 0:7] <- u64(at)
        omitted[length_start + 8L * section + 0:7] <- u64(n)
        padding <- ceiling(n / 64) * 64 - n
        streams[[section + 1L]] <- c(streams[[section + 1L]], raw(padding))
        at <- at + n + padding
      }
      omitted[33:40] <- u64(at)
      omitted <- c(head(omitted, ceiling(number(raw, 25L) / 64) * 64), unlist(streams))
      bad <- file.path(directory, "required.snap")
      writeBin(omitted, bad)
      stopifnot(system2(binary, shQuote(bad)) == 0L)
      error <- tryCatch({ restore("retry", bad); NULL }, error = identity)
      message <- if (slot == 15L) "cDNA flag does not match its arrays" else "lacks a required array"
      stopifnot(inherits(error, "error"), grepl(message, conditionMessage(error), fixed = TRUE))
      restore("retry", path)
      drop("retry")
    }
    for (version in if (v1) 2L else 1L) {
      bad <- file.path(directory, "shape.snap")
      damaged <- raw
      damaged[9:16] <- u64(version)
      writeBin(damaged, bad)
      error <- tryCatch({ restore("retry", bad); NULL }, error = identity)
      stopifnot(inherits(error, "error"), grepl("incompatible", conditionMessage(error), fixed = TRUE))
    }
    changes <- list(version = 9L, endian = 17L, header = 25L, file = 33L,
                    checksum = 41L, count = 57L,
                    flags = 8L * (7L + counts) + 1L,
                    presence = 8L * (8L + counts) + 6L,
                    width = width_start, offset = offset_start, length = length_start)
    for (label in names(changes)) {
      damaged <- raw
      at <- changes[[label]]
      damaged[[at]] <- as.raw(bitwXor(as.integer(damaged[[at]]), 128L))
      bad <- file.path(directory, "bad.snap")
      writeBin(damaged, bad)
      error <- tryCatch({ restore("retry", bad); NULL }, error = identity)
      stopifnot(inherits(error, "error"))
      restore("retry", path)
      drop("retry")
    }
    for (at in c(65L, offset_start, length_start)) {
      damaged <- raw
      damaged[at + 0:7] <- rep(as.raw(255), 8L)
      bad <- file.path(directory, "overflow.snap")
      writeBin(damaged, bad)
      error <- tryCatch({ restore("retry", bad); NULL }, error = identity)
      stopifnot(inherits(error, "error"))
      restore("retry", path)
      drop("retry")
    }
    for (length in unique(c(0L, 8L, 55L, 1095L, 1096L, 1199L, length(raw) - 1L))) {
      bad <- file.path(directory, "short.snap")
      writeBin(head(raw, length), bad)
      error <- tryCatch({ restore("retry", bad); NULL }, error = identity)
      stopifnot(inherits(error, "error"))
      restore("retry", path)
      drop("retry")
    }
  }
  for (model in c("plain", "curated", "combined", "combined_copy")) drop(model)
  invisible(gc())
  retained <- budget()
  stopifnot(all(retained$current_bytes[retained$owner %in% c("model", "index", "reference")] == 0))
  for (i in seq_len(8L)) {
    restore("cleanup", combined)
    drop("cleanup")
    invisible(gc())
    stopifnot(identical(retained, budget()))
  }
  cat("PASS: real v1 raw annotation and fingerprints; combined v2 roundtrip; authenticated required arrays; corruption and truncation; zero model/index/reference bytes and stable cleanup budget\n")
}

if (legacy_mode) {
  run(args[[2L]], Sys.getenv("DUCKVEP_LEGACY_EXTENSION"), TRUE)
} else {
  directory <- tempfile("duckvep-snapshot-")
  dir.create(directory)
  legacy <- Sys.getenv("DUCKVEP_LEGACY_EXTENSION")
  if (!nzchar(legacy)) stop("Set DUCKVEP_LEGACY_EXTENSION to a published v1 snapshot writer artifact")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  tryCatch({
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), "write-v1", shQuote(directory)), env = paste0("DUCKVEP_LEGACY_EXTENSION=", legacy))
    stopifnot(status == 0L)
    run(directory, Sys.getenv("DUCKVEP_EXTENSION", "build/release/duckvep.duckdb_extension"), FALSE)
  }, finally = unlink(directory, recursive = TRUE))
}
