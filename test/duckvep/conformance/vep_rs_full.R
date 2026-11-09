#!/usr/bin/env Rscript
# Complete-corpus independent-event diagnostic against captured VEP 116.2 JSON.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || !args[1L] %in% c("GRCh37", "GRCh38")) {
  stop("Usage: Rscript vep_rs_full.R GRCh37|GRCh38 GOLDEN_ROOT MODEL EMPTY_TMP_OUTPUT", call. = FALSE)
}
suppressPackageStartupMessages({ library(DBI); library(Rduckvep) })
start <- proc.time()[["elapsed"]]
assembly <- args[1L]
source_revision <- "70df135effe286bc96c8dc93381bb06f2e2979cf"
paths <- c(vcf = file.path(args[2L], assembly, "variants.vcf"),
           provenance = file.path(args[2L], assembly, "provenance.json"),
           expected = file.path(args[2L], assembly, "expected/json.jsonl.gz"))
model_path <- normalizePath(args[3L], mustWork = TRUE)
out <- args[4L]
parent <- normalizePath(dirname(out), mustWork = TRUE)
tmp <- normalizePath("/tmp", mustWork = TRUE)
if (!(identical(parent, tmp) || startsWith(parent, paste0(tmp, "/"))) ||
    (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))) {
  stop("Output must be an empty directory below /tmp", call. = FALSE)
}
if (file.exists(out) && !dir.exists(out)) stop("Output path is not a directory", call. = FALSE)
dir.create(out, showWarnings = FALSE)
if (dirname(normalizePath(out, mustWork = TRUE)) != parent) stop("Output directory must not be a symlink")
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
hashes <- vapply(paths, sha256, "")
expected_hashes <- switch(assembly,
  GRCh37 = c(vcf = "3f879f671744354bb672900a35240ea132caf0411bc4f5bdf3b9d10bb684d849", provenance = "60a8e5c51f2ef6e5626edab1db4469fbdb6cc1fdb979080cf5a9d85733dac171", expected = "9d20a5735562fbe0e13448a1b73c8fc07998a51d4d4f7dad9ee845572c48d371"),
  GRCh38 = c(vcf = "ca155257893e0ea3036147298c86581179444715027b916454dcf8cfdb7a6167", provenance = "d828a11af6a185e22faa9d10aa8e7f08fc5b7cd58e3e3feedae006c48953b08b", expected = "3e5fef8f2a63ac3fa7e83b715792a513a859f49e58804762553a702fb8ad68c0"))
if (!identical(hashes, expected_hashes)) stop("Golden input hash mismatch", call. = FALSE)
provenance <- jsonlite::fromJSON(paths[["provenance"]], simplifyVector = FALSE)
if (!identical(provenance$perl_cache, paste0("homo_sapiens/116_", assembly)) ||
    !grepl("116.2", provenance$vep_version_reported, fixed = TRUE)) stop("Wrong VEP provenance")

lines <- readLines(paths[["vcf"]], warn = FALSE)
lines <- lines[!startsWith(lines, "#")]
f <- strsplit(lines, "\t", fixed = TRUE)
if (any(lengths(f) < 8L) || anyDuplicated(lines)) stop("VCF records must be complete and unique")
column <- function(i) vapply(f, `[[`, "", i)
info_value <- function(text, key) {
  x <- strsplit(text, ";", fixed = TRUE)[[1L]]
  x <- x[startsWith(x, paste0(key, "="))]
  if (length(x) > 0L) substring(x[1L], nchar(key) + 2L) else NA_character_
}
records <- data.frame(record_index = seq_along(lines), chrom = column(1L),
  pos = as.numeric(column(2L)), original_id = column(3L), ref = column(4L),
  alt = column(5L), info = column(8L), input = lines, stringsAsFactors = FALSE)
records$info_end <- vapply(records$info, info_value, "", key = "END")
records$info_svtype <- vapply(records$info, info_value, "", key = "SVTYPE")
records$info_svlen <- vapply(records$info, info_value, "", key = "SVLEN")
alts <- strsplit(records$alt, ",", fixed = TRUE)
alleles <- data.frame(record_index = rep(records$record_index, lengths(alts)),
  alt_index = sequence(lengths(alts)), source_alt = unlist(alts, use.names = FALSE))
alleles$allele_id <- seq_len(nrow(alleles))
# The builder receives each original ALT list, not separately split structural records.
con <- rduckvep_connect()
dbExecute(con, "SET threads = 2")
dbExecute(con, "SET memory_limit = '4GB'")
dbExecute(con, paste("ATTACH", dbQuoteString(con, model_path), "AS golden_model (READ_ONLY)"))
model_receipt <- dbGetQuery(con, "SELECT source_version, assembly, model_sha256, reference_sha256 FROM golden_model.model_receipt")
if (nrow(model_receipt) != 1L || model_receipt$assembly[1L] != assembly || model_receipt$source_version[1L] != "116") stop("Wrong model receipt")
tables <- dbGetQuery(con, "SELECT table_name FROM information_schema.tables WHERE table_catalog = 'golden_model'")$table_name
if (all(c("duckvep_sequence_regions", "duckvep_transcripts", "duckvep_exons", "duckvep_transcript_names") %in% tables)) {
  regions <- "SELECT seq_region, sequence_length FROM golden_model.duckvep_sequence_regions ORDER BY seq_region"
  region_names <- "SELECT seq_region, name AS seq_region_name FROM golden_model.duckvep_sequence_regions"
  transcripts <- "SELECT * FROM golden_model.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index"
  exons <- "SELECT * FROM golden_model.duckvep_exons ORDER BY transcript_index, exon_start"
  names <- "SELECT transcript_index, transcript_id FROM golden_model.duckvep_transcript_names"
} else if (all(c("model_regions", "model_transcripts") %in% tables)) {
  regions <- "SELECT seq_region, sequence_length FROM golden_model.model_regions ORDER BY seq_region"
  region_names <- "SELECT seq_region, seq_region_name FROM golden_model.model_regions"
  transcripts <- paste("SELECT transcript_index, seq_region, transcript_start, transcript_end, strand,",
    "gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table,",
    "pre_cds_sequence, post_cds_sequence FROM golden_model.model_transcripts ORDER BY seq_region, transcript_start, transcript_index")
  exons <- paste("SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start,",
    "e.exon_cdna_end, e.phase, e.end_phase FROM golden_model.model_transcripts,",
    "unnest(exons) u(e) ORDER BY transcript_index, e.exon_cdna_start")
  names <- "SELECT transcript_index, transcript_stable_id AS transcript_id FROM golden_model.model_transcripts"
} else stop("Unsupported model schema")
columns <- dbGetQuery(con, "SELECT column_name FROM information_schema.columns WHERE table_catalog = 'golden_model' AND table_name = 'model_transcripts'")$column_name
mirna <- if ("duckvep_mature_mirna" %in% tables) {
  "SELECT * FROM golden_model.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start"
} else if ("mature_mirna_regions" %in% columns) {
  "SELECT transcript_index, m.mature_mirna_start, m.mature_mirna_end FROM golden_model.model_transcripts, unnest(mature_mirna_regions) u(m) ORDER BY transcript_index, m.mature_mirna_start"
} else NULL
peptide <- if ("duckvep_peptide_edits" %in% tables) {
  "SELECT * FROM golden_model.duckvep_peptide_edits ORDER BY transcript_index, protein_position"
} else if ("peptide_edits" %in% columns) {
  "SELECT transcript_index, p.protein_position, p.alternate_amino_acid FROM golden_model.model_transcripts, unnest(peptide_edits) u(p) ORDER BY transcript_index, p.protein_position"
} else NULL
load_start <- proc.time()[["elapsed"]]
rduckvep_load_model(con, "golden", regions, transcripts, exons,
                    mature_mirna_query = mirna, peptide_edit_query = peptide)
load_seconds <- proc.time()[["elapsed"]] - load_start
dbExecute(con, paste("CREATE TEMP VIEW golden_regions AS", region_names))
dbExecute(con, paste("CREATE TEMP VIEW golden_names AS", names))
dbWriteTable(con, "golden_records", records[, c("record_index", "chrom", "pos", "ref", "alt", "info_end", "info_svtype", "info_svlen")])
dbExecute(con, "CREATE TEMP VIEW golden_vcf AS SELECT chrom, pos, record_index::VARCHAR AS id, ref, string_split(alt, ',') AS alt, info_end, info_svtype, info_svlen FROM golden_records")
events_sql <- rduckvep_vcf_events_sql(con, "golden_vcf", "golden_regions")
prepare_start <- proc.time()[["elapsed"]]
dbExecute(con, paste("CREATE TEMP TABLE golden_events AS", events_sql))
events <- dbGetQuery(con, "SELECT * FROM golden_events ORDER BY source_id::INTEGER, alt_index")
events$record_index <- as.integer(events$source_id)
if (nrow(events) != nrow(alleles) || anyDuplicated(events[c("record_index", "alt_index")]) ||
    !identical(events$source_alt, alleles$source_alt) ||
    !identical(events$record_index, alleles$record_index) ||
    !all(events$alt_index == alleles$alt_index)) stop("Event-to-source allele accounting failed")
prepare_seconds <- proc.time()[["elapsed"]] - prepare_start
dbExecute(con, "CREATE TEMP VIEW golden_supported AS SELECT * FROM golden_events WHERE event_status = 'ok'")
annotate_sql <- rduckvep_annotate_sql(con, "golden_supported", "golden", rich = TRUE)
annotate_start <- proc.time()[["elapsed"]]
observed <- dbGetQuery(con, paste("SELECT e.source_id::INTEGER AS record_index, e.alt_index,",
  "e.source_alt, n.transcript_id, a.* FROM (", annotate_sql, ") a",
  "JOIN golden_events e USING (event_index) LEFT JOIN golden_names n USING (transcript_index)"))
annotate_seconds <- proc.time()[["elapsed"]] - annotate_start
if (any(!observed$event_index %in% events$event_index)) stop("Annotation event identity failed")
dbDisconnect(con, shutdown = TRUE)

json_con <- gzfile(paths[["expected"]], "rt")
json_lines <- readLines(json_con, warn = FALSE)
close(json_con)
golden <- lapply(json_lines, jsonlite::fromJSON, simplifyVector = FALSE)
golden_input <- vapply(golden, `[[`, "", "input")
if (anyDuplicated(golden_input) || any(!golden_input %in% lines)) stop("VEP JSON contains duplicate or unknown input records")
records$vep_output_present <- lines %in% golden_input
golden <- golden[match(lines, golden_input)]
scalar <- function(x, name) if (is.null(x[[name]])) NA_character_ else as.character(x[[name]])
normal_terms <- function(x) paste(sort(unlist(x, use.names = FALSE)), collapse = "&")
reference <- do.call(rbind, lapply(seq_along(golden), function(i) {
  tx <- golden[[i]]$transcript_consequences
  if (length(tx) == 0L) return(NULL)
  data.frame(record_index = i, transcript_id = vapply(tx, scalar, "", name = "transcript_id"),
    variant_allele = vapply(tx, scalar, "", name = "variant_allele"),
    terms = vapply(tx, function(x) normal_terms(x$consequence_terms), ""), stringsAsFactors = FALSE)
}))
if (is.null(reference)) reference <- data.frame(record_index = integer(), transcript_id = character(), variant_allele = character(), terms = character())
# Single ALT establishes provenance without guessing VEP's normalized allele spelling.
# For multiallelic SNPs only, an exact A/C/G/T variant_allele uniquely identifies ALT.
reference$alt_index <- NA_integer_
for (i in seq_along(alts)) {
  ix <- which(reference$record_index == i)
  if (length(alts[[i]]) == 1L) reference$alt_index[ix] <- 1L
  else if (grepl("^[ACGT]$", records$ref[i]) && all(grepl("^[ACGT]$", alts[[i]])) && !anyDuplicated(alts[[i]])) {
    reference$alt_index[ix] <- match(reference$variant_allele[ix], alts[[i]])
  }
}
reference$allele_mapping <- ifelse(is.na(reference$alt_index), "uncomparable_multiallelic", "source_allele_identified")
observed$terms <- vapply(strsplit(observed$consequence, "&", fixed = TRUE), normal_terms, "")
# No-feature and unmapped rows are retained in the full observed artifact, not made into transcript tuples.
obs_tx <- observed[!is.na(observed$transcript_id), c("record_index", "alt_index", "transcript_id", "terms"), drop = FALSE]
ref_tx <- reference[!is.na(reference$alt_index) & !is.na(reference$transcript_id), c("record_index", "alt_index", "transcript_id", "terms"), drop = FALSE]
# Unsupported events and unidentified VEP alleles remain explicit, not differences or matches.
event_key <- paste(events$record_index, events$alt_index, sep = ":")
ref_key <- paste(ref_tx$record_index, ref_tx$alt_index, sep = ":")
ref_tx <- ref_tx[events$event_status[match(ref_key, event_key)] == "ok", , drop = FALSE]
# A multiallelic record with unresolved VEP tuples cannot support allele-specific comparison.
unresolved_records <- unique(reference$record_index[is.na(reference$alt_index) | is.na(reference$transcript_id)])
ref_tx <- ref_tx[!ref_tx$record_index %in% unresolved_records, , drop = FALSE]
obs_tx <- obs_tx[!obs_tx$record_index %in% unresolved_records & records$vep_output_present[obs_tx$record_index], , drop = FALSE]
alleles$event_status <- events$event_status
alleles$event_reason <- events$event_reason
alleles$comparison_scope <- ifelse(!records$vep_output_present[alleles$record_index], "missing_vep_record",
  ifelse(events$event_status != "ok", "unsupported_event",
    ifelse(alleles$record_index %in% unresolved_records, "unresolved_vep_allele_or_transcript", "allele_agnostic_transcript_projection")))
alleles$vep_mapped_transcripts <- tabulate(match(paste(reference$record_index, reference$alt_index, sep = ":"), event_key), nrow(alleles))
alleles$duck_annotations <- tabulate(match(paste(observed$record_index, observed$alt_index, sep = ":"), event_key), nrow(alleles))
key <- function(x) paste(x$record_index, x$alt_index, x$transcript_id, x$terms, sep = "\t")
rk <- key(ref_tx); ok <- key(obs_tx)
keys <- sort(unique(c(rk, ok)))
differences <- data.frame(tuple = keys, vep_count = tabulate(match(rk, keys), length(keys)), duckvep_count = tabulate(match(ok, keys), length(keys)))
differences <- differences[differences$vep_count != differences$duckvep_count, , drop = FALSE]
# Geometry equality is meaningful here only for unambiguous literal SNPs.
vep_records <- data.frame(record_index = seq_along(golden),
  vep_seq_region = vapply(golden, scalar, "", name = "seq_region_name"),
  vep_start = vapply(golden, scalar, "", name = "start"),
  vep_end = vapply(golden, scalar, "", name = "end"),
  vep_allele_string = vapply(golden, scalar, "", name = "allele_string"),
  vep_most_severe = vapply(golden, scalar, "", name = "most_severe_consequence"),
  vep_output_present = records$vep_output_present,
  transcript_rows = ifelse(records$vep_output_present, lengths(lapply(golden, `[[`, "transcript_consequences")), NA_integer_),
  other_features_json = vapply(golden, function(x) as.character(jsonlite::toJSON(x[grepl("_consequences$", names(x)) & names(x) != "transcript_consequences"], auto_unbox = TRUE, null = "null")), ""))
geometry <- merge(events, vep_records, by = "record_index", sort = FALSE)
geometry$comparison <- "uncomparable_geometry_semantics"
snp <- geometry$vep_output_present & grepl("^[ACGT]$", records$ref[geometry$record_index]) & grepl("^[ACGT]$", geometry$source_alt) & geometry$event_status == "ok"
geometry$comparison[snp] <- ifelse(geometry$source_chrom[snp] == geometry$vep_seq_region[snp] &
  as.character(geometry$position[snp]) == geometry$vep_start[snp] & geometry$vep_start[snp] == geometry$vep_end[snp], "equal_snp", "different_snp")
reference$variant_allele_comparison <- "uncomparable_normalization"
snp_ref <- !is.na(reference$alt_index) & grepl("^[ACGT]$", records$ref[reference$record_index])
idx <- match(paste(reference$record_index, reference$alt_index, sep = ":"), event_key)
snp_ref <- snp_ref & !is.na(idx) & grepl("^[ACGT]$", events$source_alt[idx])
reference$variant_allele_comparison[snp_ref] <- ifelse(reference$variant_allele[snp_ref] == events$source_alt[idx[snp_ref]], "equal_snp", "different_snp")
bnd_idx <- which(!is.na(idx) & events$event_kind[idx] == "breakend" & events$event_status[idx] == "ok")
reference$vep_breakend_role <- NA_character_
reference$vep_breakend_role[bnd_idx] <- ifelse(reference$variant_allele[bnd_idx] == events$source_alt[idx[bnd_idx]], "source_alt",
  ifelse(grepl("^\\.|\\.$", reference$variant_allele[bnd_idx]), "single_breakend", "other"))
breakend_alleles <- reference[bnd_idx, c("record_index", "alt_index", "transcript_id", "terms", "variant_allele", "vep_breakend_role"), drop = FALSE]
breakend_alleles$source_alt <- events$source_alt[idx[bnd_idx]]
breakend_alleles$seq_region <- events$seq_region[idx[bnd_idx]]
breakend_alleles$position <- events$position[idx[bnd_idx]]
breakend_alleles$mate_seq_region <- events$mate_seq_region[idx[bnd_idx]]
breakend_alleles$mate_position <- events$mate_position[idx[bnd_idx]]
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
write_tsv(records, "records.tsv"); write_tsv(alleles, "alleles.tsv"); write_tsv(events, "events.tsv")
write_tsv(vep_records, "vep-records.tsv"); write_tsv(reference, "vep-transcripts.tsv")
write_tsv(observed, "duckvep-annotations.tsv"); write_tsv(geometry, "geometry.tsv"); write_tsv(differences, "differences.tsv")
write_tsv(breakend_alleles, "breakend-alleles.tsv")
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]), mustWork = TRUE)
git_args <- c("-C", shQuote(dirname(script_path)))
git_head <- system2("git", c(git_args, "rev-parse", "HEAD"), stdout = TRUE)
git_status <- system2("git", c(git_args, "status", "--porcelain"), stdout = TRUE)
git_diff <- system2("git", c(git_args, "diff", "--binary", "HEAD"), stdout = TRUE)
extension_path <- system.file("duckvep_extension", "artifacts", paste0("v", as.character(packageVersion("duckdb"))), "duckvep.duckdb_extension", package = "Rduckvep", mustWork = TRUE)
rss <- if (file.exists("/proc/self/status")) grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE) else character()
counts <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "<NA>"
  as.list(table(x))
}
receipt <- list(scope = "complete VCF accounting; allele-agnostic transcript-term projection, not full VEP parity",
  assembly = assembly, vep_rs_revision = source_revision, golden_sha256 = as.list(hashes),
  vep_image = provenance$vep_image, vep_version = provenance$vep_version_reported, vep_command = provenance$vep_command,
  model_receipt = model_receipt, model_file_sha256 = sha256(model_path), extension_path = extension_path, extension_sha256 = sha256(extension_path),
  script_sha256 = sha256(script_path), git_head = git_head, git_status = git_status,
  git_diff_sha256 = digest::digest(paste(git_diff, collapse = "\n"), algo = "sha256", serialize = FALSE),
  host = as.list(Sys.info()), r_version = R.version.string, duckdb_version = as.character(packageVersion("duckdb")), duckvep_version = as.character(packageVersion("Rduckvep")),
  limits = list(duckdb_threads = 2, duckdb_memory_limit = "4GB", external_timeout_required = TRUE,
    transcript_comparison = "mapped input alleles and supported events only; multiallelic records with any unresolved VEP tuple excluded as a whole; VEP breakend endpoint alleles are collapsed; multiplicity retained",
    variant_allele = "literal SNPs only; other representations uncomparable", geometry = "literal SNPs only; other event coordinate conventions uncomparable",
    other_features = "VEP feature arrays retained but not compared; interval features not loaded", vep = "captured output; VEP not run locally",
    mature_mirna_loaded = !is.null(mirna), peptide_edits_loaded = !is.null(peptide), reference_fasta_loaded = FALSE),
  counts = list(records = nrow(records), alleles = nrow(alleles), event_status = counts(events$event_status), unsupported_reasons = counts(events$event_reason[events$event_status != "ok"]),
    vep_present_records = sum(records$vep_output_present), vep_missing_records = sum(!records$vep_output_present),
    vep_transcript_rows = nrow(reference), vep_no_transcript_records = sum(vep_records$transcript_rows == 0L, na.rm = TRUE),
    vep_allele_mapping = counts(reference$allele_mapping), skipped_unresolved_records = length(unresolved_records),
    duck_annotations = nrow(observed), duck_status = counts(observed$duckvep_status), duck_reason = counts(observed$duckvep_reason),
    duck_no_or_unmapped_transcript_rows = sum(is.na(observed$transcript_id)), projected_vep_tuples = nrow(ref_tx), projected_duck_tuples = nrow(obs_tx),
    tuple_intersection = sum(pmin(tabulate(match(rk, keys), length(keys)), tabulate(match(ok, keys), length(keys)))),
    allele_comparison_scope = counts(alleles$comparison_scope),
    vep_unmapped_or_unsupported_transcript_rows = nrow(reference) - nrow(ref_tx),
    vep_only = sum(pmax(differences$vep_count - differences$duckvep_count, 0L)), duck_only = sum(pmax(differences$duckvep_count - differences$vep_count, 0L)),
    geometry = counts(geometry$comparison), variant_allele = counts(reference$variant_allele_comparison),
    vep_breakend_roles = counts(reference$vep_breakend_role[bnd_idx])),
  timings_seconds = list(model_load = unname(load_seconds), event_preparation = unname(prepare_seconds), annotation = unname(annotate_seconds), elapsed_before_receipt = unname(proc.time()[["elapsed"]] - start)),
  peak_rss = if (length(rss) > 0L) rss[1L] else NA_character_)
jsonlite::write_json(receipt, file.path(out, "receipt.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
cat(assembly, ":", nrow(records), "records,", nrow(alleles), "alleles,", nrow(differences), "differing allele-agnostic tuple keys; see receipt\n")
