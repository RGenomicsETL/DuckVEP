#!/usr/bin/env Rscript
# Compare lengthless symbolic insertion subtypes and interchromosomal BNDs
# with the VEP 116.2 output captured by vep-rs at the revision below.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || !args[[1L]] %in% c("GRCh37", "GRCh38")) {
  stop("Usage: Rscript vep_rs_sv_boundary.R GRCh37|GRCh38 GOLDEN_ROOT MODEL OUTPUT_DIR",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(DBI)
  library(Rduckvep)
})

assembly <- args[[1L]]
golden_dir <- file.path(args[[2L]], assembly)
model_path <- normalizePath(args[[3L]], mustWork = TRUE)
out_dir <- args[[4L]]
paths <- c(vcf = file.path(golden_dir, "variants.vcf"),
           provenance = file.path(golden_dir, "provenance.json"),
           expected = file.path(golden_dir, "expected", "json.jsonl.gz"))
source_revision <- "70df135effe286bc96c8dc93381bb06f2e2979cf"
expected_hashes <- switch(assembly,
  GRCh37 = c(vcf = "3f879f671744354bb672900a35240ea132caf0411bc4f5bdf3b9d10bb684d849",
             provenance = "60a8e5c51f2ef6e5626edab1db4469fbdb6cc1fdb979080cf5a9d85733dac171",
             expected = "9d20a5735562fbe0e13448a1b73c8fc07998a51d4d4f7dad9ee845572c48d371"),
  GRCh38 = c(vcf = "ca155257893e0ea3036147298c86581179444715027b916454dcf8cfdb7a6167",
             provenance = "d828a11af6a185e22faa9d10aa8e7f08fc5b7cd58e3e3feedae006c48953b08b",
             expected = "3e5fef8f2a63ac3fa7e83b715792a513a859f49e58804762553a702fb8ad68c0"))
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
actual_hashes <- vapply(paths, sha256, "")
if (!identical(actual_hashes, expected_hashes)) {
  stop("Golden input hashes do not match vep-rs revision ", source_revision, call. = FALSE)
}
provenance <- jsonlite::fromJSON(paths[["provenance"]], simplifyVector = FALSE)
if (!identical(provenance$perl_cache, paste0("homo_sapiens/116_", assembly)) ||
    !grepl("116.2", provenance$vep_version_reported, fixed = TRUE)) {
  stop("Golden provenance has the wrong assembly or VEP version", call. = FALSE)
}

start <- proc.time()[["elapsed"]]
vcf <- readLines(paths[["vcf"]], warn = FALSE)
vcf <- vcf[!startsWith(vcf, "#")]
fields <- strsplit(vcf, "\t", fixed = TRUE)
if (any(lengths(fields) < 8L) || anyDuplicated(vcf)) {
  stop("Golden VCF records must be complete and unique", call. = FALSE)
}
chrom <- vapply(fields, `[[`, "", 1L)
alt <- vapply(fields, `[[`, "", 5L)
info <- vapply(fields, `[[`, "", 8L)
subtype <- grepl("^<INS:[A-Za-z0-9_.-]+(:[A-Za-z0-9_.-]+)*>$", alt) &
  !grepl("(^|;)(END|SVLEN)=", info)
matches <- regexpr("[\\[\\]][^\\[\\]]+:[0-9]+[\\[\\]]", alt, perl = TRUE)
mate <- rep(NA_character_, length(alt))
has_mate <- matches > 0L
mate[has_mate] <- substring(alt[has_mate], matches[has_mate] + 1L,
                            matches[has_mate] + attr(matches, "match.length")[has_mate] - 2L)
mate[has_mate] <- sub(":[0-9]+$", "", mate[has_mate])
cross_contig <- has_mate & mate != chrom
selected <- which(subtype | cross_contig)
if (length(selected) == 0L || any(grepl(",", alt[selected], fixed = TRUE))) {
  stop("Boundary corpus is empty or contains a multiallelic selected record", call. = FALSE)
}
info_value <- function(text, key) {
  entry <- strsplit(text, ";", fixed = TRUE)[[1L]]
  entry <- entry[startsWith(entry, paste0(key, "="))]
  if (length(entry) == 0L) return(NA_character_)
  substring(entry[[1L]], nchar(key) + 2L)
}
records <- data.frame(
  record_index = selected,
  chrom = chrom[selected],
  pos = as.numeric(vapply(fields[selected], `[[`, "", 2L)),
  id = as.character(selected),
  original_id = vapply(fields[selected], `[[`, "", 3L),
  ref = vapply(fields[selected], `[[`, "", 4L),
  alt = alt[selected],
  info_end = vapply(info[selected], info_value, "", key = "END"),
  info_svtype = vapply(info[selected], info_value, "", key = "SVTYPE"),
  info_svlen = vapply(info[selected], info_value, "", key = "SVLEN"),
  input = vcf[selected],
  stringsAsFactors = FALSE
)

con <- rduckvep_connect()
dbExecute(con, paste("ATTACH", dbQuoteString(con, model_path), "AS golden_model (READ_ONLY)"))
model_receipt <- dbGetQuery(con, "SELECT source_version, assembly, model_sha256, reference_sha256 FROM golden_model.model_receipt")
if (nrow(model_receipt) != 1L || model_receipt$assembly[[1L]] != assembly ||
    model_receipt$source_version[[1L]] != "116") {
  stop("Model receipt must identify the selected assembly and Ensembl release 116", call. = FALSE)
}
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
    "pre_cds_sequence, post_cds_sequence FROM golden_model.model_transcripts",
    "ORDER BY seq_region, transcript_start, transcript_index")
  exons <- paste("SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start,",
    "e.exon_cdna_end, e.phase, e.end_phase FROM golden_model.model_transcripts,",
    "unnest(exons) u(e) ORDER BY transcript_index, e.exon_cdna_start")
  names <- "SELECT transcript_index, transcript_stable_id AS transcript_id FROM golden_model.model_transcripts"
} else {
  stop("Model does not expose a supported release-116 transcript schema", call. = FALSE)
}
rduckvep_load_model(con, "golden", regions, transcripts, exons)
dbExecute(con, paste("CREATE TEMP VIEW golden_regions AS", region_names))
dbExecute(con, paste("CREATE TEMP VIEW golden_names AS", names))
dbWriteTable(con, "golden_records", records)
dbExecute(con, paste("CREATE TEMP VIEW golden_vcf AS SELECT chrom, pos, id, ref, [alt] AS alt,",
  "info_end, info_svtype, info_svlen FROM golden_records"))
events_sql <- rduckvep_vcf_events_sql(con, "golden_vcf", "golden_regions")
dbExecute(con, paste("CREATE TEMP TABLE golden_events AS", events_sql))
events <- dbGetQuery(con, "SELECT * FROM golden_events ORDER BY source_id::INTEGER")
if (nrow(events) != nrow(records) || any(events$event_status != "ok")) {
  stop("At least one selected record was not prepared as a supported event", call. = FALSE)
}
annotate_sql <- rduckvep_annotate_sql(con, "golden_events", "golden", rich = TRUE)
observed <- dbGetQuery(con, paste("SELECT e.source_id::INTEGER AS record_index,",
  "n.transcript_id, a.consequence FROM (", annotate_sql, ") a",
  "JOIN golden_events e USING (event_index)",
  "JOIN golden_names n USING (transcript_index)"))
normalize_terms <- function(x) paste(sort(unlist(x, use.names = FALSE)), collapse = "&")
observed$terms <- vapply(strsplit(observed$consequence, "&", fixed = TRUE), normalize_terms, "")
observed$consequence <- NULL

json_con <- gzfile(paths[["expected"]], "rt")
json_lines <- readLines(json_con, warn = FALSE)
close(json_con)
golden <- lapply(json_lines, jsonlite::fromJSON, simplifyVector = FALSE)
golden_input <- vapply(golden, `[[`, "", "input")
if (anyDuplicated(golden_input) || any(!golden_input %in% vcf)) {
  stop("VEP JSON input records do not uniquely match the golden VCF", call. = FALSE)
}
if (any(!records$input %in% golden_input)) {
  stop("Selected VCF records are missing from the VEP JSON", call. = FALSE)
}
reference <- lapply(seq_len(nrow(records)), function(i) {
  index <- match(records$input[[i]], golden_input)
  tx <- golden[[index]]$transcript_consequences
  if (is.null(tx) || length(tx) == 0L) return(NULL)
  data.frame(record_index = records$record_index[[i]],
    transcript_id = vapply(tx, `[[`, "", "transcript_id"),
    variant_allele = vapply(tx, `[[`, "", "variant_allele"),
    terms = vapply(tx, function(x) normalize_terms(x$consequence_terms), ""),
    stringsAsFactors = FALSE)
})
reference <- do.call(rbind, reference)
if (is.null(reference)) reference <- data.frame(record_index = integer(), transcript_id = character(),
                                                variant_allele = character(), terms = character())
key <- function(x) paste(x$record_index, x$transcript_id, x$terms, sep = "\t")
reference_key <- key(reference)
observed_key <- key(observed)
keys <- sort(unique(c(reference_key, observed_key)))
reference_count <- tabulate(match(reference_key, keys), length(keys))
observed_count <- tabulate(match(observed_key, keys), length(keys))
differences <- data.frame(tuple = keys, vep_count = reference_count, duckvep_count = observed_count)
differences <- differences[differences$vep_count != differences$duckvep_count, , drop = FALSE]

if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Output directory must be empty", call. = FALSE)
}
dbDisconnect(con, shutdown = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
git <- Sys.which("git")
git_head <- NA_character_
git_status <- character()
git_diff_sha256 <- NA_character_
if (nzchar(git)) {
  git_args <- c("-C", shQuote(dirname(script_path)))
  head <- suppressWarnings(system2(git, c(git_args, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE))
  if (length(head) == 1L && grepl("^[0-9a-f]{40}$", head)) {
    git_head <- head
    git_status <- system2(git, c(git_args, "status", "--porcelain"), stdout = TRUE)
    git_diff <- system2(git, c(git_args, "diff", "--binary", "HEAD"), stdout = TRUE)
    git_diff_sha256 <- digest::digest(paste(git_diff, collapse = "\n"), algo = "sha256", serialize = FALSE)
  }
}
extension_path <- system.file("duckvep_extension", "artifacts",
  paste0("v", as.character(packageVersion("duckdb"))), "duckvep.duckdb_extension",
  package = "Rduckvep", mustWork = TRUE)
status <- if (file.exists("/proc/self/status")) readLines("/proc/self/status", warn = FALSE) else character()
rss <- grep("^VmHWM:", status, value = TRUE)
receipt <- list(scope = "lengthless insertion subtypes and cross-contig breakends",
  assembly = assembly, vep_rs_revision = source_revision, golden_sha256 = as.list(actual_hashes),
  vep_image = provenance$vep_image, vep_version = provenance$vep_version_reported,
  vep_command = provenance$vep_command, model_receipt = model_receipt,
  model_file_sha256 = sha256(model_path), extension_sha256 = sha256(extension_path),
  script_sha256 = sha256(script_path), git_head = git_head, git_status = git_status,
  git_diff_sha256 = git_diff_sha256, duckdb_version = as.character(packageVersion("duckdb")),
  duckvep_version = as.character(packageVersion("Rduckvep")),
  record_count = nrow(records), vep_tuple_count = nrow(reference),
  duckvep_tuple_count = nrow(observed), differing_tuple_count = nrow(differences),
  elapsed_seconds = unname(proc.time()[["elapsed"]] - start),
  peak_rss = if (length(rss) > 0L) rss[[1L]] else NA_character_)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(records, file.path(out_dir, "records.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(events, file.path(out_dir, "events.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(reference, file.path(out_dir, "vep-tuples.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(observed, file.path(out_dir, "duckvep-tuples.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(differences, file.path(out_dir, "differences.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
jsonlite::write_json(receipt, file.path(out_dir, "receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(assembly, ":", nrow(records), "records;", nrow(reference), "VEP tuples;",
    nrow(observed), "DuckVEP tuples;", nrow(differences), "differences\n")
if (nrow(differences)) quit(status = 1L)
