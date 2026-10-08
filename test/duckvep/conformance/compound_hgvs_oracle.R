#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(jsonlite))

fixture_dir <- "test/duckvep/conformance/compound_hgvs_oracle_fixtures"
reference_path <- file.path(fixture_dir, "reference.fa")
model_path <- file.path(fixture_dir, "model.gff3")
cases_path <- file.path(fixture_dir, "cases.jsonl")
adapter_path <- "test/duckvep/conformance/mutalyzer_oracle/normalize_jsonl.py"
lock_path <- "test/duckvep/conformance/mutalyzer_oracle/requirements.lock"
python <- Sys.which(Sys.getenv("MUTALYZER_PYTHON", "python3"))
retriever <- Sys.which(Sys.getenv("MUTALYZER_RETRIEVER", "mutalyzer_retriever"))
artifact_root <- Sys.getenv("DUCKVEP_ARTIFACTS_DIR", "test/duckvep/conformance/results/mutalyzer")

stopifnot(all(file.exists(c(reference_path, model_path, cases_path, adapter_path, lock_path))),
  file.access(python, 1L) == 0L, file.access(retriever, 1L) == 0L)
Sys.setenv(PYTHONDONTWRITEBYTECODE = "1")

cases <- lapply(readLines(cases_path), jsonlite::fromJSON, simplifyVector = FALSE)
cases <- Filter(function(case) !is.null(case$expected_normalized_description), cases)
case_ids <- vapply(cases, `[[`, "", "case_id")
reference_ids <- vapply(cases, `[[`, "", "reference_id")
transcript_ids <- vapply(cases, `[[`, "", "transcript_id")
stopifnot(length(cases) == 6L, !anyDuplicated(case_ids),
  all(vapply(cases, function(x) length(x$edits) == 2L, logical(1L))),
  all(vapply(cases, function(x) x$strand %in% c("+", "-"), logical(1L))))
edits <- unlist(lapply(cases, `[[`, "edits"), recursive = FALSE)
stopifnot(all(vapply(edits, function(edit)
  identical(sort(names(edit)), c("alternate", "position", "reference")), logical(1L))))
positions <- vapply(edits, `[[`, numeric(1L), "position")
references <- vapply(edits, `[[`, "", "reference")
alternates <- vapply(edits, `[[`, "", "alternate")
stopifnot(all(is.finite(positions)), all(positions > 0L),
  all(positions == floor(positions)), all(nzchar(references)),
  all(vapply(cases, function(x) all(vapply(x$edits, function(edit)
    nzchar(edit$reference) || nzchar(edit$alternate), logical(1L))), logical(1L))))

descriptions <- vapply(cases, function(x) {
  if (!is.null(x$input_description)) return(x$input_description)
  edits <- vapply(x$edits, function(edit) paste0(edit$position, edit$reference,
    ">", edit$alternate), "")
  paste0(x$reference_id, "(", x$transcript_id, "):c.[", paste(edits, collapse = ";"), "]")
}, "")

fasta <- readLines(reference_path)
fasta_headers <- which(startsWith(fasta, ">"))
fasta_ends <- c(fasta_headers[-1L] - 1L, length(fasta))
fasta_ids <- substring(fasta[fasta_headers], 2L)
sequences <- vapply(seq_along(fasta_headers), function(i) {
  paste0(fasta[seq.int(fasta_headers[i] + 1L, fasta_ends[i])], collapse = "")
}, "")
names(sequences) <- fasta_ids
stopifnot(!anyDuplicated(fasta_ids), all(nchar(sequences) == 104L),
  all(unique(reference_ids) %in% fasta_ids))

gff <- readLines(model_path)
gff_rows <- gff[nzchar(gff) & !startsWith(gff, "#")]
gff_fields <- strsplit(gff_rows, "\t", fixed = TRUE)
stopifnot(all(lengths(gff_fields) == 9L))
gff_references <- vapply(gff_fields, `[[`, "", 1L)

if (!dir.exists(artifact_root)) dir.create(artifact_root, recursive = TRUE)
run_dir <- file.path(artifact_root, paste0("run-", format(Sys.time(), "%Y%m%dT%H%M%S"),
  "-", Sys.getpid()))
stopifnot(dir.create(run_dir))
cache_dir <- file.path(run_dir, "cache")
stopifnot(dir.create(cache_dir))
settings_path <- file.path(run_dir, "settings.txt")
writeLines(c(paste("MUTALYZER_CACHE_DIR =", cache_dir),
  "MUTALYZER_FILE_CACHE_ADD = false",
  paste("MUTALYZER_LOG_DIR =", file.path(run_dir, "mutalyzer.log"))), settings_path)
Sys.setenv(MUTALYZER_SETTINGS = settings_path,
  MUTALYZER_EXPECTED_CACHE_DIR = cache_dir)
Sys.unsetenv("MUTALYZER_API_URL")
retriever_version <- system2(retriever, "-v", stdout = TRUE, stderr = TRUE)
stopifnot(any(grepl("^mutalyzer_retriever version 0\\.6\\.0$", retriever_version)))

for (reference_id in unique(reference_ids)) {
  selected <- which(reference_ids == reference_id)
  stopifnot(length(unique(transcript_ids[selected])) == 1L,
    length(unique(vapply(cases[selected], `[[`, "", "strand"))) == 1L)
  transcript_id <- transcript_ids[selected[1L]]
  strand <- cases[[selected[1L]]]$strand
  gff_path <- file.path(run_dir, paste0(reference_id, ".gff3"))
  fasta_path <- file.path(run_dir, paste0(reference_id, ".fa"))
  reference_gff <- gff_rows[gff_references == reference_id]
  stopifnot(length(reference_gff) == 4L)
  writeLines(c("##gff-version 3", paste("##sequence-region", reference_id, 1L,
    nchar(sequences[[reference_id]])), reference_gff), gff_path)
  writeLines(c(paste0(">", reference_id), sequences[[reference_id]]), fasta_path)

  stdout_path <- file.path(run_dir, paste0(reference_id, ".retriever.stdout"))
  stderr_path <- file.path(run_dir, paste0(reference_id, ".retriever.stderr"))
  status <- system2(retriever, shQuote(c("--id", reference_id, "--output", cache_dir,
    "--split", "from_file", "--paths", gff_path, fasta_path)),
    stdout = stdout_path, stderr = stderr_path)
  if (status != 0L) stop("Mutalyzer reference-model build failed for ", reference_id,
    "; see ", stderr_path, call. = FALSE)

  annotations_path <- file.path(cache_dir, paste0(reference_id, ".annotations"))
  sequence_path <- file.path(cache_dir, paste0(reference_id, ".sequence"))
  stopifnot(file.exists(annotations_path), file.exists(sequence_path),
    identical(readChar(sequence_path, file.info(sequence_path)$size, useBytes = TRUE),
      sequences[[reference_id]]))
  annotations <- jsonlite::fromJSON(annotations_path, simplifyVector = FALSE)
  mrna <- unlist(lapply(annotations$features, function(gene) {
    Filter(function(feature) identical(feature$type, "mRNA"), gene$features)
  }), recursive = FALSE)
  stopifnot(identical(annotations$id, reference_id), length(mrna) == 1L,
    identical(mrna[[1L]]$id, transcript_id),
    identical(as.integer(mrna[[1L]]$location$strand), if (strand == "+") 1L else -1L))
}

input_path <- file.path(run_dir, "input.jsonl")
writeLines(vapply(seq_along(cases), function(i) jsonlite::toJSON(list(
  case_id = case_ids[i], description = descriptions[i]), auto_unbox = TRUE), ""), input_path)
output_path <- file.path(run_dir, "output.jsonl")
stdout_path <- file.path(run_dir, "normalizer.stdout")
stderr_path <- file.path(run_dir, "normalizer.stderr")
status <- system2(python, shQuote(c(adapter_path, "--input", input_path, "--output",
  output_path, "--lock", lock_path)), stdout = stdout_path, stderr = stderr_path)
if (status != 0L) stop("Mutalyzer normalization failed; see ", stderr_path,
  call. = FALSE)

outputs <- lapply(readLines(output_path), jsonlite::fromJSON, simplifyVector = FALSE)
output_fields <- c("schema", "case_id", "input_description", "normalized_description",
  "protein_description", "protein_reference", "protein_predicted", "errors",
  "mutalyzer_version", "retriever_version", "python_version")
stopifnot(length(outputs) == length(cases), all(vapply(outputs, function(x)
  setequal(names(x), output_fields), logical(1L))))
for (i in seq_along(cases)) {
  result <- outputs[[i]]
  expected <- cases[[i]]
  comparisons <- c(
    schema = identical(result$schema, "mutalyzer-oracle-result/v1"),
    case_id = identical(result$case_id, case_ids[i]),
    input_description = identical(result$input_description, descriptions[i]),
    normalized_description = identical(result$normalized_description,
      expected$expected_normalized_description),
    protein_description = identical(result$protein_description,
      expected$expected_protein_description),
    protein_reference = identical(result$protein_reference,
      expected$expected_protein_reference),
    protein_predicted = identical(result$protein_predicted,
      expected$expected_protein_predicted),
    errors = identical(result$errors, expected$expected_errors),
    mutalyzer_version = identical(result$mutalyzer_version, "3.1.1"),
    retriever_version = identical(result$retriever_version, "0.6.0"),
    python_version = identical(result$python_version, "3.13.12")
  )
  if (!all(comparisons)) stop("Mutalyzer mismatch for ", case_ids[i], ": ",
    paste(names(comparisons)[!comparisons], collapse = ", "), call. = FALSE)
}

fixture_files <- c(reference_path, model_path, cases_path, adapter_path, lock_path)
model_files <- list.files(cache_dir, full.names = TRUE)
hash_files <- c(fixture_files, model_files)
hashes <- system2("sha256sum", shQuote(normalizePath(hash_files)), stdout = TRUE)
stopifnot(length(hashes) == length(hash_files),
  all(is.na(attr(hashes, "status")) | attr(hashes, "status") == 0L))
fixture_hashes <- as.list(setNames(substring(hashes[seq_along(fixture_files)], 1L, 64L),
  basename(fixture_files)))
model_hashes <- as.list(setNames(substring(hashes[length(fixture_files) + seq_along(model_files)],
  1L, 64L), basename(model_files)))
receipt <- list(schema = "duckvep.compound-hgvs-oracle-receipt/v1",
  oracle = "Mutalyzer normalizer API", mutalyzer_version = "3.1.1",
  retriever_version = "0.6.0", python_version = "3.13.12",
  reference_provenance = "Synthetic 104-nt plus/minus references; 72-nt coding region; no assembly accession or external reference cache.",
  fixture_sha256 = fixture_hashes, reference_model_sha256 = model_hashes,
  network_access = "blocked in adapter; API cache disabled",
  file_cache_writes = FALSE, cases = length(cases), run_directory = normalizePath(run_dir))
jsonlite::write_json(receipt, file.path(run_dir, "receipt.json"), pretty = TRUE,
  auto_unbox = TRUE)
message("Mutalyzer compound HGVS oracle passed ", length(cases), " cases: ", run_dir)
