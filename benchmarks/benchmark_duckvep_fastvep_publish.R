#!/usr/bin/env Rscript

# Mechanical field-campaign packaging; comparison values and verdicts are copied.
duckvep_fastvep_publish <- function(source, output) {
  source <- normalizePath(source, mustWork = TRUE)
  output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
  if (!dir.exists(source) || file.exists(output) || startsWith(output, paste0(source, "/"))) {
    stop("publication requires an existing campaign and a new external output directory")
  }
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  read_table <- function(path) utils::read.delim(path, colClasses = "character", quote = "", comment.char = "")
  write_table <- function(value, path) utils::write.table(value, path,
    sep = "\t", quote = FALSE, row.names = FALSE)
  manifest_path <- file.path(source, "artifacts.tsv")
  manifest_hash <- hash(manifest_path)
  manifest <- read_table(manifest_path)
  if (!identical(names(manifest), c("path", "sha256")) || anyNA(manifest) ||
      anyDuplicated(manifest$path) || any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("invalid execution artifact manifest")
  }
  internal <- manifest[startsWith(manifest$path, paste0(source, "/")), ]
  relative <- substring(internal$path, nchar(source) + 2L)
  required <- c("receipt.tsv", "summary.csv", "cases.csv", "errors.csv", "commands.csv",
    "environment.txt", "inputs_before.tsv", "inputs_after.tsv")
  if (!all(required %in% relative) || any(grepl("(^|/)\\.\\.(/|$)|[\r\n\t]", relative)) ||
      !setequal(relative, setdiff(list.files(source, recursive = TRUE), "artifacts.tsv")) ||
      any(!file.exists(internal$path)) || any(nzchar(Sys.readlink(internal$path)))) {
    stop("execution manifest does not cover the campaign's regular files")
  }
  if (!identical(unname(vapply(internal$path, hash, character(1L))), internal$sha256)) {
    stop("campaign artifact differs from its execution digest")
  }
  inputs <- read_table(file.path(source, "inputs_before.tsv"))
  ids <- c("projection_reference", "projection_reference_fai", "projection_model_gff")
  if (!all(c("id", "path", "sha256") %in% names(inputs)) || anyDuplicated(inputs$id) ||
      !all(ids %in% inputs$id)) stop("campaign lacks the three registered fixture inputs")
  inputs <- inputs[match(ids, inputs$id), c("id", "path", "sha256")]
  input_names <- file.path("inputs", basename(inputs$path))
  index <- match(inputs$path, manifest$path)
  if (anyNA(inputs) || anyNA(index) || anyDuplicated(input_names) ||
      !identical(inputs$sha256, manifest$sha256[index]) ||
      !identical(unname(vapply(inputs$path, hash, character(1L))), inputs$sha256)) {
    stop("registered fixture input differs from its execution digest")
  }
  keep <- !grepl("(^|/)(model[.]duckdb|fastvep[.]cache|vep_home)(/|$)", relative)
  duplicates <- data.frame(file = character(), retained_file = character(), sha256 = character())
  for (i in which(basename(relative) == "generated.vcf")) {
    counterpart <- file.path(dirname(relative[[i]]), "input.vcf")
    j <- match(counterpart, relative)
    if (!is.na(j) && identical(internal$sha256[[i]], internal$sha256[[j]])) {
      keep[[i]] <- FALSE
      duplicates <- rbind(duplicates, data.frame(file = relative[[i]],
        retained_file = counterpart, sha256 = internal$sha256[[i]]))
    }
  }
  internal <- internal[keep, ]
  relative <- relative[keep]
  compress <- basename(relative) %in% c("duckvep_native_tab17.tsv", "duckvep_vep_csq.tsv",
    "fastvep_native_tab17.tsv", "fastvep_vep_csq.vcf", "vep.vcf")
  packed <- paste0(relative, ifelse(compress, ".gz", ""))
  reserved <- c("artifacts.tsv", "execution_artifacts.tsv", "compressed_outputs.tsv", "omitted_duplicates.tsv")
  if (anyDuplicated(c(packed, input_names, reserved))) stop("portable artifact names collide")
  if (any(compress) && !nzchar(Sys.which("gzip"))) stop("gzip is required for portable raw outputs")
  staging <- tempfile(".field-pack-", tmpdir = dirname(output))
  if (!dir.create(staging)) stop("could not create publication staging directory")
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  copy_checked <- function(input, name, expected) {
    target <- file.path(staging, name)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(input, target) || !identical(hash(target), expected)) stop("copied artifact digest mismatch: ", name)
  }
  for (i in seq_len(nrow(internal))) {
    if (!compress[[i]]) {
      copy_checked(internal$path[[i]], packed[[i]], internal$sha256[[i]])
      next
    }
    target <- file.path(staging, packed[[i]])
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (system2("gzip", c("-n", "-c", "--", shQuote(internal$path[[i]])), stdout = target) != 0L) {
      stop("gzip failed: ", relative[[i]])
    }
    verified <- tempfile(".unpacked-", tmpdir = staging)
    if (system2("gzip", c("-d", "-c", "--", shQuote(target)), stdout = verified) != 0L ||
        !identical(hash(verified), internal$sha256[[i]])) stop("gzip changed raw bytes: ", relative[[i]])
    unlink(verified)
  }
  for (i in seq_len(nrow(inputs))) copy_checked(inputs$path[[i]], input_names[[i]], inputs$sha256[[i]])
  copy_checked(manifest_path, "execution_artifacts.tsv", manifest_hash)
  write_table(data.frame(file = relative[compress], original_sha256 = internal$sha256[compress]),
    file.path(staging, "compressed_outputs.tsv"))
  write_table(duplicates, file.path(staging, "omitted_duplicates.tsv"))
  files <- list.files(staging, recursive = TRUE)
  write_table(data.frame(path = files,
    sha256 = unname(vapply(file.path(staging, files), hash, character(1L)))), file.path(staging, "artifacts.tsv"))
  if (!identical(hash(manifest_path), manifest_hash) ||
      !identical(unname(vapply(internal$path, hash, character(1L))), internal$sha256)) {
    stop("campaign changed during publication")
  }
  if (file.exists(output) || !file.rename(staging, output)) stop("could not publish without replacing existing data")
  invisible(output)
}

main <- function() {
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--source"), optparse::make_option("--output"))))
  if (is.null(options$source) || is.null(options$output)) stop("--source and --output are required")
  duckvep_fastvep_publish(options$source, options$output)
  message("Portable field evidence retained: ", options$output)
}

if (sys.nframe() == 0L) main()
