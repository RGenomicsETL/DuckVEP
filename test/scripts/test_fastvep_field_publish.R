#!/usr/bin/env Rscript

# Offline packaging checks use synthetic observations, never an annotation run.
main <- function() {
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_publish.R"))
  work <- tempfile("fastvep-field-publish-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  table <- function(x, path) utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  read_table <- function(path) utils::read.delim(path, colClasses = "character")
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  exact_failure <- "forward/comparison_duckvep_vep_csq/field_failures.parquet"
  body_failure <- "forward/comparison_duckvep_vep_csq/hgvs_body_failures.parquet"
  body_summary <- "forward/comparison_duckvep_vep_csq/hgvs_bodies.csv"
  fixture <- function(name) {
    path <- file.path(work, name)
    dir.create(path)
    input_dir <- paste0(path, "-inputs")
    dir.create(input_dir)
    input_names <- c("minimal.fa", "minimal.fa.fai", "minimal.gff3")
    for (file in input_names) writeLines(paste("fixture", file), file.path(input_dir, file))
    inputs <- data.frame(id = c("projection_reference", "projection_reference_fai", "projection_model_gff"),
      path = file.path(input_dir, input_names))
    inputs$sha256 <- unname(vapply(inputs$path, hash, character(1L)))
    table(inputs, file.path(path, "inputs_before.tsv"))
    table(inputs, file.path(path, "inputs_after.tsv"))
    table(data.frame(field = "binding", value = "diagnostic_unbound"), file.path(path, "receipt.tsv"))
    writeLines(c('"field_failures","passed"', '1,FALSE'), file.path(path, "summary.csv"))
    for (file in c("cases.csv", "errors.csv", "commands.csv", "environment.txt", "build.log")) {
      writeLines(paste("retained", file), file.path(path, file))
    }
    for (profile in c("forward", "reverse")) {
      dir.create(file.path(path, profile))
      for (file in c("input.vcf", "generated.vcf", "model.gff3", "model.duckdb", "fastvep.cache", "vep.log")) {
        writeLines(paste(profile, file), file.path(path, profile, file))
      }
      for (file in c("duckvep_native_tab17.tsv", "duckvep_vep_csq.tsv", "fastvep_native_tab17.tsv",
          "fastvep_vep_csq.vcf", "vep.vcf")) {
        writeLines(c("#header", paste(profile, file, "c.2C>A", "c.1C>A", sep = "\t")), file.path(path, profile, file))
      }
      for (comparison in c("native_tab17", "duckvep_vep_csq", "fastvep_vep_csq")) {
        comparison_dir <- file.path(path, profile, paste0("comparison_", comparison))
        dir.create(comparison_dir)
        failure <- file.path(comparison_dir, "field_failures.parquet")
        DBI::dbExecute(con, paste0("COPY (SELECT 1::UBIGINT record_index,
          1::UBIGINT alt_index, 'ENST1' Feature, 'HGVSc' field,
          'ENST1.1:c.2C>A' actual, 'ENST1.1:c.1C>A' expected) TO ",
          DBI::dbQuoteString(con, failure), " (FORMAT PARQUET)"))
        if (comparison == "native_tab17") next
        body <- file.path(comparison_dir, "hgvs_body_failures.parquet")
        DBI::dbExecute(con, paste0("COPY (SELECT 1::UBIGINT record_index,
          1::UBIGINT alt_index, 'ENST1' Feature, 'HGVSc' field,
          'c.2C>A' actual, 'c.1C>A' expected) TO ",
          DBI::dbQuoteString(con, body), " (FORMAT PARQUET)"))
        utils::write.csv(data.frame(field = "HGVSc", failures = 1L),
          file.path(comparison_dir, "hgvs_bodies.csv"), row.names = FALSE)
      }
    }
    stopifnot(file.copy(file.path(path, "reverse/input.vcf"), file.path(path, "reverse/generated.vcf"), overwrite = TRUE))
    dir.create(file.path(path, "forward/vep_home"))
    writeLines("private cache", file.path(path, "forward/vep_home/private.dat"))
    seal(path)
    path
  }
  seal <- function(path) {
    files <- setdiff(list.files(path, recursive = TRUE, full.names = TRUE), file.path(path, "artifacts.tsv"))
    inputs <- read_table(file.path(path, "inputs_before.tsv"))
    files <- c(inputs$path, files)
    table(data.frame(path = files, sha256 = unname(vapply(files, hash, character(1L)))), file.path(path, "artifacts.tsv"))
  }
  validate <- function(path) {
    manifest <- read_table(file.path(path, "artifacts.tsv"))
    stopifnot(!anyDuplicated(manifest$path), !any(startsWith(manifest$path, "/")),
      !any(grepl("(^|/)\\.\\.(/|$)", manifest$path)),
      identical(unname(vapply(file.path(path, manifest$path), hash, character(1L))), manifest$sha256))
    manifest
  }
  original <- fixture("original")
  first <- file.path(work, "pack-one")
  second <- file.path(work, "pack-two")
  duckvep_fastvep_publish(original, first)
  duckvep_fastvep_publish(original, second)
  stopifnot(identical(validate(first), validate(second)),
    identical(hash(file.path(original, "artifacts.tsv")), hash(file.path(first, "execution_artifacts.tsv"))),
    identical(hash(file.path(original, "summary.csv")), hash(file.path(first, "summary.csv"))),
    identical(hash(file.path(original, exact_failure)), hash(file.path(first, exact_failure))),
    identical(hash(file.path(original, body_failure)), hash(file.path(first, body_failure))),
    identical(hash(file.path(original, body_summary)), hash(file.path(first, body_summary))),
    identical(hash(file.path(original, "forward/generated.vcf")), hash(file.path(first, "forward/generated.vcf"))),
    !file.exists(file.path(first, "reverse/generated.vcf")),
    !any(grepl("model[.]duckdb|fastvep[.]cache|vep_home", list.files(first, recursive = TRUE))))
  body_files <- list.files(first, pattern = "^hgvs_(body_failures[.]parquet|bodies[.]csv)$",
    recursive = TRUE)
  stopifnot(length(body_files) == 8L,
    all(grepl("comparison_(duckvep|fastvep)_vep_csq/", body_files)),
    !any(grepl("comparison_native_tab17/", body_files)))
  omitted <- read_table(file.path(first, "omitted_duplicates.tsv"))
  stopifnot(nrow(omitted) == 1L, omitted$file == "reverse/generated.vcf", omitted$retained_file == "reverse/input.vcf")
  compressed <- read_table(file.path(first, "compressed_outputs.tsv"))
  stopifnot(nrow(compressed) == 10L)
  for (file in compressed$file) {
    connection <- gzfile(file.path(first, paste0(file, ".gz")), "rb")
    actual <- readBin(connection, "raw", n = file.info(file.path(original, file))$size + 1L)
    close(connection)
    expected <- readBin(file.path(original, file), "raw", n = length(actual) + 1L)
    stopifnot(identical(actual, expected))
  }
  inputs <- read_table(file.path(original, "inputs_before.tsv"))
  stopifnot(identical(unname(vapply(file.path(first, "inputs", basename(inputs$path)), hash, character(1L))), inputs$sha256))
  relocated <- file.path(work, "relocated")
  stopifnot(file.rename(first, relocated))
  validate(relocated)
  rejected <- character()
  check <- function(label, mutate, reseal = FALSE) {
    path <- fixture(label)
    mutate(path)
    if (reseal) seal(path)
    output <- file.path(work, paste0(label, "-pack"))
    error <- tryCatch({duckvep_fastvep_publish(path, output); NULL}, error = function(error) error)
    stopifnot(inherits(error, "error"), !file.exists(output))
    rejected <<- c(rejected, label)
  }
  check("changed_failure", function(path) writeLines("changed", file.path(path, exact_failure)))
  check("changed_hgvs_body_failure", function(path) writeLines("changed", file.path(path, body_failure)))
  check("changed_hgvs_body_summary", function(path) writeLines("changed", file.path(path, body_summary)))
  check("missing_hgvs_body_failure", function(path) unlink(file.path(path, body_failure)))
  check("missing_hgvs_body_summary", function(path) unlink(file.path(path, body_summary)))
  check("changed_excluded_model", function(path) writeLines("changed", file.path(path, "forward/model.duckdb")))
  check("changed_input", function(path) {
    input <- read_table(file.path(path, "inputs_before.tsv"))$path[[1L]]
    writeLines("changed", input)
  })
  check("missing_manifest_row", function(path) {
    file <- file.path(path, "artifacts.tsv")
    value <- read_table(file)
    table(value[basename(value$path) != "field_failures.parquet", ], file)
  })
  check("missing_hgvs_body_manifest_row", function(path) {
    file <- file.path(path, "artifacts.tsv")
    value <- read_table(file)
    table(value[basename(value$path) != "hgvs_body_failures.parquet", ], file)
  })
  check("unmanifested_file", function(path) writeLines("extra", file.path(path, "extra.tsv")))
  check("duplicate_manifest_row", function(path) {
    file <- file.path(path, "artifacts.tsv")
    value <- read_table(file)
    table(rbind(value, value[1L, ]), file)
  })
  check("compressed_name_collision", function(path) {
    writeLines("collision", file.path(path, "forward/vep.vcf.gz"))
  }, reseal = TRUE)
  existing <- tryCatch({duckvep_fastvep_publish(original, second); NULL}, error = function(error) error)
  nested <- tryCatch({duckvep_fastvep_publish(original, file.path(original, "nested")); NULL}, error = function(error) error)
  stopifnot(inherits(existing, "error"), inherits(nested, "error"),
    !file.exists(file.path(original, "nested")), identical(validate(relocated), validate(second)))
  cat("Portable field pack: deterministic gzip, relocation, full records and failure bytes retained;",
    length(rejected), "corruption controls and two destination guards rejected\n")
}

main()
