#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "readme-evidence.R"))
source(file.path(root, "pipelines", "duckvep", "R", "root.R"))

restore_environment <- function(names, values) {
  for (i in seq_along(names)) {
    if (is.na(values[[i]])) {
      Sys.unsetenv(names[[i]])
    } else {
      do.call(Sys.setenv, as.list(setNames(values[[i]], names[[i]])))
    }
  }
}

with_environment <- function(names, code) {
  previous <- Sys.getenv(names, unset = NA_character_)
  on.exit(restore_environment(names, previous), add = TRUE)
  force(code)
}

path_function <- function(path, name) {
  expressions <- parse(path)
  definitions <- vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(as.character(expression[[2L]]), name)
  }, logical(1))
  stopifnot(sum(definitions) == 1L)
  environment <- new.env(parent = baseenv())
  eval(expressions[[which(definitions)]], environment)
  get(name, envir = environment, inherits = FALSE)
}

map <- c(old_reachable = "new_reachable", old_unreachable = "new_unreachable")
ancestry <- c("head", "new_reachable", "older")
stopifnot(identical(
  duckvep_newest_tested(c("old_reachable", "old_unreachable"), ancestry, map),
  c(TRUE, FALSE)
))
stopifnot(identical(
  duckvep_newest_tested(c("unmapped", "old_unreachable"), ancestry, map),
  c(FALSE, FALSE)
))

with_environment(c("DUCKVEP_ROOT", "DUCKHTS_ROOT"), {
  extracted <- tempfile("duckvep-extracted-")
  script <- file.path(extracted, "pipelines", "duckvep", "_targets.R")
  dir.create(dirname(script), recursive = TRUE)
  dir.create(file.path(extracted, "pipelines", "duckvep", "R"), recursive = TRUE)
  dir.create(file.path(extracted, "scripts"), recursive = TRUE)
  file.create(script)
  file.create(file.path(extracted, "Makefile"))
  file.create(file.path(extracted, "pipelines", "duckvep", "R", "pipeline.R"))
  file.create(file.path(extracted, "scripts", "duckvep_evidence.R"))

  Sys.unsetenv(c("DUCKVEP_ROOT", "DUCKHTS_ROOT"))
  stopifnot(identical(duckvep_targets_root(script), normalizePath(extracted)))
  Sys.setenv(DUCKVEP_ROOT = extracted)
  stopifnot(identical(duckvep_targets_root(script), normalizePath(extracted)))
  Sys.unsetenv("DUCKVEP_ROOT")
  Sys.setenv(DUCKHTS_ROOT = extracted)
  stopifnot(identical(duckvep_targets_root(script), normalizePath(extracted)))

  manifest_root <- tempfile("duckvep-manifest-")
  dir.create(file.path(manifest_root, "pipelines"), recursive = TRUE)
  dir.create(file.path(manifest_root, "scripts"), recursive = TRUE)
  file.create(file.path(manifest_root, "Makefile"))
  stopifnot(file.copy(
    file.path(root, "pipelines", "duckvep"),
    file.path(manifest_root, "pipelines"),
    recursive = TRUE
  ))
  stopifnot(file.copy(
    file.path(root, "scripts", "duckvep_evidence.R"),
    file.path(manifest_root, "scripts")
  ))
  Sys.unsetenv(c("DUCKVEP_ROOT", "DUCKHTS_ROOT"))
  manifest_script <- file.path(manifest_root, "pipelines", "duckvep", "_targets.R")
  expression <- paste0(
    "setwd(", encodeString(tempdir(), quote = '\"'), "); ",
    "manifest <- targets::tar_manifest(script = ",
    encodeString(manifest_script, quote = '\"'), "); ",
    "stopifnot(nrow(manifest) > 0L)"
  )
  manifest <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(expression)),
    stdout = TRUE,
    stderr = TRUE
  )
  stopifnot(is.null(attr(manifest, "status")))

  incomplete <- tempfile("duckvep-incomplete-")
  dir.create(incomplete)
  Sys.setenv(DUCKVEP_ROOT = incomplete)
  message <- tryCatch(
    duckvep_targets_root(script),
    error = function(error) conditionMessage(error)
  )
  stopifnot(grepl("complete DuckVEP checkout", message, fixed = TRUE))
})

with_environment(c("DUCKVEP_DATA", "PR66_TEST_PATH"), {
  for (file in c(
    "breakend_fusion_corpus_differential.R",
    "expansionhunter_v5_differential.R",
    "structural_hgvs_vep116_differential.R",
    "sv_payload_vep116_differential.R"
  )) {
    resolver <- path_function(
      file.path(root, "test", "duckvep", "conformance", file),
      "duckvep_data"
    )
    Sys.unsetenv(c("DUCKVEP_DATA", "PR66_TEST_PATH"))
    message <- tryCatch(
      resolver("models/model.duckdb", "PR66_TEST_PATH"),
      error = function(error) conditionMessage(error)
    )
    stopifnot(grepl("PR66_TEST_PATH or DUCKVEP_DATA", message, fixed = TRUE))
    data_root <- tempfile("duckvep-data-")
    dir.create(data_root)
    Sys.setenv(DUCKVEP_DATA = data_root)
    stopifnot(identical(
      resolver("models/model.duckdb", "PR66_TEST_PATH"),
      file.path(data_root, "models", "model.duckdb")
    ))
    Sys.setenv(PR66_TEST_PATH = "/configured/model.duckdb")
    stopifnot(identical(
      resolver("models/model.duckdb", "PR66_TEST_PATH"),
      "/configured/model.duckdb"
    ))
  }
})

with_environment(c("DUCKHTS_CACHE_DIR", "PR66_TEST_CACHE"), {
  for (file in c(
    "structural_hgvs_vep116_differential.R",
    "sv_payload_vep116_differential.R"
  )) {
    resolver <- path_function(
      file.path(root, "test", "duckvep", "conformance", file),
      "duckhts_cache"
    )
    Sys.unsetenv(c("DUCKHTS_CACHE_DIR", "PR66_TEST_CACHE"))
    message <- tryCatch(
      resolver("cache/leaf", "PR66_TEST_CACHE"),
      error = function(error) conditionMessage(error)
    )
    stopifnot(grepl("PR66_TEST_CACHE or DUCKHTS_CACHE_DIR", message, fixed = TRUE))
    cache_root <- tempfile("duckhts-cache-")
    dir.create(cache_root)
    Sys.setenv(DUCKHTS_CACHE_DIR = cache_root)
    stopifnot(identical(
      resolver("cache/leaf", "PR66_TEST_CACHE"),
      file.path(cache_root, "cache", "leaf")
    ))
    Sys.setenv(PR66_TEST_CACHE = "/configured/cache")
    stopifnot(identical(
      resolver("cache/leaf", "PR66_TEST_CACHE"),
      "/configured/cache"
    ))
  }
})

for (file in c(
  "ambiguous_codon_differential.R",
  "corpus_differential.R",
  "haplotype_differential.R",
  "haplotype_grouped_flags.R",
  "haplotype_model_differential.R",
  "haplotype_observer_contract.R",
  "haplotype_phase_differential.R",
  "haplotype_record_differential.R",
  "haplotype_reference_differential.R",
  "haplotype_sql_differential.R",
  "peptide_edit_differential.R",
  "reference_translation_differential.R",
  "sequence_diff_differential.R"
)) {
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    file.path(root, "test", "duckvep", "conformance", file),
    stdout = TRUE,
    stderr = TRUE,
    env = "VEP_PREFIX="
  ))
  status <- attr(output, "status")
  stopifnot(!is.null(status), status != 0L)
  stopifnot(any(grepl("VEP_PREFIX.*required", output)))
}

with_environment("DUCKVEP_CORPUS_STAGER", {
  wrapper <- file.path(root, "scripts", "stage_duckvep_conformance_corpora.sh")
  stager <- tempfile(fileext = ".R")
  Sys.setenv(DUCKVEP_CORPUS_STAGER = stager)
  output <- suppressWarnings(system2("bash", shQuote(wrapper), stdout = TRUE, stderr = TRUE))
  stopifnot(identical(attr(output, "status"), 2L),
    any(grepl("Set DUCKVEP_CORPUS_STAGER", output, fixed = TRUE)))
  writeLines("cat(paste(commandArgs(trailingOnly = TRUE), collapse = '\\n'))", stager)
  output <- system2("bash", c(shQuote(wrapper), shQuote("output path"), "dbvar-chr22"),
    stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(output, "status")),
    identical(output, c("--corpus", "dbvar-chr22", "--output-dir", "output path")))
  unlink(stager)
})

evidence <- readme_evidence(root)
stopifnot(identical(evidence$newest_conformance_run, "2026-07-16"))
stopifnot(identical(evidence$newest_hgvs_run, "2026-09-07"))

# A corpus with no reachable revision is omitted from the complete evidence table.
probe <- readme_evidence
environment(probe) <- list2env(list(read.csv = function(file, ...) {
  rows <- utils::read.csv(file, ...)
  if (basename(file) == "conformance_history.csv") {
    rows$source_revision[rows$corpus == "final_dbsnp157_windows_seed29"] <- "unmapped"
  }
  rows
}), parent = environment(readme_evidence))
partial <- probe(root)
stopifnot(nrow(partial$conformance) == nrow(evidence$conformance) - 1L,
  !"dbSNP 157 windows" %in% partial$conformance$label)
