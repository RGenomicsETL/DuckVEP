#!/usr/bin/env Rscript

run_file <- local({
  argument_files <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  )
  candidates <- argument_files[basename(argument_files) == "run.R"]
  if (length(candidates) == 0L) {
    stop("run this workflow with Rscript pipelines/duckvep/run.R")
  }
  normalizePath(candidates[[length(candidates)]], mustWork = TRUE)
})
source(file.path(dirname(run_file), "R", "root.R"))
root <- duckvep_targets_root(run_file)
script <- file.path(root, "pipelines", "duckvep", "_targets.R")
store <- Sys.getenv(
  "DUCKVEP_TARGETS_STORE",
  unset = file.path(root, ".tmp", "duckvep_targets", "store")
)

targets::tar_make(script = script, store = store)
errored <- targets::tar_errored(store = store)
if (length(errored) > 0L) {
  stop(
    "DuckVEP targets failed: ",
    paste(errored, collapse = ", "),
    call. = FALSE
  )
}
