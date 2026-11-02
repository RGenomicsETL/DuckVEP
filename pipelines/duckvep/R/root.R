duckvep_targets_root <- function(script_file) {
  configured <- c(Sys.getenv("DUCKVEP_ROOT"), Sys.getenv("DUCKHTS_ROOT"))
  configured <- configured[nzchar(configured)]
  candidate <- if (length(configured) > 0L) {
    configured[[1L]]
  } else {
    file.path(dirname(script_file), "..", "..")
  }
  if (!dir.exists(candidate)) {
    stop("DuckVEP root does not exist: ", candidate)
  }
  root <- normalizePath(candidate, mustWork = TRUE)
  required <- c(
    "Makefile",
    "scripts/duckvep_evidence.R",
    "pipelines/duckvep/_targets.R",
    "pipelines/duckvep/R/pipeline.R"
  )
  missing <- required[!file.exists(file.path(root, required))]
  if (length(missing) > 0L) {
    stop(
      "DuckVEP root is missing ", paste(missing, collapse = ", "),
      "; set DUCKVEP_ROOT to a complete DuckVEP checkout"
    )
  }
  root
}
