#!/usr/bin/env Rscript
root <- normalizePath(Sys.getenv("DUCKVEP_REPO_ROOT", getwd()), mustWork = TRUE)
reports <- c(throughput = "duckvep_throughput.Rmd",
  vep_rs = "benchmark_duckvep_vep_rs.Rmd",
  haplotypes = "duckvep_haplotypes.Rmd",
  indels = "duckvep_haplotypes_indel.Rmd",
  projection = "benchmark_duckvep_projection.Rmd",
  conformance = "duckvep_conformance.Rmd")
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- names(reports)
if (!all(requested %in% names(reports))) {
  stop("Unknown report: ", paste(setdiff(requested, names(reports)), collapse = ", "))
}
for (report in requested) {
  rmarkdown::render(file.path(root, "benchmarks", reports[[report]]),
    knit_root_dir = root, output_options = list(html_preview = FALSE),
    envir = new.env(parent = globalenv()), quiet = TRUE)
}
