#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript r/Rduckvep/bootstrap.R <repo>")
repo <- normalizePath(args[[1L]], mustWork = TRUE)
package <- file.path(repo, "r", "Rduckvep")
destination <- file.path(package, "inst", "duckvep_extension")
if (dir.exists(destination)) unlink(destination, recursive = TRUE)
tracked <- system2("git", c("-C", shQuote(repo), "ls-files", "--",
                            "src", "cmake", "duckdb_capi", "third_party/htslib",
                            "third_party/cgranges", "CMakeLists.txt"), stdout = TRUE)
if (!length(tracked)) stop("No extension sources found in ", repo)
for (path in tracked) {
  if (startsWith(path, "src/duckvep/")) stop("Sources must be at the repository root")
  target <- file.path(destination, path)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(repo, path), target, overwrite = TRUE))
    stop("Cannot copy ", path)
}
cat("Bundled", length(tracked), "tracked extension source files\n")
