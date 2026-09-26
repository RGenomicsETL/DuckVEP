#!/usr/bin/env Rscript
cli <- Sys.getenv("DUCKDB_CLI", unset = Sys.which("duckdb"))
if (!nzchar(cli)) stop("A DuckDB CLI is required (set DUCKDB_CLI)")
Sys.setenv(DUCKDB_CLI = normalizePath(cli, mustWork = TRUE))
rmarkdown::render("README.Rmd", output_file = "README.md", quiet = TRUE,
                  envir = new.env(parent = globalenv()))
lines <- readLines("README.md", warn = FALSE)
lines <- sub("[[:blank:]]+$", "", lines)
lines <- sub("^``` (sql|sh|r)$", "```\\1", lines)
writeLines(lines, "README.md")
