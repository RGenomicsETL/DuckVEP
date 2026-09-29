#!/usr/bin/env Rscript
source("r/Rduckvep/R/expansionhunter.R")
source("r/Rduckvep/R/structural_geometry.R")
row <- strsplit(tail(readLines(
  "test/duckvep/conformance/data/expansionhunter_v5_documented_als.vcf"), 1L),
  "\t", fixed = TRUE)[[1L]]
for (i in seq_len(100L)) {
  rduckvep_prepare_expansionhunter(row[8L], row[9L], row[10L], row[4L],
    row[5L], strrep("GGCCCC", 3L), alt_index = 1L)
  rduckvep_prepare_sv_geometry(13546123, "G", "<INS>",
    "END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG")
}
profile <- function(expr, n = 10000L) {
  gc(reset = TRUE)
  for (i in seq_len(n)) eval(expr)
  peak <- gc()
  c(peak_Ncells_MB = peak[1L, 6L],
    peak_Vcells_MB = peak[2L, 6L])
}
for (name in c("str_exact", "str_summary", "sv_confidence", "sv_literal")) {
  expr <- switch(name,
    str_exact = quote(rduckvep_prepare_expansionhunter(row[8L], row[9L],
      row[10L], row[4L], row[5L], strrep("GGCCCC", 3L), alt_index = 1L)),
    str_summary = quote(rduckvep_prepare_expansionhunter(row[8L], row[9L],
      row[10L], row[4L], row[5L], strrep("GGCCCC", 3L), alt_index = 2L)),
    sv_confidence = quote(rduckvep_prepare_sv_geometry(13546123, "G", "<INS>",
      "END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG")),
    sv_literal = quote(rduckvep_prepare_sv_geometry(13546123, "G", "GATG", ".")))
  print(c(path = name, profile(expr)))
}
results <- bench::mark(
  str_exact = rduckvep_prepare_expansionhunter(row[8L], row[9L], row[10L],
    row[4L], row[5L], strrep("GGCCCC", 3L), alt_index = 1L),
  str_summary = rduckvep_prepare_expansionhunter(row[8L], row[9L], row[10L],
    row[4L], row[5L], strrep("GGCCCC", 3L), alt_index = 2L),
  sv_confidence = rduckvep_prepare_sv_geometry(13546123, "G", "<INS>",
    "END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG"),
  sv_literal = rduckvep_prepare_sv_geometry(13546123, "G", "GATG", "."),
  iterations = 2000L, check = FALSE, memory = TRUE)
print(results[, c("expression", "min", "median", "n_itr", "n_gc")])
