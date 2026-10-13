# Focused R process boundary for the compiled opt-in product probe.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
source_dir <- file.path(root, "src/kernel/src")
bin <- tempfile("duckvep-products-")
on.exit(unlink(bin), add = TRUE)
compiler <- Sys.getenv("CC", "cc")
`%||%` <- function(x, y) if (is.null(x)) y else x
files <- file.path(root, c("test/duckvep/conformance/products_native_probe.c",
                           "src/kernel/src/duckvep_product_execute.c",
                           "src/kernel/src/duckvep_product_profile.c",
                           "src/kernel/src/duckvep_haplotype.c",
                           "src/kernel/src/duckvep_codon.c"))
output <- system2(compiler, c("-std=gnu11", paste0("-I", source_dir),
                              paste0("-I", file.path(root, "src/kernel/include")), files,
                              "-o", bin), stdout = TRUE, stderr = TRUE)
if (!identical(attr(output, "status") %||% 0L, 0L)) stop(paste(output, collapse = "\n"))

input <- list(dna = "ATGAAATAAGGGCCCGGGATGAAATAA", ref = "TAA", alt = "CAA")
before <- input
run <- function(mapping, rna, ref = input$ref) {
  result <- system2(bin, c(input$dna, mapping, rna, ref, input$alt),
                    stdout = TRUE, stderr = TRUE)
  stopifnot(identical(attr(result, "status") %||% 0L, 0L), length(result) == 1L)
  result
}
skip <- run("skip", "ATGAAATAAATGAAATAA")
include <- run("include", "ATGAAATAACCCATGAAATAA")
stopifnot(grepl("status=0 product=skip ", skip, fixed = TRUE),
          grepl("status=0 product=include ", include, fixed = TRUE),
          grepl("alt_rna=ATGAAACAAATGAAATAA", skip, fixed = TRUE),
          grepl("alt_rna=ATGAAACAACCCATGAAATAA", include, fixed = TRUE),
          grepl("ref_peptide=MK*MK* alt_peptide=MKQMK*", skip, fixed = TRUE),
          grepl("ref_peptide=MK*PMK* alt_peptide=MKQPMK*", include, fixed = TRUE),
          grepl("source=11/rna7", skip, fixed = TRUE),
          grepl("source=11/rna7", include, fixed = TRUE),
          grepl("status=2", run("skip", "ATGAAAGAAATGAAATAA"), fixed = TRUE),
          grepl("status=2", run("skip", "ATGAAATAAATGAAATAA", "AAA"), fixed = TRUE),
          identical(input, before))
message("Compiled product probe: paired RNA/peptide and independent REF boundaries passed")
