#!/usr/bin/env Rscript
source("benchmarks/benchmark_duckvep_fastvep_extract.R")
source("benchmarks/benchmark_duckvep_fastvep_fields.R")

main <- function() {
  directory <- tempfile("fastvep-extract-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- file.path(directory, "input.vcf")
  header <- c("##fileformat=VCFv4.2",
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature|Codons|HGVSp">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tX=3;CSQ=C|tx1|-|p.%3D,C|tx2||",
    "1\t4\te2\tG\tT\t.\tPASS\tCSQ=T|tx1|A/B|p.(a&b);X=4"), input)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fields <- c("Uploaded_variation", "Allele", "Feature", "Codons", "HGVSp")
  duckvep_fastvep_extract_csq(con, input, "actual", fields)
  actual <- DBI::dbGetQuery(con, "SELECT * FROM actual ORDER BY Uploaded_variation, Feature")
  expected <- data.frame(Uploaded_variation = c("e1", "e1", "e2"), Allele = c("C", "C", "T"),
    Feature = c("tx1", "tx2", "tx1"), Codons = c("-", "", "A/B"),
    HGVSp = c("p.%3D", "", "p.(a&b)"))
  stopifnot(identical(actual, expected))
  expect_error <- function(expression) stopifnot(inherits(tryCatch({
    force(expression)
    NULL
  }, error = identity), "error"))
  expect_error(duckvep_fastvep_extract_csq(con, input, "unsupported", c(fields, "GENCODE_PRIMARY")))
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tCSQ=C|tx1|-"), input)
  expect_error(duckvep_fastvep_extract_csq(con, input, "short", fields))
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tX=1"), input)
  expect_error(duckvep_fastvep_extract_csq(con, input, "missing", fields))
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tCSQ=C|tx1||;CSQ=C|tx2||"), input)
  expect_error(duckvep_fastvep_extract_csq(con, input, "duplicate", fields))
  tab <- file.path(directory, "output.tsv")
  native <- duckvep_fastvep_fields("native_tab17")
  writeLines(c("## fastvep version=0.3.0", paste0("#", paste(native, collapse = "\t"))), tab)
  stopifnot(duckvep_fastvep_tab_header(tab, native) == 1L)
  writeLines(paste(native, collapse = "\t"), tab)
  stopifnot(duckvep_fastvep_tab_header(tab, native) == 0L)
  renamed <- native
  renamed[native == "Gene"] <- "NotGene"
  duplicate <- native
  duplicate[native == "FLAGS"] <- "Gene"
  for (invalid in list(renamed, native[-length(native)], c(native, "Extra"), duplicate, rev(native))) {
    writeLines(paste(invalid, collapse = "\t"), tab)
    expect_error(duckvep_fastvep_tab_header(tab, native))
  }
  writeLines(character(), tab)
  expect_error(duckvep_fastvep_tab_header(tab, native))
  cat("CSQ extraction: transport spelling, missing fields and width controls passed\n")
}

main()
