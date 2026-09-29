#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("r/Rduckvep/R/expansionhunter.R")
vcf <- readLines("test/duckvep/conformance/data/expansionhunter_v5_example.vcf")
record <- strsplit(tail(vcf, 1L), "\t", fixed = TRUE)[[1L]]
stopifnot(identical(record[1L], "chr1"), identical(record[2L], "2005"))
fasta <- paste(readLines("test/duckvep/conformance/data/expansionhunter_v5_reference.fa")[-1L],
               collapse = "")
stopifnot(identical(substr(fasta, 2005L, 2005L), record[4L]))
info <- record[8L]
stopifnot(grepl("END=2008", info, fixed = TRUE))
reference <- substr(fasta, 2006L, 2008L)
prepare <- function(info, i) rduckvep_prepare_expansionhunter(info, record[9L],
  record[10L], record[4L], record[5L], reference, alt_index = i)
# The upstream fixture reports RL=1 despite a three-base literal interval.
stopifnot(identical(reference, "CAG"),
          identical(prepare(info, 1L)$reason, "reference_mismatch"),
          identical(prepare(info, 2L)$reason, "reference_mismatch"))
# Controlled correction of that field exercises the exact reconstruction path.
corrected <- sub("RL=1;", "RL=3;", info, fixed = TRUE)
stopifnot(!identical(corrected, info))
local({
directory <- tempfile("duckvep-receipt-")
dir.create(directory)
binary <- file.path(directory, "duckvep.duckdb_extension")
stopifnot(file.copy("build/release/duckvep.duckdb_extension", binary))
on.exit(unlink(directory, recursive = TRUE), add = TRUE)
con <- dbConnect(duckdb(shared_home = FALSE,
                        config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
dbExecute(con, paste("LOAD", as.character(dbQuoteString(con, binary))))
for (i in 1:2) {
  prepared <- prepare(corrected, i)
  stopifnot(identical(prepared$status, "ok"), prepared$sequence_exact)
  ref <- prepared$reference_components
  alt <- prepared$alternate_components
  query <- sprintf("SELECT r.* FROM (SELECT duckvep_repeat_alleles(
    [{unit:%s,count:%s}],[{unit:%s,count:%s}],true) r)",
    as.character(dbQuoteString(con, ref$unit)), ref$count,
    as.character(dbQuoteString(con, alt$unit)), alt$count)
  observed <- dbGetQuery(con, query)
  # The independent literal oracle uses the upstream FASTA and STRn definition.
  expected <- rep(reference, 1L)
  alternate <- strrep("CAG", c(2L, 10L)[i])
  stopifnot(identical(observed$status, "ok"),
            identical(observed$reference, expected),
            identical(observed$alternate, alternate),
            observed$length_change == nchar(alternate) - nchar(reference))
}
})
cat("ExpansionHunter v5 example: 2 raw reference mismatches, 2 controlled exact alleles\n")
