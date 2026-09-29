#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("r/Rduckvep/R/builders.R")
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
directory <- tempfile("duckvep-receipt-")
dir.create(directory)
binary <- file.path(directory, "duckvep.duckdb_extension")
stopifnot(file.copy(Sys.getenv("DUCKVEP_EXTENSION_FILE", "build/release/duckvep.duckdb_extension"), binary))
con <- dbConnect(duckdb(shared_home = FALSE,
                        config = list(allow_unsigned_extensions = "true")))
dbExecute(con, paste("LOAD", as.character(dbQuoteString(con, binary))))
prepare <- function(info, i, row = record, literal = reference) {
  dbWriteTable(con, "eh_input", data.frame(event_index = 1L, info,
    format = row[9L], sample = row[10L], ref = row[4L], alt = row[5L],
    alt_index = i), temporary = TRUE, overwrite = TRUE)
  dbWriteTable(con, "eh_reference", data.frame(event_index = 1L,
    reference_sequence = literal), temporary = TRUE, overwrite = TRUE)
  as.list(rduckvep_prepare_expansionhunter(con, "eh_input", "eh_reference")[1L, ])
}
# The upstream fixture reports RL=1 despite a three-base literal interval.
stopifnot(identical(reference, "CAG"),
          identical(prepare(info, 1L)$reason, "reference_mismatch"),
          identical(prepare(info, 2L)$reason, "reference_mismatch"))
# Controlled correction of that field exercises the exact reconstruction path.
corrected <- sub("RL=1;", "RL=3;", info, fixed = TRUE)
stopifnot(!identical(corrected, info))
local({
on.exit({dbDisconnect(con, shutdown = TRUE); unlink(directory, recursive = TRUE)}, add = TRUE)
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
# The documented ALS example has a consistent, reference-backed exact allele.
als <- strsplit(tail(readLines(
  "test/duckvep/conformance/data/expansionhunter_v5_documented_als.vcf"), 1L),
  "\t", fixed = TRUE)[[1L]]
fasta37 <- Sys.getenv("DUCKVEP_GRCH37_FASTA", "/root/duckvep/data/mane-grch37/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
anchor <- tail(system2("samtools", c("faidx", fasta37,
  "9:27573526-27573526"), stdout = TRUE), 1L)
literal <- tail(system2("samtools", c("faidx", fasta37,
  "9:27573527-27573544"), stdout = TRUE), 1L)
stopifnot(identical(anchor, als[4L]), identical(literal, strrep("GGCCCC", 3L)))
call <- function(i) prepare(als[8L], i, als, literal)
exact <- call(1L)
stopifnot(identical(exact$status, "ok"),
  identical(call(2L)$reason, "estimated_count"))
ref <- exact$reference_components
alt <- exact$alternate_components
query <- sprintf("SELECT r.* FROM (SELECT duckvep_repeat_alleles(
  [{unit:%s,count:%s}],[{unit:%s,count:%s}],true) r)",
  as.character(dbQuoteString(con, ref$unit)), ref$count,
  as.character(dbQuoteString(con, alt$unit)), alt$count)
observed <- dbGetQuery(con, query)
stopifnot(identical(observed$status, "ok"),
          identical(observed$reference, literal),
          identical(observed$alternate, strrep("GGCCCC", 2L)),
          observed$length_change == -6)
})
cat("ExpansionHunter v5: 1 raw exact, 1 raw summary, 2 raw mismatches, 2 controlled exact alleles\n")
