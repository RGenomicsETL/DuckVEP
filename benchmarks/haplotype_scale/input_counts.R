#!/usr/bin/env Rscript
# Physical records and ALT alleles of bgzipped single-sample VCF inputs (untimed):
#   Rscript benchmarks/haplotype_scale/input_counts.R FILE.vcf.gz ... > input_counts.tsv
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
con <- dbConnect(duckdb(config = list(threads = "4")))
cat("input\trecords\talt_alleles\n")
for (path in commandArgs(TRUE)) {
  n <- as.integer(system2("sh", c("-c", shQuote(paste0("zcat ", shQuote(path), " | head -n 5000 | grep -c '^#'"))), stdout = TRUE))
  r <- dbGetQuery(con, sprintf("SELECT count(*) AS records, coalesce(sum(len(string_split(alt, ','))) FILTER (WHERE alt <> '.'), 0) AS alts
    FROM read_csv(%s, delim='\\t', header=false, skip=%d, auto_detect=false, quote='', escape='', strict_mode=false,
    columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR',
    'info':'VARCHAR','fmt':'VARCHAR','sample':'VARCHAR'}, compression='gzip')", dbQuoteString(con, path), n))
  cat(basename(path), "\t", format(r$records, scientific = FALSE), "\t", format(r$alts, scientific = FALSE), "\n", sep = "")
}
