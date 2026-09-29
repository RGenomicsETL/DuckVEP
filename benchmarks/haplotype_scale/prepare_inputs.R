#!/usr/bin/env Rscript
# Builds the untimed inputs of the #2 slice 7 scale qualification from the staged HG002 VCF and gnomAD v4.1 exomes.
# Nothing here is measured and nothing it writes is committed (inputs are large; their hashes are in the README).
#
#   Rscript benchmarks/haplotype_scale/prepare_inputs.R WORKDIR
#
# Writes into WORKDIR:
#   hg002_calls.parquet, hg002_calls_inside.parquet, hg002_calls_outside.parquet
#       the ordered calls relation of the identical HG002 input on the full Ensembl 116 model (mode A input) and its
#       partition by the committed csq accounting: records whose category is 'compared' (23,866) versus the rest
#   qual5m.vcf.gz, qual5m_calls.parquet
#       the 5,000,000-record single-sample job (HG002 4,023,088 records plus 976,912 seeded gnomAD exome records) and its
#       calls relation on the MANE Select model
#   dense.vcf.gz / low_sharing.vcf.gz and their calls parquet: stress controls (see the README)
#
# Top-up rule (deterministic): candidate sites are the gnomAD v4.1 exome ALT rows with status 'literal' and FILTER PASS on
# chromosomes 1-22, X and Y whose (contig, position) is absent from HG002; one ALT per site is kept by the smallest
# hash(seed, contig, position, ref, alt); the first N sites by that hash are taken, N = 5,000,000 - 4,023,088. The seed is
# 20260929. Genotype is a function of the same hash: 0|1 for 45%, 1|0 for 45%, 1|1 for 10%; PS is PATMAT
# for heterozygotes and HOMVAR for homozygotes: the HG002 file's string PS labels, which are no phase set, so every phased
# record of the sample is one phase domain per transcript and the top-up records compose with HG002's own (a numeric PS would put
# them in a different domain from HG002's and make most transcripts unresolved_cross_ps_phase).
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript prepare_inputs.R WORKDIR")
work <- normalizePath(args[[1L]], mustWork = TRUE)
here <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
scale_dir <- dirname(here)
root <- normalizePath(file.path(scale_dir, "..", ".."))
extension <- Sys.getenv("DUCKVEP_SCALE_EXTENSION", file.path(root, "build", "release", "duckvep.duckdb_extension"))
data <- Sys.getenv("DUCKVEP_DATA", "/root/duckvep/data")
hg002 <- file.path(data, "hg002-csq", "hg002.ens.vcf.gz")
ledger <- Sys.getenv("DUCKVEP_HG002_LEDGER", file.path(work, "hg002_domain.counts.tsv.records.tsv.gz"))
seed <- 20260929L
target <- 5000000L
stopifnot(file.exists(hg002))
worker <- function(...) {
  status <- system2("Rscript", c(file.path(scale_dir, "worker.R"), "--extension", shQuote(extension), ...))
  if (status != 0L) stop("worker failed")
}
con <- dbConnect(duckdb(shared_home = FALSE, config = list(threads = "4", memory_limit = "16GB")))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
sql <- function(s) invisible(dbExecute(con, s))
get <- function(s) dbGetQuery(con, s)
q <- function(x) as.character(dbQuoteString(con, x))
header_n <- as.integer(system2("sh", c("-c", shQuote(paste0("zcat ", shQuote(hg002), " | head -n 5000 | grep -c '^#'"))), stdout = TRUE))
vcf_cols <- "columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR','fmt':'VARCHAR','sample':'VARCHAR'}"
read_vcf <- function(path, skip) sprintf("read_csv(%s, delim='\\t', header=false, skip=%d, auto_detect=false, quote='', escape='', strict_mode=false, %s, compression='gzip')",
  q(path), skip, vcf_cols)
header_file <- file.path(work, "header.txt")
system2("sh", c("-c", shQuote(paste0("zcat ", shQuote(hg002), " | head -n ", header_n, " > ", shQuote(header_file)))))
contigs <- system2("sh", c("-c", shQuote(paste0("grep '^##contig' ", shQuote(header_file), " | sed 's/.*ID=\\([^,]*\\).*/\\1/'"))), stdout = TRUE)
write_vcf <- function(query, path) {
  body <- tempfile("body", tmpdir = work, fileext = ".tsv")
  sql(sprintf("COPY (%s) TO %s (FORMAT csv, DELIMITER '\\t', HEADER false, QUOTE '', ESCAPE '')", query, q(body)))
  status <- system2("sh", c("-c", shQuote(sprintf("cat %s %s | bgzip -c > %s && tabix -p vcf %s && rm %s", shQuote(header_file), shQuote(body), shQuote(path), shQuote(path), shQuote(body)))))
  stopifnot(status == 0L)
}
contig_rank <- function(column) sprintf("list_position(%s, %s) - 1", paste0("[", paste(q(contigs), collapse = ","), "]"), column)

built <- function(name) file.exists(file.path(work, name))   # steps whose output exists are skipped, so an interrupted run resumes
if (!built("qual5m.vcf.gz")) {
# gnomAD exome candidates. The model's contig names have no chr prefix.
message("gnomAD exome candidates ...")
sql(sprintf("CREATE TABLE hg002_sites AS SELECT DISTINCT chrom, pos FROM %s", read_vcf(hg002, header_n)))
sql(sprintf("CREATE TABLE candidates AS SELECT replace(chrom, 'chr', '') AS chrom, position AS pos, reference AS ref, alternate AS alt,
  hash(%d::VARCHAR || ':' || chrom || ':' || position || ':' || reference || ':' || alternate) AS h
  FROM read_parquet(%s) WHERE status = 'literal' AND filter = 'PASS' AND replace(chrom, 'chr', '') IN (SELECT unnest(%s))",
  seed, q(file.path(data, "gnomad-v4.1", "exomes-*", "part-*.parquet")), paste0("[", paste(q(contigs[!grepl("^(MT|GL|KI)", contigs)]), collapse = ","), "]")))
sql("CREATE TABLE sites AS SELECT * FROM (SELECT c.*, row_number() OVER (PARTITION BY c.chrom, c.pos ORDER BY c.h) AS k
  FROM candidates c ANTI JOIN hg002_sites s ON s.chrom = c.chrom AND s.pos = c.pos) WHERE k = 1")
n_hg002 <- get(sprintf("SELECT count(*) AS n FROM %s", read_vcf(hg002, header_n)))$n
n_topup <- target - n_hg002
message(sprintf("HG002 records %d, top-up %d", n_hg002, n_topup))
sql(sprintf("CREATE TABLE topup AS SELECT chrom, pos, '.' AS id, ref, alt, '.' AS qual, 'PASS' AS filter, '.' AS info, 'GT:PS' AS fmt,
  CASE WHEN h %% 100 < 45 THEN '0|1' WHEN h %% 100 < 90 THEN '1|0' ELSE '1|1' END || ':' ||
  CASE WHEN h %% 100 < 90 THEN 'PATMAT' ELSE 'HOMVAR' END AS sample
  FROM sites ORDER BY h LIMIT %d", n_topup))
merged <- sprintf("SELECT chrom, pos, id, ref, alt, qual, filter, info, fmt, sample FROM (
   SELECT chrom, pos, id, ref, alt, qual, filter, info, fmt, sample FROM %s UNION ALL
   SELECT chrom, pos, id, ref, alt, qual, filter, info, fmt, sample FROM topup) ORDER BY %s, pos, ref, alt", read_vcf(hg002, header_n), contig_rank("chrom"))
write_vcf(merged, file.path(work, "qual5m.vcf.gz"))
message(get(sprintf("SELECT count(*) AS n FROM %s", read_vcf(file.path(work, "qual5m.vcf.gz"), header_n)))$n, " records in qual5m.vcf.gz")
}

# HG002 calls on the full model, then the csq-domain partition from the committed accounting's regenerated ledger.
if (!built("hg002_calls.parquet"))
  worker("--model", "full", "--mode", "stage", "--vcf", shQuote(hg002), "--calls", shQuote(file.path(work, "hg002_calls.parquet")), "--label", "stage-hg002")
if (file.exists(ledger)) {
  sql(sprintf("CREATE TABLE ledger AS SELECT record_index, category FROM read_csv(%s, delim='\\t', header=true, columns={'record_index':'BIGINT','chrom':'VARCHAR','pos':'BIGINT','ref':'VARCHAR','alt':'VARCHAR','category':'VARCHAR','reason':'VARCHAR'})", q(ledger)))
  for (part in c("inside", "outside")) {
    op <- if (part == "inside") "=" else "<>"
    sql(sprintf("COPY (SELECT c.* FROM read_parquet(%s) c JOIN ledger l ON l.record_index = (c.event_index >> 6) WHERE l.category %s 'compared' ORDER BY seq_region, position, event_index, transcript_index)
      TO %s (FORMAT parquet)", q(file.path(work, "hg002_calls.parquet")), op, q(file.path(work, sprintf("hg002_calls_%s.parquet", part)))))
  }
} else message("no csq ledger at ", ledger, ": the inside/outside partition is skipped")

# MANE Select calls for the 5M job.
if (!built("qual5m_calls.parquet"))
  worker("--model", "mane", "--mode", "stage", "--vcf", shQuote(file.path(work, "qual5m.vcf.gz")), "--calls", shQuote(file.path(work, "qual5m_calls.parquet")), "--label", "stage-5m")

# Controls on the MANE model, built from CDS-exon bases only (the coding part of every exon, so every record reaches the classifier).
# Dense/long-transcript: gnomAD exome variants with REF and ALT of at most 3 bases inside the CDS of the MANE transcript with the
# longest coding sequence, at most one per 16-base block and only in the first 8 bases of a block (so no two edits overlap), phased
# in one phase set: thousands of edits on one haplotype of a 100 kb coding sequence, indels included. Low-sharing: one variant
# per MANE transcript, so no two calls share a transcript and nothing shares an edit path.
sql("ATTACH '/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb' AS m (READ_ONLY)")
sql("CREATE TABLE mane_exons AS SELECT t.transcript_stable_id, t.seq_region_name AS chrom, greatest(e.exon_start, t.cds_start) AS s, least(e.exon_end, t.cds_end) AS en,
  octet_length(t.cds_sequence) AS cds_bytes FROM m.model_transcripts t, unnest(t.exons) u(e)
  WHERE t.mane_select_refseq IS NOT NULL AND t.cds_sequence IS NOT NULL AND greatest(e.exon_start, t.cds_start) <= least(e.exon_end, t.cds_end)")
sql("CREATE TABLE longest AS SELECT transcript_stable_id, chrom, cds_bytes FROM mane_exons GROUP BY ALL ORDER BY cds_bytes DESC, transcript_stable_id LIMIT 1")
print(get("SELECT * FROM longest"))
exomes <- q(file.path(data, "gnomad-v4.1", "exomes-*", "part-*.parquet"))
dense_n <- 4000L
sql(sprintf("CREATE TABLE dense_sites AS SELECT * FROM (SELECT * FROM (SELECT replace(g.chrom, 'chr', '') AS chrom, g.position AS pos, g.reference AS ref, g.alternate AS alt,
      hash(%d::VARCHAR || ':' || g.chrom || ':' || g.position || ':' || g.reference || ':' || g.alternate) AS h
    FROM read_parquet(%s) g JOIN longest l ON l.chrom = replace(g.chrom, 'chr', '')
    JOIN mane_exons x ON x.transcript_stable_id = l.transcript_stable_id AND g.position >= x.s AND g.position + length(g.reference) - 1 <= x.en
    WHERE g.status = 'literal' AND g.filter = 'PASS' AND length(g.reference) <= 3 AND length(g.alternate) <= 3 AND g.position %% 16 < 8)
  QUALIFY row_number() OVER (PARTITION BY pos // 16 ORDER BY h) = 1) ORDER BY h LIMIT %d", seed, exomes, dense_n))
sql("CREATE TABLE dense AS SELECT chrom, pos, '.' AS id, ref, alt, '.' AS qual, 'PASS' AS filter, '.' AS info, 'GT:PS' AS fmt,
  CASE WHEN h % 2 = 0 THEN '0|1' ELSE '1|0' END || ':1' AS sample FROM dense_sites")
print(get("SELECT count(*) AS dense_records, count(*) FILTER (WHERE length(ref) <> length(alt)) AS indels FROM dense"))
write_vcf(sprintf("SELECT * FROM dense ORDER BY %s, pos, ref, alt", contig_rank("chrom")), file.path(work, "dense.vcf.gz"))
sql(sprintf("CREATE TABLE low_sites AS SELECT * FROM (SELECT replace(g.chrom, 'chr', '') AS chrom, g.position AS pos, g.reference AS ref, g.alternate AS alt,
    hash(%d::VARCHAR || ':' || g.chrom || ':' || g.position || ':' || g.reference || ':' || g.alternate) AS h, x.transcript_stable_id
  FROM read_parquet(%s) g JOIN mane_exons x ON x.chrom = replace(g.chrom, 'chr', '') AND g.position >= x.s AND g.position + length(g.reference) - 1 <= x.en
  WHERE g.status = 'literal' AND g.filter = 'PASS' AND length(g.reference) <= 50 AND length(g.alternate) <= 50)
  QUALIFY row_number() OVER (PARTITION BY transcript_stable_id ORDER BY h) = 1", seed, exomes))
sql("CREATE TABLE low AS SELECT chrom, pos, '.' AS id, ref, alt, '.' AS qual, 'PASS' AS filter, '.' AS info, 'GT:PS' AS fmt,
  CASE WHEN h % 2 = 0 THEN '0|1' ELSE '1|0' END || ':1' AS sample FROM (SELECT * FROM low_sites QUALIFY row_number() OVER (PARTITION BY chrom, pos ORDER BY h) = 1)")
print(get("SELECT count(*) AS low_sharing_records FROM low"))
write_vcf(sprintf("SELECT * FROM low ORDER BY %s, pos, ref, alt", contig_rank("chrom")), file.path(work, "low_sharing.vcf.gz"))
for (control in c("dense", "low_sharing")) {
  file.remove(file.path(work, paste0(control, "_calls.parquet"))[file.exists(file.path(work, paste0(control, "_calls.parquet")))])
  worker("--model", "mane", "--mode", "stage", "--vcf", shQuote(file.path(work, paste0(control, ".vcf.gz"))),
    "--calls", shQuote(file.path(work, paste0(control, "_calls.parquet"))), "--label", paste0("stage-", control))
}
message("done")
