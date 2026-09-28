library(DBI)
library(duckdb)
root <- Sys.getenv("DUCKVEP_GNOMAD_ROOT", "/root/duckvep/data/gnomad-v4.1")
model <- Sys.getenv("DUCKVEP_SCALE_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
shards <- list.dirs(root, recursive = FALSE, full.names = TRUE)
shards <- shards[grepl("/exomes-chr(1[0-9]|2[0-2]|[1-9]|X|Y)$", shards)]
if (length(shards) == 0L || any(!file.exists(file.path(shards, "source.tsv")))) stop("no verified exome shards")
files <- unlist(lapply(sort(shards), function(dir) {
  parts <- list.files(dir, "^part-[0-9]+\\.parquet$", full.names = TRUE)
  if (length(parts) == 0L) stop("empty verified exome shard")
  sort(parts)
}), use.names = FALSE)
sha256 <- function(path) substr(system2("sha256sum", path, stdout = TRUE), 1L, 64L)
source_receipts <- data.frame(file = files, bytes = file.size(files),
  sha256 = vapply(files, sha256, character(1L)))
tmp <- tempfile("exome-sources-", tmpdir = root)
write.table(source_receipts, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
out <- file.path(root, "panels", paste0("exomes-v1-", substr(sha256(tmp), 1L, 16L)))
unlink(tmp)
if (dir.exists(out)) stop("exome panel artifact already exists")
dir.create(out, recursive = TRUE)
write.table(source_receipts, file.path(out, "sources.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
con <- dbConnect(duckdb(shared_home = FALSE))
q <- function(x) as.character(dbQuoteString(con, x))
execute <- function(sql) invisible(dbExecute(con, sql))
get <- function(sql) dbGetQuery(con, sql)
execute("SET threads=2")
execute("SET memory_limit='4GB'")
execute(paste0("SET temp_directory=", q(file.path(out, "spill"))))
execute("SET max_temp_directory_size='12GiB'")
execute(paste0("ATTACH ", q(normalizePath(model)), " AS model (READ_ONLY)"))
execute(paste0("CREATE TEMP VIEW source_alleles AS SELECT * FROM read_parquet([",
  paste(vapply(files, q, character(1L)), collapse = ","), "])"))
if (get("SELECT count(*) AS n FROM source_alleles a
  LEFT JOIN model.duckvep_sequence_regions r ON r.name = regexp_replace(a.chrom, '^chr', '')
  WHERE r.seq_region IS NULL OR a.seq_region IS DISTINCT FROM r.seq_region")$n != 0) {
  stop("source region index differs from model")
}
execute("CREATE TEMP VIEW distinct_alleles AS SELECT * EXCLUDE(source_rank) FROM (
  SELECT *, row_number() OVER (PARTITION BY chrom, position, reference, alternate
    ORDER BY source_object, record_index, alt_index) AS source_rank FROM source_alleles
  ) WHERE source_rank = 1")
quotas <- data.frame(bin = c("cds_snv", "cds_indel", "cds_mnv", "splice_remainder"),
  required = c(1000000L, 500000L, 100000L, 400000L))
upper <- get("SELECT
  count(*) FILTER (WHERE status = 'literal' AND length(reference) = 1
    AND length(alternate) = 1) AS snv,
  count(*) FILTER (WHERE status = 'literal'
    AND length(reference) != length(alternate)) AS indel,
  count(*) FILTER (WHERE status = 'literal' AND length(reference) = length(alternate)
    AND length(reference) > 1) AS mnv,
  count(*) FILTER (WHERE status = 'literal') AS total FROM distinct_alleles")
bounds <- data.frame(bin = quotas$bin, required = quotas$required,
  source_upper_bound = as.numeric(unlist(upper[1L, c("snv", "indel", "mnv", "total")],
    use.names = FALSE)))
write.table(bounds, file.path(out, "bounds.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (any(bounds$source_upper_bound < bounds$required)) {
  write.table(data.frame(rows = 0L, sha256 = "unavailable",
    status = "insufficient_source_upper_bounds"), file.path(out, "panel.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  dbDisconnect(con, shutdown = TRUE)
  cat(out, "\n")
  quit(save = "no", status = 0L)
}
execute("CREATE TEMP TABLE intervals AS SELECT t.seq_region,
  CASE WHEN t.cds_start IS NOT NULL AND t.cds_end IS NOT NULL
    THEN greatest(e.exon_start, t.cds_start) END AS cds_first,
  CASE WHEN t.cds_start IS NOT NULL AND t.cds_end IS NOT NULL
    THEN least(e.exon_end, t.cds_end) END AS cds_last,
  CASE WHEN e.exon_cdna_start > 1 THEN e.exon_start END AS exon_first,
  CASE WHEN e.exon_cdna_end < max(e.exon_cdna_end) OVER (PARTITION BY e.transcript_index)
    THEN e.exon_end END AS exon_last
  FROM model.duckvep_exons e JOIN model.model_transcripts t USING(transcript_index)
  WHERE t.seq_region IN (SELECT DISTINCT seq_region FROM source_alleles)")
execute("CREATE TEMP VIEW classified AS WITH annotated AS (
 SELECT v.*, EXISTS (
   SELECT 1 FROM intervals i WHERE i.seq_region = v.seq_region
     AND i.cds_first IS NOT NULL AND i.cds_first <= i.cds_last AND v.position <= i.cds_last
     AND v.position + length(v.reference) - 1 >= i.cds_first) AS in_cds,
   EXISTS (SELECT 1 FROM intervals i WHERE i.seq_region = v.seq_region AND
     (v.position BETWEEN greatest(1, i.exon_first - 8) AND i.exon_first + 8
      OR v.position BETWEEN greatest(1, i.exon_last - 8) AND i.exon_last + 8)) AS near_splice
 FROM distinct_alleles v WHERE v.status = 'literal'
 ) SELECT *, CASE
 WHEN in_cds AND length(reference) = 1 AND length(alternate) = 1 THEN 'cds_snv'
 WHEN in_cds AND length(reference) != length(alternate) THEN 'cds_indel'
 WHEN in_cds AND length(reference) = length(alternate) AND length(reference) > 1 THEN 'cds_mnv'
 WHEN near_splice THEN 'splice_remainder'
 ELSE 'other' END AS bin FROM annotated")
availability <- get("SELECT bin, count(*) AS available FROM classified GROUP BY bin ORDER BY bin")
availability <- merge(quotas, availability, all.x = TRUE, sort = FALSE)
availability$available[is.na(availability$available)] <- 0
write.table(availability, file.path(out, "availability.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (all(availability$available >= availability$required)) {
  execute("CREATE TEMP VIEW ranked AS SELECT *, row_number() OVER (PARTITION BY bin ORDER BY
    sha256('gnomad-v4.1-exome-v1|GRCh38|' || length(chrom)::VARCHAR || ':' || chrom ||
      position::VARCHAR || ':' || length(reference)::VARCHAR || ':' || reference ||
      length(alternate)::VARCHAR || ':' || alternate), source_object, record_index, alt_index
    ) AS rank FROM classified")
  target <- file.path(out, "exome-2m.parquet")
  execute(paste0("COPY (SELECT * EXCLUDE(in_cds, near_splice, bin, rank) FROM ranked
    WHERE (bin = 'cds_snv' AND rank <= 1000000)
       OR (bin = 'cds_indel' AND rank <= 500000)
       OR (bin = 'cds_mnv' AND rank <= 100000)
       OR (bin = 'splice_remainder' AND rank <= 400000)
    ORDER BY seq_region, position, source_object, record_index, alt_index) TO ",
    q(target), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  rows <- get(paste0("SELECT count(*) AS n FROM read_parquet(", q(target), ")"))$n
  if (rows != 2000000) stop("exome panel did not meet its quota")
  receipt <- data.frame(rows = rows, sha256 = sha256(target), status = "complete")
} else receipt <- data.frame(rows = 0L, sha256 = "unavailable", status = "insufficient_bin_alleles")
write.table(receipt, file.path(out, "panel.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
cat(out, "\n")
