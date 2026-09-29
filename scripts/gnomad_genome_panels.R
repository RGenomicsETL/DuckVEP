library(DBI)
library(duckdb)
root <- Sys.getenv("DUCKVEP_GNOMAD_ROOT", "/root/duckvep/data/gnomad-v4.1")
shards <- list.dirs(root, recursive = FALSE, full.names = TRUE)
shards <- shards[grepl("/genomes-chr(1[0-9]|2[0-2]|[1-9]|X|Y)$", shards)]
if (length(shards) == 0L || any(!file.exists(file.path(shards, "source.tsv")))) {
  stop("no verified genome shards")
}
files <- unlist(lapply(sort(shards), function(dir) {
  parts <- list.files(dir, "^part-[0-9]+\\.parquet$", full.names = TRUE)
  if (length(parts) == 0L) stop("empty verified genome shard: ", dir)
  sort(parts)
}), use.names = FALSE)
sha256 <- function(path) substr(system2("sha256sum", path, stdout = TRUE), 1L, 64L)
source_receipts <- data.frame(file = files, bytes = file.size(files),
  sha256 = vapply(files, sha256, character(1L)))
source_key <- tempfile("gnomad-source-", tmpdir = root)
write.table(source_receipts, source_key, row.names = FALSE, sep = "\t", quote = FALSE)
source_sha <- sha256(source_key)
unlink(source_key)
out <- file.path(root, "panels", paste0("genomes-v1-", substr(source_sha, 1L, 16L)))
if (dir.exists(out)) {
  # Same sources give the same panels: reuse a finished artifact, refuse a half-built one.
  if (!file.exists(file.path(out, "panels.tsv"))) stop("incomplete panel artifact, remove it: ", out)
  cat(out, "\n")
  quit(save = "no", status = 0L)
}
dir.create(out, recursive = TRUE)
write.table(source_receipts, file.path(out, "sources.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
con <- dbConnect(duckdb(shared_home = FALSE))
q <- function(x) as.character(dbQuoteString(con, x))
execute <- function(sql) invisible(dbExecute(con, sql))
get <- function(sql) dbGetQuery(con, sql)
execute(paste0("SET threads=", Sys.getenv("DUCKVEP_PANEL_THREADS", "2")))
execute(paste0("SET memory_limit=", q(Sys.getenv("DUCKVEP_PANEL_MEMORY", "4GB"))))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "gnomad_panel_tmp.R"))
panel_spill_setup(execute, q)
paths <- paste(vapply(files, q, character(1L)), collapse = ",")
execute(paste0("CREATE TEMP VIEW source_alleles AS SELECT * FROM read_parquet([", paths, "])"))
model <- Sys.getenv("DUCKVEP_SCALE_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
execute(paste0("ATTACH ", q(normalizePath(model)), " AS model (READ_ONLY)"))
if (get("SELECT count(*) AS n FROM source_alleles a
  LEFT JOIN model.duckvep_sequence_regions r ON r.name = regexp_replace(a.chrom, '^chr', '')
  WHERE r.seq_region IS NULL OR a.seq_region IS DISTINCT FROM r.seq_region")$n != 0) {
  stop("source region index differs from model")
}
# Ranking runs on narrow keys only; the wide rows are joined back for the selected ranks.
execute("CREATE TEMP VIEW source_keys AS SELECT chrom, position, reference, alternate, status,
  source_object, record_index, alt_index FROM source_alleles")
execute("CREATE TEMP VIEW distinct_alleles AS
  SELECT * EXCLUDE(source_rank) FROM (
    SELECT *, row_number() OVER (PARTITION BY chrom, position, reference, alternate
      ORDER BY source_object, record_index, alt_index) AS source_rank FROM source_keys
  ) WHERE source_rank = 1")
counts <- get("SELECT (SELECT count(*) FROM source_alleles) AS source_alts,
  (SELECT count(*) FROM distinct_alleles) AS distinct_alts")
if (counts$source_alts < counts$distinct_alts) stop("deduplication increased rows")
map <- file.path(out, "dedup-map.parquet")
execute(paste0("COPY (SELECT source_object, record_index, alt_index, chrom, position,
  reference, alternate, kept_object, kept_record, kept_alt FROM (
    SELECT source_object, record_index, alt_index, chrom, position, reference, alternate,
      row_number() OVER (PARTITION BY chrom, position, reference, alternate
        ORDER BY source_object, record_index, alt_index) AS source_rank,
      first_value(source_object) OVER (PARTITION BY chrom, position, reference, alternate
        ORDER BY source_object, record_index, alt_index) AS kept_object,
      first_value(record_index) OVER (PARTITION BY chrom, position, reference, alternate
        ORDER BY source_object, record_index, alt_index) AS kept_record,
      first_value(alt_index) OVER (PARTITION BY chrom, position, reference, alternate
        ORDER BY source_object, record_index, alt_index) AS kept_alt
    FROM source_alleles) WHERE source_rank > 1
    ORDER BY source_object, record_index, alt_index) TO ", q(map),
    " (FORMAT PARQUET, COMPRESSION ZSTD)"))
if (get(paste0("SELECT count(*) AS n FROM read_parquet(", q(map), ")"))$n !=
    counts$source_alts - counts$distinct_alts) stop("dedup map does not conserve rows")
write.table(data.frame(rows = counts$source_alts - counts$distinct_alts,
  sha256 = sha256(map)), file.path(out, "dedup-map.tsv"), row.names = FALSE,
  sep = "\t", quote = FALSE)
quotas <- c(1000000L, 5000000L, 25000000L, 100000000L)
execute(paste0("CREATE TEMP TABLE panel_keys AS SELECT source_object, record_index, alt_index, panel_rank FROM (WITH classified AS (
  SELECT *, CASE WHEN status <> 'literal' THEN status
    WHEN length(reference) = 1 AND length(alternate) = 1 THEN 'snv'
    WHEN length(reference) <> length(alternate) THEN 'indel'
    ELSE 'mnv' END AS allele_class FROM distinct_alleles
 ), ranked AS (
  SELECT *, sha256('gnomad-v4.1-panel-v1|GRCh38|' || length(chrom)::VARCHAR || ':' || chrom ||
    position::VARCHAR || ':' || length(reference)::VARCHAR || ':' || reference ||
    length(alternate)::VARCHAR || ':' || alternate) AS allele_rank,
    row_number() OVER (PARTITION BY chrom, allele_class ORDER BY
      sha256('gnomad-v4.1-panel-v1|GRCh38|' || length(chrom)::VARCHAR || ':' || chrom ||
      position::VARCHAR || ':' || length(reference)::VARCHAR || ':' || reference ||
      length(alternate)::VARCHAR || ':' || alternate)) AS stratum_rank,
    count(*) OVER (PARTITION BY chrom, allele_class) AS stratum_size FROM classified
  ) SELECT *, row_number() OVER (ORDER BY (stratum_rank - 0.5) / stratum_size,
    allele_rank, source_object, record_index, alt_index) AS panel_rank FROM ranked) WHERE panel_rank <= ", max(quotas)))
results <- lapply(quotas, function(n) {
  if (counts$distinct_alts < n) return(data.frame(quota = n, rows = 0L,
    sha256 = "unavailable", status = "insufficient_distinct_alleles"))
  path <- file.path(out, paste0("panel-", n, ".parquet"))
  execute(paste0("COPY (SELECT s.* FROM source_alleles s JOIN panel_keys k USING (source_object, record_index, alt_index) ",
    "WHERE k.panel_rank <= ", n,
    " ORDER BY s.seq_region, s.position, s.source_object, s.record_index, s.alt_index) TO ",
    q(path), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  rows <- get(paste0("SELECT count(*) AS n FROM read_parquet(", q(path), ")"))$n
  if (rows != n) stop("panel row count differs from quota")
  size <- as.numeric(strsplit(system2("du", c("-sb", "--exclude=spill", root), stdout = TRUE),
    "\t", fixed = TRUE)[[1L]][1L])
  if (size > 25e9) stop("25 GB staged-output limit reached")
  data.frame(quota = n, rows = format(rows, scientific = FALSE),
    sha256 = sha256(path), status = "complete")
})
write.table(do.call(rbind, results), file.path(out, "panels.tsv"), row.names = FALSE,
  sep = "\t", quote = FALSE)
write.table(counts, file.path(out, "counts.tsv"), row.names = FALSE,
  sep = "\t", quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
cat(out, "\n")
