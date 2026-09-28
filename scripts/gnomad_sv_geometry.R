library(DBI)
library(duckdb)
root <- "/root/duckvep/data/gnomad-v4.1/sv-sites"
if (!file.exists(file.path(root, "source.tsv"))) stop("SV source is not verified")
output <- file.path(root, "geometry.parquet")
if (file.exists(output)) stop("SV geometry already exists")
con <- dbConnect(duckdb(shared_home = FALSE))
q <- function(x) as.character(dbQuoteString(con, x))
execute <- function(sql) invisible(dbExecute(con, sql))
source <- paste0("read_parquet(", q(file.path(root, "part-*.parquet")), ")")
execute(paste0("CREATE TEMP VIEW geometry AS SELECT a.record_index, a.alt_index,
  CASE WHEN a.status = 'symbolic' AND a.svtype IN ('BND', 'CTX')
         AND contains(a.filter, 'UNRESOLVED') THEN 'unresolved_symbolic_breakpoint'
       WHEN a.status = 'symbolic' AND a.svtype IN ('BND', 'CTX')
         THEN 'symbolic_breakpoint'
       WHEN a.status = 'symbolic' AND (a.end_position IS NULL OR a.end_position < a.position)
         THEN 'missing_or_invalid_end'
       WHEN a.status = 'symbolic' AND (a.cipos IS NOT NULL OR a.ciend IS NOT NULL)
         THEN 'imprecise'
       WHEN a.status = 'breakend' AND (a.mateid IS NULL OR a.mateid = '')
         THEN 'missing_mate_id'
       WHEN a.status = 'breakend' AND NOT EXISTS (
         SELECT 1 FROM ", source, " b WHERE b.id = a.mateid AND b.mateid = a.id
           AND b.record_index <> a.record_index)
         THEN 'unpaired_mate'
       WHEN a.status = 'breakend' THEN 'paired'
       WHEN a.status = 'symbolic' THEN 'bounded'
       ELSE a.status END AS geometry_status
 FROM ", source, " a"))
execute(paste0("COPY (SELECT * FROM geometry ORDER BY record_index, alt_index) TO ",
  q(output), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
counts <- dbGetQuery(con, paste0("SELECT geometry_status, count(*) AS alts FROM read_parquet(",
  q(output), ") GROUP BY geometry_status ORDER BY geometry_status"))
all_rows <- dbGetQuery(con, paste0("SELECT count(*) AS n FROM ", source))$n
if (sum(counts$alts) != all_rows) stop("SV geometry did not conserve alleles")
write.table(counts, file.path(root, "geometry.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
cat(output, "\n")
