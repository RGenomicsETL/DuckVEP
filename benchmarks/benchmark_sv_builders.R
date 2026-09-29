#!/usr/bin/env Rscript
# Time and peak resident memory of the five native SV/STR/BND preparation builders.
#
#   Rscript benchmarks/benchmark_sv_builders.R EXTENSION_COPY OUTPUT_CSV [N] [SNIFFLES_VCF_GZ]
#
# EXTENSION_COPY must be an immutable copy outside build/release. Every
# (builder, corpus) pair runs in its own R process (one DuckDB thread) so that
# VmHWM is that pair's peak resident set. The timed step is
#   CREATE TEMP TABLE out AS SELECT * FROM query(<builder SQL>)
# on an already-materialized input table; input generation is outside it. Each
# pair is timed three times and the median is reported.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(trailingOnly = TRUE)
rss <- function(field) {
  line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", line)) / 1024
}

## ---- child: one builder on one corpus -------------------------------------------------
if (identical(args[[1L]], "--child")) {
  builder <- args[[2L]]; corpus <- args[[3L]]; n <- as.integer(args[[4L]])
  extension <- args[[5L]]; sniffles <- args[[6L]]
  con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  run <- function(sql) invisible(dbExecute(con, sql))
  run(paste0("LOAD '", extension, "'")); run("SET threads=1")
  seq_expr <- function(len_expr, salt) sprintf(
    "array_to_string(list_transform(range((%s)::BIGINT), lambda x: substr('ACGT', 1 + (hash(event_index, x, %d) %% 4)::INTEGER, 1)), '')",
    len_expr, salt)
  if (corpus == "sniffles2_chr22") {
    run(sprintf("CREATE TEMP TABLE raw AS SELECT * FROM read_csv('%s', delim = '\t', quote = '', escape = '', comment = '#', header = false, auto_detect = false, strict_mode = false,
      columns = {'chrom': 'VARCHAR', 'pos': 'BIGINT', 'id': 'VARCHAR', 'ref': 'VARCHAR', 'alt': 'VARCHAR', 'qual': 'VARCHAR',
      'filter': 'VARCHAR', 'info': 'VARCHAR'})", sniffles))
  }
  input <- switch(paste(builder, corpus),
    "sv_geometry sniffles2_chr22" = "CREATE TEMP TABLE input AS SELECT (row_number() OVER () - 1)::BIGINT AS event_index, pos, ref, alt, info FROM raw",
    "breakend_pairs sniffles2_chr22" = "CREATE TEMP TABLE input AS SELECT (row_number() OVER () - 1)::BIGINT AS event_index, chrom, pos, id, ref, alt, info FROM raw",
    "sv_geometry synthetic" = sprintf("CREATE TEMP TABLE input AS SELECT event_index, 1000 + hash(event_index) %% 200000000 AS pos,
      CASE event_index %% 6 WHEN 5 THEN 'A' ELSE 'N' END AS ref,
      CASE event_index %% 6 WHEN 0 THEN '<DEL>' WHEN 1 THEN '<DUP>' WHEN 2 THEN '<INV>' WHEN 3 THEN '<INS>' WHEN 4 THEN '<DUP:TANDEM>' ELSE 'AACGT' END AS alt,
      'END=' || (1000 + hash(event_index) %% 200000000 + 1 + hash(event_index, 1) %% 100000) ||
        CASE WHEN event_index %% 5 = 0 THEN ';CIPOS=-10,12;CIEND=-5,7' ELSE '' END ||
        CASE WHEN event_index %% 6 = 3 THEN ';SEQ=' || %s ELSE '' END AS info
      FROM range(%d) t(event_index)", seq_expr("1 + hash(event_index, 2) % 50", 3), n),
    "breakend_pairs synthetic" = sprintf("CREATE TEMP TABLE input AS WITH p AS (SELECT p AS pair, 1 + hash(p) %% 22 AS c1, 1 + hash(p, 1) %% 22 AS c2,
      1000 + hash(p, 2) %% 100000000 AS p1, 1000 + hash(p, 3) %% 100000000 AS p2, p %% 4 AS shape, p %% 50 AS defect FROM range(%d) t(p)),
      recs AS (SELECT pair * 2 AS event_index, c1::VARCHAR AS chrom, p1 AS pos, 'p' || pair || 'a' AS id, 'N' AS ref,
        CASE shape WHEN 0 THEN 'N[' || c2 || ':' || p2 || '[' WHEN 1 THEN ']' || c2 || ':' || p2 || ']N'
          WHEN 2 THEN 'N]' || c2 || ':' || p2 || ']' ELSE '[' || c2 || ':' || p2 || '[N' END AS alt,
        'SVTYPE=BND;MATEID=p' || pair || 'b;EVENT=e' || pair AS info FROM p
      UNION ALL SELECT pair * 2 + 1, c2::VARCHAR, p2, 'p' || pair || 'b', 'N',
        CASE shape WHEN 0 THEN ']' || c1 || ':' || p1 || ']N' WHEN 1 THEN 'N[' || c1 || ':' || p1 || '['
          WHEN 2 THEN 'N]' || c1 || ':' || p1 || ']' ELSE '[' || c1 || ':' || p1 || '[N' END,
        'SVTYPE=BND;MATEID=p' || pair || 'a;EVENT=' || CASE WHEN defect = 7 THEN 'other' ELSE 'e' || pair END FROM p)
      SELECT event_index, chrom, pos, id, ref, CASE WHEN event_index %% 97 = 5 THEN '.N' ELSE alt END AS alt, info FROM recs", n %/% 2L),
    "structural_hgvs synthetic" = sprintf("CREATE TEMP TABLE input AS SELECT event_index, '1' AS chrom, 1000 + hash(event_index) %% 200000000 AS pos,
      'N' AS ref, CASE event_index %% 5 WHEN 0 THEN '<DEL>' WHEN 1 THEN '<DUP>' WHEN 2 THEN '<INV>' WHEN 3 THEN '<DUP:TANDEM>' ELSE '<CNV>' END AS alt,
      'END=' || (1000 + hash(event_index) %% 200000000 + 1 + hash(event_index, 1) %% 40) AS info
      FROM range(%d) t(event_index)", n),
    "expansionhunter synthetic" = sprintf("CREATE TEMP TABLE input AS SELECT event_index,
      'END=' || (2000 + event_index) || ';REF=' || (1 + hash(event_index) %% 20) || ';RL=' || 3 * (1 + hash(event_index) %% 20) || ';RU=CAG' AS info,
      CASE event_index %% 2 WHEN 0 THEN 'GT:SO:REPCN:REPCI' ELSE 'GT:SO:CN:CI' END AS format,
      CASE event_index %% 2 WHEN 0 THEN '1/2:SPANNING/SPANNING:2/10:2-2/10-10' ELSE '1/2:SPANNING/INREPEAT:2/349:2-2/323-376' END AS \"sample\",
      'C' AS ref, CASE event_index %% 2 WHEN 0 THEN '<STR2>,<STR10>' ELSE '<STR2>,<STR349>' END AS alt,
      (1 + event_index %% 2)::DOUBLE AS alt_index FROM range(%d) t(event_index)", n),
    "breakend_fusion synthetic" = NULL, stop("unknown pair ", builder, " ", corpus))
  reference <- NULL
  if (builder == "structural_hgvs") {
    run(input); run("CREATE TEMP TABLE reference AS SELECT event_index,
      substr(array_to_string(list_transform(range((2 * (try_cast(regexp_extract(info, 'END=([0-9]+)', 1) AS BIGINT) - pos) + 1)::BIGINT),
      lambda x: substr('ACGT', 1 + (hash(event_index, x, 9) % 4)::INTEGER, 1)), ''), 1) AS reference_sequence FROM input");
    reference <- "reference"
  } else if (builder == "expansionhunter") {
    run(input); run("CREATE TEMP TABLE reference AS SELECT event_index,
      CASE WHEN event_index % 10 = 3 THEN repeat('CAG', (1 + hash(event_index) % 20 + 1)::BIGINT) ELSE repeat('CAG', (1 + hash(event_index) % 20)::BIGINT) END AS reference_sequence FROM input")
    reference <- "reference"
  } else if (builder == "breakend_fusion") {
    # The fusion builder's input is the pair-identity table plus endpoint genes.
    ids <- sprintf("CREATE TEMP TABLE raw_pairs AS WITH p AS (SELECT p AS pair, 1 + hash(p) %% 22 AS c1, 1 + hash(p, 1) %% 22 AS c2,
      1000 + hash(p, 2) %% 100000000 AS p1, 1000 + hash(p, 3) %% 100000000 AS p2, p %% 4 AS shape FROM range(%d) t(p))
      SELECT pair * 2 AS event_index, c1::VARCHAR AS chrom, p1 AS pos, 'p' || pair || 'a' AS id, 'N' AS ref,
        CASE shape WHEN 0 THEN 'N[' || c2 || ':' || p2 || '[' WHEN 1 THEN ']' || c2 || ':' || p2 || ']N'
          WHEN 2 THEN 'N]' || c2 || ':' || p2 || ']' ELSE '[' || c2 || ':' || p2 || '[N' END AS alt,
        'SVTYPE=BND;MATEID=p' || pair || 'b;EVENT=e' || pair AS info FROM p
      UNION ALL SELECT pair * 2 + 1, c2::VARCHAR, p2, 'p' || pair || 'b', 'N',
        CASE shape WHEN 0 THEN ']' || c1 || ':' || p1 || ']N' WHEN 1 THEN 'N[' || c1 || ':' || p1 || '['
          WHEN 2 THEN 'N]' || c1 || ':' || p1 || ']' ELSE '[' || c1 || ':' || p1 || '[N' END,
        'SVTYPE=BND;MATEID=p' || pair || 'a;EVENT=e' || pair FROM p", n %/% 2L)
    run(ids)
    run("CREATE TEMP TABLE input AS SELECT * FROM query(duckvep_prepare_breakend_pairs_sql('raw_pairs'))")
    run(sprintf("CREATE TEMP TABLE genes AS SELECT event_index, 'G' || (hash(event_index, k) %% 40000) AS gene_id
      FROM range(%d) t(event_index), range(0, 3) u(k) WHERE hash(event_index, k, 5) %% 3 <> 0", n))
    reference <- "genes"
  }
  if (is.null(reference) && builder != "breakend_fusion") run(input)
  table_name <- "input"
  sql_call <- switch(builder,
    sv_geometry = "duckvep_prepare_sv_geometry_sql('input')",
    breakend_pairs = "duckvep_prepare_breakend_pairs_sql('input')",
    breakend_fusion = "duckvep_prepare_breakend_fusion_sql('input', 'genes')",
    structural_hgvs = "duckvep_prepare_structural_hgvs_sql('input', 'reference')",
    expansionhunter = "duckvep_prepare_expansionhunter_sql('input', 'reference')")
  rows_in <- dbGetQuery(con, "SELECT count(*) n FROM input")$n
  before <- rss("VmHWM")
  sql_time <- system.time(dbGetQuery(con, paste("SELECT", sql_call)))[["elapsed"]]
  times <- vapply(1:3, function(i) {
    run("DROP TABLE IF EXISTS out")
    system.time(run(paste0("CREATE TEMP TABLE out AS SELECT * FROM query(", sql_call, ")")))[["elapsed"]]
  }, 0)
  rows_out <- dbGetQuery(con, "SELECT count(*) n FROM out")$n
  status_col <- if (builder == "structural_hgvs") "hgvs_status" else "status"
  status <- dbGetQuery(con, sprintf("SELECT %s AS s, count(*) n FROM out GROUP BY 1 ORDER BY 1", status_col))
  stopifnot(rows_out == rows_in)
  cat(paste(builder, corpus, format(rows_in, scientific = FALSE, trim = TRUE), format(rows_out, scientific = FALSE, trim = TRUE), sprintf("%.4f", sql_time), sprintf("%.3f", median(times)),
    sprintf("%.3f", min(times)), sprintf("%.3f", max(times)), sprintf("%.0f", before), sprintf("%.0f", rss("VmHWM")),
    paste(status$s, format(status$n, scientific = FALSE, trim = TRUE), sep = ":", collapse = ";"), sep = "\t"), "\n", sep = "")
  quit(status = 0L)
}

## ---- parent -------------------------------------------------------------------------------
extension <- normalizePath(args[[1L]], mustWork = TRUE)
if (grepl("/build/release/", extension, fixed = TRUE)) {
  stop("EXTENSION_COPY must be an immutable copy outside build/release", call. = FALSE)
}
output <- args[[2L]]
n <- if (length(args) >= 3L) as.integer(args[[3L]]) else 100000L
sniffles <- if (length(args) >= 4L) normalizePath(args[[4L]], mustWork = TRUE) else
  "/root/.cache/duckhts/corpora/duckvep/sniffles2_1kgp/sniffles2_joint_chr22_sites.vcf.gz"
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
plain <- tempfile(fileext = ".vcf"); writeLines(system2("gzip", c("-dc", shQuote(sniffles)), stdout = TRUE), plain)
plan <- list(c("sv_geometry", "sniffles2_chr22"), c("breakend_pairs", "sniffles2_chr22"),
  c("sv_geometry", "synthetic"), c("breakend_pairs", "synthetic"), c("breakend_fusion", "synthetic"),
  c("structural_hgvs", "synthetic"), c("expansionhunter", "synthetic"))
lines <- character()
for (item in plan) {
  cat("running", item, "\n")
  out <- system2("Rscript", c("--vanilla", shQuote(script), "--child", item[[1L]], item[[2L]], n, shQuote(extension), shQuote(plain)),
    stdout = TRUE, stderr = FALSE)
  row <- tail(grep("\t", out, value = TRUE), 1L)
  stopifnot(length(row) == 1L)
  lines <- c(lines, row)
}
header <- c("builder", "corpus", "rows_in", "rows_out", "sql_build_s", "median_s", "min_s", "max_s",
            "peak_rss_before_mb", "peak_rss_mb", "status_counts")
table <- read.delim(text = paste(c(paste(header, collapse = "\t"), lines), collapse = "\n"), sep = "\t", stringsAsFactors = FALSE)
table$rows_per_second <- round(table$rows_in / table$median_s)
table$extension_sha256 <- tolower(strsplit(system2("sha256sum", shQuote(extension), stdout = TRUE), " ")[[1L]][1L])
table$duckdb_r_package <- as.character(packageVersion("duckdb"))
table$synthetic_events <- n
write.csv(table, output, row.names = FALSE)
print(table[c("builder", "corpus", "rows_in", "rows_out", "median_s", "rows_per_second", "peak_rss_before_mb", "peak_rss_mb")])
