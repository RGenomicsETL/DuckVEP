#!/usr/bin/env Rscript

# Compare diagnostic relations with physical record/ALT ordinals. Empty strings,
# SQL NULL, dashes, ordered lists and percent escapes are distinct field values.
duckvep_fastvep_compare <- function(con, actual, expected, source, fields, directory) {
  allele_keys <- c("record_index", "alt_index")
  keys <- c("record_index", "alt_index", "Feature")
  columns <- unique(c(keys, fields))
  if (anyDuplicated(fields) || !all(keys %in% columns) ||
      !all(columns %in% DBI::dbListFields(con, actual)) ||
      !all(columns %in% DBI::dbListFields(con, expected))) {
    stop("comparison relations must contain the declared fields and physical keys")
  }
  if (!all(allele_keys %in% DBI::dbListFields(con, source))) {
    stop("independent source relation must contain record_index and alt_index")
  }
  if (dir.exists(directory) || !dir.create(directory, recursive = TRUE)) {
    stop("comparison artifacts require a new output directory")
  }
  qi <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  qs <- function(x) as.character(DBI::dbQuoteString(con, x))
  execute <- function(sql) invisible(DBI::dbExecute(con, sql))
  select <- function(alias, names) paste0(alias, ".", qi(names), collapse = ", ")
  same_keys <- function(left, right, names = keys) paste(paste0(left, ".", qi(names),
    " IS NOT DISTINCT FROM ", right, ".", qi(names)), collapse = " AND ")
  artifact <- function(sql, name) execute(paste0("COPY (", sql, ") TO ",
    qs(file.path(directory, paste0(name, ".parquet"))), " (FORMAT PARQUET)"))
  grouped <- function(table) paste0("SELECT ", paste(qi(keys), collapse = ", "),
    ", count(*) AS copies FROM ", qi(table), " GROUP BY ALL")
  table <- "_fastvep_comparison_keys"
  coverage <- "_fastvep_source_coverage"
  if (DBI::dbExistsTable(con, table) || DBI::dbExistsTable(con, coverage)) {
    stop("comparison workspace is already in use")
  }
  execute(paste0("CREATE TEMP TABLE ", table, " AS SELECT ",
    paste(paste0("coalesce(a.", qi(keys), ", e.", qi(keys), ") AS ", qi(keys)), collapse = ", "),
    ", coalesce(a.copies, 0) AS actual_rows, coalesce(e.copies, 0) AS expected_rows
    FROM (", grouped(actual), ") a FULL OUTER JOIN (", grouped(expected),
    ") e ON ", same_keys("a", "e")))
  on.exit(DBI::dbRemoveTable(con, table), add = TRUE)
  # Coverage is anchored to the independent source, including alleles with no
  # output on either side. Transcript multiplicity does not multiply this count.
  execute(paste0("CREATE TEMP TABLE ", coverage, " AS
    WITH declared AS (
      SELECT record_index::VARCHAR AS record_index, alt_index::VARCHAR AS alt_index,
        count(*) AS source_rows FROM ", qi(source), " GROUP BY ALL
    ), emitted AS (
      SELECT record_index::VARCHAR AS record_index, alt_index::VARCHAR AS alt_index,
        sum(actual_rows) AS actual_rows, sum(expected_rows) AS expected_rows
      FROM ", table, " GROUP BY ALL
    ) SELECT coalesce(s.record_index, o.record_index) AS record_index,
      coalesce(s.alt_index, o.alt_index) AS alt_index,
      coalesce(s.source_rows, 0) AS source_rows,
      coalesce(o.actual_rows, 0) AS actual_rows, coalesce(o.expected_rows, 0) AS expected_rows
    FROM declared s FULL OUTER JOIN emitted o ON ", same_keys("s", "o", allele_keys)))
  on.exit(DBI::dbRemoveTable(con, coverage), add = TRUE)
  valid_ordinal <- function(name) paste0("(", name, " IS NOT NULL AND regexp_full_match(",
    name, ", '[1-9][0-9]*') AND try_cast(", name, " AS UBIGINT) IS NOT NULL)")
  valid_allele <- paste(vapply(allele_keys, valid_ordinal, character(1L)), collapse = " AND ")
  artifact(paste0("SELECT *, (", valid_allele, ") AS valid_allele FROM ", coverage,
    " WHERE source_rows != 1 OR actual_rows = 0 OR expected_rows = 0 OR NOT (",
    valid_allele, ")"), "source_coverage_failures")
  artifact(paste0("SELECT * FROM ", table,
    " WHERE actual_rows != 1 OR expected_rows != 1
       OR record_index IS NULL OR alt_index IS NULL"), "key_failures")
  # Keep every duplicate, not one arbitrary representative or a Cartesian join.
  artifact(paste0("SELECT 'actual' AS side, ", select("a", columns), " FROM ", qi(actual),
    " a JOIN ", table, " k ON ", same_keys("a", "k"), " WHERE k.actual_rows > 1
    UNION ALL SELECT 'expected' AS side, ", select("e", columns), " FROM ", qi(expected),
    " e JOIN ", table, " k ON ", same_keys("e", "k"), " WHERE k.expected_rows > 1"),
    "duplicate_rows")
  compared <- setdiff(fields, keys)
  if (!length(compared)) stop("at least one payload field must be compared")
  values <- paste(vapply(compared, function(field) paste0("(", qs(field),
    ", a.", qi(field), "::VARCHAR, e.", qi(field), "::VARCHAR)"), character(1L)), collapse = ", ")
  joined <- paste0(" FROM ", table, " k JOIN ", qi(actual), " a ON ", same_keys("k", "a"),
    " JOIN ", qi(expected), " e ON ", same_keys("k", "e"),
    ", LATERAL (VALUES ", values, ") v(field, actual, expected)
    WHERE k.actual_rows = 1 AND k.expected_rows = 1
      AND k.record_index IS NOT NULL AND k.alt_index IS NOT NULL")
  artifact(paste0("SELECT ", select("k", keys), ", v.*", joined,
    " AND v.actual IS DISTINCT FROM v.expected"), "field_failures")
  summary <- DBI::dbGetQuery(con, paste0("SELECT
    coalesce(sum(actual_rows), 0) AS actual_rows,
    coalesce(sum(expected_rows), 0) AS expected_rows,
    count(*) AS union_keys,
    count(*) FILTER (WHERE actual_rows = 0) AS missing_keys,
    count(*) FILTER (WHERE expected_rows = 0) AS extra_keys,
    count(*) FILTER (WHERE actual_rows > 1) AS actual_duplicate_keys,
    count(*) FILTER (WHERE expected_rows > 1) AS expected_duplicate_keys,
    count(*) FILTER (WHERE record_index IS NULL OR alt_index IS NULL) AS invalid_keys,
    count(*) FILTER (WHERE actual_rows = 1 AND expected_rows = 1
      AND record_index IS NOT NULL AND alt_index IS NOT NULL) AS compared_keys
    FROM ", table))
  source_summary <- DBI::dbGetQuery(con, paste0("SELECT
    coalesce(sum(source_rows), 0) AS source_rows,
    count(*) FILTER (WHERE source_rows > 0) AS source_alleles,
    count(*) FILTER (WHERE source_rows > 1) AS source_duplicate_alleles,
    count(*) FILTER (WHERE source_rows > 0 AND NOT (", valid_allele, ")) AS source_invalid_alleles,
    count(*) FILTER (WHERE source_rows > 0 AND actual_rows = 0) AS actual_missing_source_alleles,
    count(*) FILTER (WHERE source_rows > 0 AND expected_rows = 0) AS expected_missing_source_alleles,
    count(*) FILTER (WHERE source_rows = 0 AND actual_rows > 0) AS actual_unknown_alleles,
    count(*) FILTER (WHERE source_rows = 0 AND expected_rows > 0) AS expected_unknown_alleles
    FROM ", coverage))
  summary <- cbind(summary, source_summary)
  summary$compared_fields <- length(compared)
  failures <- DBI::dbGetQuery(con, paste0("SELECT field, count(*) AS failures FROM read_parquet(",
    qs(file.path(directory, "field_failures.parquet")), ") GROUP BY field ORDER BY field"))
  summary$field_failures <- sum(failures$failures)
  error_columns <- c("missing_keys", "extra_keys", "actual_duplicate_keys",
    "expected_duplicate_keys", "invalid_keys", "field_failures", "source_duplicate_alleles",
    "source_invalid_alleles", "actual_missing_source_alleles", "expected_missing_source_alleles",
    "actual_unknown_alleles", "expected_unknown_alleles")
  summary$passed <- all(unlist(summary[error_columns]) == 0)
  utils::write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  utils::write.csv(failures, file.path(directory, "fields.csv"), row.names = FALSE)
  summary
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--actual"), optparse::make_option("--expected"), optparse::make_option("--source"),
    optparse::make_option("--contract"), optparse::make_option("--output"))))
  if (any(!vapply(opt[c("actual", "expected", "source", "contract", "output")],
      function(x) length(x) == 1L && nzchar(x), logical(1L)))) stop("all five options are required")
  root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"))
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  for (name in c("actual", "expected", "source")) {
    path <- normalizePath(opt[[name]], mustWork = TRUE)
    if (!endsWith(path, ".parquet")) stop("diagnostic inputs must be Parquet to preserve NULL values")
    DBI::dbExecute(con, paste0("CREATE VIEW ", name, " AS SELECT * FROM read_parquet(",
      DBI::dbQuoteString(con, path), ")"))
  }
  result <- duckvep_fastvep_compare(con, "actual", "expected", "source",
    duckvep_fastvep_fields(opt$contract), opt$output)
  print(result)
  if (!result$passed) quit(status = 1L)
}

if (sys.nframe() == 0L) main()
