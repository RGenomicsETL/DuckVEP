#!/usr/bin/env Rscript

# Read the declared VCF transport without CSQ type inference or unescaping.
duckvep_fastvep_vcf_header <- function(path) {
  connection <- if (grepl("\\.(gz|bgz)$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection), add = TRUE)
  csq <- character()
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) stop("VCF has no #CHROM header")
    if (startsWith(line, "##INFO=<ID=CSQ,")) csq <- c(csq, line)
    if (startsWith(line, "#CHROM\t")) break
    if (!startsWith(line, "#")) stop("VCF record precedes #CHROM header")
  }
  columns <- strsplit(sub("^#", "", line), "\t", fixed = TRUE)[[1L]]
  if (anyDuplicated(columns) || !identical(head(columns, 8L),
      c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))) {
    stop("invalid VCF columns")
  }
  fields <- character()
  if (length(csq)) {
    if (length(csq) != 1L || !grepl("Format: [^\"]+\"", csq)) stop("ambiguous CSQ schema")
    fields <- strsplit(sub('.*Format: ([^"]+)".*', "\\1", csq), "|", fixed = TRUE)[[1L]]
    if (anyDuplicated(fields) || any(!nzchar(fields))) stop("invalid CSQ fields")
  }
  list(columns = columns, fields = fields)
}

duckvep_fastvep_vcf_relation <- function(con, path, header) {
  qs <- function(x) as.character(DBI::dbQuoteString(con, x))
  schema <- paste(paste0(qs(header$columns), ": 'VARCHAR'"), collapse = ", ")
  paste0("read_csv(", qs(normalizePath(path, mustWork = TRUE)),
    ", delim = '\t', header = false, comment = '#', quote = '', escape = '',",
    " force_not_null = [", paste(qs(header$columns), collapse = ", "),
    "], columns = {", schema, "}, auto_detect = false)")
}

duckvep_fastvep_extract_csq <- function(con, input, table, fields, source_map = "") {
  header <- duckvep_fastvep_vcf_header(input)
  payload <- setdiff(fields, "Uploaded_variation")
  absent <- setdiff(payload, header$fields)
  if (length(absent)) stop("CSQ schema does not emit: ", paste(absent, collapse = ", "))
  qi <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  raw <- duckvep_fastvep_vcf_relation(con, input, header)
  if (DBI::dbExistsTable(con, table)) stop("extraction table already exists")
  expressions <- c(Uploaded_variation = "ID")
  for (field in payload) expressions[[field]] <- paste0("values[", match(field, header$fields), "]")
  values <- paste(paste0(expressions[fields], " AS ", qi(fields)), collapse = ", ")
  identity <- ""
  source_join <- ""
  if (nzchar(source_map)) {
    if (!"Allele" %in% fields) stop("identity extraction requires the CSQ Allele field")
    q <- function(x) as.character(DBI::dbQuoteString(con, x))
    DBI::dbExecute(con, paste0("CREATE TEMP VIEW fastvep_extract_source AS
      SELECT chrom, position, coalesce(variant_id, '.') AS variant_id,
        reference, raw_alternates, native_allele, count(*) AS source_matches,
        min(record_index) AS record_index, min(alt_index) AS alt_index
      FROM read_parquet(", q(source_map), ") GROUP BY ALL"))
    on.exit(DBI::dbExecute(con, "DROP VIEW fastvep_extract_source"), add = TRUE)
    identity <- "s.record_index, s.alt_index, coalesce(s.source_matches, 0) AS _source_matches, "
    source_join <- paste0(" LEFT JOIN fastvep_extract_source s
      ON a.CHROM = s.chrom AND try_cast(a.POS AS UBIGINT) = s.position
        AND a.ID = s.variant_id AND a.REF = s.reference AND a.ALT = s.raw_alternates
        AND a.values[", match("Allele", header$fields), "] = s.native_allele")
  }
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", qi(table), " AS WITH records AS (
    SELECT *, regexp_extract(INFO, '(?:^|;)CSQ=([^;]*)', 1) AS csq,
      length(regexp_extract_all(INFO, '(?:^|;)CSQ=')) AS csq_keys FROM ", raw,
    "), annotations AS (
      SELECT CHROM, POS, ID, REF, ALT, csq_keys, string_split(annotation, '|') AS values FROM records,
        UNNEST(string_split(csq, ',')) a(annotation)
    ) SELECT ", identity, values, ", length(values) AS _csq_width, csq_keys AS _csq_keys
      FROM annotations a", source_join))
  width <- DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM ", qi(table),
    " WHERE _csq_keys != 1 OR _csq_width != ", length(header$fields)))$n
  if (width != 0) stop("missing/duplicate CSQ or field count differs from its header: ", width, " rows")
  if (nzchar(source_map)) {
    invalid <- DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM ", qi(table),
      " WHERE _source_matches != 1"))$n
    if (invalid != 0) stop("CSQ output has unknown or ambiguous physical ALT identity: ", invalid, " rows")
    DBI::dbExecute(con, paste0("ALTER TABLE ", qi(table), " DROP COLUMN _source_matches"))
  }
  DBI::dbExecute(con, paste0("ALTER TABLE ", qi(table), " DROP COLUMN _csq_width"))
  DBI::dbExecute(con, paste0("ALTER TABLE ", qi(table), " DROP COLUMN _csq_keys"))
  invisible(table)
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--input"), optparse::make_option("--output"),
    optparse::make_option("--source-map", dest = "source_map", default = ""),
    optparse::make_option("--memory-limit", dest = "memory_limit", default = "4GB"),
    optparse::make_option("--max-spill", dest = "max_spill", default = "8GB"),
    optparse::make_option("--threads", type = "integer", default = 1L))))
  if (is.null(opt$input) || is.null(opt$output) || opt$threads < 1L) stop("input, output and positive threads required")
  root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"))
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, paste0("SET threads=", opt$threads))
  DBI::dbExecute(con, paste("SET memory_limit =", DBI::dbQuoteString(con, opt$memory_limit)))
  DBI::dbExecute(con, paste("SET max_temp_directory_size =", DBI::dbQuoteString(con, opt$max_spill)))
  DBI::dbExecute(con, paste("SET temp_directory =", DBI::dbQuoteString(con, file.path(tempdir(), "duckdb"))))
  duckvep_fastvep_extract_csq(con, opt$input, "csq", duckvep_fastvep_fields("vep_csq"), opt$source_map)
  DBI::dbExecute(con, paste0("COPY csq TO ", DBI::dbQuoteString(con, opt$output),
    " (FORMAT CSV, DELIMITER '\t', HEADER TRUE, QUOTE '', ESCAPE '')"))
}

if (sys.nframe() == 0L) main()
