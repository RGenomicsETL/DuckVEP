#!/usr/bin/env Rscript
source("benchmarks/benchmark_duckvep_fastvep_extract.R")
source("benchmarks/benchmark_duckvep_fastvep_fields.R")

main <- function() {
  directory <- tempfile("fastvep-extract-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- file.path(directory, "input.vcf")
  header <- c("##fileformat=VCFv4.2",
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature|Codons|HGVSp">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tX=3;CSQ=C|tx1|-|p.%3D,C|tx2||",
    "1\t4\te2\tG\tT\t.\tPASS\tCSQ=T|tx1|A/B|p.(a&b);X=4"), input)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fields <- c("Uploaded_variation", "Allele", "Feature", "Codons", "HGVSp")
  duckvep_fastvep_extract_csq(con, input, "actual", fields)
  actual <- DBI::dbGetQuery(con, "SELECT * FROM actual ORDER BY Uploaded_variation, Feature")
  expected <- data.frame(Uploaded_variation = c("e1", "e1", "e2"), Allele = c("C", "C", "T"),
    Feature = c("tx1", "tx2", "tx1"), Codons = c("-", "", "A/B"),
    HGVSp = c("p.%3D", "", "p.(a&b)"))
  stopifnot(identical(actual, expected))
  correlated_query <- function(connection, path, fields) {
    declared <- duckvep_fastvep_vcf_header(path)
    raw <- duckvep_fastvep_vcf_relation(connection, path, declared)
    payload <- setdiff(fields, "Uploaded_variation")
    columns <- paste0("values[", match(payload, declared$fields), "] AS ",
      DBI::dbQuoteIdentifier(connection, payload), collapse = ", ")
    paste0("WITH records AS (SELECT ID, regexp_extract(INFO, '(?:^|;)CSQ=([^;]*)', 1) AS csq FROM ", raw,
      "), annotations AS (SELECT ID, string_split(annotation, '|') AS values
        FROM records, UNNEST(string_split(csq, ',')) a(annotation))
      SELECT ID AS Uploaded_variation, ", columns, " FROM annotations")
  }
  DBI::dbExecute(con, paste("CREATE TEMP TABLE correlated AS", correlated_query(con, input, fields)))
  stopifnot(DBI::dbGetQuery(con, "SELECT count(*) n FROM (
    (SELECT * FROM actual EXCEPT ALL SELECT * FROM correlated)
    UNION ALL (SELECT * FROM correlated EXCEPT ALL SELECT * FROM actual))")$n == 0L)
  expect_error <- function(expression) stopifnot(inherits(tryCatch({
    force(expression)
    NULL
  }, error = identity), "error"))
  expect_error(duckvep_fastvep_extract_csq(con, input, "unsupported", c(fields, "GENCODE_PRIMARY")))
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tCSQ=C|tx1|-"), input)
  expect_error(duckvep_fastvep_extract_csq(con, input, "short", fields))
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tX=1"), input)
  expect_error(duckvep_fastvep_extract_csq(con, input, "missing", fields))
  writeLines(c(header, "1\t3\te1\tA\tC\t.\tPASS\tCSQ=C|tx1||;CSQ=C|tx2||"), input)
  expect_error(duckvep_fastvep_extract_csq(con, input, "duplicate", fields))
  tab <- file.path(directory, "output.tsv")
  native <- duckvep_fastvep_fields("native_tab17")
  writeLines(c("## fastvep version=0.3.0", paste0("#", paste(native, collapse = "\t"))), tab)
  stopifnot(duckvep_fastvep_tab_header(tab, native) == 1L)
  writeLines(paste(native, collapse = "\t"), tab)
  stopifnot(duckvep_fastvep_tab_header(tab, native) == 0L)
  renamed <- native
  renamed[native == "Gene"] <- "NotGene"
  duplicate <- native
  duplicate[native == "FLAGS"] <- "Gene"
  for (invalid in list(renamed, native[-length(native)], c(native, "Extra"), duplicate, rev(native))) {
    writeLines(paste(invalid, collapse = "\t"), tab)
    expect_error(duckvep_fastvep_tab_header(tab, native))
  }
  writeLines(character(), tab)
  expect_error(duckvep_fastvep_tab_header(tab, native))
  # The conformance driver imports helpers inside main, not into the global
  # environment. Its reader must resolve the same locally imported header check.
  helpers <- new.env(parent = baseenv())
  source("benchmarks/benchmark_duckvep_fastvep_fields.R", local = helpers)
  writeLines(c("Feature\tCodons\tHGVSp", "tx1\t\t-"), tab)
  helpers$duckvep_fastvep_read_field_tab(con, tab, "local_tab", c("Feature", "Codons", "HGVSp"))
  stopifnot(identical(DBI::dbGetQuery(con, "SELECT * FROM local_tab"),
    data.frame(Feature = "tx1", Codons = "", HGVSp = "-")))

  # Physical geometry and the upstream native allele identify both ALTs even
  # when IDs are missing. Full ALT lists and known star output remain intact.
  source_map <- file.path(directory, "source.parquet")
  source <- data.frame(record_index = c(1L, 1L, 2L, 2L), alt_index = c(1L, 2L, 1L, 2L),
    chrom = "1", position = c(10L, 10L, 20L, 20L), variant_id = NA_character_,
    reference = "TAA", raw_alternates = c("TA,T", "TA,T", "T,*", "T,*"),
    native_allele = c("A", "-", "-", "*"))
  DBI::dbWriteTable(con, "source_map", source)
  DBI::dbExecute(con, paste0("COPY source_map TO ", DBI::dbQuoteString(con, source_map), " (FORMAT PARQUET)"))
  identity_input <- file.path(directory, "identity.vcf")
  identity_records <- c(
    "1\t20\t.\tTAA\tT,*\t.\tPASS\tCSQ=-|-||,*|-||",
    "1\t10\t.\tTAA\tTA,T\t.\tPASS\tCSQ=A|tx1||p.%3D,-|tx1||,A|tx2||")
  writeLines(c(header, identity_records), identity_input)
  duckvep_fastvep_extract_csq(con, identity_input, "plain", fields)
  duckvep_fastvep_extract_csq(con, identity_input, "identified", fields, source_map)
  observed <- DBI::dbGetQuery(con, "SELECT * FROM identified ORDER BY record_index, alt_index, Feature")
  stopifnot(identical(names(observed), c(duckvep_fastvep_identity_fields, fields)),
    identical(as.integer(observed$record_index), c(1L, 1L, 1L, 2L, 2L)),
    identical(as.integer(observed$alt_index), c(1L, 1L, 2L, 1L, 2L)),
    identical(observed$Allele, c("A", "A", "-", "-", "*")),
    all(observed$Uploaded_variation == "."), all(observed$Feature[4:5] == "-"))
  payload <- paste(DBI::dbQuoteIdentifier(con, fields), collapse = ", ")
  differences <- DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM (
    (SELECT ", payload, " FROM identified EXCEPT ALL SELECT * FROM plain)
    UNION ALL (SELECT * FROM plain EXCEPT ALL SELECT ", payload, " FROM identified))"))$n
  stopifnot(differences == 0)
  writeLines(c(header, sub("\t20\t", "\t21\t", identity_records)), identity_input)
  expect_error(duckvep_fastvep_extract_csq(con, identity_input, "unknown_identity", fields, source_map))
  writeLines(c(header, identity_records), identity_input)
  duplicate_source <- rbind(source, transform(source[1L, ], record_index = 3L))
  DBI::dbWriteTable(con, "duplicate_source", duplicate_source)
  duplicate_map <- file.path(directory, "duplicate_source.parquet")
  DBI::dbExecute(con, paste0("COPY duplicate_source TO ", DBI::dbQuoteString(con, duplicate_map), " (FORMAT PARQUET)"))
  expect_error(duckvep_fastvep_extract_csq(con, identity_input, "ambiguous_identity", fields, duplicate_map))
  # A correlated CSQ split retains the full CSQ value as a delimiter-join key.
  # Both queries face the same memory cap and cannot spill. Duplicate annotations
  # are intentional: each physical record must retain all 1,024 copies.
  bounded_input <- file.path(directory, "multiplicity.vcf")
  wide_fields <- duckvep_fastvep_fields("vep_csq")
  payload <- setdiff(wide_fields, "Uploaded_variation")
  bounded_header <- c(header[[1L]], paste0('##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ',
    paste(payload, collapse = "|"), '">'), header[[3L]])
  records <- vapply(seq_len(128L), function(i) {
    values <- setNames(rep("", length(payload)), payload)
    values[c("Allele", "Feature", "Codons", "HGVSp")] <- c("C", paste0("tx", i), "-", "p.%3D")
    paste("1", i, paste0("e", i), "A", "C", ".", "PASS",
      paste0("CSQ=", paste(rep(paste(values, collapse = "|"), 1024L), collapse = ",")), sep = "\t")
  }, character(1L))
  writeLines(c(bounded_header, records), bounded_input)
  limited <- function() DBI::dbConnect(duckdb::duckdb(config = list(
    threads = "1", memory_limit = "64MB", max_temp_directory_size = "0B")))
  check_correlated <- function() {
    connection <- limited()
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    query <- correlated_query(connection, bounded_input, wide_fields)
    error <- tryCatch({
      DBI::dbExecute(connection, paste("CREATE TEMP TABLE correlated AS", query))
      NULL
    }, error = identity)
    stopifnot(inherits(error, "error"), grepl("Out of Memory", conditionMessage(error), fixed = TRUE))
  }
  check_projected <- function() {
    connection <- limited()
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    duckvep_fastvep_extract_csq(connection, bounded_input, "bounded", wide_fields)
    empty_fields <- setdiff(payload, c("Allele", "Feature", "Codons", "HGVSp"))
    empty <- paste0(DBI::dbQuoteIdentifier(connection, empty_fields), " = ''", collapse = " AND ")
    counts <- DBI::dbGetQuery(connection, paste0("SELECT count(*) AS records, sum(copies) AS annotations,
      count(*) FILTER (WHERE copies != 1024 OR NOT exact_fields) AS invalid FROM (
      SELECT Uploaded_variation, count(*) AS copies,
        bool_and(Feature = 'tx' || substr(Uploaded_variation, 2) AND Allele = 'C'
          AND Codons = '-' AND HGVSp = 'p.%3D' AND ", empty, ") AS exact_fields
      FROM bounded GROUP BY Uploaded_variation)"))
    stopifnot(counts$records == 128L, counts$annotations == 131072L, counts$invalid == 0L)
  }
  check_correlated()
  check_projected()
  # Exercise the file-writing process, not only its materializing R helper.
  # Malformed rows follow valid rows so a filtering implementation cannot pass.
  check_cli <- function(label, path, map = "", expected_rows = NULL, expected_error = NULL) {
    api_table <- paste0("api_", label)
    api_error <- tryCatch({
      duckvep_fastvep_extract_csq(con, path, api_table, wide_fields, map)
      NULL
    }, error = identity)
    output <- file.path(directory, paste0(label, ".tsv"))
    log <- file.path(directory, paste0(label, ".log"))
    args <- c("benchmarks/benchmark_duckvep_fastvep_extract.R", "--input", path,
      "--output", output, "--memory-limit", "64MB", "--max-spill", "0B", "--threads", "1")
    if (nzchar(map)) args <- c(args, "--source-map", map)
    status <- suppressWarnings(system2("Rscript", shQuote(args), stdout = log, stderr = log))
    if (!is.null(expected_error)) {
      stopifnot(inherits(api_error, "error"), grepl(expected_error, conditionMessage(api_error), fixed = TRUE),
        status != 0L, any(grepl(expected_error, readLines(log), fixed = TRUE)))
      return(invisible(NULL))
    }
    stopifnot(is.null(api_error), status == 0L)
    cli_table <- paste0("cli_", label)
    columns <- c(if (nzchar(map)) duckvep_fastvep_identity_fields, wide_fields)
    duckvep_fastvep_read_field_tab(con, output, cli_table, columns)
    differences <- DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM (
      (SELECT * FROM ", api_table, " EXCEPT ALL SELECT * FROM ", cli_table, ")
      UNION ALL (SELECT * FROM ", cli_table, " EXCEPT ALL SELECT * FROM ", api_table, "))"))$n
    stopifnot(differences == 0L, DBI::dbGetQuery(con, paste("SELECT count(*) n FROM", cli_table))$n == expected_rows)
  }
  check_cli("bounded", bounded_input, expected_rows = 131072L)
  annotation <- function(allele, feature, codons = "", hgvsp = "") {
    values <- setNames(rep("", length(payload)), payload)
    values[c("Allele", "Feature", "Codons", "HGVSp")] <- c(allele, feature, codons, hgvsp)
    paste(values, collapse = "|")
  }
  a <- annotation("A", "tx1", "-", "p.%3D")
  full_records <- c(
    paste0("1\t20\t.\tTAA\tT,*\t.\tPASS\tCSQ=", annotation("-", "-"), ",", annotation("*", "-")),
    paste0("1\t10\t.\tTAA\tTA,T\t.\tPASS\tCSQ=", a, ",", annotation("-", "tx1"), ",", a))
  cli_input <- file.path(directory, "cli.vcf")
  writeLines(c(bounded_header, full_records), cli_input)
  check_cli("identity", cli_input, source_map, expected_rows = 5L)
  check_cli("ambiguous", cli_input, duplicate_map,
    expected_error = "CSQ output has unknown or ambiguous physical ALT identity")
  writeLines(c(bounded_header, sub("\t20\t", "\t21\t", full_records)), cli_input)
  check_cli("unknown", cli_input, source_map,
    expected_error = "CSQ output has unknown or ambiguous physical ALT identity")
  invalid_info <- c(missing_key = "X=1", duplicate_key = paste0("CSQ=", a, ";CSQ=", a),
    short_width = paste0("CSQ=", substr(a, 1L, nchar(a) - 1L)), long_width = paste0("CSQ=", a, "|"))
  for (label in names(invalid_info)) {
    invalid <- paste0("1\t10\t.\tTAA\tTA,T\t.\tPASS\t", invalid_info[[label]])
    writeLines(c(bounded_header, full_records[[1L]], invalid), cli_input)
    check_cli(label, cli_input, expected_error = "missing/duplicate CSQ or field count differs from its header")
  }
  cat("CSQ extraction: transport spelling, missing fields and width controls passed\n")
  cat("CSQ extraction: exact lateral/projection multiset; 131,072 rows at 64MB with zero spill; correlated negative control failed OOM\n")
  cat("CSQ extraction: direct CLI COPY matches materialized rows; API and CLI reject all malformed/identity controls\n")
}

main()
