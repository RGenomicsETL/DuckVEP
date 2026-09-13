#!/usr/bin/env Rscript
source("benchmarks/benchmark_duckvep_fastvep_fields.R")

local({
  work <- tempfile("fastvep-source-coverage-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "source.vcf")
  source_map <- file.path(work, "source.parquet")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\t.\tTAA\tTA,T\t.\tPASS\t.",
    "chr1\t20\t.\tTAA\tT,*\t.\tPASS\t.",
    "chr1\t30\tsame_id\tAC\tG,A\t.\tPASS\t.",
    "chr1\t40\tsame_id\tA\tG\t.\tPASS\t.",
    "chr1\t50\t.\tT\tT\t.\tPASS\t."), input)
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con,
    normalizePath("build/release/duckhts.duckdb_extension"))))
  DBI::dbExecute(con, "SET preserve_insertion_order = false")
  duckvep_fastvep_write_source_map(con, input, source_map)
  stopifnot(!DBI::dbGetQuery(con,
    "SELECT current_setting('preserve_insertion_order') AS preserve_order")$preserve_order,
    !DBI::dbExistsTable(con, "fastvep_source"))
  map <- DBI::dbGetQuery(con, paste0("SELECT * REPLACE(record_index::VARCHAR AS record_index,
    alt_index::VARCHAR AS alt_index) FROM read_parquet(", DBI::dbQuoteString(con, source_map),
    ") ORDER BY record_index, alt_index"))
  stopifnot(nrow(map) == 8L, sum(map$eligible) == 6L,
    identical(map$native_allele[map$record_index == "1"], c("A", "-")),
    identical(map$native_allele[map$record_index == "2"], c("-", "*")),
    all(map$native_uploaded[map$record_index == "1"] == "chr1:11-12_AA/A/-"),
    identical(map$eligibility[!map$eligible], c("spanning_deletion", "reference_equal")))

  payload <- function(contract, source_rows = map) {
    fields <- duckvep_fastvep_fields(contract)
    rows <- as.data.frame(setNames(lapply(fields, function(field) rep("", nrow(source_rows))), fields))
    prefix <- if (contract == "operational17") "operational_" else "native_"
    rows$Uploaded_variation <- source_rows[[paste0(prefix, "uploaded")]]
    rows$Allele <- source_rows[[paste0(prefix, "allele")]]
    rows$Feature <- ifelse(source_rows$record_index == "4", "-", "tx1")
    rows$Consequence <- ifelse(source_rows$record_index == "4", "intergenic_variant", "intron_variant")
    if (contract == "vep_csq") {
      # Both physical deletion ALTs can have the same DuckVEP display allele.
      rows$Allele[source_rows$record_index == "1"] <- "-"
      rows <- cbind(source_rows[c("record_index", "alt_index")], rows)
    } else rows$Location <- source_rows[[paste0(prefix, "location")]]
    rows
  }
  checks <- 0L
  check <- function(rows, contract, passed, source = source_map) {
    checks <<- checks + 1L
    DBI::dbWriteTable(con, "output", rows, overwrite = TRUE)
    failures <- file.path(work, paste0("failures", checks, ".parquet"))
    result <- duckvep_fastvep_source_coverage(con, "output", source, contract, failures)
    stopifnot(identical(result$passed, passed), as.numeric(result$output_rows) == nrow(rows))
    n <- DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM read_parquet(",
      DBI::dbQuoteString(con, failures), ")"))$n
    stopifnot(if (passed) n == 0 else n > 0)
    result
  }
  for (contract in c("native_tab17", "operational17", "vep_csq")) {
    rows <- payload(contract)
    second_transcript <- rows[1L, ]
    second_transcript$Feature <- "tx2"
    result <- check(rbind(rows, second_transcript), contract, TRUE)
    stopifnot(result$source_alleles == 6, result$covered_alleles == 6,
      result$emitted_excluded_alleles == 2, result$excluded_output_rows == 2)
    # An intergenic source allele needs final output just like a transcript pair.
    result <- check(rows[map$record_index != "4", ], contract, FALSE)
    stopifnot(result$missing_alleles == 1)
    # Equal row count cannot hide a missing ALT replaced with another ALT's row.
    changed <- rows
    changed[2L, ] <- rows[1L, ]
    stopifnot(check(changed, contract, FALSE)$missing_alleles == 1)
    changed <- rows
    if (contract == "vep_csq") changed$record_index[2L] <- "99" else changed$Allele[2L] <- "UNKNOWN"
    stopifnot(check(changed, contract, FALSE)$unknown_alleles == 1)
  }
  rows <- payload("vep_csq")
  for (key in c("record_index", "alt_index")) for (invalid in
      c("0", "-1", "1.0", "01", "1e0", " 1", "", NA, "18446744073709551616")) {
    changed <- rows
    changed[[key]][1L] <- invalid
    stopifnot(check(changed, "vep_csq", FALSE)$unknown_alleles == 1)
  }
  DBI::dbWriteTable(con, "output", rows, overwrite = TRUE)
  plan <- DBI::dbGetQuery(con, paste0("EXPLAIN SELECT * FROM output o LEFT JOIN read_parquet(",
    DBI::dbQuoteString(con, source_map), ") s ON ",
    duckvep_fastvep_ordinal_predicate(c("record_index", "alt_index"))))$explain_value
  stopifnot(any(grepl("HASH_JOIN", plan, fixed = TRUE)),
    !any(grepl("BLOCKWISE_NL_JOIN", plan, fixed = TRUE)))
  # A display-key collision cannot be resolved by choosing one source record.
  duplicate_map <- rbind(map, transform(map[1L, ], record_index = "6"))
  DBI::dbWriteTable(con, "duplicate_source", duplicate_map)
  ambiguous_map <- file.path(work, "ambiguous.parquet")
  DBI::dbExecute(con, paste0("COPY duplicate_source TO ", DBI::dbQuoteString(con, ambiguous_map),
    " (FORMAT PARQUET)"))
  stopifnot(check(payload("native_tab17"), "native_tab17", FALSE, ambiguous_map)$ambiguous_alleles == 1)
  for (kind in c("duplicate_ordinal", "negative_ordinal", "false_eligibility")) {
    invalid_map <- map
    if (kind == "duplicate_ordinal") invalid_map <- rbind(map, map[1L, ])
    if (kind == "negative_ordinal") invalid_map$alt_index[1L] <- "-1"
    if (kind == "false_eligibility") invalid_map$eligible[1L] <- FALSE
    DBI::dbWriteTable(con, "invalid_source", invalid_map, overwrite = TRUE)
    path <- file.path(work, paste0(kind, ".parquet"))
    DBI::dbExecute(con, paste0("COPY invalid_source TO ", DBI::dbQuoteString(con, path), " (FORMAT PARQUET)"))
    error <- tryCatch(check(payload("native_tab17"), "native_tab17", FALSE, path), error = identity)
    stopifnot(inherits(error, "error"))
  }

  # Exercise the publication-facing receipt against a physical final file.
  output <- file.path(work, "output.tsv")
  receipt <- file.path(work, "receipt.csv")
  coverage <- file.path(work, "coverage.csv")
  input_hash <- strsplit(system2("sha256sum", shQuote(input), stdout = TRUE), " +")[[1L]][1L]
  arguments <- c("benchmarks/benchmark_duckvep_fastvep_receipt.R", "--input", output,
    "--tool", "duckvep", "--output-contract", "vep_csq", "--output", receipt,
    "--source-map", source_map, "--source-sha256", input_hash, "--coverage-output", coverage)
  write.table(rows, output, sep = "\t", quote = FALSE, row.names = FALSE)
  log <- file.path(work, "receipt.log")
  stopifnot(system2("Rscript", shQuote(arguments), stdout = log, stderr = log) == 0L)
  result <- read.csv(coverage, colClasses = "character")
  observed <- read.csv(receipt, colClasses = "character")
  stopifnot(result$passed == "TRUE", result$scope == "final_output",
    result$input_sha256 == input_hash, result$output_sha256 == observed$sha256,
    result$source_alleles == "6", result$covered_alleles == "6")
  write.table(rows[-1L, ], output, sep = "\t", quote = FALSE, row.names = FALSE)
  stopifnot(system2("Rscript", shQuote(arguments), stdout = log, stderr = log) != 0L)
  result <- read.csv(coverage, colClasses = "character")
  stopifnot(result$passed == "FALSE", result$missing_alleles == "1")

  # Run the real flattening CLI with its complete 32-field CSQ payload. These
  # transport values are explicit test inputs, not biological expectations.
  fields <- duckvep_fastvep_fields("vep_csq")
  payload_fields <- setdiff(fields, "Uploaded_variation")
  fast <- payload("vep_csq")[fields]
  fast$Allele <- map$native_allele
  records <- vapply(unique(map$record_index), function(index) {
    selected <- which(map$record_index == index)
    first <- selected[[1L]]
    annotations <- apply(fast[selected, payload_fields, drop = FALSE], 1L, paste, collapse = "|")
    id <- if (is.na(map$variant_id[first])) "." else map$variant_id[first]
    paste(c(map$chrom[first], map$position[first], id, map$reference[first],
      map$raw_alternates[first], ".", "PASS", paste0("CSQ=", paste(annotations, collapse = ","))),
      collapse = "\t")
  }, character(1L))
  annotated <- file.path(work, "annotated.vcf")
  writeLines(c("##fileformat=VCFv4.2", paste0(
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ',
    paste(payload_fields, collapse = "|"), '">'),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", rev(records)), annotated)
  extract_args <- c("benchmarks/benchmark_duckvep_fastvep_extract.R", "--input", annotated,
    "--source-map", source_map, "--output", output, "--threads", "4")
  stopifnot(system2("Rscript", shQuote(extract_args), stdout = log, stderr = log) == 0L,
    duckvep_fastvep_tab_header(output, duckvep_fastvep_transport_fields("vep_csq")) == 0L,
    system2("Rscript", shQuote(arguments), stdout = log, stderr = log) == 0L)
  result <- read.csv(coverage, colClasses = "character")
  stopifnot(result$passed == "TRUE", result$source_alleles == "6",
    result$emitted_excluded_alleles == "2", result$output_rows == "8")
  cat("Final source-ALT coverage:", checks, "mapping controls and receipt pass/fail passed\n")
})
