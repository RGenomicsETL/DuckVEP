#!/usr/bin/env Rscript
source("benchmarks/benchmark_duckvep_fastvep_compare.R")

main <- function() {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  directory <- tempfile("fastvep-compare-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  expected <- data.frame(record_index = as.character(1:6), alt_index = "1",
    Feature = c(rep("ENST1", 5), ""), payload = c(NA, "", "-", "a&b", "%3D", ""),
    stringsAsFactors = FALSE)
  DBI::dbWriteTable(con, "expected", expected)
  source <- data.frame(record_index = as.character(1:6), alt_index = "1")
  compare <- function(actual, label, oracle = expected, declared = source) {
    DBI::dbWriteTable(con, "actual", actual, overwrite = TRUE)
    DBI::dbWriteTable(con, "expected", oracle, overwrite = TRUE)
    DBI::dbWriteTable(con, "source", declared, overwrite = TRUE)
    duckvep_fastvep_compare(con, "actual", "expected", "source", c("Feature", "payload"),
      file.path(directory, label))
  }
  result <- compare(expected[6:1, ], "equal")
  stopifnot(result$passed, result$union_keys == 6, result$compared_keys == 6,
    result$source_rows == 6, result$source_alleles == 6)
  for (i in 1:5) {
    changed <- expected
    changed$payload[[i]] <- c("", "-", NA, "b&a", "=")[[i]]
    result <- compare(changed, paste0("field", i))
    stopifnot(!result$passed, result$field_failures == 1, result$union_keys == 6)
  }
  result <- compare(expected[-1, ], "missing")
  stopifnot(!result$passed, result$missing_keys == 1, result$union_keys == 6)
  extra <- expected[1, ]
  extra$record_index <- "7"
  result <- compare(rbind(expected, extra), "extra")
  stopifnot(!result$passed, result$extra_keys == 1, result$union_keys == 7)
  duplicate <- rbind(expected, expected[1, ])
  result <- compare(duplicate, "actual_duplicate")
  stopifnot(!result$passed, result$actual_rows == 7, result$actual_duplicate_keys == 1,
    result$compared_keys == 5)
  result <- compare(expected, "expected_duplicate", duplicate)
  stopifnot(!result$passed, result$expected_rows == 7, result$expected_duplicate_keys == 1)
  result <- compare(duplicate, "both_duplicate", duplicate)
  stopifnot(!result$passed, result$actual_duplicate_keys == 1, result$expected_duplicate_keys == 1)
  retained <- DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
    DBI::dbQuoteString(con, file.path(directory, "both_duplicate/duplicate_rows.parquet")), ")"))
  stopifnot(nrow(retained) == 4L, sum(retained$side == "actual") == 2L)
  missing_key <- expected
  missing_key$record_index[[1L]] <- NA_character_
  result <- compare(missing_key, "invalid", missing_key)
  stopifnot(!result$passed, result$invalid_keys == 1)
  other_alt <- expected
  other_alt$alt_index[[1L]] <- "2"
  result <- compare(other_alt, "alt_identity")
  stopifnot(!result$passed, result$missing_keys == 1, result$extra_keys == 1)
  result <- compare(expected[FALSE, ], "empty", expected[FALSE, ], source[FALSE, ])
  stopifnot(result$passed, result$union_keys == 0, result$compared_keys == 0)
  result <- compare(expected[-1, ], "shared_omission", expected[-1, ])
  stopifnot(!result$passed, result$union_keys == 5, result$compared_keys == 5,
    result$missing_keys == 0, result$extra_keys == 0,
    result$source_alleles == 6, result$actual_missing_source_alleles == 1,
    result$expected_missing_source_alleles == 1)
  retained <- DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
    DBI::dbQuoteString(con, file.path(directory, "shared_omission/source_coverage_failures.parquet")), ")"))
  stopifnot(nrow(retained) == 1L, retained$record_index == "1", retained$alt_index == "1",
    retained$source_rows == 1, retained$actual_rows == 0, retained$expected_rows == 0)
  result <- compare(expected[-6, ], "shared_intergenic_omission", expected[-6, ])
  stopifnot(!result$passed, result$actual_missing_source_alleles == 1,
    result$expected_missing_source_alleles == 1)
  result <- compare(expected[FALSE, ], "both_empty", expected[FALSE, ])
  stopifnot(!result$passed, result$union_keys == 0, result$source_alleles == 6,
    result$actual_missing_source_alleles == 6, result$expected_missing_source_alleles == 6)
  result <- compare(rbind(expected, extra), "shared_unexpected", rbind(expected, extra))
  stopifnot(!result$passed, result$union_keys == 7, result$missing_keys == 0,
    result$extra_keys == 0, result$actual_unknown_alleles == 1, result$expected_unknown_alleles == 1)
  result <- compare(expected, "empty_source", expected, source[FALSE, ])
  stopifnot(!result$passed, result$source_alleles == 0,
    result$actual_unknown_alleles == 6, result$expected_unknown_alleles == 6)
  result <- compare(expected, "duplicate_source", declared = rbind(source, source[1L, ]))
  stopifnot(!result$passed, result$source_rows == 7, result$source_alleles == 6,
    result$source_duplicate_alleles == 1)
  for (axis in c("record_index", "alt_index")) {
    invalid <- c(NA_character_, "0", "-1", "1.5", "18446744073709551616", "not_an_ordinal")
    for (i in seq_along(invalid)) {
      declared <- source
      declared[[axis]][[1L]] <- invalid[[i]]
      result <- compare(expected, paste0("invalid_source_", axis, "_", i), declared = declared)
      stopifnot(!result$passed, result$source_invalid_alleles == 1)
    }
  }
  another_transcript <- expected[1L, ]
  another_transcript$Feature <- "ENST2"
  result <- compare(rbind(expected, another_transcript), "transcript_multiplicity",
    rbind(expected, another_transcript))
  stopifnot(result$passed, result$union_keys == 7, result$source_alleles == 6)
  source("benchmarks/benchmark_duckvep_fastvep_fields.R", local = TRUE)
  cli_rows <- source
  for (field in duckvep_fastvep_fields("native_tab17")) cli_rows[[field]] <- ""
  cli_rows$Feature <- expected$Feature
  DBI::dbWriteTable(con, "cli_rows", cli_rows)
  output_path <- file.path(directory, "output.parquet")
  source_path <- file.path(directory, "source.parquet")
  DBI::dbExecute(con, paste0("COPY cli_rows TO ", DBI::dbQuoteString(con, output_path), " (FORMAT PARQUET)"))
  DBI::dbExecute(con, paste0('COPY "source" TO ', DBI::dbQuoteString(con, source_path), " (FORMAT PARQUET)"))
  arguments <- c("benchmarks/benchmark_duckvep_fastvep_compare.R", "--actual", output_path,
    "--expected", output_path, "--contract", "native_tab17", "--output", file.path(directory, "cli"))
  status <- system2("Rscript", shQuote(c(arguments, "--source", source_path)),
    stdout = file.path(directory, "cli.log"), stderr = file.path(directory, "cli.stderr"))
  stopifnot(status == 0L)
  status <- system2("Rscript", shQuote(arguments),
    stdout = file.path(directory, "cli_missing_source.log"), stderr = file.path(directory, "cli_missing_source.stderr"))
  stopifnot(status != 0L,
    any(grepl("all five options are required", readLines(file.path(directory, "cli_missing_source.stderr")), fixed = TRUE)))
  cat("Exact field comparator: all corruption controls passed\n")
}

main()
