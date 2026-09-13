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
  compare <- function(actual, label, oracle = expected) {
    DBI::dbWriteTable(con, "actual", actual, overwrite = TRUE)
    DBI::dbWriteTable(con, "expected", oracle, overwrite = TRUE)
    duckvep_fastvep_compare(con, "actual", "expected", c("Feature", "payload"),
      file.path(directory, label))
  }
  result <- compare(expected[6:1, ], "equal")
  stopifnot(result$passed, result$union_keys == 6, result$compared_keys == 6)
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
  result <- compare(expected[FALSE, ], "empty", expected[FALSE, ])
  stopifnot(result$passed, result$union_keys == 0, result$compared_keys == 0)
  cat("Exact field comparator: all corruption controls passed\n")
}

main()
