library(DBI)
library(duckdb)

check_stage <- function() {
  root <- "/root/duckvep/data/gnomad-v4.1"
  output <- tempfile("stage-contract-", tmpdir = root)
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  script <- normalizePath("scripts/stage_gnomad_v4.1.R")
  command <- paste("Rscript", shQuote(script), "sv", "fixture", "fixture-md5", shQuote(output))
  stream <- pipe(command, "w")
  writeLines(c(
    "chrY\t123\tfirst\tA\tC,G\t.\tPASS\t.",
    "chr22\t345\tsecond\tA\t<DEL>,*\t.\tq10\tEND=600;SVTYPE=DEL",
    "chr1\t789\tthird\tN\tN[chr2:456[\t.\tPASS\tMATEID=partner",
    "chrUnknown\t10\tfourth\tA\tC\t.\tPASS\t.",
    "chr1\t900\tfifth\tA\tR\t.\tPASS\t."
  ), stream)
  status <- close(stream)
  if (!is.null(status) && status != 0L) stop("streaming staging fixture failed")
  con <- dbConnect(duckdb(shared_home = FALSE))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  path <- as.character(dbQuoteString(con, file.path(output, "part-*.parquet")))
  rows <- dbGetQuery(con, paste0("SELECT record_index, alt_index, chrom, seq_region,
    alternate, status, end_position FROM read_parquet(", path,
    ") ORDER BY record_index, alt_index"))
  observed <- list(record_index = as.integer(rows$record_index),
    alt_index = as.integer(rows$alt_index), seq_region = as.integer(rows$seq_region),
    alternate = rows$alternate, status = rows$status,
    end_position = as.numeric(rows$end_position))
  expected <- list(record_index = c(1L, 1L, 2L, 2L, 3L, 4L, 5L),
    alt_index = c(1L, 2L, 1L, 2L, 1L, 1L, 1L),
    seq_region = c(16L, 16L, 20L, 20L, 13L, NA_integer_, 13L),
    alternate = c("C", "G", "<DEL>", "*", "N[chr2:456[", "C", "R"),
    status = c("literal", "literal", "symbolic", "star", "breakend",
      "unsupported", "unsupported"),
    end_position = c(NA_real_, NA_real_, 600, 600, NA_real_, NA_real_, NA_real_))
  if (!identical(observed, expected)) {
    print(rows)
    stop("physical record, ALT, region or eligibility contract differs")
  }
  counts <- read.delim(file.path(output, "counts.tsv"))
  if (counts$records != 5 || counts$alts != 7 || counts$unsupported != 2) {
    stop("staged input was not conserved")
  }
  cat("staging identity and eligibility: passed\n")
}
check_stage()
