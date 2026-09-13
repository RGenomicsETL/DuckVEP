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
  cat("CSQ extraction: transport spelling, missing fields and width controls passed\n")
}

main()
