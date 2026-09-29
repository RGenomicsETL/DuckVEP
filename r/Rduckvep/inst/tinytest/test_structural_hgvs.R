library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
events <- data.frame(event_index = 0:7, chrom = "21", pos = 7467462, ref = "T",
  alt = c("<DEL>", "<DUP>", "<INV>", "<DEL>", "N[21:7467900[", "<INS>", "<DEL>", "<CNV>"),
  info = c("END=7467465", "END=7467464", "END=7467468", "END=7467465;CIPOS=-1,1",
           "SVTYPE=BND", "END=7467462", "SVTYPE=DEL", "END=7467465"),
  stringsAsFactors = FALSE)
refs <- data.frame(event_index = c(0L, 1L, 2L, 3L),
  reference_sequence = c("TCAG", "TCAGC", "TCAGTAGCAGTAG", "TCAG"))
dbWriteTable(con, "shg_events", events, temporary = TRUE)
dbWriteTable(con, "shg_refs", refs, temporary = TRUE)
rows <- rduckvep_prepare_structural_hgvs(con, "shg_events", "shg_refs")
expect_equal(rows$event_index, 0:7)
expect_identical(rows$hgvs_status,
  c("supported", "supported", "supported", "unsupported", "unsupported", "unsupported",
    "unavailable", "unsupported"))
expect_identical(rows$hgvs_g[1:3], c("21:g.7467463_7467465del", "21:g.7467463_7467464dup",
                                     "21:g.7467463_7467468inv"))
expect_true(all(is.na(rows$hgvs_g[4:8])))
expect_identical(rows$hgvs_reason[4:8], c("imprecise", "breakend", "symbolic_insertion",
                                          "missing_end", "copy_number"))
expect_identical(rows$literal_reference[1], "TCAG")
expect_identical(rows$literal_alternate[1], "T")
expect_identical(rows$normalization[1:3], rep("none", 3L))
expect_identical(rows$transcript_hgvs_route[c(1, 4, 7)], c("literal_equivalent", "unsupported", "unavailable"))
capped <- rduckvep_prepare_structural_hgvs(con, "shg_events", "shg_refs", max_span = 2L)
expect_identical(capped$hgvs_reason[c(1, 3)], c("span_capacity", "span_capacity"))
expect_error(rduckvep_prepare_structural_hgvs(con, "shg_events", "shg_refs", max_span = 0L),
             "between 1 and 60000")
dbDisconnect(con, shutdown = TRUE)
