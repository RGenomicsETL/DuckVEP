library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
# One plus-strand transcript: exons 101-200, 215-314, 330-429, 446-545, 566-665 with the stop codon
# at 600, introns of 14, 15, 16 and 20 nt. Variant rows carry the columns of the projected annotation.
dbExecute(con, "CREATE TABLE lof_ref AS SELECT '1' AS chrom, 0::BIGINT AS \"start\", 1000::BIGINT AS \"end\",
  repeat('C', 200) || 'GT' || repeat('T', 10) || 'AG' || repeat('C', 100) || 'GT' || repeat('T', 8) || 'AGTAG' ||
  repeat('C', 100) || 'AT' || repeat('T', 12) || 'AC' || repeat('C', 100) || 'GC' || repeat('T', 16) || 'AG' ||
  repeat('C', 435) AS seq")
dbExecute(con, "CREATE TABLE lof_tx AS SELECT 0 AS transcript_index, '1' AS seq_region_name, 1 AS strand,
  220 AS cds_start, 600 AS cds_end, 'protein_coding' AS transcript_biotype,
  [{'exon_start': 101, 'exon_end': 200}, {'exon_start': 215, 'exon_end': 314}, {'exon_start': 330, 'exon_end': 429},
   {'exon_start': 446, 'exon_end': 545}, {'exon_start': 566, 'exon_end': 665}] AS exons")
dbExecute(con, "CREATE TABLE lof_ann AS SELECT event_index::UBIGINT AS event_index, 0::UINTEGER AS transcript_index,
  consequence, ref AS reference, alt AS alternate,
  struct_pack(feature_start0 := (vs - 1)::UBIGINT, feature_end0 := ve::UBIGINT, reference_difference_offset := 0::USMALLINT,
    reference_difference_length := length(ref)::USMALLINT, alternate_difference_offset := 0::USMALLINT,
    alternate_difference_length := length(alt)::USMALLINT) AS geometry,
  exon_first::UINTEGER AS exon_first, exon_first::UINTEGER AS exon_last, 5::UINTEGER AS exon_total,
  intron_first::UINTEGER AS intron_first, intron_first::UINTEGER AS intron_last, 4::UINTEGER AS intron_total,
  cds_end::UINTEGER AS cds_start, cds_end::UINTEGER AS cds_end, false AS interbase,
  false AS cds_start_nf, false AS cds_end_nf
  FROM (VALUES
  (1, 'frameshift_variant', 'C', '', 495, 495, 4, NULL, 245),
  (2, 'frameshift_variant', 'C', '', 494, 494, 4, NULL, 244),
  (3, 'splice_donor_variant', 'G', 'A', 201, 201, NULL, 1, NULL),
  (4, 'splice_acceptor_variant', 'G', 'A', 329, 329, NULL, 2, NULL),
  (5, 'splice_donor_variant', 'A', 'G', 430, 430, NULL, 3, NULL),
  (6, 'missense_variant', 'C', 'A', 400, 400, 3, NULL, 166)
  ) t(event_index, consequence, ref, alt, vs, ve, exon_first, intron_first, cds_end)")

sql <- rduckvep_lof_sql(con, "lof_ann", "lof_tx", "lof_ref")
rows <- dbGetQuery(con, paste0("SELECT * FROM query(", as.character(dbQuoteString(con, sql)),
                               ") ORDER BY event_index"))
expect_identical(rows$lof, c("LC", "HC", "LC", "HC", "HC", NA))
expect_identical(rows$lof_filter, c("END_TRUNC", NA, "SMALL_INTRON,5UTR_SPLICE", NA, NA, NA))
expect_identical(rows$lof_flags, c(NA, NA, NA, "NAGNAG_SITE", "NON_CAN_SPLICE", NA))
expect_identical(rows$lof_info[1], "PERCENTILE:0.742424242424242,BP_DIST:84,DIST_FROM_LAST_EXON:50,50_BP_RULE:FAIL")

small <- rduckvep_lof_sql(con, "lof_ann", "lof_tx", "lof_ref", min_intron_size = 0L)
rows <- dbGetQuery(con, paste0("SELECT * FROM query(", as.character(dbQuoteString(con, small)),
                               ") WHERE event_index = 3"))
expect_identical(rows$lof_filter, "5UTR_SPLICE")
expect_error(rduckvep_lof_sql(con, "lof_ann", "lof_tx", "lof_ref", min_intron_size = -1L),
             "min_intron_size")
dbDisconnect(con, shutdown = TRUE)
