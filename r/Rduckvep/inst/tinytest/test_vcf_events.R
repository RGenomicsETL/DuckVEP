library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
dbExecute(con, "CREATE TEMP TABLE vcf_regions AS SELECT * FROM (VALUES
  (0::UINTEGER, 'chr1')) t(seq_region, seq_region_name)")
# read_bcf columns; VEP 116 starts a symbolic record at POS + 1 and skips a
# record whose ALT list mixes a literal allele with <NON_REF>.
dbExecute(con, "CREATE TEMP TABLE vcf_records AS SELECT * FROM (VALUES
  ('chr1', 100::BIGINT, 'snv', 'A', ['G'], NULL::INTEGER, NULL::VARCHAR, NULL::INTEGER[]),
  ('chr1', 200, 'del', 'C', ['<DEL>'], 260, 'DEL', [-60]),
  ('chr1', 300, 'gvcf', 'T', ['C', '<NON_REF>'], NULL, NULL, NULL),
  ('chr2', 5, 'elsewhere', 'A', ['G'], NULL, NULL, NULL)
) t(\"CHROM\", \"POS\", \"ID\", \"REF\", \"ALT\", \"INFO_END\", \"INFO_SVTYPE\", \"INFO_SVLEN\")")
sql <- rduckvep_vcf_events_sql(con, "vcf_records", "vcf_regions")
expect_true(is.character(sql) && length(sql) == 1L)
rows <- dbGetQuery(con, paste0("SELECT * FROM (", sql, ") ORDER BY event_index"))
# Unsupported rows follow the ok rows, by source contig and position.
expect_identical(rows$source_id, c("snv", "del", "gvcf", "gvcf", "elsewhere"))
expect_identical(rows$event_status, c("ok", "ok", rep("unsupported", 3)))
expect_equal(rows$position[1:2], c(100, 201))
expect_equal(rows$end_position[2], 260)
expect_identical(rows$structural_type[2], "DEL")
expect_identical(rows$event_reason[3:5],
                 c(rep("multiallelic_structural_record", 2), "contig_not_in_model"))
expect_error(rduckvep_vcf_events_sql(con, "vcf_records", NA_character_))

rduckvep_load_model(
  con, "r-point",
  "SELECT 0::UINTEGER AS seq_region, 100::UBIGINT AS sequence_length",
  paste("SELECT 0::UINTEGER AS transcript_index, 0::UINTEGER AS seq_region,",
    "11::UBIGINT AS transcript_start, 34::UBIGINT AS transcript_end,",
    "1::TINYINT AS strand, 0::UINTEGER AS gene_index, 0::UBIGINT AS transcript_flags,",
    "11::UBIGINT AS cds_start, 34::UBIGINT AS cds_end,",
    "'ATGGCCGCTGAATTTGGTCCCTAA'::BLOB AS cds_sequence,",
    "1::UTINYINT AS codon_table, ''::BLOB AS pre_cds_sequence,",
    "''::BLOB AS post_cds_sequence"),
  paste("SELECT 0::UINTEGER AS transcript_index, 11::UBIGINT AS exon_start,",
    "34::UBIGINT AS exon_end, 1::UBIGINT AS exon_cdna_start,",
    "24::UBIGINT AS exon_cdna_end, 0::TINYINT AS phase, -1::TINYINT AS end_phase")
)
dbExecute(con, "CREATE TEMP TABLE point_records AS SELECT * FROM (VALUES
  ('chr1', 33::BIGINT, 'mobile', 'A', ['<INS:ME>'], NULL::INTEGER),
  ('chr1', 34::BIGINT, 'interbase', 'A', ['<INS>'], 34::INTEGER)
) t(CHROM, POS, ID, REF, ALT, INFO_END)")
point_sql <- rduckvep_vcf_events_sql(con, "point_records", "vcf_regions")
dbExecute(con, paste("CREATE TEMP TABLE point_events AS", point_sql))
points <- dbGetQuery(con, "SELECT source_id, source_alt, position, end_position,
  structural_geometry FROM point_events ORDER BY source_id")
expect_identical(points$source_alt, c("<INS>", "<INS:ME>"))
expect_equal(points$position, c(34, 34))
expect_equal(points$end_position, c(34, 34))
expect_identical(points$structural_geometry, c(NA_character_, "POINT"))

annotation_sql <- rduckvep_annotate_sql(con, "point_events", "r-point",
  rich = TRUE, upstream_distance = 0, downstream_distance = 0)
annotations <- dbGetQuery(con, paste("SELECT e.source_id, a.consequence, a.duckvep_reason",
  "FROM (", annotation_sql, ") a JOIN point_events e USING (event_index)",
  "ORDER BY e.source_id"))
expect_identical(annotations$source_id, c("interbase", "mobile"))
expect_identical(annotations$consequence,
  c("sequence_variant", "feature_elongation&coding_sequence_variant"))
expect_identical(annotations$duckvep_reason, c("no_feature_in_loaded_model", NA_character_))
dbDisconnect(con, shutdown = TRUE)
