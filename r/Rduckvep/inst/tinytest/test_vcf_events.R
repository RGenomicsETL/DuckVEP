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
dbDisconnect(con, shutdown = TRUE)
