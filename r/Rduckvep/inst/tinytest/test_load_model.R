library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)

dbExecute(con, "CREATE TABLE lm_regions AS SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER)) t(seq_region)")
dbExecute(con, "CREATE TABLE lm_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region,
  100::UBIGINT transcript_start, 250::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index,
  3::UBIGINT transcript_flags, 120::UBIGINT cds_start, 240::UBIGINT cds_end,
  'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence,
  1::UTINYINT codon_table, 'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence, 'ACGTACGTAC'::BLOB post_cds_sequence")
dbExecute(con, "CREATE TABLE lm_exons AS SELECT * FROM (VALUES
  (0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT),
  (0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT))
  t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase)")

expect_true(rduckvep_load_model(con, "lm", "SELECT * FROM lm_regions ORDER BY seq_region",
  "SELECT * FROM lm_transcripts ORDER BY seq_region, transcript_start",
  "SELECT * FROM lm_exons ORDER BY transcript_index, exon_start"))
expect_error(rduckvep_load_model(con, "lm", "SELECT * FROM lm_regions ORDER BY seq_region",
  "SELECT * FROM lm_transcripts ORDER BY seq_region, transcript_start",
  "SELECT * FROM lm_exons ORDER BY transcript_index, exon_start"), "already exists")
expect_error(rduckvep_load_model(con, "", "a", "b", "c"), "non-empty")
expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('lm') AS dropped")$dropped)
dbDisconnect(con, shutdown = TRUE)
