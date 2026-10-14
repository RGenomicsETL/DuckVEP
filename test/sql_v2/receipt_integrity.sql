CREATE TABLE receipt_regions AS SELECT * FROM query(duckvep_ensembl_regions_sql('duckvep_core', 'duckvep_reference_chunks', 'GRCh38'));
CREATE TABLE receipt_model AS SELECT * FROM query(duckvep_ensembl_transcripts_sql('duckvep_core', 'duckvep_reference_chunks', 'GRCh38'));
CREATE TABLE receipt_model_before AS SELECT * FROM receipt_model;
CREATE TABLE receipt_regions_before AS SELECT * FROM receipt_regions;
CREATE VIEW receipt_legacy AS SELECT * EXCLUDE (rna_cds_cdna_start, rna_cds_cdna_end, rna_codon_table) FROM receipt_model;
CREATE VIEW receipt_null AS SELECT * REPLACE (NULL::UBIGINT AS rna_cds_cdna_start, NULL::UBIGINT AS rna_cds_cdna_end, NULL::UTINYINT AS rna_codon_table) FROM receipt_model;
CREATE VIEW receipt_extra AS SELECT *, 'ignored' AS model, 'ignored' AS authority, 'ignored' AS legacy_row, 'ignored' AS row_payloads FROM receipt_model;
CREATE VIEW receipt_extra_legacy AS SELECT * EXCLUDE (rna_cds_cdna_start, rna_cds_cdna_end, rna_codon_table) FROM receipt_extra;
CREATE VIEW receipt_reordered AS SELECT rna_codon_table, rna_cds_cdna_end, rna_cds_cdna_start, * EXCLUDE (rna_cds_cdna_start, rna_cds_cdna_end, rna_codon_table) FROM receipt_model;
CREATE VIEW receipt_start AS SELECT * REPLACE (CASE WHEN source_transcript_id = 32 THEN rna_cds_cdna_start + 1 ELSE rna_cds_cdna_start END AS rna_cds_cdna_start) FROM receipt_model;
CREATE VIEW receipt_end AS SELECT * REPLACE (CASE WHEN source_transcript_id = 32 THEN rna_cds_cdna_end + 1 ELSE rna_cds_cdna_end END AS rna_cds_cdna_end) FROM receipt_model;
CREATE VIEW receipt_code AS SELECT * REPLACE (CASE WHEN source_transcript_id = 32 THEN 1::UTINYINT ELSE rna_codon_table END AS rna_codon_table) FROM receipt_model;
CREATE VIEW receipt_flags AS SELECT * REPLACE (CASE WHEN source_transcript_id = 32 THEN transcript_flags + 1 ELSE transcript_flags END AS transcript_flags) FROM receipt_model;
CREATE VIEW receipt_wide AS SELECT * REPLACE (CASE WHEN source_transcript_id = 32 THEN rna_cds_cdna_start + 4294967296::UBIGINT ELSE rna_cds_cdna_start END AS rna_cds_cdna_start) FROM receipt_model;
CREATE VIEW receipt_transcript_source AS SELECT * REPLACE (CASE WHEN source_transcript_id = 32 THEN source_transcript_id + 1 ELSE source_transcript_id END AS source_transcript_id) FROM receipt_model;
CREATE VIEW receipt_regions_source AS SELECT * REPLACE (source_seq_region_id + 1 AS source_seq_region_id) FROM receipt_regions;
CREATE VIEW receipt_regions_topology AS SELECT * REPLACE (true AS circular) FROM receipt_regions;
CREATE MACRO receipt_for(t, r) AS TABLE SELECT * FROM query(duckvep_model_receipt_sql(r, t, 'Ensembl', '116', 'GRCh38', repeat('a', 64), repeat('b', 64), 'fixture transcripts'));
CREATE MACRO receipt_hash(t) AS (SELECT model_sha256 FROM receipt_for(t, 'receipt_regions'));
CREATE TABLE receipt_before AS SELECT * FROM receipt_for('receipt_model', 'receipt_regions');
CREATE TABLE receipt_after AS SELECT * FROM receipt_for('receipt_model', 'receipt_regions');
SELECT CASE WHEN (SELECT count(*) FROM receipt_model WHERE source_transcript_id = 32 AND sequence_withheld_reason = 'rna_edit' AND cds_sequence IS NULL AND rna_cds_cdna_start IS NOT NULL AND rna_cds_cdna_end IS NOT NULL AND rna_codon_table IS NOT NULL) = 1
  AND (SELECT count(*) FROM receipt_model) = 3
  AND receipt_hash('receipt_legacy') = receipt_hash('receipt_null')
  AND receipt_hash('receipt_legacy') <> receipt_hash('receipt_model')
  AND receipt_hash('receipt_legacy') = receipt_hash('receipt_extra_legacy')
  AND receipt_hash('receipt_model') = receipt_hash('receipt_extra')
  AND receipt_hash('receipt_model') = receipt_hash('receipt_reordered')
  AND receipt_hash('receipt_model') <> receipt_hash('receipt_start')
  AND receipt_hash('receipt_model') <> receipt_hash('receipt_end')
  AND receipt_hash('receipt_model') <> receipt_hash('receipt_code')
  AND receipt_hash('receipt_model') <> receipt_hash('receipt_flags')
  AND receipt_hash('receipt_model') <> receipt_hash('receipt_wide')
  AND receipt_hash('receipt_model') <> receipt_hash('receipt_transcript_source')
  AND (SELECT model_sha256 FROM receipt_before) <> (SELECT model_sha256 FROM receipt_for('receipt_model', 'receipt_regions_source'))
  AND (SELECT model_sha256 FROM receipt_before) <> (SELECT model_sha256 FROM query(duckvep_model_receipt_sql('receipt_regions', 'receipt_model', 'Ensembl', '117', 'GRCh38', repeat('a', 64), repeat('b', 64), 'fixture transcripts')))
  AND (SELECT topology_sha256 FROM receipt_before) = (SELECT topology_sha256 FROM receipt_for('receipt_code', 'receipt_regions'))
  AND (SELECT model_sha256 FROM receipt_before) = (SELECT model_sha256 FROM receipt_for('receipt_model', 'receipt_regions_topology'))
  AND (SELECT topology_sha256 FROM receipt_before) <> (SELECT topology_sha256 FROM receipt_for('receipt_model', 'receipt_regions_topology'))
  AND NOT EXISTS ((SELECT * FROM receipt_before EXCEPT ALL SELECT * FROM receipt_after) UNION ALL (SELECT * FROM receipt_after EXCEPT ALL SELECT * FROM receipt_before))
  AND NOT EXISTS ((SELECT * FROM receipt_model EXCEPT ALL SELECT * FROM receipt_model_before) UNION ALL (SELECT * FROM receipt_model_before EXCEPT ALL SELECT * FROM receipt_model))
  AND NOT EXISTS ((SELECT * FROM receipt_regions EXCEPT ALL SELECT * FROM receipt_regions_before) UNION ALL (SELECT * FROM receipt_regions_before EXCEPT ALL SELECT * FROM receipt_regions))
  AND (SELECT count(*) FROM receipt_model WHERE source_transcript_id = 32 AND rna_codon_table = 2) = 1
  AND (SELECT count(*) FROM receipt_regions WHERE circular) = 0
  THEN 'PASS' ELSE error('receipt integrity regression') END AS receipt_integrity;
SELECT 'legacy', receipt_hash('receipt_legacy') UNION ALL SELECT 'null', receipt_hash('receipt_null') UNION ALL SELECT 'authority', receipt_hash('receipt_model') UNION ALL SELECT 'start', receipt_hash('receipt_start') UNION ALL SELECT 'end', receipt_hash('receipt_end') UNION ALL SELECT 'code', receipt_hash('receipt_code') UNION ALL SELECT 'flags', receipt_hash('receipt_flags') UNION ALL SELECT 'wide', receipt_hash('receipt_wide') ORDER BY 1;
SELECT 'complete_receipt', to_json(r) FROM receipt_before r;
