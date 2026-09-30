-- v2 model sink scenarios, run in ONE session (the staging area is per database) by
-- test/scripts/run_v2_tests.py without -bail. A line "-- expect error: <text>" before a
-- statement says it must fail with that text; every other statement must succeed, and a
-- SELECT that returns true is an assertion (error() otherwise).

CREATE TEMP TABLE m_regions AS SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER)) t(seq_region);
CREATE TEMP TABLE m_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, 100::UBIGINT transcript_start, 250::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 120::UBIGINT cds_start, 240::UBIGINT cds_end, 'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence, 1::UTINYINT codon_table, 'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence, 'ACGTACGTAC'::BLOB post_cds_sequence;
CREATE TEMP TABLE m_exons AS SELECT * FROM (VALUES (0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT), (0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)) t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase);

-- TEMP tables are visible to the COPY (they run in the caller's session and transaction).
COPY (SELECT * FROM m_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'temp-model', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_transcripts ORDER BY seq_region, transcript_start) TO 'x' (FORMAT duckvep_stage, MODEL 'temp-model', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'temp-model', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN duckvep_model_publish('temp-model') THEN true ELSE error('publish of a TEMP-table model') END;
SELECT CASE WHEN _duckvep_model_fingerprint('temp-model') = 16162758230738989510 THEN true ELSE error('TEMP model equals the README model fingerprint') END;

-- Uncommitted rows are visible: an extra region inserted in an open transaction changes the model.
CREATE TABLE u_regions AS SELECT * FROM m_regions;
BEGIN;
INSERT INTO u_regions VALUES (2::UINTEGER);
COPY (SELECT * FROM u_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'uncommitted', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
ROLLBACK;
COPY (SELECT * FROM m_transcripts ORDER BY seq_region, transcript_start) TO 'x' (FORMAT duckvep_stage, MODEL 'uncommitted', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'uncommitted', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN duckvep_model_publish('uncommitted') THEN true ELSE error('publish with staged uncommitted rows') END;
SELECT CASE WHEN _duckvep_model_fingerprint('uncommitted') <> _duckvep_model_fingerprint('temp-model') THEN true ELSE error('the uncommitted region was staged') END;

-- A COPY whose query fails stages nothing.
-- expect error: boom
COPY (SELECT CASE WHEN seq_region = 1 THEN error('boom') ELSE seq_region END AS seq_region FROM m_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'failed-copy', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_transcripts ORDER BY seq_region, transcript_start) TO 'x' (FORMAT duckvep_stage, MODEL 'failed-copy', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'failed-copy', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: relation 'regions' is not staged
SELECT duckvep_model_publish('failed-copy');
SELECT CASE WHEN _duckvep_model_fingerprint('failed-copy') IS NULL THEN true ELSE error('a failed COPY published a model') END;

-- Invalid staged rows fail the publish, publish nothing, and consume the staging.
COPY (SELECT * FROM m_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'bad-rows', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (5::UINTEGER AS transcript_index) FROM m_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'bad-rows', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'bad-rows', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: transcript_index is duplicated or outside the dense zero-based range
SELECT duckvep_model_publish('bad-rows');
SELECT CASE WHEN _duckvep_model_fingerprint('bad-rows') IS NULL AND NOT duckvep_model_drop('bad-rows') THEN true ELSE error('a failed publish installed a model') END;
-- expect error: relation 'regions' is not staged
SELECT duckvep_model_publish('bad-rows');

-- A wrong column type is the same error as on v1.
COPY (SELECT seq_region::BIGINT AS seq_region FROM m_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'bad-type', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'bad-type', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM m_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'bad-type', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: query column seq_region has the wrong type
SELECT duckvep_model_publish('bad-type');

-- Dropping a model name also discards its staged rows.
COPY (SELECT * FROM m_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'discarded', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN NOT duckvep_model_drop('discarded') THEN true ELSE error('drop of a staged-only name reports no model') END;
-- expect error: relation 'regions' is not staged
SELECT duckvep_model_publish('discarded');

-- Names and options are checked.
-- expect error: model name already exists
COPY (SELECT * FROM m_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'temp-model', RELATION 'regions', USE_TMP_FILE FALSE);
-- expect error: RELATION must be
COPY (SELECT * FROM m_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'n', RELATION 'nope', USE_TMP_FILE FALSE);
-- expect error: MODEL and RELATION options are required
COPY (SELECT * FROM m_regions) TO 'x' (FORMAT duckvep_stage, USE_TMP_FILE FALSE);
-- expect error: name must be a non-empty string
SELECT duckvep_model_publish('');
-- expect error: reference_fasta must be a non-empty string
SELECT duckvep_model_publish('temp-model', {'reference_fasta': ''});

-- Restaging a relation replaces its rows; the drop releases the model.
SELECT CASE WHEN duckvep_model_drop('temp-model') AND duckvep_model_drop('uncommitted') THEN true ELSE error('drop of published models') END;
SELECT CASE WHEN _duckvep_model_fingerprint('temp-model') IS NULL THEN true ELSE error('dropped model is gone') END;

-- The statements duckvep_model_load_sql returns load a model.
SELECT CASE WHEN len(duckvep_model_load_sql('m1', 'SELECT 1', 'SELECT 2', 'SELECT 3')) = 4 AND len(duckvep_model_load_sql('m1', 'SELECT 1', 'SELECT 2', 'SELECT 3', {'mature_mirna_query': 'SELECT 4', 'reference_fasta': '/x', 'transcript_coverage_complete': true})) = 5 THEN true ELSE error('statement counts') END;
SELECT CASE WHEN duckvep_model_load_sql('it''s', 'SELECT 1 -- c', 'SELECT 2', 'SELECT 3')[4] = 'SELECT duckvep_model_publish(''it''''s'')' THEN true ELSE error('publish statement text') END;
SELECT CASE WHEN contains(duckvep_model_load_sql('m', 'SELECT 1 -- c', 'SELECT 2', 'SELECT 3')[1], E'-- c\n) TO ''duckvep_stage''') THEN true ELSE error('trailing comment stays inside the parentheses') END;
-- expect error: arguments must be non-empty strings
SELECT duckvep_model_load_sql('m', '', 'SELECT 2', 'SELECT 3');
-- expect error: unknown option
SELECT duckvep_model_load_sql('m', 'SELECT 1', 'SELECT 2', 'SELECT 3', {'nope': 1});
