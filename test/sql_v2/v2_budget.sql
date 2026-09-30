-- v2 budget scenarios (issue #8 slice 6), run in ONE session like v2_model.sql. They mirror test/sql/duckvep_budget.test
-- on the v2 paths: a model publish, an annotation, a haplotype scan and a worker lease over a tiny budget give the explicit
-- capacity error, publish nothing and leave the connection and the loaded model usable.

CREATE TEMP TABLE b_regions AS SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER)) t(seq_region);
CREATE TEMP TABLE b_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, 100::UBIGINT transcript_start, 250::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 120::UBIGINT cds_start, 240::UBIGINT cds_end, 'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence, 1::UTINYINT codon_table, 'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence, 'ACGTACGTAC'::BLOB post_cds_sequence;
CREATE TEMP TABLE b_exons AS SELECT * FROM (VALUES (0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT), (0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)) t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase);
CREATE TEMP TABLE b_events AS SELECT 1::UBIGINT event_index, 1::UINTEGER seq_region, 130::UBIGINT AS position, 'A' reference, 'C' alternate, NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change, NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position;
SELECT duckvep_worker_limits_set(6, 134217728, 67108864, 0);

COPY (SELECT * FROM b_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-ok', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM b_transcripts ORDER BY seq_region, transcript_start) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-ok', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM b_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-ok', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN duckvep_model_publish('budget-ok') THEN true ELSE error('publish under the default budget') END;
SELECT CASE WHEN (SELECT max(current_bytes) FILTER (WHERE owner = 'model') > 0 AND max(current_bytes) FILTER (WHERE owner = 'index') >= 0 FROM duckvep_native_budget()) THEN true ELSE error('the published model is charged to the model owner') END;
CREATE TEMP TABLE b_golden AS SELECT * FROM query(duckvep_annotate_sql('b_events', 'budget-ok', {'upstream_distance': 0, 'downstream_distance': 0}));
SELECT CASE WHEN count(*) > 0 THEN true ELSE error('golden annotation is empty') END FROM b_golden;

-- A publish over budget: explicit capacity error, nothing published, nothing left charged.
COPY (SELECT * FROM b_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-over', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM b_transcripts ORDER BY seq_region, transcript_start) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-over', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM b_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-over', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SET VARIABLE b_before = (SELECT current_bytes FROM duckvep_native_budget() WHERE owner = 'total');
SELECT CASE WHEN duckvep_native_budget_set(getvariable('b_before') + 65536) - getvariable('b_before') = 65536 THEN true ELSE error('budget set') END;
-- expect error: capacity error
SELECT duckvep_model_publish('budget-over');
SELECT CASE WHEN duckvep_native_budget_set(4294967296) = 4294967296 THEN true ELSE error('budget restore') END;
SELECT CASE WHEN (SELECT current_bytes = getvariable('b_before') FROM duckvep_native_budget() WHERE owner = 'total') THEN true ELSE error('the refused publish left a charge behind') END;
SELECT CASE WHEN _duckvep_model_fingerprint('budget-over') IS NULL THEN true ELSE error('a refused publish published a model') END;
-- expect error: unknown model name
SELECT count(*) FROM query(duckvep_annotate_sql('b_events', 'budget-over'));
-- The refused publish consumed its staging: publishing again needs a new staging.
-- expect error: is not staged for this model
SELECT duckvep_model_publish('budget-over');

-- An annotation over budget (the pooled worker is warm, a budget at the current charge leaves it no room).
SELECT count(*) FROM query(duckvep_annotate_sql('b_events', 'budget-ok', {'upstream_distance': 0, 'downstream_distance': 0}));
SET VARIABLE b_before = (SELECT current_bytes FROM duckvep_native_budget() WHERE owner = 'total');
SELECT CASE WHEN duckvep_native_budget_set(getvariable('b_before')) = getvariable('b_before') THEN true ELSE error('budget set') END;
-- expect error: capacity error
SELECT count(*) FROM query(duckvep_annotate_sql('b_events', 'budget-ok'));
SELECT CASE WHEN duckvep_native_budget_set(4294967296) = 4294967296 THEN true ELSE error('budget restore') END;
SELECT CASE WHEN (SELECT count(*) FROM ((SELECT * FROM b_golden EXCEPT SELECT * FROM query(duckvep_annotate_sql('b_events', 'budget-ok', {'upstream_distance': 0, 'downstream_distance': 0}))) UNION ALL (SELECT * FROM query(duckvep_annotate_sql('b_events', 'budget-ok', {'upstream_distance': 0, 'downstream_distance': 0})) EXCEPT SELECT * FROM b_golden))) = 0 THEN true ELSE error('annotation after a refusal differs') END;

-- Per-worker leases: an over-limit vector is an explicit capacity error.
SELECT duckvep_worker_limits_set(6, 1024, 67108864, 0);
-- expect error: capacity error: per-worker scratch lease
SELECT count(*) FROM query(duckvep_annotate_sql('b_events', 'budget-ok'));
SELECT duckvep_worker_limits_set(6, 134217728, 1024, 0);
-- expect error: capacity error: per-worker emitted-output lease
SELECT count(*) FROM query(duckvep_annotate_sql('b_events', 'budget-ok'));
SELECT duckvep_worker_limits_set(6, 134217728, 67108864, 0);

-- A haplotype scan over budget: the workspace is refused, the job is released, the connection is usable.
COPY (SELECT 1::UBIGINT, 1::UINTEGER, 130::UBIGINT, 'A', 'C', 1::UINTEGER, 0::UINTEGER, 0::UINTEGER, [0, 1]::INTEGER[], [false, true]::BOOLEAN[], NULL::BIGINT, [NULL]::BIGINT[], 1::BIGINT, 1::BIGINT, 1::BIGINT) TO 'x' (FORMAT duckvep_stage, JOB 'budget-job', MODEL 'budget-ok', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SET VARIABLE b_before = (SELECT current_bytes FROM duckvep_native_budget() WHERE owner = 'total');
SELECT CASE WHEN duckvep_native_budget_set(getvariable('b_before') + 4096) > 0 THEN true ELSE error('budget set') END;
-- expect error: duckvep_haplotypes:
SELECT count(*) FROM duckvep_haplotype_scan('budget-job');
SELECT CASE WHEN duckvep_native_budget_set(4294967296) = 4294967296 THEN true ELSE error('budget restore') END;
SELECT CASE WHEN (SELECT current_bytes = getvariable('b_before') FROM duckvep_native_budget() WHERE owner = 'total') THEN true ELSE error('the refused scan left a charge behind') END;
-- expect error: is not staged
SELECT count(*) FROM duckvep_haplotype_scan('budget-job');

-- The coding_calls reader reserves htslib's buffers from the same budget.
COPY (SELECT seq_region, 1000::UBIGINT sequence_length, ('c' || seq_region) seq_region_name FROM b_regions ORDER BY seq_region) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-cc', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM b_transcripts ORDER BY seq_region, transcript_start) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-cc', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM b_exons ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'budget-cc', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN duckvep_model_publish('budget-cc') THEN true ELSE error('publish') END;
SELECT CASE WHEN (SELECT count(*) >= 0 FROM duckvep_coding_calls('budget-cc', 'test/data/coding_calls/cases.vcf')) THEN true ELSE error('coding_calls under the default budget') END;
SET VARIABLE b_before = (SELECT current_bytes FROM duckvep_native_budget() WHERE owner = 'total');
SELECT CASE WHEN duckvep_native_budget_set(getvariable('b_before') + 4096) > 0 THEN true ELSE error('budget set') END;
-- expect error: capacity error
SELECT count(*) FROM duckvep_coding_calls('budget-cc', 'test/data/coding_calls/cases.vcf');
SELECT CASE WHEN duckvep_native_budget_set(4294967296) = 4294967296 THEN true ELSE error('budget restore') END;
SELECT CASE WHEN (SELECT current_bytes = getvariable('b_before') FROM duckvep_native_budget() WHERE owner = 'total') THEN true ELSE error('the refused reader left a charge behind') END;
SELECT CASE WHEN duckvep_model_drop('budget-cc') THEN true ELSE error('drop') END;

-- The function surface: bad arguments are errors, not clipped values.
-- expect error: the budget must be a positive byte count
SELECT duckvep_native_budget_set(0);
-- expect error: capacity error: 1 bytes is below
SELECT duckvep_native_budget_set(1);
-- expect error: workers must be 1..1024
SELECT duckvep_worker_limits_set(0, 1, 1, 1);

-- Dropping the model returns every model and index byte.
SELECT CASE WHEN duckvep_model_drop('budget-ok') THEN true ELSE error('drop') END;
SELECT CASE WHEN (SELECT sum(current_bytes) FROM duckvep_native_budget() WHERE owner IN ('model', 'index')) = 0 THEN true ELSE error('model and index bytes returned') END;
