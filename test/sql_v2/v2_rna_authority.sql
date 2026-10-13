CREATE TEMP TABLE a_regions AS SELECT 0::UINTEGER seq_region, 12::UBIGINT sequence_length;
CREATE TEMP TABLE a_transcripts AS SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region, 1::UBIGINT transcript_start, 12::UBIGINT transcript_end, -1::TINYINT strand, 0::UINTEGER gene_index, 0::UBIGINT transcript_flags, 4::UBIGINT cds_start, 12::UBIGINT cds_end, 'ATGCCCTAA'::BLOB cds_sequence, 2::UTINYINT codon_table, ''::BLOB pre_cds_sequence, 'AGG'::BLOB post_cds_sequence, 'ATGCCCTAAAGG'::BLOB cdna_sequence, 1::UBIGINT rna_cds_cdna_start, 9::UBIGINT rna_cds_cdna_end, 2::UTINYINT rna_codon_table;
CREATE TEMP TABLE a_exons AS SELECT 0::UINTEGER transcript_index, 1::UBIGINT exon_start, 12::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase;
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-v2', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-v2', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-v2', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN duckvep_model_publish('authority-v2') THEN true ELSE error('17-column v2 transcript model') END;
SELECT CASE WHEN duckvep_model_save('authority-v2', 'build/v2_authority.dvsnap') THEN true ELSE error('v2 authority save') END;
SELECT CASE WHEN duckvep_model_restore('authority-v2-copy', 'build/v2_authority.dvsnap') THEN true ELSE error('v2 authority restore') END;
SELECT CASE WHEN _duckvep_model_fingerprint('authority-v2') = _duckvep_model_fingerprint('authority-v2-copy') THEN true ELSE error('v2 authority fingerprint') END;
SELECT CASE WHEN duckvep_model_drop('authority-v2') AND duckvep_model_drop('authority-v2-copy') THEN true ELSE error('v2 authority cleanup') END;
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-start-overflow', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (4294967296::UBIGINT AS rna_cds_cdna_start) FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-start-overflow', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-start-overflow', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: original RNA coding authority
SELECT duckvep_model_publish('authority-start-overflow');
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-end-overflow', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (4294967296::UBIGINT AS rna_cds_cdna_end) FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-end-overflow', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-end-overflow', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: original RNA coding authority
SELECT duckvep_model_publish('authority-end-overflow');
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-wrap-valid', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (4294967297::UBIGINT AS rna_cds_cdna_start) FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-wrap-valid', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-wrap-valid', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: original RNA coding authority
SELECT duckvep_model_publish('authority-wrap-valid');
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-end-wrap-valid', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (4294967305::UBIGINT AS rna_cds_cdna_end) FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-end-wrap-valid', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-end-wrap-valid', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: original RNA coding authority
SELECT duckvep_model_publish('authority-end-wrap-valid');
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-partial', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (NULL::UBIGINT AS rna_cds_cdna_end) FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-partial', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-partial', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: original RNA coding authority
SELECT duckvep_model_publish('authority-partial');
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-code', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * REPLACE (99::UTINYINT AS rna_codon_table) FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-code', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-code', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
-- expect error: original RNA coding authority
SELECT duckvep_model_publish('authority-code');
COPY (SELECT * FROM a_regions) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-reuse', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_transcripts) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-reuse', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
COPY (SELECT * FROM a_exons) TO 'x' (FORMAT duckvep_stage, MODEL 'authority-reuse', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE);
SELECT CASE WHEN duckvep_model_publish('authority-reuse') THEN true ELSE error('v2 authority reuse') END;
SELECT CASE WHEN duckvep_model_drop('authority-reuse') THEN true ELSE error('v2 authority reuse cleanup') END;
