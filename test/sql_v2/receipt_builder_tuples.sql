WITH checks(label, sql_text, expected_md5, expected_length) AS (VALUES
 ('defaults', duckvep_ensembl_transcripts_sql('core', 'reference_chunks', 'GRCh38'), '7acafe02a20f054e56926e767a306d93', 31672),
 ('qualified', duckvep_ensembl_transcripts_sql('core_1', 'sch.chunks', 'GRCh37'), '3869ec30614101ba2094618f1ac3442c', 31704),
 ('quoted', duckvep_ensembl_transcripts_sql('co"re', 'chunks"', 'asm''bly'), 'c00a4cd8b68369f8f4dc99a67b3df313', 31698),
 ('species1', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 1}), 'e591f1adff1b0786bb08922543c734f6', 31643),
 ('species0', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 0}), 'b0f328fcdaceabb387d89815cc39a67e', 31643),
 ('species7', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 7}), 'b44ede8034900cf11f742752822b4b6d', 31643),
 ('species_minus2', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': -2::BIGINT}), '99f1c0d5fca105052766ed9cf5c0daf5', 31644),
 ('species255', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 255::UTINYINT}), '183377b1ade28277a2f85cbc5101bd4d', 31645),
 ('species_max', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 9223372036854775807::UBIGINT}), 'f5408908b0f709eb25197f77518e1a6a', 31661),
 ('species_null_integer', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': NULL::INTEGER}), '0dd36f0ca8a1f30276931da02571713f', 31646),
 ('species_null', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': NULL}), '0dd36f0ca8a1f30276931da02571713f', 31646),
 ('null_options', duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', NULL), 'e591f1adff1b0786bb08922543c734f6', 31643)
) SELECT label, md5(sql_text), length(sql_text), CASE WHEN md5(sql_text) = expected_md5 AND length(sql_text) = expected_length THEN 'PASS' ELSE error('transcript builder tuple mismatch: ' || label) END FROM checks ORDER BY label;
SELECT count(*), md5(string_agg(sql_text, '' ORDER BY i)), CASE WHEN count(*) = 1200 AND md5(string_agg(sql_text, '' ORDER BY i)) = 'dad4c327842bd8a669377da5eeb7baf0' THEN 'PASS' ELSE error('transcript builder batch tuple mismatch') END FROM (SELECT i, duckvep_ensembl_transcripts_sql('c' || (i % 4), 'r' || (i % 3), 'a', {'species_id': i % 5}) AS sql_text FROM range(1200) t(i));
