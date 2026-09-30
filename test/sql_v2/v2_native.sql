-- v2-only assertions, run in one session by test/scripts/run_v2_tests.py after
-- LOAD. Each statement returns true or raises error(...).

-- LOAD registered exactly the three ported functions and no catalog object.
SELECT CASE WHEN (SELECT count(DISTINCT function_name) FROM duckdb_functions() WHERE function_name IN ('duckvep_so_terms', 'duckvep_allele_geometry', 'duckvep_breakend_geometry')) = 3 THEN true ELSE error('ported functions missing') END;
SELECT CASE WHEN (SELECT count(*) FROM duckdb_functions() WHERE function_type = 'macro' AND function_name LIKE '%duckvep%') = 0 THEN true ELSE error('LOAD persisted a macro') END;
SELECT CASE WHEN (SELECT count(*) FROM duckdb_views() WHERE NOT internal) = 0 AND (SELECT count(*) FROM duckdb_tables()) = 0 THEN true ELSE error('LOAD created a table or view') END;

-- Table function: 41 rows, resumable across several output chunks (16 rows each).
SELECT CASE WHEN count(*) = 41 AND count(DISTINCT bit_index) = 41 AND min(bit_index) = 0 AND max(bit_index) = 40 THEN true ELSE error('so_terms rows') END FROM duckvep_so_terms();
SELECT CASE WHEN bool_and(consequence_mask = (1::UBIGINT << bit_index)) THEN true ELSE error('so_terms mask') END FROM duckvep_so_terms();
SELECT CASE WHEN bool_and(consequence <> '' AND impact IN ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')) THEN true ELSE error('so_terms strings') END FROM duckvep_so_terms();
SELECT CASE WHEN (SELECT list(bit_index) FROM duckvep_so_terms()) = list_transform(range(41), lambda x: x::UTINYINT) THEN true ELSE error('so_terms order') END;
SELECT CASE WHEN count(*) = 1681 THEN true ELSE error('so_terms self cross join') END FROM duckvep_so_terms() a, duckvep_so_terms() b;
SELECT CASE WHEN (SELECT count(*) FROM (SELECT * FROM duckvep_so_terms() LIMIT 5)) = 5 THEN true ELSE error('so_terms limit') END;
SELECT CASE WHEN (SELECT count(*) FROM duckvep_so_terms()) = 41 THEN true ELSE error('so_terms rescan after limit') END;

-- Struct NULL invariants: a NULL struct row has NULL children; valid rows keep valid children.
SELECT CASE WHEN g IS NULL AND g.kind_code IS NULL AND g.insertion_boundary0 IS NULL THEN true ELSE error('allele NULL row') END FROM (SELECT duckvep_allele_geometry(NULL, 'A', 'G') AS g);
SELECT CASE WHEN g IS NOT NULL AND g.insertion_boundary0 IS NULL AND g.interbase = false THEN true ELSE error('allele SNV boundary') END FROM (SELECT duckvep_allele_geometry(5, 'A', 'G') AS g);
SELECT CASE WHEN g.interbase AND g.insertion_boundary0 IS NOT NULL THEN true ELSE error('allele insertion boundary') END FROM (SELECT duckvep_allele_geometry(5, 'A', 'AC') AS g);
SELECT CASE WHEN g IS NULL AND g.mate_chrom IS NULL AND g.replacement_sequence IS NULL THEN true ELSE error('breakend non-BND row') END FROM (SELECT duckvep_breakend_geometry('ACGT') AS g);
SELECT CASE WHEN g IS NOT NULL AND g.mate_chrom IS NULL AND g.mate_position IS NULL AND g.mate_extends_right IS NULL AND g.replacement_sequence = 'A' AND g.local_join_after IS NOT NULL THEN true ELSE error('single breakend fields') END FROM (SELECT duckvep_breakend_geometry('.A') AS g);

-- NULL-heavy vectors copy cleanly (TRY, projection, aggregation read every child).
SELECT CASE WHEN count(*) = 1000 AND count(g) = 500 THEN true ELSE error('breakend null-heavy copy') END FROM (SELECT duckvep_breakend_geometry(CASE WHEN i % 2 = 0 THEN 'A]1:' || i || ']' ELSE NULL END) AS g FROM range(1000) t(i));
SELECT CASE WHEN count(*) = 1000 AND count(g) = 600 THEN true ELSE error('allele null-heavy copy') END FROM (SELECT duckvep_allele_geometry(CASE WHEN i % 5 < 2 THEN NULL ELSE (i + 1)::UBIGINT END, 'A', 'C') AS g FROM range(1000) t(i));
SELECT CASE WHEN count(DISTINCT g) = 10 THEN true ELSE error('dictionary-shaped duplicates') END FROM (SELECT duckvep_allele_geometry((1 + (i % 10))::UBIGINT, 'A', 'G') AS g FROM range(10000) t(i));

-- An error in one row fails the statement, leaves the connection usable, and TRY absorbs it.
SELECT CASE WHEN TRY(duckvep_allele_geometry(0, 'A', 'G')) IS NULL THEN true ELSE error('TRY absorbs') END;
SELECT CASE WHEN TRY(duckvep_breakend_geometry('G]1:123[')) IS NULL THEN true ELSE error('TRY absorbs breakend') END;

-- Registered kinds.
SELECT CASE WHEN (SELECT count(*) FROM duckdb_functions() WHERE function_name = 'duckvep_so_terms' AND function_type = 'table') = 1 THEN true ELSE error('so_terms registered as a table function') END;
SELECT CASE WHEN (SELECT count(*) FROM duckdb_functions() WHERE function_name IN ('duckvep_allele_geometry', 'duckvep_breakend_geometry') AND function_type = 'scalar') = 2 THEN true ELSE error('geometry registered as scalar functions') END;

-- Nested output beyond one vector (2,048) and beyond the list-child capacity of one chunk.
SELECT CASE WHEN count(*) = 18000 AND sum(u.input_slot::BIGINT) = 27000 AND count(u.haplotype_lane) = 9000 THEN true ELSE error('phase output of 18,000 records') END FROM (SELECT unnest(duckvep_phase_call([i % 3, (i + 1) % 3], [false, i % 2 = 0])) AS u FROM range(9000) t(i));
SELECT CASE WHEN count(*) = 5000 AND max(u.ploidy) = 5000 AND min(u.input_slot) = 1 AND max(u.input_slot) = 5000 THEN true ELSE error('phase single 5,000-slot list') END FROM (SELECT unnest(duckvep_phase_call(list_transform(range(5000), lambda x: x % 3), list_transform(range(5000), lambda x: true))) AS u);
SELECT CASE WHEN length(g.reference) = 9998 AND g.reference_length = 9998 AND g.alternate_length = 3000 AND g.status = 'ok' THEN true ELSE error('repeat from 5,000-element lists') END FROM (SELECT duckvep_repeat_alleles(list_transform(range(5000), lambda x: {'unit': 'AC', 'count': x % 3}), list_transform(range(3000), lambda x: {'unit': 'G', 'count': 1}), true, {'max_allele_bases': 100000}) AS g);
SELECT CASE WHEN length(g.alternate) = 2000007 AND g.length_change = 7 THEN true ELSE error('repeat long output') END FROM (SELECT duckvep_repeat_alleles([{'unit': 'ACGTN', 'count': 400000}], [{'unit': 'ACGTN', 'count': 400001}, {'unit': 'RY', 'count': 1}], true, {'max_allele_bases': 2100000}) AS g);

-- NULLs inside lists and structs.
SELECT CASE WHEN g[2].allele_index IS NULL AND g[1].allele_index = 0 AND g[3].allele_index = 1 AND len(g) = 3 THEN true ELSE error('phase NULL allele element') END FROM (SELECT duckvep_phase_call([0, NULL, 1], [NULL, true, false]) AS g);
SELECT CASE WHEN g.status = 'incomplete_input' AND g.reference IS NULL AND g.length_change IS NULL THEN true ELSE error('repeat NULL unit') END FROM (SELECT duckvep_repeat_alleles([{'unit': NULL::VARCHAR, 'count': 3}], [{'unit': 'A', 'count': 4}], true) AS g);
SELECT CASE WHEN g.status = 'incomplete_input' THEN true ELSE error('repeat NULL record') END FROM (SELECT duckvep_repeat_alleles([NULL::STRUCT(unit VARCHAR, count INTEGER)], [{'unit': 'A', 'count': 4}], true) AS g);
SELECT CASE WHEN g IS NOT NULL AND g.status = 'incomplete_input' THEN true ELSE error('repeat untyped NULLs') END FROM (SELECT duckvep_repeat_alleles(NULL, NULL, true) AS g);
SELECT CASE WHEN count(*) = 6000 AND count(r) = 5454 AND count(r.status) = 5454 AND count(r.allele0) = 5454 THEN true ELSE error('raw_gt NULL structs have NULL children') END FROM (SELECT _duckvep_raw_gt(CASE WHEN i % 11 = 0 THEN NULL ELSE '0/1' END, 1::UINTEGER) AS r FROM range(6000) t(i));
SELECT CASE WHEN count(*) = 4000 AND count(s) = 3000 THEN true ELSE error('revcomp NULL rows') END FROM (SELECT _duckvep_revcomp(CASE WHEN i % 4 = 0 THEN NULL ELSE 'ACGT' END) AS s FROM range(4000) t(i));

-- Constant, dictionary and selected inputs agree with row-wise evaluation.
SELECT CASE WHEN (SELECT count(DISTINCT g) FROM (SELECT duckvep_phase_call([0, 1], [false, true]) AS g FROM range(5000))) = 1 THEN true ELSE error('phase constant input') END;
SELECT CASE WHEN bool_and(a = b) THEN true ELSE error('phase filtered equals row-wise') END FROM (SELECT duckvep_phase_call([i % 2, 1], [false, i % 4 = 0]) AS a, (SELECT duckvep_phase_call([i % 2, 1], [false, i % 4 = 0])) AS b FROM range(6000) t(i) WHERE i % 13 = 3);
SELECT CASE WHEN bool_and(a = b) THEN true ELSE error('repeat dictionary equals row-wise') END FROM (SELECT duckvep_repeat_alleles(d.r, d.a, true) AS a, (SELECT duckvep_repeat_alleles(d.r, d.a, true)) AS b FROM range(3000) k(i) JOIN (VALUES (0, [{'unit': 'CAG', 'count': 3}], [{'unit': 'CAG', 'count': 4}]), (1, [{'unit': 'A', 'count': 2}], [{'unit': 'T', 'count': 2}])) d(id, r, a) ON d.id = k.i % 2);

-- A failed call leaves the connection usable.
SELECT CASE WHEN TRY(duckvep_repeat_alleles([{'unit': 'CXG', 'count': 3}], [{'unit': 'A', 'count': 4}], true)) IS NULL THEN true ELSE error('TRY absorbs repeat error') END;
SELECT CASE WHEN TRY(duckvep_phase_call([0, -1], [false, true])) IS NULL THEN true ELSE error('TRY absorbs phase error') END;
SELECT CASE WHEN (SELECT count(*) FROM duckdb_functions() WHERE function_name IN ('duckvep_repeat_alleles', 'duckvep_phase_call', '_duckvep_revcomp', '_duckvep_raw_gt', '_duckvep_record_order' ) AND function_type = 'scalar') >= 5 THEN true ELSE error('nested scalars registered') END;

-- The twelve SQL builders are registered with and without the options argument.
SELECT CASE WHEN (SELECT count(DISTINCT function_name) FROM duckdb_functions() WHERE function_name IN ('duckvep_ensembl_regions_sql', 'duckvep_ensembl_transcripts_sql', 'duckvep_ensembl_regulation_features_sql', 'duckvep_model_receipt_sql', 'duckvep_annotate_sql', 'duckvep_annotate_projected_sql', 'duckvep_transcript_projection_sql', 'duckvep_prepare_sv_geometry_sql', 'duckvep_prepare_expansionhunter_sql', 'duckvep_prepare_breakend_pairs_sql', 'duckvep_prepare_breakend_fusion_sql', 'duckvep_prepare_structural_hgvs_sql')) = 12 THEN true ELSE error('builders registered') END;
SELECT CASE WHEN (SELECT count(*) FROM duckdb_functions() WHERE function_name LIKE 'duckvep_%_sql' AND function_type = 'scalar') = 26 THEN true ELSE error('two overloads per builder, plus load_sql') END;
