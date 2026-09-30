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
