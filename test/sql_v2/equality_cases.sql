-- Queries run on BOTH hosts (v1 stable C API and v2). Each case is wrapped by the
-- runner as SELECT to_json(t) FROM (<query>) t, so results are compared as JSON
-- row text, in order; an error is compared by its duckvep_*: message.
-- Cases start at a "-- case: <name>" line and end at the next one.

-- case: so_terms ordered
SELECT * FROM duckvep_so_terms() ORDER BY bit_index

-- case: so_terms scan order
SELECT * FROM duckvep_so_terms()

-- case: so_terms filter and projection
SELECT consequence, impact FROM duckvep_so_terms() WHERE bit_index > 30 AND impact_code >= 0 ORDER BY 1

-- case: so_terms aggregate
SELECT count(*) AS n, sum(consequence_mask::HUGEINT) AS masks, max(severity_rank) AS worst, count(DISTINCT impact) AS impacts FROM duckvep_so_terms()

-- case: so_terms self join
SELECT count(*) AS n, sum(a.bit_index::INTEGER * b.bit_index::INTEGER) AS s FROM duckvep_so_terms() a, duckvep_so_terms() b

-- case: so_terms limit
SELECT bit_index FROM duckvep_so_terms() LIMIT 3

-- case: so_terms describe types
SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM duckvep_so_terms())

-- case: allele SNV
SELECT duckvep_allele_geometry(100, 'A', 'G') AS g

-- case: allele MNV
SELECT duckvep_allele_geometry(100, 'AC', 'GT') AS g

-- case: allele insertion (interbase)
SELECT duckvep_allele_geometry(100, 'A', 'AGT') AS g

-- case: allele deletion
SELECT duckvep_allele_geometry(100, 'ACG', 'A') AS g

-- case: allele indel
SELECT duckvep_allele_geometry(100, 'ACG', 'TTA') AS g

-- case: allele lowercase and N
SELECT duckvep_allele_geometry(7, 'acgtn', 'ACGTA') AS g

-- case: allele position one
SELECT duckvep_allele_geometry(1, 'A', 'AC') AS g

-- case: allele position max uinteger
SELECT duckvep_allele_geometry(4294967295, 'A', 'C') AS g

-- case: allele long non-inlined strings
SELECT duckvep_allele_geometry(5000, 'ACGTACGTACGTACGTACGTACGT', 'ACGTACGTACGTAAGTACGTACGT') AS g

-- case: allele 65535 bases ok
SELECT g.kind_code, g.reference_difference_length FROM (SELECT duckvep_allele_geometry(10, repeat('A', 65535), repeat('A', 65534) || 'C') AS g)

-- case: allele 65536 bases error
SELECT duckvep_allele_geometry(10, repeat('A', 65536), 'A') AS g

-- case: allele null arguments
SELECT duckvep_allele_geometry(NULL, 'A', 'G') AS a, duckvep_allele_geometry(1, NULL, 'G') AS b, duckvep_allele_geometry(1, 'A', NULL) AS c, duckvep_allele_geometry(NULL, NULL, NULL) AS d

-- case: allele null is a struct null
SELECT duckvep_allele_geometry(NULL, 'A', 'G') IS NULL AS is_null, (duckvep_allele_geometry(NULL, 'A', 'G')).kind_code IS NULL AS field_null

-- case: allele position zero error
SELECT duckvep_allele_geometry(0, 'A', 'G') AS g

-- case: allele position overflow error
SELECT duckvep_allele_geometry(4294967296, 'A', 'G') AS g

-- case: allele identical error
SELECT duckvep_allele_geometry(5, 'ACG', 'ACG') AS g

-- case: allele empty reference error
SELECT duckvep_allele_geometry(5, '', 'A') AS g

-- case: allele empty alternate error
SELECT duckvep_allele_geometry(5, 'A', '') AS g

-- case: allele bad base error
SELECT duckvep_allele_geometry(5, 'A', 'X') AS g

-- case: allele error under TRY
SELECT TRY(duckvep_allele_geometry(5, 'A', 'X')) AS g, TRY(duckvep_allele_geometry(5, 'A', 'G')).kind_code AS kind

-- case: allele constant vector
SELECT i, duckvep_allele_geometry(100, 'AC', 'A') AS g FROM range(5) t(i)

-- case: allele constant null
SELECT i, duckvep_allele_geometry(NULL, 'AC', 'A') AS g FROM range(3) t(i)

-- case: allele mixed null and valid rows
SELECT i, duckvep_allele_geometry(CASE WHEN i % 4 = 0 THEN NULL ELSE (100 + i)::UBIGINT END, 'AC', CASE WHEN i % 5 = 0 THEN NULL ELSE 'A' END) AS g FROM range(40) t(i)

-- case: allele filtered selection
SELECT i, duckvep_allele_geometry((i + 1)::UBIGINT, 'ACGT', 'A') AS g FROM range(6000) t(i) WHERE i % 97 = 3 ORDER BY i

-- case: allele dictionary join
SELECT k.i, duckvep_allele_geometry(d.pos, d.ref, d.alt) AS g FROM range(5000) k(i) JOIN (VALUES (1, 100::UBIGINT, 'A', 'AT'), (2, 200::UBIGINT, 'GG', 'G'), (3, 300::UBIGINT, 'C', 'T')) d(id, pos, ref, alt) ON d.id = k.i % 3 + 1 ORDER BY k.i

-- case: allele over several chunks checksum
SELECT count(*) AS rows, count(g) AS valid, sum(g.raw_start0) AS s0, sum(g.edit_end0) AS e, sum(g.kind_code::INTEGER) AS k, sum(g.alternate_difference_length::INTEGER) AS adl, count(g.insertion_boundary0) AS ib FROM (SELECT TRY(duckvep_allele_geometry(((i * 7919) % 1000000 + 1)::UBIGINT, substr('ACGTNACGTTGCAAGCTTAGGCTAACGGTAC', i % 7 + 1, i % 5 + 1), substr('TTGCAACCGGAATTCGATCGACGGCTA', i % 6 + 1, i % 4 + 1))) AS g FROM range(100000) t(i))

-- case: allele table input with NULLs and repeats
SELECT r.*, duckvep_allele_geometry(r.pos, r.ref, r.alt) AS g FROM (SELECT * FROM (VALUES (10::UBIGINT, 'A', 'C'), (NULL, 'A', 'C'), (10, 'AT', 'A'), (10, 'AT', 'A'), (11, NULL, 'C'), (12, 'TTTTTTTTTTTTTTTTT', 'TTTTTTTTTTTTTTTTTA')) v(pos, ref, alt)) r

-- case: breakend four bracket forms
SELECT a, duckvep_breakend_geometry(a) AS g FROM (VALUES ('G]1:123]'), (']1:123]G'), ('G[1:123['), ('[1:123[G')) v(a)

-- case: breakend single forms
SELECT a, duckvep_breakend_geometry(a) AS g FROM (VALUES ('.A'), ('A.'), ('.ACGTACGTACGTACGTACGT'), ('ACGTACGTACGTACGTACGT.')) v(a)

-- case: breakend non-BND is null
SELECT a, duckvep_breakend_geometry(a) AS g FROM (VALUES ('A'), ('<DEL>'), ('ACGT'), ('C')) v(a)

-- case: breakend null
SELECT duckvep_breakend_geometry(NULL) AS g, duckvep_breakend_geometry(NULL) IS NULL AS is_null

-- case: breakend long strings non-inlined
SELECT duckvep_breakend_geometry('GGACGTACGTACGTAC]very_long_contig_name_1:12345]') AS g

-- case: breakend telomeric mate
SELECT duckvep_breakend_geometry('G]1:0]') AS g

-- case: breakend chromosome with colons
SELECT duckvep_breakend_geometry('T]HLA-A*01:01:01:1:5]') AS g

-- case: breakend mismatched brackets error
SELECT duckvep_breakend_geometry('G]1:123[') AS g

-- case: breakend empty string error
SELECT duckvep_breakend_geometry('') AS g

-- case: breakend star allele is null
SELECT duckvep_breakend_geometry('*') AS g

-- case: breakend bad position error
SELECT duckvep_breakend_geometry('G]1:x]') AS g

-- case: breakend position overflow error
SELECT duckvep_breakend_geometry('G]1:99999999999999999999]') AS g

-- case: breakend bad replacement error
SELECT duckvep_breakend_geometry('X]1:5]') AS g

-- case: breakend error under TRY
SELECT TRY(duckvep_breakend_geometry('G]1:123[')) AS g, TRY(duckvep_breakend_geometry('G]1:123]')).mate_position AS p

-- case: breakend constant vector
SELECT i, duckvep_breakend_geometry('G]2:99]') AS g FROM range(5) t(i)

-- case: breakend mixed rows over several chunks
SELECT count(*) AS rows, count(g) AS valid, count(g.mate_chrom) AS mates, sum(g.mate_position) AS p, sum(g.local_join_after::INTEGER) AS j, sum(g.mate_extends_right::INTEGER) AS x, sum(length(g.replacement_sequence)) AS l FROM (SELECT duckvep_breakend_geometry(CASE i % 6 WHEN 0 THEN 'G]' || (i % 23) || ':' || i || ']' WHEN 1 THEN '[' || (i % 5) || ':' || i || '[AC' WHEN 2 THEN '.A' WHEN 3 THEN 'ACGT' WHEN 4 THEN NULL ELSE 'T.' END) AS g FROM range(50000) t(i))

-- case: breakend filtered selection
SELECT i, duckvep_breakend_geometry('C[' || i || ':' || (i * 3) || '[') AS g FROM range(5000) t(i) WHERE i % 211 = 7 ORDER BY i

-- case: geometry functions combined
SELECT duckvep_allele_geometry(10, 'A', 'C').kind_code AS k, duckvep_breakend_geometry('A]1:5]').mate_position AS p, (SELECT count(*) FROM duckvep_so_terms()) AS n
