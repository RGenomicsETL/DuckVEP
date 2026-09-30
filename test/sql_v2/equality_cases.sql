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

-- case: phase diploid phased and unphased
SELECT duckvep_phase_call([0, 1], [false, true]) AS a, duckvep_phase_call([0, 1], [false, false]) AS b, duckvep_phase_call([1, 1], [false, true]) AS c

-- case: phase haploid and polyploid
SELECT duckvep_phase_call([1], [false]) AS h, duckvep_phase_call([0, 1, 2], [false, true, true]) AS t, duckvep_phase_call([0, 1, 1, 0], [false, true, false, true]) AS q

-- case: phase missing alleles and NULL flags
SELECT duckvep_phase_call([0, NULL, 1], [NULL, true, false]) AS a, duckvep_phase_call([NULL, NULL], [NULL, NULL]) AS b, duckvep_phase_call([0, 1], [NULL, NULL]) AS c

-- case: phase null lists
SELECT duckvep_phase_call(NULL, NULL) AS a, duckvep_phase_call(NULL::INTEGER[], NULL::BOOLEAN[]) AS b, duckvep_phase_call([0, 1], NULL::BOOLEAN[]) AS c, duckvep_phase_call([0, 1], NULL) AS d

-- case: phase null of lists typed NULL elements
SELECT duckvep_phase_call([NULL, NULL], [true, false]) AS a, duckvep_phase_call([0, 1], [NULL, NULL]) AS b

-- case: phase allele and flag types
SELECT duckvep_phase_call([0::TINYINT, 1::TINYINT], [0, 1]) AS a, duckvep_phase_call([0::UBIGINT, 1::UBIGINT], ['true', 'FALSE']) AS b, duckvep_phase_call(['0', '1'], ['1', '0']) AS c, duckvep_phase_call([0.0::DOUBLE, 1.0::DOUBLE], [1::BIGINT, 0::BIGINT]) AS d, duckvep_phase_call([0.5::DECIMAL(4,1), 1.0::DECIMAL(4,1)], [true, true]) AS e, duckvep_phase_call([0::HUGEINT, 1::HUGEINT], [true, false]) AS f

-- case: phase options phase_set and policy
SELECT duckvep_phase_call([0, 1], [false, true], {'phase_set': 5}) AS a, duckvep_phase_call([0, 1], [false, true], {'phase_set': 9::UBIGINT, 'phase_policy': 'strict'}) AS b, duckvep_phase_call([0, 1], [false, true], {'phase_policy': 'vep116_compat'}) AS c, duckvep_phase_call([0, 1], [true, true], {'phase_set': NULL::INTEGER, 'phase_policy': NULL::VARCHAR}) AS d, duckvep_phase_call([0, 1], [false, true], NULL) AS e

-- case: phase options with NULL struct rows
SELECT i, duckvep_phase_call([0, 1], [false, true], CASE WHEN i % 2 = 0 THEN NULL ELSE {'phase_set': i, 'phase_policy': 'vep116_compat'} END) AS g FROM range(6) t(i)

-- case: phase option phase_set too large error
SELECT duckvep_phase_call([0, 1], [false, true], {'phase_set': 18446744073709551615::UBIGINT}) AS g

-- case: phase bad policy error
SELECT duckvep_phase_call([0, 1], [false, true], {'phase_policy': 'loose'}) AS g

-- case: phase unknown option error
SELECT duckvep_phase_call([0, 1], [false, true], {'other': 1}) AS g

-- case: phase wrong option type error
SELECT duckvep_phase_call([0, 1], [false, true], {'phase_set': 'x'}) AS g

-- case: phase options not a struct error
SELECT duckvep_phase_call([0, 1], [false, true], 5) AS g

-- case: phase empty list error
SELECT duckvep_phase_call([]::INTEGER[], []::BOOLEAN[]) AS g

-- case: phase length mismatch error
SELECT duckvep_phase_call([0, 1], [true]) AS g

-- case: phase flags without alleles error
SELECT duckvep_phase_call(NULL::INTEGER[], [true]) AS g

-- case: phase negative allele error
SELECT duckvep_phase_call([0, -1], [false, true]) AS g

-- case: phase bad allele text error
SELECT duckvep_phase_call(['0', 'x'], [false, true]) AS g

-- case: phase bad flag text error
SELECT duckvep_phase_call([0, 1], ['maybe', 'true']) AS g

-- case: phase decimal flag error
SELECT duckvep_phase_call([0, 1], [1.0::DECIMAL(3,1), 0.0::DECIMAL(3,1)]) AS g

-- case: phase ploidy over limit error
SELECT duckvep_phase_call(list_transform(range(65536), lambda x: 0), list_transform(range(65536), lambda x: false)) AS g

-- case: phase error under TRY
SELECT TRY(duckvep_phase_call([0, -1], [false, true])) AS a, len(duckvep_phase_call([0, 1], [false, true])) AS b

-- case: phase one row beyond one vector
SELECT count(*) AS n, sum(u.input_slot::BIGINT) AS slots, count(u.allele_index) AS called, count(u.haplotype_lane) AS lanes, max(u.ploidy) AS ploidy, count(DISTINCT u.phase_scope) AS scopes FROM (SELECT unnest(duckvep_phase_call(list_transform(range(5000), lambda x: x % 3), list_transform(range(5000), lambda x: x % 2 = 0))) AS u)

-- case: phase many rows beyond one vector
SELECT count(*) AS n, sum(u.input_slot::BIGINT) AS slots, sum(u.allele_index::BIGINT) AS alleles, count(u.haplotype_lane) AS lanes, sum(u.phase_set) AS sets, count(DISTINCT u.status) AS statuses FROM (SELECT unnest(duckvep_phase_call([i % 3, (i + 1) % 3], [false, i % 2 = 0], {'phase_set': i})) AS u FROM range(9000) t(i))

-- case: phase rows with mixed ploidy NULLs and errors avoided
SELECT count(*) AS rows, count(g) AS valid, sum(len(g)) AS slots, sum(len(list_filter(g, lambda x: x.allele_index IS NULL))) AS missing FROM (SELECT duckvep_phase_call(CASE WHEN i % 7 = 0 THEN NULL WHEN i % 3 = 0 THEN [i % 2] WHEN i % 3 = 1 THEN [i % 2, NULL] ELSE [0, 1, i % 2, NULL] END, CASE WHEN i % 5 = 0 OR i % 7 = 0 THEN NULL WHEN i % 3 = 0 THEN [false] WHEN i % 3 = 1 THEN [true, NULL] ELSE [false, true, true, false] END) AS g FROM range(7000) t(i))

-- case: phase constant lists
SELECT i, duckvep_phase_call([0, 1], [false, true]) AS g FROM range(4) t(i)

-- case: phase filtered selection
SELECT i, duckvep_phase_call([i % 2, 1], [false, i % 4 = 0]) AS g FROM range(6000) t(i) WHERE i % 499 = 3 ORDER BY i

-- case: phase dictionary join
SELECT k.i, duckvep_phase_call(d.gt, d.ph) AS g FROM range(3000) k(i) JOIN (VALUES (0, [0, 1], [false, true]), (1, [1, 1], [false, false]), (2, [0], [false])) d(id, gt, ph) ON d.id = k.i % 3 ORDER BY k.i LIMIT 30

-- case: repeat basic gain loss neutral
SELECT duckvep_repeat_alleles([{'unit': 'CAG', 'count': 3}], [{'unit': 'CAG', 'count': 5}], true) AS gain, duckvep_repeat_alleles([{'unit': 'CAG', 'count': 5}], [{'unit': 'CAG', 'count': 3}], true) AS loss, duckvep_repeat_alleles([{'unit': 'CAG', 'count': 3}], [{'unit': 'CAG', 'count': 3}], true) AS neutral

-- case: repeat multiple units and field order
SELECT duckvep_repeat_alleles([{'unit': 'CAG', 'count': 2}, {'unit': 'CCG', 'count': 1}], [{'count': 4, 'unit': 'ac'}], true) AS g

-- case: repeat empty lists
SELECT duckvep_repeat_alleles([]::STRUCT(unit VARCHAR, count INTEGER)[], [{'unit': 'A', 'count': 2}], true) AS g

-- case: repeat count types
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3::TINYINT}], [{'unit': 'A', 'count': 4::UBIGINT}], true) AS a, duckvep_repeat_alleles([{'unit': 'A', 'count': 3::HUGEINT}], [{'unit': 'A', 'count': 4.0::DOUBLE}], true) AS b, duckvep_repeat_alleles([{'unit': 'A', 'count': 3::DECIMAL(6,2)}], [{'unit': 'A', 'count': 4::DECIMAL(30,3)}], true) AS c, duckvep_repeat_alleles([{'unit': 'A', 'count': 3::FLOAT}], [{'unit': 'A', 'count': 4::DECIMAL(2,0)}], true) AS d

-- case: repeat nonintegral counts
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3.5::DOUBLE}], [{'unit': 'A', 'count': 4}], true) AS a, duckvep_repeat_alleles([{'unit': 'A', 'count': 3.5::DECIMAL(5,1)}], [{'unit': 'A', 'count': 4}], true) AS b, duckvep_repeat_alleles([{'unit': 'A', 'count': 7::DECIMAL(30,2) / 2}], [{'unit': 'A', 'count': 4}], true) AS c

-- case: repeat summary only and incomplete
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], false) AS a, duckvep_repeat_alleles([{'unit': 'A', 'count': NULL::INTEGER}], [{'unit': 'A', 'count': 4}], true) AS b, duckvep_repeat_alleles([{'unit': NULL::VARCHAR, 'count': 3}], [{'unit': 'A', 'count': 4}], true) AS c, duckvep_repeat_alleles([NULL::STRUCT(unit VARCHAR, count INTEGER)], [{'unit': 'A', 'count': 4}], true) AS d, duckvep_repeat_alleles(NULL::STRUCT(unit VARCHAR, count INTEGER)[], [{'unit': 'A', 'count': 4}], true) AS e

-- case: repeat untyped NULL lists
SELECT duckvep_repeat_alleles(NULL, NULL, true) AS a, duckvep_repeat_alleles([NULL, NULL], [NULL], true) AS b, duckvep_repeat_alleles(NULL, [{'unit': 'A', 'count': 1}], false) AS c

-- case: repeat options
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 6000}], [{'unit': 'A', 'count': 6001}], true, {'max_allele_bases': 7000}) AS a, duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, {'max_allele_bases': 4.0::DECIMAL(5,1)}) AS b, duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, NULL) AS c, duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], false, {'max_allele_bases': 1::UTINYINT}) AS d

-- case: repeat null struct option rows
SELECT i, duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, CASE WHEN i % 2 = 0 THEN NULL ELSE {'max_allele_bases': 10} END) AS g FROM range(4) t(i)

-- case: repeat exceeds default cap error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 5001}], true) AS g

-- case: repeat exceeds custom cap error
SELECT duckvep_repeat_alleles([{'unit': 'CAG', 'count': 40}], [{'unit': 'A', 'count': 1}], true, {'max_allele_bases': 100}) AS g

-- case: repeat bad cap fractional error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, {'max_allele_bases': 1.5}) AS g

-- case: repeat bad cap negative error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, {'max_allele_bases': -1}) AS g

-- case: repeat bad cap null error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, {'max_allele_bases': NULL}) AS g

-- case: repeat unknown option error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, {'limit': 10}) AS g

-- case: repeat wrong option type error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, {'max_allele_bases': 'ten'}) AS g

-- case: repeat options not a struct error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true, 10) AS g

-- case: repeat sequence exact required error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], NULL::BOOLEAN) AS g

-- case: repeat invalid unit error
SELECT duckvep_repeat_alleles([{'unit': 'CXG', 'count': 3}], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat empty unit error
SELECT duckvep_repeat_alleles([{'unit': '', 'count': 3}], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat negative count error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': -1}], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat infinite count error
SELECT duckvep_repeat_alleles([{'unit': 'A', 'count': 'Infinity'::DOUBLE}], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat not lists error
SELECT duckvep_repeat_alleles([1, 2], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat wrong struct names error
SELECT duckvep_repeat_alleles([{'base': 'A', 'count': 3}], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat scalar not a list error
SELECT duckvep_repeat_alleles('ACG', [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat unit not text error
SELECT duckvep_repeat_alleles([{'unit': 1, 'count': 3}], [{'unit': 'A', 'count': 4}], true) AS g

-- case: repeat error under TRY
SELECT TRY(duckvep_repeat_alleles([{'unit': 'CXG', 'count': 3}], [{'unit': 'A', 'count': 4}], true)) AS a, duckvep_repeat_alleles([{'unit': 'A', 'count': 1}], [{'unit': 'A', 'count': 2}], true).status AS b

-- case: repeat long output beyond inline strings
SELECT length(g.reference) AS r, length(g.alternate) AS a, g.length_change AS c, md5(g.alternate) AS h FROM (SELECT duckvep_repeat_alleles([{'unit': 'ACGTN', 'count': 400000}], [{'unit': 'ACGTN', 'count': 400001}, {'unit': 'RY', 'count': 7}], true, {'max_allele_bases': 2100000}) AS g)

-- case: repeat one row beyond one vector of elements
SELECT length(g.reference) AS r, g.reference_length AS n, md5(g.reference) AS h, g.status AS s FROM (SELECT duckvep_repeat_alleles(list_transform(range(5000), lambda x: {'unit': 'AC', 'count': x % 3}), list_transform(range(3000), lambda x: {'unit': 'G', 'count': 1}), true, {'max_allele_bases': 100000}) AS g)

-- case: repeat many rows beyond one vector
SELECT count(*) AS rows, count(g.reference) AS ok, sum(g.reference_length) AS r, sum(g.alternate_length) AS a, sum(g.length_change) AS c, count(DISTINCT g.length_direction) AS d, count(DISTINCT g.status) AS s FROM (SELECT duckvep_repeat_alleles([{'unit': 'CAG', 'count': i % 11}], [{'unit': 'CAG', 'count': i % 13}, {'unit': 'T', 'count': i % 2}], i % 17 <> 0) AS g FROM range(9000) t(i))

-- case: repeat mixed null and incomplete rows
SELECT i, duckvep_repeat_alleles(CASE WHEN i % 4 = 0 THEN NULL ELSE [{'unit': 'A', 'count': i}] END, [{'unit': CASE WHEN i % 5 = 0 THEN NULL ELSE 'T' END, 'count': 2}], i % 3 <> 0) AS g FROM range(30) t(i)

-- case: repeat constant lists
SELECT i, duckvep_repeat_alleles([{'unit': 'CAG', 'count': 3}], [{'unit': 'CAG', 'count': 5}], true) AS g FROM range(4) t(i)

-- case: repeat filtered selection
SELECT i, duckvep_repeat_alleles([{'unit': 'CAG', 'count': i % 7}], [{'unit': 'CAG', 'count': i % 5}], true).status AS s FROM range(6000) t(i) WHERE i % 613 = 5 ORDER BY i

-- case: repeat dictionary join
SELECT k.i, duckvep_repeat_alleles(d.r, d.a, true) AS g FROM range(3000) k(i) JOIN (VALUES (0, [{'unit': 'CAG', 'count': 3}], [{'unit': 'CAG', 'count': 4}]), (1, [{'unit': 'A', 'count': 2}], [{'unit': 'T', 'count': 2}]), (2, [{'unit': 'G', 'count': 1}], [{'unit': 'G', 'count': 1}])) d(id, r, a) ON d.id = k.i % 3 ORDER BY k.i LIMIT 30

-- case: revcomp bases and case
SELECT s, _duckvep_revcomp(s) AS r FROM (VALUES ('ACGT'), ('acgtn'), ('RYSWKMBDHVN'), ('ryswkmbdhvn'), ('ACGU'), (''), ('A'), (NULL)) v(s)

-- case: revcomp utf8 kept
SELECT _duckvep_revcomp('ACé€TG') AS a, _duckvep_revcomp('日本AC') AS b, _duckvep_revcomp('AC' || chr(128512) || 'GT') AS c

-- case: revcomp long inputs
SELECT length(r) AS n, md5(r) AS h FROM (SELECT _duckvep_revcomp(repeat('ACGTNRY', 30000) || 'AAC') AS r)

-- case: revcomp many rows
SELECT count(*) AS n, count(r) AS ok, sum(length(r)) AS total, md5(string_agg(r, '' ORDER BY i)) AS h FROM (SELECT i, _duckvep_revcomp(CASE WHEN i % 9 = 0 THEN NULL ELSE repeat('ACGT', i % 7) || substr('NRYSW', i % 5 + 1, 1) END) AS r FROM range(6000) t(i))

-- case: revcomp constant and filtered
SELECT i, _duckvep_revcomp('GATTACA') AS c, _duckvep_revcomp(repeat('AC', i % 10)) AS f FROM range(3000) t(i) WHERE i % 401 = 7 ORDER BY i

-- case: raw_gt forms
SELECT g, n, _duckvep_raw_gt(g, n) AS r FROM (VALUES ('0/1', 1::UINTEGER), ('1|0', 1), ('0|1', 2), ('.', 1), ('./.', 1), ('1/2', 2), ('|1', 1), ('0/0', 1), ('1', 3), ('2/2', 2), ('1/1/1', 1)) v(g, n)

-- case: raw_gt null inputs
SELECT _duckvep_raw_gt(NULL, 1::UINTEGER) AS a, _duckvep_raw_gt('0/1', NULL) AS b, _duckvep_raw_gt('0/1', NULL) IS NULL AS c, (_duckvep_raw_gt(NULL, 1::UINTEGER)).status IS NULL AS d

-- case: raw_gt many rows
SELECT count(*) AS n, count(r) AS ok, sum(r.status::BIGINT) AS s, sum(r.allele0::BIGINT) AS a0, sum(r.allele1::BIGINT) AS a1, sum(r.parsed_slots::BIGINT) AS p, sum(r.disposition::BIGINT) AS d FROM (SELECT _duckvep_raw_gt(CASE WHEN i % 11 = 0 THEN NULL ELSE (i % 3)::VARCHAR || CASE i % 4 WHEN 0 THEN '/' WHEN 1 THEN '|' ELSE '/' END || (i % 2)::VARCHAR END, (1 + i % 3)::UINTEGER) AS r FROM range(6000) t(i))

-- case: record_order values
SELECT c, o, _duckvep_record_order(c, o) AS r FROM (VALUES (1::UBIGINT, 1::UBIGINT), (2, 1), (2, 2), (5, 3), (100, 50), (NULL, 1), (3, NULL)) v(c, o)

-- case: record_order invalid ordinal error
SELECT _duckvep_record_order(2::UBIGINT, 3::UBIGINT) AS r

-- case: record_order zero ordinal error
SELECT _duckvep_record_order(2::UBIGINT, 0::UBIGINT) AS r

-- case: record_order many rows
SELECT count(*) AS n, sum(r) AS s, max(r) AS m FROM (SELECT _duckvep_record_order((i % 50 + 1)::UBIGINT, (i % 50 + 1)::UBIGINT) AS r FROM range(5000) t(i))

-- case: builder annotate defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('events', 'model') AS x)

-- case: builder annotate schema qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('main.events', 'model') AS x)

-- case: builder annotate quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('ev"il', 'mo''del') AS x)

-- case: builder annotate quoting schema
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('s"1.t"2', 'm') AS x)

-- case: builder annotate extra dots
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('a.b.c', 'm') AS x)

-- case: builder annotate unicode
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('événements', 'modèle ✓') AS x)

-- case: builder annotate empty events
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('', 'm') AS x)

-- case: builder annotate long names
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql(repeat('x', 300), repeat('y', 300)) AS x)

-- case: builder annotate distances
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'upstream_distance': 10, 'downstream_distance': 20}) AS x)

-- case: builder annotate distance types
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'upstream_distance': 3::TINYINT, 'downstream_distance': 18446744073709551615::UBIGINT}) AS x)

-- case: builder annotate distance negative
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'upstream_distance': -1::BIGINT}) AS x)

-- case: builder annotate distance null
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'upstream_distance': NULL::INTEGER, 'downstream_distance': NULL}) AS x)

-- case: builder annotate null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', NULL) AS x)

-- case: builder annotate reversed keys
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'downstream_distance': 2, 'upstream_distance': 1}) AS x)

-- case: builder annotate unknown option error
SELECT duckvep_annotate_sql('e', 'm', {'other': 1}) AS x

-- case: builder annotate distance wrong type error
SELECT duckvep_annotate_sql('e', 'm', {'upstream_distance': 'x'}) AS x

-- case: builder annotate distance double type error
SELECT duckvep_annotate_sql('e', 'm', {'upstream_distance': 1.5}) AS x

-- case: builder annotate options not a struct error
SELECT duckvep_annotate_sql('e', 'm', 5) AS x

-- case: builder annotate null events error
SELECT duckvep_annotate_sql(NULL, 'm') AS x

-- case: builder annotate null model error
SELECT duckvep_annotate_sql('e', NULL) AS x

-- case: builder annotate null events with options error
SELECT duckvep_annotate_sql(NULL, 'm', {'nope': 1}) AS x

-- case: builder annotate embedded nul error
SELECT duckvep_annotate_sql('e' || chr(0) || 'x', 'm') AS x

-- case: builder projected defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('events', 'model') AS x)

-- case: builder projected schema qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('main.events', 'model') AS x)

-- case: builder projected quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('ev"il', 'mo''del') AS x)

-- case: builder projected quoting schema
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('s"1.t"2', 'm') AS x)

-- case: builder projected extra dots
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('a.b.c', 'm') AS x)

-- case: builder projected unicode
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('événements', 'modèle ✓') AS x)

-- case: builder projected empty events
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('', 'm') AS x)

-- case: builder projected long names
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql(repeat('x', 300), repeat('y', 300)) AS x)

-- case: builder projected distances
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('e', 'm', {'upstream_distance': 10, 'downstream_distance': 20}) AS x)

-- case: builder projected distance types
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('e', 'm', {'upstream_distance': 3::TINYINT, 'downstream_distance': 18446744073709551615::UBIGINT}) AS x)

-- case: builder projected distance negative
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('e', 'm', {'upstream_distance': -1::BIGINT}) AS x)

-- case: builder projected distance null
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('e', 'm', {'upstream_distance': NULL::INTEGER, 'downstream_distance': NULL}) AS x)

-- case: builder projected null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('e', 'm', NULL) AS x)

-- case: builder projected reversed keys
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_projected_sql('e', 'm', {'downstream_distance': 2, 'upstream_distance': 1}) AS x)

-- case: builder projected unknown option error
SELECT duckvep_annotate_projected_sql('e', 'm', {'other': 1}) AS x

-- case: builder projected distance wrong type error
SELECT duckvep_annotate_projected_sql('e', 'm', {'upstream_distance': 'x'}) AS x

-- case: builder projected distance double type error
SELECT duckvep_annotate_projected_sql('e', 'm', {'upstream_distance': 1.5}) AS x

-- case: builder projected options not a struct error
SELECT duckvep_annotate_projected_sql('e', 'm', 5) AS x

-- case: builder projected null events error
SELECT duckvep_annotate_projected_sql(NULL, 'm') AS x

-- case: builder projected null model error
SELECT duckvep_annotate_projected_sql('e', NULL) AS x

-- case: builder projected null events with options error
SELECT duckvep_annotate_projected_sql(NULL, 'm', {'nope': 1}) AS x

-- case: builder projected embedded nul error
SELECT duckvep_annotate_projected_sql('e' || chr(0) || 'x', 'm') AS x

-- case: builder annotate hgvs and rich
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'hgvs': true, 'rich': true}) AS x)

-- case: builder annotate all options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'hgvs': false, 'upstream_distance': 1, 'downstream_distance': 2, 'rich': false}) AS x)

-- case: builder annotate null booleans
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', 'm', {'hgvs': NULL::BOOLEAN, 'rich': NULL}) AS x)

-- case: builder annotate hgvs wrong type error
SELECT duckvep_annotate_sql('e', 'm', {'hgvs': 1}) AS x

-- case: builder annotate rich wrong type error
SELECT duckvep_annotate_sql('e', 'm', {'rich': 'yes'}) AS x

-- case: builder projected hgvs unsupported error
SELECT duckvep_annotate_projected_sql('e', 'm', {'hgvs': true}) AS x

-- case: builder projected empty model error
SELECT duckvep_annotate_projected_sql('e', '') AS x

-- case: builder annotate empty model
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_annotate_sql('e', '') AS x)

-- case: builder annotate many rows
SELECT count(*) AS n, count(DISTINCT x) AS d, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_annotate_sql('t' || (i % 37), 'm' || (i % 5), {'upstream_distance': i % 11, 'rich': i % 2 = 0}) AS x FROM range(3000) t(i))

-- case: builder projected many rows
SELECT count(*) AS n, count(DISTINCT x) AS d, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_annotate_projected_sql('t' || (i % 37), 'm' || (i % 5), {'downstream_distance': i % 13}) AS x FROM range(3000) t(i))

-- case: builder annotate filtered rows
SELECT i, md5(duckvep_annotate_sql('t' || i, 'm')) AS h FROM range(5000) t(i) WHERE i % 401 = 9 ORDER BY i

-- case: builder annotate constant rows
SELECT i, md5(duckvep_annotate_sql('t', 'm')) AS h FROM range(4) t(i)

-- case: builder annotate null row among many error
SELECT duckvep_annotate_sql(CASE WHEN i = 1500 THEN NULL ELSE 't' END, 'm') AS x FROM range(3000) t(i)

-- case: builder projection defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_transcript_projection_sql('events', 'annotations', 'transcripts') AS x)

-- case: builder projection schema qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_transcript_projection_sql('a.events', 'b.annotations', 'c.transcripts') AS x)

-- case: builder projection quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_transcript_projection_sql('e"1', 'a"2', 't"3') AS x)

-- case: builder projection extra dots
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_transcript_projection_sql('a.b.c', 'x', 'y') AS x)

-- case: builder projection null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_transcript_projection_sql('e', 'a', 't', NULL) AS x)

-- case: builder projection empty events error
SELECT duckvep_transcript_projection_sql('', 'a', 't') AS x

-- case: builder projection empty transcripts error
SELECT duckvep_transcript_projection_sql('e', 'a', '') AS x

-- case: builder projection null annotations error
SELECT duckvep_transcript_projection_sql('e', NULL, 't') AS x

-- case: builder projection option not supported error
SELECT duckvep_transcript_projection_sql('e', 'a', 't', {'x': 1}) AS x

-- case: builder projection options not a struct error
SELECT duckvep_transcript_projection_sql('e', 'a', 't', 'x') AS x

-- case: builder projection many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_transcript_projection_sql('e' || (i % 9), 'a' || (i % 4), 't') AS x FROM range(2500) t(i))

-- case: builder prepare sv defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_sv_geometry_sql('events') AS x)

-- case: builder prepare sv schema
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_sv_geometry_sql('main.events') AS x)

-- case: builder prepare sv quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_sv_geometry_sql('e"v') AS x)

-- case: builder prepare sv null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_sv_geometry_sql('events', NULL) AS x)

-- case: builder prepare sv null relation error
SELECT duckvep_prepare_sv_geometry_sql(NULL) AS x

-- case: builder prepare sv empty relation error
SELECT duckvep_prepare_sv_geometry_sql('') AS x

-- case: builder prepare sv extra dots error
SELECT duckvep_prepare_sv_geometry_sql('a.b.c') AS x

-- case: builder prepare sv leading dot error
SELECT duckvep_prepare_sv_geometry_sql('.x') AS x

-- case: builder prepare sv trailing dot error
SELECT duckvep_prepare_sv_geometry_sql('x.') AS x

-- case: builder prepare sv option unsupported error
SELECT duckvep_prepare_sv_geometry_sql('e', {'x': 1}) AS x

-- case: builder prepare sv options not a struct error
SELECT duckvep_prepare_sv_geometry_sql('e', 1) AS x

-- case: builder prepare expansionhunter defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_expansionhunter_sql('events', 'reference') AS x)

-- case: builder prepare expansionhunter schemas
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_expansionhunter_sql('a.events', 'b.reference') AS x)

-- case: builder prepare expansionhunter null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_expansionhunter_sql('e', 'r', NULL) AS x)

-- case: builder prepare expansionhunter null reference error
SELECT duckvep_prepare_expansionhunter_sql('e', NULL) AS x

-- case: builder prepare expansionhunter empty events error
SELECT duckvep_prepare_expansionhunter_sql('', 'r') AS x

-- case: builder prepare expansionhunter extra dots error
SELECT duckvep_prepare_expansionhunter_sql('e', 'a.b.c') AS x

-- case: builder prepare expansionhunter option unsupported error
SELECT duckvep_prepare_expansionhunter_sql('e', 'r', {'x': 1}) AS x

-- case: builder prepare sv many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_prepare_sv_geometry_sql('t' || (i % 7)) AS x FROM range(2500) t(i))

-- case: builder prepare expansionhunter many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_prepare_expansionhunter_sql('e' || i, 'r' || (i % 3)) AS x FROM range(2500) t(i))

-- case: builder breakend_pairs defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_pairs_sql('rel1') AS x)

-- case: builder breakend_pairs schema
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_pairs_sql('s.r1') AS x)

-- case: builder breakend_pairs quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_pairs_sql('q"1') AS x)

-- case: builder breakend_pairs null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_pairs_sql('rel1', NULL) AS x)

-- case: builder breakend_pairs null relation error
SELECT duckvep_prepare_breakend_pairs_sql(NULL) AS x

-- case: builder breakend_pairs empty relation error
SELECT duckvep_prepare_breakend_pairs_sql('') AS x

-- case: builder breakend_pairs extra dots error
SELECT duckvep_prepare_breakend_pairs_sql('a.b.c') AS x

-- case: builder breakend_pairs option unsupported error
SELECT duckvep_prepare_breakend_pairs_sql('rel1', {'max_span': 10}) AS x

-- case: builder breakend_pairs options not a struct error
SELECT duckvep_prepare_breakend_pairs_sql('rel1', 3) AS x

-- case: builder breakend_fusion defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_fusion_sql('rel1', 'rel2') AS x)

-- case: builder breakend_fusion schema
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_fusion_sql('s.r1', 's.r2') AS x)

-- case: builder breakend_fusion quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_fusion_sql('q"1', 'q"2') AS x)

-- case: builder breakend_fusion null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_breakend_fusion_sql('rel1', 'rel2', NULL) AS x)

-- case: builder breakend_fusion null relation error
SELECT duckvep_prepare_breakend_fusion_sql(NULL, 'r') AS x

-- case: builder breakend_fusion empty relation error
SELECT duckvep_prepare_breakend_fusion_sql('r', '') AS x

-- case: builder breakend_fusion extra dots error
SELECT duckvep_prepare_breakend_fusion_sql('a.b.c', 'r') AS x

-- case: builder breakend_fusion option unsupported error
SELECT duckvep_prepare_breakend_fusion_sql('rel1', 'rel2', {'max_span': 10}) AS x

-- case: builder breakend_fusion options not a struct error
SELECT duckvep_prepare_breakend_fusion_sql('rel1', 'rel2', 3) AS x

-- case: builder structural_hgvs defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('rel1', 'rel2') AS x)

-- case: builder structural_hgvs schema
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('s.r1', 's.r2') AS x)

-- case: builder structural_hgvs quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('q"1', 'q"2') AS x)

-- case: builder structural_hgvs null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('rel1', 'rel2', NULL) AS x)

-- case: builder structural_hgvs null relation error
SELECT duckvep_prepare_structural_hgvs_sql(NULL, 'r') AS x

-- case: builder structural_hgvs empty relation error
SELECT duckvep_prepare_structural_hgvs_sql('r', '') AS x

-- case: builder structural_hgvs extra dots error
SELECT duckvep_prepare_structural_hgvs_sql('a.b.c', 'r') AS x

-- case: builder structural_hgvs max_span 1
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 1}) AS x)

-- case: builder structural_hgvs max_span 5000
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 5000}) AS x)

-- case: builder structural_hgvs max_span 60000
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 60000}) AS x)

-- case: builder structural_hgvs max_span types
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 77::UTINYINT}) AS x)

-- case: builder structural_hgvs max_span null
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': NULL::INTEGER}) AS x)

-- case: builder structural_hgvs max_span zero error
SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 0}) AS x

-- case: builder structural_hgvs max_span too large error
SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 60001}) AS x

-- case: builder structural_hgvs max_span huge unsigned error
SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 18446744073709551615::UBIGINT}) AS x

-- case: builder structural_hgvs max_span negative error
SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': -3}) AS x

-- case: builder structural_hgvs max_span wrong type error
SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'max_span': 'x'}) AS x

-- case: builder structural_hgvs unknown option error
SELECT duckvep_prepare_structural_hgvs_sql('a', 'b', {'other': 1}) AS x

-- case: builder structural many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_prepare_structural_hgvs_sql('a' || (i % 5), 'b' || (i % 3), {'max_span': 1 + i % 60000}) AS x FROM range(2600) t(i))

-- case: builder ensembl regions defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'reference_chunks', 'GRCh38') AS x)

-- case: builder ensembl regions schema qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core_1', 'sch.chunks', 'GRCh37') AS x)

-- case: builder ensembl regions quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('co"re', 'chunks"', 'asm''bly') AS x)

-- case: builder ensembl regions species 1
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 1}) AS x)

-- case: builder ensembl regions species 0
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 0}) AS x)

-- case: builder ensembl regions species 7
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 7}) AS x)

-- case: builder ensembl regions species m2BIGINT
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': -2::BIGINT}) AS x)

-- case: builder ensembl regions species 255UTINYINT
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 255::UTINYINT}) AS x)

-- case: builder ensembl regions species 9223372036854775807UBIGI
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 9223372036854775807::UBIGINT}) AS x)

-- case: builder ensembl regions species NULLINTEGER
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': NULL::INTEGER}) AS x)

-- case: builder ensembl regions species NULL
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': NULL}) AS x)

-- case: builder ensembl regions null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', NULL) AS x)

-- case: builder ensembl regions species too large error
SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 18446744073709551615::UBIGINT}) AS x

-- case: builder ensembl regions species wrong type error
SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species_id': 'x'}) AS x

-- case: builder ensembl regions unknown option error
SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', {'species': 1}) AS x

-- case: builder ensembl regions options not a struct error
SELECT duckvep_ensembl_regions_sql('core', 'ref', 'asm', 1) AS x

-- case: builder ensembl regions null core error
SELECT duckvep_ensembl_regions_sql(NULL, 'ref', 'asm') AS x

-- case: builder ensembl regions empty assembly error
SELECT duckvep_ensembl_regions_sql('core', 'ref', '') AS x

-- case: builder ensembl regions empty reference error
SELECT duckvep_ensembl_regions_sql('core', '', 'asm') AS x

-- case: builder ensembl regions bad qualified reference error
SELECT duckvep_ensembl_regions_sql('core', '.x', 'asm') AS x

-- case: builder ensembl regions trailing dot reference error
SELECT duckvep_ensembl_regions_sql('core', 'x.', 'asm') AS x

-- case: builder ensembl regions many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_ensembl_regions_sql('c' || (i % 4), 'r' || (i % 3), 'a', {'species_id': i % 5}) AS x FROM range(1200) t(i))

-- case: builder ensembl transcripts defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'reference_chunks', 'GRCh38') AS x)

-- case: builder ensembl transcripts schema qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core_1', 'sch.chunks', 'GRCh37') AS x)

-- case: builder ensembl transcripts quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('co"re', 'chunks"', 'asm''bly') AS x)

-- case: builder ensembl transcripts species 1
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 1}) AS x)

-- case: builder ensembl transcripts species 0
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 0}) AS x)

-- case: builder ensembl transcripts species 7
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 7}) AS x)

-- case: builder ensembl transcripts species m2BIGINT
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': -2::BIGINT}) AS x)

-- case: builder ensembl transcripts species 255UTINYINT
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 255::UTINYINT}) AS x)

-- case: builder ensembl transcripts species 9223372036854775807UBIGI
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 9223372036854775807::UBIGINT}) AS x)

-- case: builder ensembl transcripts species NULLINTEGER
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': NULL::INTEGER}) AS x)

-- case: builder ensembl transcripts species NULL
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': NULL}) AS x)

-- case: builder ensembl transcripts null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', NULL) AS x)

-- case: builder ensembl transcripts species too large error
SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 18446744073709551615::UBIGINT}) AS x

-- case: builder ensembl transcripts species wrong type error
SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species_id': 'x'}) AS x

-- case: builder ensembl transcripts unknown option error
SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', {'species': 1}) AS x

-- case: builder ensembl transcripts options not a struct error
SELECT duckvep_ensembl_transcripts_sql('core', 'ref', 'asm', 1) AS x

-- case: builder ensembl transcripts null core error
SELECT duckvep_ensembl_transcripts_sql(NULL, 'ref', 'asm') AS x

-- case: builder ensembl transcripts empty assembly error
SELECT duckvep_ensembl_transcripts_sql('core', 'ref', '') AS x

-- case: builder ensembl transcripts empty reference error
SELECT duckvep_ensembl_transcripts_sql('core', '', 'asm') AS x

-- case: builder ensembl transcripts bad qualified reference error
SELECT duckvep_ensembl_transcripts_sql('core', '.x', 'asm') AS x

-- case: builder ensembl transcripts trailing dot reference error
SELECT duckvep_ensembl_transcripts_sql('core', 'x.', 'asm') AS x

-- case: builder ensembl transcripts many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_ensembl_transcripts_sql('c' || (i % 4), 'r' || (i % 3), 'a', {'species_id': i % 5}) AS x FROM range(1200) t(i))

-- case: builder ensembl regulation defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regulation_features_sql('funcgen', 'regions') AS x)

-- case: builder ensembl regulation schema regions
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regulation_features_sql('funcgen_1', 'sch.regions') AS x)

-- case: builder ensembl regulation quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_ensembl_regulation_features_sql('fun"c', 'reg"') AS x)

-- case: builder ensembl regulation null schema error
SELECT duckvep_ensembl_regulation_features_sql(NULL, 'regions') AS x

-- case: builder ensembl regulation empty regions error
SELECT duckvep_ensembl_regulation_features_sql('f', '') AS x

-- case: builder ensembl regulation options unsupported error
SELECT duckvep_ensembl_regulation_features_sql('f', 'r', {'species_id': 1}) AS x

-- case: builder ensembl regulation null options unsupported error
SELECT duckvep_ensembl_regulation_features_sql('f', 'r', NULL) AS x

-- case: builder ensembl regulation many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_ensembl_regulation_features_sql('f' || (i % 4), 'r' || (i % 3)) AS x FROM range(1200) t(i))

-- case: builder receipt defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none') AS x)

-- case: builder receipt schema tables
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('s.regions', 't.model', 'src', 'v1', 'asm', 'a', 'b', 'c') AS x)

-- case: builder receipt quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('re"g', 'mo"d', 'sr''c', 'v''1', 'as"m', 'a', 'b', 'c') AS x)

-- case: builder receipt null parameters
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('regions', 'model', NULL, NULL, NULL, NULL, NULL, NULL) AS x)

-- case: builder receipt regulation table
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', {'regulation_features_table': 'reg'}) AS x)

-- case: builder receipt regulation qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', {'regulation_features_table': 's.reg"x'}) AS x)

-- case: builder receipt regulation null
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', {'regulation_features_table': NULL::VARCHAR}) AS x)

-- case: builder receipt null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', NULL) AS x)

-- case: builder receipt null regions error
SELECT duckvep_model_receipt_sql(NULL, 'model', 'a', 'b', 'c', 'd', 'e', 'f') AS x

-- case: builder receipt null model error
SELECT duckvep_model_receipt_sql('regions', NULL, 'a', 'b', 'c', 'd', 'e', 'f') AS x

-- case: builder receipt unknown option error
SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', {'other': 'x'}) AS x

-- case: builder receipt option wrong type error
SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', {'regulation_features_table': 5}) AS x

-- case: builder receipt options not a struct error
SELECT duckvep_model_receipt_sql('regions', 'model', 'src', 'v1', 'GRCh38', 'mhash', 'rhash', 'none', 'x') AS x

-- case: builder receipt nul in parameter error
SELECT duckvep_model_receipt_sql('regions', 'model', 'a' || chr(0), 'b', 'c', 'd', 'e', 'f') AS x

-- case: builder receipt many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_model_receipt_sql('r' || (i % 3), 'm' || (i % 5), 's' || i, NULL, 'a', 'b', 'c', 'd', CASE WHEN i % 2 = 0 THEN {'regulation_features_table': 'rg' || (i % 7)} ELSE NULL END) AS x FROM range(1500) t(i))

-- fixture: temp tables for the builders whose SQL runs on both hosts
CREATE TEMP TABLE sv_records AS SELECT * FROM (VALUES (0::BIGINT, 13546123::BIGINT, 'G', '<INS>', 'SVTYPE=INS;END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG'), (1, 100, 'A', '<DEL>', 'END=220;CIPOS=-5,5'), (2, 500, 'G', 'GATG', '.'), (3, 9, 'A', '<DEL>', 'END=x')) t(event_index, pos, ref, alt, info)

-- fixture: breakend records
CREATE TEMP TABLE bnd_records AS SELECT * FROM (VALUES (0::INTEGER, '1', 100::BIGINT, 'a1', 'A', 'ATT[2:200[', 'MATEID=a2;EVENT=e1'), (1, '2', 200, 'a2', 'C', ']1:100]AAC', 'MATEID=a1;EVENT=e1'), (2, '3', 10, 's1', 'G', '.G', 'SVTYPE=BND'), (3, '4', 5, 'x1', 'T', 'T[9:9]', 'MATEID=zz'), (4, '5', 7, 'n1', 'A', 'C', '.')) t(event_index, chrom, pos, id, ref, alt, info)

-- fixture: breakend pairs
CREATE TEMP TABLE bnd_pairs AS SELECT * FROM query(duckvep_prepare_breakend_pairs_sql('bnd_records'))

-- fixture: breakend genes
CREATE TEMP TABLE bnd_genes AS SELECT * FROM (VALUES (0, 'G1'), (1, 'G2'), (2, 'G1')) t(event_index, gene_id)

-- fixture: structural HGVS events
CREATE TEMP TABLE shgvs_events AS SELECT * FROM (VALUES (0::INTEGER, '21', 7467462::DOUBLE, 'T', '<DEL>', 'END=7467465'), (1, '21', 7467462, 'T', 'N[21:7467900[', 'SVTYPE=BND'), (2, '21', 100, 'A', '<DUP>', 'END=110'), (3, '21', 200, 'C', '<INV>', 'END=203')) t(event_index, chrom, pos, ref, alt, info)

-- fixture: structural HGVS reference
CREATE TEMP TABLE shgvs_reference AS SELECT * FROM (VALUES (0::INTEGER, 'TCAG'), (2, 'AAAAAAAAAAA'), (3, 'CGTA')) t(event_index, reference_sequence)

-- fixture: ExpansionHunter records
CREATE TEMP TABLE eh_records AS SELECT * FROM (VALUES (0::INTEGER, 'END=2008;REF=1;RL=3;RU=CAG', 'GT:SO:REPCN:REPCI', '1/2:SPANNING/SPANNING:2/10:2-2/10-10', 'C', '<STR2>,<STR10>', 1::DOUBLE), (1, 'END=2008;REF=1;RL=3;RU=CAG', 'GT:SO:REPCN:REPCI', '1/2:SPANNING/SPANNING:2/10:2-2/10-10', 'C', '<STR2>,<STR10>', 2), (2, 'END=2008;REF=1;RL=3;RU=CAG', 'GT:SO:REPCN:REPCI', '1/1:SPANNING/SPANNING:2/2:2-2/2-2', 'C', '<STR2>', 1), (3, NULL, 'GT', '1', 'C', '<STR2>', 1)) t(event_index, info, format, "sample", ref, alt, alt_index)

-- fixture: ExpansionHunter reference
CREATE TEMP TABLE eh_reference AS SELECT event_index, 'CAG' AS reference_sequence FROM eh_records

-- case: builder executes prepare sv geometry
SELECT * FROM query(duckvep_prepare_sv_geometry_sql('sv_records')) ORDER BY event_index

-- case: builder executes prepare sv geometry qualified
SELECT * FROM query(duckvep_prepare_sv_geometry_sql('temp.sv_records')) ORDER BY event_index

-- case: builder executes breakend pairs
SELECT * FROM query(duckvep_prepare_breakend_pairs_sql('bnd_records')) ORDER BY event_index

-- case: builder executes breakend fusion
SELECT * FROM query(duckvep_prepare_breakend_fusion_sql('bnd_pairs', 'bnd_genes')) ORDER BY event_index

-- case: builder executes structural hgvs
SELECT * FROM query(duckvep_prepare_structural_hgvs_sql('shgvs_events', 'shgvs_reference')) ORDER BY event_index

-- case: builder executes structural hgvs with max_span
SELECT event_index, hgvs_status, hgvs_reason FROM query(duckvep_prepare_structural_hgvs_sql('shgvs_events', 'shgvs_reference', {'max_span': 5})) ORDER BY event_index

-- case: builder executes expansionhunter
SELECT * FROM query(duckvep_prepare_expansionhunter_sql('eh_records', 'eh_reference')) ORDER BY event_index

-- case: builder executes expansionhunter into repeat alleles
SELECT event_index, status, duckvep_repeat_alleles([reference_components], [alternate_components], sequence_exact) AS r FROM query(duckvep_prepare_expansionhunter_sql('eh_records', 'eh_reference')) ORDER BY event_index

-- fixture: model relations (README model)
CREATE TABLE readme_regions AS SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER)) t(seq_region)

-- fixture: model README transcripts
CREATE TABLE readme_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, 100::UBIGINT transcript_start, 250::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 120::UBIGINT cds_start, 240::UBIGINT cds_end, 'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence, 1::UTINYINT codon_table, 'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence, 'ACGTACGTAC'::BLOB post_cds_sequence

-- fixture: model README exons
CREATE TABLE readme_exons AS SELECT * FROM (VALUES (0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT), (0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)) t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase)

-- fixture: model rich regions
CREATE TABLE rich_regions AS SELECT * FROM (VALUES (0::UINTEGER, 100000::UBIGINT, 'chr0', false), (1::UINTEGER, 100000::UBIGINT, 'chr1', false)) t(seq_region, sequence_length, seq_region_name, circular)

-- fixture: model rich transcripts
CREATE TABLE rich_transcripts AS SELECT * FROM (VALUES (0::UINTEGER, 0::UINTEGER, 50::UBIGINT, 200::UBIGINT, 1::TINYINT, 0::UINTEGER, 3::UBIGINT, 61::UBIGINT, 190::UBIGINT, ('ATG' || repeat('GCA', 25) || 'TAA')::BLOB, 1::UTINYINT, 'ACGTACGTACG'::BLOB, 'TTGACCAGTA'::BLOB), (1::UINTEGER, 1::UINTEGER, 100::UBIGINT, 250::UBIGINT, 1::TINYINT, 1::UINTEGER, 3::UBIGINT, 120::UBIGINT, 240::UBIGINT, 'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB, 1::UTINYINT, 'TACGTACGTACGTACGTACG'::BLOB, 'ACGTACGTAC'::BLOB), (2::UINTEGER, 1::UINTEGER, 300::UBIGINT, 400::UBIGINT, -1::TINYINT, 2::UINTEGER, 3::UBIGINT, 310::UBIGINT, 390::UBIGINT, ('ATG' || repeat('CCA', 25) || 'TGA')::BLOB, 2::UTINYINT, 'GGGGGGGGGG'::BLOB, 'CCCCCCCCCC'::BLOB), (3::UINTEGER, 1::UINTEGER, 500::UBIGINT, 600::UBIGINT, 1::TINYINT, 3::UINTEGER, 8::UBIGINT, NULL::UBIGINT, NULL::UBIGINT, NULL::BLOB, NULL::UTINYINT, NULL::BLOB, NULL::BLOB)) t(transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence)

-- fixture: model rich exons
CREATE TABLE rich_exons AS SELECT * FROM (VALUES (0::UINTEGER, 50::UBIGINT, 100::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT), (0::UINTEGER, 150::UBIGINT, 200::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT), (1::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT), (1::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT), (2::UINTEGER, 300::UBIGINT, 400::UBIGINT, 1::UBIGINT, 101::UBIGINT, 0::TINYINT, 0::TINYINT), (3::UINTEGER, 500::UBIGINT, 600::UBIGINT, 1::UBIGINT, 101::UBIGINT, -1::TINYINT, -1::TINYINT)) t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase)

-- fixture: model rich mature miRNA
CREATE TABLE rich_mirna AS SELECT * FROM (VALUES (3::UINTEGER, 510::UBIGINT, 530::UBIGINT), (3::UINTEGER, 550::UBIGINT, 570::UBIGINT)) t(transcript_index, mature_mirna_start, mature_mirna_end)

-- fixture: model rich peptide edits
CREATE TABLE rich_peptides AS SELECT * FROM (VALUES (1::UINTEGER, 3::UINTEGER, 'W'), (1::UINTEGER, 9::UINTEGER, 'K')) t(transcript_index, protein_position, alternate_amino_acid)

-- fixture: model rich regulation
CREATE TABLE rich_regulation AS SELECT * FROM (VALUES (0::UINTEGER, 1::UINTEGER, 1000::UINTEGER, 1020::UINTEGER, 1::UTINYINT), (1::UINTEGER, 1::UINTEGER, 1000::UINTEGER, 1020::UINTEGER, 2::UTINYINT), (2::UINTEGER, 1::UINTEGER, 5000::UINTEGER, 5100::UINTEGER, 1::UTINYINT)) t(regulation_feature_index, seq_region, feature_start, feature_end, feature_kind)

-- fixture: model reference regions
CREATE TABLE ref_regions AS SELECT * FROM (VALUES (0::UINTEGER, 1::UBIGINT, 'chrEmpty'), (1::UINTEGER, 260::UBIGINT, 'chrDuck')) t(seq_region, sequence_length, seq_region_name)

-- fixture: model reference transcripts
CREATE TABLE ref_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, 11::UBIGINT transcript_start, 22::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 11::UBIGINT cds_start, 22::UBIGINT cds_end, 'ATGAAACCCGGG'::BLOB cds_sequence, 1::UTINYINT codon_table, ''::BLOB pre_cds_sequence, ''::BLOB post_cds_sequence

-- fixture: model reference exons
CREATE TABLE ref_exons AS SELECT 0::UINTEGER transcript_index, 11::UBIGINT exon_start, 22::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase

-- fixture-v1: load README model on the v1 host (private connection: permanent tables only)
SELECT loaded FROM duckvep_model_load('readme', 'SELECT * FROM readme_regions ORDER BY seq_region', 'SELECT * FROM readme_transcripts ORDER BY seq_region, transcript_start', 'SELECT * FROM readme_exons ORDER BY transcript_index, exon_start')

-- fixture-v1: load rich model on the v1 host
SELECT loaded FROM duckvep_model_load('rich', 'SELECT * FROM rich_regions ORDER BY seq_region', 'SELECT * FROM rich_transcripts ORDER BY seq_region, transcript_start', 'SELECT * FROM rich_exons ORDER BY transcript_index, exon_start', mature_mirna_query := 'SELECT * FROM rich_mirna ORDER BY transcript_index, mature_mirna_start', peptide_edit_query := 'SELECT * FROM rich_peptides ORDER BY transcript_index, protein_position', interval_feature_query := 'SELECT * FROM rich_regulation ORDER BY seq_region, feature_start, regulation_feature_index', transcript_coverage_complete := true)

-- fixture-v1: load reference model on the v1 host
SELECT loaded FROM duckvep_model_load('withref', 'SELECT * FROM ref_regions ORDER BY seq_region', 'SELECT * FROM ref_transcripts ORDER BY seq_region, transcript_start', 'SELECT * FROM ref_exons ORDER BY transcript_index, exon_start', reference_fasta := 'test/data/duckvep/minimal.fa')

-- fixture-v2: stage readme regions
COPY (
SELECT * FROM readme_regions ORDER BY seq_region
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'readme', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage readme transcripts
COPY (
SELECT * FROM readme_transcripts ORDER BY seq_region, transcript_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'readme', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage readme exons
COPY (
SELECT * FROM readme_exons ORDER BY transcript_index, exon_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'readme', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: publish readme
SELECT duckvep_model_publish('readme')

-- fixture-v2: stage rich regions
COPY (
SELECT * FROM rich_regions ORDER BY seq_region
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'rich', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage rich transcripts
COPY (
SELECT * FROM rich_transcripts ORDER BY seq_region, transcript_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'rich', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage rich exons
COPY (
SELECT * FROM rich_exons ORDER BY transcript_index, exon_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'rich', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage rich mature_mirna
COPY (
SELECT * FROM rich_mirna ORDER BY transcript_index, mature_mirna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'rich', RELATION 'mature_mirna', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage rich peptide_edits
COPY (
SELECT * FROM rich_peptides ORDER BY transcript_index, protein_position
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'rich', RELATION 'peptide_edits', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage rich interval_features
COPY (
SELECT * FROM rich_regulation ORDER BY seq_region, feature_start, regulation_feature_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'rich', RELATION 'interval_features', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: publish rich
SELECT duckvep_model_publish('rich', {'transcript_coverage_complete': true})

-- fixture-v2: stage withref regions
COPY (
SELECT * FROM ref_regions ORDER BY seq_region
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'withref', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage withref transcripts
COPY (
SELECT * FROM ref_transcripts ORDER BY seq_region, transcript_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'withref', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage withref exons
COPY (
SELECT * FROM ref_exons ORDER BY transcript_index, exon_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'withref', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: publish withref
SELECT duckvep_model_publish('withref', {'reference_fasta': 'test/data/duckvep/minimal.fa'})

-- case: model fingerprint readme
SELECT _duckvep_model_fingerprint('readme') AS f

-- case: model fingerprint rich
SELECT _duckvep_model_fingerprint('rich') AS f

-- case: model fingerprint with reference
SELECT _duckvep_model_fingerprint('withref') AS f

-- case: model fingerprints differ
SELECT count(DISTINCT f) AS n FROM (SELECT _duckvep_model_fingerprint(m) AS f FROM (VALUES ('readme'), ('rich'), ('withref')) v(m))

-- case: model fingerprint of an unknown model is null
SELECT _duckvep_model_fingerprint('no-such-model') IS NULL AS missing, _duckvep_model_fingerprint(NULL) IS NULL AS null_name

-- case: model drop of an unknown model
SELECT duckvep_model_drop('no-such-model') AS dropped

-- case: model drop empty name error
SELECT duckvep_model_drop('') AS dropped

-- fixture: hgvs model relations
CREATE TABLE hg_regions AS SELECT * FROM (VALUES (0::UINTEGER, 1::UBIGINT, 'chrEmpty'), (1::UINTEGER, 260::UBIGINT, 'chrDuck')) t(seq_region, sequence_length, seq_region_name)

-- fixture: hgvs model transcripts
CREATE TABLE hg_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, 21::UBIGINT transcript_start, 100::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 21::UBIGINT cds_start, 98::UBIGINT cds_end, (repeat('ACGT', 19) || 'AC')::BLOB cds_sequence, 1::UTINYINT codon_table, ''::BLOB pre_cds_sequence, 'GT'::BLOB post_cds_sequence

-- fixture: hgvs model exons
CREATE TABLE hg_exons AS SELECT 0::UINTEGER transcript_index, 21::UBIGINT exon_start, 100::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 80::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase

-- fixture: annotation events on the README model
CREATE TABLE ann_readme AS SELECT row_number() OVER (ORDER BY seq_region, position, alternate)::UBIGINT AS event_index, seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM (VALUES (0, 50, 'A', 'C'), (1, 90, 'A', 'C'), (1, 99, 'C', 'G'), (1, 100, 'A', 'T'), (1, 110, 'G', 'A'), (1, 121, 'T', 'C'), (1, 122, 'G', 'A'), (1, 123, 'A', 'GG'), (1, 124, 'T', 'G'), (1, 125, 'AC', 'A'), (1, 126, 'A', 'ACG'), (1, 130, 'AAC', 'A'), (1, 135, 'C', 'CAT'), (1, 149, 'G', 'T'), (1, 150, 'A', 'C'), (1, 151, 'A', 'G'), (1, 175, 'T', 'C'), (1, 199, 'A', 'G'), (1, 200, 'C', 'T'), (1, 239, 'A', 'T'), (1, 240, 'G', 'A'), (1, 241, 'C', 'G'), (1, 250, 'A', 'T'), (1, 255, 'A', 'C'), (1, 121, 'T', '<*>'), (1, 300, 'A', 'G')) v(seq_region, position, reference, alternate)

-- fixture: annotation events on the rich model
CREATE TABLE ann_rich AS SELECT row_number() OVER (ORDER BY seq_region, position, alternate)::UBIGINT AS event_index, seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM (VALUES (0, 49, 'A', 'T'), (0, 55, 'C', 'G'), (0, 61, 'A', 'C'), (0, 62, 'T', 'TA'), (0, 75, 'G', 'C'), (0, 99, 'A', 'G'), (0, 100, 'G', 'A'), (0, 101, 'A', 'C'), (0, 120, 'A', 'T'), (0, 149, 'C', 'T'), (0, 150, 'A', 'G'), (0, 151, 'T', 'C'), (0, 190, 'A', 'T'), (0, 191, 'A', 'G'), (0, 199, 'C', 'T'), (0, 205, 'A', 'C'), (1, 95, 'A', 'C'), (1, 110, 'G', 'A'), (1, 125, 'C', 'CA'), (1, 149, 'G', 'T'), (1, 160, 'A', 'T'), (1, 205, 'A', 'G'), (1, 239, 'A', 'T'), (1, 245, 'T', 'C'), (1, 290, 'A', 'G'), (1, 305, 'C', 'T'), (1, 330, 'CAA', 'C'), (1, 350, 'A', 'T'), (1, 391, 'A', 'G'), (1, 395, 'G', 'T'), (1, 405, 'A', 'C'), (1, 512, 'A', 'T'), (1, 520, 'G', 'C'), (1, 555, 'A', 'T'), (1, 580, 'C', 'G'), (1, 1005, 'A', 'T'), (1, 1010, 'C', 'CA'), (1, 5050, 'G', 'A'), (1, 9999, 'A', 'C')) v(seq_region, position, reference, alternate)

-- fixture: annotation events with a reference FASTA
CREATE TABLE ann_hgvs AS SELECT row_number() OVER (ORDER BY seq_region, position, alternate)::UBIGINT AS event_index, seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM (VALUES (1, 15, 'A', 'C'), (1, 21, 'A', 'G'), (1, 22, 'C', 'T'), (1, 23, 'G', 'A'), (1, 24, 'T', 'C'), (1, 25, 'A', 'T'), (1, 26, 'C', 'G'), (1, 27, 'G', 'T'), (1, 30, 'C', 'CA'), (1, 31, 'G', 'GTT'), (1, 34, 'ACG', 'A'), (1, 40, 'GTA', 'G'), (1, 50, 'A', 'ACGT'), (1, 60, 'T', 'TACG'), (1, 61, 'A', 'T'), (1, 80, 'C', 'A'), (1, 96, 'C', 'T'), (1, 97, 'G', 'C'), (1, 98, 'T', 'A'), (1, 99, 'A', 'T'), (1, 100, 'C', 'G'), (1, 101, 'G', 'C'), (1, 150, 'C', 'T')) v(seq_region, position, reference, alternate)

-- fixture: annotation events: mixed families
CREATE TABLE ann_mixed AS SELECT * FROM (VALUES (1::UBIGINT, 1::UINTEGER, 90::UBIGINT, NULL::VARCHAR, '<DUP>'::VARCHAR, 260::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, NULL::UINTEGER, NULL::UBIGINT), (2::UBIGINT, 1::UINTEGER, 100::UBIGINT, NULL::VARCHAR, '<DEL>'::VARCHAR, 180::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, NULL::UINTEGER, NULL::UBIGINT), (3::UBIGINT, 1::UINTEGER, 124::UBIGINT, 'T'::VARCHAR, 'C'::VARCHAR, NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, NULL::UINTEGER, NULL::UBIGINT), (4::UBIGINT, 1::UINTEGER, 130::UBIGINT, NULL::VARCHAR, '<INV>'::VARCHAR, 210::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, NULL::UINTEGER, NULL::UBIGINT), (5::UBIGINT, 1::UINTEGER, 159::UBIGINT, NULL::VARCHAR, 'N]1:170]'::VARCHAR, NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, 1::UINTEGER, 170::UBIGINT), (6::UBIGINT, 1::UINTEGER, 210::UBIGINT, NULL::VARCHAR, 'N[1:230['::VARCHAR, NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, 1::UINTEGER, 230::UBIGINT), (7::UBIGINT, 1::UINTEGER, 300::UBIGINT, 'A'::VARCHAR, 'G'::VARCHAR, NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR, NULL::UINTEGER, NULL::UBIGINT)) t(event_index, seq_region, position, reference, alternate, end_position, structural_type, copy_change, mate_seq_region, mate_position)

-- fixture: annotation events: literal only
CREATE TABLE ann_lit AS SELECT * FROM ann_readme WHERE alternate <> '<*>' ORDER BY event_index

-- fixture: annotation events: generated
CREATE TABLE ann_big AS SELECT (row_number() OVER (ORDER BY seq_region, position, alternate) - 1)::UBIGINT AS event_index, seq_region::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM (SELECT (i % 2)::INTEGER AS seq_region, (40 + (i * 7919) % 650 + (i % 2) * 100)::BIGINT AS position, substr('ACGT', 1 + i % 4, 1) AS reference, CASE i % 5 WHEN 0 THEN substr('ACGT', 1 + (i + 1) % 4, 1) WHEN 1 THEN substr('ACGT', 1 + (i + 2) % 4, 1) WHEN 2 THEN substr('ACGT', 1 + i % 4, 1) || 'GA' WHEN 3 THEN substr('ACGT', 1 + (i + 3) % 4, 1) ELSE 'T' END AS alternate FROM range(6000) t(i)) v WHERE alternate <> reference

-- fixture-v1: load hgvs model on the v1 host
SELECT loaded FROM duckvep_model_load('hgvsref', 'SELECT * FROM hg_regions ORDER BY seq_region', 'SELECT * FROM hg_transcripts ORDER BY seq_region, transcript_start', 'SELECT * FROM hg_exons ORDER BY transcript_index, exon_start', reference_fasta := 'test/data/duckvep/minimal.fa')

-- fixture-v2: stage hgvsref regions
COPY (
SELECT * FROM hg_regions ORDER BY seq_region
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'hgvsref', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage hgvsref transcripts
COPY (
SELECT * FROM hg_transcripts ORDER BY seq_region, transcript_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'hgvsref', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: stage hgvsref exons
COPY (
SELECT * FROM hg_exons ORDER BY transcript_index, exon_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'hgvsref', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: publish hgvsref
SELECT duckvep_model_publish('hgvsref', {'reference_fasta': 'test/data/duckvep/minimal.fa'})

-- fixture: projection transcripts of the README model
CREATE TABLE readme_projection_transcripts AS SELECT t.*, (SELECT list(struct_pack(exon_start := e.exon_start, exon_end := e.exon_end, exon_cdna_start := e.exon_cdna_start, exon_cdna_end := e.exon_cdna_end, phase := e.phase, end_phase := e.end_phase) ORDER BY e.exon_start) FROM readme_exons e WHERE e.transcript_index = t.transcript_index) AS exons, []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] AS peptide_edits FROM readme_transcripts t

-- fixture: projection annotations of the README model
CREATE TABLE readme_annotations AS SELECT * FROM query(duckvep_annotate_sql('ann_lit', 'readme'))

-- case: annotate readme compact
SELECT * FROM query(duckvep_annotate_sql('ann_readme', 'readme')) ORDER BY ALL

-- case: annotate readme compact hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_readme', 'readme')) t)

-- case: annotate readme rich
SELECT * FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'rich': true})) ORDER BY ALL

-- case: annotate readme rich hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'rich': true})) t)

-- case: annotate readme hgvs without reference
SELECT * FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'hgvs': true})) ORDER BY ALL

-- case: annotate readme hgvs without reference hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'hgvs': true})) t)

-- case: annotate readme rich hgvs
SELECT * FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'hgvs': true, 'rich': true})) ORDER BY ALL

-- case: annotate readme rich hgvs hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'hgvs': true, 'rich': true})) t)

-- case: annotate readme distances zero
SELECT * FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'upstream_distance': 0, 'downstream_distance': 0, 'rich': true})) ORDER BY ALL

-- case: annotate readme distances zero hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'upstream_distance': 0, 'downstream_distance': 0, 'rich': true})) t)

-- case: annotate readme distances small
SELECT * FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'upstream_distance': 10, 'downstream_distance': 20})) ORDER BY ALL

-- case: annotate readme distances small hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_readme', 'readme', {'upstream_distance': 10, 'downstream_distance': 20})) t)

-- case: annotate readme mixed families compact
SELECT * FROM query(duckvep_annotate_sql('ann_mixed', 'readme')) ORDER BY ALL

-- case: annotate readme mixed families compact hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_mixed', 'readme')) t)

-- case: annotate readme mixed families rich
SELECT * FROM query(duckvep_annotate_sql('ann_mixed', 'readme', {'rich': true})) ORDER BY ALL

-- case: annotate readme mixed families rich hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_mixed', 'readme', {'rich': true})) t)

-- case: annotate rich model compact
SELECT * FROM query(duckvep_annotate_sql('ann_rich', 'rich')) ORDER BY ALL

-- case: annotate rich model compact hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_rich', 'rich')) t)

-- case: annotate rich model rich
SELECT * FROM query(duckvep_annotate_sql('ann_rich', 'rich', {'rich': true, 'upstream_distance': 100})) ORDER BY ALL

-- case: annotate rich model rich hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_rich', 'rich', {'rich': true, 'upstream_distance': 100})) t)

-- case: annotate rich model hgvs
SELECT * FROM query(duckvep_annotate_sql('ann_rich', 'rich', {'hgvs': true, 'rich': true})) ORDER BY ALL

-- case: annotate rich model hgvs hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_rich', 'rich', {'hgvs': true, 'rich': true})) t)

-- case: annotate rich model distances
SELECT * FROM query(duckvep_annotate_sql('ann_rich', 'rich', {'upstream_distance': 0, 'downstream_distance': 5000})) ORDER BY ALL

-- case: annotate rich model distances hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_rich', 'rich', {'upstream_distance': 0, 'downstream_distance': 5000})) t)

-- case: annotate reference model compact
SELECT * FROM query(duckvep_annotate_sql('ann_hgvs', 'hgvsref')) ORDER BY ALL

-- case: annotate reference model compact hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_hgvs', 'hgvsref')) t)

-- case: annotate reference model hgvs
SELECT * FROM query(duckvep_annotate_sql('ann_hgvs', 'hgvsref', {'hgvs': true})) ORDER BY ALL

-- case: annotate reference model hgvs hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_hgvs', 'hgvsref', {'hgvs': true})) t)

-- case: annotate reference model rich hgvs
SELECT * FROM query(duckvep_annotate_sql('ann_hgvs', 'hgvsref', {'hgvs': true, 'rich': true})) ORDER BY ALL

-- case: annotate reference model rich hgvs hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_hgvs', 'hgvsref', {'hgvs': true, 'rich': true})) t)

-- case: annotate readme projected
SELECT * FROM query(duckvep_annotate_projected_sql('ann_lit', 'readme')) ORDER BY ALL

-- case: annotate readme projected hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_projected_sql('ann_lit', 'readme')) t)

-- case: annotate readme projected distances
SELECT * FROM query(duckvep_annotate_projected_sql('ann_lit', 'readme', {'upstream_distance': 0, 'downstream_distance': 0})) ORDER BY ALL

-- case: annotate readme projected distances hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_projected_sql('ann_lit', 'readme', {'upstream_distance': 0, 'downstream_distance': 0})) t)

-- case: annotate rich model projected
SELECT * FROM query(duckvep_annotate_projected_sql('ann_rich', 'rich')) ORDER BY ALL

-- case: annotate rich model projected hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_projected_sql('ann_rich', 'rich')) t)

-- case: annotate reference model projected
SELECT * FROM query(duckvep_annotate_projected_sql('ann_hgvs', 'hgvsref')) ORDER BY ALL

-- case: annotate reference model projected hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_projected_sql('ann_hgvs', 'hgvsref')) t)

-- case: annotate transcript projection
SELECT * FROM query(duckvep_transcript_projection_sql('ann_lit', 'readme_annotations', 'readme_projection_transcripts')) ORDER BY ALL

-- case: annotate transcript projection hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_transcript_projection_sql('ann_lit', 'readme_annotations', 'readme_projection_transcripts')) t)

-- case: annotate generated events compact
SELECT * FROM query(duckvep_annotate_sql('ann_big', 'rich')) ORDER BY ALL

-- case: annotate generated events compact hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_big', 'rich')) t)

-- case: annotate generated events rich
SELECT * FROM query(duckvep_annotate_sql('ann_big', 'rich', {'rich': true})) ORDER BY ALL

-- case: annotate generated events rich hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_big', 'rich', {'rich': true})) t)

-- case: annotate generated events hgvs
SELECT * FROM query(duckvep_annotate_sql('ann_big', 'rich', {'hgvs': true, 'rich': true, 'upstream_distance': 200})) ORDER BY ALL

-- case: annotate generated events hgvs hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_sql('ann_big', 'rich', {'hgvs': true, 'rich': true, 'upstream_distance': 200})) t)

-- case: annotate generated events projected
SELECT * FROM query(duckvep_annotate_projected_sql('ann_big', 'rich')) ORDER BY ALL

-- case: annotate generated events projected hashes
SELECT count(*) AS n, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h FROM query(duckvep_annotate_projected_sql('ann_big', 'rich')) t)

-- case: annotate unknown model error
SELECT count(*) FROM query(duckvep_annotate_sql('ann_readme', 'no-such-model'))

-- case: annotate natives types
SELECT typeof(_duckvep_annotate_small_rich('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C')) AS a, typeof(_duckvep_annotate_small_compact('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C')) AS b, typeof(_duckvep_annotate_small_hgvs('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C')) AS c, typeof(_duckvep_annotate_small_rich_hgvs('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C')) AS d, typeof(_duckvep_annotate_small_projected('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C')) AS e, typeof(_duckvep_annotate_small_projected_hgvs('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C')) AS f, typeof(_duckvep_annotate_structural_rich('readme', 1::UINTEGER, 90::UBIGINT, 200::UBIGINT, 'DEL', 'LOSS')) AS g, typeof(_duckvep_annotate_structural_compact('readme', 1::UINTEGER, 90::UBIGINT, 200::UBIGINT, 'DEL', 'LOSS')) AS h, typeof(_duckvep_annotate_breakend_rich('readme', 1::UINTEGER, 159::UBIGINT, 1::UINTEGER, 170::UBIGINT)) AS i, typeof(_duckvep_annotate_breakend_compact('readme', 1::UINTEGER, 159::UBIGINT, 1::UINTEGER, 170::UBIGINT)) AS j

-- case: annotate natives direct
SELECT u.* FROM (SELECT unnest(_duckvep_annotate_small_rich('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'C', 0::UBIGINT, 0::UBIGINT)) AS u) ORDER BY ALL

-- case: annotate natives distance overloads
SELECT len(_duckvep_annotate_small_compact('readme', 1::UINTEGER, 90::UBIGINT, 'A', 'C')) AS a, len(_duckvep_annotate_small_compact('readme', 1::UINTEGER, 90::UBIGINT, 'A', 'C', 10::UBIGINT)) AS b, len(_duckvep_annotate_small_compact('readme', 1::UINTEGER, 90::UBIGINT, 'A', 'C', 10::UBIGINT, 10::UBIGINT)) AS c

-- case: annotate natives null allele gives null
SELECT _duckvep_annotate_small_compact('readme', 1::UINTEGER, 124::UBIGINT, NULL::VARCHAR, 'C') AS r

-- case: annotate natives bad allele error
SELECT _duckvep_annotate_small_compact('readme', 1::UINTEGER, 124::UBIGINT, 'T', 'X') AS r

-- case: annotate natives empty model error
SELECT _duckvep_annotate_small_compact('', 1::UINTEGER, 124::UBIGINT, 'T', 'C') AS r

-- case: annotate natives long alleles
SELECT len(_duckvep_annotate_small_rich('readme', 1::UINTEGER, 90::UBIGINT, repeat('A', 300), repeat('A', 100) || 'C' || repeat('A', 199))) AS n, md5(to_json(_duckvep_annotate_small_rich('readme', 1::UINTEGER, 90::UBIGINT, repeat('A', 300), repeat('A', 100) || 'C' || repeat('A', 199)))::VARCHAR) AS h

-- case: annotate natives many rows
SELECT count(*) AS n, sum(hash(u)::HUGEINT) AS s, bit_xor(hash(u)) AS x FROM (SELECT unnest(_duckvep_annotate_small_rich('readme', 1::UINTEGER, (90 + i % 200)::UBIGINT, substr('ACGT', 1 + i % 4, 1), substr('ACGT', 1 + (i + 1) % 4, 1))) AS u FROM range(9000) t(i))

-- case: projection code
SELECT string_agg(__duckvep_projection_code(c::UTINYINT), '|' ORDER BY c) AS codes FROM (VALUES (1), (2), (3)) v(c)

-- case: projection code null
SELECT __duckvep_projection_code(NULL::UTINYINT) IS NULL AS n

-- case: projection code unsupported error
SELECT __duckvep_projection_code(200::UTINYINT) AS c

-- fixture: lof README annotations
CREATE TABLE lof_ann_readme AS SELECT * FROM query(duckvep_annotate_projected_sql('ann_lit', 'readme'))

-- fixture: lof README transcripts
CREATE TABLE lof_tx_readme AS SELECT t.transcript_index, '1' AS seq_region_name, t.strand, t.cds_start, t.cds_end, 'protein_coding' AS transcript_biotype, 'ENST_README' AS transcript_stable_id, (SELECT list(struct_pack(exon_start := e.exon_start, exon_end := e.exon_end) ORDER BY e.exon_start) FROM readme_exons e WHERE e.transcript_index = t.transcript_index) AS exons FROM readme_transcripts t

-- fixture: lof README reference
CREATE TABLE lof_ref_readme AS SELECT '1' AS chrom, 0 AS "start", 300 AS "end", repeat('A', 150) || 'GT' || repeat('T', 45) || 'AG' || repeat('A', 101) AS seq

-- fixture: lof gerp
CREATE TABLE lof_gerp AS SELECT '1' AS chrom, i * 10 AS "start", i * 10 + 10 AS "end", (i % 7) - 3.5 AS score FROM range(30) t(i)

-- fixture: lof ancestor
CREATE TABLE lof_anc AS SELECT '1' AS chrom, 0 AS "start", 300 AS "end", repeat('ACGT', 75) AS seq

-- fixture: lof phylocsf
CREATE TABLE lof_pcsf AS SELECT * FROM (VALUES ('ENST_README', 1, -1.5, 2.0), ('ENST_README', 2, 3.0, 4.0)) t(transcript, exon, corresponding_orf_score, max_score)

-- fixture: lof rich annotations
CREATE TABLE lof_ann_rich AS SELECT * FROM query(duckvep_annotate_projected_sql('ann_rich', 'rich'))

-- fixture: lof rich transcripts
CREATE TABLE lof_tx_rich AS SELECT t.transcript_index, 'chr' || t.seq_region AS seq_region_name, t.strand, t.cds_start, t.cds_end, CASE WHEN t.transcript_index = 3 THEN 'miRNA' ELSE 'protein_coding' END AS transcript_biotype, 'ENST_RICH' || t.transcript_index AS transcript_stable_id, (SELECT list(struct_pack(exon_start := e.exon_start, exon_end := e.exon_end) ORDER BY e.exon_start) FROM rich_exons e WHERE e.transcript_index = t.transcript_index) AS exons FROM rich_transcripts t

-- fixture: lof rich reference
CREATE TABLE lof_ref_rich AS SELECT * FROM (SELECT 'chr0' AS chrom, 0 AS "start", 300 AS "end", repeat('ACGTTGCA', 37) || 'GTAG' AS seq UNION ALL SELECT 'chr1', 0, 700, repeat('GATTACA', 100))

-- case: builder lof defaults
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a', 't', 'r') AS x)

-- case: builder lof schema qualified
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('s.a', 's.t', 's.r') AS x)

-- case: builder lof quoting
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a"1', 't"2', 'r"3') AS x)

-- case: builder lof gerp
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp': 'g'}) AS x)

-- case: builder lof ancestor
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'ancestor': 's.anc'}) AS x)

-- case: builder lof phylocsf
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'phylocsf': 'p'}) AS x)

-- case: builder lof all relations
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp': 'g', 'ancestor': 'n', 'phylocsf': 'p'}) AS x)

-- case: builder lof min intron
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': 30}) AS x)

-- case: builder lof min intron types
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': 7::UTINYINT}) AS x)

-- case: builder lof min intron bounds
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': 1000000000}) AS x)

-- case: builder lof cutoff double
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': -12.25::DOUBLE}) AS x)

-- case: builder lof cutoff integer
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': 4}) AS x)

-- case: builder lof cutoff float
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': 1.5::FLOAT}) AS x)

-- case: builder lof cutoff decimal error
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': 1.5::DECIMAL(4,1)}) AS x)

-- case: builder lof check cds
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'check_complete_cds': true}) AS x)

-- case: builder lof all options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp': 'g', 'ancestor': 'n', 'phylocsf': 'p', 'min_intron_size': 20, 'gerp_end_trunc_cutoff': -40.5::DOUBLE, 'check_complete_cds': true}) AS x)

-- case: builder lof null options
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',NULL) AS x)

-- case: builder lof null option values
SELECT md5(x) AS h, length(x) AS n FROM (SELECT duckvep_lof_sql('a','t','r',{'gerp': NULL::VARCHAR, 'min_intron_size': NULL::INTEGER, 'gerp_end_trunc_cutoff': NULL, 'check_complete_cds': NULL::BOOLEAN}) AS x)

-- case: builder lof unknown option error
SELECT duckvep_lof_sql('a','t','r',{'other': 1}) AS x

-- case: builder lof options not a struct error
SELECT duckvep_lof_sql('a','t','r',5) AS x

-- case: builder lof min intron wrong type error
SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': 'x'}) AS x

-- case: builder lof min intron too large error
SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': 1000000001}) AS x

-- case: builder lof min intron negative error
SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': -1}) AS x

-- case: builder lof min intron huge unsigned error
SELECT duckvep_lof_sql('a','t','r',{'min_intron_size': 18446744073709551615::UBIGINT}) AS x

-- case: builder lof cutoff wrong type error
SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': 'x'}) AS x

-- case: builder lof cutoff not finite error
SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': 'Infinity'::DOUBLE}) AS x

-- case: builder lof cutoff nan error
SELECT duckvep_lof_sql('a','t','r',{'gerp_end_trunc_cutoff': 'NaN'::DOUBLE}) AS x

-- case: builder lof check cds wrong type error
SELECT duckvep_lof_sql('a','t','r',{'check_complete_cds': 1}) AS x

-- case: builder lof gerp wrong type error
SELECT duckvep_lof_sql('a','t','r',{'gerp': 5}) AS x

-- case: builder lof null annotations error
SELECT duckvep_lof_sql(NULL,'t','r') AS x

-- case: builder lof null reference error
SELECT duckvep_lof_sql('a','t',NULL) AS x

-- case: builder lof extra dots error
SELECT duckvep_lof_sql('a.b.c','t','r') AS x

-- case: builder lof empty transcripts error
SELECT duckvep_lof_sql('a','','r') AS x

-- case: builder lof many rows
SELECT count(*) AS n, md5(string_agg(x, '' ORDER BY i)) AS h FROM (SELECT i, duckvep_lof_sql('a' || (i % 5), 't' || (i % 3), 'r', {'min_intron_size': i % 50, 'gerp': CASE WHEN i % 2 = 0 THEN 'g' ELSE NULL END}) AS x FROM range(1200) t(i))

-- case: lof executes README defaults
SELECT * FROM query(duckvep_lof_sql('lof_ann_readme', 'lof_tx_readme', 'lof_ref_readme')) ORDER BY ALL

-- case: lof executes README defaults hashes
SELECT count(*) AS n, count(lof) AS calls, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h, t.lof FROM query(duckvep_lof_sql('lof_ann_readme', 'lof_tx_readme', 'lof_ref_readme')) t)

-- case: lof executes README all options
SELECT * FROM query(duckvep_lof_sql('lof_ann_readme', 'lof_tx_readme', 'lof_ref_readme', {'gerp': 'lof_gerp', 'ancestor': 'lof_anc', 'phylocsf': 'lof_pcsf', 'min_intron_size': 40, 'check_complete_cds': true})) ORDER BY ALL

-- case: lof executes README all options hashes
SELECT count(*) AS n, count(lof) AS calls, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h, t.lof FROM query(duckvep_lof_sql('lof_ann_readme', 'lof_tx_readme', 'lof_ref_readme', {'gerp': 'lof_gerp', 'ancestor': 'lof_anc', 'phylocsf': 'lof_pcsf', 'min_intron_size': 40, 'check_complete_cds': true})) t)

-- case: lof executes README gerp cutoff
SELECT * FROM query(duckvep_lof_sql('lof_ann_readme', 'lof_tx_readme', 'lof_ref_readme', {'gerp': 'lof_gerp', 'gerp_end_trunc_cutoff': -2.5::DOUBLE})) ORDER BY ALL

-- case: lof executes README gerp cutoff hashes
SELECT count(*) AS n, count(lof) AS calls, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h, t.lof FROM query(duckvep_lof_sql('lof_ann_readme', 'lof_tx_readme', 'lof_ref_readme', {'gerp': 'lof_gerp', 'gerp_end_trunc_cutoff': -2.5::DOUBLE})) t)

-- case: lof executes rich model defaults
SELECT * FROM query(duckvep_lof_sql('lof_ann_rich', 'lof_tx_rich', 'lof_ref_rich')) ORDER BY ALL

-- case: lof executes rich model defaults hashes
SELECT count(*) AS n, count(lof) AS calls, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h, t.lof FROM query(duckvep_lof_sql('lof_ann_rich', 'lof_tx_rich', 'lof_ref_rich')) t)

-- case: lof executes rich model check cds
SELECT * FROM query(duckvep_lof_sql('lof_ann_rich', 'lof_tx_rich', 'lof_ref_rich', {'check_complete_cds': true, 'min_intron_size': 100})) ORDER BY ALL

-- case: lof executes rich model check cds hashes
SELECT count(*) AS n, count(lof) AS calls, sum(h::HUGEINT) AS s, bit_xor(h) AS x FROM (SELECT hash(t) AS h, t.lof FROM query(duckvep_lof_sql('lof_ann_rich', 'lof_tx_rich', 'lof_ref_rich', {'check_complete_cds': true, 'min_intron_size': 100})) t)

-- Haplotype fixtures: the model, transcripts, exons and phased calls of the haplotype suites (test/sql/duckvep_haplotype_*.test),
-- read from test/data/haplotype. "-- job:" blocks are haplotype inputs: v1 reads the query itself (duckvep_haplotypes), v2 stages it with
-- the statements of duckvep_haplotype_load_sql and scans it (duckvep_haplotype_scan); both are reached through the table macro of the job name.

-- fixture: haplotype vt tx
CREATE TABLE vt_tx AS SELECT seq_region::UINTEGER transcript_index,seq_region::UINTEGER seq_region, transcript_start::UBIGINT transcript_start,transcript_end::UBIGINT transcript_end,strand::TINYINT strand, seq_region::UINTEGER gene_index,3::UBIGINT transcript_flags,transcript_start::UBIGINT cds_start, transcript_end::UBIGINT cds_end,cds_sequence::BLOB cds_sequence,1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence,"case" AS case_name FROM read_csv('test/data/haplotype/vertical_transcripts.tsv',delim='\t',header=true)

-- fixture: haplotype vt exons
CREATE TABLE vt_exons AS SELECT seq_region::UINTEGER transcript_index,exon_start::UBIGINT exon_start, exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,exon_cdna_end::UBIGINT exon_cdna_end, phase::TINYINT phase,end_phase::TINYINT end_phase FROM read_csv('test/data/haplotype/vertical_exons.tsv',delim='\t',header=true)

-- fixture: haplotype vt vcf
CREATE TABLE vt_vcf AS SELECT string_split(line,chr(9)) f FROM (SELECT unnest(string_split(content,chr(10))) line FROM read_text('test/data/haplotype/vertical.vcf')) WHERE line<>'' AND NOT starts_with(line,'#')

-- fixture: haplotype vt calls
CREATE TABLE vt_calls AS SELECT (row_number() OVER ())::BIGINT event_index,t.seq_region::INT seq_region, v.f[2]::BIGINT AS position,v.f[4] AS reference,v.f[5] AS alternate,1 AS alt_index,t.transcript_index::INT transcript_index, 0 sample_index,(CASE v.f[10] WHEN '1|0' THEN [1,0] WHEN '0|1' THEN [0,1] ELSE [1,1] END)::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set FROM vt_vcf v JOIN vt_tx t ON t.case_name=v.f[1]

-- fixture-v1: haplotype vt model
SELECT loaded FROM duckvep_model_load('vt', 'SELECT i::UINTEGER seq_region FROM range(13) t(i)', 'SELECT * EXCLUDE(case_name) FROM vt_tx ORDER BY transcript_index', 'SELECT * FROM vt_exons ORDER BY transcript_index,exon_cdna_start')

-- fixture-v2: haplotype vt stage regions
COPY (
SELECT i::UINTEGER seq_region FROM range(13) t(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'vt', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype vt stage transcripts
COPY (
SELECT * EXCLUDE(case_name) FROM vt_tx ORDER BY transcript_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'vt', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype vt stage exons
COPY (
SELECT * FROM vt_exons ORDER BY transcript_index,exon_cdna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'vt', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype vt publish
SELECT duckvep_model_publish('vt')

-- job: hap_vt vt
SELECT * FROM vt_calls

-- fixture: haplotype sc tx
CREATE TABLE sc_tx AS SELECT seq_region::UINTEGER transcript_index,seq_region::UINTEGER seq_region, transcript_start::UBIGINT transcript_start,transcript_end::UBIGINT transcript_end,strand::TINYINT strand, seq_region::UINTEGER gene_index,3::UBIGINT transcript_flags,transcript_start::UBIGINT cds_start, transcript_end::UBIGINT cds_end,cds_sequence::BLOB cds_sequence,1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence,"case" AS case_name FROM read_csv('test/data/haplotype/same_codon_transcripts.tsv',delim='\t',header=true)

-- fixture: haplotype sc exons
CREATE TABLE sc_exons AS SELECT seq_region::UINTEGER transcript_index,exon_start::UBIGINT exon_start, exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,exon_cdna_end::UBIGINT exon_cdna_end, phase::TINYINT phase,end_phase::TINYINT end_phase FROM read_csv('test/data/haplotype/same_codon_exons.tsv',delim='\t',header=true)

-- fixture: haplotype sc vcf
CREATE TABLE sc_vcf AS SELECT string_split(line,chr(9)) f FROM (SELECT unnest(string_split(content,chr(10))) line FROM read_text('test/data/haplotype/same_codon.vcf')) WHERE line<>'' AND NOT starts_with(line,'#')

-- fixture: haplotype sc calls
CREATE TABLE sc_calls AS SELECT (row_number() OVER ())::BIGINT event_index,t.seq_region::INT seq_region, v.f[2]::BIGINT AS position,v.f[4] AS reference,v.f[5] AS alternate,1 AS alt_index,t.transcript_index::INT transcript_index, 0 sample_index,(CASE v.f[10] WHEN '1|0' THEN [1,0] WHEN '0|1' THEN [0,1] ELSE [1,1] END)::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set FROM sc_vcf v JOIN sc_tx t ON t.case_name=v.f[1]

-- fixture-v1: haplotype sc model
SELECT loaded FROM duckvep_model_load('sc', 'SELECT i::UINTEGER seq_region FROM range(138) t(i)', 'SELECT * EXCLUDE(case_name) FROM sc_tx ORDER BY transcript_index', 'SELECT * FROM sc_exons ORDER BY transcript_index,exon_cdna_start')

-- fixture-v2: haplotype sc stage regions
COPY (
SELECT i::UINTEGER seq_region FROM range(138) t(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'sc', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype sc stage transcripts
COPY (
SELECT * EXCLUDE(case_name) FROM sc_tx ORDER BY transcript_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'sc', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype sc stage exons
COPY (
SELECT * FROM sc_exons ORDER BY transcript_index,exon_cdna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'sc', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype sc publish
SELECT duckvep_model_publish('sc')

-- job: hap_sc sc
SELECT * FROM sc_calls

-- fixture: haplotype fr tx
CREATE TABLE fr_tx AS SELECT seq_region::UINTEGER transcript_index,seq_region::UINTEGER seq_region, transcript_start::UBIGINT transcript_start,transcript_end::UBIGINT transcript_end,strand::TINYINT strand, seq_region::UINTEGER gene_index,3::UBIGINT transcript_flags,transcript_start::UBIGINT cds_start, transcript_end::UBIGINT cds_end,cds_sequence::BLOB cds_sequence,1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence,"case" AS case_name FROM read_csv('test/data/haplotype/frame_transcripts.tsv',delim='\t',header=true)

-- fixture: haplotype fr exons
CREATE TABLE fr_exons AS SELECT seq_region::UINTEGER transcript_index,exon_start::UBIGINT exon_start, exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,exon_cdna_end::UBIGINT exon_cdna_end, phase::TINYINT phase,end_phase::TINYINT end_phase FROM read_csv('test/data/haplotype/frame_exons.tsv',delim='\t',header=true)

-- fixture: haplotype fr vcf
CREATE TABLE fr_vcf AS SELECT string_split(line,chr(9)) f FROM (SELECT unnest(string_split(content,chr(10))) line FROM read_text('test/data/haplotype/frame.vcf')) WHERE line<>'' AND NOT starts_with(line,'#')

-- fixture: haplotype fr calls
CREATE TABLE fr_calls AS SELECT (row_number() OVER ())::BIGINT event_index,t.seq_region::INT seq_region, v.f[2]::BIGINT AS position,v.f[4] AS reference,v.f[5] AS alternate,1 AS alt_index,t.transcript_index::INT transcript_index, 0 sample_index,(CASE v.f[10] WHEN '1|0' THEN [1,0] WHEN '0|1' THEN [0,1] ELSE [1,1] END)::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set FROM fr_vcf v JOIN fr_tx t ON t.case_name=v.f[1]

-- fixture-v1: haplotype fr model
SELECT loaded FROM duckvep_model_load('fr', 'SELECT i::UINTEGER seq_region FROM range(180) t(i)', 'SELECT * EXCLUDE(case_name) FROM fr_tx ORDER BY transcript_index', 'SELECT * FROM fr_exons ORDER BY transcript_index,exon_cdna_start')

-- fixture-v2: haplotype fr stage regions
COPY (
SELECT i::UINTEGER seq_region FROM range(180) t(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'fr', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype fr stage transcripts
COPY (
SELECT * EXCLUDE(case_name) FROM fr_tx ORDER BY transcript_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'fr', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype fr stage exons
COPY (
SELECT * FROM fr_exons ORDER BY transcript_index,exon_cdna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'fr', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype fr publish
SELECT duckvep_model_publish('fr')

-- job: hap_fr fr
SELECT * FROM fr_calls

-- fixture: haplotype ss tx
CREATE TABLE ss_tx AS SELECT seq_region::UINTEGER transcript_index,seq_region::UINTEGER seq_region, transcript_start::UBIGINT transcript_start,transcript_end::UBIGINT transcript_end,strand::TINYINT strand, seq_region::UINTEGER gene_index,3::UBIGINT transcript_flags,transcript_start::UBIGINT cds_start, transcript_end::UBIGINT cds_end,cds_sequence::BLOB cds_sequence,1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence,"case" AS case_name FROM read_csv('test/data/haplotype/startstop_transcripts.tsv',delim='\t',header=true)

-- fixture: haplotype ss exons
CREATE TABLE ss_exons AS SELECT seq_region::UINTEGER transcript_index,exon_start::UBIGINT exon_start, exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,exon_cdna_end::UBIGINT exon_cdna_end, phase::TINYINT phase,end_phase::TINYINT end_phase FROM read_csv('test/data/haplotype/startstop_exons.tsv',delim='\t',header=true)

-- fixture: haplotype ss vcf
CREATE TABLE ss_vcf AS SELECT string_split(line,chr(9)) f FROM (SELECT unnest(string_split(content,chr(10))) line FROM read_text('test/data/haplotype/startstop.vcf')) WHERE line<>'' AND NOT starts_with(line,'#')

-- fixture: haplotype ss calls
CREATE TABLE ss_calls AS SELECT (row_number() OVER ())::BIGINT event_index,t.seq_region::INT seq_region, v.f[2]::BIGINT AS position,v.f[4] AS reference,v.f[5] AS alternate,1 AS alt_index,t.transcript_index::INT transcript_index, 0 sample_index,(CASE v.f[10] WHEN '1|0' THEN [1,0] WHEN '0|1' THEN [0,1] ELSE [1,1] END)::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set FROM ss_vcf v JOIN ss_tx t ON t.case_name=v.f[1]

-- fixture-v1: haplotype ss model
SELECT loaded FROM duckvep_model_load('ss', 'SELECT i::UINTEGER seq_region FROM range(348) t(i)', 'SELECT * EXCLUDE(case_name) FROM ss_tx ORDER BY transcript_index', 'SELECT * FROM ss_exons ORDER BY transcript_index,exon_cdna_start')

-- fixture-v2: haplotype ss stage regions
COPY (
SELECT i::UINTEGER seq_region FROM range(348) t(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'ss', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype ss stage transcripts
COPY (
SELECT * EXCLUDE(case_name) FROM ss_tx ORDER BY transcript_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'ss', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype ss stage exons
COPY (
SELECT * FROM ss_exons ORDER BY transcript_index,exon_cdna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'ss', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype ss publish
SELECT duckvep_model_publish('ss')

-- job: hap_ss ss
SELECT * FROM ss_calls

-- fixture: haplotype nm tx
CREATE TABLE nm_tx AS SELECT seq_region::UINTEGER transcript_index,seq_region::UINTEGER seq_region, transcript_start::UBIGINT transcript_start,transcript_end::UBIGINT transcript_end,strand::TINYINT strand, seq_region::UINTEGER gene_index,3::UBIGINT transcript_flags,cds_start::UBIGINT cds_start, cds_end::UBIGINT cds_end,cds_sequence::BLOB cds_sequence,1::UTINYINT codon_table, coalesce(nullif(pre_cds_sequence,'.'),'')::BLOB pre_cds_sequence,coalesce(nullif(post_cds_sequence,'.'),'')::BLOB post_cds_sequence,"case" AS case_name FROM read_csv('test/data/haplotype/nmd_transcripts.tsv',delim='\t',header=true)

-- fixture: haplotype nm exons
CREATE TABLE nm_exons AS SELECT seq_region::UINTEGER transcript_index,exon_start::UBIGINT exon_start, exon_end::UBIGINT exon_end,exon_cdna_start::UBIGINT exon_cdna_start,exon_cdna_end::UBIGINT exon_cdna_end, phase::TINYINT phase,end_phase::TINYINT end_phase FROM read_csv('test/data/haplotype/nmd_exons.tsv',delim='\t',header=true)

-- fixture: haplotype nm vcf
CREATE TABLE nm_vcf AS SELECT string_split(line,chr(9)) f FROM (SELECT unnest(string_split(content,chr(10))) line FROM read_text('test/data/haplotype/nmd.vcf')) WHERE line<>'' AND NOT starts_with(line,'#')

-- fixture: haplotype nm calls
CREATE TABLE nm_calls AS SELECT (row_number() OVER ())::BIGINT event_index,t.seq_region::INT seq_region, v.f[2]::BIGINT AS position,v.f[4] AS reference,v.f[5] AS alternate,1 AS alt_index,t.transcript_index::INT transcript_index, 0 sample_index,(CASE v.f[10] WHEN '1|0' THEN [1,0] WHEN '0|1' THEN [0,1] ELSE [1,1] END)::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set FROM nm_vcf v JOIN nm_tx t ON t.case_name=v.f[1]

-- fixture-v1: haplotype nm model
SELECT loaded FROM duckvep_model_load('nm', 'SELECT i::UINTEGER seq_region FROM range(238) t(i)', 'SELECT * EXCLUDE(case_name) FROM nm_tx ORDER BY transcript_index', 'SELECT * FROM nm_exons ORDER BY transcript_index,exon_cdna_start')

-- fixture-v2: haplotype nm stage regions
COPY (
SELECT i::UINTEGER seq_region FROM range(238) t(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'nm', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype nm stage transcripts
COPY (
SELECT * EXCLUDE(case_name) FROM nm_tx ORDER BY transcript_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'nm', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype nm stage exons
COPY (
SELECT * FROM nm_exons ORDER BY transcript_index,exon_cdna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'nm', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: haplotype nm publish
SELECT duckvep_model_publish('nm')

-- job: hap_nm nm
SELECT * FROM nm_calls


-- Haplotype cases. Each table macro is scanned once per case.
-- case: haplotypes vt full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x, count(DISTINCT transcript_index) AS transcripts FROM (SELECT * FROM hap_vt()) h

-- case: haplotypes vt first rows
SELECT * FROM hap_vt() LIMIT 3

-- case: haplotypes vt status counts
SELECT prediction_status, prediction_reason, haplotype_impact, count(*) AS n FROM hap_vt() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt nested carriers
SELECT transcript_index, c.sample_index, c.haplotype_lane, c.prediction_status, c.haplotype_impact, c.haplotype_consequences, c.nmd_prediction, c.nmd_stop_position, c.nmd_junction_position FROM (SELECT transcript_index, unnest(carrier_predictions) AS c FROM hap_vt()) ORDER BY ALL

-- case: haplotypes sc full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x, count(DISTINCT transcript_index) AS transcripts FROM (SELECT * FROM hap_sc()) h

-- case: haplotypes sc first rows
SELECT * FROM hap_sc() LIMIT 3

-- case: haplotypes sc status counts
SELECT prediction_status, prediction_reason, haplotype_impact, count(*) AS n FROM hap_sc() GROUP BY ALL ORDER BY ALL

-- case: haplotypes sc nested carriers
SELECT transcript_index, c.sample_index, c.haplotype_lane, c.prediction_status, c.haplotype_impact, c.haplotype_consequences, c.nmd_prediction, c.nmd_stop_position, c.nmd_junction_position FROM (SELECT transcript_index, unnest(carrier_predictions) AS c FROM hap_sc()) ORDER BY ALL

-- case: haplotypes fr full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x, count(DISTINCT transcript_index) AS transcripts FROM (SELECT * FROM hap_fr()) h

-- case: haplotypes fr first rows
SELECT * FROM hap_fr() LIMIT 3

-- case: haplotypes fr status counts
SELECT prediction_status, prediction_reason, haplotype_impact, count(*) AS n FROM hap_fr() GROUP BY ALL ORDER BY ALL

-- case: haplotypes fr nested carriers
SELECT transcript_index, c.sample_index, c.haplotype_lane, c.prediction_status, c.haplotype_impact, c.haplotype_consequences, c.nmd_prediction, c.nmd_stop_position, c.nmd_junction_position FROM (SELECT transcript_index, unnest(carrier_predictions) AS c FROM hap_fr()) ORDER BY ALL

-- case: haplotypes ss full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x, count(DISTINCT transcript_index) AS transcripts FROM (SELECT * FROM hap_ss()) h

-- case: haplotypes ss first rows
SELECT * FROM hap_ss() LIMIT 3

-- case: haplotypes ss status counts
SELECT prediction_status, prediction_reason, haplotype_impact, count(*) AS n FROM hap_ss() GROUP BY ALL ORDER BY ALL

-- case: haplotypes ss nested carriers
SELECT transcript_index, c.sample_index, c.haplotype_lane, c.prediction_status, c.haplotype_impact, c.haplotype_consequences, c.nmd_prediction, c.nmd_stop_position, c.nmd_junction_position FROM (SELECT transcript_index, unnest(carrier_predictions) AS c FROM hap_ss()) ORDER BY ALL

-- case: haplotypes nm full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x, count(DISTINCT transcript_index) AS transcripts FROM (SELECT * FROM hap_nm()) h

-- case: haplotypes nm first rows
SELECT * FROM hap_nm() LIMIT 3

-- case: haplotypes nm status counts
SELECT prediction_status, prediction_reason, haplotype_impact, count(*) AS n FROM hap_nm() GROUP BY ALL ORDER BY ALL

-- case: haplotypes nm nested carriers
SELECT transcript_index, c.sample_index, c.haplotype_lane, c.prediction_status, c.haplotype_impact, c.haplotype_consequences, c.nmd_prediction, c.nmd_stop_position, c.nmd_junction_position FROM (SELECT transcript_index, unnest(carrier_predictions) AS c FROM hap_nm()) ORDER BY ALL

-- case: haplotypes vt schema
SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM hap_vt())

-- case: haplotypes vt projected list columns
SELECT transcript_index, len(coding_blocks) AS blocks, len(contributors) AS contributors, len(cds_differences) AS cds_diffs, len(protein_differences) AS protein_diffs, len(normalized_edits) AS edits, nmd_contributors FROM hap_vt() ORDER BY transcript_index

-- case: haplotypes vt filter after the scan
SELECT count(*) AS n FROM hap_vt() WHERE prediction_status = 'predicted'

-- case: haplotypes vt limit one
SELECT transcript_index FROM hap_vt() LIMIT 1

-- case: haplotypes nmd distinct reasons
SELECT DISTINCT nmd_prediction, nmd_rule FROM hap_nm() ORDER BY ALL

-- Discovery and coding_calls fixtures: the model of test/sql/duckvep_coding_calls.test (region names c0 to c2).

-- fixture: discovery tx
CREATE TABLE cc_tx AS SELECT * FROM (VALUES (0::UINTEGER, 0::UINTEGER, 100::UBIGINT, 132::UBIGINT, 1::TINYINT, 0::UINTEGER, 3::UBIGINT, 103::UBIGINT, 129::UBIGINT, 'ATGAAACCCGGGTTTTAA'::BLOB, 1::UTINYINT, 'AAA'::BLOB, 'CCC'::BLOB), (1::UINTEGER, 0::UINTEGER, 104::UBIGINT, 131::UBIGINT, 1::TINYINT, 1::UINTEGER, 3::UBIGINT, 104::UBIGINT, 131::UBIGINT, 'ATGAAACCCTAA'::BLOB, 1::UTINYINT, ''::BLOB, ''::BLOB), (2::UINTEGER, 0::UINTEGER, 105::UBIGINT, 115::UBIGINT, 1::TINYINT, 2::UINTEGER, 0::UBIGINT, NULL::UBIGINT, NULL::UBIGINT, NULL::BLOB, NULL::UTINYINT, NULL::BLOB, NULL::BLOB), (3::UINTEGER, 1::UINTEGER, 200::UBIGINT, 232::UBIGINT, -1::TINYINT, 3::UINTEGER, 3::UBIGINT, 203::UBIGINT, 229::UBIGINT, 'ATGAAACCCGGGTTTTAA'::BLOB, 1::UTINYINT, 'GGG'::BLOB, 'TTT'::BLOB), (4::UINTEGER, 2::UINTEGER, 300::UBIGINT, 349::UBIGINT, 1::TINYINT, 4::UINTEGER, 3::UBIGINT, 300::UBIGINT, 349::UBIGINT, 'ATGAAACCCGGGTTTAAACC'::BLOB, 1::UTINYINT, ''::BLOB, ''::BLOB), (5::UINTEGER, 2::UINTEGER, 400::UBIGINT, 449::UBIGINT, -1::TINYINT, 5::UINTEGER, 3::UBIGINT, 400::UBIGINT, 449::UBIGINT, 'ATGAAACCCGGGTTTAAACC'::BLOB, 1::UTINYINT, ''::BLOB, ''::BLOB) ) t(transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence)

-- fixture: discovery exons
CREATE TABLE cc_ex AS SELECT * FROM (VALUES (0::UINTEGER, 100::UBIGINT, 110::UBIGINT, 1::UBIGINT, 11::UBIGINT, -1::TINYINT, 2::TINYINT), (0::UINTEGER, 120::UBIGINT, 132::UBIGINT, 12::UBIGINT, 24::UBIGINT, 2::TINYINT, -1::TINYINT), (1::UINTEGER, 104::UBIGINT, 108::UBIGINT, 1::UBIGINT, 5::UBIGINT, 0::TINYINT, 2::TINYINT), (1::UINTEGER, 125::UBIGINT, 131::UBIGINT, 6::UBIGINT, 12::UBIGINT, 2::TINYINT, 0::TINYINT), (2::UINTEGER, 105::UBIGINT, 115::UBIGINT, 1::UBIGINT, 11::UBIGINT, -1::TINYINT, -1::TINYINT), (3::UINTEGER, 220::UBIGINT, 232::UBIGINT, 1::UBIGINT, 13::UBIGINT, -1::TINYINT, 1::TINYINT), (3::UINTEGER, 200::UBIGINT, 210::UBIGINT, 14::UBIGINT, 24::UBIGINT, 1::TINYINT, -1::TINYINT), (4::UINTEGER, 300::UBIGINT, 309::UBIGINT, 1::UBIGINT, 10::UBIGINT, 0::TINYINT, 1::TINYINT), (4::UINTEGER, 340::UBIGINT, 349::UBIGINT, 11::UBIGINT, 20::UBIGINT, 1::TINYINT, 2::TINYINT), (5::UINTEGER, 440::UBIGINT, 449::UBIGINT, 1::UBIGINT, 10::UBIGINT, 0::TINYINT, 1::TINYINT), (5::UINTEGER, 400::UBIGINT, 409::UBIGINT, 11::UBIGINT, 20::UBIGINT, 1::TINYINT, 2::TINYINT) ) e(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase)

-- fixture: discovery alleles
CREATE TABLE dc_alleles AS SELECT * FROM (VALUES ('snv', 'A', 'C'), ('mnv', 'AG', 'CT'), ('del1', 'AC', 'A'), ('del2', 'ACG', 'A'), ('del3', 'ACGT', 'A'), ('del5', 'ACGTAC', 'A'), ('ins1', 'A', 'AT'), ('ins2', 'A', 'ATG'), ('indel', 'ACG', 'AT'), ('delins', 'ACG', 'TT')) v(shape, reference, alternate)

-- fixture: discovery events
CREATE TABLE dc_ev AS SELECT (row_number() OVER ())::UBIGINT event_index, r::UINTEGER seq_region, p::UBIGINT AS position, reference, alternate, NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change, NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position FROM range(0, 3) a(r), range(90, 440) b(p), dc_alleles

-- fixture-v1: discovery model
SELECT loaded FROM duckvep_model_load('cc', 'SELECT i::UINTEGER seq_region, 1000::UBIGINT sequence_length, (''c'' || i) seq_region_name FROM range(3) t(i)', 'SELECT * FROM cc_tx ORDER BY seq_region, transcript_start, transcript_index', 'SELECT * FROM cc_ex ORDER BY transcript_index, exon_cdna_start', transcript_coverage_complete := TRUE)

-- fixture-v2: discovery stage regions
COPY (
SELECT i::UINTEGER seq_region, 1000::UBIGINT sequence_length, ('c' || i) seq_region_name FROM range(3) t(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'cc', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: discovery stage transcripts
COPY (
SELECT * FROM cc_tx ORDER BY seq_region, transcript_start, transcript_index
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'cc', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: discovery stage exons
COPY (
SELECT * FROM cc_ex ORDER BY transcript_index, exon_cdna_start
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'cc', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: discovery publish
SELECT duckvep_model_publish('cc', {'transcript_coverage_complete': true})

-- job: hap_cc_cases cc
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf')

-- job: hap_cc_multi cc
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/multi.vcf')

-- job: hap_cc_demo cc
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/demo.vcf')

-- job: hap_cc_bcf cc
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.bcf')


-- Policy, limit and error cases of test/sql/duckvep_haplotypes.test on its model "hap" (calls tables hap_calls, hap_raw_overlap, hap_raw).

-- fixture: policy hap_tx
CREATE TABLE hap_tx AS SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region, 100::UBIGINT transcript_start, 111::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 100::UBIGINT cds_start, 111::UBIGINT cds_end, 'AAAAAAAAAAAA'::BLOB cds_sequence, 1::UTINYINT codon_table

-- fixture: policy hap_exons
CREATE TABLE hap_exons AS SELECT 0::UINTEGER transcript_index, 100::UBIGINT exon_start, 111::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase

-- fixture: policy hap_calls
CREATE TABLE hap_calls AS SELECT event_index, 0::UINTEGER seq_region, position, 'A' reference, alternate, 1::UINTEGER alt_index, 0::UINTEGER transcript_index, 0::UINTEGER sample_index, alleles, [true,true] phase_before, phase_set FROM (VALUES (1::UBIGINT,100::UBIGINT,'C',[1,1],NULL::BIGINT), (2,101,'G',[1,0],10), (3,102,'C',[0,1],20)) v(event_index,position,alternate,alleles,phase_set)

-- fixture: policy hap_raw_overlap
CREATE TABLE hap_raw_overlap AS SELECT event_index,0 seq_region,position,reference, alternates,0 transcript_index,0 sample_index,gt FROM (VALUES (1,100,'AAA',['CAA'],'0|1'),(2,101,'A',['G'],'1|1')) v(event_index,position,reference,alternates,gt)

-- fixture: policy hap_raw
CREATE TABLE hap_raw AS SELECT i event_index,0 seq_region,99+i AS position,'A' AS reference, CASE i WHEN 1 THEN ['C','T'] ELSE ['G'] END alternates,0 transcript_index,0 sample_index, CASE i WHEN 1 THEN '.|1' ELSE '1|1' END gt FROM range(1,3) r(i)

-- fixture-v1: policy model
SELECT loaded FROM duckvep_model_load('hap', 'SELECT i::UINTEGER seq_region FROM range(2) r(i)', 'SELECT * FROM hap_tx', 'SELECT * FROM hap_exons')

-- fixture-v2: policy stage regions
COPY (
SELECT i::UINTEGER seq_region FROM range(2) r(i)
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'hap', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: policy stage transcripts
COPY (
SELECT * FROM hap_tx
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'hap', RELATION 'transcripts', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: policy stage exons
COPY (
SELECT * FROM hap_exons
) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL 'hap', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)

-- fixture-v2: policy publish
SELECT duckvep_model_publish('hap')

-- job: hap_p01 hap
SELECT * FROM hap_calls ORDER BY event_index DESC

-- job: hap_p02 hap
SELECT * FROM hap_calls

-- job: hap_p03 hap max_alignment_cells=39
SELECT 1 event_index,0 seq_region,105 AS position,'A' AS reference,'AA' alternate,
 1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set

-- job: hap_p04 hap max_alignment_cells=38
SELECT 1 event_index,0 seq_region,105 AS position,'A' AS reference,'AA' alternate,
 1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set

-- job: hap_p05 hap max_leaf_differences=1
SELECT * FROM hap_calls

-- job: hap_p06 hap max_alignment_cells=1 max_leaf_differences=1
SELECT event_index,0 seq_region,position,reference,alternate,
 1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set
 FROM (VALUES (1,101,'A','AA'),(2,108,'AA','A'))
 v(event_index,position,reference,alternate)

-- job: hap_p07 hap phase_policy='vep116_compat'
SELECT * FROM hap_calls

-- job: hap_p08 hap
SELECT * REPLACE(CASE WHEN event_index=3 THEN [false,false] ELSE phase_before END AS phase_before) FROM hap_calls

-- job: hap_p09 hap
SELECT * REPLACE([NULL,NULL]::INTEGER[] AS alleles) FROM hap_calls

-- job: hap_p10 hap
SELECT * FROM hap_calls UNION ALL SELECT * FROM hap_calls

-- job: hap_p11 hap
SELECT * REPLACE([0,0] AS alleles) FROM hap_calls UNION ALL SELECT * REPLACE([0,0] AS alleles) FROM hap_calls

-- job: hap_p12 hap
SELECT * FROM hap_calls UNION ALL SELECT * REPLACE(1 AS sample_index,position+1 AS position) FROM hap_calls

-- job: hap_p13 hap phase_policy='vep116_compat'
SELECT * FROM hap_calls UNION ALL SELECT * REPLACE(1 AS sample_index,2 AS alt_index) FROM hap_calls

-- job: hap_p14 hap
SELECT * REPLACE(NULL::UINTEGER AS sample_index) FROM hap_calls

-- job: hap_p15 hap
SELECT * REPLACE(NULL::UBIGINT AS event_index) FROM hap_calls

-- job: hap_p16 hap
SELECT * REPLACE(CASE WHEN event_index=3 THEN [1] ELSE alleles END AS alleles, NULL::BOOLEAN[] AS phase_before) FROM hap_calls

-- job: hap_p17 hap
SELECT * FROM hap_calls UNION ALL SELECT * REPLACE(1 AS sample_index,[1] AS alleles,[true] AS phase_before) FROM hap_calls

-- job: hap_p18 hap workspace_limit=1
SELECT * FROM hap_calls

-- job: hap_p19 hap max_ploidy=1
SELECT * FROM hap_calls

-- job: hap_p20 hap
SELECT * FROM hap_calls WHERE false

-- job: hap_p21 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * FROM hap_raw_overlap ORDER BY event_index DESC

-- job: hap_p22 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * FROM hap_raw_overlap

-- job: hap_p23 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(CASE event_index WHEN 1 THEN '.|.' ELSE gt END AS gt)
 FROM hap_raw_overlap

-- job: hap_p24 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(100 AS position,'A' AS reference,'1|1' AS gt,
 CASE event_index WHEN 1 THEN ['C'] ELSE ['G'] END AS alternates) FROM hap_raw_overlap

-- job: hap_p25 hap input_mode='source_records' phase_policy='vep116_compat' max_leaf_edits=1
SELECT * FROM hap_raw_overlap

-- job: hap_p26 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT event_index,0 seq_region,position,reference,alternates,0 transcript_index,0 sample_index,gt
 FROM (VALUES (1,100,'A',['C'],'1|1'),(2,100,'A',['G'],'1|1'),
 (3,109,'A',['T'],'0|0')) r(event_index,position,reference,alternates,gt) ORDER BY event_index DESC

-- job: hap_p27 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT event_index,0 seq_region,100 AS position,'A' AS reference,['C'] alternates,
 0 transcript_index,0 sample_index,gt FROM (VALUES (1,'1|0'),(2,'0|1')) r(event_index,gt)

-- job: hap_p28 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT event_index,0 seq_region,100 AS position,'A' AS reference,alternates,
 0 transcript_index,0 sample_index,'1|1' gt
 FROM (VALUES (1,['C','G']),(2,['C','T'])) r(event_index,alternates)

-- job: hap_p29 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * FROM hap_raw ORDER BY event_index DESC

-- job: hap_p30 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * FROM hap_raw

-- job: hap_p31 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(CASE event_index WHEN 1 THEN '|0|1' ELSE gt END AS gt) FROM hap_raw

-- job: hap_p32 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(CASE event_index WHEN 1 THEN '0|1' ELSE gt END AS gt) FROM hap_raw

-- job: hap_p33 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(CASE event_index WHEN 1 THEN '.' ELSE gt END AS gt) FROM hap_raw

-- job: hap_p34 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE('T' AS reference,'.' AS gt) FROM hap_raw WHERE event_index=1

-- job: hap_p35 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * FROM hap_raw UNION ALL SELECT * FROM hap_raw

-- job: hap_p36 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(['C',NULL] AS alternates) FROM hap_raw

-- job: hap_p37 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(NULL AS alternates) FROM hap_raw

-- job: hap_p38 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE('3|1' AS gt) FROM hap_raw

-- job: hap_p39 hap max_ploidy=1 input_mode='source_records' phase_policy='vep116_compat'
SELECT * FROM hap_raw

-- job: hap_p40 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * EXCLUDE(i) REPLACE(s.i AS sample_index) FROM hap_raw CROSS JOIN range(4097) s(i)

-- job: hap_p41 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(NULL AS event_index) FROM hap_raw

-- job: hap_p42 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(NULL AS sample_index) FROM hap_raw

-- job: hap_p43 hap input_mode='source_records' phase_policy='vep116_compat'
SELECT * REPLACE(NULL AS transcript_index) FROM hap_raw

-- job: hap_p44 hap max_active_events=1
SELECT * FROM hap_calls

-- job: hap_p45 hap max_sequence_bases=11
SELECT * FROM hap_calls

-- job: hap_p46 hap
SELECT event_index,0 seq_region,100 AS position,'A' AS reference,alternate,
  1 alt_index,0 transcript_index,sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set
  FROM (VALUES (1,0,'TA'),(2,1,'AT')) v(event_index,sample_index,alternate)

-- job: hap_p47 hap
SELECT event_index,0 seq_region,100 AS position,'A' AS reference,alternate,
  1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set
  FROM (VALUES (1,'C'),(2,'G')) v(event_index,alternate)

-- case: haplotypes policy 01
SELECT * FROM hap_p01()

-- case: haplotypes policy 02
SELECT * FROM hap_p02()

-- case: haplotypes policy 03
SELECT * FROM hap_p03()

-- case: haplotypes policy 04 error
SELECT * FROM hap_p04()

-- case: haplotypes policy 05 error
SELECT * FROM hap_p05()

-- case: haplotypes policy 06
SELECT * FROM hap_p06()

-- case: haplotypes policy 07
SELECT * FROM hap_p07()

-- case: haplotypes policy 08
SELECT * FROM hap_p08()

-- case: haplotypes policy 09
SELECT * FROM hap_p09()

-- case: haplotypes policy 10 error
SELECT * FROM hap_p10()

-- case: haplotypes policy 11 error
SELECT * FROM hap_p11()

-- case: haplotypes policy 12 error
SELECT * FROM hap_p12()

-- case: haplotypes policy 13 error
SELECT * FROM hap_p13()

-- case: haplotypes policy 14 error
SELECT * FROM hap_p14()

-- case: haplotypes policy 15 error
SELECT * FROM hap_p15()

-- case: haplotypes policy 16 error
SELECT * FROM hap_p16()

-- case: haplotypes policy 17
SELECT * FROM hap_p17()

-- case: haplotypes policy 18 error
SELECT * FROM hap_p18()

-- case: haplotypes policy 19 error
SELECT * FROM hap_p19()

-- case: haplotypes policy 20
SELECT * FROM hap_p20()

-- case: haplotypes policy 21
SELECT * FROM hap_p21()

-- case: haplotypes policy 22
SELECT * FROM hap_p22()

-- case: haplotypes policy 23
SELECT * FROM hap_p23()

-- case: haplotypes policy 24
SELECT * FROM hap_p24()

-- case: haplotypes policy 25 error
SELECT * FROM hap_p25()

-- case: haplotypes policy 26
SELECT * FROM hap_p26()

-- case: haplotypes policy 27
SELECT * FROM hap_p27()

-- case: haplotypes policy 28
SELECT * FROM hap_p28()

-- case: haplotypes policy 29
SELECT * FROM hap_p29()

-- case: haplotypes policy 30
SELECT * FROM hap_p30()

-- case: haplotypes policy 31
SELECT * FROM hap_p31()

-- case: haplotypes policy 32
SELECT * FROM hap_p32()

-- case: haplotypes policy 33
SELECT * FROM hap_p33()

-- case: haplotypes policy 34
SELECT * FROM hap_p34()

-- case: haplotypes policy 35 error
SELECT * FROM hap_p35()

-- case: haplotypes policy 36 error
SELECT * FROM hap_p36()

-- case: haplotypes policy 37 error
SELECT * FROM hap_p37()

-- case: haplotypes policy 38 error
SELECT * FROM hap_p38()

-- case: haplotypes policy 39 error
SELECT * FROM hap_p39()

-- case: haplotypes policy 40
SELECT * FROM hap_p40()

-- case: haplotypes policy 41 error
SELECT * FROM hap_p41()

-- case: haplotypes policy 42 error
SELECT * FROM hap_p42()

-- case: haplotypes policy 43 error
SELECT * FROM hap_p43()

-- case: haplotypes policy 44 error
SELECT * FROM hap_p44()

-- case: haplotypes policy 45 error
SELECT * FROM hap_p45()

-- case: haplotypes policy 46
SELECT * FROM hap_p46()

-- case: haplotypes policy 47
SELECT * FROM hap_p47()

-- job: hap_vt_compat vt phase_policy='vep116_compat'
SELECT * FROM vt_calls

-- job: hap_vt_hgvs vt hgvs=true
SELECT * FROM vt_calls

-- job: hap_vt_limit vt max_active_events=1
SELECT * FROM vt_calls

-- job: hap_vt_dup vt
SELECT * FROM vt_calls UNION ALL SELECT * FROM vt_calls

-- job: hap_vt_null vt
SELECT * REPLACE(NULL::UINTEGER AS sample_index) FROM vt_calls

-- job: hap_vt_ploidy vt max_ploidy=1
SELECT * FROM vt_calls

-- job: hap_vt_sets vt
SELECT * REPLACE(CASE WHEN event_index % 2 = 0 THEN 7 ELSE 9 END::BIGINT AS phase_set) FROM vt_calls

-- job: hap_vt_unphased vt
SELECT * REPLACE([false, false] AS phase_before) FROM vt_calls

-- job: hap_vt_missing vt
SELECT * REPLACE(CASE WHEN event_index % 3 = 0 THEN [NULL, 1]::INTEGER[] ELSE alleles END AS alleles) FROM vt_calls

-- job: hap_vt_triploid vt
SELECT * REPLACE([1, 0, 1] AS alleles, [false, true, true] AS phase_before) FROM vt_calls

-- job: hap_vt_shuffled vt
SELECT * FROM vt_calls ORDER BY hash(event_index)

-- job: hap_nm_hgvs nm hgvs=true
SELECT * FROM nm_calls

-- job: hap_nm_compat nm phase_policy='vep116_compat'
SELECT * FROM nm_calls

-- job: hap_nm_events nm max_leaf_events=1
SELECT * FROM nm_calls

-- job: hap_nm_empty nm
SELECT * FROM nm_calls WHERE false

-- case: haplotypes vt compat full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_compat()) h

-- case: haplotypes vt compat status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_compat() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt hgvs full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_hgvs()) h

-- case: haplotypes vt hgvs status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_hgvs() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt shuffled full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_shuffled()) h

-- case: haplotypes vt shuffled status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_shuffled() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt sets full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_sets()) h

-- case: haplotypes vt sets status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_sets() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt unphased full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_unphased()) h

-- case: haplotypes vt unphased status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_unphased() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt missing full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_missing()) h

-- case: haplotypes vt missing status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_missing() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt triploid full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_vt_triploid()) h

-- case: haplotypes vt triploid status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_vt_triploid() GROUP BY ALL ORDER BY ALL

-- case: haplotypes nm hgvs full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_nm_hgvs()) h

-- case: haplotypes nm hgvs status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_nm_hgvs() GROUP BY ALL ORDER BY ALL

-- case: haplotypes nm compat full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_nm_compat()) h

-- case: haplotypes nm compat status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_nm_compat() GROUP BY ALL ORDER BY ALL

-- case: haplotypes nm empty full-row hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_nm_empty()) h

-- case: haplotypes nm empty status counts
SELECT prediction_status, prediction_reason, hgvsp_status, count(*) AS n FROM hap_nm_empty() GROUP BY ALL ORDER BY ALL

-- case: haplotypes vt hgvs rows
SELECT transcript_index, hgvsp, hgvsp_status FROM hap_vt_hgvs() ORDER BY transcript_index

-- case: haplotypes vt limit error
SELECT count(*) FROM hap_vt_limit()

-- case: haplotypes vt duplicate calls error
SELECT count(*) FROM hap_vt_dup()

-- case: haplotypes vt null sample error
SELECT count(*) FROM hap_vt_null()

-- case: haplotypes vt ploidy error
SELECT count(*) FROM hap_vt_ploidy()

-- case: haplotypes nm leaf events error
SELECT count(*) FROM hap_nm_events()

-- case: discovery all events
SELECT count(*) AS events, sum(len(t)) AS pairs, sum(hash(t)::HUGEINT) AS h FROM (SELECT duckvep_coding_transcripts('cc', seq_region, position, reference, alternate) AS t FROM dc_ev)

-- case: discovery unnested pairs
SELECT event_index, unnest(duckvep_coding_transcripts('cc', seq_region, position, reference, alternate)) AS transcript_index FROM dc_ev WHERE position BETWEEN 100 AND 112

-- case: discovery named cases
SELECT duckvep_coding_transcripts('cc', 0, 129, 'AGG', 'A') AS a, duckvep_coding_transcripts('cc', 0, 102, 'TAA', 'T') AS b, duckvep_coding_transcripts('cc', 0, 110, 'AGG', 'A') AS c, duckvep_coding_transcripts('cc', 0, 119, 'AGG', 'A') AS d, duckvep_coding_transcripts('cc', 2, 449, 'AGG', 'A') AS e, duckvep_coding_transcripts('cc', 1, 229, 'T', 'TA') AS f

-- case: discovery non literal alleles
SELECT duckvep_coding_transcripts('cc', 0, 103, 'A', 'C') AS snv, duckvep_coding_transcripts('cc', 0, 103, 'A', 'A') AS same, duckvep_coding_transcripts('cc', 0, 104, 'A', '<DEL>') AS symbolic, duckvep_coding_transcripts('cc', 0, 104, 'A', '*') AS star, duckvep_coding_transcripts('cc', 7, 104, 'A', 'C') AS unknown_region

-- case: discovery null arguments
SELECT duckvep_coding_transcripts('cc', NULL, 103, 'A', 'C') AS a, duckvep_coding_transcripts(NULL, 0, 103, 'A', 'C') AS b, duckvep_coding_transcripts('cc', 0, NULL, 'A', 'C') AS c, duckvep_coding_transcripts('cc', 0, 103, NULL, 'C') AS d, duckvep_coding_transcripts('cc', 0, 103, 'A', NULL) AS e

-- case: discovery constant and varying model columns
SELECT sum(len(duckvep_coding_transcripts(m, 0, 100 + i, 'A', 'C'))) AS n FROM (SELECT 'cc' AS m, i FROM range(40) t(i))

-- case: discovery unknown model error
SELECT duckvep_coding_transcripts('nope', 0, 103, 'A', 'C')

-- case: discovery bad position error
SELECT duckvep_coding_transcripts('cc', 0, 0, 'A', 'C')

-- case: discovery long allele error
SELECT duckvep_coding_transcripts('cc', 0, 103, repeat('A', 65536), 'C')

-- case: coding_calls cases.vcf rows
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf')

-- case: coding_calls multi.vcf rows
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/multi.vcf')

-- case: coding_calls cases.bcf rows
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.bcf')

-- case: coding_calls cases.vcf.gz rows
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf.gz')

-- case: coding_calls demo.vcf rows
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/demo.vcf')

-- case: coding_calls many alts error
SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/many_alts.vcf')

-- case: coding_calls schema
SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf'))

-- case: coding_calls hash
SELECT count(*) AS n, sum(hash(c)::HUGEINT) AS total, bit_xor(hash(c)) AS x FROM (SELECT * FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf')) c

-- case: coding_calls filter and limit
SELECT event_index, transcript_index FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf') WHERE sample_index = 0 LIMIT 5

-- case: coding_calls self join
SELECT count(*) AS n FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf') a JOIN duckvep_coding_calls('cc', 'test/data/coding_calls/multi.vcf') b USING (event_index)

-- case: coding_calls plain versus gzip versus bcf
SELECT (SELECT sum(hash(c)::HUGEINT) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf') c) = (SELECT sum(hash(c)::HUGEINT) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf.gz') c) AS gz_equal, (SELECT sum(hash(c)::HUGEINT) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.vcf') c) = (SELECT sum(hash(c)::HUGEINT) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/cases.bcf') c) AS bcf_equal

-- case: coding_calls bad position error
SELECT count(*) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/bad_pos.vcf')

-- case: coding_calls no GT header error
SELECT count(*) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/no_gt_header.vcf')

-- case: coding_calls no GT record error
SELECT count(*) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/no_gt_record.vcf')

-- case: coding_calls no samples error
SELECT count(*) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/no_samples.vcf')

-- case: coding_calls not a vcf error
SELECT count(*) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/not_a_vcf.txt')

-- case: coding_calls missing file error
SELECT count(*) FROM duckvep_coding_calls('cc', 'test/data/coding_calls/missing.vcf')

-- case: coding_calls unknown model error
SELECT count(*) FROM duckvep_coding_calls('nope', 'test/data/coding_calls/cases.vcf')

-- case: coding_calls empty path error
SELECT count(*) FROM duckvep_coding_calls('cc', '')

-- case: haplotypes from coding_calls cases error
SELECT * FROM hap_cc_cases()

-- case: haplotypes from coding_calls multi rows
SELECT * FROM hap_cc_multi()

-- case: haplotypes from coding_calls bcf error
SELECT * FROM hap_cc_bcf()

-- case: haplotypes from coding_calls demo hash
SELECT count(*) AS n, sum(hash(h)::HUGEINT) AS total, bit_xor(hash(h)) AS x FROM (SELECT * FROM hap_cc_demo()) h

-- case: haplotypes from coding_calls demo rows
SELECT * FROM hap_cc_demo()
