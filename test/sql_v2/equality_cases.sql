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
