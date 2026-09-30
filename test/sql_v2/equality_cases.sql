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
