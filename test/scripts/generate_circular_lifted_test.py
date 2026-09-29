#!/usr/bin/env python3
"""Generate test/sql/duckvep_circular_lifted.test.

A 300 bp circular region carries plus- and minus-strand coding and non-coding
transcripts and regulation features, several of which cross the origin. The
same world is rotated by several offsets in SQL: reference, model and events
move together. Annotation through the lifted model must agree for every
offset once rows are keyed by stable identifier. Regenerate with:

    python3 test/scripts/generate_circular_lifted_test.py > test/sql/duckvep_circular_lifted.test
"""
import random
import sys

L = 3000
ROTATIONS = [0, 611, 1498, 2203, 2999]
SHIFT = L - 900  # the layout below was designed on a 900 bp circle with the origin after 900
rng = random.Random(20260929)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


ref = [rng.choice("ACGT") for _ in range(L)]


def pos(u):
    return (u - 1) % L + 1



def base(u):
    return ref[(u - 1) % L]


def setbase(u, b):
    ref[(u - 1) % L] = b


# Transcripts in unrolled coordinates: exons ascending genomic, listed as
# (u_start, u_end); rank order is ascending on + and descending on -.
TX = [
    # Exon lists are unrolled and ascending; cds is the 1-based cDNA CDS range.
    # T_SW: one wrapped exon, CDS starting at its first and ending at its last base.
    dict(name="T_SW", gene="G_SW", strand=1, biotype="protein_coding",
         exons=[(885, 923)], cds=(1, 39)),
    # T_PW: exon ends on base L, the intron begins after the origin.
    dict(name="T_PW", gene="G_PW", strand=1, biotype="protein_coding",
         exons=[(860, 900), (960, 1010)], cds=(45, 89)),
    # T_MW: wrapped 5' UTR exon on the minus strand, CDS in the low exon.
    dict(name="T_MW", gene="G_MW", strand=-1, biotype="protein_coding",
         exons=[(700, 760), (870, 930)], cds=(65, 121)),
    dict(name="T_NC", gene="G_NC", strand=-1, biotype="lncRNA",
         exons=[(890, 910), (950, 980)], cds=None),
    # T_UW: wrapped 3' UTR exon.
    dict(name="T_UW", gene="G_UW", strand=1, biotype="protein_coding",
         exons=[(200, 260), (280, 310), (880, 905)], cds=(10, 90)),
    # T_NB and T_NF are unwrapped, but their 5' flanks cross the origin.
    dict(name="T_NB", gene="G_NB", strand=1, biotype="protein_coding",
         exons=[(26, 58)], cds=(3, 32)),
    dict(name="T_NF", gene="G_NF", strand=-1, biotype="protein_coding",
         exons=[(852, 880)], cds=(2, 28)),
    # T_MI: minus strand, low exon starts on base 1 so the intron ends on base L.
    dict(name="T_MI", gene="G_MI", strand=-1, biotype="protein_coding",
         exons=[(600, 640), (901, 940)], cds=(44, 79)),
]
FEATURES = [
    dict(name="R_WRAP", kind="reg", u=(880, 915)),
    dict(name="M_WRAP", kind="motif", u=(897, 905)),
    dict(name="R_MID", kind="reg", u=(300, 330)),
    dict(name="R_NEAR", kind="reg", u=(10, 20)),
]


# T_NB and R_NEAR sit just after the origin and keep their small coordinates.
AFTER_ORIGIN = {"T_NB", "R_NEAR"}
for tx in TX:
    d = 0 if tx["name"] in AFTER_ORIGIN else SHIFT
    tx["exons"] = [(a + d, b + d) for a, b in tx["exons"]]
for f in FEATURES:
    d = 0 if f["name"] in AFTER_ORIGIN else SHIFT
    f["u"] = (f["u"][0] + d, f["u"][1] + d)
# A coding transcript and a feature far from the origin, and the same shapes
# elsewhere, so the linear oracle also covers ordinary positions.
TX.append(dict(name="T_MID", gene="G_MID", strand=1, biotype="protein_coding",
               exons=[(1400, 1440), (1500, 1560)], cds=(6, 92)))
FEATURES[2]["u"] = (1600, 1630)  # R_MID
# A regulatory region over all but one base wraps in every rotation except the one
# that puts the origin next to it, so each rotated model executes lifted even where
# no transcript wraps.
FEATURES.append(dict(name="R_GUARD", kind="reg", u=(1000, 1000 + L - 2)))


def cdna_positions(tx):
    """Unrolled coordinate of every cDNA base of a transcript."""
    order = tx["exons"] if tx["strand"] > 0 else list(reversed(tx["exons"]))
    out = []
    for s, e in order:
        rng_ = range(s, e + 1) if tx["strand"] > 0 else range(e, s - 1, -1)
        out.extend(rng_)
    return out


NONSTOP = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
           if a + b + c not in ("TAA", "TAG", "TGA")]
for tx in TX:
    tx["cdna_u"] = cdna_positions(tx)
    if tx["cds"]:
        a, b = tx["cds"]
        codons = (b - a + 1) // 3
        seq = "ATG" + "".join(rng.choice(NONSTOP) for _ in range(codons - 2)) + "TAA"
        tx["cds_seq"] = seq
        for i, ch in enumerate(seq):
            u = tx["cdna_u"][a - 1 + i]
            setbase(u, ch if tx["strand"] > 0 else COMP[ch])
seen = {}
for tx in TX:
    if tx["cds"]:
        a, b = tx["cds"]
        for c in range(a, b + 1):
            u = pos(tx["cdna_u"][c - 1])
            assert u not in seen, (tx["name"], seen[u], u)
            seen[u] = tx["name"]
# A homopolymer across the origin inside the T_SW coding exon would need its
# codons rewritten; the designed CDS is left as generated and repeats arise by chance.
# Splice consensus: GT..AG on the transcript strand (intron ends), where the
# introns are not shared with coding sequence.
for tx in TX:
    exons = tx["exons"]
    for i in range(len(exons) - 1):
        lo = exons[i][1] + 1
        hi = exons[i + 1][0] - 1
        if hi - lo + 1 < 6:
            continue
        if tx["strand"] > 0:
            setbase(lo, "G"); setbase(lo + 1, "T"); setbase(hi - 1, "A"); setbase(hi, "G")
        else:
            setbase(lo, "C"); setbase(lo + 1, "T"); setbase(hi - 1, "A"); setbase(hi, "C")


def sqlstr(s):
    return "'" + s.replace("'", "''") + "'"


lines = []


def emit(kind, sql, expect=None, cols=None):
    if kind == "ok":
        lines.append("statement ok\n" + sql.strip() + "\n")
    elif kind == "error":
        lines.append("statement error\n" + sql.strip() + "\n----\n" + expect + "\n")
    else:
        lines.append("query " + cols + "\n" + sql.strip() + "\n----\n" + expect.strip() + "\n")


lines.append("""# name: test/sql/duckvep_circular_lifted.test
# description: Lifted-interval execution gives identical annotation under rotation of a circular region
# group: [duckvep]
# Generated by test/scripts/generate_circular_lifted_test.py; do not edit by hand.

require duckvep
""")

emit("ok", "SET threads = 1;")
emit("ok", f"CREATE MACRO rot(x, k) AS ((x - 1 + k) % {L}) + 1;")
emit("ok", "CREATE SCHEMA core;")
emit("ok", "CREATE TABLE core.coord_system(coord_system_id BIGINT, species_id BIGINT, name VARCHAR, version VARCHAR, rank BIGINT);")
emit("ok", "INSERT INTO core.coord_system VALUES (1, 1, 'chromosome', 'synthetic', 1);")
emit("ok", "CREATE TABLE core.seq_region(seq_region_id BIGINT, name VARCHAR, coord_system_id BIGINT, length BIGINT);")
emit("ok", f"INSERT INTO core.seq_region VALUES (10, 'C', 1, {L});")
emit("ok", "CREATE TABLE core.attrib_type(attrib_type_id BIGINT, code VARCHAR);")
emit("ok", "INSERT INTO core.attrib_type VALUES (1, 'circular_seq');")
emit("ok", "CREATE TABLE core.seq_region_attrib(seq_region_id BIGINT, attrib_type_id BIGINT, value BIGINT);")
emit("ok", "INSERT INTO core.seq_region_attrib VALUES (10, 1, 1);")
emit("ok", "CREATE TABLE core.gene(gene_id BIGINT, biotype VARCHAR, is_current BIGINT, stable_id VARCHAR, version BIGINT);")
emit("ok", "INSERT INTO core.gene VALUES " + ", ".join(
    f"({20 + i}, {sqlstr(t['biotype'])}, 1, {sqlstr(t['gene'])}, 1)" for i, t in enumerate(TX)) + ";")
emit("ok", """CREATE TABLE core.transcript(transcript_id BIGINT, gene_id BIGINT, seq_region_id BIGINT,
  seq_region_start BIGINT, seq_region_end BIGINT, seq_region_strand BIGINT,
  biotype VARCHAR, is_current BIGINT, stable_id VARCHAR, version BIGINT);""")
emit("ok", """CREATE TABLE core.exon(exon_id BIGINT, seq_region_id BIGINT, seq_region_start BIGINT,
  seq_region_end BIGINT, seq_region_strand BIGINT, phase BIGINT, end_phase BIGINT,
  is_current BIGINT, stable_id VARCHAR);""")
emit("ok", "CREATE TABLE core.exon_transcript(exon_id BIGINT, transcript_id BIGINT, rank BIGINT);")
emit("ok", """CREATE TABLE core.translation(translation_id BIGINT, transcript_id BIGINT,
  seq_start BIGINT, start_exon_id BIGINT, seq_end BIGINT, end_exon_id BIGINT,
  stable_id VARCHAR, version BIGINT);""")
emit("ok", "CREATE TABLE core.transcript_attrib(transcript_id BIGINT, attrib_type_id BIGINT, value VARCHAR);")
emit("ok", "CREATE TABLE core.translation_attrib(translation_id BIGINT, attrib_type_id BIGINT, value VARCHAR);")

tx_rows, exon_rows, et_rows, tr_rows = [], [], [], []
exon_id = 100
for i, tx in enumerate(TX):
    tid = 30 + i
    exons = tx["exons"]
    span_s = exons[0][0]
    span_e = exons[-1][1]
    tx_rows.append(f"({tid}, {20 + i}, 10, {pos(span_s)}, {pos(span_e)}, {tx['strand']}, "
                   f"{sqlstr(tx['biotype'])}, 1, {sqlstr(tx['name'])}, 1)")
    order = exons if tx["strand"] > 0 else list(reversed(exons))
    cum = 0
    ids = []
    a, b = tx["cds"] if tx["cds"] else (0, 0)
    before = 0
    start_exon = end_exon = None
    seq_start = seq_end = 0
    for rank, (s, e) in enumerate(order, 1):
        exon_id += 1
        ids.append(exon_id)
        length = e - s + 1
        cs, ce = cum + 1, cum + length
        phase = end_phase = -1
        if tx["cds"] and ce >= a and cs <= b:
            first, last = max(cs, a), min(ce, b)
            if cs >= a:
                phase = before % 3
            if ce < b:
                end_phase = (before + last - first + 1) % 3
            if cs <= a <= ce:
                start_exon, seq_start = exon_id, a - cs + 1
            if cs <= b <= ce:
                end_exon, seq_end = exon_id, b - cs + 1
            before += last - first + 1
        exon_rows.append(f"({exon_id}, 10, {pos(s)}, {pos(e)}, {tx['strand']}, {phase}, {end_phase}, 1, "
                         f"{sqlstr(tx['name'] + '_E' + str(rank))})")
        et_rows.append(f"({exon_id}, {tid}, {rank})")
        cum = ce
    if tx["cds"]:
        tr_rows.append(f"({50 + i}, {tid}, {seq_start}, {start_exon}, {seq_end}, {end_exon}, "
                       f"{sqlstr('P_' + tx['name'][2:])}, 1)")
emit("ok", "INSERT INTO core.transcript VALUES\n  " + ",\n  ".join(tx_rows) + ";")
emit("ok", "INSERT INTO core.exon VALUES\n  " + ",\n  ".join(exon_rows) + ";")
emit("ok", "INSERT INTO core.exon_transcript VALUES " + ", ".join(et_rows) + ";")
emit("ok", "INSERT INTO core.translation VALUES " + ", ".join(tr_rows) + ";")

emit("ok", "CREATE SCHEMA funcgen;")
emit("ok", "CREATE TABLE funcgen.feature_type(feature_type_id BIGINT, name VARCHAR, so_accession VARCHAR, so_term VARCHAR);")
emit("ok", "INSERT INTO funcgen.feature_type VALUES (70, 'Promoter', 'SO:0000167', 'promoter');")
emit("ok", """CREATE TABLE funcgen.regulatory_feature(regulatory_feature_id BIGINT, feature_type_id BIGINT,
  seq_region_id BIGINT, seq_region_strand BIGINT, seq_region_start BIGINT,
  seq_region_end BIGINT, stable_id VARCHAR, regulatory_build_id BIGINT);""")
emit("ok", """CREATE TABLE funcgen.motif_feature(motif_feature_id BIGINT, binding_matrix_id BIGINT,
  seq_region_id BIGINT, seq_region_start BIGINT, seq_region_end BIGINT,
  seq_region_strand BIGINT, score DOUBLE, stable_id VARCHAR);""")
regs = [(81 + i, f) for i, f in enumerate(FEATURES) if f["kind"] == "reg"]
mots = [(91 + i, f) for i, f in enumerate(FEATURES) if f["kind"] == "motif"]
emit("ok", "INSERT INTO funcgen.regulatory_feature VALUES " + ", ".join(
    f"({i}, 70, 10, 0, {pos(f['u'][0])}, {pos(f['u'][1])}, {sqlstr(f['name'])}, 1)" for i, f in regs) + ";")
emit("ok", "INSERT INTO funcgen.motif_feature VALUES " + ", ".join(
    f"({i}, 90, 10, {pos(f['u'][0])}, {pos(f['u'][1])}, 1, 1.0, {sqlstr(f['name'])})" for i, f in mots) + ";")

emit("ok", "CREATE TABLE ref(chrom VARCHAR, \"start\" BIGINT, \"end\" BIGINT, seq VARCHAR);")
emit("ok", f"INSERT INTO ref VALUES ('C', 0, {L}, {sqlstr(''.join(ref))});")

# Events, defined once in unrolled base coordinates: every position, several
# shapes; long alleles; and one nearly circumferential deletion.
INS = "GTC"
_long = random.Random(7)
LONG_INS = "".join(_long.choice("ACGT") for _ in range(61))
shapes = [
    ("snv", 1, lambda r: COMP[r[0]]),
    ("del1", 2, lambda r: r[0]),
    ("del2", 3, lambda r: r[0]),
    ("del5", 6, lambda r: r[0]),
    ("insA", 1, lambda r: r[0] + INS),
    ("insB", 1, lambda r: INS + r[0]),
    ("mnv3", 3, lambda r: revcomp(r)),
    ("del45", 46, lambda r: r[0]),
    ("ins61", 1, lambda r: r[0] + LONG_INS),
    ("delins", 12, lambda r: "ACGTACGTACGTACGTACGTACGTA"),
]
events = []
# Every position within 300 bases of the origin, every 25th elsewhere.
POSITIONS = sorted(set(list(range(L - 299, L + 1)) + list(range(1, 301)) + list(range(301, L - 299, 25))))
for p in POSITIONS:
    for k, (name, rlen, altf) in enumerate(shapes):
        r = "".join(base(p + i) for i in range(rlen))
        alt = altf(r)
        if alt == r:
            continue
        events.append((p * 100 + k, p, rlen, alt))
events.append((L * 100 + 90, 100, 200, None))  # huge deletion, alt = anchor
# The same three bases inserted between L and 1, anchored on the left (after base L)
# and on the right (before base 1). The inserted text starts with a base that
# differs from base 1 so that VEP's leading-base trim cannot move the anchor.
PAIR_INS = ("A" if base(1) != "A" else "C") + "GT"
events.append((L * 100 + 91, L, 1, base(L) + PAIR_INS))
events.append((L * 100 + 92, 1, 1, PAIR_INS + base(1)))
emit("ok", "CREATE TABLE ev0(event_index UBIGINT, u_pos BIGINT, ref_length BIGINT, alt VARCHAR);")
vals = []
for eid, p, rlen, alt in events:
    if alt is None:
        alt = base(p)
    vals.append(f"({eid}, {p}, {rlen}, {sqlstr(alt)})")
for i in range(0, len(vals), 400):
    emit("ok", "INSERT INTO ev0 VALUES\n  " + ",\n  ".join(vals[i:i + 400]) + ";")

TABLES_SHARED = ["coord_system", "seq_region", "attrib_type", "seq_region_attrib", "gene",
                 "exon_transcript", "translation", "transcript_attrib", "translation_attrib"]
for k in ROTATIONS:
    s = f"core_{k}"
    f = f"funcgen_{k}"
    emit("ok", f"CREATE SCHEMA {s};")
    emit("ok", f"CREATE SCHEMA {f};")
    for t in TABLES_SHARED:
        emit("ok", f"CREATE VIEW {s}.{t} AS SELECT * FROM core.{t};")
    emit("ok", f"""CREATE TABLE {s}.transcript AS SELECT transcript_id, gene_id, seq_region_id,
  rot(seq_region_start, {k}) AS seq_region_start, rot(seq_region_end, {k}) AS seq_region_end,
  seq_region_strand, biotype, is_current, stable_id, version FROM core.transcript;""")
    emit("ok", f"""CREATE TABLE {s}.exon AS SELECT exon_id, seq_region_id,
  rot(seq_region_start, {k}) AS seq_region_start, rot(seq_region_end, {k}) AS seq_region_end,
  seq_region_strand, phase, end_phase, is_current, stable_id FROM core.exon;""")
    emit("ok", f"CREATE VIEW {f}.feature_type AS SELECT * FROM funcgen.feature_type;")
    emit("ok", f"""CREATE TABLE {f}.regulatory_feature AS SELECT regulatory_feature_id, feature_type_id,
  seq_region_id, seq_region_strand, rot(seq_region_start, {k}) AS seq_region_start,
  rot(seq_region_end, {k}) AS seq_region_end, stable_id, regulatory_build_id FROM funcgen.regulatory_feature;""")
    emit("ok", f"""CREATE TABLE {f}.motif_feature AS SELECT motif_feature_id, binding_matrix_id,
  seq_region_id, rot(seq_region_start, {k}) AS seq_region_start, rot(seq_region_end, {k}) AS seq_region_end,
  seq_region_strand, score, stable_id FROM funcgen.motif_feature;""")
    # Reference rotated so base position p sits at rot(p, k).
    emit("ok", f"""CREATE TABLE ref_{k} AS SELECT chrom, "start", "end",
  substr(seq, {L - k} + 1) || substr(seq, 1, {L - k}) AS seq FROM ref;""" if k else
         f"CREATE TABLE ref_{k} AS SELECT * FROM ref;")
    emit("ok", f"""COPY (SELECT line FROM (SELECT 0 AS o, '>' || chrom AS line FROM ref_{k}
  UNION ALL SELECT 1, seq FROM ref_{k}) ORDER BY o) TO '__TEST_DIR__/circular_{k}.fa'
  (FORMAT csv, HEADER false, QUOTE '', ESCAPE '');""")
    emit("ok", f"""COPY (SELECT chrom, {L}, length('>' || chrom) + 1, {L}, {L} + 1 FROM ref_{k})
  TO '__TEST_DIR__/circular_{k}.fa.fai' (FORMAT csv, HEADER false, DELIMITER '\t', QUOTE '', ESCAPE '');""")
    emit("ok", f"CREATE TABLE regions_{k} AS FROM query(duckvep_ensembl_regions_sql('{s}', 'ref_{k}', 'synthetic'));")
    emit("ok", f"CREATE TABLE transcripts_{k} AS FROM query(duckvep_ensembl_transcripts_sql('{s}', 'ref_{k}', 'synthetic'));")
    emit("ok", f"CREATE TABLE features_{k} AS FROM query(duckvep_ensembl_regulation_features_sql('{f}', 'regions_{k}'));")
    emit("query", f"""SELECT loaded FROM duckvep_model_load(
  'circ_{k}',
  'SELECT seq_region, sequence_length, seq_region_name, circular FROM regions_{k} ORDER BY seq_region',
  'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM transcripts_{k} ORDER BY seq_region, transcript_start, transcript_index',
  'SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM transcripts_{k}, LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start',
  interval_feature_query := 'SELECT regulation_feature_index, seq_region, feature_start, feature_end, feature_kind FROM features_{k} ORDER BY regulation_feature_index',
  reference_fasta := '__TEST_DIR__/circular_{k}.fa'
);""", "1", "I")
    emit("ok", f"""CREATE TABLE ev_{k} AS SELECT e.event_index, (SELECT seq_region FROM regions_{k}) AS seq_region,
  rot(e.u_pos, {k})::UBIGINT AS position,
  substr(r.seq || r.seq, rot(e.u_pos, {k}), e.ref_length) AS reference,
  e.alt AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
  NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
  FROM ev0 e, ref_{k} r;""")
    # Stable keys replace frame-specific ordinals; positions never enter.
    emit("ok", f"""CREATE TABLE res_{k} AS
SELECT a.event_index, a.duckvep_event_kind,
  coalesce(t.transcript_stable_id, f.feature_class || ':' || f.stable_id, 'none') AS object,
  a.consequence_mask, a.region_mask, a.impact_code, a.status_code, a.reason_code,
  a.cdna_position, a.cds_position, a.protein_position,
  a.reference_amino_acid_code, a.alternate_amino_acid_code, a.nmd_prediction_code,
  a.nmd_escape_reasons, a.overlap_object_code,
  a.transcript_hgvs, a.protein_hgvs, a.hgvs_shift, a.transcript_hgvs_status,
  a.transcript_hgvs_reason, a.protein_hgvs_status, a.protein_hgvs_reason
FROM query(duckvep_annotate_sql('ev_{k}', 'circ_{k}',
  struct_pack(hgvs := true, upstream_distance := 40, downstream_distance := 40))) a
LEFT JOIN transcripts_{k} t ON t.transcript_index = a.transcript_index
LEFT JOIN features_{k} f ON f.regulation_feature_index = a.regulation_feature_index;""")



OBJ_COLS = """a.consequence_mask, a.region_mask, a.impact_code, a.status_code, a.reason_code,
  a.cdna_position, a.cds_position, a.protein_position,
  a.reference_amino_acid_code, a.alternate_amino_acid_code, a.nmd_prediction_code,
  a.nmd_escape_reasons, a.overlap_object_code,
  a.transcript_hgvs, a.protein_hgvs, a.hgvs_shift, a.transcript_hgvs_status,
  a.transcript_hgvs_reason, a.protein_hgvs_status, a.protein_hgvs_reason"""

# ---------------------------------------------------------------- checks --
EXPECTED = {
    'shape': '28349\t6963\t2015\t1815\t2078\t15\t2200',
    'oracle_counts': 'T_MI\t3808\nT_MID\t93\nT_MW\t3169\nT_NB\t1190\nT_NC\t1773\nT_NF\t1150\nT_PW\t2373\nT_SW\t1251\nT_UW\t3618\nregulatory_region:R_GUARD\t6962\nregulatory_region:R_MID\t21\nregulatory_region:R_NEAR\t168\nregulatory_region:R_WRAP\t421\ntranscription_factor_binding_site:M_WRAP\t151',
    'interbase': '0\tT_MI,T_PW\t11',
    'named_snv': '100\tT_MI\t5_prime_UTR_variant&splice_region_variant\t-\tc.-4C>G\t-\n100\tT_PW\tsplice_donor_variant\t-\tc.-4+1G>C\t-\n100\tT_SW\tmissense_variant\t17\tc.17G>C\tp.Arg6Pro\n2300\tT_SW\tstop_lost\t39\tc.39A>T\tp.Ter13TyrextTer?\n5900\tT_PW\tsplice_acceptor_variant\t-\tc.-3-1G>C\t-\n298500\tT_SW\tstart_lost\t1\tc.1A>T\tp.Met1?\n300000\tT_MI\tsplice_donor_variant\t-\tc.-4+1G>C\t-\n300000\tT_PW\t5_prime_UTR_variant&splice_region_variant\t-\tc.-4C>G\t-\n300000\tT_SW\tmissense_variant\t16\tc.16C>G\tp.Arg6Gly',
    'named_indel': '298007\tT_MW\t5_prime_UTR_variant\t-\tc.-59_-15del\t-\n298007\tT_SW\ttranscript_ablation\t-\tc.1_39del\t-\n299803\tT_PW\t5_prime_UTR_variant&intron_variant&splice_donor_region_variant&splice_donor_variant\t-\tc.-5_-4+3del\t-\n299803\tT_SW\tframeshift_variant\t-\tc.15_19del\tp.Arg6Ter\n299803\tT_UW\t3_prime_UTR_variant\t-\tc.*22_*26del\t-\n299906\tT_PW\t5_prime_UTR_variant&splice_donor_variant\t-\tc.-5_-4+1inv\t-\n299906\tT_SW\tmissense_variant\t-\tc.15_17inv\tp.Arg6Val\n299906\tT_UW\t3_prime_UTR_variant\t-\tc.*22_*24inv\t-\n299906\tregulatory_region:R_WRAP\tregulatory_region_variant\t-\t-\t-\n299906\ttranscription_factor_binding_site:M_WRAP\tTF_binding_site_variant\t-\t-\t-\n300008\tT_SW\tframeshift_variant&stop_gained\t-\tc.17_18insCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCG\tp.Thr7Ter\n300008\tT_UW\t3_prime_UTR_variant\t-\tc.*24_*25insCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCG\t-',
    'flank': 'T_NB\t40\t1\t3000\nT_NF\t40\t1\t3000',
    'default': '10043\t0\t0\t0',
    'proj_shape': '28349\t2332\t2053',
    'mixed': '2474\t1976\t498\t0\t0\t0\t0\t3',
}


def exp(key):
    return EXPECTED.get(key, "@@" + key + "@@")


emit("ok", """CREATE MACRO terms(m) AS (SELECT string_agg(t.consequence, '&' ORDER BY t.consequence)
  FROM duckvep_so_terms() t WHERE (m & t.consequence_mask) <> 0);""")

# Non-vacuity of the frame-0 result the rotations are compared against.
emit("query", """SELECT count(*), count(DISTINCT event_index), count(*) FILTER (WHERE protein_hgvs IS NOT NULL),
  count(*) FILTER (WHERE nmd_prediction_code <> 0), count(*) FILTER (WHERE hgvs_shift > 0),
  count(DISTINCT object), count(*) FILTER (WHERE object = 'none')
FROM res_0;""", exp("shape"), "IIIIIII")

# Rotation equivariance: every rotation is identical to frame 0 in both directions.
for k in ROTATIONS[1:]:
    emit("query", f"""SELECT (SELECT count(*) FROM (SELECT * FROM res_0 EXCEPT SELECT * FROM res_{k})),
  (SELECT count(*) FROM (SELECT * FROM res_{k} EXCEPT SELECT * FROM res_0)),
  (SELECT count(*) FROM res_{k}) - (SELECT count(*) FROM res_0);""", "0\t0\t0", "III")

# A lifted duplicate never produces two rows for one event/object pair.
emit("query", " UNION ALL ".join(
    f"SELECT {k}, count(*) - count(DISTINCT (event_index, object)) FROM res_{k}" for k in ROTATIONS) + " ORDER BY 1;",
     "\n".join(f"{k}\t0" for k in ROTATIONS), "II")

# Independent oracle, one object at a time. Each object is placed in the middle
# of its own ordinary linear (non-circular) model by rotating reference, object and
# events; no object wraps there, so the linear kernel and linear HGVS windows give
# the answer the lifted circular execution must reproduce, for wrapped objects too.
OBJECTS = [(t["name"], t["name"], t["exons"][0][0], t["exons"][-1][1], "tx") for t in TX] + \
          [(f["name"], ("regulatory_region:" if f["kind"] == "reg" else "transcription_factor_binding_site:") + f["name"],
            f["u"][0], f["u"][1], f["kind"]) for f in FEATURES]
oracle_rows = []
for i, (short, key, us, ue, kind) in enumerate(OBJECTS):
    center = pos((us + ue) // 2)
    k = (1500 - center) % L
    c, f = f"lin{i}", f"linf{i}"
    emit("ok", f"CREATE SCHEMA {c};")
    emit("ok", f"CREATE SCHEMA {f};")
    for t in ["coord_system", "seq_region", "attrib_type", "gene", "translation", "transcript_attrib", "translation_attrib"]:
        emit("ok", f"CREATE VIEW {c}.{t} AS SELECT * FROM core.{t};")
    emit("ok", f"CREATE TABLE {c}.seq_region_attrib AS SELECT * FROM core.seq_region_attrib WHERE false;")
    if kind == "tx":
        emit("ok", f"""CREATE TABLE {c}.transcript AS SELECT transcript_id, gene_id, seq_region_id,
  rot(seq_region_start, {k}) AS seq_region_start, rot(seq_region_end, {k}) AS seq_region_end,
  seq_region_strand, biotype, is_current, stable_id, version FROM core.transcript WHERE stable_id = '{short}';""")
        emit("ok", f"CREATE TABLE {c}.exon_transcript AS SELECT * FROM core.exon_transcript WHERE transcript_id IN (SELECT transcript_id FROM {c}.transcript);")
        emit("ok", f"""CREATE TABLE {c}.exon AS SELECT exon_id, seq_region_id, rot(seq_region_start, {k}) AS seq_region_start,
  rot(seq_region_end, {k}) AS seq_region_end, seq_region_strand, phase, end_phase, is_current, stable_id
  FROM core.exon WHERE exon_id IN (SELECT exon_id FROM {c}.exon_transcript);""")
    else:
        # A feature-only model still needs one transcript; a non-coding filler
        # sits at 10..20 and its rows are never compared.
        emit("ok", f"""CREATE TABLE {c}.transcript AS SELECT 900 AS transcript_id, 20 AS gene_id, 10 AS seq_region_id,
  10 AS seq_region_start, 20 AS seq_region_end, 1 AS seq_region_strand, 'lncRNA' AS biotype, 1 AS is_current,
  'FILLER' AS stable_id, 1 AS version;""")
        emit("ok", f"CREATE TABLE {c}.exon_transcript AS SELECT 900 AS exon_id, 900 AS transcript_id, 1 AS rank;")
        emit("ok", f"""CREATE TABLE {c}.exon AS SELECT 900 AS exon_id, 10 AS seq_region_id, 10 AS seq_region_start,
  20 AS seq_region_end, 1 AS seq_region_strand, -1 AS phase, -1 AS end_phase, 1 AS is_current, 'FILLER_E1' AS stable_id;""")
    emit("ok", f"CREATE VIEW {f}.feature_type AS SELECT * FROM funcgen.feature_type;")
    emit("ok", f"""CREATE TABLE {f}.regulatory_feature AS SELECT regulatory_feature_id, feature_type_id, seq_region_id,
  seq_region_strand, rot(seq_region_start, {k}) AS seq_region_start, rot(seq_region_end, {k}) AS seq_region_end,
  stable_id, regulatory_build_id FROM funcgen.regulatory_feature WHERE stable_id = '{short}';""")
    emit("ok", f"""CREATE TABLE {f}.motif_feature AS SELECT motif_feature_id, binding_matrix_id, seq_region_id,
  rot(seq_region_start, {k}) AS seq_region_start, rot(seq_region_end, {k}) AS seq_region_end,
  seq_region_strand, score, stable_id FROM funcgen.motif_feature WHERE stable_id = '{short}';""")
    emit("ok", f"""CREATE TABLE ref_lin{i} AS SELECT chrom, "start", "end",
  substr(seq, {L - k} + 1) || substr(seq, 1, {L - k}) AS seq FROM ref;""" if k else
         f"CREATE TABLE ref_lin{i} AS SELECT * FROM ref;")
    emit("ok", f"""COPY (SELECT line FROM (SELECT 0 AS o, '>' || chrom AS line FROM ref_lin{i}
  UNION ALL SELECT 1, seq FROM ref_lin{i}) ORDER BY o) TO '__TEST_DIR__/linear_{i}.fa'
  (FORMAT csv, HEADER false, QUOTE '', ESCAPE '');""")
    emit("ok", f"""COPY (SELECT chrom, {L}, length('>' || chrom) + 1, {L}, {L} + 1 FROM ref_lin{i})
  TO '__TEST_DIR__/linear_{i}.fa.fai' (FORMAT csv, HEADER false, DELIMITER '\t', QUOTE '', ESCAPE '');""")
    emit("ok", f"CREATE TABLE regions_lin{i} AS FROM query(duckvep_ensembl_regions_sql('{c}', 'ref_lin{i}', 'synthetic'));")
    emit("ok", f"CREATE TABLE transcripts_lin{i} AS FROM query(duckvep_ensembl_transcripts_sql('{c}', 'ref_lin{i}', 'synthetic'));")
    emit("ok", f"CREATE TABLE features_lin{i} AS FROM query(duckvep_ensembl_regulation_features_sql('{f}', 'regions_lin{i}'));")
    emit("query", f"""SELECT loaded FROM duckvep_model_load(
  'lin_{i}',
  'SELECT seq_region, sequence_length, seq_region_name, circular FROM regions_lin{i} ORDER BY seq_region',
  'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM transcripts_lin{i} ORDER BY seq_region, transcript_start, transcript_index',
  'SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM transcripts_lin{i}, LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start',
  interval_feature_query := 'SELECT regulation_feature_index, seq_region, feature_start, feature_end, feature_kind FROM features_lin{i} ORDER BY regulation_feature_index',
  reference_fasta := '__TEST_DIR__/linear_{i}.fa'
);""", "1", "I")
    emit("ok", f"""CREATE TABLE ev_lin{i} AS SELECT e.event_index, (SELECT seq_region FROM regions_lin{i}) AS seq_region,
  rot(e.u_pos, {k})::UBIGINT AS position, substr(r.seq || r.seq, rot(e.u_pos, {k}), e.ref_length) AS reference,
  e.alt AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
  NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
  FROM ev0 e, ref_lin{i} r WHERE rot(e.u_pos, {k}) + e.ref_length - 1 <= {L};""")
    emit("ok", f"""CREATE TABLE res_lin{i} AS
SELECT a.event_index, a.duckvep_event_kind,
  coalesce(t.transcript_stable_id, fe.feature_class || ':' || fe.stable_id, 'none') AS object,
  {OBJ_COLS}
FROM query(duckvep_annotate_sql('ev_lin{i}', 'lin_{i}',
  struct_pack(hgvs := true, upstream_distance := 40, downstream_distance := 40))) a
LEFT JOIN transcripts_lin{i} t ON t.transcript_index = a.transcript_index
LEFT JOIN features_lin{i} fe ON fe.regulation_feature_index = a.regulation_feature_index;""")
    oracle_rows.append((i, key))
sel = ",\n  ".join(
    f"(SELECT count(*) FROM res_0 WHERE object = '{key}' AND event_index IN (SELECT event_index FROM ev_lin{i}))"
    for i, key in oracle_rows)
diffs = " + ".join(
    f"(SELECT count(*) FROM (SELECT * EXCLUDE (duckvep_event_kind) FROM res_0 WHERE object = '{key}' AND event_index IN (SELECT event_index FROM ev_lin{i}) "
    f"EXCEPT SELECT * EXCLUDE (duckvep_event_kind) FROM res_lin{i} WHERE object = '{key}')) + "
    f"(SELECT count(*) FROM (SELECT * EXCLUDE (duckvep_event_kind) FROM res_lin{i} WHERE object = '{key}' "
    f"EXCEPT SELECT * EXCLUDE (duckvep_event_kind) FROM res_0 WHERE object = '{key}' AND event_index IN (SELECT event_index FROM ev_lin{i})))"
    for i, key in oracle_rows)
emit("query", "SELECT " + diffs + ";", "0", "I")
emit("query", "SELECT object, count(*) FROM (" + " UNION ALL ".join(
    f"SELECT object FROM res_lin{i} WHERE object = '{key}'" for i, key in oracle_rows) + ") GROUP BY object ORDER BY object;",
     exp("oracle_counts"), "TI")

# The two insertion interbase orientations at the origin, one anchored on the last
# base and one on the first, insert the same bases into the same gap. Their HGVS
# and shift agree; only VEP's anchor-side region for the exon-boundary transcripts
# may differ, exactly as the linear oracle above reproduces mid-sequence.
emit("query", """SELECT count(*) FILTER (WHERE a.transcript_hgvs IS DISTINCT FROM b.transcript_hgvs
    OR a.hgvs_shift IS DISTINCT FROM b.hgvs_shift OR a.protein_hgvs IS DISTINCT FROM b.protein_hgvs),
  string_agg(object, ',' ORDER BY object) FILTER (WHERE a.region_mask <> b.region_mask), count(*)
FROM (SELECT * FROM res_0 WHERE event_index = 300091) a
JOIN (SELECT * FROM res_0 WHERE event_index = 300092) b USING (object);""", exp("interbase"), "ITI")

# Named origin cases: SNVs on the last and first base, the exon-intron junction
# of both strands at the origin, CDS start and end of the wrapped single-exon CDS,
# an MNV, deletions and a long insertion across L/1.
named = [(100, "T_PW"), (100, "T_MI"), (300000, "T_PW"), (300000, "T_MI"), (300000, "T_SW"),
         (5900, "T_PW"), (298500, "T_SW"), (2300, "T_SW"), (100, "T_SW")]
emit("query", "SELECT event_index, object, terms(consequence_mask), coalesce(cds_position::VARCHAR, '-'), coalesce(transcript_hgvs, '-'), coalesce(protein_hgvs, '-') FROM res_0 WHERE (event_index, object) IN (VALUES " +
     ", ".join(f"({e}, '{o}')" for e, o in named) + ") ORDER BY 1, 2;", exp("named_snv"), "TTTTTT")
named2 = [(299906, "T_SW"), (299906, "T_PW"), (299906, "T_UW"), (299803, "T_SW"), (299803, "T_PW"), (299803, "T_UW"),
          (298007, "T_SW"), (298007, "T_MW"), (300008, "T_SW"), (300008, "T_UW"), (299906, "regulatory_region:R_WRAP"),
          (299906, "transcription_factor_binding_site:M_WRAP")]
emit("query", "SELECT event_index, object, terms(consequence_mask), coalesce(cds_position::VARCHAR, '-'), coalesce(transcript_hgvs, '-'), coalesce(protein_hgvs, '-') FROM res_0 WHERE (event_index, object) IN (VALUES " +
     ", ".join(f"({e}, '{o}')" for e, o in named2) + ") ORDER BY 1, 2;", exp("named_indel"), "TTTTTT")

# Flanks cross the origin for objects that do not: T_NB (plus, starts at 26) is
# upstream of the last 15 bases before the origin and of bases 1..25, and T_NF
# (minus, ends at 2980) is upstream of 2981..3000 and 1..20.
emit("query", """WITH r AS (SELECT object, event_index, terms(consequence_mask) AS t FROM res_0
  WHERE object IN ('T_NB', 'T_NF') AND event_index % 100 = 0)
SELECT object, count(*), min(event_index // 100), max(event_index // 100)
FROM r WHERE t = 'upstream_gene_variant' GROUP BY 1 ORDER BY 1;""", exp("flank"), "TIII")

# Default 5000 bp flank windows on a 900 bp circle see every object at several laps
# and keep the nearest image; the result is still rotation-invariant.
for k in (0, 611):
    emit("ok", f"""CREATE TABLE res_default_{k} AS
SELECT a.event_index, coalesce(t.transcript_stable_id, f.feature_class || ':' || f.stable_id, 'none') AS object,
  a.consequence_mask, a.region_mask, a.impact_code, a.cdna_position, a.cds_position, a.protein_position,
  a.transcript_hgvs, a.protein_hgvs, a.nmd_prediction_code
FROM query(duckvep_annotate_sql('ev_{k}', 'circ_{k}', struct_pack(hgvs := true))) a
LEFT JOIN transcripts_{k} t ON t.transcript_index = a.transcript_index
LEFT JOIN features_{k} f ON f.regulation_feature_index = a.regulation_feature_index
WHERE a.event_index % 7 = 0;""")
emit("query", """SELECT (SELECT count(*) FROM res_default_0),
  (SELECT count(*) FROM (SELECT * FROM res_default_0 EXCEPT SELECT * FROM res_default_611)),
  (SELECT count(*) FROM (SELECT * FROM res_default_611 EXCEPT SELECT * FROM res_default_0)),
  (SELECT count(*) - count(DISTINCT (event_index, object)) FROM res_default_611);""", exp("default"), "IIII")

# Worker splits: the same events on four threads and 9001 events across several
# DuckDB vectors give the frame-311 result again.
emit("ok", "SET threads = 4;")
emit("ok", """CREATE TABLE res_parallel AS
SELECT a.event_index, a.duckvep_event_kind,
  coalesce(t.transcript_stable_id, f.feature_class || ':' || f.stable_id, 'none') AS object,
  a.consequence_mask, a.region_mask, a.impact_code, a.status_code, a.reason_code,
  a.cdna_position, a.cds_position, a.protein_position,
  a.reference_amino_acid_code, a.alternate_amino_acid_code, a.nmd_prediction_code,
  a.nmd_escape_reasons, a.overlap_object_code,
  a.transcript_hgvs, a.protein_hgvs, a.hgvs_shift, a.transcript_hgvs_status,
  a.transcript_hgvs_reason, a.protein_hgvs_status, a.protein_hgvs_reason
FROM query(duckvep_annotate_sql('ev_611', 'circ_611',
  struct_pack(hgvs := true, upstream_distance := 40, downstream_distance := 40))) a
LEFT JOIN transcripts_611 t ON t.transcript_index = a.transcript_index
LEFT JOIN features_611 f ON f.regulation_feature_index = a.regulation_feature_index;""")
emit("query", """SELECT (SELECT count(*) FROM (SELECT * FROM res_parallel EXCEPT SELECT * FROM res_611)),
  (SELECT count(*) FROM (SELECT * FROM res_611 EXCEPT SELECT * FROM res_parallel)),
  (SELECT count(*) FROM res_parallel) = (SELECT count(*) FROM res_0);""", "0\t0\ttrue", "IIT")
emit("ok", "SET threads = 1;")

# Contract limits.
emit("error", f"""SELECT count(*) FROM query(duckvep_annotate_sql(
  (SELECT 'ev_long' WHERE false), 'circ_0'));""", "") if False else None
emit("ok", f"""CREATE TABLE ev_long AS SELECT 1::UBIGINT AS event_index, (SELECT seq_region FROM regions_0) AS seq_region,
  2::UBIGINT AS position, repeat('A', {L + 1})::VARCHAR AS reference, 'C'::VARCHAR AS alternate,
  NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
  NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position;""")
emit("error", "SELECT count(*) FROM query(duckvep_annotate_sql('ev_long', 'circ_0'));",
     "variant span exceeds sequence-region length")
emit("ok", f"""CREATE TABLE ev_past AS SELECT 1::UBIGINT AS event_index, (SELECT seq_region FROM regions_0) AS seq_region,
  {L + 1}::UBIGINT AS position, 'A'::VARCHAR AS reference, 'C'::VARCHAR AS alternate,
  NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
  NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position;""")
emit("error", "SELECT count(*) FROM query(duckvep_annotate_sql('ev_past', 'circ_0'));",
     "variant span exceeds sequence-region length")
emit("ok", f"""CREATE TABLE ev_sv AS SELECT 1::UBIGINT AS event_index, (SELECT seq_region FROM regions_0) AS seq_region,
  100::UBIGINT AS position, 'N'::VARCHAR AS reference, '<DEL>'::VARCHAR AS alternate,
  200::UBIGINT AS end_position, 'DEL'::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
  NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position;""")
emit("error", "SELECT count(*) FROM query(duckvep_annotate_sql('ev_sv', 'circ_0'));",
     "structural and breakend annotation is not supported for models with wrapped circular objects")
emit("error", "SELECT count(*) FROM duckvep_haplotypes('SELECT 1', 'circ_0');",
     "phased edit sets are not supported for models with wrapped circular objects")


# Projected edits (typed transcript presentation: alleles, positions, exon/intron
# ordinals, codons, amino acids) are rotation-invariant as well.
PROJ_COLS = """a.event_index, coalesce(t.transcript_stable_id, 'none') AS object, a.output_allele, a.interbase,
  a.cdna_start, a.cdna_end, a.cds_start, a.cds_end, a.protein_start, a.protein_end,
  a.exon_first, a.exon_last, a.exon_total, a.intron_first, a.intron_last, a.intron_total,
  a.transcript_distance, a.cds_start_nf, a.cds_end_nf, a.reference_amino_acids,
  a.alternate_amino_acids, a.reference_codons, a.alternate_codons, a.consequence_mask, a.region_mask"""
for k in ROTATIONS:
    emit("ok", f"""CREATE TABLE proj_{k} AS
SELECT {PROJ_COLS}
FROM query(duckvep_annotate_projected_sql('ev_{k}', 'circ_{k}',
  struct_pack(upstream_distance := 40, downstream_distance := 40))) a
LEFT JOIN transcripts_{k} t ON t.transcript_index = a.transcript_index;""")
emit("query", "SELECT (SELECT count(*) FROM proj_0), (SELECT count(*) FILTER (WHERE cds_start IS NOT NULL) FROM proj_0), " +
     "(SELECT count(*) FILTER (WHERE reference_codons IS NOT NULL) FROM proj_0);", exp("proj_shape"), "III")
for k in ROTATIONS[1:]:
    emit("query", f"""SELECT (SELECT count(*) FROM (SELECT * FROM proj_0 EXCEPT SELECT * FROM proj_{k})),
  (SELECT count(*) FROM (SELECT * FROM proj_{k} EXCEPT SELECT * FROM proj_0));""", "0\t0", "II")

# ---------------------------------------------------- mixed-region model --
# A lifted circular region and an ordinary linear region in one model: the linear
# region must annotate exactly as it does in a model of its own, the circular one
# as it does alone, with events from both interleaved in the same DuckDB vectors.
mx = random.Random(4242)
CM, LM = 400, 600
mix_ref = {"CM": [mx.choice("ACGT") for _ in range(CM)], "LM": [mx.choice("ACGT") for _ in range(LM)]}


def design_cds(contig, positions, strand):
    codons = (len(positions)) // 3
    seq = "ATG" + "".join(mx.choice(NONSTOP) for _ in range(codons - 2)) + "TAA"
    for u, ch in zip(positions, seq):
        mix_ref[contig][(u - 1) % len(mix_ref[contig])] = ch if strand > 0 else COMP[ch]


design_cds("CM", list(range(381, 381 + 39)), 1)   # wrapped single exon, CDS = whole exon
design_cds("LM", list(range(103, 103 + 54)), 1)   # exon 100..160, CDS cDNA 4..57
mix_seq = {k: "".join(v) for k, v in mix_ref.items()}
emit("ok", "CREATE SCHEMA core_mix;")
emit("ok", "CREATE VIEW core_mix.coord_system AS SELECT * FROM core.coord_system;")
emit("ok", "CREATE VIEW core_mix.attrib_type AS SELECT * FROM core.attrib_type;")
emit("ok", "CREATE VIEW core_mix.transcript_attrib AS SELECT * FROM core.transcript_attrib;")
emit("ok", "CREATE VIEW core_mix.translation_attrib AS SELECT * FROM core.translation_attrib;")
emit("ok", f"""CREATE TABLE core_mix.seq_region AS SELECT * FROM (VALUES (20, 'CM', 1, {CM}), (21, 'LM', 1, {LM}))
  t(seq_region_id, name, coord_system_id, length);""")
emit("ok", "CREATE TABLE core_mix.seq_region_attrib AS SELECT * FROM (VALUES (20, 1, 1)) t(seq_region_id, attrib_type_id, value);")
emit("ok", "CREATE TABLE core_mix.gene AS SELECT * FROM (VALUES (60, 'protein_coding', 1, 'G_CM', 1), (61, 'protein_coding', 1, 'G_LM', 1)) t(gene_id, biotype, is_current, stable_id, version);")
emit("ok", f"""CREATE TABLE core_mix.transcript AS SELECT * FROM (VALUES
  (70, 60, 20, 381, 20, 1, 'protein_coding', 1, 'TM_CM', 1),
  (71, 61, 21, 100, 160, 1, 'protein_coding', 1, 'TM_LM', 1))
  t(transcript_id, gene_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, biotype, is_current, stable_id, version);""")
emit("ok", """CREATE TABLE core_mix.exon AS SELECT * FROM (VALUES
  (700, 20, 381, 20, 1, 0, -1, 1, 'TM_CM_E1'), (701, 21, 100, 160, 1, -1, -1, 1, 'TM_LM_E1'))
  t(exon_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, phase, end_phase, is_current, stable_id);""")
emit("ok", "CREATE TABLE core_mix.exon_transcript AS SELECT * FROM (VALUES (700, 70, 1), (701, 71, 1)) t(exon_id, transcript_id, rank);")
emit("ok", """CREATE TABLE core_mix.translation AS SELECT * FROM (VALUES
  (80, 70, 1, 700, 39, 700, 'P_CM', 1), (81, 71, 4, 701, 57, 701, 'P_LM', 1))
  t(translation_id, transcript_id, seq_start, start_exon_id, seq_end, end_exon_id, stable_id, version);""")
emit("ok", f"""CREATE TABLE ref_mix AS SELECT * FROM (VALUES
  ('CM', 0, {CM}, {sqlstr(mix_seq['CM'])}), ('LM', 0, {LM}, {sqlstr(mix_seq['LM'])}))
  t(chrom, "start", "end", seq);""")
emit("ok", f"""COPY (SELECT line FROM (SELECT 0 AS o, '>CM' AS line UNION ALL SELECT 1, seq FROM ref_mix WHERE chrom = 'CM'
  UNION ALL SELECT 2, '>LM' UNION ALL SELECT 3, seq FROM ref_mix WHERE chrom = 'LM') ORDER BY o)
  TO '__TEST_DIR__/circular_mixed.fa' (FORMAT csv, HEADER false, QUOTE '', ESCAPE '');""")
emit("ok", f"""COPY (SELECT * FROM (VALUES ('CM', {CM}, 4, {CM}, {CM + 1}), ('LM', {LM}, {4 + CM + 1 + 4}, {LM}, {LM + 1})))
  TO '__TEST_DIR__/circular_mixed.fa.fai' (FORMAT csv, HEADER false, DELIMITER '\\t', QUOTE '', ESCAPE '');""")
emit("ok", "CREATE TABLE regions_mix AS FROM query(duckvep_ensembl_regions_sql('core_mix', 'ref_mix', 'synthetic'));")
emit("ok", "CREATE TABLE transcripts_mix AS FROM query(duckvep_ensembl_transcripts_sql('core_mix', 'ref_mix', 'synthetic'));")
MIX_MODELS = (("mix", "true"), ("mix_c", "seq_region_name = 'CM'"), ("mix_l", "seq_region_name = 'LM'"))
for name, where in MIX_MODELS:
    if name != "mix":
        emit("ok", f"""CREATE TABLE transcripts_{name} AS SELECT * REPLACE ((row_number() OVER (ORDER BY seq_region, transcript_start) - 1)::UINTEGER AS transcript_index)
  FROM transcripts_mix WHERE {where};""")
    else:
        pass
    t = "transcripts_mix" if name == "mix" else f"transcripts_{name}"
    where_q = where.replace("'", "''")
    emit("query", f"""SELECT loaded FROM duckvep_model_load(
  '{name}',
  'SELECT seq_region, sequence_length, seq_region_name, circular FROM regions_mix WHERE {where_q} ORDER BY seq_region',
  'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM {t} ORDER BY seq_region, transcript_start, transcript_index',
  'SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM {t}, LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start',
  reference_fasta := '__TEST_DIR__/circular_mixed.fa'
);""", "1", "I")
mix_events = []
for contig, rng_ in (("CM", range(1, CM + 1)), ("LM", range(80, 181))):
    seq = mix_seq[contig]
    for p in rng_:
        for k, (rlen, altf) in enumerate([(1, lambda r: COMP[r[0]]), (3, lambda r: r[0]), (1, lambda r: r[0] + "GTC"),
                                          (1, lambda r: "TGC" + r[0]), (4, lambda r: revcomp(r))]):
            r = "".join(seq[(p - 1 + i) % len(seq)] for i in range(rlen))
            if contig == "LM" and p + rlen - 1 > LM:
                continue
            alt = altf(r)
            if alt != r:
                mix_events.append((p * 10 + k + (0 if contig == "CM" else 10_000_000), contig, p, r, alt))
emit("ok", "CREATE TABLE ev_mix_raw(event_index UBIGINT, contig VARCHAR, position UBIGINT, reference VARCHAR, alternate VARCHAR);")
mv = [f"({e}, '{c}', {p}, {sqlstr(r)}, {sqlstr(a)})" for e, c, p, r, a in mix_events]
for i in range(0, len(mv), 500):
    emit("ok", "INSERT INTO ev_mix_raw VALUES\n  " + ",\n  ".join(mv[i:i + 500]) + ";")
emit("ok", """CREATE TABLE ev_mix AS SELECT r.event_index, g.seq_region, r.position, r.reference, r.alternate,
  NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
  NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
  FROM ev_mix_raw r JOIN regions_mix g ON g.seq_region_name = r.contig ORDER BY r.position, r.event_index;""")
for name, where in MIX_MODELS:
    t = "transcripts_mix" if name == "mix" else f"transcripts_{name}"
    if name != "mix":
        emit("ok", f"CREATE TABLE ev_{name} AS SELECT * FROM ev_mix WHERE seq_region IN (SELECT seq_region FROM regions_mix WHERE {where});")
    emit("ok", f"""CREATE TABLE res_{name} AS
SELECT a.event_index, a.duckvep_event_kind, coalesce(t.transcript_stable_id, 'none') AS object,
  {OBJ_COLS}
FROM query(duckvep_annotate_sql('ev_{name}', '{name}', struct_pack(hgvs := true, upstream_distance := 40, downstream_distance := 40))) a
LEFT JOIN {t} t ON t.transcript_index = a.transcript_index;""")
emit("query", """SELECT (SELECT count(*) FROM res_mix), (SELECT count(*) FROM res_mix_c), (SELECT count(*) FROM res_mix_l),
  (SELECT count(*) FROM (SELECT * FROM res_mix WHERE event_index < 10000000 EXCEPT SELECT * FROM res_mix_c)),
  (SELECT count(*) FROM (SELECT * FROM res_mix_c EXCEPT SELECT * FROM res_mix WHERE event_index < 10000000)),
  (SELECT count(*) FROM (SELECT * FROM res_mix WHERE event_index >= 10000000 EXCEPT SELECT * FROM res_mix_l)),
  (SELECT count(*) FROM (SELECT * FROM res_mix_l EXCEPT SELECT * FROM res_mix WHERE event_index >= 10000000)),
  (SELECT count(DISTINCT object) FROM res_mix);""", exp("mixed"), "IIIIIIII")

open("/dev/stderr", "w").write("events=%d\n" % len(events))
sys.stdout.write("\n".join(lines).rstrip("\n") + "\n")
