#!/usr/bin/env python3
"""Differential fuzz of compound haplotype composition and HGVS against independent oracles.

Each case is a random single-exon transcript (5' UTR, a complete CDS, 3' UTR) on either strand and two to four
literal cis records inside its CDS: substitutions, deletions, insertions and unequal delins. Short CDSs make
overlapping records common. Three checks per case:

  composition  an independent literal applier (below) decides conflict versus joint sequence; DuckVEP must report
               edit_conflict exactly when the applier finds one, and otherwise the same edited CDS.
  DNA HGVS     Mutalyzer normalizes the applier's joint allele; DuckVEP's hgvsc must name the same allele.
  protein      Mutalyzer's predicted protein (read through the 3' UTR when the stop is lost) must equal DuckVEP's,
               and the two protein HGVS strings are compared after known notation differences.

The literal rule: every record applies at its own coordinates. A reference base is kept, substituted or removed;
an interbase site receives at most one inserted sequence. Two records that claim one base or site differently
conflict. An opaque delins (unequal nonzero REF/ALT) is a unit: an identical copy composes, an insertion at its
first site precedes it, any other record that changes or removes one of its bases conflicts. Unchanged bases
inside a substitution are context, not a claim.

Usage:
  compound_hgvs_fuzz.py --extension build/release/duckvep.duckdb_extension \
      --mutalyzer-python /path/to/venv/bin/python --cases 2000 --seed 1 --work DIR
Writes DIR/cases.jsonl (every case with all three verdicts) and prints a summary; exits 1 on any mismatch
that is not a documented notation difference.
"""
import argparse, json, os, random, subprocess, sys, tempfile

import duckdb

STOPS = {"TAA", "TAG", "TGA"}
BASES = "ACGT"
CODE = {}
for i, a in enumerate("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    CODE["TCAG"[i // 16] + "TCAG"[i // 4 % 4] + "TCAG"[i % 4]] = a
COMP = str.maketrans("ACGT", "TGCA")
U5, U3 = 20, 30


def rc(s):
    return s.translate(COMP)[::-1]


def translate(seq):
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODE[seq[i:i + 3]]
        out.append(aa)
        if aa == "*":
            break
    return "".join(out)


def random_case(rng, case_id):
    n = rng.randint(4, 12)
    codons = ["ATG"]
    while len(codons) < n - 1:
        c = "".join(rng.choice(BASES) for _ in range(3))
        if c not in STOPS:
            codons.append(c)
    codons.append(rng.choice(sorted(STOPS)))
    cds = "".join(codons)
    u5 = "".join(rng.choice(BASES) for _ in range(U5))
    u3 = "".join(rng.choice(BASES) for _ in range(U3))
    tx = u5 + cds + u3
    edits = []
    for _ in range(rng.randint(2, 4)):
        kind = rng.choice(["sub", "sub", "del", "ins", "delins"])
        if kind == "ins":
            gap = rng.randint(U5 + 1, U5 + len(cds) - 1)  # before transcript base `gap` (0-based), inside the CDS
            alt = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 6)))
            edits.append((gap, "", alt))
            continue
        length = rng.randint(1, 3) if kind == "sub" else rng.randint(1, 6) if kind == "del" else rng.randint(1, 4)
        start = rng.randint(U5, U5 + len(cds) - length)
        ref = tx[start:start + length]
        if kind == "del":
            edits.append((start, ref, ""))
        elif kind == "sub":
            alt = [rng.choice([b for b in BASES if b != r]) for r in ref]
            for j in range(1, len(ref) - 1):  # interior bases may stay unchanged
                if rng.random() < 0.3:
                    alt[j] = ref[j]
            edits.append((start, ref, "".join(alt)))
        else:
            alt_len = rng.choice([k for k in range(1, 6) if k != length])
            alt = [rng.choice(BASES) for _ in range(alt_len)]
            alt[0] = rng.choice([b for b in BASES if b != ref[0]])
            alt[-1] = rng.choice([b for b in BASES if b != ref[-1] and (alt_len > 1 or b != ref[0])])
            if alt_len == 1 and alt[0] == ref[-1]:
                alt[0] = rng.choice([b for b in BASES if b not in (ref[0], ref[-1])])
            edits.append((start, ref, "".join(alt)))
    return {"case_id": case_id, "strand": rng.choice([1, -1]), "u5": u5, "cds": cds, "u3": u3, "edits": edits}


def apply_literal(tx, edits):
    """Independent literal applier on transcript-orientation islands (start0, ref, alt). Returns (conflict, seq).

    Claims are changed bases (substituted or removed), insertion sites and opaque spans. Unchanged bases inside
    a substitution are context and claim nothing."""
    sub, dele, gaps, opaque = {}, set(), {}, []
    inside_opaque = lambda at: any(s <= at < e for s, e, _ in opaque)
    for start, ref, alt in edits:
        end = start + len(ref)
        if ref and alt and len(ref) != len(alt):
            if (start, end, alt) in opaque:
                continue
            if any(s < end and start < e for s, e, _ in opaque) or any(start <= at < end for at in list(sub) + list(dele)) \
                    or any(start < g < end for g in gaps):
                return True, None
            opaque.append((start, end, alt))
            continue
        if not ref:
            if gaps.get(start, alt) != alt or any(s < start < e for s, e, _ in opaque):
                return True, None
            gaps[start] = alt
            continue
        for j in range(len(ref)):
            at = start + j
            if not alt:
                if at in sub or inside_opaque(at):
                    return True, None
                dele.add(at)
            elif alt[j] != ref[j]:
                if at in dele or sub.get(at, alt[j]) != alt[j] or inside_opaque(at):
                    return True, None
                sub[at] = alt[j]
    out, i = [], 0
    starts = {s: (e, a) for s, e, a in opaque}
    while i <= len(tx):
        out.append(gaps.get(i, ""))
        if i == len(tx):
            break
        if i in starts:
            e, a = starts[i]
            out.append(a)
            i = e
            continue
        if i not in dele:
            out.append(sub.get(i, tx[i]))
        i += 1
    return False, "".join(out)


def vcf_record(case, start, ref, alt):
    """Genomic VCF record for a transcript-orientation island: 1-based POS, REF, ALT on the + strand genome."""
    tx = case["u5"] + case["cds"] + case["u3"]
    total = len(tx)
    genome = tx if case["strand"] == 1 else rc(tx)
    if case["strand"] == 1:
        g0, gref, galt = start, ref, alt
    else:
        g0, gref, galt = total - (start + len(ref)), rc(ref), rc(alt)
    if gref and galt:
        return 1000 + g0 + 1, gref, galt
    anchor = g0 - 1  # pure deletion or insertion: the base before, genomic orientation
    return 1000 + anchor + 1, genome[anchor] + gref, genome[anchor] + galt


def run_duckvep(extension, cases):
    con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    con.execute(f"LOAD '{extension}'")
    tx_rows, ex_rows, calls = [], [], []
    for i, c in enumerate(cases):
        total = len(c["u5"]) + len(c["cds"]) + len(c["u3"])
        cds_lo = 1001 + (U5 if c["strand"] == 1 else U3)
        tx_rows.append((i, i, 1001, 1000 + total, c["strand"], i, 3, cds_lo, cds_lo + len(c["cds"]) - 1,
                        c["cds"], 1, c["u5"], c["u3"]))
        ex_rows.append((i, 1001, 1000 + total, 1, total, 0, 0))
        for k, (start, ref, alt) in enumerate(c["edits"]):
            pos, r, a = vcf_record(c, start, ref, alt)
            calls.append((i * 10 + k + 1, i, pos, r, a, 1, i, i))
    con.execute("CREATE TABLE tx(transcript_index UINTEGER, seq_region UINTEGER, transcript_start UBIGINT, transcript_end UBIGINT,"
                " strand TINYINT, gene_index UINTEGER, transcript_flags UBIGINT, cds_start UBIGINT, cds_end UBIGINT,"
                " cds_sequence VARCHAR, codon_table UTINYINT, pre_cds_sequence VARCHAR, post_cds_sequence VARCHAR)")
    con.executemany("INSERT INTO tx VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", tx_rows)
    con.execute("CREATE TABLE ex(transcript_index UINTEGER, exon_start UBIGINT, exon_end UBIGINT, exon_cdna_start UBIGINT,"
                " exon_cdna_end UBIGINT, phase TINYINT, end_phase TINYINT)")
    con.executemany("INSERT INTO ex VALUES (?,?,?,?,?,?,?)", ex_rows)
    con.execute("CREATE TABLE calls(event_index BIGINT, seq_region INTEGER, position BIGINT, reference VARCHAR, alternate VARCHAR,"
                " alt_index INTEGER, transcript_index INTEGER, sample_index INTEGER)")
    con.executemany("INSERT INTO calls VALUES (?,?,?,?,?,?,?,?)", calls)
    con.execute("SELECT loaded FROM duckvep_model_load('fuzz', 'SELECT DISTINCT seq_region FROM tx ORDER BY 1',"
                " 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags,"
                " cds_start, cds_end, cds_sequence::BLOB cds_sequence, codon_table, pre_cds_sequence::BLOB pre_cds_sequence,"
                " post_cds_sequence::BLOB post_cds_sequence FROM tx ORDER BY transcript_index',"
                " 'SELECT * FROM ex ORDER BY transcript_index')").fetchall()
    rows = con.execute("SELECT carriers[1].sample_index, prediction_status, prediction_reason, cds, protein, hgvsc, hgvsp,"
                       " haplotype_consequences, hgvsc_status, hgvsp_status FROM duckvep_haplotypes('SELECT *, [1, 0] alleles, [false, true] phase_before,"
                       " 7::BIGINT phase_set FROM calls ORDER BY seq_region, position, event_index', 'fuzz',"
                       " phase_policy := 'strict', hgvs := true)").fetchall()
    out = {}
    for sample, status, reason, cds, protein, hgvsc, hgvsp, terms, hgvsc_status, hgvsp_status in rows:
        out[sample] = dict(status=status, reason=reason, cds=cds, protein=protein, hgvsc=hgvsc, hgvsp=hgvsp, terms=terms,
                           hgvsc_status=hgvsc_status, hgvsp_status=hgvsp_status)
    return out


MUTALYZER = r'''
import json, sys
from mutalyzer.normalizer import normalize
for line in sys.stdin:
    q = json.loads(line)
    try:
        r = normalize(q["description"])
        p = r.get("protein") or {}
        print(json.dumps({"id": q["id"], "dna": r.get("normalized_description"), "p": p.get("description"),
                          "predicted": p.get("predicted"), "errors": [e.get("code") for e in r.get("errors", [])]}))
    except Exception as e:
        print(json.dumps({"id": q["id"], "errors": ["EXC:" + type(e).__name__]}))
'''


def annotation(name, total, cds_start0, cds_end0):
    def rng(s, e, strand=True):
        loc = {"type": "range", "start": {"type": "point", "position": s}, "end": {"type": "point", "position": e}}
        if strand:
            loc["strand"] = 1
        return loc
    return {"id": name, "type": "record", "location": rng(0, total, False), "features": [
        {"type": "gene", "location": rng(0, total), "id": "G_" + name, "qualifiers": {"name": "G_" + name}, "features": [
            {"type": "mRNA", "location": rng(0, total), "id": "T_" + name, "features": [
                {"type": "exon", "location": rng(0, total), "id": "T_" + name + ".e1"},
                {"type": "CDS", "location": rng(cds_start0, cds_end0), "id": "P_T_" + name}]}]}]}


def joint_description(tx, seq, cds_start0, cds_end0):
    """Applier's joint allele as one delins over the differing span, in c. coordinates of the reference CDS."""
    p = 0
    while p < len(tx) and p < len(seq) and tx[p] == seq[p]:
        p += 1
    s = 0
    while s < len(tx) - p and s < len(seq) - p and tx[len(tx) - 1 - s] == seq[len(seq) - 1 - s]:
        s += 1
    ref_span, alt = (p, len(tx) - s), seq[p:len(seq) - s]

    def c(i):  # 0-based transcript index -> c. position text
        if i >= cds_end0:
            return f"*{i - cds_end0 + 1}"
        return str(i - cds_start0 + 1) if i >= cds_start0 else str(i - cds_start0)
    if ref_span[0] == ref_span[1]:
        return f"c.{c(p - 1)}_{c(p)}ins{alt}"
    return f"c.{c(ref_span[0])}_{c(ref_span[1] - 1)}delins{alt}" if alt else f"c.{c(ref_span[0])}_{c(ref_span[1] - 1)}del"


def run_mutalyzer(python, cases, joint, work):
    cache = os.path.join(work, "mutalyzer-cache")
    os.makedirs(cache, exist_ok=True)
    queries = []
    for i, c in enumerate(cases):
        if joint[i] is None:
            continue
        name = f"f{i}"
        tx = c["u5"] + c["cds"] + c["u3"]
        json.dump(annotation(name, len(tx), U5, U5 + len(c["cds"])), open(os.path.join(cache, name + ".annotations"), "w"))
        open(os.path.join(cache, name + ".sequence"), "w").write(tx)
        queries.append({"id": i, "description": f"{name}(T_{name}):" + joint_description(tx, joint[i], U5, U5 + len(c["cds"]))})
    settings = os.path.join(work, "mutalyzer.settings")
    open(settings, "w").write(f"MUTALYZER_CACHE_DIR = {cache}\nMUTALYZER_FILE_CACHE_ADD = false\n")
    proc = subprocess.run([python, "-c", MUTALYZER], input="\n".join(json.dumps(q) for q in queries) + "\n",
                          capture_output=True, text=True, env=dict(os.environ, MUTALYZER_SETTINGS=settings))
    if proc.returncode:
        sys.exit("mutalyzer failed: " + proc.stderr[-2000:])
    out = {}
    for line in proc.stdout.splitlines():
        r = json.loads(line)
        out[r["id"]] = r
    return out, {q["id"]: q["description"] for q in queries}


def run_vep(prefix, cases, joint, work):
    """Executable VEP on one merged record per composed case: the joint allele over its differing span."""
    d = os.path.join(work, "vep")
    os.makedirs(d, exist_ok=True)
    fasta, gff, vcf = (os.path.join(d, n) for n in ("ref.fa", "model.gff3", "merged.vcf"))
    records, gff_rows = [], []
    with open(fasta, "w") as fa:
        for i, c in enumerate(cases):
            tx = c["u5"] + c["cds"] + c["u3"]
            genome = tx if c["strand"] == 1 else rc(tx)
            fa.write(f">f{i}\n{'N' * 1000}{genome}{'N' * 100}\n")
            total, lo = len(tx), 1001 + (U5 if c["strand"] == 1 else U3)
            st = "+" if c["strand"] == 1 else "-"
            gff_rows += [f"f{i}\tfuzz\tgene\t1001\t{1000 + total}\t.\t{st}\t.\tID=gene:G{i};biotype=protein_coding",
                         f"f{i}\tfuzz\tmRNA\t1001\t{1000 + total}\t.\t{st}\t.\tID=transcript:T{i};Parent=gene:G{i};biotype=protein_coding",
                         f"f{i}\tfuzz\texon\t1001\t{1000 + total}\t.\t{st}\t.\tParent=transcript:T{i};rank=1",
                         f"f{i}\tfuzz\tCDS\t{lo}\t{lo + len(c['cds']) - 1}\t.\t{st}\t0\tID=CDS:P{i};Parent=transcript:T{i};protein_id=P{i}"]
            if joint[i] is None:
                continue
            seq = joint[i]
            p = 0
            while p < len(tx) and p < len(seq) and tx[p] == seq[p]:
                p += 1
            q = 0
            while q < len(tx) - p and q < len(seq) - p and tx[len(tx) - 1 - q] == seq[len(seq) - 1 - q]:
                q += 1
            if p == len(tx) and p == len(seq):
                continue
            ref_t, alt_t = tx[p:len(tx) - q], seq[p:len(seq) - q]
            g0 = p if c["strand"] == 1 else total - (len(tx) - q)
            ref_g, alt_g = (ref_t, alt_t) if c["strand"] == 1 else (rc(ref_t), rc(alt_t))
            pos = 1000 + g0 + 1
            if not ref_g or not alt_g:
                anchor = genome[g0 - 1]
                pos, ref_g, alt_g = pos - 1, anchor + ref_g, anchor + alt_g
            records.append((f"f{i}", pos, i, ref_g, alt_g))
    with open(gff, "w") as g:
        g.write("##gff-version 3\n" + "\n".join(sorted(gff_rows, key=lambda r: (r.split("\t")[0], int(r.split("\t")[3])))) + "\n")
    with open(vcf, "w") as v:
        v.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID=f{i}>\n" for i in range(len(cases))))
        v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, i, r, a in sorted(records, key=lambda x: (x[0], x[1])):
            v.write(f"{chrom}\t{pos}\tcase{i}\t{r}\t{a}\t.\tPASS\t.\n")
    binp = os.path.join(prefix, "bin")
    env = dict(os.environ, PATH=binp + os.pathsep + os.environ["PATH"])
    subprocess.run(["bgzip", "-f", gff], check=True, env=env)
    subprocess.run(["tabix", "-f", "-p", "gff", gff + ".gz"], check=True, env=env)
    out = os.path.join(d, "vep.json")
    r = subprocess.run(["vep", "-i", vcf, "--gff", gff + ".gz", "--fasta", fasta, "--hgvs", "--json", "-o", out,
                        "--force_overwrite", "--no_stats", "--distance", "0", "--fork", "4"],
                       capture_output=True, text=True, env=env)
    if r.returncode:
        sys.exit("vep failed: " + r.stderr[-2000:])
    result = {}
    for line in open(out):
        j = json.loads(line)
        i = int(j["id"][4:])
        for t in j.get("transcript_consequences", []):
            if t.get("transcript_id") == f"T{i}":
                result[i] = {"hgvsp": strip(t.get("hgvsp")), "hgvsc": strip(t.get("hgvsc")),
                             "terms": sorted(t.get("consequence_terms", []))}
    return result


def strip(desc):
    return desc.split(":", 1)[1] if desc and ":" in desc else desc


def compare_vep(duck, vep):
    """DuckVEP's single-change protein name against VEP's on the merged record."""
    if vep is None or duck.get("status") != "predicted" or not duck.get("hgvsp"):
        return "n/a"
    if duck["hgvsp"].startswith("p.[") or duck["hgvsp"].startswith("p.("):
        d = duck["hgvsp"].replace("p.(", "p.", 1).rstrip(")") if duck["hgvsp"].startswith("p.(") else duck["hgvsp"]
        if duck["hgvsp"].startswith("p.["):
            return "n/a separated"
    else:
        d = duck["hgvsp"]
    v = vep.get("hgvsp")
    if v is None:
        return f"vep none (duckvep {duck['hgvsp']})"
    return "ok" if d == v else f"DIFF duckvep {duck['hgvsp']} vep {v}"


def compare(case, duck, joint, mut):
    """Verdicts for one case: (composition, dna, protein); each 'ok', a documented difference, or 'MISMATCH ...'."""
    if joint is None:
        comp = "ok" if duck["status"] == "edit_conflict" else f"MISMATCH applier conflict, duckvep {duck['status']}"
        return comp, "n/a", "n/a"
    if duck["status"] != "predicted":
        return f"MISMATCH applier composes, duckvep {duck['status']}/{duck['reason']}", "n/a", "n/a"
    cds = joint[U5:len(joint) - U3]
    comp = "ok" if duck["cds"] == cds else f"MISMATCH cds {duck['cds']} != {cds}"
    if mut is None or mut.get("errors"):
        return comp, f"oracle error {mut and mut.get('errors')}", "n/a"
    dna_m, dna_d = strip(mut["dna"]), duck["hgvsc"]
    if dna_d == dna_m:
        dna = "ok"
    elif dna_d is None:
        dna = f"MISMATCH duckvep hgvsc NULL ({duck['hgvsc_status']}), mutalyzer {dna_m}"
    elif "inv" in dna_m:
        dna = "notation: mutalyzer inversion"
    elif "[" in dna_d or "[" in dna_m:
        dna = f"notation: grouping duckvep {dna_d} mutalyzer {dna_m}"
    else:
        dna = f"MISMATCH duckvep {dna_d} mutalyzer {dna_m}"
    if not cds.startswith("ATG"):  # start lost: neither tool predicts a product
        return comp, dna, "n/a start lost"
    joint_protein = translate(joint[U5:])  # the joint allele's product, read through the 3' UTR when needed
    predicted = (mut.get("predicted") or "").rstrip("X")
    prot = "ok"
    if duck["protein"] != joint_protein:
        prot = f"MISMATCH duckvep protein {duck['protein']} vs joint translation {joint_protein}"
    elif predicted and not joint_protein.startswith(predicted.rstrip("*")):
        prot = f"MISMATCH mutalyzer predicted {predicted} vs joint translation {joint_protein}"
    else:
        p_m, p_d = strip(mut.get("p")), duck["hgvsp"]
        if p_d != p_m:
            if p_m == "p.?":
                prot = "notation: mutalyzer p.?"
            elif p_d and p_d.startswith("p.[("):
                prot = f"notation: separated duckvep {p_d} mutalyzer {p_m}"
            else:
                prot = f"DIFF hgvsp duckvep {p_d} mutalyzer {p_m}"
    return comp, dna, prot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extension", required=True)
    ap.add_argument("--mutalyzer-python", required=True)
    ap.add_argument("--cases", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--work", default=None)
    ap.add_argument("--vep-prefix", default=None, help="conda prefix with vep, bgzip and tabix; omit to skip the VEP leg")
    a = ap.parse_args()
    work = a.work or tempfile.mkdtemp(prefix="compound-hgvs-fuzz-")
    os.makedirs(work, exist_ok=True)
    rng = random.Random(a.seed)
    cases = [random_case(rng, i) for i in range(a.cases)]
    joint = []
    for c in cases:
        conflict, seq = apply_literal(c["u5"] + c["cds"] + c["u3"], c["edits"])
        joint.append(None if conflict else seq)
    duck = run_duckvep(os.path.abspath(a.extension), cases)
    mut, described = run_mutalyzer(a.mutalyzer_python, cases, joint, work)
    vep = run_vep(a.vep_prefix, cases, joint, work) if a.vep_prefix else {}
    tally, bad = {}, 0
    with open(os.path.join(work, "cases.jsonl"), "w") as f:
        for i, c in enumerate(cases):
            d = duck.get(i, {"status": "missing", "reason": None})
            verdict = compare(c, d, joint[i], mut.get(i)) + ((compare_vep(d, vep.get(i)),) if a.vep_prefix else ())
            for name, v in zip(("composition", "dna", "protein", "vep_protein"), verdict):
                key = (name, v.split(" ")[0] if v.startswith("MISMATCH") else v.split(":")[0] if ":" in v else v)
                tally[key] = tally.get(key, 0) + 1
                bad += v.startswith("MISMATCH")
            f.write(json.dumps({**c, "joint": joint[i], "oracle_input": described.get(i), "duckvep": d,
                                "mutalyzer": mut.get(i), "vep": vep.get(i), "verdict": verdict}, default=str) + "\n")
    conflicts = sum(j is None for j in joint)
    print(f"cases {len(cases)} (seed {a.seed}); applier conflicts {conflicts}; work {work}")
    for (name, v), n in sorted(tally.items()):
        print(f"  {name:12s} {v:40s} {n}")
    print("MISMATCHES", bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
