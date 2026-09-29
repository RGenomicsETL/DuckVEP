#!/usr/bin/env python3
"""Origin-focused circular workload: throughput, output equality and peak memory.

A synthetic circular region carries many transcripts and regulation features
packed around its origin, with most events falling within a few kilobases of it.
The same world is executed in three frames:

  lifted      origin inside the object cloud, so most objects wrap;
  lifted_rot  the same world rotated by a small offset, a different wrap pattern;
  linear      the world rotated so the origin sits in an empty gap; nothing wraps,
              the model runs on the ordinary linear kernel and is the control.

Each frame is annotated with the fused consequence + HGVS path on 1 and N threads.
Every run is a fresh child process that loads an immutable, hashed copy of the
extension, so peak resident memory is that run's alone. The output checksum keys
rows by transcript stable identifier, not by frame-specific ordinal, so it must
be identical across thread counts and across frames.

    python3 benchmarks/duckvep_circular_origin.py \\
        --extension build/release/extension/duckvep/duckvep.duckdb_extension \\
        --out benchmarks/data/duckvep_circular_origin.csv
"""
import argparse
import csv
import hashlib
import json
import os
import random
import resource
import shutil
import subprocess
import sys
import tempfile
import time

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def build_world(length, n_tx, n_features, seed):
    """Objects in unrolled coordinates around u = length (the origin)."""
    rng = random.Random(seed)
    ref = "".join(rng.choice("ACGT") for _ in range(length))
    txs = []
    for i in range(n_tx):
        strand = rng.choice((1, -1))
        n_ex = rng.randint(1, 5)
        cursor = length + rng.randint(-15000, 15000)
        exons = []
        for _ in range(n_ex):
            ln = rng.randint(40, 400)
            exons.append((cursor, cursor + ln - 1))
            cursor += ln + rng.randint(80, 900)
        cdna_len = sum(b - a + 1 for a, b in exons)
        coding = rng.random() < 0.7 and cdna_len >= 60
        cds = None
        if coding:
            codons = rng.randint(10, (cdna_len - 6) // 3)
            a = rng.randint(1, cdna_len - codons * 3 + 1)
            cds = (a, a + codons * 3 - 1)
        txs.append(dict(name=f"T{i:04d}", strand=strand, exons=exons, cds=cds))
    feats = []
    for i in range(n_features):
        s = length + rng.randint(-15000, 15000)
        feats.append(dict(name=f"R{i:03d}", start=s, end=s + rng.randint(30, 900),
                          kind=rng.choice((1, 2))))
    return ref, txs, feats


def frame_tables(length, ref, txs, feats, k):
    """Model rows for the world rotated forward by k bases."""
    def pos(u):
        return (u - 1 + k) % length + 1

    def base(u):
        return ref[(u - 1) % length]

    rows = []
    for t in txs:
        order = t["exons"] if t["strand"] > 0 else list(reversed(t["exons"]))
        cdna_u = []
        for a, b in order:
            cdna_u.extend(range(a, b + 1) if t["strand"] > 0 else range(b, a - 1, -1))
        seq = "".join(base(u) if t["strand"] > 0 else COMP[base(u)] for u in cdna_u)
        ex_rows, cum, before = [], 0, 0
        for a, b in order:
            ln = b - a + 1
            cs, ce = cum + 1, cum + ln
            phase = end_phase = -1
            if t["cds"]:
                ca, cb = t["cds"]
                if ce >= ca and cs <= cb:
                    first, last = max(cs, ca), min(ce, cb)
                    if cs >= ca:
                        phase = before % 3
                    if ce < cb:
                        end_phase = (before + last - first + 1) % 3
                    before += last - first + 1
            ex_rows.append((pos(a), pos(b), cs, ce, phase, end_phase))
            cum = ce
        lo, hi = t["exons"][0][0], t["exons"][-1][1]
        row = dict(name=t["name"], strand=t["strand"], start=pos(lo), end=pos(hi), exons=ex_rows)
        if t["cds"]:
            ca, cb = t["cds"]
            row["cds_start"] = pos(cdna_u[(ca if t["strand"] > 0 else cb) - 1])
            row["cds_end"] = pos(cdna_u[(cb if t["strand"] > 0 else ca) - 1])
            row["cds_seq"] = seq[ca - 1:cb]
            row["pre"] = seq[:ca - 1]
            row["post"] = seq[cb:]
        rows.append(row)
    rows.sort(key=lambda r: (r["start"], r["name"]))
    frows = []
    for f in feats:
        frows.append((pos(f["start"]), pos(f["end"]), f["kind"], f["name"]))
    frows.sort()
    return rows, frows


def build_events(length, ref, n_events, seed, seam_exclusion):
    """Events keyed by unrolled position; 70% within 4 kb of the origin."""
    rng = random.Random(seed ^ 0x5bd1e995)
    events = []
    for i in range(n_events):
        if rng.random() < 0.7:
            u = length + rng.randint(-4000, 4000)
        else:
            u = rng.randint(1, length)
        kind = rng.random()
        base = ref[(u - 1) % length]
        if kind < 0.6:
            rl, alt = 1, rng.choice([c for c in "ACGT" if c != base])
        elif kind < 0.75:
            rl = rng.randint(2, 12)
            alt = base
        elif kind < 0.9:
            rl = 1
            alt = base + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 9)))
        else:
            rl = rng.randint(2, 4)
            alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 5)))
        refseq = "".join(ref[(u - 1 + j) % length] for j in range(rl))
        if refseq == alt:
            continue
        events.append((i, u, refseq, alt))
    return events


def write_frame(directory, tag, length, ref, txs, feats, events, k, seam_margin):
    import duckdb
    rows, frows = frame_tables(length, ref, txs, feats, k)
    rotated_ref = ref[length - k:] + ref[:length - k] if k else ref
    with open(os.path.join(directory, f"{tag}.fa"), "w") as h:
        h.write(">C\n" + rotated_ref + "\n")
    with open(os.path.join(directory, f"{tag}.fa.fai"), "w") as h:
        h.write(f"C\t{length}\t3\t{length}\t{length + 1}\n")
    con = duckdb.connect()
    con.execute("CREATE TABLE tx(transcript_index UINTEGER, stable_id VARCHAR, seq_region UINTEGER, transcript_start UBIGINT, "
                "transcript_end UBIGINT, strand TINYINT, gene_index UINTEGER, transcript_flags UBIGINT, cds_start UBIGINT, "
                "cds_end UBIGINT, cds_sequence BLOB, codon_table UTINYINT, pre_cds_sequence BLOB, post_cds_sequence BLOB)")
    con.execute("CREATE TABLE ex(transcript_index UINTEGER, exon_start UBIGINT, exon_end UBIGINT, exon_cdna_start UBIGINT, "
                "exon_cdna_end UBIGINT, phase TINYINT, end_phase TINYINT)")
    tx_batch, ex_batch = [], []
    for i, r in enumerate(rows):
        coding = "cds_seq" in r
        tx_batch.append((i, r["name"], 0, r["start"], r["end"], r["strand"], i, 3 if coding else 0,
                         r.get("cds_start"), r.get("cds_end"), r["cds_seq"].encode() if coding else None,
                         1 if coding else None, r["pre"].encode() if coding else None,
                         r["post"].encode() if coding else None))
        ex_batch.extend((i, *e) for e in r["exons"])
    con.executemany("INSERT INTO tx VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", tx_batch)
    con.executemany("INSERT INTO ex VALUES (?,?,?,?,?,?,?)", ex_batch)
    con.execute("CREATE TABLE ft(regulation_feature_index UINTEGER, seq_region UINTEGER, feature_start UINTEGER, "
                "feature_end UINTEGER, feature_kind UTINYINT)")
    con.executemany("INSERT INTO ft VALUES (?,0,?,?,?)",
                    [(i, s, e, kind) for i, (s, e, kind, _name) in enumerate(frows)])
    con.execute("CREATE TABLE ev(event_index UBIGINT, seq_region UINTEGER, position UBIGINT, reference VARCHAR, "
                "alternate VARCHAR, end_position UBIGINT, structural_type VARCHAR, copy_change VARCHAR, "
                "mate_seq_region UINTEGER, mate_position UBIGINT)")
    csv_path = os.path.join(directory, tag + "_ev.csv")
    with open(csv_path, "w") as handle:
        for i, u, refseq, alt in events:
            handle.write("%d,0,%d,%s,%s,,,,,\n" % (i, (u - 1 + k) % length + 1, refseq, alt))
    con.execute("COPY ev FROM '%s' (FORMAT csv, HEADER false)" % csv_path)
    kept = len(events)
    os.remove(csv_path)
    for name in ("tx", "ex", "ft", "ev"):
        con.execute(f"COPY {name} TO '{os.path.join(directory, tag + '_' + name)}.parquet' (FORMAT parquet)")
    con.execute("CREATE TABLE names AS SELECT transcript_index AS ordinal, 't' AS kind, stable_id FROM tx")
    con.executemany("INSERT INTO names VALUES (?, 'f', ?)", [(i, name) for i, (_s, _e, _k, name) in enumerate(frows)])
    con.execute(f"COPY names TO '{os.path.join(directory, tag + '_names')}.parquet' (FORMAT parquet)")
    return kept


CHILD_SQL = """
CREATE TABLE tx AS SELECT * EXCLUDE (stable_id) FROM read_parquet('{d}/{tag}_tx.parquet');
CREATE TABLE ex AS SELECT * FROM read_parquet('{d}/{tag}_ex.parquet');
CREATE TABLE ft AS SELECT * FROM read_parquet('{d}/{tag}_ft.parquet');
CREATE TABLE ev AS SELECT * FROM read_parquet('{d}/{tag}_ev.parquet');
CREATE TABLE names AS SELECT * FROM read_parquet('{d}/{tag}_names.parquet');
CREATE TABLE regions AS SELECT 0::UINTEGER AS seq_region, {length}::UBIGINT AS sequence_length, 'C' AS seq_region_name, true AS circular;
"""


def vm_hwm_mb():
    """Peak resident set of this process image, from the kernel high-water mark."""
    with open("/proc/self/status") as h:
        for line in h:
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) / 1024.0
    return float("nan")


def child(args):
    import duckdb
    con = duckdb.connect(config={"allow_unsigned_extensions": "true", "threads": str(args.threads[0] if isinstance(args.threads, list) else args.threads)})
    con.execute(f"LOAD '{args.extension}'")
    con.execute(CHILD_SQL.format(d=args.directory, tag=args.tag, length=args.length))
    con.execute("""SELECT loaded FROM duckvep_model_load('bench',
      'SELECT * FROM regions',
      'SELECT * FROM tx ORDER BY seq_region, transcript_start, transcript_index',
      'SELECT * FROM ex ORDER BY transcript_index, exon_cdna_start',
      interval_feature_query := 'SELECT * FROM ft ORDER BY regulation_feature_index',
      reference_fasta := '""" + os.path.join(args.directory, args.tag + ".fa") + "')").fetchall()
    n_events = con.execute("SELECT count(*) FROM ev").fetchone()[0]
    loaded_rss = vm_hwm_mb()
    sql = """SELECT count(*), sum(hash(a.event_index, coalesce(n.stable_id, m.stable_id),
        a.consequence_mask, a.region_mask, a.impact_code, coalesce(a.cdna_position, 0), coalesce(a.cds_position, 0),
        coalesce(a.protein_position, 0), coalesce(a.transcript_hgvs, ''), coalesce(a.protein_hgvs, ''),
        coalesce(a.hgvs_shift, 0), a.nmd_prediction_code)::HUGEINT)
      FROM query(duckvep_annotate_sql('ev', 'bench', struct_pack(hgvs := {hgvs}, upstream_distance := 1000, downstream_distance := 1000))) a
      LEFT JOIN names n ON n.kind = 't' AND n.ordinal = a.transcript_index
      LEFT JOIN names m ON m.kind = 'f' AND m.ordinal = a.regulation_feature_index""".format(hgvs="true" if args.hgvs else "false")
    times = []
    rows = checksum = None
    for _ in range(args.passes):
        t0 = time.perf_counter()
        rows, checksum = con.execute(sql).fetchone()
        times.append(time.perf_counter() - t0)
    print(json.dumps(dict(events=n_events, rows=rows, checksum=str(checksum), seconds=times,
                          peak_rss_mb=vm_hwm_mb(), loaded_rss_mb=loaded_rss)))


def run_child(args, tag, threads, hgvs):
    cmd = [sys.executable, os.path.abspath(__file__), "--child", "--extension", args.extension_copy,
           "--directory", args.directory, "--tag", tag, "--threads", str(threads), "--length", str(args.length),
           "--passes", str(args.passes)] + (["--hgvs"] if hgvs else [])
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(err[-2000:])
    return json.loads(out.strip().splitlines()[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extension", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--length", type=int, default=100000)
    ap.add_argument("--transcripts", type=int, default=400)
    ap.add_argument("--features", type=int, default=120)
    ap.add_argument("--events", type=int, default=400000)
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 8])
    ap.add_argument("--passes", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--child", action="store_true")
    ap.add_argument("--directory")
    ap.add_argument("--tag")
    ap.add_argument("--hgvs", action="store_true")
    args = ap.parse_args()
    if args.child:
        child(args)
        return
    work = tempfile.mkdtemp(prefix="duckvep_circular_origin_")
    args.directory = work
    # Immutable copy: the benchmark never loads a file a later build could overwrite.
    args.extension_copy = os.path.join(work, "duckvep.duckdb_extension")
    shutil.copyfile(args.extension, args.extension_copy)
    os.chmod(args.extension_copy, 0o444)
    sha = hashlib.sha256(open(args.extension_copy, "rb").read()).hexdigest()
    ref, txs, feats = build_world(args.length, args.transcripts, args.features, args.seed)
    events = build_events(args.length, ref, args.events, args.seed, 0)
    frames = {"lifted": 0, "lifted_rot": 777, "linear": args.length // 2}
    # The linear control cannot hold an event that crosses its seam, and the two
    # lifted frames can. All frames annotate the same events: drop any event within
    # 300 bases of a seam or across it in any frame.
    def clear_of_seams(u, refseq):
        return all(300 < (u - 1 + k) % args.length + 1 and
                   (u - 1 + k) % args.length + len(refseq) <= args.length - 300 for k in frames.values())
    events = [e for e in events if clear_of_seams(e[1], e[2])]
    kept = {tag: write_frame(work, tag, args.length, ref, txs, feats, events, k, 300)
            for tag, k in frames.items()}
    import datetime
    import duckdb
    revision = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                              cwd=os.path.dirname(os.path.abspath(__file__))).stdout.strip()
    cpu = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), "unknown")
    meta = dict(run_date=datetime.date.today().isoformat(), source_revision=revision, duckdb_version=duckdb.__version__,
                cpu=cpu, length=args.length, transcripts=args.transcripts, features=args.features)
    results = []
    for tag in frames:
        for hgvs in (False, True):
            for threads in args.threads:
                r = run_child(args, tag, threads, hgvs)
                secs = sorted(r["seconds"])
                median = secs[len(secs) // 2]
                results.append(dict(**meta, frame=tag, hgvs=hgvs, threads=threads, events=r["events"], rows=r["rows"],
                                    seconds_median=round(median, 3), events_per_second=round(r["events"] / median),
                                    peak_rss_mb=round(r["peak_rss_mb"], 1),
                                    loaded_rss_mb=round(r["loaded_rss_mb"], 1), checksum=r["checksum"],
                                    extension_sha256=sha))
                print(results[-1], flush=True)
    if args.out:
        with open(args.out, "w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(results[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(results)
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
