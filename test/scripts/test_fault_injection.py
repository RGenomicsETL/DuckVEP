#!/usr/bin/env python3
"""Fail the Nth native allocation during model load and annotation, under ASan.

The extension is built with -DDUCKVEP_FAULT_INJECTION (make test_fault_injection),
which adds duckvep_fault_arm(n): the nth budget allocation after arming returns
NULL. Every allocation the fixture performs is failed once, in order, in a single
DuckDB CLI session so the registry survives from one failure to the next. After
each failure the harness checks that

  * the statement failed with an error (or, when the failed allocation sits on a
    tolerated path, produced exactly the golden result);
  * nothing was published: the model name is unknown after a failed load, and a
    loaded model still returns the golden result after a failed annotation;
  * the native budget is back to its baseline (no charged byte leaked);
  * ASan and LeakSanitizer stay silent for the whole session.

Set DUCKVEP_FAULT_TRACE=1 (the harness does) so each failure prints the calling
stack; the harness reports how many distinct call sites were failed.
"""
import argparse
import collections
import concurrent.futures
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

FIXTURE = """
CREATE TABLE f_src AS SELECT * FROM (VALUES
  (0, 4, 1, 13, 'G', 'GT'), (1, 17, 1, 14, 'A', 'AG'), (2, 18, 1, 14, 'A', 'AT'),
  (3, 73, 2, 14, 'A', 'AG'), (4, 74, 2, 14, 'A', 'AT'), (5, 76, 2, 14, 'A', 'AGCC'),
  (6, 284, 6, 13, 'G', 'GT'), (7, 353, 9, 14, 'A', 'AG'), (8, 354, 9, 14, 'A', 'AT'),
  (9, 356, 9, 14, 'A', 'AGCC'), (10, 5400, 1, 14, 'NAAG', 'N'), (11, 32261, 1, 13, 'G', 'GAC')
) s(i, event_index, codon_table, position, reference, alternate);
CREATE TABLE f_tr AS SELECT i::UINTEGER transcript_index, i::UINTEGER seq_region,
  11::UBIGINT transcript_start, 22::UBIGINT transcript_end, 1::TINYINT strand,
  i::UINTEGER gene_index, (CASE WHEN i = 3 THEN 8 ELSE 3 END)::UBIGINT transcript_flags,
  (CASE WHEN i = 3 THEN NULL ELSE transcript_start END)::UBIGINT cds_start,
  (CASE WHEN i = 3 THEN NULL ELSE transcript_end END)::UBIGINT cds_end,
  (CASE i WHEN 3 THEN NULL WHEN 10 THEN 'ATGNAAGCCTAA' WHEN 11 THEN 'ATGNNAGCCTAA' ELSE 'ATGAAAGCCTAA' END)::BLOB cds_sequence,
  (CASE WHEN i = 3 THEN NULL ELSE codon_table END)::UTINYINT codon_table, (CASE WHEN i = 3 THEN NULL ELSE '' END)::BLOB pre_cds_sequence,
  (CASE WHEN i = 3 THEN NULL ELSE '' END)::BLOB post_cds_sequence FROM f_src;
CREATE TABLE f_pe AS SELECT 0::UINTEGER transcript_index, 2::UINTEGER protein_position,
  'W'::VARCHAR alternate_amino_acid;
CREATE TABLE f_mirna AS SELECT 3::UINTEGER transcript_index, 12::UBIGINT mature_mirna_start,
  15::UBIGINT mature_mirna_end;
CREATE TABLE f_feat AS SELECT * FROM (VALUES
  (0::UINTEGER, 1::UINTEGER, 5::UINTEGER, 30::UINTEGER, 1::UTINYINT),
  (1::UINTEGER, 2::UINTEGER, 5::UINTEGER, 30::UINTEGER, 2::UTINYINT))
  t(regulation_feature_index, seq_region, feature_start, feature_end, feature_kind);
CREATE TABLE f_ev AS SELECT event_index::UBIGINT event_index, i::UINTEGER seq_region,
  position::UBIGINT AS position, reference, alternate, NULL::UBIGINT end_position,
  NULL::VARCHAR structural_type, NULL::VARCHAR copy_change, NULL::UINTEGER mate_seq_region,
  NULL::UBIGINT mate_position FROM f_src
  UNION ALL SELECT 900::UBIGINT, 1::UINTEGER, 12::UBIGINT, NULL, '<DEL>', 20::UBIGINT, NULL, NULL, NULL, NULL
  UNION ALL SELECT 901::UBIGINT, 1::UINTEGER, 12::UBIGINT, NULL, 'N]2:20]', NULL, NULL, NULL, 2::UINTEGER, 20::UBIGINT
  UNION ALL SELECT 902::UBIGINT, 2::UINTEGER, 14::UBIGINT, NULL, '<DUP>', 25::UBIGINT, NULL, NULL, NULL, NULL;
-- A circular region with a transcript, an exon and a regulatory feature that cross the
-- origin: the model is executed on a lifted linear copy (duckvep_lift.c). Region 1 is
-- table1 (32 bp), circular; the transcript covers 27..32 and 1..6.
CREATE TABLE f_ctr AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region,
  27::UBIGINT transcript_start, 6::UBIGINT transcript_end, 1::TINYINT strand,
  0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 27::UBIGINT cds_start, 6::UBIGINT cds_end,
  'ATGAAAGCCTAA'::BLOB cds_sequence, 1::UTINYINT codon_table, ''::BLOB pre_cds_sequence,
  ''::BLOB post_cds_sequence;
CREATE TABLE f_cfeat AS SELECT * FROM (VALUES
  (0::UINTEGER, 1::UINTEGER, 30::UINTEGER, 3::UINTEGER, 1::UTINYINT))
  t(regulation_feature_index, seq_region, feature_start, feature_end, feature_kind);
CREATE TABLE f_cev AS SELECT * FROM (VALUES
  (0::UBIGINT, 1::UINTEGER, 29::UBIGINT, 'A', 'G'), (1::UBIGINT, 1::UINTEGER, 31::UBIGINT, 'A', 'T'),
  (2::UBIGINT, 1::UINTEGER, 2::UBIGINT, 'A', 'C'), (3::UBIGINT, 1::UINTEGER, 5::UBIGINT, 'A', 'AGCC'),
  (4::UBIGINT, 1::UINTEGER, 30::UBIGINT, 'AA', 'A'), (5::UBIGINT, 1::UINTEGER, 8::UBIGINT, 'A', 'G'),
  (6::UBIGINT, 1::UINTEGER, 28::UBIGINT, 'A', 'AGCC'), (7::UBIGINT, 1::UINTEGER, 3::UBIGINT, 'AAA', 'A'))
  t(event_index, seq_region, position, reference, alternate);
CREATE VIEW f_cev_full AS SELECT event_index, seq_region, position, reference, alternate,
  NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change,
  NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position FROM f_cev;
-- Phased calls for duckvep_haplotypes on the first model: a same-codon pair, a frame-opening
-- insertion restored by a downstream deletion, and a stop gain, on three transcripts (the
-- frame classifier and the per-carrier consequence and impact lists).
CREATE TABLE f_hcalls AS SELECT (row_number() OVER ())::BIGINT event_index, seq_region::INT seq_region,
  position::BIGINT AS position, reference, alternate, 1 alt_index, seq_region::INT transcript_index,
  0 sample_index, [1,1]::INTEGER[] alleles, [false,true]::BOOLEAN[] phase_before, NULL::BIGINT phase_set
  FROM (VALUES (4, 14, 'A', 'G'), (4, 15, 'A', 'T'), (1, 13, 'G', 'GT'), (1, 16, 'AG', 'A'),
   (2, 14, 'A', 'T'), (2, 15, 'A', 'G')) t(seq_region, position, reference, alternate);
CREATE VIEW f_small AS SELECT * FROM f_ev WHERE reference IS NOT NULL
  ORDER BY seq_region, position, event_index;
"""

LOAD = """SELECT 'LOADED ' || loaded FROM duckvep_model_load('fm',
 'SELECT seq_region, 32::UBIGINT sequence_length, CASE seq_region WHEN 10 THEN ''residual_naa'' WHEN 11 THEN ''residual_nna'' ELSE ''table'' || seq_region END seq_region_name FROM f_tr ORDER BY seq_region',
 'SELECT * FROM f_tr ORDER BY transcript_index',
 'SELECT transcript_index, 11::UBIGINT exon_start, 22::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase FROM f_tr ORDER BY transcript_index',
 peptide_edit_query := 'SELECT * FROM f_pe',
 mature_mirna_query := 'SELECT * FROM f_mirna',
 interval_feature_query := 'SELECT * FROM f_feat ORDER BY regulation_feature_index',
 reference_fasta := '{fasta}');"""

LOAD_CIRCULAR = """SELECT 'LOADED ' || loaded FROM duckvep_model_load('fc',
 'SELECT 1::UINTEGER seq_region, 32::UBIGINT sequence_length, ''table1'' seq_region_name, true circular',
 'SELECT * FROM f_ctr',
 'SELECT 0::UINTEGER transcript_index, 27::UBIGINT exon_start, 6::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase',
 interval_feature_query := 'SELECT * FROM f_cfeat',
 reference_fasta := '{fasta}');"""

# Each query prints one line: 'RESULT <name> <rows> <hash>'.
QUERIES = {
    "hgvs": """SELECT 'RESULT hgvs ' || count(*) || ' ' || hash(list(a ORDER BY hash(a))) FROM
 query(duckvep_annotate_sql('f_ev', 'fm', struct_pack(hgvs := true, upstream_distance := 0, downstream_distance := 0))) a;""",
    "regulation": """SELECT 'RESULT regulation ' || count(*) || ' ' || hash(list(a ORDER BY hash(a))) FROM
 query(duckvep_annotate_sql('f_ev', 'fm')) a;""",
    "projected": """SELECT 'RESULT projected ' || count(*) || ' ' || hash(list(a ORDER BY hash(a))) FROM
 query(duckvep_annotate_projected_sql('f_small', 'fm')) a;""",
    # Whole-haplotype classifier over phased calls, with per-carrier consequence and impact lists.
    "haplotype": """SELECT 'RESULT haplotype ' || count(*) || ' ' || hash(list(h ORDER BY hash(h))) FROM
 duckvep_haplotypes('SELECT * FROM f_hcalls', 'fm') h;""",
    # Transcript discovery for phased calls: model-name copy, hit buffer and result list.
    "discovery": """SELECT 'RESULT discovery ' || count(*) || ' ' || hash(list(d ORDER BY hash(d))) FROM
 (SELECT event_index, duckvep_coding_transcripts('fm', seq_region::BIGINT, position::BIGINT,
  reference, alternate) AS transcripts FROM f_ev WHERE reference IS NOT NULL) d;""",
    # Circular model executed on a lifted copy: lift_resolve and its HGVS and projection copies.
    "lifted_hgvs": """SELECT 'RESULT lifted_hgvs ' || count(*) || ' ' || hash(list(a ORDER BY hash(a))) FROM
 query(duckvep_annotate_sql('f_cev_full', 'fc', struct_pack(hgvs := true, upstream_distance := 0, downstream_distance := 0))) a;""",
    "lifted": """SELECT 'RESULT lifted ' || count(*) || ' ' || hash(list(a ORDER BY hash(a))) FROM
 query(duckvep_annotate_sql('f_cev_full', 'fc')) a;""",
    "lifted_projected": """SELECT 'RESULT lifted_projected ' || count(*) || ' ' || hash(list(a ORDER BY hash(a))) FROM
 query(duckvep_annotate_projected_sql('f_cev', 'fc')) a;""",
}
# The model each annotation phase runs on, and the query that proves a load published nothing.
MODEL_OF = {"haplotype": "fm", "discovery": "fm", "hgvs": "fm", "regulation": "fm", "projected": "fm",
            "lifted_hgvs": "fc", "lifted": "fc", "lifted_projected": "fc"}
LOAD_PHASES = {"load": ("fm", "hgvs"), "load_circular": ("fc", "lifted_hgvs")}

BASELINE_SQL = "SELECT 'TOTAL ' || current_bytes FROM duckvep_native_budget() WHERE owner = 'total';"


def run(extension, script, duckdb, env_extra=None, timeout=3600):
    asan = subprocess.check_output(["gcc", "-print-file-name=libasan.so"], text=True).strip()
    env = dict(os.environ)
    env["LD_PRELOAD"] = asan
    env["ASAN_OPTIONS"] = "detect_leaks=1:halt_on_error=1:abort_on_error=0:exitcode=86"
    env["LSAN_OPTIONS"] = "print_suppressions=0:exitcode=87"
    env["DUCKVEP_FAULT_TRACE"] = "1"
    env.update(env_extra or {})
    text = ".bail off\nLOAD '%s';\n%s" % (extension, script)
    with tempfile.NamedTemporaryFile("w", suffix=".sql", delete=False) as handle:
        handle.write(text)
        path = handle.name
    try:
        proc = subprocess.run([duckdb, "-unsigned", "-csv", "-noheader", "-f", path],
                              env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, timeout=timeout)
    finally:
        os.unlink(path)
    return proc.returncode, proc.stdout


def value(output, prefix, index=-1):
    lines = [line[len(prefix):].strip() for line in output.splitlines() if line.startswith(prefix)]
    if not lines:
        raise SystemExit("missing %r in output:\n%s" % (prefix, output[-3000:]))
    return lines[index]


def iteration(args, phase, index, loads, base, allowed, golden):
    """Fail allocation `index` of `phase` in a fresh DuckDB process; return problems."""
    script = [FIXTURE, "SELECT duckvep_worker_limits_set(6, 134217728, 67108864, 0);"]
    if phase in LOAD_PHASES:
        model, probe = LOAD_PHASES[phase]
        script += ["SELECT duckvep_fault_arm(%d);" % index, loads[model],
                   "SELECT 'FIRED ' || duckvep_fault_fired();", "SELECT duckvep_fault_arm(0);",
                   BASELINE_SQL.replace("TOTAL", "LEAK"),
                   QUERIES[probe],  # nothing may have been published
                   "SELECT 'DROP ' || duckvep_model_drop('%s');" % model]
    else:
        model = MODEL_OF[phase]
        script += [loads[model], "SELECT duckvep_fault_arm(%d);" % index, QUERIES[phase],
                   "SELECT 'FIRED ' || duckvep_fault_fired();", "SELECT duckvep_fault_arm(0);",
                   QUERIES[phase].replace("RESULT ", "AGAIN "),  # the model must still answer
                   "SELECT 'DROP ' || duckvep_model_drop('%s');" % model]
    script += [BASELINE_SQL.replace("TOTAL", "END")]
    rc, out = run(args.extension, "\n".join(script), args.duckdb)
    tag = "%s#%d" % (phase, index)
    problems = []
    # The DuckDB shell exits 1 after any statement error; sanitizers use 86/87.
    if rc not in (0, 1):
        problems.append("%s: exit status %d" % (tag, rc))
    for marker in ("AddressSanitizer", "LeakSanitizer", "runtime error"):
        if marker in out:
            problems.append("%s: sanitizer report (%s)" % (tag, marker))
    fired = re.search(r"^FIRED (\d+)", out, re.M)
    if not fired or int(fired.group(1)) != 1:
        problems.append("%s: allocation was not failed" % tag)
        return problems, None, out
    site = None
    frames = re.findall(r"#\d+ 0x[0-9a-f]+ in (\S+) (\S+?):(\d+)", out.split("DUCKVEP_FAULT_FIRED", 1)[-1])
    skip = {"duckvep_sql_resize", "duckvep_scalar_grow", "duckvep_budget_strdup",
            "duckvep_budget_strndup", "cr_grow_", "__sanitizer_print_stack_trace", "fault_fires"}
    for function, path, line in frames:
        if "duckvep_budget.c" not in path and function not in skip:
            site = (os.path.basename(path), function, int(line))
            break
    end = int(value(out, "END ")) - base
    if end not in allowed:
        problems.append("%s: budget not restored after drop (%+d bytes)" % (tag, end))
    if phase in LOAD_PHASES:
        model, probe = LOAD_PHASES[phase]
        errored = "LOADED true" not in out
        if int(value(out, "LEAK ")) - base != 0:
            problems.append("%s: bytes still charged after the failed load" % tag)
        if errored:
            if "DROP false" not in out:
                problems.append("%s: failed load published a model" % tag)
            if "RESULT %s" % probe in out:
                problems.append("%s: failed load left an annotatable model" % tag)
            if not re.search(r"(?i)error", out):
                problems.append("%s: failed load gave no explicit error" % tag)
        elif "RESULT %s %s" % (probe, golden[probe]) not in out:
            problems.append("%s: tolerated failure produced a different model" % tag)
        return problems, (site, "load error" if errored else "load tolerated"), out
    first = re.search(r"^RESULT %s (.*)$" % phase, out, re.M)
    again = re.search(r"^AGAIN %s (.*)$" % phase, out, re.M)
    if not again or again.group(1) != golden[phase]:
        problems.append("%s: model unusable after the failure: %s" % (tag, again and again.group(0)))
    if first and first.group(1) != golden[phase]:
        problems.append("%s: tolerated failure changed the result" % tag)
    if not first and not re.search(r"(?i)(out of memory|capacity error|error)", out):
        problems.append("%s: failed annotation gave no explicit error" % tag)
    return problems, (site, "annotation tolerated" if first else "annotation error"), out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extension", default=os.path.join(ROOT, "build/fault/duckvep.duckdb_extension"))
    parser.add_argument("--duckdb", default="duckdb")
    parser.add_argument("--fasta", default=os.path.join(ROOT, "test/data/duckvep/indel_translation.fa"))
    parser.add_argument("--limit", type=int, default=0, help="stop after this many failures per phase (debugging)")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    loads = {"fm": LOAD.format(fasta=args.fasta), "fc": LOAD_CIRCULAR.format(fasta=args.fasta)}

    # Count pass: allocations per phase, golden results, and the bytes a pooled
    # annotation worker legitimately retains (idle limit 0 trims its buffers).
    script = [FIXTURE, "SELECT duckvep_worker_limits_set(6, 134217728, 67108864, 0);",
              BASELINE_SQL.replace("TOTAL", "BASE")]
    for model, tag in (("fm", "NLOAD"), ("fc", "NLOADC")):
        script += ["SELECT duckvep_fault_arm(0);", loads[model],
                   "SELECT '%s ' || duckvep_fault_allocations();" % tag]
    script += ["SELECT duckvep_fault_arm(0);"]
    for name, query in QUERIES.items():
        script += [query, "SELECT 'NQUERY %s ' || duckvep_fault_allocations();" % name,
                   "SELECT duckvep_fault_arm(0);"]
    script += ["SELECT duckvep_model_drop('fm');", "SELECT duckvep_model_drop('fc');",
               BASELINE_SQL.replace("TOTAL", "USED")]
    rc, out = run(args.extension, "\n".join(script), args.duckdb)
    if rc not in (0, 1) or out.count("LOADED true") != 2:
        print(out[-4000:])
        raise SystemExit("count pass failed (exit %d)" % rc)
    base = int(value(out, "BASE "))
    n_load = {"load": int(value(out, "NLOAD ")), "load_circular": int(value(out, "NLOADC "))}
    golden = {name: value(out, "RESULT %s " % name) for name in QUERIES}
    for name, text in golden.items():
        if int(text.split()[0]) == 0:
            raise SystemExit("golden %s is empty; the fixture does not exercise it" % name)
    n_query = {name: int(value(out, "NQUERY %s " % name)) for name in QUERIES}
    used = int(value(out, "USED ")) - base
    rc, out = run(args.extension, "\n".join([FIXTURE,
        "SELECT duckvep_worker_limits_set(6, 134217728, 67108864, 0);", QUERIES["hgvs"],
        BASELINE_SQL.replace("TOTAL", "UNKNOWN")]), args.duckdb)
    unknown = int(value(out, "UNKNOWN ")) - base
    allowed = {0, used, unknown}
    print("allocations: %s %s; baseline=%d retained-worker bytes=%s" % (
        " ".join("%s=%d" % kv for kv in n_load.items()), " ".join("%s=%d" % kv for kv in n_query.items()), base, sorted(allowed)))
    for name, text in golden.items():
        print("  golden %s: %s" % (name, text))

    plan = []
    for name, count in n_load.items():
        plan += [(name, i) for i in range(1, (min(count, args.limit) if args.limit else count) + 1)]
    for name, count in n_query.items():
        plan += [(name, i) for i in range(1, (min(count, args.limit) if args.limit else count) + 1)]
    failures, sites, outcomes = [], collections.Counter(), collections.Counter()
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(iteration, args, phase, index, loads, base, allowed, golden)
                   for phase, index in plan]
        for future in futures:
            problems, info, out = future.result()
            failures += problems
            if info:
                sites[info[0]] += 1
                outcomes[info[1]] += 1
            if problems:
                failures.append(out[-1500:])
    print("failed %d allocations (%s); %d distinct call sites" % (
        len(plan), ", ".join("%s=%d" % kv for kv in sorted(outcomes.items())), len(sites)))
    for site, count in sorted(sites.items(), key=lambda kv: str(kv[0])):
        print("  site %s x%d" % (site, count))
    # The new lifted-execution sites must actually have been failed.
    for needed in ("duckvep_lift.c", "duckvep_scalar_lift_resolve", "duckvep_scalar_lift_reserve"):
        if not args.limit and not any(needed in "%s %s" % (site[0], site[1]) for site in sites if site):
            failures.append("no failed allocation at %s" % needed)
    if failures:
        print("FAILURES (%d):" % len(failures))
        for item in failures[:40]:
            print("  " + item)
        sys.exit(1)
    print("fault injection: all %d failures were clean" % len(plan))


if __name__ == "__main__":
    main()
