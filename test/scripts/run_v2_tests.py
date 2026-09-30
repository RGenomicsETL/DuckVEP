#!/usr/bin/env python3
"""Run the DuckVEP v2-host tests.

The v2 extension is built against a pinned preview SDK whose ABI is not frozen,
so it must run on a DuckDB built at that revision (duckvep-package.json,
v2_host.duckdb_sdk_revision). Other 2.0 preview snapshots, including the PyPI
`duckdb>=2.0.0.dev0` wheels, fail to load it (segfault or garbage function
table), so this runner refuses an engine whose source id differs.

  DUCKVEP_V2_DUCKDB   DuckDB CLI built at the pinned revision (scripts/stage-v2-duckdb.py)
  --v1-extension      the v1 build (default build/release/duckvep.duckdb_extension); used with
                      the `duckdb` Python package (a stable DuckDB >= 1.4) to re-verify the
                      golden file live. Skipped when absent.

Checks: 1. LOAD (read-only, twice, no DDL, no database change); 2. v2-only SQL
assertions; 3. every equality case against test/sql_v2/equality_golden.json, the
recorded v1-host results; 4. with --v1-extension, the golden against the live v1 host.
Use --record to regenerate the golden file from the v1 host.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "test/sql_v2/equality_cases.sql"
NATIVE = ROOT / "test/sql_v2/v2_native.sql"
GOLDEN = ROOT / "test/sql_v2/equality_golden.json"
PIN = json.loads((ROOT / "duckvep-package.json").read_text())["v2_host"]
MESSAGE = re.compile(r"duckvep_\w+: [^\n]*")


def quote(path):
    return str(path).replace("'", "''")


def parse_cases():
    cases, name, lines = [], None, []
    for line in CASES.read_text().splitlines():
        if line.startswith("-- case:"):
            if name:
                cases.append((name, " ".join(lines).strip()))
            name, lines = line[len("-- case:"):].strip(), []
        elif name and not line.startswith("--") and line.strip():
            lines.append(line.strip())
    cases.append((name, " ".join(lines).strip()))
    return cases


def wrapped(query):
    return f"SELECT CAST(to_json(t) AS VARCHAR) FROM ({query}) t"


def outcome_error(text):
    match = MESSAGE.search(text)
    return {"error": match.group(0).strip() if match else text.strip().splitlines()[0]}


def check_case_shape(outcomes):
    """A case named *error* must fail with a duckvep message; no other case may fail.
    This keeps a typo in a case from being recorded as an 'equal' error."""
    for name, outcome in outcomes.items():
        failed = "error" in outcome
        if name.endswith("error") != failed or (failed and not MESSAGE.match(outcome["error"])):
            raise SystemExit(f"case {name!r} has an unexpected outcome: {outcome}")


class V2Host:
    def __init__(self, cli, extension):
        self.cli = str(cli)
        self.extension = Path(extension).resolve()

    # -column: the pinned snapshot's CLI swallows statement errors (exit 0) in its streaming
    # output modes (-list, -csv, -json, -line), but reports them in -column, -table and -box.
    # -column pads lines, so callers strip trailing whitespace.
    def run(self, sql, database=":memory:", readonly=False, check=True):
        command = [self.cli, "-unsigned", "-no-init", "-batch", "-bail", "-noheader", "-column"]
        if readonly:
            command.append("-readonly")
        command.append(str(database))
        preamble = "SET autoinstall_known_extensions=false; SET autoload_known_extensions=false;\n"
        result = subprocess.run(command, input=preamble + sql, text=True, capture_output=True,
                                timeout=600, cwd=ROOT)
        if check and result.returncode != 0:
            raise AssertionError(f"{result.stderr}\n{result.stdout}")
        return result

    def load(self):
        return f"LOAD '{quote(self.extension)}';\n"

    def case(self, query):
        result = self.run(self.load() + wrapped(query) + ";", check=False)
        if result.returncode != 0:
            return outcome_error(result.stderr)
        return {"rows": [line.rstrip() for line in result.stdout.splitlines()]}


def v1_outcomes(extension, cases):
    import duckdb
    connection = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    connection.execute(f"LOAD '{quote(Path(extension).resolve())}'")
    outcomes = {}
    for name, query in cases:
        try:
            outcomes[name] = {"rows": [row[0] for row in connection.execute(wrapped(query)).fetchall()]}
        except duckdb.Error as error:
            outcomes[name] = outcome_error(str(error))
    return outcomes


def check_engine(host):
    source = host.run("SELECT source_id FROM pragma_version();").stdout.strip()
    revision = PIN["duckdb_sdk_revision"]
    if not revision.startswith(source) or len(source) < 10:
        raise SystemExit(f"DuckDB engine {source!r} is not the pinned revision {revision}; "
                        "the preview v2 ABI is not frozen across snapshots. "
                        "Build one with scripts/stage-v2-duckdb.py.")
    print(f"engine: DuckDB {source} (pinned revision)")


def test_load(host):
    with tempfile.TemporaryDirectory(prefix="duckvep-v2-") as directory:
        database = Path(directory) / "sentinel.duckdb"
        host.run("CREATE TABLE sentinel AS SELECT 42 AS value; CHECKPOINT;", database)
        digest = hashlib.sha256(database.read_bytes()).hexdigest()
        absent = ("SELECT CASE WHEN count(*) = 0 THEN true ELSE error('persisted duckvep catalog entry') END "
                  "FROM duckdb_functions() WHERE function_name LIKE '%duckvep%' AND function_type = 'macro';\n")
        for readonly in (False, True):
            host.run(absent, database, readonly)
            # Twice: a repeated LOAD is a no-op, not an error.
            host.run(host.load() + host.load() + "SELECT * FROM sentinel;\n"
                     "SELECT count(*) FROM duckvep_so_terms();\n"
                     "SELECT duckvep_allele_geometry(1, 'A', 'C');\n"
                     "SELECT duckvep_breakend_geometry('A]1:2]');\n", database, readonly)
            host.run(absent, database, readonly)
            assert hashlib.sha256(database.read_bytes()).hexdigest() == digest, "LOAD changed the database file"
    print("load: read-only and writable primary, repeated LOAD, no DDL, database bytes unchanged")


def test_native(host):
    host.run(host.load() + NATIVE.read_text())
    print("native: v2-only assertions passed")


def test_equality(host, golden):
    outcomes = {name: host.case(query) for name, query in parse_cases()}
    check_case_shape(outcomes)
    bad = [name for name, outcome in outcomes.items() if golden.get(name) != outcome]
    missing = sorted(set(golden) - set(outcomes))
    if bad or missing:
        for name in bad:
            print(f"MISMATCH {name}\n  v1: {golden.get(name)}\n  v2: {outcomes[name]}", file=sys.stderr)
        raise SystemExit(f"{len(bad)} case(s) differ from the v1 host; golden-only cases: {missing}")
    errors = sum("error" in outcome for outcome in outcomes.values())
    print(f"equality: {len(outcomes)} cases identical to the v1 host ({errors} error cases)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v2-extension", default=ROOT / "build/release_v2/duckvep.duckdb_extension")
    parser.add_argument("--v2-duckdb", default=os.environ.get("DUCKVEP_V2_DUCKDB"))
    parser.add_argument("--v1-extension", default=ROOT / "build/release/duckvep.duckdb_extension")
    parser.add_argument("--record", action="store_true", help="write the golden file from the v1 host and exit")
    args = parser.parse_args()
    if args.record:
        outcomes = v1_outcomes(args.v1_extension, parse_cases())
        check_case_shape(outcomes)
        GOLDEN.write_text(json.dumps(outcomes, indent=1, sort_keys=True) + "\n")
        print(f"recorded {len(outcomes)} cases from the v1 host into {GOLDEN.relative_to(ROOT)}")
        return
    if not args.v2_duckdb:
        raise SystemExit("set DUCKVEP_V2_DUCKDB or --v2-duckdb to a DuckDB CLI built at the pinned revision")
    host = V2Host(args.v2_duckdb, args.v2_extension)
    check_engine(host)
    golden = json.loads(GOLDEN.read_text())
    test_load(host)
    test_native(host)
    test_equality(host, golden)
    if Path(args.v1_extension).exists():
        live = v1_outcomes(args.v1_extension, parse_cases())
        if live != golden:
            diff = [name for name in live if live[name] != golden.get(name)]
            raise SystemExit(f"golden file is stale against the live v1 host: {diff}")
        print("golden: matches the live v1 host")
    else:
        print("golden: live v1 check skipped (no v1 extension build)")


if __name__ == "__main__":
    main()
