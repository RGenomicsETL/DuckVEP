#!/usr/bin/env python3
"""Check docs/functions.md against the built extension.

1. Every public function registered by the extension (duckdb_functions(),
   names starting with duckvep_, so the _duckvep_* internals are excluded)
   has a `### name` section, and no section names a function that is not
   registered.
2. Every section contains at least one ```sql example that calls its function.
3. Every ```sql example runs, in document order on one connection, after the
   README fixture model (test/data/duckvep/readme.sql) and the small Ensembl
   core fixture (test/data/duckvep/functions_docs_ensembl.sql). Blocks fenced
   as ```sql no-run are illustrative and skipped.

Usage: scripts/check-function-docs.py [extension_file]
Run from the repository root; `make check-function-docs` does that with the
DuckDB Python module of the configure venv. DUCKVEP_EXTENSION_FILE overrides
the default build/release/duckvep.duckdb_extension.
"""
import os
import re
import sys

import duckdb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs", "functions.md")
FIXTURES = [
    os.path.join(ROOT, "test", "data", "duckvep", "readme.sql"),
    os.path.join(ROOT, "test", "data", "duckvep", "functions_docs_ensembl.sql"),
]


def main():
    extension = (sys.argv[1] if len(sys.argv) > 1
                 else os.environ.get("DUCKVEP_EXTENSION_FILE")
                 or os.path.join(ROOT, "build", "release", "duckvep.duckdb_extension"))
    if not os.path.isfile(extension):
        sys.exit("check-function-docs: %s not found; run `make release` first" % extension)
    os.chdir(ROOT)  # fixtures read test/data/... relative to the repository root
    text = open(DOCS, encoding="utf-8").read()

    con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    con.execute("LOAD '%s'" % extension.replace("'", "''"))

    registered = {r[0] for r in con.execute(
        "SELECT DISTINCT function_name FROM duckdb_functions() "
        "WHERE starts_with(function_name, 'duckvep_')").fetchall()}
    documented = set(re.findall(r"^### (duckvep_[a-z0-9_]+)\s*$", text, re.M))
    failures = []
    for name in sorted(registered - documented):
        failures.append("function not documented: " + name)
    for name in sorted(documented - registered):
        failures.append("documented function is not registered: " + name)

    # Split into sections, then collect the runnable sql blocks in order.
    sections = re.split(r"^### ", text, flags=re.M)[1:]
    for section in sections:
        name = section.split("\n", 1)[0].strip()
        blocks = re.findall(r"^```sql[ \t]*\n(.*?)^```", section, re.M | re.S)
        if not any(name in b for b in blocks):
            failures.append("no runnable example calls " + name)

    for path in FIXTURES:
        with open(path, encoding="utf-8") as handle:
            con.execute(handle.read())

    ran = 0
    for match in re.finditer(r"^```sql[ \t]*\n(.*?)^```", text, re.M | re.S):
        sql = match.group(1)
        line = text.count("\n", 0, match.start()) + 2
        try:
            con.execute(sql)
            con.fetchall() if con.description else None
            ran += 1
        except Exception as error:  # report every failing example
            failures.append("example at docs/functions.md:%d failed: %s\n%s"
                            % (line, str(error).splitlines()[0], sql))
    skipped = len(re.findall(r"^```sql no-run", text, re.M))

    if failures:
        print("\n".join(failures), file=sys.stderr)
        sys.exit("check-function-docs: FAILED (%d problem(s))" % len(failures))
    print("check-function-docs: %d public functions documented, %d examples ran, "
          "%d illustrative example(s) skipped" % (len(registered), ran, skipped))


if __name__ == "__main__":
    main()
