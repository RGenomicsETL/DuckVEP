#!/usr/bin/env python3
"""Build a DuckDB CLI at the pinned v2 SDK revision into .deps-v2/duckdb-build/duckdb.

The preview v2 extension ABI is not frozen across snapshots, so the v2 host must run on a
DuckDB built at the revision its headers came from (duckvep-package.json). Prints the path.
"""
from pathlib import Path
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REVISION = json.loads((ROOT / "duckvep-package.json").read_text())["v2_host"]["duckdb_sdk_revision"]
SOURCE = ROOT / ".deps-v2/duckdb"
BUILD = ROOT / ".deps-v2/duckdb-build"
BINARY = BUILD / "duckdb"


def run(*command, cwd=None):
    subprocess.run(command, check=True, cwd=cwd)


def main():
    jobs = sys.argv[1] if len(sys.argv) > 1 else str(os.cpu_count() or 2)
    if BINARY.exists():
        version = subprocess.check_output([str(BINARY), "--version"], text=True)
        if REVISION[:10] not in version:
            raise SystemExit(f"cached DuckDB {version.strip()} does not match {REVISION}")
        print(BINARY)
        return
    if not (SOURCE / ".git").exists():
        SOURCE.mkdir(parents=True, exist_ok=True)
        run("git", "init", "-q", cwd=SOURCE)
        run("git", "remote", "add", "origin", "https://github.com/duckdb/duckdb.git", cwd=SOURCE)
    run("git", "fetch", "-q", "--depth", "1", "origin", REVISION, cwd=SOURCE)
    run("git", "checkout", "-q", "--detach", REVISION, cwd=SOURCE)
    # json provides to_json, which the equality tests serialize rows with.
    run("cmake", "-S", str(SOURCE), "-B", str(BUILD), "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_UNITTESTS=OFF", "-DBUILD_SHELL=ON", "-DBUILD_EXTENSIONS=json")
    run("cmake", "--build", str(BUILD), "--target", "shell", f"-j{jobs}")
    print(BINARY)


if __name__ == "__main__":
    main()
