#!/usr/bin/env python3
"""Verify (default) or refresh the pinned DuckDB C API v2 headers in duckdb_capi_v2/."""
from pathlib import Path
import hashlib
import json
import sys
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
PIN = json.loads((ROOT / "duckvep-package.json").read_text())["v2_host"]
REVISION = PIN["duckdb_sdk_revision"]
HEADERS = PIN["duckdb_sdk_sha256"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    refresh = "--refresh" in sys.argv[1:]
    target = ROOT / PIN["headers"]
    target.mkdir(exist_ok=True)
    failed = False
    for name, expected in HEADERS.items():
        path = target / name
        if path.exists() and digest(path) == expected:
            continue
        if not refresh:
            print(f"SDK checksum mismatch or missing: {path}", file=sys.stderr)
            failed = True
            continue
        url = f"https://raw.githubusercontent.com/duckdb/duckdb/{REVISION}/src/include/{name}"
        with urlopen(url, timeout=60) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != expected:
            raise RuntimeError(f"SDK checksum mismatch for {name} at {REVISION}")
        path.write_bytes(data)
    revision = target / "REVISION"
    if not revision.exists() or revision.read_text().strip() != REVISION:
        if refresh:
            revision.write_text(REVISION + "\n")
        else:
            print("REVISION does not match duckvep-package.json", file=sys.stderr)
            failed = True
    if failed:
        raise SystemExit(1)
    print(f"DuckDB v2 SDK {REVISION}: headers verified")


if __name__ == "__main__":
    main()
