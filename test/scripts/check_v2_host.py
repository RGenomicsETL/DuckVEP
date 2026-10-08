#!/usr/bin/env python3
"""Static and binary gates for the DuckDB C API v2 host (host_v2).

  symbols  the binary exports exactly duckvep_init_c_api_v2 and imports no duckdb_* symbol
           (a v2 extension reaches DuckDB only through the function table it is handed, so
           an imported DuckDB symbol would mean a v1-style or unstable linkage); mirrors
           `make test-extension-symbols` for the v1 binary
  footer   the appended metadata is ABI C_STRUCT (stable), extension API v2.0.0
  static   the host_v2 sources use only stable, non-deprecated duckdb_v2_* functions from the pinned
           headers, defines the unstable/deprecated opt-ins as 0, never includes the v1
           headers, and keeps src/core and src/kernel free of DuckDB
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("DUCKVEP_V2_EXTENSION", ROOT / "build/release_v2/duckvep.duckdb_extension"))
PIN = json.loads((ROOT / "duckvep-package.json").read_text())["v2_host"]


def fail(message):
    raise SystemExit("FAIL: " + message)


def symbols():
    if not BINARY.exists():
        fail(f"{BINARY} missing; run make release_v2")
    if sys.platform == "darwin":
        defined = subprocess.check_output(["nm", "-gU", str(BINARY)], text=True)
        exported = {line.split()[-1].lstrip("_") for line in defined.splitlines() if " T " in line}
        undefined = subprocess.check_output(["nm", "-gu", str(BINARY)], text=True)
        imports = {line.split()[-1].lstrip("_") for line in undefined.splitlines()}
    else:
        defined = subprocess.check_output(["nm", "-D", "--defined-only", str(BINARY)], text=True)
        exported = {line.split()[2] for line in defined.splitlines()
                    if len(line.split()) == 3 and line.split()[1] in "TDB" and line.split()[2] not in ("_init", "_fini")}
        undefined = subprocess.check_output(["nm", "-D", "-u", str(BINARY)], text=True)
        imports = {line.split()[-1] for line in undefined.splitlines() if line.strip()}
    if exported != {"duckvep_init_c_api_v2"}:
        fail(f"unexpected exports: {sorted(exported)}")
    bad = sorted(name for name in imports if name.startswith("duckdb"))
    if bad:
        fail(f"unexpected DuckDB imports: {bad}")
    print("symbols: exports exactly duckvep_init_c_api_v2; no duckdb_* imports")


def footer():
    data = BINARY.read_bytes()
    # append_extension_metadata.py writes FIELD8 first and FIELD1 (signature) last, then 256 signature bytes.
    fields = [data[-512 + 32 * i:-512 + 32 * (i + 1)].rstrip(b"\0").decode() for i in range(8)][::-1]
    signature, platform, engine, version, abi = fields[:5]
    if signature != "4":
        fail(f"footer signature {signature!r}")
    if abi != "C_STRUCT":
        fail(f"footer ABI type {abi!r}; expected C_STRUCT (C_STRUCT_UNSTABLE pins one DuckDB build)")
    if engine != PIN["extension_api_version"]:
        fail(f"footer API version {engine!r}; expected {PIN['extension_api_version']}")
    print(f"footer: abi={abi} api={engine} platform={platform} extension={version}")


def static():
    files = sorted((ROOT / "host_v2").glob("*.[ch]"))
    source = "\n".join(path.read_text() for path in files)
    header = (ROOT / PIN["headers"] / "duckdb_v2.h").read_text()
    for name, expected in PIN["duckdb_sdk_sha256"].items():
        if hashlib.sha256((ROOT / PIN["headers"] / name).read_bytes()).hexdigest() != expected:
            fail(f"{name} does not match its pinned checksum")
    for needle in ('#include "duckdb_extension_v2.h"', "#define DUCKDB_V2_API_ALLOW_UNSTABLE 0",
                   "#define DUCKDB_V2_API_ALLOW_DEPRECATED 0"):
        if needle not in source:
            fail(f"host_v2 sources lack {needle}")
    if re.search(r'#include\s+[<"]duckdb(_extension)?\.h[>"]', source) or "DUCKDB_EXTENSION_API_VERSION_UNSTABLE" in source:
        fail("host_v2 sources reference the v1 headers or the v1 unstable opt-in")
    code = re.sub(r"/\*.*?\*/|//[^\n]*", "", source, flags=re.S)
    called = set(re.findall(r"\b(duckdb_[a-z0-9_]+)\s*\(", code))
    stray = sorted(name for name in called if not name.startswith("duckdb_v2_"))
    if stray:
        fail(f"non-v2 DuckDB calls in host_v2: {stray}")
    for name in sorted(called):
        declaration = re.search(rf"DUCKDB_C_API[^;{{}}]*?\b{name}\s*\(", header, re.S)
        if not declaration:
            fail(f"{name} is not declared in the pinned header")
        doc = header[header.rfind("/*!", 0, declaration.start()):declaration.start()]
        history = doc.split("history:")[-1].split("@param")[0]
        if "- stable: v2.0.0" not in history or re.search(r"unstable|deprecated", history, re.I):
            fail(f"{name} is not stable v2.0.0 in the pinned header: {history.strip()}")
    # duckdb_ext_api is the function table; no other DuckDB name may appear in host-neutral code.
    for directory in ("src/core", "src/kernel"):
        for path in (ROOT / directory).rglob("*"):
            if path.suffix in (".c", ".h", ".inc") and re.search(r'duckdb(_extension)?\.h|\bduckdb_[a-z]+_?[a-z_]*\s*\(', path.read_text()):
                fail(f"{path.relative_to(ROOT)} mentions the DuckDB API; it must stay host-neutral")
    print(f"static: {len(called)} duckdb_v2_* functions, all stable v2.0.0; unstable/deprecated opt-ins off; "
          "src/core and src/kernel host-neutral")


if __name__ == "__main__":
    {"symbols": symbols, "footer": footer, "static": static}[sys.argv[1]]()
