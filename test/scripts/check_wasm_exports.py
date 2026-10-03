#!/usr/bin/env python3
"""Check that a webR/Emscripten side module of DuckVEP keeps its bundled libraries private.

Another side module in the same process (Rduckhts) carries its own HTSlib, zlib and cgranges.
The module must export only its DuckDB entry point, and must not import any bundled symbol
through the GOT, where the dynamic loader could bind it to the other module's definition.

Usage: check_wasm_exports.py <module.so|module.wasm|module.duckdb_extension>
"""
import re
import sys

ALLOWED_EXPORTS = re.compile(r"^(duckvep_init_c_api|__wasm_call_ctors|__wasm_apply_data_relocs|dynCall_\w+)$")
BUNDLED = re.compile(
    r"^(hts_|bgzf_|fai_|faidx_|bcf_|vcf_|sam_|bam_|cram_|tbx_|seq_nt16|hfile_|hopen|hclose|hts|kputs|ks_|kh_"
    r"|cr_|z_|zc|inflate|deflate|crc32|adler32|_length_code|_dist_code|duckvep_)")


def leb(data, pos):
    value = shift = 0
    while True:
        byte = data[pos]
        pos += 1
        value |= (byte & 0x7F) << shift
        shift += 7
        if not byte & 0x80:
            return value, pos


def name(data, pos):
    size, pos = leb(data, pos)
    return data[pos:pos + size].decode("utf-8"), pos + size


def limits(data, pos):
    flags, pos = leb(data, pos)
    _, pos = leb(data, pos)
    if flags & 1:
        _, pos = leb(data, pos)
    return pos


def sections(data):
    if data[:4] != b"\0asm":
        sys.exit("not a WebAssembly module")
    pos = 8
    while pos < len(data):
        section = data[pos]
        size, pos = leb(data, pos + 1)
        yield section, pos, pos + size
        pos += size


def imports(data, pos):
    count, pos = leb(data, pos)
    for _ in range(count):
        module, pos = name(data, pos)
        field, pos = name(data, pos)
        kind = data[pos]
        pos += 1
        if kind == 0:  # function: type index
            _, pos = leb(data, pos)
        elif kind == 1:  # table: reference type, limits
            pos = limits(data, pos + 1)
        elif kind == 2:  # memory: limits
            pos = limits(data, pos)
        elif kind == 3:  # global: value type, mutability
            pos += 2
        elif kind == 4:  # tag: attribute, type index
            _, pos = leb(data, pos + 1)
        else:
            sys.exit(f"unknown import kind {kind}")
        yield module, field


def exports(data, pos):
    count, pos = leb(data, pos)
    for _ in range(count):
        field, pos = name(data, pos)
        _, pos = leb(data, pos + 1)
        yield field


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    data = open(sys.argv[1], "rb").read()
    bad_exports, bad_imports, seen = [], [], False
    for section, start, _ in sections(data):
        if section == 2:
            bad_imports = sorted(f"{module}.{field}" for module, field in imports(data, start)
                                 if module.startswith("GOT.") and BUNDLED.match(field))
        elif section == 7:
            found = list(exports(data, start))
            seen = "duckvep_init_c_api" in found
            bad_exports = sorted(field for field in found if not ALLOWED_EXPORTS.match(field))
    if not seen:
        sys.exit("duckvep_init_c_api is not exported")
    if bad_exports or bad_imports:
        for field in bad_exports[:20]:
            print(f"Unexpected export: {field}")
        for field in bad_imports[:20]:
            print(f"Interposable import of a bundled symbol: {field}")
        sys.exit(f"{len(bad_exports)} unexpected exports, {len(bad_imports)} interposable imports")
    print("wasm side module exports only duckvep_init_c_api")


if __name__ == "__main__":
    main()
