#!/usr/bin/env python3
"""Call the pinned Mutalyzer normalizer on one description per JSONL record."""

import argparse
import importlib.metadata
import json
import os
import re
import socket
import sys


def package_name(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def read_pins(path):
    pins = {}
    with open(path, encoding="utf-8") as lock_file:
        for line in lock_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            name, separator, version = line.partition("==")
            if not separator or not name or not version:
                raise ValueError(f"invalid lock entry: {line!r}")
            key = package_name(name)
            if key in pins:
                raise ValueError(f"duplicate lock entry: {name}")
            pins[key] = version
    return pins


def verify_environment(pins):
    installed = {
        package_name(distribution.metadata["Name"]): distribution.version
        for distribution in importlib.metadata.distributions()
        if distribution.metadata.get("Name")
    }
    mismatches = [
        f"{name}: expected {version}, found {installed.get(name, 'missing')}"
        for name, version in pins.items()
        if installed.get(name) != version
    ]
    extras = sorted(set(installed) - set(pins))
    if mismatches or extras:
        details = mismatches
        if extras:
            details.append(f"unlocked distributions: {', '.join(extras)}")
        raise RuntimeError(
            "Mutalyzer environment does not match requirements.lock: "
            + "; ".join(details)
        )


def disable_network():
    attempts = []

    def blocked(*args, **kwargs):
        address = args[-1] if args else kwargs
        attempts.append(repr(address))
        raise OSError("network access disabled by compound HGVS oracle")

    socket.socket.connect = blocked
    socket.socket.connect_ex = blocked
    socket.socket.sendto = blocked
    if hasattr(socket.socket, "sendmsg"):
        socket.socket.sendmsg = blocked
    socket.create_connection = blocked
    socket.getaddrinfo = blocked
    socket.gethostbyname = blocked
    socket.gethostbyname_ex = blocked
    return attempts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--lock", required=True)
    args = parser.parse_args()

    pins = read_pins(args.lock)
    verify_environment(pins)
    if sys.implementation.name != "cpython" or sys.version.split()[0] != "3.13.12":
        raise RuntimeError("oracle requires CPython 3.13.12")
    attempts = disable_network()

    import mutalyzer_retriever.configuration as configuration
    from mutalyzer.normalizer import normalize

    expected_cache = os.environ.get("MUTALYZER_EXPECTED_CACHE_DIR")
    actual_cache = os.path.realpath(configuration.cache_dir() or "")
    if not expected_cache or actual_cache != os.path.realpath(expected_cache):
        raise RuntimeError(
            "Mutalyzer cache does not match MUTALYZER_EXPECTED_CACHE_DIR"
        )
    if configuration.cache_url() is not None:
        raise RuntimeError(
            "Mutalyzer API cache is enabled; offline oracle requires file cache only"
        )
    if configuration.cache_add() is not False:
        raise RuntimeError("Mutalyzer file-cache writes must be disabled")

    seen = set()
    with open(args.input, encoding="utf-8") as source, open(
        args.output, "w", encoding="utf-8"
    ) as target:
        for line_number, line in enumerate(source, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            if not isinstance(record, dict) or set(record) != {"case_id", "description"}:
                raise ValueError(f"input line {line_number}: expected case_id and description")
            case_id = record["case_id"]
            description = record["description"]
            if not isinstance(case_id, str) or not case_id or case_id in seen:
                raise ValueError(f"input line {line_number}: case_id must be unique and non-empty")
            if not isinstance(description, str) or not description:
                raise ValueError(f"input line {line_number}: description must be non-empty")
            seen.add(case_id)

            result = normalize(description)
            if attempts:
                raise RuntimeError(
                    "Mutalyzer attempted network access: " + ", ".join(attempts)
                )
            protein = result.get("protein") or {}
            errors = result.get("errors") or []
            output = {
                "schema": "mutalyzer-oracle-result/v1",
                "case_id": case_id,
                "input_description": description,
                "normalized_description": result.get("normalized_description"),
                "protein_description": protein.get("description"),
                "protein_reference": protein.get("reference"),
                "protein_predicted": protein.get("predicted"),
                "errors": [error.get("code", "?") for error in errors],
                "mutalyzer_version": pins["mutalyzer"],
                "retriever_version": pins["mutalyzer-retriever"],
                "python_version": sys.version.split()[0],
            }
            target.write(json.dumps(output, ensure_ascii=False, separators=(",", ":")) + "\n")

    if not seen:
        raise ValueError("input contains no descriptions")


if __name__ == "__main__":
    main()
