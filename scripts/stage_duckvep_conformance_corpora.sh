#!/usr/bin/env bash
# Compatibility entry point for a separately supplied registry staging runner.
#
# Usage:
#   scripts/stage_duckvep_conformance_corpora.sh [OUTPUT_DIR] [CORPUS]
#
# Set DUCKVEP_CORPUS_STAGER to the registry's R staging script.
# OUTPUT_DIR defaults to the registry paths below $DUCKHTS_CACHE_DIR. CORPUS is
# one of: all (default), hprc-african4-chr22, sniffles2-chr22, dbvar-chr22.
set -euo pipefail

if [[ $# -gt 2 ]]; then
  echo "usage: $0 [OUTPUT_DIR] [CORPUS]" >&2
  exit 2
fi
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGER="${DUCKVEP_CORPUS_STAGER:-$ROOT_DIR/r/duckhtsbench/scripts/stage_duckvep_conformance_corpora.R}"
if [[ ! -f "$STAGER" ]]; then
  echo "Set DUCKVEP_CORPUS_STAGER to an existing registry R staging script; this checkout does not supply one." >&2
  exit 2
fi
OUTPUT_DIR="${1:-}"
CORPUS="${2:-all}"
ARGS=(--corpus "$CORPUS")
if [[ -n "$OUTPUT_DIR" ]]; then
  ARGS+=(--output-dir "$OUTPUT_DIR")
fi
exec Rscript "$STAGER" "${ARGS[@]}"
