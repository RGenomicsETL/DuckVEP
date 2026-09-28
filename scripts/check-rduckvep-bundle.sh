#!/usr/bin/env bash
# Fails if the bundled extension sources committed at
# r/Rduckvep/inst/duckvep_extension differ from what bootstrap.R would
# regenerate from the tracked sources (src, cmake, duckdb_capi,
# third_party/htslib, third_party/cgranges, CMakeLists.txt). r-universe
# cannot run bootstrap.R itself, so the bundle is committed and must be kept
# in sync by hand.
#
# Usage:
#   scripts/check-rduckvep-bundle.sh
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUNDLE_DIR="r/Rduckvep/inst/duckvep_extension"

cd "$ROOT_DIR"
Rscript r/Rduckvep/bootstrap.R "$ROOT_DIR"

CHANGES="$(git status --porcelain -- "$BUNDLE_DIR")"
if [[ -n "$CHANGES" ]]; then
  echo "error: $BUNDLE_DIR is out of date with its tracked sources." >&2
  echo "The following paths differ from what bootstrap.R regenerates:" >&2
  echo "$CHANGES" >&2
  echo >&2
  echo "Fix: run 'Rscript r/Rduckvep/bootstrap.R .' and commit the regenerated bundle." >&2
  exit 1
fi

echo "$BUNDLE_DIR is in sync with its tracked sources."
