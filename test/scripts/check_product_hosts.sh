#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
: "${DUCKVEP_V2_DUCKDB:?Set DUCKVEP_V2_DUCKDB to the pinned CLI}"
exec bash "$root/test/scripts/check_products_typed.sh" \
  "${DUCKVEP_V1_DUCKDB:-duckdb}" "$DUCKVEP_V2_DUCKDB" \
  "${DUCKVEP_PRODUCTS_OUTPUT:-$root/build/products-typed}"
