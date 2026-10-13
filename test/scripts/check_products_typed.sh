#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
v1="${1:-${DUCKVEP_V1_DUCKDB:-duckdb}}"
v2="${2:-${DUCKVEP_V2_DUCKDB:?Supply the pinned v2 CLI as argument 2 or DUCKVEP_V2_DUCKDB}}"
output="${3:-$root/build/products-typed}"
mkdir -p "$output"
for host in v1 v2; do
  cli="$v1"
  ext="$root/build/release/duckvep.duckdb_extension"
  if [[ "$host" == v2 ]]; then
    cli="$v2"
    ext="$root/build/release_v2/duckvep.duckdb_extension"
  fi
  {
    printf "LOAD '%s';\n" "$ext"
    awk '{print}' "$root/test/sql_v2/products_typed_cases.sql"
  } | "$cli" -unsigned -no-init -batch -bail -csv >"$output/$host-full.csv"
  [[ "$(wc -l <"$output/$host-full.csv")" == 9 ]]
  {
    printf "LOAD '%s';\n" "$ext"
    awk '{print}' "$root/test/sql_v2/products_translation_cases.sql"
  } | "$cli" -unsigned -no-init -batch -bail -csv >"$output/$host-translation.csv"
  [[ "$(wc -l <"$output/$host-translation.csv")" == 6 ]]
  {
    printf "LOAD '%s';\n" "$ext"
    awk '{print}' "$root/test/sql_v2/products_257.sql"
  } | "$cli" -unsigned -no-init -batch -bail -csv >"$output/$host-257.csv"
  grep -qx '257,257,257,257,1' "$output/$host-257.csv"
  {
    printf "LOAD '%s';\n" "$ext"
    awk '{print}' "$root/test/sql_v2/products_uint64_edge.sql"
  } | "$cli" -unsigned -no-init -batch -bail -csv >"$output/$host-uint64.csv"
  grep -qx '18446744073709551610,18446744073709551610,18446744073709551607' \
    "$output/$host-uint64.csv"
  {
    printf "LOAD '%s';\n" "$ext"
    awk '{print}' "$root/test/sql_v2/products_validation.sql"
  } | "$cli" -unsigned -no-init >"$output/$host-validation.out" 2>"$output/$host-validation.err" || true
  [[ "$(grep -c 'Invalid Input Error:' "$output/$host-validation.err")" == 10 ]]
  grep -q valid_after_errors "$output/$host-validation.out"
  grep -q two_sources "$output/$host-validation.out"
done
cmp "$output/v1-full.csv" "$output/v2-full.csv"
cmp "$output/v1-translation.csv" "$output/v2-translation.csv"
cmp "$output/v1-257.csv" "$output/v2-257.csv"
cmp "$output/v1-uint64.csv" "$output/v2-uint64.csv"
cmp "$output/v1-validation.err" "$output/v2-validation.err"
cmp "$output/v1-validation.out" "$output/v2-validation.out"
sha256sum "$output"/{v1,v2}-{full,translation,257,uint64}.csv
