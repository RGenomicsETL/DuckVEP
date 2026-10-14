#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
out=${1:?Provide a writable evidence directory}
mkdir -p "$out"
out=$(realpath "$out")
cli_v1=${DUCKVEP_V1_DUCKDB:?Set DUCKVEP_V1_DUCKDB to the pinned v1 CLI}
cli_v2=${DUCKVEP_V2_DUCKDB:?Set DUCKVEP_V2_DUCKDB to the pinned v2 CLI}
awk '/^CREATE SCHEMA duckvep_core;$/ {start=1} start && /^statement ok$/ {next} start && /^----$/ {exit} start && /^# / {next} start && /^query / {exit} start {print}' test/sql/duckvep_ensembl.test > "$out/seed.sql"
for version in v1 v2; do
  cli=$cli_v1
  ext=build/release/duckvep.duckdb_extension
  if [[ $version == v2 ]]; then cli=$cli_v2; ext=build/release_v2/duckvep.duckdb_extension; fi
  mkdir -p "$out/$version"
  cp "$ext" "$out/$version/duckvep.duckdb_extension"
  sha256sum "$out/$version/duckvep.duckdb_extension" > "$out/$version.sha256"
  { printf "LOAD '%s/%s/duckvep.duckdb_extension';\n" "$out" "$version"; printf ".read %s/seed.sql\n" "$out"; printf ".read %s/test/sql_v2/receipt_integrity.sql\n" "$PWD"; printf ".read %s/test/sql_v2/receipt_builder_tuples.sql\n" "$PWD"; } > "$out/$version.sql"
  "$cli" -unsigned -no-init -batch -bail -csv -noheader -f "$out/$version.sql" > "$out/$version.log" 2>&1
  grep -q '^PASS$' "$out/$version.log"
  "$cli" -version > "$out/$version.version"
done
cmp "$out/v1.log" "$out/v2.log"
sha256sum "$out/v1.log" "$out/v2.log" "$out/seed.sql" "$out/v1.sql" "$out/v2.sql"
