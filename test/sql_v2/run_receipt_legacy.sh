#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
out=$(realpath "${1:?Provide the receipt-integrity evidence directory}")
cli=${DUCKVEP_V1_DUCKDB:?Set DUCKVEP_V1_DUCKDB to the pinned v1 CLI}
baseline=${DUCKVEP_RECEIPT_BASELINE:?Set DUCKVEP_RECEIPT_BASELINE to the preserved extension}
mkdir -p "$out/baseline"
cp "$baseline" "$out/baseline/duckvep.duckdb_extension"
sha256sum "$out/baseline/duckvep.duckdb_extension" > "$out/baseline.sha256"
{
  printf "LOAD '%s/v1/duckvep.duckdb_extension';\n" "$out"
  printf ".read %s/seed.sql\n" "$out"
  printf "CREATE TABLE receipt_regions AS SELECT * FROM query(duckvep_ensembl_regions_sql('duckvep_core', 'duckvep_reference_chunks', 'GRCh38'));\n"
  printf "CREATE TABLE receipt_model AS SELECT * FROM query(duckvep_ensembl_transcripts_sql('duckvep_core', 'duckvep_reference_chunks', 'GRCh38'));\n"
  printf "CREATE TABLE receipt_legacy AS SELECT * EXCLUDE (rna_cds_cdna_start, rna_cds_cdna_end, rna_codon_table) FROM receipt_model;\n"
  printf "COPY receipt_legacy TO '%s/legacy.parquet' (FORMAT PARQUET);\n" "$out"
  printf "COPY receipt_regions TO '%s/regions.parquet' (FORMAT PARQUET);\n" "$out"
} > "$out/export-legacy.sql"
"$cli" -unsigned -no-init -batch -bail -f "$out/export-legacy.sql" > "$out/export-legacy.log" 2>&1
for version in baseline v1 v2; do
  ext=$out/$version/duckvep.duckdb_extension
  client=$cli
  if [[ $version == v2 ]]; then client=${DUCKVEP_V2_DUCKDB:?Set pinned v2 CLI}; fi
  {
    printf "LOAD '%s';\n" "$ext"
    printf "CREATE VIEW receipt_regions AS SELECT * FROM read_parquet('%s/regions.parquet');\n" "$out"
    printf "CREATE VIEW receipt_legacy AS SELECT * FROM read_parquet('%s/legacy.parquet');\n" "$out"
    printf "SELECT to_json(r) FROM query(duckvep_model_receipt_sql('receipt_regions', 'receipt_legacy', 'Ensembl', '116', 'GRCh38', repeat('a', 64), repeat('b', 64), 'fixture transcripts')) r;\n"
  } > "$out/$version-legacy.sql"
  "$client" -unsigned -no-init -batch -bail -csv -noheader -f "$out/$version-legacy.sql" > "$out/$version-legacy.log" 2>&1
done
cmp "$out/baseline-legacy.log" "$out/v1-legacy.log"
cmp "$out/baseline-legacy.log" "$out/v2-legacy.log"
sha256sum "$out"/{baseline,v1,v2}-legacy.log "$out"/{legacy,regions}.parquet
