#!/bin/sh
set -eu
if [ -z "${DUCKHTS_EXTENSION:-}" ]; then
  echo 'SKIP DuckHTS integration (DUCKHTS_EXTENSION unset)'
  exit 0
fi
if [ ! -f "$DUCKHTS_EXTENSION" ]; then
  echo 'DUCKHTS_EXTENSION does not name an extension file' >&2
  exit 1
fi
"${DUCKDB_CLI:-duckdb}" -unsigned \
  -cmd "LOAD '$DUCKHTS_EXTENSION'" \
  -cmd "LOAD '$(pwd)/build/release/duckvep.duckdb_extension'" \
  -cmd '.read test/data/duckvep/readme.sql' \
  -cmd "CREATE TABLE integration_bcf_events AS SELECT row_number() OVER ()::UBIGINT AS event_index, 1::UINTEGER AS seq_region, POS::UBIGINT AS position, REF AS reference, ALT[1] AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM read_bcf('test/data/duckvep/minimal_bcsq.vcf') WHERE POS = 124" \
  -c "SELECT event_index, transcript_hgvs FROM query(duckvep_annotate_sql('integration_bcf_events', 'readme', struct_pack(hgvs := true)))" \
  | grep 'c.5T>C' >/dev/null
