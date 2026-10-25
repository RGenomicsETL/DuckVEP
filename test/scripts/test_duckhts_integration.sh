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
  -cmd "CREATE TABLE integration_regions AS SELECT * FROM (VALUES (0::UINTEGER, 'chrEmpty'), (1::UINTEGER, 'chrDuck')) t(seq_region, seq_region_name)" \
  -cmd "CREATE VIEW integration_records AS SELECT * FROM read_bcf('test/data/duckvep/minimal_bcsq.vcf')" \
  -cmd "CREATE TABLE integration_bcf_events AS SELECT * FROM query(duckvep_vcf_events_sql('integration_records', 'integration_regions')) WHERE event_status = 'ok' ORDER BY event_index" \
  -c "SELECT source_id, transcript_hgvs FROM query(duckvep_annotate_sql('integration_bcf_events', 'readme', struct_pack(hgvs := true))) JOIN integration_bcf_events USING (event_index)" \
  | grep 'duck_missense.*c.5T>C' >/dev/null
