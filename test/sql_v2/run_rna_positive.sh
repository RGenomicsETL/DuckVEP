#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
cli=${DUCKVEP_V2_DUCKDB:?Set DUCKVEP_V2_DUCKDB to the pinned v2 CLI}
ext=${DUCKVEP_V2_EXTENSION:-build/release_v2/duckvep.duckdb_extension}
tmp=$(mktemp -d)
trap 'rc=$?; if ((rc)); then echo "fixture failed ($rc): $tmp" >&2; for log in "$tmp"/*.log; do [[ -f $log ]] && { echo "$log" >&2; tail -50 "$log" >&2; }; done; else rm -rf "$tmp"; fi' EXIT
cp "$ext" "$tmp/duckvep.duckdb_extension"
sha256sum "$tmp/duckvep.duckdb_extension"
# Use the actual source-owned seed, without a second transcription of Ensembl rows.
awk '/^CREATE SCHEMA duckvep_core;$/ {start=1} start && /^statement ok$/ {next} start && /^----$/ {exit} start && /^# / {next} start && /^query / {exit} start {print}' test/sql/duckvep_ensembl.test > "$tmp/seed.sql"
# The seed stops before the first query assertion; source authority is built by the fixture.
for phase in v2_rna_positive v2_rna_replay; do
  { printf "LOAD '%s/duckvep.duckdb_extension';\n" "$tmp"; if [[ $phase == v2_rna_positive ]]; then cat "$tmp/seed.sql"; fi; sed "s|__RNA_TMP__|$tmp|g" "test/sql_v2/$phase.sql"; } > "$tmp/$phase.sql"
done
"$cli" -unsigned -bail -batch -csv -noheader -f "$tmp/v2_rna_positive.sql" > "$tmp/live.log" 2>&1
for mode in mmap no-mmap; do
  if [[ $mode == mmap ]]; then DUCKVEP_SNAPSHOT_MMAP=1 "$cli" -unsigned -bail -batch -csv -noheader -f "$tmp/v2_rna_replay.sql" > "$tmp/$mode.log" 2>&1; else DUCKVEP_SNAPSHOT_MMAP=0 "$cli" -unsigned -bail -batch -csv -noheader -f "$tmp/v2_rna_replay.sql" > "$tmp/$mode.log" 2>&1; fi
done
cmp "$tmp/mmap.log" "$tmp/no-mmap.log"
for mode in live mmap no-mmap; do grep -E '^(fingerprint-|full-)' "$tmp/$mode.log" > "$tmp/$mode.full"; done
cmp "$tmp/live.full" "$tmp/mmap.full"
cmp "$tmp/live.full" "$tmp/no-mmap.full"
sha256sum "$tmp/reverse.dvsnap" "$tmp/forward.dvsnap" "$tmp/live.log" "$tmp/mmap.log" "$tmp/no-mmap.log" "$tmp/live.full" "$tmp/mmap.full" "$tmp/no-mmap.full"
printf '\nLIVE\n'; cat "$tmp/live.log"
printf '\nREPLAY mmap / no-mmap identical\n'; cat "$tmp/mmap.log"
