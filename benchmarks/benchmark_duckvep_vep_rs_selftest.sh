#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
runner=$root/benchmarks/benchmark_duckvep_vep_rs.sh
real_duckdb=$(command -v duckdb)
scratch_root=${BENCHMARK_SELFTEST_TMPDIR:-$(cd "$root/.." && pwd)/artifacts/benchmark}
mkdir -p "$scratch_root"
tmp=$(mktemp -d "$scratch_root/selftest.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin" "$tmp/prefix/bin" "$tmp/work dir's" "$tmp/output dir's" "$tmp/resources"

cat > "$tmp/bin/duckdb" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
sql=$(mktemp)
cat > "$sql"
cat "$sql" >> "$SQL_CAPTURE"
printf '\n-- invocation --\n' >> "$SQL_CAPTURE"
if grep -q 'CREATE TABLE sites AS' "$sql"; then
  printf 'prepared\n' > sites.vcf
  printf 'snapshot\n' > model.dvsnap
  printf 'ids\n' > transcript_ids.parquet
  printf 'regions\n' > regions.parquet
  printf 'prepare\n' >> "$PREPARE_LOG"
elif grep -q '^\.timer on' "$sql"; then
  for stage in compact_count rich_count rich_touch rich_joined_parquet compact_parquet; do
    printf 'STAGE %s\nRun Time (s): real 1.25 user 0.01 sys 0.01\n' "$stage"
  done
fi
rm -f "$sql"
MOCK
cat > "$tmp/bin/time" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
format= output=
while (($#)); do
  case $1 in
    -f) format=$2; shift 2;;
    -a) shift;;
    -o) output=$2; shift 2;;
    *) break;;
  esac
done
printf 'time' >> "$TIME_LOG"
for arg in "$@"; do printf '\t%s' "$arg" >> "$TIME_LOG"; done
printf '\n' >> "$TIME_LOG"
"$@"
format=${format//%e/1.25}
format=${format//%U/1.00}
format=${format//%S/0.25}
format=${format//%M/42}
printf '%b\n' "$format" >> "$output"
MOCK
cat > "$tmp/bin/taskset" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
[[ $1 == -c && $2 == 3 ]] || exit 2
shift 2
printf '%s\n' "$*" >> "$AFFINITY_LOG"
exec "$@"
MOCK
cat > "$tmp/bin/vep-rs" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf 'vep-rs' >> "$ARG_LOG"
for arg in "$@"; do printf '\t%s' "$arg" >> "$ARG_LOG"; done
printf '\n' >> "$ARG_LOG"
while (($#)); do
  if [[ $1 == -o ]]; then
    if [[ ${MOCK_EMPTY_ROWS:-0} == 1 ]]; then printf '#header\n' > "$2";
    else printf 'vep-rs output\n' > "$2"; fi
    break
  fi
  shift
done
MOCK
cat > "$tmp/prefix/bin/vep" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf 'vep' >> "$ARG_LOG"
for arg in "$@"; do printf '\t%s' "$arg" >> "$ARG_LOG"; done
printf '\n' >> "$ARG_LOG"
MOCK
chmod +x "$tmp/bin/duckdb" "$tmp/bin/time" "$tmp/bin/taskset" "$tmp/bin/vep-rs" "$tmp/prefix/bin/vep"

vcf=$tmp/input.vcf.gz
printf 'mock input\n' > "$vcf"
for resource in "reader's.duckdb_extension" "model's.db" "cache.json" "genes.gtf.gz" "reference.fa" "duck'vep.duckdb_extension"; do
  printf 'mock resource\n' > "$tmp/resources/$resource"
done
export PATH="$tmp/bin:$PATH"
export DUCKDB_BIN=duckdb TIME_BIN=$tmp/bin/time TASKSET_BIN=$tmp/bin/taskset
export VEP_RS_BIN=$tmp/bin/vep-rs VEP_BIN=$tmp/prefix/bin/vep VEP_PREFIX=$tmp/prefix
export BENCH_OUT="$tmp/output dir's"
export DUCKVEP_READER_EXT="$tmp/resources/reader's.duckdb_extension"
export DUCKVEP_MODEL="$tmp/resources/model's.db" DUCKVEP_EXT="$tmp/resources/duck'vep.duckdb_extension"
mkdir -p "$tmp/resources/cache/transcripts"
export VEP_RS_CACHE=$tmp/resources/cache VEP_GTF=$tmp/resources/genes.gtf.gz VEP_FASTA=$tmp/resources/reference.fa
export SQL_CAPTURE=$tmp/sql.log PREPARE_LOG=$tmp/prepare.log TIME_LOG=$tmp/time.log
export AFFINITY_LOG=$tmp/affinity.log ARG_LOG=$tmp/args.log
: > "$SQL_CAPTURE"; : > "$PREPARE_LOG"; : > "$TIME_LOG"; : > "$AFFINITY_LOG"; : > "$ARG_LOG"
run_benchmark() {
  if ! "$runner" --vcf "$vcf" --work "$tmp/work dir's" --threads 1 --runs 1 --cpu-affinity 3 \
      > "$tmp/run.log" 2>&1; then
    cat "$tmp/run.log" >&2
    exit 1
  fi
}
run_benchmark
[[ $(wc -l < "$tmp/work dir's/timings.tsv") -eq 3 ]]
[[ $(wc -l < "$tmp/work dir's/stage_timings.tsv") -eq 6 ]]
[[ $(wc -l < "$PREPARE_LOG") -eq 1 ]]
awk -F '\t' 'NR > 1 { if ($1 != "benchmark" || $3 != 1 || $4 != 3 || $6 != 1.25) exit 1; n++ } END { if (n != 2) exit 1 }' "$tmp/work dir's/timings.tsv"
awk -F '\t' 'NR > 1 { if ($1 != "stage" || $2 != "duckvep" || $3 != 1 || $4 != 3 || $6 != 1.25) exit 1; n++ } END { if (n != 5) exit 1 }' "$tmp/work dir's/stage_timings.tsv"
grep -Fq -- $'--fork\t1' "$ARG_LOG"
grep -Fq -- '--force_overwrite' "$ARG_LOG"
grep -Fq -- "$tmp/resources/reader''s.duckdb_extension" "$SQL_CAPTURE"
grep -Fq -- "$tmp/resources/model''s.db" "$SQL_CAPTURE"
grep -Fq -- "$tmp/output dir''s/stage-rich.parquet" "$SQL_CAPTURE"
grep -Fq -- 'SET threads=1' "$SQL_CAPTURE"
[[ $(wc -l < "$AFFINITY_LOG") -eq 3 ]]
grep -Fq "$tmp/bin/vep-rs" "$AFFINITY_LOG"
grep -Fq "$tmp/bin/duckdb -unsigned" "$AFFINITY_LOG"
grep -Fq $'cpu_affinity\t3' "$tmp/work dir's/run_receipt.tsv"

# Preparation is content-bound, including same-size edits with unchanged timestamps.
run_benchmark
[[ $(wc -l < "$PREPARE_LOG") -eq 1 ]]
cp -p "$vcf" "$tmp/original.vcf"
printf 'edit input\n' > "$vcf"
touch -r "$tmp/original.vcf" "$vcf"
run_benchmark
[[ $(wc -l < "$PREPARE_LOG") -eq 2 ]]
printf 'corrupted\n' > "$tmp/work dir's/model.dvsnap"
run_benchmark
[[ $(wc -l < "$PREPARE_LOG") -eq 3 ]]

# Execute the generated replay and adjudication SQL on one-sided and paired mismatches.
mkdir "$tmp/adjudication"
cd "$tmp/adjudication"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > sites.vcf
for id in d_only v_only paired absent; do
  printf '1\t10\t%s\tA\tC\t.\t.\t.\n' "$id" >> sites.vcf
done
for row in 'd_only d' 'v_only v' 'paired other'; do
  read -r id cons <<< "$row"
  printf '%s\t1:10\tC\tg\ttx\tTranscript\t%s' "$id" "$cons" >> disagree.vep.tab
  for ((i=7; i<21; i++)); do printf '\t-' >> disagree.vep.tab; done
  printf '\n' >> disagree.vep.tab
done
cat > check.sql <<'SQL'
CREATE TABLE only_d(id VARCHAR,feature VARCHAR,cons VARCHAR);
INSERT INTO only_d VALUES ('d_only','tx','d'),('paired','tx','d'),('absent','tx','d');
CREATE TABLE only_v(id VARCHAR,feature VARCHAR,cons VARCHAR);
INSERT INTO only_v VALUES ('v_only','tx','v'),('paired','tx','v'),('absent','tx','v');
CREATE TABLE d AS FROM only_d;
CREATE TABLE v AS FROM only_v;
SQL
awk '/^COPY \(SELECT line FROM \(SELECT 0 o, 0::BIGINT k/ { copying=1 }
     /^-- invocation --/ && copying { exit } copying { print }' "$SQL_CAPTURE" >> check.sql
awk '/^CREATE OR REPLACE TABLE p/ { copying=1 }
     /^-- invocation --/ && copying { exit } copying { print }' "$SQL_CAPTURE" >> check.sql
"$real_duckdb" -bail -csv < check.sql > verdicts.csv
[[ $(wc -l < disagree.vcf) -eq 6 ]]
for verdict in 'VEP agrees with DuckVEP' 'VEP agrees with vep-rs' 'VEP differs from both' 'VEP has no row'; do
  grep -Fxq "$verdict,1" verdicts.csv
done
if VEP_RS_CACHE="$tmp/resources" "$runner" --vcf "$vcf" --work "$tmp/work dir's" --runs 1 > "$tmp/bad-cache.log" 2>&1; then
  echo 'invalid cache root was accepted' >&2; exit 1
fi
if MOCK_EMPTY_ROWS=1 "$runner" --vcf "$vcf" --work "$tmp/work dir's" --threads 1 --runs 1 > "$tmp/empty.log" 2>&1; then
  echo 'empty annotation output was accepted' >&2; exit 1
fi
printf 'benchmark runner selftest: PASS\n'
