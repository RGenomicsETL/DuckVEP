#!/usr/bin/env bash
# Compare DuckVEP with vep-rs on the same sites-only VCF; executable Ensembl VEP adjudicates disagreements.
# One ALT allele per record and unique VCF IDs make output rows joinable on (ID, transcript).
#
#   benchmark_duckvep_vep_rs.sh --vcf INPUT.vcf.gz --work DIR [--threads N] [--runs N] [--cpu-affinity LIST]
#
# Required environment:
#   DUCKVEP_READER_EXT  DuckHTS extension used to prepare the VCF
#   DUCKVEP_MODEL       DuckVEP model database
#   VEP_RS_CACHE        vep-rs JSON cache
#   VEP_GTF             GTF for executable Ensembl VEP
#   VEP_FASTA           indexed reference FASTA
#   VEP_RS_DIR          directory with `vep`, unless VEP_RS_BIN is set
#
# Executables can be selected with DUCKDB_BIN, TIME_BIN, TASKSET_BIN, VEP_RS_BIN and VEP_BIN.
# CPU_AFFINITY sets a default CPU list; --cpu-affinity overrides it and uses taskset -c.
# TIME_BIN must accept GNU time's -f/-a/-o options.
# Optional VEP_PREFIX adds its bin directory to PATH for the adjudicating VEP run.
# BENCH_OUT stores output files (default /dev/shm/duckvep-vep-rs).
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
caller_pwd=$(pwd -P)
vcf=; work=; threads=16; runs=3; cpu_affinity=${CPU_AFFINITY:-}
while [[ $# -gt 0 ]]; do
  case $1 in
    --vcf|--work|--threads|--runs|--cpu-affinity)
      [[ $# -ge 2 ]] || { echo "missing value for $1" >&2; exit 2; }
      case $1 in
        --vcf) vcf=$2;; --work) work=$2;; --threads) threads=$2;; --runs) runs=$2;; --cpu-affinity) cpu_affinity=$2;;
      esac
      shift 2;;
    --help) sed -n '1,24p' "$0"; exit 0;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
done
[[ -n $vcf && -n $work ]] || { echo "need --vcf and --work" >&2; exit 2; }
[[ $threads =~ ^[1-9][0-9]*$ && $runs =~ ^[1-9][0-9]*$ ]] || { echo "--threads and --runs must be positive integers" >&2; exit 2; }

absolute_path() {
  case $1 in /*) printf '%s\n' "$1";; *) printf '%s/%s\n' "$caller_pwd" "$1";; esac
}
absolute_command_path() {
  case $1 in */*) absolute_path "$1";; *) command -v "$1";; esac
}
sql_literal() {
  local escaped=${1//\'/\'\'}
  printf "'%s'" "$escaped"
}
require_command() {
  command -v "$1" >/dev/null 2>&1 || { echo "executable not found: $1" >&2; exit 2; }
}

vcf=$(absolute_path "$vcf")
work=$(absolute_path "$work")
[[ -f $vcf ]] || { echo "VCF not found: $vcf" >&2; exit 2; }
: "${DUCKVEP_READER_EXT:?}" "${DUCKVEP_MODEL:?}" "${VEP_RS_CACHE:?}" "${VEP_GTF:?}" "${VEP_FASTA:?}"
reader_ext=$(absolute_path "$DUCKVEP_READER_EXT")
model=$(absolute_path "$DUCKVEP_MODEL")
vep_cache=$(absolute_path "$VEP_RS_CACHE")
[[ -d $vep_cache/transcripts ]] || { echo "VEP_RS_CACHE must contain transcripts/: $vep_cache" >&2; exit 2; }
vep_gtf=$(absolute_path "$VEP_GTF")
vep_fasta=$(absolute_path "$VEP_FASTA")
ext=$(absolute_path "${DUCKVEP_EXT:-$root/build/release/duckvep.duckdb_extension}")
if [[ -n ${VEP_PREFIX:-} ]]; then
  vep_prefix=$(absolute_path "$VEP_PREFIX")
  export PATH="$vep_prefix/bin:$PATH"
fi
duckdb_bin=$(absolute_command_path "${DUCKDB_BIN:-duckdb}")
time_bin=$(absolute_command_path "${TIME_BIN:-/usr/bin/time}")
if [[ -n ${VEP_RS_BIN:-} ]]; then
  vep_rs_bin=$(absolute_command_path "$VEP_RS_BIN")
else
  : "${VEP_RS_DIR:?set VEP_RS_DIR or VEP_RS_BIN}"
  vep_rs_bin=$(absolute_path "$VEP_RS_DIR/vep")
fi
vep_bin=$(absolute_command_path "${VEP_BIN:-vep}")
out=$(absolute_path "${BENCH_OUT:-/dev/shm/duckvep-vep-rs}")
require_command sha256sum
require_command "$duckdb_bin"
require_command "$time_bin"
require_command "$vep_rs_bin"
require_command "$vep_bin"
taskset_bin=none
if [[ -n $cpu_affinity ]]; then
  taskset_bin=$(absolute_command_path "${TASKSET_BIN:-taskset}")
  require_command "$taskset_bin"
  affinity_cmd=("$taskset_bin" -c "$cpu_affinity")
  affinity_label=$cpu_affinity
else
  affinity_cmd=()
  affinity_label=none
fi
mkdir -p "$work" "$out"
work=$(cd "$work" && pwd -P)
out=$(cd "$out" && pwd -P)
cd "$work"

sql_reader_ext=$(sql_literal "$reader_ext")
sql_ext=$(sql_literal "$ext")
sql_model=$(sql_literal "$model")
sql_vcf=$(sql_literal "$vcf")
sql_duckvep_output=$(sql_literal "$out/duckvep.parquet")
sql_stage_rich=$(sql_literal "$out/stage-rich.parquet")
sql_stage_compact=$(sql_literal "$out/stage-compact.parquet")

# Reuse preparation only when all generated files and their input identities still match.
prep_expected=$work/preparation.expected.tsv
prep_receipt=$work/preparation.receipt.tsv
for input in "$vcf" "$reader_ext" "$model" "$ext" "$duckdb_bin" "$root/benchmarks/benchmark_duckvep_vep_rs.sh"; do
  digest=$(sha256sum -- "$input")
  printf '%s\t%s\n' "$input" "${digest%% *}"
done > "$prep_expected"
if ! cmp -s "$prep_expected" "$prep_receipt" ||
   ! sha256sum --check --status preparation_outputs.sha256 2>/dev/null; then
  rm -f sites.vcf model.dvsnap transcript_ids.parquet regions.parquet "$prep_receipt"
  "$duckdb_bin" -unsigned -bail <<SQL
LOAD $sql_reader_ext;
CREATE TABLE sites AS SELECT CHROM, POS, REF, unnest(ALT) AS ALT FROM read_bcf($sql_vcf);
COPY (SELECT line FROM (
  SELECT 0 o, 0 r, '##fileformat=VCFv4.2' line UNION ALL SELECT 1, 0,
    '#CHROM' || chr(9) || 'POS' || chr(9) || 'ID' || chr(9) || 'REF' || chr(9) || 'ALT' || chr(9) || 'QUAL' || chr(9) || 'FILTER' || chr(9) || 'INFO'
  UNION ALL SELECT 2, row_number() OVER (), CHROM || chr(9) || POS || chr(9) || 'e' || row_number() OVER () || chr(9) ||
    REF || chr(9) || ALT || chr(9) || '.' || chr(9) || '.' || chr(9) || '.' FROM sites) ORDER BY o, r)
 TO 'sites.vcf' (FORMAT csv, HEADER false, QUOTE '', ESCAPE '', DELIMITER '\x01');
LOAD $sql_ext;
ATTACH $sql_model AS m (READ_ONLY);
SELECT loaded FROM duckvep_model_load('g', 'SELECT seq_region, sequence_length FROM m.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM m.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM m.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM m.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM m.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 transcript_coverage_complete := TRUE);
SELECT duckvep_model_save('g', 'model.dvsnap');
COPY (SELECT transcript_index, transcript_stable_id FROM m.model_transcripts) TO 'transcript_ids.parquet' (FORMAT parquet);
COPY (SELECT seq_region, seq_region_name FROM m.model_regions) TO 'regions.parquet' (FORMAT parquet);
SQL
  for prepared in sites.vcf model.dvsnap transcript_ids.parquet regions.parquet; do
    [[ -s $prepared ]] || { echo "preparation did not produce $work/$prepared" >&2; exit 1; }
  done
  sha256sum sites.vcf model.dvsnap transcript_ids.parquet regions.parquet > preparation_outputs.sha256
  cp "$prep_expected" "$prep_receipt"
fi
rm -f "$prep_expected"

duckvep_setup_sql() {
  cat <<SQL
SET threads=$1; SET memory_limit='8GB';
LOAD $sql_ext;
SELECT duckvep_model_restore('g', 'model.dvsnap');
CREATE TEMP TABLE ev AS
 SELECT row_number() OVER (ORDER BY r.seq_region, v.pos, v.id)::UBIGINT AS event_index, v.id AS vcf_id,
        r.seq_region::UINTEGER AS seq_region, v.pos::UBIGINT AS "position", v.ref AS reference, v.alt AS alternate
 FROM read_csv('sites.vcf', delim='\t', header=false, skip=2, quote='', escape='', auto_detect=false,
   columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR'}) v
 JOIN read_parquet('regions.parquet') r ON r.seq_region_name = v.chrom;
CREATE TEMP TABLE alleles AS SELECT event_index, seq_region, "position", reference, alternate, NULL::UBIGINT AS end_position,
  NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
 FROM ev ORDER BY event_index;
SQL
}
duckvep_sql() {
  duckvep_setup_sql "$1"
  cat <<SQL
COPY (SELECT e.vcf_id, t.transcript_stable_id, a.*
      FROM query(duckvep_annotate_sql('alleles', 'g', struct_pack(rich := true))) a
      JOIN ev e USING (event_index) LEFT JOIN read_parquet('transcript_ids.parquet') t USING (transcript_index))
 TO $sql_duckvep_output (FORMAT parquet);
SQL
}
duckvep_stage_sql() {
  duckvep_setup_sql 1
  cat <<SQL
.timer on
.print STAGE compact_count
SELECT count(*) FROM query(duckvep_annotate_sql('alleles', 'g'));
.print STAGE rich_count
SELECT count(*) FROM query(duckvep_annotate_sql('alleles', 'g', struct_pack(rich := true)));
.print STAGE rich_touch
SELECT sum(hash(a)::HUGEINT) IS NOT NULL FROM query(duckvep_annotate_sql('alleles', 'g', struct_pack(rich := true))) a;
.print STAGE rich_joined_parquet
COPY (SELECT e.vcf_id, t.transcript_stable_id, a.* FROM query(duckvep_annotate_sql('alleles', 'g', struct_pack(rich := true))) a
      JOIN ev e USING (event_index) LEFT JOIN read_parquet('transcript_ids.parquet') t USING (transcript_index))
 TO $sql_stage_rich (FORMAT parquet);
.print STAGE compact_parquet
COPY (SELECT * FROM query(duckvep_annotate_sql('alleles', 'g'))) TO $sql_stage_compact (FORMAT parquet);
SQL
}

# End-to-end measurements include process startup, extension load, restore, VCF read/sort, annotation and output.
printf 'kind\ttool\tthreads\tcpu_affinity\trun\twall_s\tuser_s\tsys_s\tmaxrss_kib\n' > timings.tsv
{
  printf 'key\tvalue\n'
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'vcf\t%s\nmodel\t%s\nduckvep_extension\t%s\nduckhts_reader_extension\t%s\n' "$vcf" "$model" "$ext" "$reader_ext"
  printf 'vep_rs_cache\t%s\nvep_gtf\t%s\nvep_fasta\t%s\nbench_out\t%s\n' "$vep_cache" "$vep_gtf" "$vep_fasta" "$out"
  printf 'duckdb_bin\t%s\nvep_rs_bin\t%s\nvep_bin\t%s\ntime_bin\t%s\n' "$duckdb_bin" "$vep_rs_bin" "$vep_bin" "$time_bin"
  printf 'taskset_bin\t%s\n' "$taskset_bin"
  printf 'configured_threads\t%s\nruns_per_thread_setting\t%s\ncpu_affinity\t%s\n' "$threads" "$runs" "$affinity_label"
  printf 'vep_rs_arguments\tassembly=GRCh38;tab;fork={threads};force_overwrite\n'
  printf 'duckvep_settings\tmemory_limit=8GB;rich=true;parquet\n'
  printf 'source_revision\t%s\n' "$(git -C "$root" rev-parse HEAD)"
  for artifact in "$duckdb_bin" "$vep_rs_bin" "$vep_bin" "$root/benchmarks/benchmark_duckvep_vep_rs.sh"; do
    digest=$(sha256sum -- "$artifact")
    printf 'sha256:%s\t%s\n' "$artifact" "${digest%% *}"
  done
} > run_receipt.tsv
thread_values=("$threads")
[[ $threads == 1 ]] || thread_values+=(1)
for t in "${thread_values[@]}"; do
  for run in $(seq 1 "$runs"); do
    if ! "$time_bin" -f "benchmark\tvep-rs\t$t\t$affinity_label\t$run\t%e\t%U\t%S\t%M" -a -o timings.tsv \
      "${affinity_cmd[@]}" "$vep_rs_bin" --json_cache "$vep_cache" --assembly GRCh38 \
      -i sites.vcf -o "$out/veprs.tab" --tab --fork "$t" --force_overwrite > veprs.log 2>&1; then
      cat veprs.log >&2
      echo "vep-rs measurement failed (threads=$t run=$run)" >&2
      exit 1
    fi
    awk '!/^#/ && NF { found=1; exit } END { exit !found }' "$out/veprs.tab" || {
      echo "vep-rs produced no annotation rows; check its cache and input" >&2
      exit 1
    }
    duckvep_sql "$t" > run.sql
    "$time_bin" -f "benchmark\tduckvep\t$t\t$affinity_label\t$run\t%e\t%U\t%S\t%M" -a -o timings.tsv \
      "${affinity_cmd[@]}" "$duckdb_bin" -unsigned -bail < run.sql > /dev/null
  done
done

# Query-only stage measurements follow the published probe: setup/restore is outside .timer; timings go to a separate receipt.
printf 'kind\ttool\tthreads\tcpu_affinity\tstage\tseconds\n' > stage_timings.tsv
rm -f "$out/stage-rich.parquet" "$out/stage-compact.parquet"
duckvep_stage_sql > stages.sql
if ! "${affinity_cmd[@]}" "$duckdb_bin" -unsigned -bail < stages.sql > stages.log 2>&1; then
  cat stages.log >&2
  echo "DuckVEP stage probe failed" >&2
  exit 1
fi
awk -v affinity="$affinity_label" '
  $1 == "STAGE" { stage = $2; next }
  $1 == "Run" && $2 == "Time" && $4 == "real" {
    if (stage == "") exit 1
    printf "stage\tduckvep\t1\t%s\t%s\t%s\n", affinity, stage, $5
    stage = ""
    count++
  }
  END { if (count != 5) exit 1 }
' stages.log > stage_timings.rows || { echo "could not parse all five DuckVEP stages; see $work/stages.log" >&2; exit 1; }
cat stage_timings.rows >> stage_timings.tsv
rm -f stage_timings.rows "$out/stage-rich.parquet" "$out/stage-compact.parquet"
printf '%s\n' 'End-to-end benchmark timings:'
cat timings.tsv
printf '%s\n' 'DuckVEP query-only stage timings (not included above):'
cat stage_timings.tsv

# ---- Tuples: (ID, transcript, sorted consequence set), on the transcripts both tools emit ------------------
tab_columns="columns={'c00':'VARCHAR','c01':'VARCHAR','c02':'VARCHAR','c03':'VARCHAR','c04':'VARCHAR','c05':'VARCHAR','c06':'VARCHAR','c07':'VARCHAR','c08':'VARCHAR','c09':'VARCHAR','c10':'VARCHAR','c11':'VARCHAR','c12':'VARCHAR','c13':'VARCHAR','c14':'VARCHAR','c15':'VARCHAR','c16':'VARCHAR','c17':'VARCHAR','c18':'VARCHAR','c19':'VARCHAR','c20':'VARCHAR'}, null_padding=true"
sql_duckvep_parquet=$(sql_literal "$out/duckvep.parquet")
sql_veprs_tab=$(sql_literal "$out/veprs.tab")
rm -f compare.duckdb
"$duckdb_bin" -bail compare.duckdb <<SQL
SET threads=$threads; SET preserve_insertion_order=false;
CREATE TABLE d AS SELECT vcf_id AS id, coalesce(transcript_stable_id, '-') AS feature,
  array_to_string(list_sort(string_split_regex(consequence, '[,&]')), ',') AS cons FROM read_parquet($sql_duckvep_parquet);
CREATE TABLE v AS SELECT c00 AS id, c04 AS feature, array_to_string(list_sort(string_split(c06, ',')), ',') AS cons
  FROM read_csv($sql_veprs_tab, delim='\t', header=false, comment='#', quote='', escape='', auto_detect=false, $tab_columns);
CREATE TABLE shared AS SELECT feature FROM (SELECT DISTINCT feature FROM d) JOIN (SELECT DISTINCT feature FROM v) USING (feature) WHERE feature <> '-';
CREATE TABLE only_d AS SELECT d.* FROM d SEMI JOIN shared USING (feature) ANTI JOIN v USING (id, feature, cons);
CREATE TABLE only_v AS SELECT v.* FROM v SEMI JOIN shared USING (feature) ANTI JOIN d USING (id, feature, cons);
SELECT (SELECT count(*) FROM d) AS duckvep_rows, (SELECT count(*) FROM v) AS veprs_rows,
  (SELECT count(DISTINCT feature) FROM d WHERE feature <> '-') AS duckvep_transcripts,
  (SELECT count(DISTINCT feature) FROM v WHERE feature <> '-') AS veprs_transcripts, (SELECT count(*) FROM shared) AS shared_transcripts;
SELECT (SELECT count(*) FROM d SEMI JOIN shared USING (feature)) AS duckvep_tuples,
  (SELECT count(*) FROM v SEMI JOIN shared USING (feature)) AS veprs_tuples,
  (SELECT count(*) FROM d SEMI JOIN shared USING (feature) SEMI JOIN v USING (id, feature, cons)) AS matching,
  (SELECT count(*) FROM only_d) AS only_duckvep, (SELECT count(*) FROM only_v) AS only_veprs;
COPY (SELECT line FROM (SELECT 0 o, 0::BIGINT k, '##fileformat=VCFv4.2' line UNION ALL SELECT 1, 0,
  '#CHROM' || chr(9) || 'POS' || chr(9) || 'ID' || chr(9) || 'REF' || chr(9) || 'ALT' || chr(9) || 'QUAL' || chr(9) || 'FILTER' || chr(9) || 'INFO'
  UNION ALL SELECT 2, row_number() OVER (ORDER BY chrom, pos::BIGINT), chrom || chr(9) || pos || chr(9) || id || chr(9) || ref || chr(9) || alt || chr(9) || '.' || chr(9) || '.' || chr(9) || '.'
  FROM (SELECT column2 AS id, column0 AS chrom, column1 AS pos, column3 AS ref, column4 AS alt
        FROM read_csv('sites.vcf', delim='\t', header=false, skip=2, quote='', all_varchar=true)) s
        SEMI JOIN (SELECT id FROM only_d UNION SELECT id FROM only_v) disagreements USING (id)) ORDER BY o, k)
 TO 'disagree.vcf' (FORMAT csv, HEADER false, QUOTE '', DELIMITER '\x01');
SQL

# ---- Executable VEP on every site where the two tools disagree --------------------------------------------
if ! "$vep_bin" --gtf "$vep_gtf" --fasta "$vep_fasta" -i disagree.vcf -o disagree.vep.tab --tab \
  --force_overwrite --no_stats --fork 8 > vep.log 2>&1; then
  cat vep.log >&2
  echo "executable VEP adjudication failed" >&2
  exit 1
fi
"$duckdb_bin" -bail compare.duckdb <<SQL
CREATE OR REPLACE TABLE p AS SELECT c00 AS id, c04 AS feature, array_to_string(list_sort(string_split(c06, ',')), ',') AS cons
  FROM read_csv('disagree.vep.tab', delim='\t', header=false, comment='#', quote='', escape='', auto_detect=false, $tab_columns);
SELECT CASE WHEN p.cons IS NULL THEN 'VEP has no row' WHEN p.cons = d.cons THEN 'VEP agrees with DuckVEP'
            WHEN p.cons = v.cons THEN 'VEP agrees with vep-rs' ELSE 'VEP differs from both' END AS verdict, count(*) AS n
FROM only_d d FULL JOIN only_v v USING (id, feature)
LEFT JOIN p ON p.id=coalesce(d.id,v.id) AND p.feature=coalesce(d.feature,v.feature)
GROUP BY ALL ORDER BY n DESC;
CREATE OR REPLACE TEMP TABLE pid AS SELECT DISTINCT id FROM p;
SELECT (SELECT count(*) FROM p) AS vep_tuples,
  (SELECT count(*) FROM p SEMI JOIN d USING (id, feature, cons)) AS duckvep_matches,
  (SELECT count(*) FROM p SEMI JOIN v USING (id, feature, cons)) AS veprs_matches,
  (SELECT count(*) FROM p ANTI JOIN d USING (id, feature)) AS vep_rows_absent_from_duckvep,
  (SELECT count(*) FROM p ANTI JOIN v USING (id, feature)) AS vep_rows_absent_from_veprs,
  (SELECT count(*) FROM (SELECT * FROM d SEMI JOIN pid USING (id)) x ANTI JOIN p USING (id, feature)) AS duckvep_rows_absent_from_vep,
  (SELECT count(*) FROM (SELECT * FROM v SEMI JOIN pid USING (id)) x ANTI JOIN p USING (id, feature)) AS veprs_rows_absent_from_vep;
SQL
