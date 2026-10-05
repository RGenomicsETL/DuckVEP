#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
V1_CLI=${DUCKVEP_V1_DUCKDB:-duckdb}
V2_CLI=${DUCKVEP_V2_DUCKDB:-}
V1_EXTENSION=${DUCKVEP_V1_EXTENSION:-$ROOT/build/release/duckvep.duckdb_extension}
V2_EXTENSION=${DUCKVEP_V2_EXTENSION:-$ROOT/build/release_v2/duckvep.duckdb_extension}
VCF=${DUCKVEP_HG002_VCF:-}
MODEL=${DUCKVEP_HG002_MODEL:-}
V1_REFERENCE=${DUCKVEP_HG002_V1_PARQUET:-}
OUT_DIR=
CPU_LIST=${DUCKVEP_CPU_LIST:-2}
THREADS=${DUCKVEP_THREADS:-1}
MEMORY_LIMIT=${DUCKVEP_MEMORY_LIMIT:-8GB}
TIME_BIN=${DUCKVEP_TIME:-/usr/bin/time}

usage() {
    cat <<'USAGE'
Usage: check_v2_hg002.sh --v2-cli PATH --vcf PATH --model PATH --out-dir DIR [options]

Runs the current v1 and pinned-v2 HG002 haplotype paths, writes Parquet outputs,
and compares the complete rows as an unordered multiset. The v1 reference is an
optional additional exact comparison, not a replacement for the v1 run.

Required options:
  --v2-cli PATH          DuckDB CLI built at the pinned v2 SDK revision
  --vcf PATH             HG002 VCF/BCF input consumed by duckvep_coding_calls
  --model PATH           Ensembl model DuckDB database
  --out-dir DIR          New or empty directory for SQL, logs, metrics and Parquet

Optional options:
  --v1-cli PATH          v1-compatible DuckDB CLI (default: duckdb on PATH)
  --v1-extension PATH    v1 extension (default: build/release/duckvep.duckdb_extension)
  --v2-extension PATH    v2 extension (default: build/release_v2/duckvep.duckdb_extension)
  --v1-reference PATH    saved v1 Parquet for a further exact comparison
  --cpu-list LIST        taskset CPU list (default: 2)
  --threads N            DuckDB threads per process (default: 1)
  --memory-limit VALUE   DuckDB memory limit (default: 8GB)

The v2 CLI must match v2_host.duckdb_sdk_revision in duckvep-package.json. The
runner requires taskset, GNU time, jq, and DuckDB's -column error mode.
USAGE
}

fail() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

need_value() {
    (($# >= 2)) || fail "missing value for $1"
}

while (($#)); do
    case $1 in
        -h|--help) usage; exit 0 ;;
        --v1-cli) need_value "$@"; V1_CLI=$2; shift 2 ;;
        --v2-cli) need_value "$@"; V2_CLI=$2; shift 2 ;;
        --v1-extension) need_value "$@"; V1_EXTENSION=$2; shift 2 ;;
        --v2-extension) need_value "$@"; V2_EXTENSION=$2; shift 2 ;;
        --vcf) need_value "$@"; VCF=$2; shift 2 ;;
        --model) need_value "$@"; MODEL=$2; shift 2 ;;
        --v1-reference) need_value "$@"; V1_REFERENCE=$2; shift 2 ;;
        --out-dir) need_value "$@"; OUT_DIR=$2; shift 2 ;;
        --cpu-list) need_value "$@"; CPU_LIST=$2; shift 2 ;;
        --threads) need_value "$@"; THREADS=$2; shift 2 ;;
        --memory-limit) need_value "$@"; MEMORY_LIMIT=$2; shift 2 ;;
        *) fail "unknown option: $1 (use --help)" ;;
    esac
done

[[ $THREADS =~ ^[1-9][0-9]*$ ]] || fail '--threads must be a positive integer'
[[ -n $CPU_LIST ]] || fail '--cpu-list must not be empty'
[[ -n $V2_CLI ]] || fail 'set --v2-cli or DUCKVEP_V2_DUCKDB'
[[ -n $VCF ]] || fail 'set --vcf or DUCKVEP_HG002_VCF'
[[ -n $MODEL ]] || fail 'set --model or DUCKVEP_HG002_MODEL'
[[ -n $OUT_DIR ]] || fail 'set --out-dir'

resolve_command() {
    local result
    result=$(command -v "$1") || fail "cannot find executable: $1"
    [[ -x $result ]] || fail "not executable: $result"
    printf '%s\n' "$result"
}

V1_CLI=$(resolve_command "$V1_CLI")
V2_CLI=$(resolve_command "$V2_CLI")
TIME_BIN=$(resolve_command "$TIME_BIN")
TASKSET=$(resolve_command taskset)
JQ=$(resolve_command jq)
SHA256=$(resolve_command sha256sum)

[[ -f $V1_EXTENSION && -r $V1_EXTENSION ]] || fail "v1 extension is not readable: $V1_EXTENSION"
[[ -f $V2_EXTENSION && -r $V2_EXTENSION ]] || fail "v2 extension is not readable: $V2_EXTENSION"
[[ -f $VCF && -r $VCF ]] || fail "VCF/BCF is not readable: $VCF"
[[ -f $MODEL && -r $MODEL ]] || fail "model database is not readable: $MODEL"
if [[ -n $V1_REFERENCE ]]; then
    [[ -f $V1_REFERENCE && -r $V1_REFERENCE ]] || fail "v1 reference is not readable: $V1_REFERENCE"
fi
"$TASKSET" -c "$CPU_LIST" true >/dev/null 2>&1 || fail "CPU list is not available: $CPU_LIST"

mkdir -p -- "$OUT_DIR" || fail "cannot create output directory: $OUT_DIR"
OUT_DIR=$(cd -- "$OUT_DIR" && pwd -P)
[[ -z $(find "$OUT_DIR" -mindepth 1 -maxdepth 1 -print -quit) ]] || fail "output directory is not empty: $OUT_DIR"
mkdir "$OUT_DIR/v1-tmp" "$OUT_DIR/v2-tmp"
"$SHA256" -- "$V1_CLI" "$V2_CLI" "$V1_EXTENSION" "$V2_EXTENSION" "$MODEL" "$VCF" \
    "$ROOT/scripts/check_v2_hg002.sh" "$ROOT/duckvep-package.json" > "$OUT_DIR/inputs.sha256"
if [[ -n $V1_REFERENCE ]]; then "$SHA256" -- "$V1_REFERENCE" >> "$OUT_DIR/inputs.sha256"; fi

V2_PIN=$("$JQ" -er '.v2_host.duckdb_sdk_revision' "$ROOT/duckvep-package.json")
[[ $V2_PIN =~ ^[0-9a-f]{40}$ ]] || fail 'could not read the pinned v2 source revision'
V2_SOURCE=$("$TASKSET" -c "$CPU_LIST" "$V2_CLI" -no-init -batch -list -noheader -c 'SELECT source_id FROM pragma_version()' 2>"$OUT_DIR/v2-version-check.log") || fail 'could not query the v2 DuckDB source revision'
V2_SOURCE=${V2_SOURCE//$'\n'/}
V2_SOURCE=${V2_SOURCE//$'\r'/}
[[ ${#V2_SOURCE} -ge 10 && ${V2_PIN:0:${#V2_SOURCE}} == "$V2_SOURCE" ]] || \
    fail "v2 DuckDB source $V2_SOURCE does not match pinned revision $V2_PIN"
V1_SOURCE=$("$TASKSET" -c "$CPU_LIST" "$V1_CLI" -no-init -batch -list -noheader -c 'SELECT source_id FROM pragma_version()' 2>"$OUT_DIR/v1-version-check.log") || fail 'could not query the v1 DuckDB source revision'
V1_SOURCE=${V1_SOURCE//$'\n'/}
V1_SOURCE=${V1_SOURCE//$'\r'/}
[[ -n $V1_SOURCE ]] || fail 'v1 DuckDB returned an empty source revision'

sql_quote() {
    local value=${1//\'/\'\'}
    printf "'%s'" "$value"
}

regions_v1='SELECT seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name FROM model_regions ORDER BY seq_region'
transcripts_v1='SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM model_transcripts'
exons_v1='SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM model_transcripts, unnest(exons) u(e)'
mirna_v1='SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM model_transcripts, unnest(mature_mirna_regions) u(x)'
peptides_v1='SELECT transcript_index, x.protein_position, x.alternate_amino_acid FROM model_transcripts, unnest(peptide_edits) u(x)'
regions_v2=${regions_v1/ FROM model_regions/ FROM m.model_regions}
transcripts_v2=${transcripts_v1/ FROM model_transcripts/ FROM m.model_transcripts}
exons_v2=${exons_v1/ FROM model_transcripts/ FROM m.model_transcripts}
mirna_v2=${mirna_v1/ FROM model_transcripts/ FROM m.model_transcripts}
peptides_v2=${peptides_v1/ FROM model_transcripts/ FROM m.model_transcripts}
calls_query="SELECT * FROM duckvep_coding_calls('hap', $(sql_quote "$VCF"))"

model_builder="duckvep_model_load_sql('hap', $(sql_quote "$regions_v2"), $(sql_quote "$transcripts_v2"), $(sql_quote "$exons_v2"), {mature_mirna_query: $(sql_quote "$mirna_v2"), peptide_edit_query: $(sql_quote "$peptides_v2"), transcript_coverage_complete: true})"
hap_builder="duckvep_haplotype_load_sql($(sql_quote "$calls_query"), 'hap', 'hg002_job', {max_alignment_cells: 268435456, workspace_limit: 1073741824})"

builder_statements() {
    local expression=$1 destination=$2 builder_log=$3 json
    json=$("$TASKSET" -c "$CPU_LIST" "$V2_CLI" -unsigned -no-init -batch -bail -list -noheader \
        -c "LOAD $(sql_quote "$V2_EXTENSION"); SELECT CAST(to_json($expression) AS VARCHAR);" \
        2>"$builder_log") || return 1
    printf '%s' "$json" | "$JQ" -er '
      if type == "array" and length > 0 and all(.[]; type == "string" and length > 0)
      then .[] + ";" else error("builder did not return SQL statements") end
    ' > "$destination"
}

builder_statements "$model_builder" "$OUT_DIR/v2-model-statements.sql" "$OUT_DIR/v2-model-builder.log" || \
    fail 'v2 model-load builder failed; see v2-model-builder.log'
builder_statements "$hap_builder" "$OUT_DIR/v2-haplotype-statements.sql" "$OUT_DIR/v2-haplotype-builder.log" || \
    fail 'v2 haplotype builder failed; see v2-haplotype-builder.log'

v1_regions_q=$(sql_quote "$regions_v1")
v1_transcripts_q=$(sql_quote "$transcripts_v1")
v1_exons_q=$(sql_quote "$exons_v1")
v1_mirna_q=$(sql_quote "$mirna_v1")
v1_peptides_q=$(sql_quote "$peptides_v1")
{
    printf 'SET threads = %s;\nSET memory_limit = %s;\nSET temp_directory = %s;\n' \
        "$THREADS" "$(sql_quote "$MEMORY_LIMIT")" "$(sql_quote "$OUT_DIR/v1-tmp")"
    printf 'LOAD %s;\n' "$(sql_quote "$V1_EXTENSION")"
    printf '.print === v1 model load ===\n.timer on\n'
    printf "SELECT loaded FROM duckvep_model_load('hap', %s, %s, %s, mature_mirna_query := %s, peptide_edit_query := %s, transcript_coverage_complete := true);\n" \
        "$v1_regions_q" "$v1_transcripts_q" "$v1_exons_q" "$v1_mirna_q" "$v1_peptides_q"
    printf "COPY (SELECT _duckvep_model_fingerprint('hap') AS fingerprint) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v1-model-fingerprint.csv")"
    printf "COPY (SELECT * FROM duckvep_native_budget() ORDER BY owner) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v1-budget-after-load.csv")"
    printf "SELECT CASE WHEN bool_and(current_bytes <= limit_bytes AND high_water_bytes <= limit_bytes) THEN true ELSE error('native budget exceeded') END FROM duckvep_native_budget();\n"
    printf '.print === v1 prediction and Parquet write ===\n'
    printf 'COPY (%s) TO %s (FORMAT PARQUET);\n' \
        "SELECT * FROM duckvep_haplotypes($(sql_quote "$calls_query"), 'hap', max_alignment_cells := 268435456, workspace_limit := 1073741824)" \
        "$(sql_quote "$OUT_DIR/v1.parquet")"
    printf "COPY (SELECT * FROM duckvep_native_budget() ORDER BY owner) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v1-budget-final.csv")"
    printf "SELECT CASE WHEN bool_and(current_bytes <= limit_bytes AND high_water_bytes <= limit_bytes) THEN true ELSE error('native budget exceeded') END FROM duckvep_native_budget();\n.timer off\n"
} > "$OUT_DIR/v1.sql"

{
    printf 'SET threads = %s;\nSET memory_limit = %s;\nSET temp_directory = %s;\n' \
        "$THREADS" "$(sql_quote "$MEMORY_LIMIT")" "$(sql_quote "$OUT_DIR/v2-tmp")"
    printf 'LOAD %s;\nATTACH %s AS m (READ_ONLY);\n' \
        "$(sql_quote "$V2_EXTENSION")" "$(sql_quote "$MODEL")"
    printf '.print === v2 model staging and publish ===\n.timer on\n'
    cat "$OUT_DIR/v2-model-statements.sql"
    printf "COPY (SELECT _duckvep_model_fingerprint('hap') AS fingerprint) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v2-model-fingerprint.csv")"
    printf "COPY (SELECT * FROM duckvep_native_budget() ORDER BY owner) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v2-budget-after-load.csv")"
    printf "SELECT CASE WHEN bool_and(current_bytes <= limit_bytes AND high_water_bytes <= limit_bytes) THEN true ELSE error('native budget exceeded') END FROM duckvep_native_budget();\n"
    printf '.print === v2 coding_calls capture ===\n'
    cat "$OUT_DIR/v2-haplotype-statements.sql"
    printf "COPY (SELECT * FROM duckvep_native_budget() ORDER BY owner) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v2-budget-after-capture.csv")"
    printf "SELECT CASE WHEN bool_and(current_bytes <= limit_bytes AND high_water_bytes <= limit_bytes) THEN true ELSE error('native budget exceeded') END FROM duckvep_native_budget();\n"
    printf '.print === v2 haplotype scan and Parquet write ===\n'
    printf 'COPY (SELECT * FROM duckvep_haplotype_scan(\047hg002_job\047)) TO %s (FORMAT PARQUET);\n' \
        "$(sql_quote "$OUT_DIR/v2.parquet")"
    printf "COPY (SELECT * FROM duckvep_native_budget() ORDER BY owner) TO %s (FORMAT CSV, HEADER);\n" \
        "$(sql_quote "$OUT_DIR/v2-budget-final.csv")"
    printf "SELECT CASE WHEN bool_and(current_bytes <= limit_bytes AND high_water_bytes <= limit_bytes) THEN true ELSE error('native budget exceeded') END FROM duckvep_native_budget();\n.timer off\n"
} > "$OUT_DIR/v2.sql"

run_timed() {
    local label=$1 cli=$2
    shift 2
    if ! "$TIME_BIN" -f 'wall_seconds=%e\nmax_rss_kib=%M' -o "$OUT_DIR/$label.resources.txt" \
        "$TASKSET" -c "$CPU_LIST" "$cli" "$@" >"$OUT_DIR/$label.log" 2>&1; then
        printf '%s failed; inspect %s\n' "$label" "$OUT_DIR/$label.log" >&2
        return 1
    fi
}

run_timed v1 "$V1_CLI" -unsigned -readonly -no-init -batch -bail -column -noheader -f "$OUT_DIR/v1.sql" "$MODEL" || exit 1
run_timed v2 "$V2_CLI" -unsigned -no-init -batch -bail -column -noheader -f "$OUT_DIR/v2.sql" :memory: || exit 1
[[ -s $OUT_DIR/v1.parquet ]] || fail 'v1 produced no Parquet output'
[[ -s $OUT_DIR/v2.parquet ]] || fail 'v2 produced no Parquet output'

compare_files() {
    local label=$1 left=$2 right=$3 left_q right_q schema_sql schema_csv rows_sql rows_csv schema_diff row_line
    left_q=$(sql_quote "$left")
    right_q=$(sql_quote "$right")
    schema_sql="$OUT_DIR/$label-schema.sql"
    schema_csv="$OUT_DIR/$label-schema.csv"
    rows_sql="$OUT_DIR/$label-rows.sql"
    rows_csv="$OUT_DIR/$label-rows.csv"
    {
        printf 'SET threads = %s;\nSET memory_limit = %s;\n' "$THREADS" "$(sql_quote "$MEMORY_LIMIT")"
        printf "COPY (WITH left_schema AS (SELECT row_number() OVER () AS ordinal, column_name, column_type FROM (DESCRIBE SELECT * FROM read_parquet(%s))), right_schema AS (SELECT row_number() OVER () AS ordinal, column_name, column_type FROM (DESCRIBE SELECT * FROM read_parquet(%s))) SELECT count(*) AS schema_differences FROM ((SELECT * FROM left_schema EXCEPT ALL SELECT * FROM right_schema) UNION ALL (SELECT * FROM right_schema EXCEPT ALL SELECT * FROM left_schema))) TO %s (FORMAT CSV, HEADER);\n" \
            "$left_q" "$right_q" "$(sql_quote "$schema_csv")"
    } > "$schema_sql"
    run_timed "$label-schema" "$V1_CLI" -no-init -batch -bail -column -noheader -f "$schema_sql" :memory: || return 1
    schema_diff=$(tail -n 1 "$schema_csv" | tr -d '\r')
    if [[ $schema_diff != 0 ]]; then
        printf '%s schema differs (%s differing schema rows)\n' "$label" "$schema_diff" >&2
        return 1
    fi
    {
        printf 'SET threads = %s;\nSET memory_limit = %s;\n' "$THREADS" "$(sql_quote "$MEMORY_LIMIT")"
        printf "COPY (SELECT (SELECT count(*) FROM read_parquet(%s)) AS left_rows, (SELECT count(*) FROM read_parquet(%s)) AS right_rows, (SELECT count(*) FROM ((SELECT * FROM read_parquet(%s) EXCEPT ALL SELECT * FROM read_parquet(%s)) UNION ALL (SELECT * FROM read_parquet(%s) EXCEPT ALL SELECT * FROM read_parquet(%s)))) AS row_differences) TO %s (FORMAT CSV, HEADER);\n" \
            "$left_q" "$right_q" "$left_q" "$right_q" "$right_q" "$left_q" "$(sql_quote "$rows_csv")"
    } > "$rows_sql"
    run_timed "$label-rows" "$V1_CLI" -no-init -batch -bail -column -noheader -f "$rows_sql" :memory: || return 1
    row_line=$(tail -n 1 "$rows_csv" | tr -d '\r')
    [[ $row_line =~ ^[0-9]+,[0-9]+,0$ ]] || {
        printf '%s exact row multiset differs: %s\n' "$label" "$row_line" >&2
        return 1
    }
}

exact_ok=1
if ! compare_files v1-v2 "$OUT_DIR/v1.parquet" "$OUT_DIR/v2.parquet"; then
    exact_ok=0
fi
reference_ok=1
if [[ -n $V1_REFERENCE ]]; then
    if ! compare_files v1-reference "$OUT_DIR/v1.parquet" "$V1_REFERENCE"; then
        reference_ok=0
    fi
fi
v1_fingerprint=$(tail -n 1 "$OUT_DIR/v1-model-fingerprint.csv" | tr -d '\r')
v2_fingerprint=$(tail -n 1 "$OUT_DIR/v2-model-fingerprint.csv" | tr -d '\r')
fingerprint_ok=1
if [[ -z $v1_fingerprint || $v1_fingerprint != "$v2_fingerprint" ]]; then
    fingerprint_ok=0
    printf 'model fingerprints differ: v1=%s v2=%s\n' "$v1_fingerprint" "$v2_fingerprint" >&2
fi

status=PASS
if (( ! exact_ok || ! reference_ok || ! fingerprint_ok )); then
    status=FAIL
fi
{
    printf 'status\t%s\n' "$status"
    printf 'v1_duckdb_source\t%s\n' "$V1_SOURCE"
    printf 'v2_duckdb_source\t%s\n' "$V2_SOURCE"
    printf 'v2_pinned_revision\t%s\n' "$V2_PIN"
    printf 'v1_extension\t%s\n' "$V1_EXTENSION"
    printf 'v2_extension\t%s\n' "$V2_EXTENSION"
    printf 'vcf\t%s\nmodel\t%s\n' "$VCF" "$MODEL"
    printf 'cpu_list\t%s\nthreads\t%s\nmemory_limit\t%s\n' "$CPU_LIST" "$THREADS" "$MEMORY_LIMIT"
    printf 'v1_parquet\t%s\nv2_parquet\t%s\n' "$OUT_DIR/v1.parquet" "$OUT_DIR/v2.parquet"
    printf 'v1_model_fingerprint\t%s\nv2_model_fingerprint\t%s\n' "$v1_fingerprint" "$v2_fingerprint"
    printf 'v1_v2_exact_multiset\t%s\n' "$([[ $exact_ok == 1 ]] && printf PASS || printf FAIL)"
    printf 'v1_reference\t%s\n' "${V1_REFERENCE:-not supplied}"
    printf 'v1_reference_exact_multiset\t%s\n' "$([[ $reference_ok == 1 ]] && printf PASS || printf FAIL)"
    printf 'v1_max_rss_kib\t%s\nv2_max_rss_kib\t%s\n' \
        "$(awk -F= '$1 == "max_rss_kib" {print $2}' "$OUT_DIR/v1.resources.txt")" \
        "$(awk -F= '$1 == "max_rss_kib" {print $2}' "$OUT_DIR/v2.resources.txt")"
    printf 'v1_output_bytes\t%s\nv2_output_bytes\t%s\n' \
        "$(wc -c < "$OUT_DIR/v1.parquet" | tr -d ' ')" "$(wc -c < "$OUT_DIR/v2.parquet" | tr -d ' ')"
} > "$OUT_DIR/receipt.tsv"

printf '%s\n' "HG002 v1-v2 check: $status" "Receipts: $OUT_DIR/receipt.tsv"
[[ $status == PASS ]]
