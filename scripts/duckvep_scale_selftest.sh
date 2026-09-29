#!/usr/bin/env bash
# shellcheck disable=SC2317  # helpers are called through check
# Capacity-failure self-test for the scale runner (started by `duckvep_scale_run.sh --self-test`).
# Four small runs, each of which must end with an explicit classification, never a success:
#   1. native budget 64 MiB          -> capacity_error, nothing published, budget back to baseline,
#                                       and the same connection annotates again after a restored budget
#   2. DuckDB memory_limit 100MB and a 1 MiB spill quota -> capacity_error with DuckDB's error text
#   3. cgroup MemoryMax just above the process baseline   -> failed (oom_kill), not ok
# Usage: duckvep_scale_selftest.sh RUNNER [--self-test] [--out DIR] [--panel P] [--extension PATH] [--model PATH]
set -euo pipefail
RUNNER="$1"; shift
OUT="" PANEL="1M" EXTRA=() TEMP_QUOTA="${DUCKVEP_SELFTEST_TEMP:-8GiB}" MIN_FREE="${DUCKVEP_MIN_FREE_GIB:-50}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --self-test) shift ;;
    --out) OUT="$2"; shift 2 ;;
    --panel) PANEL="$2"; shift 2 ;;
    --extension|--model|--cgroup) EXTRA+=("$1" "$2"); shift 2 ;;
    --allow-unenforced) EXTRA+=("$1"); shift ;;
    *) echo "self-test: unknown argument $1" >&2; exit 2 ;;
  esac
done
OUT="${OUT:-$(mktemp -d /tmp/duckvep-scale-selftest.XXXXXX)}"
mkdir -p "$OUT"
fail=0
field() { Rscript -e 'r <- read.csv(commandArgs(TRUE)[1], colClasses = "character", na.strings = ""); v <- r[[commandArgs(TRUE)[2]]][1]; cat(if (is.na(v)) "" else v)' "$1" "$2"; }
result_field() { Rscript -e 'r <- read.delim(commandArgs(TRUE)[1], colClasses = "character", na.strings = "", quote = ""); v <- r[[commandArgs(TRUE)[2]]][1]; cat(if (is.na(v)) "" else v)' "$1" "$2"; }
check() {  # name, detail shown on failure, then the command that must succeed
  local name="$1" detail="$2"
  shift 2
  if "$@"; then echo "PASS  $name"; else echo "FAIL  $name: $detail"; fail=1; fi
}
within() { [[ -n "$1" && -n "$2" && "$2" -le $(($1 + $3)) ]]; }
contains() { [[ "$1" == *"$2"* ]]; }
nonempty_positive() { [[ -n "$1" && "$1" -gt 0 ]]; }
run() {  # dir, then runner flags; the runner exits 3 when a job is not ok, which is expected
  local dir="$1"; shift
  set +e
  "$RUNNER" --jobs 1 --panel "$PANEL" --out "$dir" --threads-per-job "${DUCKVEP_SELFTEST_THREADS:-6}" --max-temp "$TEMP_QUOTA" --min-free-gib "$MIN_FREE" "${EXTRA[@]}" "$@" >"$dir.log" 2>&1
  local code=$?
  set -e
  echo "$code" >"$dir.exit"
  return 0
}
no_output() { ! find "$1" -name 'output-*' | grep -q .; }
exits() { [[ "$(<"$1.exit")" == "$2" ]]; }

# 1. Native budget of 64 MiB: the model cannot load; nothing is published; the connection recovers.
run "$OUT/native-64mib" --limit-rows 20000 --modes compact --native-budget-mib 64 --retry-after-capacity
d="$OUT/native-64mib"
o="$(field "$d/receipt.csv" outcome)"; r="$(field "$d/receipt.csv" reason)"
check "native budget 64 MiB: outcome capacity_error" "outcome=$o reason=$r" test "$o" = capacity_error
check "native budget 64 MiB: explicit capacity error text" "$r" contains "$r" "capacity error"
check "native budget 64 MiB: nothing published" "output file exists" no_output "$d"
check "native budget 64 MiB: runner exit status is 3" "exit=$(<"$d.exit")" exits "$d" 3
base="$(result_field "$d/job-1/result.tsv" native_total_baseline)"
after="$(result_field "$d/job-1/result.tsv" native_total_after_error)"
check "native budget 64 MiB: no model or index bytes left ($base -> $after; the 512-byte registry stays)" "baseline=$base after=$after" within "$base" "$after" 512
retry="$(result_field "$d/job-1/result.tsv" retry_ok)"
check "connection reusable after a capacity error" "retry_ok=$retry" test "$retry" = TRUE

# 2. DuckDB memory limit 100MB with a 1 MiB spill quota: DuckDB's own explicit error.
run "$OUT/duckdb-100mb" --modes complete17 --duckdb-memory-limit 100MB --max-temp 1MiB
d="$OUT/duckdb-100mb"
o="$(field "$d/receipt.csv" outcome)"; r="$(field "$d/receipt.csv" reason)"
check "DuckDB 100MB + 1MiB spill: outcome capacity_error" "outcome=$o reason=$r" test "$o" = capacity_error
check "DuckDB 100MB + 1MiB spill: DuckDB error text" "$r" contains "$r" "Out of Memory"
check "DuckDB 100MB + 1MiB spill: nothing published" "output file exists" no_output "$d"

# 3. cgroup MemoryMax just above the baseline (a bare R + DuckDB process), far below the model load.
if [[ " ${EXTRA[*]} " == *" none "* ]]; then
  echo "SKIP  cgroup MemoryMax: --cgroup none does not enforce a process ceiling"
else
  run "$OUT/cgroup-oom" --limit-rows 20000 --modes compact --memory-max 700M
  d="$OUT/cgroup-oom"
  o="$(field "$d/receipt.csv" outcome)"; k="$(field "$d/receipt.csv" oom_kill)"; r="$(field "$d/receipt.csv" reason)"
  check "cgroup MemoryMax 700M: classified failed, not ok" "outcome=$o reason=$r" test "$o" = failed
  check "cgroup MemoryMax 700M: oom_kill counted" "oom_kill=$k" nonempty_positive "$k"
  check "cgroup MemoryMax 700M: nothing published" "output file exists" no_output "$d"
fi
echo "self-test artifacts: $OUT"
if [[ "$fail" == 0 ]]; then echo "self-test: all capacity checks passed"; else echo "self-test: FAILED"; fi
exit "$fail"
