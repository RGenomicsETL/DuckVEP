#!/usr/bin/env bash
# DuckVEP scale runner: N concurrent annotation jobs over one gnomAD panel, each in its own
# cgroup v2 scope with the agreed ceilings, then one receipt and one summary.
#
#   scripts/duckvep_scale_run.sh --jobs N --panel 5M --out DIR [options]
#   scripts/duckvep_scale_run.sh --self-test [--out DIR]
#
# Ceilings enforced per job (never advisory):
#   process   cgroup v2 memory.max 16 GiB, memory.swap.max 0 (an OOM kill is "failed")
#   DuckDB    memory_limit 8GB, job-local temp_directory, max_temp_directory_size 32GiB
#   native    4 GiB budget, 6 workers x 128 MiB scratch + 256 MiB emit, threads 6
# Each job ends ok, capacity_error (explicit) or failed (with the reason).
set -euo pipefail

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
REPO="$(dirname "$SCRIPT_DIR")"
GIB=$((1024 * 1024 * 1024))
PANEL_ROOT="${DUCKVEP_GNOMAD_ROOT:-/root/duckvep/data/gnomad-v4.1}/panels"

usage() {
  cat <<'EOF'
usage: duckvep_scale_run.sh --jobs N --panel PANEL --out DIR [options]

  --jobs N                 concurrent jobs (default 1)
  --panel P                1M, 5M, 25M, 100M, exome2M, structural, or a Parquet path
  --out DIR                output directory (receipt, summary, per-job files); must be empty or new
  --model PATH             GRCh38 model DuckDB file
                           (default $DUCKVEP_SCALE_MODEL or /root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb)
  --extension PATH         duckvep.duckdb_extension (default build/release/extension/duckvep/duckvep.duckdb_extension)
  --threads-per-job T      DuckDB threads per job (default 6)
  --cgroup MODE            systemd | manual | none | auto (default auto: systemd, then manual)
  --allow-unenforced       required with --cgroup none; the receipt says "ceilings not enforced"
  --min-free-gib G         disk that must stay free after the spill quotas (default 10)
  --oversubscribe          run although the host cannot hold N jobs at the ceilings
  --modes LIST             compact,complete17 (default both)
  --spill-dir DIR          parent for per-job spill directories (default OUT)
  --keep-output            keep each job's output Parquet (default: checksum, then delete)
  --regulation             load the resident regulatory and motif feature intervals (default, the production configuration)
  --no-regulation          annotate without them
  --retry-after-capacity   test hook: after a load capacity error, restore the ceiling and annotate again
  --limit-rows N           annotate only the first N panel rows (tests and smoke runs)
  --self-test              run the capacity-failure tests and exit

Test overrides for the ceilings (defaults are the agreed values):
  --memory-max SIZE (16G)  --duckdb-memory-limit (8GB)  --max-temp (32GiB)
  --native-budget-mib (4096)  --workers (6)  --scratch-mib (128)  --emit-mib (256, the extension default; 64 fails complete-17 on gene-dense panels)
Exit status: 0 all jobs ok; 2 usage or refused by the pre-flight; 3 some job was not ok.
EOF
}
die() { echo "duckvep_scale_run: $*" >&2; exit 2; }

# ---------------------------------------------------------------- job wrapper (runs inside the cgroup)
# Reads the cgroup counters after the job process exits; the wrapper survives a kernel OOM kill of
# the job, so an OOM is recorded rather than lost.
if [[ "${1:-}" == "__job" ]]; then
  shift
  jobdir="$1"; shift
  mkdir -p "$jobdir"
  cgroup_mode="${SCALE_CGROUP_MODE:?}"
  if [[ "$cgroup_mode" == manual ]]; then
    echo $$ >"${SCALE_MANUAL_CGROUP:?}/cgroup.procs"
    cgdir="$SCALE_MANUAL_CGROUP"
  elif [[ "$cgroup_mode" == systemd ]]; then
    cgdir="/sys/fs/cgroup$(sed -n 's/^0:://p' /proc/self/cgroup)"
  else
    cgdir=""
  fi
  spill="$jobdir/spill"
  [[ -n "${SCALE_SPILL_DIR:-}" ]] && spill="$SCALE_SPILL_DIR"
  mkdir -p "$spill"
  spill_max_file="$jobdir/spill-peak.txt"
  echo 0 >"$spill_max_file"
  (
    peak=0
    while :; do
      now="$(du -sb "$spill" 2>/dev/null | cut -f1 || echo 0)"
      if [[ "${now:-0}" -gt "$peak" ]]; then peak="$now"; echo "$peak" >"$spill_max_file"; fi
      sleep 0.25
    done
  ) &
  sampler=$!
  start="$(date +%s.%N)"
  set +e
  Rscript "$SCRIPT_DIR/duckvep_scale_job.R" --out "$jobdir" --temp-dir "$spill" "$@" >"$jobdir/job.log" 2>&1
  code=$?
  set -e
  end="$(date +%s.%N)"
  kill "$sampler" 2>/dev/null || true
  wait "$sampler" 2>/dev/null || true
  peak="" max="" swapmax="" swappeak=""
  declare -A ev=()
  if [[ -n "$cgdir" && -r "$cgdir/memory.peak" ]]; then
    peak="$(<"$cgdir/memory.peak")"
    max="$(<"$cgdir/memory.max")"
    swapmax="$(<"$cgdir/memory.swap.max")"
    swappeak="$(cat "$cgdir/memory.swap.peak" 2>/dev/null || true)"
    while read -r k v; do ev[$k]="$v"; done <"$cgdir/memory.events"
  fi
  {
    printf 'exit_code\twall_s\tmemory_peak\tmemory_max\tmemory_swap_max\tmemory_swap_peak\toom\toom_kill\tevent_max\tevent_high\tspill_peak_bytes\tcgroup\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$code" \
      "$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.3f", b-a}')" "$peak" "$max" "$swapmax" "$swappeak" \
      "${ev[oom]:-}" "${ev[oom_kill]:-}" "${ev[max]:-}" "${ev[high]:-}" "$(<"$spill_max_file")" "$cgdir"
  } >"$jobdir/cgroup.tsv"
  exit 0
fi

# ---------------------------------------------------------------- self-test
if [[ "${1:-}" == "--self-test" ]] || [[ " $* " == *" --self-test "* ]]; then
  exec "$SCRIPT_DIR/duckvep_scale_selftest.sh" "$SCRIPT_PATH" "$@"
fi

# ---------------------------------------------------------------- arguments
JOBS=1 PANEL="" OUT="" MODEL="${DUCKVEP_SCALE_MODEL:-/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb}"
EXTENSION="$REPO/build/release/extension/duckvep/duckvep.duckdb_extension"
THREADS=6 CGROUP=auto ALLOW_UNENFORCED=0 OVERSUBSCRIBE=0 MODES="compact,complete17" SPILL_PARENT="" KEEP=0
MIN_FREE_GIB=10 REGULATION=1 RETRY=0 LIMIT_ROWS="" MEMORY_MAX=16G DUCKDB_LIMIT=8GB MAX_TEMP=32GiB BUDGET_MIB=4096 WORKERS=6 SCRATCH_MIB=128 EMIT_MIB=256
while [[ $# -gt 0 ]]; do
  need() { [[ $# -ge 2 ]] || die "missing value for $1"; }
  case "$1" in
    --jobs) need "$@"; JOBS="$2"; shift 2 ;;
    --panel) need "$@"; PANEL="$2"; shift 2 ;;
    --out) need "$@"; OUT="$2"; shift 2 ;;
    --model) need "$@"; MODEL="$2"; shift 2 ;;
    --extension) need "$@"; EXTENSION="$2"; shift 2 ;;
    --threads-per-job) need "$@"; THREADS="$2"; shift 2 ;;
    --cgroup) need "$@"; CGROUP="$2"; shift 2 ;;
    --allow-unenforced) ALLOW_UNENFORCED=1; shift ;;
    --oversubscribe) OVERSUBSCRIBE=1; shift ;;
    --min-free-gib) need "$@"; MIN_FREE_GIB="$2"; shift 2 ;;
    --modes) need "$@"; MODES="$2"; shift 2 ;;
    --spill-dir) need "$@"; SPILL_PARENT="$2"; shift 2 ;;
    --keep-output) KEEP=1; shift ;;
    --retry-after-capacity) RETRY=1; shift ;;
    --regulation) REGULATION=1; shift ;;
    --no-regulation) REGULATION=0; shift ;;
    --limit-rows) need "$@"; LIMIT_ROWS="$2"; shift 2 ;;
    --memory-max) need "$@"; MEMORY_MAX="$2"; shift 2 ;;
    --duckdb-memory-limit) need "$@"; DUCKDB_LIMIT="$2"; shift 2 ;;
    --max-temp) need "$@"; MAX_TEMP="$2"; shift 2 ;;
    --native-budget-mib) need "$@"; BUDGET_MIB="$2"; shift 2 ;;
    --workers) need "$@"; WORKERS="$2"; shift 2 ;;
    --scratch-mib) need "$@"; SCRATCH_MIB="$2"; shift 2 ;;
    --emit-mib) need "$@"; EMIT_MIB="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown argument: $1" ;;
  esac
done
[[ -n "$PANEL" && -n "$OUT" ]] || { usage >&2; die "--panel and --out are required"; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads-per-job must be a positive integer"
case "$CGROUP" in auto|systemd|manual|none) ;; *) die "--cgroup: systemd, manual, none or auto" ;; esac

size_bytes() {  # 16G, 512M, 8GB, 32GiB, plain bytes
  local v="$1" n unit
  n="${v%%[!0-9]*}"; unit="${v#"$n"}"
  [[ -n "$n" ]] || die "bad size: $v"
  case "${unit^^}" in
    ""|B) echo "$n" ;; K|KB|KIB) echo $((n * 1024)) ;; M|MB|MIB) echo $((n * 1024 * 1024)) ;;
    G|GB|GIB) echo $((n * GIB)) ;; T|TB|TIB) echo $((n * 1024 * GIB)) ;; *) die "bad size: $v" ;;
  esac
}

# Panel names resolve through the symlinks that scripts/gnomad_build_panels.sh maintains.
resolve_panel() {
  case "$1" in
    1M) echo "$PANEL_ROOT/genomes-1M.parquet" ;; 5M) echo "$PANEL_ROOT/genomes-5M.parquet" ;;
    25M) echo "$PANEL_ROOT/genomes-25M.parquet" ;; 100M) echo "$PANEL_ROOT/genomes-100M.parquet" ;;
    exome2M) echo "$PANEL_ROOT/exomes-2M.parquet" ;; structural) echo "$PANEL_ROOT/structural-controls.parquet" ;;
    *) echo "$1" ;;
  esac
}
PANEL_FILE="$(resolve_panel "$PANEL")"
[[ -e "$PANEL_FILE" ]] || die "panel $PANEL not found ($PANEL_FILE); build it with scripts/gnomad_build_panels.sh (a 100M panel needs every genome shard)"
PANEL_FILE="$(readlink -f "$PANEL_FILE")"
[[ -r "$MODEL" ]] || die "model not readable: $MODEL"
[[ -r "$EXTENSION" ]] || die "extension not found: $EXTENSION (run make release or pass --extension)"
command -v Rscript >/dev/null || die "Rscript is required"

# ---------------------------------------------------------------- pre-flight
MEMORY_MAX_BYTES="$(size_bytes "$MEMORY_MAX")"
MAX_TEMP_BYTES="$(size_bytes "$MAX_TEMP")"
mem_total=$(awk '/^MemTotal:/{print $2*1024}' /proc/meminfo)
cores=$(nproc)
mkdir -p "$OUT"
if [[ -n "$(ls -A "$OUT" 2>/dev/null)" ]]; then die "--out $OUT is not empty"; fi
SPILL_PARENT="${SPILL_PARENT:-$OUT}"
mkdir -p "$SPILL_PARENT"
disk_free=$(df -PB1 "$SPILL_PARENT" | awk 'NR==2{print $4}')
have_v2=0; [[ "$(stat -fc %T /sys/fs/cgroup 2>/dev/null)" == cgroup2fs ]] && have_v2=1
have_systemd_run=0
if command -v systemd-run >/dev/null && [[ -d /run/systemd/system ]] && [[ "$(id -u)" == 0 || -n "${XDG_RUNTIME_DIR:-}" ]]; then have_systemd_run=1; fi
if [[ "$CGROUP" == auto ]]; then
  if (( have_v2 && have_systemd_run )); then CGROUP=systemd
  elif (( have_v2 )); then CGROUP=manual
  else CGROUP=none; fi
fi
problems=()
if [[ "$CGROUP" == none ]]; then
  (( ALLOW_UNENFORCED )) || problems+=("--cgroup none needs --allow-unenforced; the process ceilings would not be enforced")
else
  (( have_v2 )) || problems+=("cgroup v2 is not mounted at /sys/fs/cgroup")
  [[ "$CGROUP" == systemd ]] && ! (( have_systemd_run )) && problems+=("systemd-run is not usable here (try --cgroup manual)")
  [[ "$CGROUP" == manual && "$(id -u)" != 0 ]] && problems+=("--cgroup manual needs root to write cgroup.procs")
fi
if (( ! OVERSUBSCRIBE )); then
  need_mem=$((JOBS * MEMORY_MAX_BYTES + 4 * GIB))
  (( mem_total >= need_mem )) || problems+=("$JOBS jobs x $MEMORY_MAX need $((need_mem / GIB)) GiB with 4 GiB for the OS; host RAM is $((mem_total / GIB)) GiB")
  (( JOBS * THREADS <= cores )) || problems+=("$JOBS jobs x $THREADS threads = $((JOBS * THREADS)) exceed $cores cores")
  need_disk=$((JOBS * MAX_TEMP_BYTES + MIN_FREE_GIB * GIB))
  (( disk_free >= need_disk )) || problems+=("spill quotas need $((need_disk / GIB)) GiB ($JOBS x $MAX_TEMP plus $MIN_FREE_GIB GiB kept free) but $SPILL_PARENT has $((disk_free / GIB)) GiB free")
fi
if ((${#problems[@]})); then
  echo "duckvep_scale_run: refusing to start:" >&2
  printf '  - %s\n' "${problems[@]}" >&2
  echo "  (--oversubscribe overrides the capacity checks; it does not override enforcement)" >&2
  exit 2
fi
ENFORCED=yes; [[ "$CGROUP" == none ]] && ENFORCED="no: ceilings not enforced"

# ---------------------------------------------------------------- run
RUN_ID="scale-$(date -u +%Y%m%dT%H%M%SZ)-$$"
{
  printf 'key\tvalue\n'
  printf 'run_id\t%s\nstarted_utc\t%s\nhost\t%s\nkernel\t%s\ncores\t%s\nmem_total_bytes\t%s\n' \
    "$RUN_ID" "$(date -u +%FT%TZ)" "$(hostname)" "$(uname -r)" "$cores" "$mem_total"
  printf 'jobs\t%s\npanel\t%s\npanel_file\t%s\npanel_bytes\t%s\npanel_sha256\t%s\n' \
    "$JOBS" "$PANEL" "$PANEL_FILE" "$(stat -c %s "$PANEL_FILE")" "$(sha256sum "$PANEL_FILE" | cut -c1-64)"
  printf 'model\t%s\nmodel_sha256\t%s\nextension\t%s\nextension_sha256\t%s\n' "$MODEL" \
    "$(sha256sum "$MODEL" | cut -c1-64)" "$EXTENSION" "$(sha256sum "$EXTENSION" | cut -c1-64)"
  printf 'git_revision\t%s\ngit_dirty\t%s\n' "$(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)" \
    "$(git -C "$REPO" status --porcelain 2>/dev/null | wc -l)"
  printf 'cgroup_mode\t%s\nenforced\t%s\noversubscribed\t%s\n' "$CGROUP" "$ENFORCED" "$OVERSUBSCRIBE"
  printf 'regulation\t%s\n' "$REGULATION"
  printf 'min_free_gib\t%s\n' "$MIN_FREE_GIB"
  printf 'threads_per_job\t%s\nmodes\t%s\nlimit_rows\t%s\n' "$THREADS" "$MODES" "${LIMIT_ROWS:-all}"
  printf 'memory_max_bytes\t%s\nmemory_swap_max\t0\nduckdb_memory_limit\t%s\nmax_temp_bytes\t%s\n' \
    "$MEMORY_MAX_BYTES" "$DUCKDB_LIMIT" "$MAX_TEMP_BYTES"
  printf 'native_budget_mib\t%s\nworkers\t%s\nscratch_mib\t%s\nemit_mib\t%s\n' "$BUDGET_MIB" "$WORKERS" "$SCRATCH_MIB" "$EMIT_MIB"
  printf 'loadavg_at_start\t%s\n' "$(cut -d' ' -f1-3 /proc/loadavg)"
} >"$OUT/run.tsv"

job_args=(--panel "$PANEL_FILE" --modes "$MODES" --threads "$THREADS" --model "$MODEL"
  --extension "$EXTENSION" --memory-limit "$DUCKDB_LIMIT" --max-temp "$MAX_TEMP"
  --native-budget-mib "$BUDGET_MIB" --workers "$WORKERS" --scratch-mib "$SCRATCH_MIB" --emit-mib "$EMIT_MIB")
(( KEEP )) && job_args+=(--keep-output)
(( RETRY )) && job_args+=(--retry-after-capacity)
if (( REGULATION )); then job_args+=(--regulation); else job_args+=(--no-regulation); fi
[[ -n "$LIMIT_ROWS" ]] && job_args+=(--limit-rows "$LIMIT_ROWS")

pids=()
manual_dirs=()
# shellcheck disable=SC2317  # invoked through the EXIT trap
cleanup() {
  local d
  for d in "${manual_dirs[@]:-}"; do
    if [[ -n "$d" && -d "$d" ]]; then rmdir "$d" 2>/dev/null || true; fi
  done
}
trap cleanup EXIT
wall_start=$(date +%s.%N)
for ((i = 1; i <= JOBS; i++)); do
  jobdir="$OUT/job-$i"
  mkdir -p "$jobdir"
  export SCALE_CGROUP_MODE="$CGROUP"
  export SCALE_SPILL_DIR="$SPILL_PARENT/spill-job-$i"
  case "$CGROUP" in
    systemd)
      systemd-run --scope -q -p "MemoryMax=$MEMORY_MAX_BYTES" -p MemorySwapMax=0 -p OOMPolicy=continue \
        --unit="duckvep-$RUN_ID-job$i" "$SCRIPT_PATH" __job "$jobdir" --job-id "$i" "${job_args[@]}" &
      ;;
    manual)
      cg="/sys/fs/cgroup/duckvep-$RUN_ID-job$i"
      echo +memory >/sys/fs/cgroup/cgroup.subtree_control 2>/dev/null || true
      mkdir "$cg"; manual_dirs+=("$cg")
      echo "$MEMORY_MAX_BYTES" >"$cg/memory.max"; echo 0 >"$cg/memory.swap.max"
      SCALE_MANUAL_CGROUP="$cg" "$SCRIPT_PATH" __job "$jobdir" --job-id "$i" "${job_args[@]}" &
      ;;
    none)
      "$SCRIPT_PATH" __job "$jobdir" --job-id "$i" "${job_args[@]}" &
      ;;
  esac
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
wall_end=$(date +%s.%N)
printf 'wall_s\t%s\nfinished_utc\t%s\n' "$(awk -v a="$wall_start" -v b="$wall_end" 'BEGIN{printf "%.3f", b-a}')" "$(date -u +%FT%TZ)" >>"$OUT/run.tsv"
for ((i = 1; i <= JOBS; i++)); do rm -rf "$SPILL_PARENT/spill-job-$i" "$OUT/job-$i/spill"; done

status=0
Rscript "$SCRIPT_DIR/duckvep_scale_aggregate.R" "$OUT" || status=$?
echo "receipt: $OUT/receipt.csv"
echo "summary: $OUT/summary.md"
exit "$status"
