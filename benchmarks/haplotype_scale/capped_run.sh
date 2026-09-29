#!/usr/bin/env bash
# Run one command as a fresh process under the qualification caps and append one TSV receipt row.
#
#   capped_run.sh CORE RECEIPT_TSV LABEL -- command [args...]
#
# Waits (60 s polls, up to 30 min) while the 1-minute load average is above DUCKVEP_LOAD_MAX (default 4), then runs
# the command in a systemd scope with MemoryMax=16G and no swap, pinned to CORE. The row records wall seconds, user
# and system seconds, GNU-time peak RSS (KiB), the cgroup memory.peak (bytes), the load average at the start and the
# busy percentage of the sibling hardware thread during the run.
# The command's stdout and stderr go to LABEL.out / LABEL.err beside RECEIPT_TSV. With SPILL_DIR set, the peak size of
# that directory (DuckDB's temp_directory) is sampled every 0.25 s and recorded as spill_peak_bytes.
set -euo pipefail
core=$1; receipt=$2; label=$3; shift 3
[[ ${1:-} == -- ]] && shift
load_max=${DUCKVEP_LOAD_MAX:-4}
waited=0
while :; do
    read -r load1 _ < /proc/loadavg
    if awk -v l="$load1" -v m="$load_max" 'BEGIN{exit !(l<=m)}'; then break; fi
    # After 30 minutes the run proceeds anyway; waited_s and load1_at_start in the receipt say so.
    (( waited >= 1800 )) && { echo "load stayed above $load_max for 30 min; running at load $load1" >&2; break; }
    sleep 60; waited=$((waited + 60))
done
dir=$(dirname "$receipt")
inner='
set -euo pipefail
cg=/sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup)
spill_peak=0
if [ -n "${SPILL_DIR:-}" ]; then
    ( peak=0; while :; do
        now=$(du -sb "$SPILL_DIR" 2>/dev/null | cut -f1); now=${now:-0}
        [ "$now" -gt "$peak" ] && { peak=$now; echo "$peak" > "$SPILL_FILE"; }
        sleep 0.25; done ) &
    sampler=$!
fi
sib_read() { awk -v c="cpu$SIBLING" '"'"'$1==c {t=0; for(i=2;i<=NF;i++) t+=$i; print t-$5-$6, t}'"'"' /proc/stat; }
read -r sb0 st0 < <(sib_read)
set +e
/usr/bin/time -f "%e\t%U\t%S\t%M" -o "$TIMEFILE" taskset -c "$CORE" "$@" >"$OUT" 2>"$ERR"
status=$?
set -e
read -r sb1 st1 < <(sib_read)
sibling=$(awk -v b=$((sb1-sb0)) -v t=$((st1-st0)) '"'"'BEGIN{printf "%.0f", (t>0 ? 100*b/t : 0)}'"'"')
[ -n "${sampler:-}" ] && kill "$sampler" 2>/dev/null || true
[ -f "${SPILL_FILE:-/nonexistent}" ] && spill_peak=$(cat "$SPILL_FILE")
peak=$(cat "$cg/memory.peak" 2>/dev/null || echo NA)
IFS=$(printf "\t") read -r wall user sys rss < <(tail -n 1 "$TIMEFILE")
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$LABEL" "$status" "$wall" "$user" "$sys" "$rss" "$peak" "$LOAD1" "$CORE" "$sibling" "$spill_peak" "$WAITED" >>"$RECEIPT"
exit $status
'
list=$(cat /sys/devices/system/cpu/cpu$core/topology/thread_siblings_list)
SIBLING=$(echo "$list" | awk -F'[,-]' -v c="$core" '{for(i=1;i<=NF;i++) if ($i!=c) {print $i; exit}}'); export SIBLING
export SPILL_DIR=${SPILL_DIR:-} SPILL_FILE=$dir/$label.spill WAITED=$waited
export TIMEFILE=$dir/$label.time OUT=$dir/$label.out ERR=$dir/$label.err CORE=$core RECEIPT=$receipt LABEL=$label LOAD1=$load1
[[ -s $receipt ]] || printf 'label\texit\twall_s\tuser_s\tsys_s\tmax_rss_kib\tcgroup_peak_bytes\tload1_at_start\tcore\tsibling_busy_pct\tspill_peak_bytes\twaited_s\n' >"$receipt"
systemd-run --scope -q -p MemoryMax=16G -p MemorySwapMax=0 --setenv=TIMEFILE --setenv=OUT --setenv=ERR --setenv=CORE --setenv=RECEIPT --setenv=LABEL --setenv=LOAD1 --setenv=SIBLING --setenv=SPILL_DIR --setenv=SPILL_FILE --setenv=WAITED bash -c "$inner" _ "$@"
