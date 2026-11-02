#!/usr/bin/env bash
# A host allocation error must release annotation ownership before the next query.
set -euo pipefail
[[ $(uname -s) == Linux ]] || { echo 'This address-space test requires Linux.' >&2; exit 2; }
root=$(cd "$(dirname "$0")/../.." && pwd)
extension=$(realpath "${1:?usage: test_host_oom.sh /path/to/duckvep.duckdb_extension}")
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
shared=${DUCKDB_SHARED:-$(Rscript -e 'cat(system.file("libs", "duckdb.so", package="duckdb"))')}
[[ -f "$shared" ]] || { echo 'Set DUCKDB_SHARED to the R DuckDB shared library.' >&2; exit 2; }
read -r -a rflags <<< "$(R CMD config --ldflags)"
"${CC:-cc}" -std=c11 -O0 -g -Wall -Wextra -Werror -I "$root/duckdb_capi" \
  "$root/test/scripts/v1_annotation_oom.c" "$shared" "${rflags[@]}" \
  -Wl,-rpath,"$(R RHOME)/lib" -Wl,-rpath,"$(dirname "$shared")" \
  -pthread -o "$work/probe"
# The standalone probe needs the shared-library paths used by the R launcher.
LD_LIBRARY_PATH="$(Rscript -e 'cat(Sys.getenv("LD_LIBRARY_PATH"))')"
export LD_LIBRARY_PATH
export MALLOC_ARENA_MAX=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ulimit -c 0
if ! timeout -k 1 10 "$work/probe" "$extension" -1 10000 recover > "$work/control.log" 2>&1; then
  printf 'Unrestricted annotation control failed:\n' >&2
  grep . "$work/control.log" >&2 || true
  exit 1
fi
grep -q '^annotation=ok' "$work/control.log"
grep -q '^recovery=ok' "$work/control.log"
grep -q '^drop=true' "$work/control.log"
host_failures=0
native_failures=0
for kb in $(seq 1024 32 8192); do
  log="$work/$kb.log"
  if ! timeout -k 1 10 "$work/probe" "$extension" "$kb" 10000 recover > "$log" 2>&1; then
    printf 'Annotation did not recover with %s KiB additional address space.\n' "$kb" >&2
    grep . "$log" >&2 || true
    exit 1
  fi
  grep -q '^recovery=ok' "$log"
  grep -q '^drop=true' "$log"
  if grep -q '^annotation=error Out of Memory Error: Allocation failure' "$log" &&
     ! grep -q '^LINE 1:' "$log"; then
    host_failures=$((host_failures + 1))
  fi
  if grep -q '^annotation=error Invalid Input Error: duckvep_annotate: out of memory' "$log"; then
    native_failures=$((native_failures + 1))
  fi
done
((host_failures > 0 && native_failures > 0)) || {
  echo "Required failure paths not exercised: host=$host_failures native=$native_failures" >&2
  exit 1
}
printf 'Host OOM recovery passed: %s host failures, %s native failures; pins released.\n' "$host_failures" "$native_failures"
