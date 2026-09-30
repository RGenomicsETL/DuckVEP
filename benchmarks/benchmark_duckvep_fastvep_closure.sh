#!/bin/bash
# Paired, alternating native_tab17 runs of DuckVEP (standalone duckvep
# extension) and FastVEP over the full GIAB HG002 v4.2.1 VCF, with the full-row
# receipts of benchmark_duckvep_fastvep_receipt.R. Usage:
#   THREADS=8 CPU=2,4,8,10,12,14,16,18 benchmark_duckvep_fastvep_closure.sh OUT_DIR [REPETITIONS]
# THREADS (default 1) is DuckVEP --threads and FastVEP RAYON_NUM_THREADS; CPU
# (default 2) is the taskset list.
# Every input is passed by environment variable; nothing is rebuilt here.
set -u
OUT=${1:?output directory}; REPS=${2:-3}
: "${DUCKVEP_EXT:?}" "${DUCKHTS_EXT:?}" "${FASTVEP_BIN:?}" "${FASTVEP_CACHE:?}" "${MODEL:?}" \
  "${INPUT_VCF:?}" "${GFF3:?}" "${SOURCE_MAP:?}" "${WORK:?}" "${TMPDIR:?}"
CPU=${CPU:-2}; THREADS=${THREADS:-1}
ROOT=$(git rev-parse --show-toplevel)
mkdir -p "$OUT" "$WORK" "$TMPDIR"
SRC_SHA=$(sha256sum "$INPUT_VCF" | cut -d' ' -f1)
waitload() { while :; do l=$(cut -d' ' -f1 /proc/loadavg); awk -v l="$l" 'BEGIN{exit !(l<=3)}' && break; sleep 60; done
  echo "$1 $(cut -d' ' -f1-3 /proc/loadavg)" >> "$OUT/load.log"; }
run_one() { # engine run
  local engine=$1 run=$2 label=${1}_native_tab17_${THREADS}_${2} tab="$WORK/output.tab"
  rm -f "$tab"; waitload "$label"
  if [ "$engine" = duckvep ]; then
    /usr/bin/time -v -o "$OUT/$label.time" taskset -c "$CPU" Rscript "$ROOT/benchmarks/benchmark_duckvep_fastvep_worker.R" \
      --extension "$DUCKVEP_EXT" --duckhts-extension "$DUCKHTS_EXT" --model "$MODEL" --input "$INPUT_VCF" \
      --output "$tab" --output-contract native_tab17 --gff3 "$GFF3" --threads "$THREADS" --distance 5000 \
      --memory-limit 16GB --max-spill 8GB > "$OUT/$label.log" 2>&1
  else
    /usr/bin/time -v -o "$OUT/$label.time" taskset -c "$CPU" env RAYON_NUM_THREADS="$THREADS" "$FASTVEP_BIN" annotate \
      --input "$INPUT_VCF" --output "$tab" --output-format tab --transcript-cache "$FASTVEP_CACHE" \
      --distance 5000 --no-progress > "$OUT/$label.log" 2>&1
  fi
  local status=$?; echo "$label exit=$status" >> "$OUT/load.log"; [ $status -eq 0 ] || return 1
  local skip; skip=$(Rscript -e 'source(file.path(commandArgs(TRUE)[1], "benchmarks/benchmark_duckvep_fastvep_fields.R")); cat(duckvep_fastvep_tab_header(commandArgs(TRUE)[2], duckvep_fastvep_transport_fields("native_tab17")))' "$ROOT" "$tab" 2>/dev/null | tail -1)
  Rscript "$ROOT/benchmarks/benchmark_duckvep_fastvep_receipt.R" --input "$tab" --tool "$engine" \
    --output-contract native_tab17 --threads "$THREADS" --run "$run" --timing-file "$OUT/$label.time" \
    --output "$OUT/$label.csv" --skip-lines "$skip" --source-map "$SOURCE_MAP" --source-sha256 "$SRC_SHA" \
    --coverage-output "$OUT/$label.coverage.csv" --memory-limit 16GB --max-spill 8GB >> "$OUT/$label.log" 2>&1
  echo "$label receipt=$?" >> "$OUT/load.log"
  rm -f "$tab"
}
for run in $(seq 1 "$REPS"); do
  if [ $((run % 2)) -eq 1 ]; then order="duckvep fastvep"; else order="fastvep duckvep"; fi
  for engine in $order; do run_one "$engine" "$run"; done
done
echo DONE >> "$OUT/load.log"
