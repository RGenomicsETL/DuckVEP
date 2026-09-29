#!/usr/bin/env bash
# The #2 slice 7 scale qualification in one command (contract section 4 of design/duckvep_haplotype_contract.md).
#
#   WORK=/path/to/inputs bash benchmarks/haplotype_scale/run_qualification.sh
#
# Inputs (built once, untimed, by `Rscript benchmarks/haplotype_scale/prepare_inputs.R WORK`): hg002_calls*.parquet,
# qual5m.vcf.gz, qual5m_calls.parquet, dense/low_sharing VCFs and calls. Environment:
#   WORK      inputs directory (required)          OUT     raw results directory (default $WORK/results)
#   EXT       extension to time (default build/release/duckvep.duckdb_extension; copied to an immutable file first)
#   CORE      the single core (default 6)          LOAD_MAX  1-minute load ceiling before a timed run (default 4)
#   ROUNDS    fresh processes per mode (default 3) VARIANTS  "name=/path/duckvep.duckdb_extension ..." for the A/B of
#             the optimizations (mode A, cold, alternating rounds)
#   STEPS     sections to run, default "hg002 domain qual5m controls variants extras failures" (a run resumes with the rest)
# Every timed run is a fresh process pinned to CORE inside a systemd scope with MemoryMax=16G and no swap (capped_run.sh),
# with a 4 GiB native budget, DuckDB memory_limit 8GB, a temp directory capped at 32GiB, and the load average checked
# (and waited on) before it starts. Rounds alternate csq and DuckVEP so both see the same host conditions.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
: "${WORK:?set WORK to the prepared inputs directory}"
OUT=${OUT:-$WORK/results}; CORE=${CORE:-6}; ROUNDS=${ROUNDS:-3}; export DUCKVEP_LOAD_MAX=${LOAD_MAX:-4}
data=/root/duckvep/data
vcf=$data/hg002-csq/hg002.ens.vcf.gz
fasta=$data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gff=$data/hg002-csq/Homo_sapiens.GRCh38.116.gff3.gz
mkdir -p "$OUT/ext" "$OUT/spill"
if [[ ! -f $OUT/ext/duckvep.duckdb_extension ]]; then   # a resumed run keeps the binary it started with
    cp "${EXT:-$root/build/release/duckvep.duckdb_extension}" "$OUT/ext/duckvep.duckdb_extension"
    chmod a-w "$OUT/ext/duckvep.duckdb_extension"
fi
ext=$OUT/ext/duckvep.duckdb_extension
sha256sum "$ext" | cut -d' ' -f1 > "$OUT/ext/sha256"
proc=$OUT/process.tsv

# Identity of every input and tool, by content hash (paths are not identity). Written once: a resumed run keeps them.
[[ -s $OUT/inputs.tsv ]] || {
    printf 'artifact\tsha256\tbytes\n'
    for f in "$vcf" "$vcf.tbi" "$gff" "$fasta" "$data/models/homo_sapiens_116_GRCh38_final.duckdb" "$ext" \
        "$WORK/qual5m.vcf.gz" "$WORK/hg002_calls.parquet" "$WORK/qual5m_calls.parquet" "$WORK/dense.vcf.gz" \
        "$WORK/low_sharing.vcf.gz" "$WORK/hg002_domain.counts.tsv.records.tsv.gz"; do
        [[ -f $f ]] && printf '%s\t%s\t%s\n' "$(basename "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c %s "$f")"
    done
    printf 'bcftools\t%s\t-\n' "$(/usr/local/bin/bcftools --version-only)"
    printf 'duckdb_r\t%s\t-\n' "$(Rscript -e 'cat(as.character(packageVersion("duckdb")))' 2>/dev/null)"
    printf 'git_head\t%s\t-\n' "$(git -C "$root" rev-parse HEAD)"
    printf 'src_tree\t%s\t-\n' "$(git -C "$root" rev-parse HEAD:src)"
    printf 'cgranges_tree\t%s\t-\n' "$(git -C "$root" rev-parse HEAD:third_party/cgranges)"
} > "$OUT/inputs.tsv"
# Physical records and ALT alleles of every VCF input (untimed).
[[ -s $OUT/input_counts.tsv ]] || Rscript "$here/input_counts.R" "$vcf" "$WORK/qual5m.vcf.gz" "$WORK/dense.vcf.gz" "$WORK/low_sharing.vcf.gz" > "$OUT/input_counts.tsv" 2>/dev/null

timed() { # label command...
    local label=$1; shift
    SPILL_DIR="$OUT/spill/$label" "$here/capped_run.sh" "$CORE" "$proc" "$label" -- "$@"
    rm -rf "$OUT/spill/$label"
}
csq() { # label output-flags...
    local label=$1; shift
    timed "$label" /usr/local/bin/bcftools csq -f "$fasta" -g "$gff" -p a "$@" "$vcf"
}
# worker.R defaults: max_alignment_cells 268,435,456 and workspace_limit 1 GiB, raised from the builder's 16,777,216 and
# 256 MiB (the default refuses HG002's longest transcripts and the 5M job, which failures.R and the README record). The dense
# control needs more: thousands of edits on one 100 kb coding sequence.
dv() { # label mode model input-flag input out extra...
    local label=$1 mode=$2 model=$3 flag=$4 input=$5; shift 5
    timed "$label" Rscript "$here/worker.R" --extension "$ext" --mode "$mode" --model "$model" "$flag" "$input" \
        --out "$OUT/$label.parquet" --receipt "$OUT/$label.worker.tsv" --label "$label" --spill "$OUT/spill/$label" "$@"
    rm -f "$OUT/$label.parquet"
}

want() { [[ " ${STEPS:-hg002 domain qual5m controls variants extras failures} " == *" $1 "* ]]; }
if want hg002; then
for round in $(seq "$ROUNDS"); do
    # Gate and diagnostic csq baselines: uncompressed output discarded (the reconstruction of the recorded baseline) and
    # BCF written to a file (the oracle scripts' output). The DuckVEP gate runs are cold-only fresh processes, so the
    # process wall clock is the cold wall clock (no audit: reading the output back for its checksum is not part of the job);
    # the *W runs add a second, warm pass in the same process and carry the audit (counts and full-output checksum).
    csq "csq_Ou_$round" -Ou -o /dev/null
    dv "hg002_B_$round" B full --vcf "$vcf" --no-audit
    dv "hg002_Bbgzip_$round" B full --vcf "$vcf" --decoder bgzip --no-audit   # diagnostic: the VCF inflated by bgzip -dc (libdeflate) instead of DuckDB's miniz
    dv "hg002_A_$round" A full --calls "$WORK/hg002_calls.parquet" --no-audit
    csq "csq_Ob_$round" -Ob -o "$OUT/csq_Ob_$round.bcf"; rm -f "$OUT/csq_Ob_$round.bcf"
    dv "hg002_BW_$round" B full --vcf "$vcf" --warm
    dv "hg002_AW_$round" A full --calls "$WORK/hg002_calls.parquet" --warm
done
fi
if want domain; then
for round in $(seq "$ROUNDS"); do
    for part in inside outside; do
        dv "hg002_${part}_A_$round" A full --calls "$WORK/hg002_calls_$part.parquet"
    done
done
fi
if want qual5m; then
for round in $(seq "$ROUNDS"); do
    dv "qual5m_B_$round" B mane --vcf "$WORK/qual5m.vcf.gz" --no-audit
    dv "qual5m_A_$round" A mane --calls "$WORK/qual5m_calls.parquet" --no-audit
    dv "qual5m_BW_$round" B mane --vcf "$WORK/qual5m.vcf.gz" --warm
    dv "qual5m_AW_$round" A mane --calls "$WORK/qual5m_calls.parquet" --warm
done
fi
if want controls; then
for control in dense low_sharing; do
    for round in $(seq "$ROUNDS"); do
        caps=()
        [[ $control == dense ]] && caps=(--max-alignment-cells 1073741824 --workspace-limit 2147483648)
        dv "${control}_B_$round" B mane --vcf "$WORK/$control.vcf.gz" "${caps[@]}"
        dv "${control}_A_$round" A mane --calls "$WORK/${control}_calls.parquet" "${caps[@]}"
    done
done
fi
if want variants && [[ -n "${VARIANTS:-}" ]]; then
    for round in $(seq "$ROUNDS"); do
        for variant in $VARIANTS; do
            name=${variant%%=*}; path=${variant#*=}
            timed "ab_${name}_$round" Rscript "$here/worker.R" --extension "$path" --mode A --model full \
                --calls "$WORK/hg002_calls.parquet" --out "$OUT/ab_${name}_$round.parquet" \
                --receipt "$OUT/ab_${name}_$round.worker.tsv" --label "ab_${name}_$round" --spill "$OUT/spill/ab_${name}_$round"
            rm -f "$OUT/ab_${name}_$round.parquet"
        done
    done
fi
if want extras; then
    # A/B of two pipeline choices, alternating, cold: ORDER BY on the model queries (load time; the compiled model is stored in load
    # order) and transcript discovery by the annotation builder against duckvep_coding_transcripts (mode B stage time).
    for round in $(seq "$ROUNDS"); do
        for order in off on; do
            flag=(); [[ $order == on ]] && flag=(--order-load)
            dv "abx_order_${order}_$round" A full --calls "$WORK/hg002_calls.parquet" --calls-filter "event_index < 0" --no-audit "${flag[@]}"
        done
        for method in annotate scalar; do
            dv "abx_discovery_${method}_$round" B full --vcf "$vcf" --discovery "$method"
        done
    done
fi
want failures && Rscript "$here/failures.R" "$ext" "$WORK/hg002_calls.parquet" "$OUT/failures" "$OUT/failures.tsv"
echo "raw receipts in $OUT; summarize with: Rscript $here/summarize.R $OUT"
