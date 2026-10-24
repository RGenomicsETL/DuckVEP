#!/usr/bin/env bash
# HG002 whole-haplotype prediction against pinned bcftools csq, one core, fresh processes.
#
#   OUT=benchmarks/data/haplotype_scale/hg002_<rev> bash benchmarks/haplotype_scale/run_hg002_csq.sh
#
# ROUNDS rounds (default 3) alternate `bcftools csq -p a -Ou -o /dev/null`, the DuckVEP CLI pipeline
# (cli_worker.sh --mode F: model load, native VCF reader, duckvep_haplotypes to Parquet) and the same pipeline
# from R (worker.R --mode F). Every run is pinned to CORE (default 6) under capped_run.sh's caps and load wait.
# Writes process.tsv (one row per run), identity.txt (input, tool and output identities) and status_counts.tsv
# (paths per prediction status and reason).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../.." && pwd)
data=${DATA:-/root/duckvep/data}
vcf=$data/hg002-csq/hg002.ens.vcf.gz
gff=$data/hg002-csq/Homo_sapiens.GRCh38.116.gff3.gz
fasta=$data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa
model=$data/models/homo_sapiens_116_GRCh38_final.duckdb
: "${OUT:?set OUT to the result directory}"
CORE=${CORE:-6}; ROUNDS=${ROUNDS:-3}; export DUCKVEP_LOAD_MAX=${LOAD_MAX:-4}
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
mkdir -p "$work/ext"; cp "${EXT:-$root/build/release/duckvep.duckdb_extension}" "$work/ext/duckvep.duckdb_extension"
ext=$work/ext/duckvep.duckdb_extension
proc=$OUT/process.tsv; rm -f "$proc"
{
    printf 'git_head %s\n' "$(git -C "$root" rev-parse HEAD)"
    printf 'extension_sha256 %s\n' "$(sha256sum "$ext" | cut -d' ' -f1)"
    printf 'vcf_sha256 %s\n' "$(sha256sum "$vcf" | cut -d' ' -f1)"
    printf 'gff3_sha256 %s\n' "$(sha256sum "$gff" | cut -d' ' -f1)"
    printf 'model_sha256 %s\n' "$(sha256sum "$model" | cut -d' ' -f1)"
    printf 'bcftools %s\n' "$(/usr/local/bin/bcftools --version-only)"
    printf 'duckdb_cli %s\n' "$(duckdb -version | head -1)"
    printf 'duckdb_r %s\n' "$(Rscript -e 'cat(as.character(packageVersion("duckdb")))' 2>/dev/null)"
    printf 'cpu %s\n' "$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
    printf 'core %s\n' "$CORE"
} > "$OUT/identity.txt"
for round in $(seq "$ROUNDS"); do
    "$here/capped_run.sh" "$CORE" "$proc" "csq_$round" -- \
        /usr/local/bin/bcftools csq -f "$fasta" -g "$gff" -p a -Ou -o /dev/null "$vcf"
    "$here/capped_run.sh" "$CORE" "$proc" "duckvep_F_$round" -- \
        bash "$here/cli_worker.sh" --extension "$ext" --mode F --vcf "$vcf" --out "$work/F_$round.parquet" --model-file "$model"
    "$here/capped_run.sh" "$CORE" "$proc" "duckvep_R_F_$round" -- \
        Rscript "$here/worker.R" --extension "$ext" --mode F --model full --vcf "$vcf" \
        --out "$work/R_$round.parquet" --receipt "$work/R_$round.tsv" --label "duckvep_R_F_$round" --no-audit
done
rows=$(duckdb -noheader -list -c "SELECT count(*) FROM read_parquet('$work/F_1.parquet')")
sum=$(duckdb -noheader -list -c "SELECT sum(hash(t)::HUGEINT) FROM read_parquet('$work/F_1.parquet') t")
printf 'output rows %s checksum %s\n' "$rows" "$sum" >> "$OUT/identity.txt"
duckdb -c "COPY (SELECT prediction_status, prediction_reason, count(*) AS paths FROM read_parquet('$work/F_1.parquet')
  GROUP BY ALL ORDER BY paths DESC) TO '$OUT/status_counts.tsv' (HEADER, DELIMITER '\t')"
rm -f "$OUT"/*.err
