#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
duckdb_bin=${DUCKDB:-duckdb}
bcftools_bin=${BCFTOOLS:-bcftools}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
mkdir -p "$work/controls"
cat > "$work/input.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10	.	A	C,G	.	PASS	.	GT:PS	0|1:7
chr1	20	.	T	G	.	PASS	.	GT:PS	1|0:7
VCF
bgzip -c "$work/input.vcf" > "$work/hg002.vcf.gz"
tabix -p vcf "$work/hg002.vcf.gz"
for target in "$work/qual5m.vcf.gz" "$work/controls/dense.vcf.gz" "$work/controls/low_sharing.vcf.gz" "$work/corpus.vcf.gz"; do
  cp "$work/hg002.vcf.gz" "$target"
  cp "$work/hg002.vcf.gz.tbi" "$target.tbi"
done
"$duckdb_bin" "$work/model.duckdb" -c "CREATE TABLE model_regions AS SELECT 0::UINTEGER seq_region, '1' seq_region_name" >/dev/null
"$duckdb_bin" "$work/corpus.duckdb" -c "CREATE TABLE bench_variants AS SELECT * FROM (VALUES (1,0,10,'A','C',1),(1,0,10,'A','G',2),(2,0,20,'T','G',1)) AS t(source_record_index,seq_region,position,reference,alternate,source_alt_index)" >/dev/null
export DUCKDB="$duckdb_bin" BCFTOOLS="$bcftools_bin"
export HG002_VCF="$work/hg002.vcf.gz" QUAL5M_VCF="$work/qual5m.vcf.gz" CONTROLS_DIR="$work/controls"
export CORPUS_DB="$work/corpus.duckdb" MODEL_DB="$work/model.duckdb" CORPUS_VCF="$work/corpus.vcf.gz"
stage=$work/stage
STAGE="$stage" bash "$root/scripts/publish-benchmark-inputs.sh" --dry-run example/benchmarks
(cd "$stage" && sha256sum -c SHA256SUMS >/dev/null)
[[ $("$duckdb_bin" -noheader -list :memory: -c "SELECT count(*) FROM read_parquet('$stage/hg002/alleles/**/*.parquet')") == 3 ]]
grep -F "| \`hg002/\` | 2 | 3 |" "$stage/README.md" >/dev/null
grep -F "| \`throughput/\` | 2 | 3 |" "$stage/README.md" >/dev/null
if grep -Eq '@[A-Z_]+@' "$stage/README.md"; then echo 'dataset card has unresolved count fields' >&2; exit 1; fi

stage=$work/occupied
mkdir -p "$stage"
printf 'keep\n' > "$stage/source"
if STAGE="$stage" bash "$root/scripts/publish-benchmark-inputs.sh" --dry-run example/benchmarks >"$work/occupied.log" 2>&1; then
  echo 'nonempty stage was accepted' >&2; exit 1
fi
grep -F 'STAGE must be an empty directory' "$work/occupied.log" >/dev/null
[[ $(<"$stage/source") == keep ]]

"$duckdb_bin" "$work/corpus.duckdb" -c "UPDATE bench_variants SET alternate = 'T' WHERE source_record_index = 2" >/dev/null
stage=$work/invalid
if STAGE="$stage" bash "$root/scripts/publish-benchmark-inputs.sh" --dry-run example/benchmarks >"$work/invalid.log" 2>&1; then
  echo 'throughput allele mismatch was accepted' >&2; exit 1
fi
grep -F 'throughput: corpus rows are missing from the GIAB VCF or differ from Parquet' "$work/invalid.log" >/dev/null
echo 'benchmark input staging tests passed'
