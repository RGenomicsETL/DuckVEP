#!/usr/bin/env bash
# Publish DuckVEP's benchmark inputs to a Hugging Face dataset repository, as the VCFs the runners read and as
# trimmed Parquet with one file per chromosome.
#
#   HF_TOKEN=... bash scripts/publish-benchmark-inputs.sh [--public] [--dry-run] REPO_ID
#
# Inputs (environment, defaults under /root/duckvep/data):
#   HG002_VCF     hg002.ens.vcf.gz from benchmarks/prepare_hg002.sh
#   QUAL5M_VCF    qual5m.vcf.gz from benchmarks/haplotype_scale/prepare_inputs.R
#   CONTROLS_DIR  dense.vcf.gz and low_sharing.vcf.gz from the same script
#   CORPUS_DB     the throughput corpus (table bench_variants, region ordinals of MODEL_DB)
#   MODEL_DB      the GRCh38 model whose model_regions name the corpus ordinals
#   CORPUS_VCF    GIAB's HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz, the corpus source, used only to check it
# Every Parquet relation has one row per ALT allele and is checked against an independent parse: bcftools query
# of the VCF it comes from. A new repository is private unless --public is given. STAGE holds the staged files.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
visibility=--private; dry_run=0
while [[ $# -gt 0 && $1 == --* ]]; do
    case $1 in --public) visibility=--public;; --dry-run) dry_run=1;; *) echo "unknown option $1" >&2; exit 2;; esac
    shift
done
repo=${1:?usage: publish-benchmark-inputs.sh [--public] [--dry-run] REPO_ID}
data=/root/duckvep/data
hg002=${HG002_VCF:-$data/hg002-csq/hg002.ens.vcf.gz}
qual5m=${QUAL5M_VCF:-$data/haplotype-scale/qual5m.vcf.gz}
controls=${CONTROLS_DIR:-$data/haplotype-scale}
corpus=${CORPUS_DB:-$data/scale5/corpus.duckdb}
model=${MODEL_DB:-$data/models/homo_sapiens_116_GRCh38_final.duckdb}
corpus_vcf=${CORPUS_VCF:-/root/.cache/duckhts/benchmarks/variantkey-providers/raw/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}
duckdb_bin=${DUCKDB:-duckdb}; bcftools=${BCFTOOLS:-bcftools}
[[ $dry_run == 1 || -n ${HF_TOKEN:-} ]] || { echo "set HF_TOKEN" >&2; exit 1; }
stage=${STAGE:-$(mktemp -d)}; mkdir -p "$stage"; work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
sql_text() { printf "'%s'" "${1//\'/\'\'}"; }
columns="{'c0':'VARCHAR','c1':'VARCHAR','c2':'VARCHAR','c3':'VARCHAR','c4':'VARCHAR','c5':'VARCHAR','c6':'VARCHAR','c7':'VARCHAR','c8':'VARCHAR','c9':'VARCHAR'}"

# One row per ALT allele: record_index (0-based record ordinal in the VCF), chrom, position, reference, alternate,
# alt_index (1-based), filter, gt and ps (raw sample text), and origin when an expression for it is given.
alleles_sql() { # vcf [origin_expression]
    local skip origin=""
    skip=$(zcat "$1" | head -n 100000 | grep -c '^#')
    [[ -n ${2:-} ]] && origin=", $2 AS origin"
    cat <<SQL
SET threads = 1; SET preserve_insertion_order = true;
CREATE TEMP TABLE records AS SELECT row_number() OVER () - 1 AS record_index, c0 AS chrom, c1 AS pos, c3 AS ref,
  string_split(c4, ',') AS alts, c6 AS "filter", c8 AS format, string_split(c8, ':') AS keys, string_split(c9, ':') AS vals
FROM read_csv($(sql_text "$1"), delim = '\t', header = false, skip = $skip, quote = '', escape = '',
  auto_detect = false, columns = $columns);
CREATE TEMP TABLE alleles AS SELECT record_index::UBIGINT AS record_index, chrom, pos::UBIGINT AS position,
  ref AS reference, alts[i] AS alternate, i::UTINYINT AS alt_index, "filter",
  vals[list_position(keys, 'GT')] AS gt, vals[list_position(keys, 'PS')] AS ps$origin
FROM records, range(1, len(alts) + 1) t(i) ORDER BY record_index, alt_index;
SQL
}
# bcftools' parse of the same VCF, one row per ALT allele, as an order-independent fingerprint.
check_sql() { # vcf
    "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t[%PS]\n' "$1" |
        awk 'BEGIN { OFS = "\t" } { print NR - 1, $0 }' > "$work/check.tsv"
    cat <<SQL
CREATE TEMP TABLE expected AS SELECT c0::UBIGINT AS record_index, c1 AS chrom, c2::UBIGINT AS position, c3 AS reference,
  string_split(c4, ',')[i] AS alternate, i::UTINYINT AS alt_index, c5 AS "filter", c6 AS gt, c7 AS ps
FROM read_csv($(sql_text "$work/check.tsv"), delim = '\t', header = false, quote = '', escape = '', auto_detect = false,
  columns = {'c0':'VARCHAR','c1':'VARCHAR','c2':'VARCHAR','c3':'VARCHAR','c4':'VARCHAR','c5':'VARCHAR','c6':'VARCHAR','c7':'VARCHAR'}),
  range(1, len(string_split(c4, ',')) + 1) t(i);
SQL
}
fingerprint="struct_pack(n := count(*), h := sum(hash(record_index, chrom, position, reference, alternate, alt_index, \"filter\", gt, coalesce(ps, '.'))::HUGEINT))"
publish_vcf() { # vcf dir [origin_expression] [single]
    local vcf=$1 dir=$stage/$2 target
    mkdir -p "$dir"; cp "$vcf" "$vcf.tbi" "$dir/"
    if [[ ${4:-} == single ]]; then target="$(sql_text "$dir/$(basename "$vcf" .vcf.gz).parquet") (FORMAT parquet, COMPRESSION zstd)"
    else target="$(sql_text "$dir/alleles") (FORMAT parquet, PARTITION_BY (chrom), WRITE_PARTITION_COLUMNS true, COMPRESSION zstd, OVERWRITE true)"; fi
    "$duckdb_bin" :memory: -bail -c "$(alleles_sql "$vcf" "${3:-}") $(check_sql "$vcf")
COPY (FROM alleles) TO $target;
CREATE TEMP TABLE published AS FROM read_parquet($(sql_text "$dir/**/*.parquet"), hive_partitioning = false);
SELECT CASE WHEN (SELECT $fingerprint FROM published) = (SELECT $fingerprint FROM expected)
  AND (SELECT count(*) FROM published) = (SELECT count(*) FROM alleles) THEN 'ok'
  ELSE error('$2: Parquet differs from bcftools parse of the VCF') END;" > /dev/null
    echo "$2: $("$duckdb_bin" -noheader -list :memory: -c "SELECT count(DISTINCT record_index) || ' records, ' || count(*) || ' ALT alleles' FROM read_parquet($(sql_text "$dir/**/*.parquet"), hive_partitioning = false)")" >&2
}

publish_vcf "$hg002" hg002
publish_vcf "$qual5m" qual5m "CASE WHEN format = 'GT:PS' THEN 'gnomad_v4.1_exomes' ELSE 'giab_hg002' END"
publish_vcf "$controls/dense.vcf.gz" controls/dense "" single
publish_vcf "$controls/low_sharing.vcf.gz" controls/low_sharing "" single

# The throughput corpus, with the model's region ordinals replaced by chromosome names and checked against the
# GIAB benchmark VCF it was taken from (chr prefix dropped, record ordinals from 0, one row per ALT allele).
mkdir -p "$stage/throughput"
"$duckdb_bin" :memory: -bail -c "$(check_sql "$corpus_vcf")
ATTACH $(sql_text "$corpus") AS c (READ_ONLY); ATTACH $(sql_text "$model") AS m (READ_ONLY);
CREATE TEMP TABLE corpus AS SELECT (v.source_record_index - 1)::UBIGINT AS record_index, r.seq_region_name AS chrom,
  v.position::UBIGINT AS position, v.reference, v.alternate, v.source_alt_index::UTINYINT AS alt_index
FROM c.bench_variants v JOIN m.model_regions r USING (seq_region) ORDER BY record_index, alt_index;
COPY (FROM corpus) TO $(sql_text "$stage/throughput/alleles") (FORMAT parquet, PARTITION_BY (chrom), WRITE_PARTITION_COLUMNS true, COMPRESSION zstd, OVERWRITE true);
SELECT CASE WHEN (SELECT count(*) FROM c.bench_variants) = (SELECT count(*) FROM read_parquet($(sql_text "$stage/throughput/alleles/**/*.parquet"), hive_partitioning = false))
  AND NOT EXISTS (SELECT record_index, chrom, position, reference, alternate, alt_index FROM corpus
    EXCEPT SELECT record_index, replace(chrom, 'chr', ''), position, reference, alternate, alt_index FROM expected)
  THEN 'ok' ELSE error('throughput: corpus rows are missing from the GIAB VCF or lost in Parquet') END;" > /dev/null
echo "throughput: $("$duckdb_bin" -noheader -list :memory: -c "SELECT count(*) || ' ALT alleles' FROM read_parquet($(sql_text "$stage/throughput/alleles/**/*.parquet"), hive_partitioning = false)")" >&2

(cd "$stage" && find . -type f ! -name README.md ! -name SHA256SUMS -printf '%P\n' | sort | xargs sha256sum > SHA256SUMS)
cp "$root/scripts/benchmark-inputs-card.md" "$stage/README.md"
sed -i "s|REPO_ID|$repo|g" "$stage/README.md"
echo "staged in $stage" >&2
[[ $dry_run == 1 ]] && exit 0

hf repos create "$repo" --repo-type dataset "$visibility" --exist-ok > /dev/null
hf upload "$repo" "$stage" . --repo-type dataset --commit-message "Publish DuckVEP benchmark inputs"
