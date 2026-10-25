#!/usr/bin/env bash
# Rebuild hg002.ens.vcf.gz, the phased GIAB HG002 input of the haplotype, vep-rs and scale benchmarks.
#
#   bash benchmarks/prepare_hg002.sh OUT_DIR
#
# 1. Download GIAB's HG002 v4.2.1 GRCh38 hifiasm phase-transfer VCF and check its SHA-256.
# 2. Drop records whose sample column has a different number of fields than FORMAT declares
#    (25,254 records in the chr6 MHC region, which bcftools cannot read).
# 3. Rename chr1..chr22, chrX and chrY to Ensembl's 1..22, X and Y, and index the result.
# The output has 4,023,088 records and 4,070,522 ALT alleles. Needs curl, bgzip and bcftools.
set -euo pipefail
out=${1:?usage: prepare_hg002.sh OUT_DIR}
mkdir -p "$out"; cd "$out"
url=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf.gz
sha=d027d7f0fed5205af713d9d477d6d3281d6e18ea3adef67dfb620c0c9c222bf6
[[ -s hg002.phased.vcf.gz ]] || curl -sfL -o hg002.phased.vcf.gz "$url"
echo "$sha  hg002.phased.vcf.gz" | sha256sum -c --quiet
zcat hg002.phased.vcf.gz | awk -F'\t' 'BEGIN { OFS = "\t" } /^#/ { print; next }
    { if (split($9, keys, ":") == split($10, values, ":")) print; else dropped++ }
    END { print dropped " records with mismatched FORMAT and sample fields dropped" > "/dev/stderr" }' |
    bgzip > hg002.fixed.vcf.gz
for i in $(seq 1 22) X Y; do echo "chr$i $i"; done > chrmap.txt
bcftools annotate --rename-chrs chrmap.txt -Oz -o hg002.ens.vcf.gz hg002.fixed.vcf.gz
bcftools index -t hg002.ens.vcf.gz
rm hg002.fixed.vcf.gz chrmap.txt
echo "records: $(bcftools view -H hg002.ens.vcf.gz | wc -l)"
