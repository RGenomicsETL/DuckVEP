#!/usr/bin/env bash
set -euo pipefail

# Physical VCF bytes only transit pipes; one verified object is committed at a time.
# Usage: scripts/stage_gnomad_v4.1.sh sv | genomes chrY | exomes chrY
lane=${1:?lane required}
chrom=${2:-}
case "$lane:$chrom" in
  sv:)
    object=release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
    shard=sv-sites ;;
  genomes:chr*|exomes:chr*)
    case "${chrom#chr}" in
      [1-9]|1[0-9]|2[0-2]|X|Y) ;;
      *) echo 'invalid chromosome' >&2; exit 2 ;;
    esac
    object="release/4.1/vcf/$lane/gnomad.$lane.v4.1.sites.$chrom.vcf.bgz"
    shard="$lane-$chrom" ;;
  *) echo 'usage: stage_gnomad_v4.1.sh sv | genomes chrN | exomes chrN' >&2; exit 2 ;;
esac
repo=$(git rev-parse --show-toplevel)
manifest="$repo/benchmarks/data/gnomad-v4.1-manifest.tsv"
root=/root/duckvep/data/gnomad-v4.1
mkdir -p "$root"
exec 9> "$root/.stage.lock"
if ! flock -n 9; then echo 'another source object is streaming' >&2; exit 1; fi
read -r path bytes generation md5 crc < <(awk -F '\t' -v key="$object" '$1 == key {print $1, $2, $3, $4, $5}' "$manifest")
if [[ "${path:-}" != "$object" || ! "$bytes" =~ ^[0-9]+$ || ! "$generation" =~ ^[0-9]+$ ]]; then
  echo 'object missing from frozen manifest' >&2; exit 2
fi
if [[ $(awk -F '\t' -v key="$object" '$1 == key {n++} END {print n+0}' "$manifest") != 1 ]]; then
  echo 'duplicate manifest object' >&2; exit 2
fi
if [[ -d "$root/$shard" ]]; then
  echo "complete: $root/$shard"; exit 0
fi
staged=$(du -sb "$root" | cut -f1)
free=$(df -PB1 "$root" | awk 'NR==2 {print $4}')
if (( staged >= 25000000000 || free < 42000000000 )); then
  echo 'staging cap or shared disk headroom reached' >&2; exit 1
fi
partial="$root/$shard.partial"
rm -rf -- "$partial"
mkdir "$partial"
cp "$repo/scripts/stage_gnomad_v4.1.R" "$partial/stager.R"
export DUCKVEP_STAGED_BYTES="$staged"
url="https://storage.googleapis.com/gcp-public-data--gnomad/$object?generation=$generation"
echo "streaming $object ($bytes bytes)" >&2
if ! curl --fail --location --silent --show-error "$url" \
    | tee >(openssl dgst -md5 -binary | base64 -w0 > "$partial/stream.md5") \
          >(wc -c > "$partial/stream.bytes") \
    | bcftools view -H - \
    | Rscript "$partial/stager.R" "$lane" "$object" "$md5" "$partial"; then
  echo "incomplete: $partial (restart this chromosome to retry)" >&2
  exit 1
fi
# tee waits for its consumers; compare the compressed byte stream to the pinned object.
read -r actual_bytes < "$partial/stream.bytes"
actual_md5=$(< "$partial/stream.md5")
if [[ "$actual_bytes" != "$bytes" || "$actual_md5" != "$md5" ]]; then
  echo "source object checksum/length mismatch: $partial" >&2; exit 1
fi
printf 'object\tgeneration\tbytes_read\tmd5_base64\tcrc32c_base64\n%s\t%s\t%s\t%s\t%s\n' \
  "$object" "$generation" "$actual_bytes" "$actual_md5" "$crc" > "$partial/source.tsv"
rm -f "$partial/stream.md5" "$partial/stream.bytes" "$partial/stager.R"
rmdir "$partial/spill" 2>/dev/null || true
staged=$(du -sb "$root" | cut -f1)
if (( staged > 25000000000 )); then
  echo 'staging cap exceeded' >&2; exit 1
fi
mv "$partial" "$root/$shard"
echo "verified: $root/$shard" >&2
