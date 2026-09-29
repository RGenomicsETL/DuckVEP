#!/usr/bin/env bash
set -euo pipefail

# Cache version is the Ensembl cache version, which can differ from VEP 116.
# Optional trailing arguments are passed to vep unchanged (for example --hgvs).
if [[ $# -lt 7 ]]; then
  echo "usage: $0 SPECIES ASSEMBLY CACHE_VERSION CACHE_ROOT FASTA INPUT_VCF OUTPUT_JSON [VEP_ARG...]" >&2
  exit 2
fi
image='ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11'
species=$1
assembly=$2
cache_version=$3
cache=$(realpath "$4")
fasta=$(realpath "$5")
input=$(realpath "$6")
output_dir=$(realpath -m "$(dirname "$7")")
output_name=$(basename "$7")
extra=("${@:8}")
test -f "$cache/$species/${cache_version}_${assembly}/info.txt"
test -f "$fasta.fai"
test -f "$input"
test -d "$output_dir"

docker run --rm --network none --user "$(id -u):$(id -g)" \
  --entrypoint vep \
  --mount "type=bind,src=$cache,dst=/cache,readonly" \
  --mount "type=bind,src=$(dirname "$fasta"),dst=/reference,readonly" \
  --mount "type=bind,src=$(dirname "$input"),dst=/input,readonly" \
  --mount "type=bind,src=$output_dir,dst=/output" \
  "$image" \
  -i "/input/$(basename "$input")" \
  --fasta "/reference/$(basename "$fasta")" \
  --cache --offline --dir_cache /cache \
  --species "$species" --assembly "$assembly" \
  --cache_version "$cache_version" --distance 5000 --buffer_size 5000 \
  --json --no_stats --force_overwrite \
  -o "/output/$output_name" ${extra[@]+"${extra[@]}"}
