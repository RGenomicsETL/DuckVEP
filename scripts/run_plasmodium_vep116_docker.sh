#!/usr/bin/env bash
set -euo pipefail

# Arguments: extracted cache root, reference FASTA, input VCF, output JSON.
if [[ $# -ne 4 ]]; then
  echo "usage: $0 CACHE_ROOT FASTA INPUT_VCF OUTPUT_JSON" >&2
  exit 2
fi
image='ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11'
cache=$(realpath "$1")
fasta=$(realpath "$2")
input=$(realpath "$3")
output_dir=$(realpath -m "$(dirname "$4")")
output_name=$(basename "$4")
test -f "$cache/plasmodium_falciparum/63_GCA000002765v3/info.txt"
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
  --species plasmodium_falciparum --assembly GCA000002765v3 \
  --cache_version 63 --distance 5000 --buffer_size 5000 \
  --json --no_stats --force_overwrite \
  -o "/output/$output_name"
