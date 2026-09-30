#!/usr/bin/env bash
set -euo pipefail

# VEP 116 with the LOFTEE plugin (konradjk/loftee at a46b502), GRCh38, offline, from
# the pinned image. The profile is LOFTEE's defaults with the human ancestor and
# PhyloCSF off; GERP is the constant-negative bigwig, which makes the GERP test
# always pass and leaves LOFTEE's unweighted 50 bp rule.
# Optional trailing arguments are appended to the LoF plugin parameters.
if [[ $# -lt 6 ]]; then
  echo "usage: $0 CACHE_ROOT FASTA LOFTEE_DIR GERP_BIGWIG INPUT_VCF OUTPUT_JSON [LOF_PARAM...]" >&2
  exit 2
fi
image='ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11'
cache=$(realpath "$1")
fasta=$(realpath "$2")
loftee=$(realpath "$3")
gerp=$(realpath "$4")
input=$(realpath "$5")
output_dir=$(realpath -m "$(dirname "$6")")
output_name=$(basename "$6")
params=""
for p in "${@:7}"; do params="$params,$p"; done
test -f "$cache/homo_sapiens/116_GRCh38/info.txt"
test -f "$fasta.fai"
test -f "$loftee/LoF.pm"
test -f "$input"
test -d "$output_dir"

docker run --rm --network none --user "$(id -u):$(id -g)" \
  --entrypoint vep \
  --mount "type=bind,src=$cache,dst=/cache,readonly" \
  --mount "type=bind,src=$(dirname "$fasta"),dst=/reference,readonly" \
  --mount "type=bind,src=$loftee,dst=/loftee,readonly" \
  --mount "type=bind,src=$(dirname "$gerp"),dst=/gerp,readonly" \
  --mount "type=bind,src=$(dirname "$input"),dst=/input,readonly" \
  --mount "type=bind,src=$output_dir,dst=/output" \
  "$image" \
  -i "/input/$(basename "$input")" \
  --fasta "/reference/$(basename "$fasta")" \
  --cache --offline --dir_cache /cache \
  --species homo_sapiens --assembly GRCh38 --cache_version 116 \
  --distance 5000 --buffer_size 5000 \
  --dir_plugins /loftee \
  --plugin "LoF,loftee_path:/loftee,gerp_bigwig:/gerp/$(basename "$gerp")$params" \
  --json --no_stats --force_overwrite \
  -o "/output/$output_name"
