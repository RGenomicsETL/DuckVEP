#!/usr/bin/env bash
# Publish DuckVEP model databases to a Hugging Face dataset repository.
#
#   HF_TOKEN=... bash scripts/publish-models.sh [--public] [--dry-run] REPO_ID MODEL.duckdb...
#
# Each model is exported to a database holding only what duckvep_model_load reads: model_regions,
# model_transcripts, duckvep_regulation_features when present, and model_receipt. The export is accepted
# only when a native snapshot saved from it is byte-identical to one saved from the source model.
# The upload holds the exports under ensembl-<release>/, SHA256SUMS and a dataset card (README.md).
# A new repository is private unless --public is given; an existing one keeps its visibility.
# --dry-run stops after the exports, checks and card. STAGE (default a temporary directory) holds them.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
visibility=--private; dry_run=0
while [[ $# -gt 0 && $1 == --* ]]; do
    case $1 in --public) visibility=--public;; --dry-run) dry_run=1;; *) echo "unknown option $1" >&2; exit 2;; esac
    shift
done
[[ $# -ge 2 ]] || { echo "usage: publish-models.sh [--public] [--dry-run] REPO_ID MODEL.duckdb..." >&2; exit 2; }
repo=$1; shift
duckdb_bin=${DUCKDB:-duckdb}
ext=${DUCKVEP_EXTENSION:-$root/build/release/duckvep.duckdb_extension}
[[ -f $ext ]] || { echo "build the extension first: $ext" >&2; exit 1; }
[[ $dry_run == 1 || -n ${HF_TOKEN:-} ]] || { echo "set HF_TOKEN" >&2; exit 1; }
stage=${STAGE:-$(mktemp -d)}; mkdir -p "$stage"
sql_text() { printf "'%s'" "${1//\'/\'\'}"; }

# The load duckvep_model_load performs from a model database attached as src.
load_sql() { # name with_regulation
    local regulation=""
    [[ $2 == 1 ]] && regulation=", interval_feature_query := 'SELECT regulation_feature_index, seq_region, feature_start, feature_end, feature_kind FROM src.duckvep_regulation_features ORDER BY seq_region, feature_start, regulation_feature_index'"
    cat <<SQL
SELECT loaded FROM duckvep_model_load('$1',
 'SELECT seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name FROM src.model_regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM src.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM src.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM src.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid, x.edit_code FROM src.model_transcripts, unnest(peptide_edits) u(x) ORDER BY 1, 2',
 transcript_coverage_complete := TRUE$regulation);
SQL
}
snapshot() { # model.duckdb with_regulation out.dvsnap
    "$duckdb_bin" -unsigned :memory: -bail -c "LOAD $(sql_text "$ext"); ATTACH $(sql_text "$1") AS src (READ_ONLY);
$(load_sql m "$2") SELECT duckvep_model_save('m', $(sql_text "$3"));" > /dev/null
}

rows=()
for source in "$@"; do
    has_regulation=$("$duckdb_bin" -readonly -noheader -list "$source" -c \
        "SELECT count(*) FROM duckdb_tables() WHERE table_name = 'duckvep_regulation_features'")
    IFS=$'\t' read -r release assembly < <("$duckdb_bin" -readonly -noheader -separator $'\t' -list "$source" -c \
        "SELECT source_version, assembly FROM model_receipt")
    name=$(basename "$source" .duckdb); name=${name%_final}
    dir=ensembl-$release; mkdir -p "$stage/$dir"
    export_path=$stage/$dir/$name.duckdb; rm -f "$export_path" "$export_path.wal"
    tables="model_regions model_transcripts model_receipt"; [[ $has_regulation == 1 ]] && tables+=" duckvep_regulation_features"
    copy="SET preserve_insertion_order = true; ATTACH $(sql_text "$source") AS src (READ_ONLY); ATTACH $(sql_text "$export_path") AS out;"
    for table in $tables; do copy+=" CREATE TABLE out.$table AS SELECT * FROM src.$table;"; done
    "$duckdb_bin" :memory: -bail -c "$copy DETACH out;" > /dev/null
    snapshot "$source" "$has_regulation" "$stage/source.dvsnap"
    snapshot "$export_path" "$has_regulation" "$stage/export.dvsnap"
    cmp -s "$stage/source.dvsnap" "$stage/export.dvsnap" || { echo "$name: exported model loads differently" >&2; exit 1; }
    rm -f "$stage/source.dvsnap" "$stage/export.dvsnap"
    IFS=$'\t' read -r transcripts regions regulation < <("$duckdb_bin" -readonly -noheader -separator $'\t' -list "$export_path" -c \
        "SELECT format('{:,}', transcript_count), format('{:,}', region_count), format('{:,}', coalesce(regulation_feature_count, 0)) FROM model_receipt")
    rows+=("| [\`$dir/$name.duckdb\`]($dir/$name.duckdb) | $assembly | $release | $transcripts | $regions | $regulation | $(du -h "$export_path" | cut -f1) |")
    echo "$name: exported and checked" >&2
done
(cd "$stage" && find . -name '*.duckdb' -printf '%P\n' | sort | xargs sha256sum > SHA256SUMS)
example=$(grep -m1 GRCh38 "$stage/SHA256SUMS" | awk '{print $2}'); example=${example:-$(awk 'NR == 1 {print $2}' "$stage/SHA256SUMS")}

{
cat <<EOF
---
license: other
license_name: ensembl
license_link: https://www.ensembl.org/info/about/legal/disclaimer.html
pretty_name: DuckVEP transcript models
tags:
- genomics
- variant-effect-prediction
- ensembl
- duckdb
---

# DuckVEP transcript models

Transcript models for [DuckVEP](https://github.com/RGenomicsETL/DuckVEP), a DuckDB extension that annotates variant consequences the way Ensembl VEP does. Each file is a DuckDB database built from an Ensembl core database and the matching primary-assembly FASTA. It holds the regions, transcripts with their exons and coding and flanking sequence, and, where present, the regulatory and motif features that \`duckvep_model_load\` reads. The \`model_receipt\` table records the source release, the transcript filter and the source and reference hashes.

| File | Assembly | Ensembl release | Transcripts | Regions | Regulatory features | Size |
|---|---|---|---:|---:|---:|---:|
EOF
printf '%s\n' "${rows[@]}"
cat <<EOF

\`SHA256SUMS\` lists the checksum of every file.

## Use

Download a model, then load it once per DuckDB session:

\`\`\`sh
hf download $repo $example --repo-type dataset --local-dir models
\`\`\`

\`\`\`sql
LOAD 'duckvep.duckdb_extension';
ATTACH 'models/$example' AS src (READ_ONLY);
EOF
load_sql grch38 "$(grep -q GRCh38 <<<"$example" && echo 1 || echo 0)"
cat <<EOF
\`\`\`

\`duckvep_model_save\` writes the loaded model to a snapshot that \`duckvep_model_restore\` maps back in under a second; snapshots belong to the extension build that wrote them, so they are not published here. The [function reference](https://github.com/RGenomicsETL/DuckVEP/blob/main/docs/functions.md) describes loading and annotation.

## Terms

The data come from [Ensembl](https://www.ensembl.org), which makes them available without restriction; see the [Ensembl disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html). When you use them, cite Ensembl and DuckVEP.
EOF
} > "$stage/README.md"
echo "staged in $stage" >&2
[[ $dry_run == 1 ]] && exit 0

hf repos create "$repo" --repo-type dataset "$visibility" --exist-ok > /dev/null
hf upload "$repo" "$stage" . --repo-type dataset --exclude '*.wal' --commit-message "Publish DuckVEP models"
