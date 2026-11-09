#!/usr/bin/env bash
# Publish DuckVEP model databases to a Hugging Face dataset repository.
#
#   HF_TOKEN=... bash scripts/publish-models.sh [--public] [--dry-run] [--complete-coverage] [--mane-grch37 DIR] REPO_ID MODEL.duckdb...
#
# Each export contains model_regions, model_transcripts, optional regulatory features and RNA edits,
# and a recomputed model_receipt. Native snapshots compare the source and export under the load profile below.
# --complete-coverage asserts that the input contains the complete transcript set for intergenic annotation.
# --mane-grch37 DIR adds the downloaded mane-grch37 data release (verified against its SHA256SUMS) beside
# the GRCh37 model, after checking that the exported model's model_sha256 is the release's join key.
# The upload holds the exports under ensembl-<release>/, SHA256SUMS and a dataset card (README.md).
# A new repository is private unless --public is given; an existing one keeps its visibility.
# --dry-run stops after the exports, checks and card. STAGE must be empty and defaults to a temporary directory.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
visibility=--private; dry_run=0; complete_coverage=0; mane=
while [[ $# -gt 0 && $1 == --* ]]; do
    case $1 in
        --public) visibility=--public;; --dry-run) dry_run=1;;
        --complete-coverage) complete_coverage=1;;
        --mane-grch37) [[ $# -ge 2 ]] || { echo "--mane-grch37 needs a directory" >&2; exit 2; }; mane=$(cd "$2" && pwd); shift;;
        *) echo "unknown option $1" >&2; exit 2;;
    esac
    shift
done
[[ $# -ge 2 ]] || { echo "usage: publish-models.sh [--public] [--dry-run] [--complete-coverage] [--mane-grch37 DIR] REPO_ID MODEL.duckdb..." >&2; exit 2; }
repo=$1; shift
duckdb_bin=${DUCKDB:-duckdb}
ext=${DUCKVEP_EXTENSION:-$root/build/release/duckvep.duckdb_extension}
[[ -f $ext ]] || { echo "build the extension first: $ext" >&2; exit 1; }
[[ $dry_run == 1 || -n ${HF_TOKEN:-} ]] || { echo "set HF_TOKEN" >&2; exit 1; }
stage=${STAGE:-$(mktemp -d)}; mkdir -p "$stage"
[[ -z $(find "$stage" -mindepth 1 -maxdepth 1 -print -quit) ]] || { echo "STAGE must be an empty directory: $stage" >&2; exit 1; }
stage=$(cd "$stage" && pwd -P)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
sql_text() { printf "'%s'" "${1//\'/\'\'}"; }

# The load duckvep_model_load performs from a model database attached as src.
load_sql() { # name regulation circular cdna rna_authority rna_edits
    local regulation="" rna="" regions="seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name" transcripts="transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence" coverage=FALSE
    [[ $2 == 1 ]] && regulation=", interval_feature_query := 'SELECT regulation_feature_index::UINTEGER AS regulation_feature_index, seq_region::UINTEGER AS seq_region, feature_start::UINTEGER AS feature_start, feature_end::UINTEGER AS feature_end, feature_kind::UTINYINT AS feature_kind FROM src.duckvep_regulation_features ORDER BY seq_region, feature_start, regulation_feature_index'"
    [[ $3 == 1 ]] && regions+=", circular"
    [[ $4 == 1 ]] && transcripts+=", cdna_sequence"
    [[ $5 == 1 ]] && transcripts+=", rna_cds_cdna_start, rna_cds_cdna_end, rna_codon_table"
    [[ $6 == 1 ]] && rna=", rna_edit_query := 'SELECT transcript_index, cdna_start, cdna_end, expected_substrate, alternate_sequence, operation_kind, programme_name, evidence_name, boundary_overrides FROM src.model_rna_edits ORDER BY transcript_index, cdna_start DESC'"
    [[ $complete_coverage == 1 ]] && coverage=TRUE
    cat <<SQL
SELECT loaded FROM duckvep_model_load('$1',
 'SELECT $regions FROM src.model_regions ORDER BY seq_region',
 'SELECT $transcripts FROM src.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM src.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM src.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid, x.edit_code FROM src.model_transcripts, unnest(peptide_edits) u(x) ORDER BY 1, 2',
 transcript_coverage_complete := $coverage$regulation$rna);
SQL
}
snapshot() { # source regulation circular cdna rna_authority rna_edits output
    "$duckdb_bin" -unsigned :memory: -bail -c "LOAD $(sql_text "$ext"); ATTACH $(sql_text "$1") AS src (READ_ONLY);
$(load_sql m "$2" "$3" "$4" "$5" "$6") SELECT duckvep_model_save('m', $(sql_text "$7"));" > /dev/null
}

rows=(); example=; example_assembly=; example_regulation=0; example_circular=0; example_cdna=0; example_authority=0; example_rna=0
for source in "$@"; do
    [[ -f $source ]] || { echo "model file not found: $source" >&2; exit 1; }
    has_regulation=$("$duckdb_bin" -readonly -noheader -list "$source" -c \
        "SELECT count(*) FROM duckdb_tables() WHERE schema_name = 'main' AND table_name = 'duckvep_regulation_features'")
    metadata=$("$duckdb_bin" -readonly -noheader -list "$source" -c "SELECT
      count(*) FILTER (WHERE table_name = 'model_regions' AND column_name = 'circular'),
      count(*) FILTER (WHERE table_name = 'model_transcripts' AND column_name = 'cdna_sequence'),
      count(*) FILTER (WHERE table_name = 'model_transcripts' AND column_name IN ('rna_cds_cdna_start', 'rna_cds_cdna_end', 'rna_codon_table')),
      (SELECT count(*) FROM duckdb_tables() WHERE schema_name = 'main' AND table_name = 'model_rna_edits')
      FROM duckdb_columns() WHERE schema_name = 'main'")
    IFS='|' read -r has_circular has_cdna authority_count has_rna <<< "$metadata"
    [[ $authority_count == 0 || ( $authority_count == 3 && $has_cdna == 1 ) ]] || { echo "$source: incomplete RNA coding authority columns" >&2; exit 1; }
    [[ $has_rna == 0 || $has_cdna == 1 ]] || { echo "$source: model_rna_edits requires cdna_sequence" >&2; exit 1; }
    if [[ $has_rna == 0 ]]; then
        # DUCKVEP_TX_RNA_EDIT is bit 8 in duckvep_kernel.h.
        rna_flags=$("$duckdb_bin" -readonly -noheader -list "$source" -c 'SELECT count(*) FROM model_transcripts WHERE (transcript_flags & 256) != 0')
        [[ $rna_flags == 0 ]] || { echo "$source: RNA-edit flags require model_rna_edits" >&2; exit 1; }
    fi
    has_authority=0; [[ $authority_count == 3 ]] && has_authority=1
    [[ $("$duckdb_bin" -readonly -noheader -list "$source" -c 'SELECT count(*) FROM model_receipt') == 1 ]] || { echo "$source: model_receipt must have one row" >&2; exit 1; }
    identity=$("$duckdb_bin" -readonly -noheader -separator $'\t' -list "$source" -c "SELECT source_version, assembly FROM model_receipt")
    IFS=$'\t' read -r release assembly <<< "$identity"
    [[ $release =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && -n $assembly ]] || { echo "$source: invalid release or assembly in model_receipt" >&2; exit 1; }
    name=$(basename "$source" .duckdb); name=${name%_final}
    dir=ensembl-$release; mkdir -p "$stage/$dir"
    export_path=$stage/$dir/$name.duckdb
    [[ ! -e $export_path && ! -L $export_path && ! -e $export_path.wal && ! -L $export_path.wal ]] || { echo "duplicate export path: $export_path" >&2; exit 1; }
    if [[ -z $example || ( $assembly == GRCh38 && $example_assembly != GRCh38 ) ]]; then
        example=$dir/$name.duckdb; example_assembly=$assembly; example_regulation=$has_regulation
        example_circular=$has_circular; example_cdna=$has_cdna; example_authority=$has_authority; example_rna=$has_rna
    fi
    tables="model_regions model_transcripts"; options=""
    if [[ $has_regulation == 1 ]]; then
        tables+=" duckvep_regulation_features"; options=", {regulation_features_table: 'duckvep_regulation_features'}"
    fi
    [[ $has_rna == 1 ]] && tables+=" model_rna_edits"
    copy="SET preserve_insertion_order = true; LOAD $(sql_text "$ext"); ATTACH $(sql_text "$source") AS src (READ_ONLY);"
    copy+=" ATTACH $(sql_text "$export_path") AS out; USE out;"
    for table in $tables; do copy+=" CREATE TABLE $table AS SELECT * FROM src.$table;"; done
    copy+=" SET VARIABLE receipt = (SELECT duckvep_model_receipt_sql('model_regions', 'model_transcripts', source_name,"
    copy+=" source_version, assembly, source_manifest_sha256, reference_sha256, transcript_filter$options) FROM src.model_receipt);"
    copy+=" CREATE TABLE model_receipt AS SELECT * FROM query(getvariable('receipt'));"
    "$duckdb_bin" -unsigned :memory: -bail -c "$copy USE memory; DETACH out;" > /dev/null
    snapshot "$source" "$has_regulation" "$has_circular" "$has_cdna" "$has_authority" "$has_rna" "$work/source.dvsnap"
    snapshot "$export_path" "$has_regulation" "$has_circular" "$has_cdna" "$has_authority" "$has_rna" "$work/export.dvsnap"
    cmp -s "$work/source.dvsnap" "$work/export.dvsnap" || { echo "$name: exported model loads differently" >&2; exit 1; }
    rm -f "$work/source.dvsnap" "$work/export.dvsnap"
    IFS=$'\t' read -r transcripts regions regulation < <("$duckdb_bin" -readonly -noheader -separator $'\t' -list "$export_path" -c \
        "SELECT format('{:,}', transcript_count), format('{:,}', region_count), format('{:,}', coalesce(regulation_feature_count, 0)) FROM model_receipt")
    rows+=("| [\`$dir/$name.duckdb\`]($dir/$name.duckdb) | $assembly | $release | $transcripts | $regions | $regulation | $(du -h "$export_path" | cut -f1) |")
    if [[ -n $mane && $assembly == GRCh37 ]]; then
        (cd "$mane" && sha256sum --quiet -c SHA256SUMS) || { echo "MANE release files fail SHA256SUMS" >&2; exit 1; }
        model_key=$("$duckdb_bin" -readonly -noheader -list "$export_path" -c "SELECT model_sha256 FROM model_receipt")
        release_key=$("$duckdb_bin" -noheader -list :memory: -c \
            "SELECT string_agg(DISTINCT model_sha256, ',') FROM read_parquet($(sql_text "$mane/mane_grch37_mapping.parquet"))")
        [[ $model_key == "$release_key" ]] || { echo "$name: model_sha256 $model_key is not the MANE join key $release_key" >&2; exit 1; }
        mane_dir=$dir/mane-grch37; mane_model=$dir/$name.duckdb; mkdir -p "$stage/$mane_dir"
        cp "$mane"/mane_grch37_mapping.parquet "$mane"/grch37_transcript_authorities.parquet "$mane"/mane_grch37_receipt.csv "$stage/$mane_dir/"
    fi
    echo "$name: exported and checked" >&2
done
[[ -z $mane || -n ${mane_dir:-} ]] || { echo "--mane-grch37 needs a GRCh37 model" >&2; exit 1; }
(cd "$stage" && find . -type f \( -name '*.duckdb' -o -path './*/mane-grch37/*' \) -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS)

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

Transcript models for [DuckVEP](https://github.com/RGenomicsETL/DuckVEP), a DuckDB extension that annotates variant consequences the way Ensembl VEP does. Each file is a DuckDB database built from an Ensembl core database and the matching primary-assembly FASTA. It holds regions, transcripts with their exons and sequences, and optional regulatory features and RNA edits. The \`model_receipt\` table records the source release, transcript filter, source and reference hashes, counts, and \`model_sha256\`. External transcript relations use that hash as a join key. Use the file checksum to identify the complete database, including optional data.

| File | Assembly | Ensembl release | Transcripts | Regions | Regulatory features | Size |
|---|---|---|---:|---:|---:|---:|
EOF
printf '%s\n' "${rows[@]}"
cat <<EOF

\`SHA256SUMS\` lists the checksums of the data files.

## Use

Download a model, then load it once per DuckDB session:

\`\`\`sh
hf download $repo $example --repo-type dataset --local-dir models
\`\`\`

\`\`\`sql
LOAD 'duckvep.duckdb_extension';
ATTACH 'models/$example' AS src (READ_ONLY);
EOF
load_sql model "$example_regulation" "$example_circular" "$example_cdna" "$example_authority" "$example_rna"
echo '```'
echo
cat <<EOF
The source and export are compared as native snapshots loaded with the queries above. This check covers that load profile. The exporter loads \`model_rna_edits\` when it is present; it does not derive edits from other source tables. It does not supply a reference FASTA. The publisher must assert complete transcript coverage separately. The [function reference](https://github.com/RGenomicsETL/DuckVEP/blob/main/docs/functions.md) describes loading and annotation.

EOF
if [[ -n ${mane_dir:-} ]]; then cat <<EOF
## MANE on GRCh37

MANE is defined only on GRCh38. [\`$mane_dir/\`]($mane_dir) holds the DuckVEP data release [MANE v1.5 mapped to GRCh37](https://github.com/RGenomicsETL/DuckVEP/blob/main/design/duckvep_mane_grch37.md), validated against \`$mane_model\`. Its rows join to annotation rows on \`model_sha256\` and \`transcript_index\`, and that model's \`model_receipt\` holds the matching \`model_sha256\`. Label joined rows "MANE mapped to GRCh37", never native MANE. Rows with \`mapping_status\` \`cds_exact_utr_differs\` speak for the coding region only. The [caller recipe](https://github.com/RGenomicsETL/DuckVEP/blob/main/scripts/mane_grch37_caller.sql) shows the join.

EOF
fi
cat <<EOF
## Terms

The data come from [Ensembl](https://www.ensembl.org), which makes them available without restriction; see the [Ensembl disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html). When you use them, cite Ensembl and DuckVEP.
EOF
} > "$stage/README.md"
echo "staged in $stage" >&2
[[ $dry_run == 1 ]] && exit 0

hf repos create "$repo" --repo-type dataset "$visibility" --exist-ok > /dev/null
hf upload "$repo" "$stage" . --repo-type dataset --exclude '*.wal' --commit-message "Publish DuckVEP models"
