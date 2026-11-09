#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
duckdb_bin=${DUCKDB:-duckdb}
extension=${DUCKVEP_EXTENSION:-$root/build/release/duckvep.duckdb_extension}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
source_db=$work/grch37.duckdb

"$duckdb_bin" -unsigned "$source_db" >/dev/null <<'SQL'
CREATE TABLE model_regions AS SELECT 0::UINTEGER seq_region, '1'::VARCHAR seq_region_name,
  1000::UBIGINT sequence_length, 1::UINTEGER source_seq_region_id,
  1::UINTEGER source_coord_system_id, 'chromosome'::VARCHAR coord_system_name, 1::INTEGER coord_system_rank;
CREATE TABLE model_transcripts AS SELECT
  0::UINTEGER transcript_index, 0::UINTEGER seq_region, 10::UBIGINT transcript_start,
  20::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index,
  0::UBIGINT transcript_flags, NULL::UBIGINT cds_start, NULL::UBIGINT cds_end,
  NULL::BLOB cds_sequence, NULL::UTINYINT codon_table, NULL::BLOB pre_cds_sequence,
  NULL::BLOB post_cds_sequence, '1'::VARCHAR seq_region_name, 1000::UBIGINT sequence_length,
  1::UINTEGER source_seq_region_id, 1::UINTEGER source_transcript_id,
  'TEST1'::VARCHAR transcript_stable_id, 1::USMALLINT transcript_version,
  'lncRNA'::VARCHAR transcript_biotype, NULL::VARCHAR mane_select_refseq,
  NULL::VARCHAR mane_plus_clinical_refseq, 1::UINTEGER source_gene_id,
  'GENE1'::VARCHAR gene_stable_id, 1::USMALLINT gene_version, 'lncRNA'::VARCHAR gene_biotype,
  NULL::UINTEGER source_translation_id, NULL::VARCHAR translation_stable_id,
  NULL::USMALLINT translation_version, NULL::VARCHAR sequence_withheld_reason,
  [struct_pack(rank := 1::UINTEGER, exon_start := 10::UBIGINT, exon_end := 20::UBIGINT,
    exon_cdna_start := 1::UBIGINT, exon_cdna_end := 11::UBIGINT, phase := -1::TINYINT,
    end_phase := -1::TINYINT, source_exon_id := 1::UINTEGER, exon_stable_id := 'EXON1'::VARCHAR)] exons,
  []::STRUCT(cdna_start UBIGINT, cdna_end UBIGINT, mature_mirna_start UBIGINT, mature_mirna_end UBIGINT)[] mature_mirna_regions,
  []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] peptide_edits;
CREATE TABLE model_receipt AS SELECT 'Ensembl'::VARCHAR source_name, '116'::VARCHAR source_version,
  'GRCh37'::VARCHAR assembly, repeat('a', 64)::VARCHAR source_manifest_sha256,
  repeat('b', 64)::VARCHAR reference_sha256, 'all'::VARCHAR transcript_filter;
SQL

export DUCKDB="$duckdb_bin" DUCKVEP_EXTENSION="$extension"
stage=$work/stage
STAGE="$stage" bash "$root/scripts/publish-models.sh" --dry-run example/models "$source_db"
(cd "$stage" && sha256sum -c SHA256SUMS >/dev/null)
grep -F 'ensembl-116/grch37.duckdb' "$stage/README.md" >/dev/null
grep -F 'transcript_coverage_complete := FALSE' "$stage/README.md" >/dev/null
[[ -z $(find "$stage" -name '*.dvsnap' -print -quit) ]]

stage=$work/complete
STAGE="$stage" bash "$root/scripts/publish-models.sh" --dry-run --complete-coverage example/models "$source_db"
grep -F 'transcript_coverage_complete := TRUE' "$stage/README.md" >/dev/null

cp "$source_db" "$work/optional.duckdb"
"$duckdb_bin" "$work/optional.duckdb" >/dev/null <<'SQL'
ALTER TABLE model_regions ADD COLUMN circular BOOLEAN;
UPDATE model_regions SET circular = true;
ALTER TABLE model_transcripts ADD COLUMN cdna_sequence BLOB;
ALTER TABLE model_transcripts ADD COLUMN rna_cds_cdna_start UBIGINT;
ALTER TABLE model_transcripts ADD COLUMN rna_cds_cdna_end UBIGINT;
ALTER TABLE model_transcripts ADD COLUMN rna_codon_table UTINYINT;
UPDATE model_transcripts SET cdna_sequence = 'AAAAAAAAAAA'::BLOB, transcript_flags = 256;
CREATE TABLE model_rna_edits AS SELECT 0::UINTEGER transcript_index, 2::UINTEGER cdna_start,
  2::UINTEGER cdna_end, 'A'::VARCHAR expected_substrate, 'G'::VARCHAR alternate_sequence,
  'substitution'::VARCHAR operation_kind, 'Ensembl116:_rna_edit'::VARCHAR programme_name,
  'transcript_attrib:SeqEdit'::VARCHAR evidence_name, 0::UTINYINT boundary_overrides;
CREATE TABLE duckvep_regulation_features (
  regulation_feature_index UINTEGER, seq_region UINTEGER, feature_start UBIGINT, feature_end UBIGINT,
  feature_kind VARCHAR, source_table VARCHAR, regulatory_feature_id UINTEGER,
  regulatory_stable_id VARCHAR, feature_type_name VARCHAR, feature_type_class VARCHAR,
  motif_feature_id UINTEGER, motif_stable_id VARCHAR, motif_name VARCHAR, motif_matrix_id UINTEGER,
  motif_strand TINYINT, motif_score DOUBLE, source_seq_region_id UINTEGER
);
INSERT INTO duckvep_regulation_features (regulation_feature_index, seq_region, feature_start,
  feature_end, feature_kind, source_table, regulatory_stable_id)
  VALUES (0, 0, 15, 18, '1', 'regulatory_feature', 'ENSR1');
SQL
stage=$work/optional
STAGE="$stage" bash "$root/scripts/publish-models.sh" --dry-run example/models "$work/optional.duckdb"
grep -F 'seq_region_name, circular FROM src.model_regions' "$stage/README.md" >/dev/null
grep -F 'cdna_sequence, rna_cds_cdna_start, rna_cds_cdna_end, rna_codon_table' "$stage/README.md" >/dev/null
grep -F 'rna_edit_query :=' "$stage/README.md" >/dev/null
grep -F 'interval_feature_query :=' "$stage/README.md" >/dev/null
[[ $("$duckdb_bin" -readonly -noheader -list "$stage/ensembl-116/optional.duckdb" -c 'SELECT count(*) FROM model_rna_edits') == 1 ]]

cp "$work/optional.duckdb" "$work/missing-rna.duckdb"
"$duckdb_bin" "$work/missing-rna.duckdb" -c 'DROP TABLE model_rna_edits' >/dev/null
if STAGE="$work/missing-rna" bash "$root/scripts/publish-models.sh" --dry-run example/models "$work/missing-rna.duckdb" >"$work/missing-rna.log" 2>&1; then
  echo 'missing RNA programme was accepted' >&2; exit 1
fi
grep -F 'RNA-edit flags require model_rna_edits' "$work/missing-rna.log" >/dev/null

stage=$work/overlap
mkdir -p "$stage/ensembl-116"
cp "$source_db" "$stage/ensembl-116/grch37.duckdb"
before=$(sha256sum "$stage/ensembl-116/grch37.duckdb")
if STAGE="$stage" bash "$root/scripts/publish-models.sh" --dry-run example/models "$stage/ensembl-116/grch37.duckdb" >"$work/overlap.log" 2>&1; then
  echo 'nonempty stage was accepted' >&2; exit 1
fi
grep -F 'STAGE must be an empty directory' "$work/overlap.log" >/dev/null
[[ $before == "$(sha256sum "$stage/ensembl-116/grch37.duckdb")" ]]

stage=$work/duplicate
if STAGE="$stage" bash "$root/scripts/publish-models.sh" --dry-run example/models "$source_db" "$source_db" >"$work/duplicate.log" 2>&1; then
  echo 'duplicate export path was accepted' >&2; exit 1
fi
grep -F 'duplicate export path' "$work/duplicate.log" >/dev/null
[[ -f $stage/ensembl-116/grch37.duckdb ]]
echo 'model publication staging tests passed'
