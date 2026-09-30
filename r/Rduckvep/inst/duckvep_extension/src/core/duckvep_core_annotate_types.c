#include "duckvep_core_annotate_types.h"

static const char *const base_names[] = {
	"transcript_index", "gene_index", "consequence", "impact",
	"region", "status", "reason", "cdna_position", "cds_position",
	"protein_position", "reference_amino_acid",
	"alternate_amino_acid", "nmd_prediction",
	"nmd_escape_intronless", "nmd_escape_early_cds",
	"nmd_escape_last_exon", "nmd_escape_penultimate_exon_end",
	"regulation_feature_index", "overlap_object",
	"consequence_mask", "region_mask", "impact_code", "status_code",
	"reason_code", "reference_amino_acid_code",
	"alternate_amino_acid_code", "nmd_prediction_code",
	"nmd_escape_reasons", "overlap_object_code"
};
static const char *const hgvs_names[] = {
	"transcript_hgvs", "protein_hgvs", "hgvs_shift",
	"transcript_hgvs_status", "transcript_hgvs_reason",
	"protein_hgvs_status", "protein_hgvs_reason"
};
static const char *const projection_names[] = {
	"output_allele", "interbase", "cdna_start", "cdna_end",
	"cds_start", "cds_end", "protein_start", "protein_end",
	"exon_first", "exon_last", "exon_total", "intron_first",
	"intron_last", "intron_total", "transcript_distance",
	"cds_start_nf", "cds_end_nf", "reference_amino_acids",
	"alternate_amino_acids", "reference_codons",
	"alternate_codons"
};
static const char *const compact_names[] = {
	"transcript_index", "gene_index", "consequence_mask",
	"region_mask", "impact_code", "status_code", "reason_code",
	"cdna_position", "cds_position", "protein_position",
	"reference_amino_acid_code", "alternate_amino_acid_code",
	"nmd_prediction_code", "nmd_escape_reasons",
	"regulation_feature_index", "overlap_object_code"
};
static const char *const compact_hgvs_names[] = {
	"transcript_hgvs", "protein_hgvs", "hgvs_shift",
	"transcript_hgvs_status", "transcript_hgvs_reason",
	"protein_hgvs_status", "protein_hgvs_reason"
};

size_t
duckvep_core_result_columns(duckvep_result_kind_t kind, int with_hgvs, int with_projection,
	duckvep_result_column_t out[DUCKVEP_RESULT_COLUMNS_MAX])
{
	size_t index, count = 0;

	if (kind == DUCKVEP_RESULT_RICH) {
		for (index = 0u; index < 29u; index++) {
			out[index].name = base_names[index];
			if (index <= 1u || (index >= 7u && index <= 9u) || index == 17u ||
			    index == 20u)
				out[index].type = DUCKVEP_COLUMN_UINTEGER;
			else if ((index >= 2u && index <= 6u) ||
			    (index >= 10u && index <= 12u) || index == 18u)
				out[index].type = DUCKVEP_COLUMN_VARCHAR;
			else if (index >= 13u && index <= 16u)
				out[index].type = DUCKVEP_COLUMN_BOOLEAN;
			else if (index == 19u)
				out[index].type = DUCKVEP_COLUMN_UBIGINT;
			else
				out[index].type = DUCKVEP_COLUMN_UTINYINT;
		}
		count = 29u;
		if (with_hgvs) {
			for (index = 0u; index < 7u; index++, count++) {
				out[count].name = hgvs_names[index];
				out[count].type = index == 2u ? DUCKVEP_COLUMN_UINTEGER : DUCKVEP_COLUMN_VARCHAR;
			}
		}
		if (with_projection) {
			for (index = 0u; index < 21u; index++, count++) {
				out[count].name = projection_names[index];
				if (index == 0u || index >= 17u)
					out[count].type = DUCKVEP_COLUMN_VARCHAR;
				else if (index == 1u || index == 15u || index == 16u)
					out[count].type = DUCKVEP_COLUMN_BOOLEAN;
				else if (index == 14u)
					out[count].type = DUCKVEP_COLUMN_UBIGINT;
				else
					out[count].type = DUCKVEP_COLUMN_UINTEGER;
			}
		}
		return count;
	}
	for (index = 0; index < 16; index++) {
		out[index].name = compact_names[index];
		out[index].type = index == 0 || index == 1 || index == 3 || (index >= 7 && index <= 9) ||
		    index == 14 ? DUCKVEP_COLUMN_UINTEGER : index == 2 ? DUCKVEP_COLUMN_UBIGINT :
		    DUCKVEP_COLUMN_UTINYINT;
	}
	count = 16;
	if (kind == DUCKVEP_RESULT_COMPACT_HGVS) {
		for (index = 0; index < 7; index++, count++) {
			out[count].name = compact_hgvs_names[index];
			out[count].type = index == 2 ? DUCKVEP_COLUMN_UINTEGER : DUCKVEP_COLUMN_VARCHAR;
		}
	}
	return count;
}

const char *
duckvep_core_column_type_name(duckvep_column_type_t type)
{
	switch (type) {
	case DUCKVEP_COLUMN_UINTEGER: return "UINTEGER";
	case DUCKVEP_COLUMN_VARCHAR: return "VARCHAR";
	case DUCKVEP_COLUMN_BOOLEAN: return "BOOLEAN";
	case DUCKVEP_COLUMN_UBIGINT: return "UBIGINT";
	default: return "UTINYINT";
	}
}
