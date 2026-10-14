#include "duckvep_core_products.h"
#include <string.h>

const char *const duckvep_product_input_type =
    "STRUCT(reference STRUCT(start1 UBIGINT, dna VARCHAR, strand TINYINT), "
    "exons STRUCT(start1 UBIGINT, end1 UBIGINT)[], "
    "blocks STRUCT(start1 UBIGINT, end1 UBIGINT)[], "
    "sources STRUCT(source_id UBIGINT, start1 UBIGINT, reference VARCHAR, alternate VARCHAR)[], "
    "programme STRUCT(id VARCHAR, evidence_id VARCHAR, assumption_id VARCHAR, "
    "expected_reference_rna VARCHAR, annotated_start1 UINTEGER, alternative_start1 UINTEGER, "
    "alternative_trigger_id VARCHAR, genetic_code UTINYINT))";
const char *const duckvep_product_output_type =
    "STRUCT(sources STRUCT(source_id UBIGINT, start1 UBIGINT, reference VARCHAR, alternate VARCHAR)[], "
    "source_effects STRUCT(source_id UBIGINT, changed_dna BOOLEAN, changed_rna_bases UINTEGER, "
    "genomic_start1 UBIGINT, genomic_end1 UBIGINT, rna_first1 UINTEGER, rna_last1 UINTEGER)[], "
    "reference_dna VARCHAR, allele_evidence VARCHAR, programme_id VARCHAR, "
    "programme_evidence_id VARCHAR, programme_assumption_id VARCHAR, product_evidence VARCHAR, "
    "alternate_dna VARCHAR, reference_rna VARCHAR, alternate_rna VARCHAR, "
    "reference_peptide VARCHAR, alternate_peptide VARCHAR, reference_stop_codon1 UINTEGER, "
    "alternate_stop_codon1 UINTEGER, annotated_start1 UINTEGER, effective_start1 UINTEGER, "
    "last_mapped_junction1 UINTEGER, verified_sources UBIGINT, conditional_nmd VARCHAR, "
    "conditional_nmd_reason VARCHAR, alternative_start1 UINTEGER, alternative_trigger_id VARCHAR, "
    "genetic_code UTINYINT, strand TINYINT, reference_start1 UBIGINT)";
const char *const duckvep_product_output_names[PRODUCT_OUTPUT_COUNT] = {
    "sources", "source_effects", "reference_dna", "allele_evidence", "programme_id",
    "programme_evidence_id", "programme_assumption_id", "product_evidence", "alternate_dna",
    "reference_rna", "alternate_rna", "reference_peptide", "alternate_peptide",
    "reference_stop_codon1", "alternate_stop_codon1", "annotated_start1", "effective_start1",
    "last_mapped_junction1", "verified_sources", "conditional_nmd", "conditional_nmd_reason",
    "alternative_start1", "alternative_trigger_id", "genetic_code", "strand", "reference_start1"
};
const char *const duckvep_product_output_types[PRODUCT_OUTPUT_COUNT] = {
    "LIST", "LIST", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR",
    "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "UINTEGER", "UINTEGER",
    "UINTEGER", "UINTEGER", "UINTEGER", "UBIGINT", "VARCHAR", "VARCHAR", "UINTEGER",
    "VARCHAR", "UTINYINT", "TINYINT", "UBIGINT"
};

int duckvep_core_product_output(const duckvep_core_product_fields_t *in,
    const duckvep_core_product_result_t *r, size_t index,
    const char **text, size_t *length, uint64_t *number) {
    const duckvep_profile_execution_result_t *o = &r->result;
    *text = NULL; *length = 0; *number = 0;
    switch (index) {
    case 2: *text = in->dna.ptr; *length = in->dna.length; break;
    case 3: *text = "hypothetical_assignment"; *length = strlen(*text); break;
    case 4: *text = r->identity; *length = strlen(*text); break;
    case 5: *text = r->evidence_id; *length = strlen(*text); break;
    case 6: *text = r->assumption_id; *length = strlen(*text); break;
    case 7: *text = "supplied_programme"; *length = strlen(*text); break;
    case 8: *text = (const char *)r->dna; *length = in->dna.length; break;
    case 9: *text = (const char *)r->reference_rna; *length = o->rna_length; break;
    case 10: *text = (const char *)r->alternate_rna; *length = o->rna_length; break;
    case 11: *text = (const char *)r->reference_peptide; *length = o->reference_peptide_length; break;
    case 12: *text = (const char *)r->alternate_peptide; *length = o->alternate_peptide_length; break;
    case 13: *number = o->reference_stop1; break;
    case 14: *number = o->alternate_stop1; break;
    case 15: *number = o->reference_start1; break;
    case 16: *number = o->alternate_start1; break;
    case 17: *number = o->last_junction1; break;
    case 18: *number = o->verified_sources; break;
    case 19: *text = "unknown"; *length = strlen(*text); break;
    case 20: *text = "shared NMD eligibility and complete CDS termination not established";
             *length = strlen(*text); break;
    case 21: *number = in->alternative_start1; break;
    case 22: *text = in->trigger_id.ptr; *length = in->trigger_id.length; break;
    case 23: *number = in->genetic_code; break;
    case 24: *number = (uint8_t)in->strand; break;
    case 25: *number = in->reference_start1; break;
    default: return 0;
    }
    return 1;
}
