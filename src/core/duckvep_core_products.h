#ifndef DUCKVEP_CORE_PRODUCTS_H
#define DUCKVEP_CORE_PRODUCTS_H
#include "kernel/src/duckvep_product_profile.h"
#include <stdint.h>

typedef struct { const char *ptr; size_t length; } duckvep_product_text_t;
typedef struct { uint64_t start1, end1; } duckvep_product_interval_t;
typedef struct {
    uint64_t source_id, start1;
    duckvep_product_text_t reference, alternate;
} duckvep_product_source_t;
typedef struct {
    uint64_t reference_start1;
    duckvep_product_text_t dna;
    int8_t strand;
    const duckvep_product_interval_t *exons, *blocks;
    size_t exon_count, block_count;
    const duckvep_product_source_t *sources;
    size_t source_count;
    duckvep_product_text_t id, evidence_id, assumption_id, expected_rna, trigger_id;
    uint32_t annotated_start1, alternative_start1;
    uint8_t genetic_code;
} duckvep_core_product_fields_t;

/* The two LIST columns are written by each host from the owned source views
 * and effects; all other fields use this ordered scalar table. */
enum { PRODUCT_OUTPUT_COUNT = 26, PRODUCT_SOURCES = 0, PRODUCT_EFFECTS = 1 };
extern const char *const duckvep_product_output_names[PRODUCT_OUTPUT_COUNT];
extern const char *const duckvep_product_output_types[PRODUCT_OUTPUT_COUNT];
extern const char *const duckvep_product_input_type;
extern const char *const duckvep_product_output_type;
typedef struct {
    uint64_t source_id, genomic_start1, genomic_end1;
    uint32_t rna_first1, rna_last1, changed_rna_bases;
    uint8_t changed_dna;
} duckvep_core_product_effect_t;
typedef struct {
    duckvep_profile_execution_result_t result;
    duckvep_core_product_effect_t *effects;
    uint8_t *dna, *reference_rna, *alternate_rna;
    uint8_t *reference_peptide, *alternate_peptide;
    char *identity, *evidence_id, *assumption_id;
} duckvep_core_product_result_t;
int duckvep_core_product_output(const duckvep_core_product_fields_t *input,
    const duckvep_core_product_result_t *result, size_t index,
    const char **text, size_t *length, uint64_t *number);
duckvep_profile_status_t duckvep_core_product_run(
    const duckvep_profile_execution_t *input, duckvep_core_product_result_t *output);
void duckvep_core_product_dispose(duckvep_core_product_result_t *output);
duckvep_profile_status_t duckvep_core_product_scenario(
    const duckvep_core_product_fields_t *fields, duckvep_core_product_result_t *output);
#endif
