#include "duckvep_core_products.h"
#include "kernel/src/duckvep_budget.h"
#include <stdint.h>
#include <string.h>

void duckvep_core_product_dispose(duckvep_core_product_result_t *out) {
    if (!out) return;
    duckvep_budget_free(out->effects);
    duckvep_budget_free(out->dna);
    duckvep_budget_free(out->reference_rna);
    duckvep_budget_free(out->alternate_rna);
    duckvep_budget_free(out->reference_peptide);
    duckvep_budget_free(out->alternate_peptide);
    duckvep_budget_free(out->identity);
    duckvep_budget_free(out->evidence_id);
    duckvep_budget_free(out->assumption_id);
    *out = (duckvep_core_product_result_t){0};
}

/* Product and workspace are transaction-local. The caller retains ownership
 * of reference, programme, exon/block and source arrays. */
duckvep_profile_status_t duckvep_core_product_run(
    const duckvep_profile_execution_t *input, duckvep_core_product_result_t *output) {
    if (!output || !input || !input->geometry || !input->programme)
        return DUCKVEP_PROFILE_INVALID;
    const size_t dna = input->geometry->dna_length, rna = input->programme->rna_length;
    const size_t sources = input->source_count;
    if (!dna || !rna || dna > UINT32_MAX || rna > UINT32_MAX ||
        dna == SIZE_MAX || rna == SIZE_MAX ||
        sources > SIZE_MAX / sizeof(duckvep_profile_source_effect_t) ||
        sources > SIZE_MAX / sizeof(duckvep_core_product_effect_t) ||
        sources > SIZE_MAX / sizeof(duckvep_haplotype_edit_t))
        return DUCKVEP_PROFILE_INVALID;
    duckvep_profile_workspace_t w = {0};
    duckvep_core_product_result_t owned = {0};
    duckvep_profile_execution_result_t result;
    size_t peptide = rna / 3u + 1u;
    w.dna_capacity = dna;
    w.rna_capacity = rna;
    w.peptide_capacity = peptide;
    w.source_effect_capacity = sources;
    w.edit_capacity = sources;
    w.dna = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, dna);
    w.reference_rna = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, rna);
    w.alternate_rna = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, rna);
    w.reference_peptide = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, peptide);
    w.alternate_peptide = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, peptide);
    if (sources) {
        w.source_effects = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,
            sources * sizeof *w.source_effects);
        w.edits = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,
            sources * sizeof *w.edits);
    }
    duckvep_profile_status_t status = DUCKVEP_PROFILE_CAPACITY;
    if (!w.dna || !w.reference_rna || !w.alternate_rna || !w.reference_peptide ||
        !w.alternate_peptide || (sources && (!w.source_effects || !w.edits))) goto done;
    status = duckvep_product_execute_v1(input, &w, &result);
    if (status != DUCKVEP_PROFILE_OK) goto done;
    owned.dna = duckvep_budget_malloc(DUCKVEP_OWNER_EMIT, dna + 1u);
    owned.reference_rna = duckvep_budget_malloc(DUCKVEP_OWNER_EMIT, rna + 1u);
    owned.alternate_rna = duckvep_budget_malloc(DUCKVEP_OWNER_EMIT, rna + 1u);
    owned.reference_peptide = duckvep_budget_malloc(DUCKVEP_OWNER_EMIT, result.reference_peptide_length + 1u);
    owned.alternate_peptide = duckvep_budget_malloc(DUCKVEP_OWNER_EMIT, result.alternate_peptide_length + 1u);
    owned.identity = duckvep_budget_strdup(DUCKVEP_OWNER_EMIT, result.identity);
    owned.evidence_id = duckvep_budget_strdup(DUCKVEP_OWNER_EMIT, result.evidence_id);
    owned.assumption_id = duckvep_budget_strdup(DUCKVEP_OWNER_EMIT, result.assumption_id);
    if (sources) owned.effects = duckvep_budget_malloc(DUCKVEP_OWNER_EMIT,
        sources * sizeof *owned.effects);
    status = DUCKVEP_PROFILE_CAPACITY;
    if (!owned.dna || !owned.reference_rna || !owned.alternate_rna ||
        !owned.reference_peptide || !owned.alternate_peptide ||
        !owned.identity || !owned.evidence_id || !owned.assumption_id ||
        (sources && !owned.effects)) goto done;
    memcpy(owned.dna, w.dna, dna); owned.dna[dna] = 0;
    memcpy(owned.reference_rna, w.reference_rna, rna); owned.reference_rna[rna] = 0;
    memcpy(owned.alternate_rna, w.alternate_rna, rna); owned.alternate_rna[rna] = 0;
    memcpy(owned.reference_peptide, w.reference_peptide, result.reference_peptide_length);
    owned.reference_peptide[result.reference_peptide_length] = 0;
    memcpy(owned.alternate_peptide, w.alternate_peptide, result.alternate_peptide_length);
    owned.alternate_peptide[result.alternate_peptide_length] = 0;
    for (size_t i = 0; i < sources; ++i) {
        const duckvep_profile_source_effect_t *e = &w.source_effects[i];
        owned.effects[i] = (duckvep_core_product_effect_t){e->source_id,
            e->genomic_start1, e->genomic_end1, e->rna_first1, e->rna_last1,
            e->changed_rna_bases, e->changed_dna};
    }
    owned.result = result;
    owned.result.identity = owned.identity;
    owned.result.evidence_id = owned.evidence_id;
    owned.result.assumption_id = owned.assumption_id;
    *output = owned;
    owned = (duckvep_core_product_result_t){0};
    status = DUCKVEP_PROFILE_OK;
done:
    duckvep_budget_free(w.dna);
    duckvep_budget_free(w.reference_rna);
    duckvep_budget_free(w.alternate_rna);
    duckvep_budget_free(w.reference_peptide);
    duckvep_budget_free(w.alternate_peptide);
    duckvep_budget_free(w.source_effects);
    duckvep_budget_free(w.edits);
    duckvep_core_product_dispose(&owned);
    return status;
}
