#include "duckvep_core_products.h"
#include "kernel/src/duckvep_budget.h"
#include <limits.h>
#include <stdbool.h>
#include <string.h>

static bool local(uint64_t position, uint64_t base, uint64_t span, uint32_t *out) {
    if (position < base || position - base >= span) return false;
    *out = (uint32_t)(position - base + 1u);
    return true;
}
static bool interval(duckvep_product_interval_t in, uint64_t base, uint64_t span,
                     duckvep_profile_interval_t *out) {
    return local(in.start1, base, span, &out->start1) &&
        local(in.end1, base, span, &out->end1) && out->start1 <= out->end1;
}
static char *copy_id(duckvep_product_text_t in) {
    if (!in.ptr || !in.length || in.length == SIZE_MAX || memchr(in.ptr, 0, in.length)) return NULL;
    char *out = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, in.length + 1u);
    if (out) { memcpy(out, in.ptr, in.length); out[in.length] = 0; }
    return out;
}

duckvep_profile_status_t duckvep_core_product_scenario(
    const duckvep_core_product_fields_t *f, duckvep_core_product_result_t *out) {
    if (!f || !out || !f->reference_start1 || !f->dna.ptr || !f->dna.length ||
        f->dna.length > UINT32_MAX || f->reference_start1 - 1u > UINT64_MAX - f->dna.length || (f->strand != 1 && f->strand != -1) ||
        !f->exon_count || !f->block_count || !f->source_count ||
        !f->exons || !f->blocks || !f->sources || !f->expected_rna.ptr ||
        !f->expected_rna.length || f->expected_rna.length > UINT32_MAX ||
        !f->annotated_start1 || f->annotated_start1 > f->expected_rna.length ||
        (f->alternative_start1 && (f->trigger_id.length != sizeof("annotated-start-lost")-1u ||
         !f->trigger_id.ptr || memcmp(f->trigger_id.ptr, "annotated-start-lost", f->trigger_id.length))) ||
        (!f->alternative_start1 && f->trigger_id.length) ||
        f->exon_count > SIZE_MAX / sizeof(duckvep_profile_interval_t) ||
        f->block_count > SIZE_MAX / sizeof(duckvep_profile_rna_block_v1_t) ||
        f->source_count > SIZE_MAX / sizeof(duckvep_haplotype_edit_t) ||
        f->source_count > SIZE_MAX / sizeof(uint64_t)) return DUCKVEP_PROFILE_INVALID;
    duckvep_profile_status_t status = DUCKVEP_PROFILE_INVALID;
    duckvep_profile_interval_t *exons = NULL;
    duckvep_profile_rna_block_v1_t *blocks = NULL;
    duckvep_haplotype_edit_t *edits = NULL;
    uint64_t *ids = NULL;
    char *names[3] = {0};
    exons = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, f->exon_count * sizeof *exons);
    blocks = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, f->block_count * sizeof *blocks);
    edits = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, f->source_count * sizeof *edits);
    ids = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, f->source_count * sizeof *ids);
    if (!exons || !blocks || !edits || !ids) { status = DUCKVEP_PROFILE_CAPACITY; goto done; }
    for (size_t i = 0; i < f->exon_count; ++i)
        if (!interval(f->exons[i], f->reference_start1, f->dna.length, &exons[i])) goto done;
    uint64_t rna_offset = 0;
    for (size_t i = 0; i < f->block_count; ++i) {
        if (!interval(f->blocks[i], f->reference_start1, f->dna.length,
                      &blocks[i].genomic_source)) goto done;
        uint64_t width = (uint64_t)blocks[i].genomic_source.end1 - blocks[i].genomic_source.start1 + 1u;
        if (width > f->expected_rna.length - rna_offset) { status = DUCKVEP_PROFILE_PARTIAL_MODEL; goto done; }
        blocks[i].rna_start1 = (uint32_t)(rna_offset + 1u);
        rna_offset += width;
    }
    if (rna_offset != f->expected_rna.length) { status = DUCKVEP_PROFILE_PARTIAL_MODEL; goto done; }
    for (size_t i = 0; i < f->source_count; ++i) {
        const duckvep_product_source_t *s = &f->sources[i];
        uint32_t start;
        if (!s->reference.ptr || !s->alternate.ptr ||
            !s->reference.length || s->reference.length > UINT32_MAX ||
            s->reference.length != s->alternate.length ||
            !local(s->start1, f->reference_start1, f->dna.length, &start) ||
            s->reference.length > f->dna.length - (size_t)(start - 1u)) goto done;
        edits[i] = (duckvep_haplotype_edit_t){.cds_start=start,
            .ref=(const uint8_t *)s->reference.ptr, .ref_len=(uint32_t)s->reference.length,
            .alt=(const uint8_t *)s->alternate.ptr, .alt_len=(uint32_t)s->alternate.length,
            .variant_strand=1};
        ids[i] = s->source_id;
    }
    duckvep_product_text_t labels[] = {f->id, f->evidence_id, f->assumption_id};
    for (size_t i = 0; i < 3; ++i)
        if (!labels[i].ptr || !labels[i].length || labels[i].length == SIZE_MAX ||
            memchr(labels[i].ptr, 0, labels[i].length)) goto done;
    for (size_t i = 0; i < 3; ++i) {
        names[i] = copy_id(labels[i]);
        if (!names[i]) { status = DUCKVEP_PROFILE_CAPACITY; goto done; }
    }
    /* The kernel consumes descending genomic edits. Retain caller order in the
     * public source/effect lists by restoring effects by unique source ID. */
    for (size_t i = 1; i < f->source_count; ++i) {
        duckvep_haplotype_edit_t edit = edits[i];
        uint64_t id = ids[i];
        size_t j = i;
        while (j && edits[j-1].cds_start < edit.cds_start) {
            edits[j] = edits[j-1]; ids[j] = ids[j-1]; --j;
        }
        edits[j] = edit; ids[j] = id;
    }
    duckvep_profile_v1_input_t geometry = {
        .reference_start1=1, .reference_dna=(const uint8_t *)f->dna.ptr,
        .dna_length=f->dna.length, .transcript={1, (uint32_t)f->dna.length},
        .strand=f->strand, .exons=exons, .exon_count=f->exon_count, .model_complete=1u
    };
    duckvep_profile_programme_v1_t programme = {
        .identity=names[0], .evidence_id=names[1], .assumption_id=names[2],
        .blocks=blocks, .block_count=f->block_count,
        .reference_rna=(const uint8_t *)f->expected_rna.ptr,
        .rna_length=f->expected_rna.length,
        .cds_start1=f->annotated_start1, .cds_stop_end1=(uint32_t)f->expected_rna.length,
        .annotated_start1=f->annotated_start1, .alternative_start1=f->alternative_start1,
        .trigger_id=f->alternative_start1 ? "annotated-start-lost" : NULL,
        .counterfactual_id=f->alternative_start1 ? "reference-annotated-start-retained" : NULL
    };
    duckvep_profile_execution_t input = {
        .geometry=&geometry, .programme=&programme, .sources=edits,
        .source_ids=ids, .source_count=f->source_count, .codon_table=f->genetic_code
    };
    status = duckvep_core_product_run(&input, out);
    if (status == DUCKVEP_PROFILE_OK) {
        for (size_t i = 0; i < f->source_count; ++i) {
            size_t j = i;
            while (j < f->source_count && out->effects[j].source_id != f->sources[i].source_id) ++j;
            if (j == f->source_count) { status = DUCKVEP_PROFILE_INVALID; break; }
            duckvep_core_product_effect_t tmp = out->effects[i];
            out->effects[i] = out->effects[j]; out->effects[j] = tmp;
            out->effects[i].genomic_start1 += f->reference_start1 - 1u;
            out->effects[i].genomic_end1 += f->reference_start1 - 1u;
        }
        if (status != DUCKVEP_PROFILE_OK) duckvep_core_product_dispose(out);
    }
done:
    for (size_t i = 0; i < 3; ++i) duckvep_budget_free(names[i]);
    duckvep_budget_free(exons); duckvep_budget_free(blocks);
    duckvep_budget_free(edits); duckvep_budget_free(ids);
    return status;
}
