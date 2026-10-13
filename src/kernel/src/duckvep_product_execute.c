#include "duckvep_product_profile.h"
#include "duckvep_haplotype_stream.h"
#include <limits.h>
#include <string.h>

static uint8_t complement(uint8_t b) {
    switch (b) {
    case 'A': return 'T'; case 'T': return 'A';
    case 'C': return 'G'; case 'G': return 'C';
    default: return 0;
    }
}

static int span_disjoint(const void *a, size_t n, const void *b, size_t m) {
    uintptr_t x = (uintptr_t)a, y = (uintptr_t)b;
    return x <= UINTPTR_MAX - n && y <= UINTPTR_MAX - m &&
        (x + n <= y || y + m <= x);
}

duckvep_profile_status_t duckvep_product_execute_v1(
    const duckvep_profile_execution_t *q,
    duckvep_profile_workspace_t *w,
    duckvep_profile_execution_result_t *out) {
    if (!q || !w || !out || !q->geometry || !q->programme)
        return DUCKVEP_PROFILE_INVALID;
    const duckvep_profile_v1_input_t *g = q->geometry;
    const duckvep_profile_programme_v1_t *p = q->programme;
    /* This bounded profile rejects overlapping source representations rather
     * than treating a repeated or aliased edit as an independent joint allele. */
    if ((q->source_count && (!q->sources || !q->source_ids)) ||
        q->source_count > SIZE_MAX / sizeof *q->sources ||
        q->source_count > SIZE_MAX / sizeof *q->source_ids ||
        !g->reference_dna || !g->reference_start1 || !g->dna_length ||
        g->dna_length > UINT32_MAX ||
        (uint64_t)g->reference_start1 + g->dna_length - 1u > UINT32_MAX ||
        !g->model_complete || !g->exons || !g->exon_count ||
        g->exon_count > SIZE_MAX / sizeof *g->exons ||
        (g->strand != 1 && g->strand != -1) ||
        !p->identity || !*p->identity || !p->evidence_id || !*p->evidence_id ||
        !p->assumption_id || !*p->assumption_id || !p->blocks || !p->block_count ||
        p->block_count > SIZE_MAX / sizeof *p->blocks ||
        !p->reference_rna || !p->rna_length || p->rna_length > UINT32_MAX ||
        !p->cds_start1 || p->cds_start1 != p->annotated_start1 ||
        p->cds_stop_end1 < p->cds_start1 || p->cds_stop_end1 > p->rna_length ||
        !w->dna || !w->reference_rna || !w->alternate_rna ||
        !w->reference_peptide || !w->alternate_peptide ||
        !duckvep_codon_table_supported((duckvep_codon_table_t)q->codon_table))
        return DUCKVEP_PROFILE_INVALID;
    if (g->transcript.start1 < g->reference_start1 ||
        g->transcript.end1 < g->transcript.start1 ||
        (uint64_t)g->transcript.end1 >= (uint64_t)g->reference_start1 + g->dna_length ||
        g->exons[0].start1 != g->transcript.start1 ||
        g->exons[g->exon_count - 1u].end1 != g->transcript.end1)
        return DUCKVEP_PROFILE_PARTIAL_MODEL;
    if (w->dna_capacity < g->dna_length || w->rna_capacity < p->rna_length ||
        w->peptide_capacity < (p->rna_length - p->cds_start1 + 1u) / 3u + 1u ||
        w->source_effect_capacity < q->source_count)
        return DUCKVEP_PROFILE_SHORT_BUFFER;
    if (q->source_count && !w->source_effects) return DUCKVEP_PROFILE_INVALID;
    if (!span_disjoint(out, sizeof *out, g, sizeof *g) ||
        !span_disjoint(out, sizeof *out, p, sizeof *p) ||
        !span_disjoint(out, sizeof *out, g->reference_dna, g->dna_length) ||
        !span_disjoint(out, sizeof *out, p->reference_rna, p->rna_length) ||
        !span_disjoint(out, sizeof *out, g->exons,
                       g->exon_count * sizeof *g->exons) ||
        !span_disjoint(out, sizeof *out, p->blocks,
                       p->block_count * sizeof *p->blocks) ||
        !span_disjoint(out, sizeof *out, q->sources,
                       q->source_count * sizeof *q->sources) ||
        !span_disjoint(out, sizeof *out, q->source_ids,
                       q->source_count * sizeof *q->source_ids))
        return DUCKVEP_PROFILE_INVALID;
    if (w->source_effect_capacity > SIZE_MAX / sizeof *w->source_effects ||
        w->edit_capacity > SIZE_MAX / sizeof *w->edits)
        return DUCKVEP_PROFILE_INVALID;
    const void *buf[] = {w->dna, w->reference_rna, w->alternate_rna,
        w->reference_peptide, w->alternate_peptide, w->source_effects};
    size_t lengths[] = {w->dna_capacity, w->rna_capacity, w->rna_capacity,
        w->peptide_capacity, w->peptide_capacity,
        w->source_effect_capacity * sizeof *w->source_effects};
    for (size_t i = 0; i < 6; i++) {
        if (!span_disjoint(buf[i], lengths[i], g->reference_dna, g->dna_length) ||
            !span_disjoint(buf[i], lengths[i], p->reference_rna, p->rna_length) ||
            !span_disjoint(buf[i], lengths[i], out, sizeof *out) ||
            !span_disjoint(buf[i], lengths[i], q->sources,
                q->source_count * sizeof *q->sources) ||
            !span_disjoint(buf[i], lengths[i], q->source_ids,
                q->source_count * sizeof *q->source_ids) ||
            !span_disjoint(buf[i], lengths[i], g, sizeof *g) ||
            !span_disjoint(buf[i], lengths[i], p, sizeof *p) ||
            !span_disjoint(buf[i], lengths[i], g->exons,
                g->exon_count * sizeof *g->exons) ||
            !span_disjoint(buf[i], lengths[i], p->blocks,
                p->block_count * sizeof *p->blocks)) return DUCKVEP_PROFILE_INVALID;
        for (size_t j = 0; j < i; j++)
            if (!span_disjoint(buf[i], lengths[i], buf[j], lengths[j]))
                return DUCKVEP_PROFILE_INVALID;
    }
    for (size_t i = 0; i < g->dna_length; i++)
        if (!complement(g->reference_dna[i])) return DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE;
    /* Complete exon topology is required even when a named product skips one. */
    for (size_t i = 0; i < g->exon_count; i++) {
        duckvep_profile_interval_t e = g->exons[i];
        if (e.start1 < g->transcript.start1 || e.end1 > g->transcript.end1 ||
            e.start1 > e.end1 ||
            (i && g->exons[i-1].end1 >= e.start1 - 1u))
            return DUCKVEP_PROFILE_PARTIAL_MODEL;
    }
    duckvep_haplotype_edit_t small_edits[256];
    duckvep_haplotype_edit_t *edits = q->source_count <= 256u ? small_edits : w->edits;
    if (q->source_count > 256u) {
        if (!edits || w->edit_capacity < q->source_count) return DUCKVEP_PROFILE_SHORT_BUFFER;
        size_t size = w->edit_capacity * sizeof *edits;
        if (!span_disjoint(edits, size, out, sizeof *out) ||
            !span_disjoint(edits, size, g->reference_dna, g->dna_length) ||
            !span_disjoint(edits, size, p->reference_rna, p->rna_length) ||
            !span_disjoint(edits, size, q->sources, q->source_count * sizeof *q->sources) ||
            !span_disjoint(edits, size, q->source_ids, q->source_count * sizeof *q->source_ids) ||
            !span_disjoint(edits, size, g->exons, g->exon_count * sizeof *g->exons) ||
            !span_disjoint(edits, size, p->blocks, p->block_count * sizeof *p->blocks))
            return DUCKVEP_PROFILE_INVALID;
        for (size_t i = 0; i < 6; i++)
            if (!span_disjoint(edits, size, buf[i], lengths[i]))
                return DUCKVEP_PROFILE_INVALID;
    }
    for (size_t i = 0; i < q->source_count; i++) {
        duckvep_haplotype_edit_t e = q->sources[i];
        if (!e.ref_len || e.ref_len != e.alt_len ||
            e.variant_strand != 1 || e.cds_start < g->reference_start1 ||
            (uint64_t)e.cds_start + e.ref_len >
                (uint64_t)g->reference_start1 + g->dna_length ||
            (i && (q->sources[i-1].cds_start <= e.cds_start ||
                (uint64_t)e.cds_start + e.ref_len > q->sources[i-1].cds_start)))
            return DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION;
        for (size_t j = 0; j < i; j++)
            if (q->source_ids[j] == q->source_ids[i])
                return DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION;
        if (!e.ref || !e.alt ||
            !span_disjoint(e.ref, e.ref_len, out, sizeof *out) ||
            !span_disjoint(e.alt, e.alt_len, out, sizeof *out))
            return DUCKVEP_PROFILE_INVALID;
        for (size_t j = 0; j < 6; j++)
            if (!span_disjoint(e.ref, e.ref_len, buf[j], lengths[j]) ||
                !span_disjoint(e.alt, e.alt_len, buf[j], lengths[j]))
                return DUCKVEP_PROFILE_INVALID;
        e.cds_start -= g->reference_start1 - 1u;
        edits[i] = e;
    }
    /* Check the original genomic REF for every source using the mutation
     * owner's independent reference-relative checks, not a supplied final DNA. */
    size_t dna_len = 0;
    duckvep_haplotype_result_t applied;
    duckvep_haplotype_status_t hs = duckvep_haplotype_apply_cds_edits(
        g->reference_dna, g->dna_length, edits, q->source_count, 1,
        w->dna, w->dna_capacity, &dna_len, &applied);
    if (hs == DUCKVEP_HAPLOTYPE_REF_MISMATCH) return DUCKVEP_PROFILE_REFERENCE_MISMATCH;
    if (hs != DUCKVEP_HAPLOTYPE_OK || dna_len != g->dna_length)
        return DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION;
    size_t offset = 0;
    uint32_t prev_lo = 0, prev_hi = 0, junction = 0;
    for (size_t i = 0; i < p->block_count; i++) {
        duckvep_profile_rna_block_v1_t b = p->blocks[i];
        uint32_t lo = b.genomic_source.start1, hi = b.genomic_source.end1;
        if (lo < g->transcript.start1 || hi > g->transcript.end1 || lo > hi ||
            b.rna_start1 != offset + 1u ||
            (size_t)(hi - lo + 1u) > p->rna_length - offset ||
            (i && (g->strand == 1 ? lo <= prev_hi : hi >= prev_lo)))
            return DUCKVEP_PROFILE_PARTIAL_MODEL;
        /* A block is a whole model exon or a whole interval strictly between
         * two adjacent model exons. No invented genomic-to-RNA interpolation. */
        int mapped = 0;
        for (size_t j = 0; j < g->exon_count; j++) {
            if (lo == g->exons[j].start1 && hi == g->exons[j].end1) mapped = 1;
            if (j && lo == g->exons[j-1].end1 + 1u &&
                hi == g->exons[j].start1 - 1u) mapped = 1;
        }
        if (!mapped) return DUCKVEP_PROFILE_PARTIAL_MODEL;
        /* Adjacent mapped blocks form one exon; an EJC requires a gap. */
        if (i && (g->strand == 1 ? lo > prev_hi + 1u : prev_lo > hi + 1u))
            junction = (uint32_t)offset;
        for (uint32_t k = 0; k <= hi - lo; k++) {
            uint32_t pos = g->strand == 1 ? lo + k : hi - k;
            size_t index = pos - g->reference_start1;
            uint8_t r = g->strand == 1 ? g->reference_dna[index] : complement(g->reference_dna[index]);
            uint8_t a = g->strand == 1 ? w->dna[index] : complement(w->dna[index]);
            if (!a) return DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE;
            if (r != p->reference_rna[offset + k]) return DUCKVEP_PROFILE_REFERENCE_MISMATCH;
            w->reference_rna[offset+k] = r;
            w->alternate_rna[offset+k] = a;
        }
        offset += hi - lo + 1u;
        prev_lo = lo; prev_hi = hi;
    }
    if (offset != p->rna_length) return DUCKVEP_PROFILE_PARTIAL_MODEL;
    uint32_t start = p->annotated_start1;
    if ((size_t)start + 2u > p->rna_length ||
        !duckvep_codon_is_start(w->reference_rna + start - 1u,
                                (duckvep_codon_table_t)q->codon_table))
        return DUCKVEP_PROFILE_PARTIAL_MODEL;
    int lost = !duckvep_codon_is_start(w->alternate_rna + start - 1u,
                                       (duckvep_codon_table_t)q->codon_table);
    uint32_t alt_start = start;
    if (p->alternative_start1 || p->trigger_id || p->counterfactual_id) {
        if (!p->alternative_start1 || !p->trigger_id || !p->counterfactual_id ||
            strcmp(p->trigger_id, "annotated-start-lost") ||
            strcmp(p->counterfactual_id, "reference-annotated-start-retained") ||
            p->alternative_start1 <= start ||
            (size_t)p->alternative_start1 + 2u > p->rna_length ||
            !duckvep_codon_is_start(w->alternate_rna + p->alternative_start1 - 1u,
                                    (duckvep_codon_table_t)q->codon_table))
            return DUCKVEP_PROFILE_UNSUPPORTED_PRODUCT;
        if (lost) alt_start = p->alternative_start1;
    } else if (lost) {
        /* No peptide is inferred from an uninitiated alternate ORF. */
        return DUCKVEP_PROFILE_UNSUPPORTED_PRODUCT;
    }
    duckvep_translation_t rt, at;
    if (duckvep_translate_cds(w->reference_rna + start - 1u,
            p->rna_length - start + 1u, (duckvep_codon_table_t)q->codon_table,
            w->reference_peptide, w->peptide_capacity, &rt) != DUCKVEP_TRANSLATION_OK ||
        duckvep_translate_cds(w->alternate_rna + alt_start - 1u,
            p->rna_length - alt_start + 1u, (duckvep_codon_table_t)q->codon_table,
            w->alternate_peptide, w->peptide_capacity, &at) != DUCKVEP_TRANSLATION_OK)
        return DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE;
    w->reference_peptide[0] = 'M';
    w->alternate_peptide[0] = 'M';
    if (rt.first_stop_position1) rt.length = rt.first_stop_position1;
    if (at.first_stop_position1) at.length = at.first_stop_position1;
    w->reference_peptide[rt.length] = 0;
    w->alternate_peptide[at.length] = 0;
    for (size_t i = 0; i < q->source_count; i++) {
        const duckvep_haplotype_edit_t *e = &q->sources[i];
        duckvep_profile_source_effect_t effect = {
            .source_id=q->source_ids[i], .genomic_start1=e->cds_start,
            .genomic_end1=e->cds_start + e->ref_len - 1u};
        for (uint32_t k = 0; k < e->ref_len; k++) {
            uint32_t pos = e->cds_start + k;
            if (g->reference_dna[pos - g->reference_start1] ==
                w->dna[pos - g->reference_start1]) continue;
            effect.changed_dna = 1;
            for (size_t j = 0; j < p->block_count; j++) {
                duckvep_profile_rna_block_v1_t b = p->blocks[j];
                if (pos < b.genomic_source.start1 || pos > b.genomic_source.end1) continue;
                uint32_t r = b.rna_start1 + (g->strand == 1
                    ? pos - b.genomic_source.start1 : b.genomic_source.end1 - pos);
                if (!effect.rna_first1 || r < effect.rna_first1) effect.rna_first1 = r;
                if (r > effect.rna_last1) effect.rna_last1 = r;
                effect.changed_rna_bases++;
            }
        }
        w->source_effects[i] = effect;
    }
    *out = (duckvep_profile_execution_result_t){
        .identity=p->identity, .evidence_id=p->evidence_id, .assumption_id=p->assumption_id,
        .rna_length=p->rna_length, .reference_peptide_length=rt.length,
        .alternate_peptide_length=at.length, .reference_stop1=rt.first_stop_position1,
        .alternate_stop1=at.first_stop_position1, .reference_start1=start,
        .alternate_start1=alt_start, .last_junction1=junction,
        .annotated_start_lost=(uint8_t)lost,
        .alternative_start_used=(uint8_t)(alt_start != start),
        .verified_sources=q->source_count};
    return DUCKVEP_PROFILE_OK;
}

duckvep_profile_status_t duckvep_product_execute_position_v1(
    const duckvep_profile_execution_t *input,
    duckvep_profile_workspace_t *workspace,
    duckvep_profile_execution_result_t *result,
    duckvep_profile_v1_fact_t *facts, size_t capacity, size_t *count) {
    if (!count) return DUCKVEP_PROFILE_INVALID;
    *count = 0;
    if (!result || (capacity && !facts) ||
        capacity > SIZE_MAX / sizeof *facts ||
        (input && input->source_count > SIZE_MAX / sizeof *input->sources) ||
        (input && input->source_count > SIZE_MAX / sizeof *input->source_ids) ||
        (workspace && workspace->source_effect_capacity >
            SIZE_MAX / sizeof *workspace->source_effects))
        return DUCKVEP_PROFILE_INVALID;
    if (capacity && input && input->geometry && input->programme && workspace) {
        const duckvep_profile_v1_input_t *g = input->geometry;
        const duckvep_profile_programme_v1_t *p = input->programme;
        size_t bytes = capacity * sizeof *facts;
        if (!span_disjoint(facts, bytes, result, sizeof *result) ||
            !span_disjoint(facts, bytes, g->reference_dna, g->dna_length) ||
            !span_disjoint(facts, bytes, p->reference_rna, p->rna_length) ||
            !span_disjoint(facts, bytes, input->sources,
                input->source_count * sizeof *input->sources) ||
            !span_disjoint(facts, bytes, input->source_ids,
                input->source_count * sizeof *input->source_ids) ||
            !span_disjoint(facts, bytes, workspace->dna, workspace->dna_capacity) ||
            !span_disjoint(facts, bytes, workspace->reference_rna, workspace->rna_capacity) ||
            !span_disjoint(facts, bytes, workspace->alternate_rna, workspace->rna_capacity) ||
            !span_disjoint(facts, bytes, workspace->reference_peptide, workspace->peptide_capacity) ||
            !span_disjoint(facts, bytes, workspace->alternate_peptide, workspace->peptide_capacity) ||
            !span_disjoint(facts, bytes, workspace->source_effects,
                workspace->source_effect_capacity * sizeof *workspace->source_effects))
            return DUCKVEP_PROFILE_INVALID;
    }
    duckvep_profile_execution_result_t pending;
    duckvep_profile_status_t status = duckvep_product_execute_v1(input, workspace, &pending);
    if (status != DUCKVEP_PROFILE_OK) return status;
    duckvep_profile_v1_input_t view = *input->geometry;
    view.reference_witness = view.reference_dna;
    view.witness_length = view.dna_length;
    view.final_dna = workspace->dna;
    view.final_length = view.dna_length;
    view.joint_status = DUCKVEP_PROFILE_OK;
    status = duckvep_product_position_v1(&view, facts, capacity, count);
    if (status == DUCKVEP_PROFILE_OK) *result = pending;
    return status;
}
