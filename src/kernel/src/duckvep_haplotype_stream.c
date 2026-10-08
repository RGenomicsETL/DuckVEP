#include "duckvep_haplotype_stream.h"
#include "duckvep_classify.h"
#include "duckvep_effect.h"
#include "duckvep_so.h"
#include "duckvep_transcript_edit.h"

#include <string.h>

int duckvep_haplotype_record_plan_init(duckvep_haplotype_record_plan_t *p,
    const duckvep_transcript_model_t *m) {
    if (!p) return 0;
    memset(p, 0, sizeof(*p));
    if (!m || (m->transcript_count && (!m->chrom_id || !m->start1 || !m->end1))) return 0;
    for (size_t i = 0u; i < m->transcript_count; i++)
        if (!m->start1[i] || m->start1[i] > m->end1[i] ||
            (i && (m->chrom_id[i] < m->chrom_id[i - 1u] ||
             (m->chrom_id[i] == m->chrom_id[i - 1u] && m->start1[i] < m->start1[i - 1u])))) return 0;
    p->model = m;
    return 1;
}

int duckvep_haplotype_record_plan_next(duckvep_haplotype_record_plan_t *p,
    uint16_t chrom, uint32_t start, uint32_t end, uint64_t *buffer, uint64_t *ordinal) {
    if (buffer) *buffer = 0u;
    if (ordinal) *ordinal = 0u;
    if (!p || !p->model || !buffer || !ordinal || !start || end < start ||
        (p->have_input && (chrom < p->chrom || (chrom == p->chrom && start < p->last_pos1)))) return 0;
    if (!p->have_input || chrom != p->chrom) p->end1 = 0u;
    p->have_input = 1u; p->chrom = chrom; p->last_pos1 = start;
    if (start > p->end1) {
        const duckvep_transcript_model_t *m = p->model;
        int overlaps = 0;
        while (p->transcript < m->transcript_count && m->chrom_id[p->transcript] < chrom)
            p->transcript++;
        while (p->transcript < m->transcript_count && m->chrom_id[p->transcript] == chrom &&
                m->start1[p->transcript] <= end) {
            uint32_t tx_end = m->end1[p->transcript++];
            if (tx_end < start) continue;
            overlaps = 1;
            if (tx_end > end) end = tx_end;
        }
        p->end1 = 0u;
        if (!overlaps) return 1;
        if (p->buffer == UINT64_MAX) return 0;
        p->buffer++; p->ordinal = 0u; p->end1 = end;
    }
    if (p->ordinal == UINT64_MAX) return 0;
    *buffer = p->buffer; *ordinal = ++p->ordinal;
    return 1;
}

uint64_t duckvep_haplotype_record_order(uint64_t count, uint64_t ordinal) {
    if (!ordinal || ordinal > count) return 0u;
    uint64_t rank = 1u;
    while (count) {
        /* Sorted red-black insertion doubles root ordinal r at 5*r-2
         * records. Its left subtree is perfect with r-1 nodes; the right
         * subtree is another sorted-insertion tree. Divide before adding
         * to keep the threshold defined through UINT64_MAX records. */
        uint64_t threshold = count / 5u + (count % 5u + 2u) / 5u;
        uint64_t root = 1u;
        while (root <= threshold) root *= 2u;
        if (ordinal == root) return rank;
        if (ordinal > root) {
            ordinal -= root; count -= root; rank += root;
            continue;
        }
        rank++;
        for (root /= 2u; root; root /= 2u) {
            if (ordinal == root) return rank;
            if (ordinal > root) { ordinal -= root; rank += root; }
            else rank++;
        }
    }
    return 0u;
}

static duckvep_haplotype_stream_status_t fail(
    duckvep_haplotype_stream_t *s, duckvep_haplotype_stream_status_t status) {
    if (s) s->error = status;
    return status;
}

static duckvep_haplotype_stream_status_t carrier_fail(
    duckvep_haplotype_stream_t *s, duckvep_carriers_status_t status) {
    s->carrier_error = status;
    return fail(s, status == DUCKVEP_CARRIERS_INPUT_ORDER
        ? DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER : DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR);
}

/* Both arguments are at most capacity; avoid overflowing an address-size sum. */
static size_t ring_add(size_t begin, size_t count, size_t capacity) {
    return count >= capacity - begin ? count - (capacity - begin) : begin + count;
}

static int valid_array(const void *p, size_t n, size_t width) {
    return p && n && n <= SIZE_MAX / width;
}

static int ranges_overlap(const void *left, size_t left_length,
    const void *right, size_t right_length)
{
    uintptr_t a = (uintptr_t)left, b = (uintptr_t)right;

    if (!left || !right || !left_length || !right_length)
        return 0;
    return a <= b ? b - a < left_length : a - b < right_length;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_init(
    duckvep_haplotype_stream_t *s, const duckvep_transcript_model_t *tx,
    const duckvep_exon_model_t *exons, const duckvep_sequence_pool_t *seq,
    const duckvep_haplotype_stream_buffers_t *b) {
    if (!s) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    memset(s, 0, sizeof(*s));
    if (!tx || !exons || !seq || !b || !tx->start1 || !tx->strand ||
        seq->transcript_count != tx->transcript_count ||
        !seq->cds_offset || !seq->cds_length || (!seq->cds_bytes && seq->cds_bytes_len) ||
        !valid_array(b->events, b->event_capacity, sizeof(*b->events)) ||
        !valid_array(b->projections, b->projection_capacity, sizeof(*b->projections)) ||
        !valid_array(b->alleles, b->allele_capacity, 1u) ||
        !valid_array(b->leaf_events, b->leaf_capacity, sizeof(*b->leaf_events)) ||
        !valid_array(b->contributors, b->leaf_capacity, sizeof(*b->contributors)) ||
        !valid_array(b->edits, b->edit_capacity, sizeof(*b->edits)) ||
        !valid_array(b->edit_event_ids, b->edit_capacity, sizeof(*b->edit_event_ids)) ||
        !valid_array(b->blocks, b->edit_capacity, sizeof(*b->blocks)) ||
        !valid_array(b->cds, b->cds_capacity, 1u) ||
        !valid_array(b->protein, b->protein_capacity, 1u) ||
        !valid_array(b->reference_protein, b->reference_protein_capacity, 1u) ||
        !valid_array(b->reference_coding_protein, b->reference_protein_capacity, 1u) ||
        (!!b->noncoding != !!b->noncoding_capacity))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    if (b->prediction_protein &&
        (ranges_overlap(b->prediction_protein, b->prediction_protein_capacity,
             b->cds, b->cds_capacity) ||
         ranges_overlap(b->prediction_protein, b->prediction_protein_capacity,
             b->protein, b->protein_capacity) ||
         ranges_overlap(b->prediction_protein, b->prediction_protein_capacity,
             b->reference_protein, b->reference_protein_capacity) ||
         ranges_overlap(b->prediction_protein, b->prediction_protein_capacity,
             b->reference_coding_protein, b->reference_protein_capacity)))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    duckvep_carriers_status_t status = duckvep_carriers_init(&s->carriers, tx, &b->carriers);
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    s->exons = exons;
    s->sequences = seq;
    s->buffers = *b;
    s->initialized = 1u;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

/* Called only after the carrier watermark has drained every earlier transcript.
 * An old event can pin younger events behind it, but nothing beyond the active
 * genomic window is retained and all three rings are reclaimed together. */
static void reclaim(duckvep_haplotype_stream_t *s, uint16_t chrom, uint32_t pos1, int eof) {
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    while (s->event_count) {
        const duckvep_haplotype_stored_event_t *event = &b->events[s->event_begin];
        if (!eof && (event->source.chrom_id > chrom ||
            (event->source.chrom_id == chrom && event->last_end1 >= pos1))) break;
        s->allele_begin = ring_add(s->allele_begin, event->allele_consumed, b->allele_capacity);
        s->allele_count -= event->allele_consumed;
        s->projection_begin = (uint32_t)ring_add(s->projection_begin,
            event->projection_count, b->projection_capacity);
        s->projection_count -= event->projection_count;
        s->event_begin = (uint32_t)ring_add(s->event_begin, 1u, b->event_capacity);
        s->event_count--;
    }
    if (!s->event_count) {
        s->event_begin = s->projection_begin = 0u;
        s->allele_begin = 0u;
    }
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_begin(
    duckvep_haplotype_stream_t *s, const duckvep_haplotype_source_t *source) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!source || !source->pos1 || !source->ref || !source->alt ||
        !source->ref_len || source->source_record > 1u ||
        (!source->source_record && (!source->alt_len || source->allele_index || source->replay_order)) ||
        (source->source_record && ((source->allele_index == UINT32_MAX) != !source->alt_len)) ||
        (source->source_record && source->allele_index != UINT32_MAX && source->allele_index > INT32_MAX) ||
        s->serial == UINT64_MAX ||
        (uint32_t)source->ref_len - 1u > UINT32_MAX - source->pos1)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    if (s->have_input && (source->chrom_id < s->last_chrom ||
        (source->chrom_id == s->last_chrom && (source->pos1 < s->last_pos1 ||
         (source->pos1 == s->last_pos1 && (source->event_id < s->last_event_id ||
          (source->event_id == s->last_event_id && (!source->source_record ||
           source->allele_index <= s->last_allele_index))))))))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER);
    if (s->have_input && source->chrom_id == s->last_chrom && source->pos1 == s->last_pos1 &&
        source->event_id == s->last_event_id) {
        const duckvep_haplotype_source_t *previous = &s->buffers.events[s->current_event].source;
        if (!s->have_current || !previous->source_record || previous->ref_len != source->ref_len ||
            memcmp(previous->ref, source->ref, source->ref_len) || previous->replay_order != source->replay_order)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    }
    if (source->source_record && !source->allele_index &&
        (source->ref_len != source->alt_len || memcmp(source->ref, source->alt, source->ref_len)))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    duckvep_event_t prepared = {0};
    int valid = source->source_record
        ? duckvep_event_prepare_replacement(source->pos1, source->ref, source->ref_len,
            source->alt, source->alt_len, &prepared)
        : duckvep_event_prepare_small(source->pos1, source->ref, source->ref_len,
            source->alt, source->alt_len, &prepared);
    if (!valid)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    prepared.chrom_id = source->chrom_id;
    uint32_t completed;
    duckvep_carriers_status_t status = duckvep_carriers_advance(&s->carriers,
        source->chrom_id, source->pos1, s->serial + 1u, &completed);
    if (status == DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
        s->closing = 1u;
        return DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY;
    }
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    reclaim(s, source->chrom_id, source->pos1, 0);
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    if (s->event_count == b->event_capacity)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL);
    size_t bytes = (size_t)source->ref_len + source->alt_len;
    size_t at = ring_add(s->allele_begin, s->allele_count, b->allele_capacity);
    size_t padding = bytes > b->allele_capacity - at ? b->allele_capacity - at : 0u;
    size_t available = b->allele_capacity - s->allele_count;
    if (padding > available || bytes > available - padding)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL);
    if (padding) at = 0u;
    uint32_t event_at = (uint32_t)ring_add(s->event_begin, s->event_count, b->event_capacity);
    duckvep_haplotype_stored_event_t *stored = &b->events[event_at];
    *stored = (duckvep_haplotype_stored_event_t){0};
    stored->source = *source;
    stored->prepared = prepared;
    stored->source.ref = b->alleles + at;
    stored->source.alt = b->alleles + at + source->ref_len;
    memcpy(b->alleles + at, source->ref, source->ref_len);
    memcpy(b->alleles + at + source->ref_len, source->alt, source->alt_len);
    stored->serial = s->serial + 1u;
    stored->allele_consumed = padding + bytes;
    stored->projection_begin = (uint32_t)ring_add(s->projection_begin,
        s->projection_count, b->projection_capacity);
    stored->last_end1 = source->pos1;
    s->current_event = event_at;
    s->event_count++;
    s->allele_count += stored->allele_consumed;
    if (s->event_count > s->peak_events) s->peak_events = s->event_count;
    if (s->allele_count > s->peak_alleles) s->peak_alleles = s->allele_count;
    s->have_current = s->have_input = 1u;
    s->serial++;
    s->last_event_id = source->event_id;
    s->last_allele_index = source->allele_index;
    s->last_pos1 = source->pos1;
    s->last_chrom = source->chrom_id;
    s->input_events++;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_project(
    duckvep_haplotype_stream_t *s, uint32_t tx) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!s->have_current || s->closing || s->carriers.pending || s->carriers.finished)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    const duckvep_transcript_model_t *model = s->carriers.model;
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    duckvep_haplotype_stored_event_t *stored = &b->events[s->current_event];
    const duckvep_event_t *prepared = &stored->prepared;
    if (tx >= model->transcript_count || model->chrom_id[tx] != stored->source.chrom_id ||
        model->end1[tx] < stored->source.pos1 ||
        (model->start1[tx] > prepared->raw_end1 &&
         model->start1[tx] > duckvep_event_feature_max1(prepared)))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    if (stored->projection_count) {
        size_t previous = ring_add(stored->projection_begin, stored->projection_count - 1u,
            b->projection_capacity);
        if (tx <= b->projections[previous].transcript_index)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER);
    }
    if (s->projection_count == b->projection_capacity)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL);
    if (s->projected_events == UINT64_MAX)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    size_t at = ring_add(stored->projection_begin, stored->projection_count, b->projection_capacity);
    duckvep_haplotype_projection_t *p = &b->projections[at];
    duckvep_prepared_cds_allele_t allele = {prepared,
        stored->source.ref + prepared->ref_diff_offset,
        stored->source.alt + prepared->alt_diff_offset,
        stored->source.ref + prepared->anchor_ref_offset,
        prepared->ref_diff_length, prepared->alt_diff_length, 1};
    p->transcript_index = tx;
    memset(&p->edit, 0, sizeof(p->edit));
    p->status = stored->source.source_record
        ? duckvep_compat_vep_source_cds_edit_build(model, s->exons, s->sequences,
            tx, model->strand[tx], &allele, &p->edit)
        : duckvep_cds_edit_build_prepared_allele(model, s->exons, s->sequences,
            tx, model->strand[tx], &allele, UINT32_MAX, &p->edit);
    p->cds_unaffected = 0u;
    p->source_selected = 1u;
    p->selection_set = 0u;
    p->source_exonic = 1u;
    if (p->status == DUCKVEP_CDS_EDIT_OUT_OF_CDS && s->sequences->cds_length[tx] &&
        model->cds_start1 && model->cds_end1 && model->cds_start1[tx] &&
        model->exon_offset && model->exon_count && model->exon_count[tx] &&
        s->exons->start1 && s->exons->end1 &&
        model->exon_offset[tx] <= s->exons->exon_count &&
        model->exon_count[tx] <= s->exons->exon_count - model->exon_offset[tx]) {
        /* OUT_OF_CDS also covers failed CDS-slice bounds and coding/noncoding
         * crossings. Only the shared topology classifier may prove absence of
         * a coding overlap; an insertion examines both reference flanks. */
        duckvep_region_state_t region = duckvep_region_classify_span(model, s->exons, tx,
            prepared->interbase ? prepared->insertion_boundary0 : prepared->start1,
            prepared->interbase ? duckvep_event_right_flank1(prepared) : prepared->end1,
            0u, 0u);
        p->source_exonic = region.overlaps_exon;
        /* Haplosaurus admits whole source spans through exons. Intronic
         * context still has provenance but cannot alter its literal CDS,
         * including the short gaps classified as frameshift introns by SO. */
        p->cds_unaffected = !region.overlaps_cds ||
            (stored->source.source_record && !region.overlaps_exon);
    }
    stored->projection_count++;
    if (model->end1[tx] > stored->last_end1) stored->last_end1 = model->end1[tx];
    s->projection_count++;
    s->projected_events++;
    if (s->projection_count > s->peak_projections) s->peak_projections = s->projection_count;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_finish(duckvep_haplotype_stream_t *s) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    uint32_t completed;
    duckvep_carriers_status_t status = duckvep_carriers_finish(&s->carriers, &completed);
    if (status == DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
        s->closing = 1u;
        return DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY;
    }
    if (status != DUCKVEP_CARRIERS_DONE) return carrier_fail(s, status);
    reclaim(s, 0u, 0u, 1);
    s->have_current = 0u;
    return DUCKVEP_HAPLOTYPE_STREAM_DONE;
}

static const duckvep_haplotype_stored_event_t *find_event(
    const duckvep_haplotype_stream_t *s, uint64_t serial) {
    size_t lo = 0u, hi = s->event_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        size_t at = ring_add(s->event_begin, mid, s->buffers.event_capacity);
        const duckvep_haplotype_stored_event_t *e = &s->buffers.events[at];
        if (e->serial == serial) return e;
        if (e->serial < serial) lo = mid + 1u;
        else hi = mid;
    }
    return NULL;
}

static duckvep_haplotype_projection_t *find_projection(
    const duckvep_haplotype_stream_t *s, const duckvep_haplotype_stored_event_t *e, uint32_t tx) {
    size_t lo = 0u, hi = e->projection_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        size_t at = ring_add(e->projection_begin, mid, s->buffers.projection_capacity);
        duckvep_haplotype_projection_t *p = &s->buffers.projections[at];
        if (p->transcript_index == tx) return p;
        if (p->transcript_index < tx) lo = mid + 1u;
        else hi = mid;
    }
    return NULL;
}

static int same_phase_set(duckvep_haplotype_phase_set_t a, duckvep_haplotype_phase_set_t b) {
    return a.present == b.present && (!a.present || a.value == b.value);
}

static duckvep_haplotype_stream_status_t push_call_lane(
    duckvep_haplotype_stream_t *s, uint32_t tx, const duckvep_haplotype_call_t *call,
    duckvep_haplotype_phase_set_t set, uint16_t lane, uint8_t evidence, uint8_t split) {
    duckvep_carrier_key_t key = {call->sample_index, set.value, lane, call->ploidy, set.present, split};
    duckvep_carriers_status_t status = duckvep_carriers_push(&s->carriers, tx, &key, evidence);
    return status == DUCKVEP_CARRIERS_OK ? DUCKVEP_HAPLOTYPE_STREAM_OK : carrier_fail(s, status);
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push_call(
    duckvep_haplotype_stream_t *s, uint32_t tx, const duckvep_haplotype_call_t *call,
    const duckvep_haplotype_phase_set_t *sets, size_t set_count) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!s->have_current || s->closing || s->carriers.pending || s->carriers.finished ||
        !call || !call->alleles || !call->ploidy || !call->alt_index ||
        call->alt_index > INT32_MAX || call->phase_set.present > 1u || call->hypothetical > 1u ||
        s->buffers.events[s->current_event].source.source_record ||
        (call->policy != DUCKVEP_PHASE_STRICT && call->policy != DUCKVEP_PHASE_VEP_COMPAT) ||
        (s->have_phase_policy && s->phase_policy != call->policy) ||
        (set_count && !sets) || set_count > SIZE_MAX / sizeof(*sets) ||
        !find_projection(s, &s->buffers.events[s->current_event], tx))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    const duckvep_haplotype_phase_set_t absent = {0, 0u};
    if (!set_count) { sets = &absent; set_count = 1u; }
    size_t declared_set = set_count;
    for (size_t i = 0u; i < set_count; i++) {
        if (sets[i].present > 1u || (i && (!sets[i].present ||
            (sets[i - 1u].present && sets[i].value <= sets[i - 1u].value))))
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
        if (same_phase_set(sets[i], call->phase_set)) declared_set = i;
    }
    if (call->policy == DUCKVEP_PHASE_VEP_COMPAT &&
        (set_count != 1u || sets[0].present))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);

    duckvep_phase_summary_t summary = {0};
    uint8_t missing = 0u, pool_missing = 0u, pool_alt = 0u;
    for (uint32_t slot = 0u; slot < call->ploidy; slot++) {
        int32_t allele = call->alleles[slot];
        uint8_t phase = call->phase_before ? call->phase_before[slot] : 0u;
        if (duckvep_phase_observe(&summary, allele, phase) != DUCKVEP_PHASE_OK)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
        if (allele < 0) missing = 1u;
        if (!phase) {
            if (allele < 0) pool_missing = 1u;
            if (allele == (int32_t)call->alt_index) pool_alt = 1u;
        }
    }
    int broadcast = summary.ploidy == 1u || summary.homozygous ||
        summary.unphased_count == summary.ploidy;
    if (call->policy == DUCKVEP_PHASE_STRICT && !broadcast && declared_set == set_count)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    s->phase_policy = call->policy;
    s->have_phase_policy = 1u;
    /* Homozygous/unphased calls are broadcast over the domain, so every carrier
     * of a sample with several heterozygous phase sets is marked, not only the
     * carriers that happen to use a second set. */
    const uint8_t split = call->policy == DUCKVEP_PHASE_STRICT && set_count > 1u;

    uint16_t called_before = 0u;
    for (uint32_t slot = 0u; slot < call->ploidy; slot++) {
        int32_t allele = call->alleles[slot];
        uint8_t phase = call->phase_before ? call->phase_before[slot] : 0u;
        duckvep_phase_assignment_t assignment;
        if (duckvep_phase_assign(&summary, (uint16_t)(slot + 1u), called_before,
            allele, phase, call->policy, &assignment) != DUCKVEP_PHASE_OK)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
        if (allele >= 0) called_before++;
        uint8_t evidence = 0u;
        if (call->policy == DUCKVEP_PHASE_VEP_COMPAT) {
            if (allele < 0) continue; /* Missing slots do not consume compacted lanes. */
            if (missing) evidence = DUCKVEP_CARRIER_MISSING;
            if (allele == (int32_t)call->alt_index) evidence |= DUCKVEP_CARRIER_CALLED;
        } else if (assignment.scope == DUCKVEP_PHASE_UNRESOLVED) {
            if (pool_alt || pool_missing) evidence = DUCKVEP_CARRIER_UNPHASED;
            if (pool_missing) evidence |= DUCKVEP_CARRIER_MISSING;
            assignment.lane = (uint16_t)(slot + 1u);
        } else if (allele < 0) {
            evidence = DUCKVEP_CARRIER_MISSING;
        } else if (allele == (int32_t)call->alt_index) {
            evidence = DUCKVEP_CARRIER_CALLED;
        }
        if (!evidence) continue;
        if (call->hypothetical) evidence |= DUCKVEP_CARRIER_HYPOTHETICAL;
        size_t first = broadcast || call->policy == DUCKVEP_PHASE_VEP_COMPAT ? 0u : declared_set;
        size_t end = broadcast ? set_count : first + 1u;
        for (size_t i = first; i < end; i++) {
            duckvep_haplotype_stream_status_t status = push_call_lane(s, tx, call, sets[i],
                assignment.lane, evidence, split);
            if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        }
    }
    /* Compaction loses the original missing slot positions, not their evidence.
     * Include the remaining lanes so an all-missing GT is never implicit REF. */
    if (call->policy == DUCKVEP_PHASE_VEP_COMPAT && missing) {
        for (uint32_t lane = (uint32_t)called_before + 1u; lane <= call->ploidy; lane++) {
            duckvep_haplotype_stream_status_t status = push_call_lane(s, tx, call, absent,
                (uint16_t)lane, DUCKVEP_CARRIER_MISSING, 0u);
            if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        }
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push_raw_call(
    duckvep_haplotype_stream_t *s, uint32_t tx, uint32_t sample, const duckvep_raw_gt_t *call,
    uint8_t source_selected) {
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!s->have_current || s->closing || s->carriers.pending || s->carriers.finished ||
        !call || !call->source_ploidy || call->source_has_missing > 1u || source_selected > 1u ||
        call->disposition < DUCKVEP_RAW_GT_OMITTED_REFERENCE ||
        call->disposition > DUCKVEP_RAW_GT_RETAINED ||
        (s->have_phase_policy && s->phase_policy != DUCKVEP_PHASE_VEP_RAW))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    const duckvep_haplotype_stored_event_t *event = &s->buffers.events[s->current_event];
    duckvep_haplotype_projection_t *projection = find_projection(s, event, tx);
    if (!event->source.source_record || !projection ||
        (projection->selection_set && projection->source_selected != source_selected))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    projection->source_selected = source_selected;
    projection->selection_set = 1u;
    if (call->disposition == DUCKVEP_RAW_GT_RETAINED) {
        if (!call->parsed_slots || call->parsed_slots > (uint32_t)call->source_ploidy + 1u ||
            call->allele_index[0] > INT32_MAX ||
            (call->parsed_slots == 1u ? call->allele_index[1] != UINT32_MAX
                                     : call->allele_index[1] > INT32_MAX))
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    } else if (call->parsed_slots || call->allele_index[0] != UINT32_MAX ||
               call->allele_index[1] != UINT32_MAX ||
               call->source_has_missing != (call->disposition == DUCKVEP_RAW_GT_OMITTED_EMPTY)) {
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    }
    s->phase_policy = DUCKVEP_PHASE_VEP_RAW;
    s->have_phase_policy = 1u;
    uint32_t allele = event->source.allele_index;
    if (call->disposition == DUCKVEP_RAW_GT_OMITTED_REFERENCE) return DUCKVEP_HAPLOTYPE_STREAM_OK;
    for (uint16_t lane = 1u; lane <= 2u; lane++) {
        uint8_t evidence = 0u;
        if (call->disposition == DUCKVEP_RAW_GT_OMITTED_EMPTY) {
            /* No upstream edit, but an omitted missing call is not proof of REF.
             * Retain its conditional no-op observation on both occupied paths. */
            if (allele) continue;
            evidence = DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_CONDITIONAL;
        } else {
            if (call->allele_index[lane - 1u] != allele) continue;
            if (allele == UINT32_MAX) evidence = DUCKVEP_CARRIER_CONDITIONAL;
            else if (allele) evidence = DUCKVEP_CARRIER_CALLED;
            else evidence = DUCKVEP_CARRIER_REFERENCE_REPLAY;
            if (call->source_has_missing)
                evidence |= DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_CONDITIONAL;
        }
        duckvep_carrier_key_t key = {sample, 0, lane, 2u, 0u, 0u};
        duckvep_carriers_status_t status = duckvep_carriers_push(&s->carriers, tx, &key, evidence);
        if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

/* Descending CDS order with ascending source ordinal for equal starts. */
static int edit_precedes_min(uint32_t start, uint64_t id, uint32_t other, uint64_t other_id,
    const duckvep_haplotype_contributor_t *contributors) {
    if (contributors) {
        if (contributors[id].source.replay_order) id = contributors[id].source.replay_order;
        if (contributors[other_id].source.replay_order) other_id = contributors[other_id].source.replay_order;
    }
    return start < other || (start == other && id > other_id);
}

static void sift_edit_min(duckvep_haplotype_edit_t *edits, uint64_t *ids,
    size_t root, size_t count, const duckvep_haplotype_contributor_t *contributors) {
    duckvep_haplotype_edit_t value = edits[root];
    uint64_t id = ids[root];
    while (root < count / 2u) {
        size_t child = root * 2u + 1u;
        if (child + 1u < count && edit_precedes_min(edits[child + 1u].cds_start,
                ids[child + 1u], edits[child].cds_start, ids[child], contributors)) child++;
        if (!edit_precedes_min(edits[child].cds_start, ids[child], value.cds_start, id, contributors)) break;
        edits[root] = edits[child];
        ids[root] = ids[child];
        root = child;
    }
    edits[root] = value;
    ids[root] = id;
}

/* Genomic upload order is not CDS edit order: trimming retained REF can move
 * an earlier upload past the next edit, and reverse transcripts invert it.
 * A typed in-place heap keeps scratch constant; libc qsort may allocate. */
static void sort_edits_descending(duckvep_haplotype_edit_t *edits, uint64_t *ids, size_t count,
    const duckvep_haplotype_contributor_t *contributors) {
    for (size_t root = count / 2u; root > 0u; root--)
        sift_edit_min(edits, ids, root - 1u, count, contributors);
    for (size_t remaining = count; remaining > 1u; remaining--) {
        duckvep_haplotype_edit_t first = edits[0];
        uint64_t id = ids[0];
        edits[0] = edits[remaining - 1u];
        edits[remaining - 1u] = first;
        ids[0] = ids[remaining - 1u];
        ids[remaining - 1u] = id;
        sift_edit_min(edits, ids, 0u, remaining - 1u, contributors);
    }
}

static duckvep_haplotype_stream_status_t append_differing_edits(
    duckvep_haplotype_stream_t *s, const duckvep_haplotype_projection_t *p,
    uint64_t event_id, duckvep_haplotype_leaf_t *leaf) {
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    duckvep_edit_set_t edits;
    duckvep_cds_edit_status_t status = duckvep_projected_cds_edit_set_build(&p->edit,
        s->carriers.model->strand[leaf->carriers.transcript_index], b->edits + leaf->edit_count,
        b->edit_capacity - leaf->edit_count, &edits);
    if (status == DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL);
    if (status == DUCKVEP_CDS_EDIT_OK) {
        for (size_t j = 0u; j < edits.count; j++)
            b->edit_event_ids[leaf->edit_count + j] = event_id;
        leaf->edit_count += edits.count;
    } else if (leaf->projection_status == DUCKVEP_CDS_EDIT_OK) leaf->projection_status = status;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

/* One preparation per closing transcript supplies curated replay/difference
 * reference and uncurated coding operands from the same translation pass. */
static duckvep_haplotype_stream_status_t prepare_reference_protein(
    duckvep_haplotype_stream_t *s, uint32_t tx) {
    if (s->have_reference_protein && s->reference_transcript == tx)
        return DUCKVEP_HAPLOTYPE_STREAM_OK;
    const duckvep_sequence_pool_t *seq = s->sequences;
    size_t length = seq->cds_length[tx];
    uint64_t offset = seq->cds_offset[tx];
    size_t begin = seq->peptide_edit_offset ? seq->peptide_edit_offset[tx] : 0u;
    size_t end = seq->peptide_edit_offset ? seq->peptide_edit_offset[tx + 1u] : 0u;
    if (offset > seq->cds_bytes_len || length > seq->cds_bytes_len - offset ||
        (seq->peptide_edit_count && !seq->peptide_edit_offset) ||
        begin > end || end > seq->peptide_edit_count ||
        (end > begin && (!seq->peptide_edit_position1 || !seq->peptide_edit_alt)))
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    duckvep_codon_table_t table = seq->codon_table
        ? (duckvep_codon_table_t)seq->codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
    duckvep_haplotype_status_t status = duckvep_haplotype_reference_proteins(
        seq->cds_bytes + (size_t)offset, length, table,
        end > begin ? seq->peptide_edit_position1 + begin : NULL,
        end > begin ? seq->peptide_edit_alt + begin : NULL, end - begin,
        s->buffers.reference_protein, s->buffers.reference_coding_protein,
        s->buffers.reference_protein_capacity, &s->reference_protein_length,
        &s->reference_coding_translation);
    if (status == DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL);
    if (status != DUCKVEP_HAPLOTYPE_OK && status != DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    s->reference_protein_known = status == DUCKVEP_HAPLOTYPE_OK;
    s->reference_transcript = tx;
    s->have_reference_protein = 1u;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

/* ---- eligibility: status/reason only, no classifier ---- */

static int stop_codon(const uint8_t *c) {
    uint8_t a = c[0] & 0xDFu, b = c[1] & 0xDFu, d = c[2] & 0xDFu;
    return a == 'T' && ((b == 'A' && (d == 'A' || d == 'G')) || (b == 'G' && d == 'A'));
}

/* Start and stop codons under the transcript's genetic code. The standard code requires the
 * canonical ATG start; another code has no single canonical start, so any start codon of that code
 * counts (vertebrate mitochondrial transcripts begin with ATA, ATT or GTG as well as ATG). */
static int start_codon_of(const uint8_t *c, duckvep_codon_table_t table) {
    if (table == DUCKVEP_CODON_TABLE_STANDARD)
        return (c[0] & 0xDFu) == 'A' && (c[1] & 0xDFu) == 'T' && (c[2] & 0xDFu) == 'G';
    return duckvep_codon_is_start(c, table);
}

static int stop_codon_of(const uint8_t *c, duckvep_codon_table_t table) {
    char codon[3];
    if (table == DUCKVEP_CODON_TABLE_STANDARD) return stop_codon(c);
    for (unsigned i = 0u; i < 3u; i++) codon[i] = (char)(c[i] & 0xDFu);
    return duckvep_translate_codon(codon, table) == '*';
}

/* A CDS that begins inside a codon (exon phase 1 or 2, only where the start is not annotated) is stored as
 * Ensembl translates it: padded in front with one N per missing base, so the stored bytes are in frame and
 * the first residue is unknown. */
static int phase_padding_ok(const uint8_t *cds, unsigned phase, int open_start) {
    if (!phase) return 1;
    if (!open_start || phase > 2u) return 0;
    for (unsigned i = 0u; i < phase; i++) if ((cds[i] & 0xDFu) != 'N') return 0;
    return 1;
}

/* CDS in a supported genetic code, stored in frame (see phase_padding_ok), with no curated RNA/peptide edit
 * or recoding. A complete CDS has a start codon, a terminal stop and no internal stop. A CDS whose start is
 * not annotated (cds_start_NF) need not begin with a start codon; one whose end is not annotated (cds_end_NF)
 * has no stop at all and may end in a partial codon. Cached per transcript. */
static duckvep_prediction_reason_t transcript_domain(duckvep_haplotype_stream_t *s, uint32_t tx) {
    if (s->have_domain && s->domain_transcript == tx) return s->domain_reason;
    const duckvep_transcript_model_t *m = s->carriers.model;
    const duckvep_sequence_pool_t *seq = s->sequences;
    duckvep_prediction_reason_t r = DUCKVEP_REASON_SUPPORTED_DOMAIN;
    size_t length = seq->cds_length ? seq->cds_length[tx] : 0u;
    uint64_t offset = seq->cds_offset ? seq->cds_offset[tx] : 0u;
    const uint64_t curated = DUCKVEP_TX_SELENOCYSTEINE | DUCKVEP_TX_STOP_CODON_READTHROUGH |
        DUCKVEP_TX_RNA_EDIT | DUCKVEP_TX_AMINO_ACID_SUB;
    uint64_t flags = m->flags ? m->flags[tx] : 0u;
    if (!length || !m->cds_start1 || !m->cds_start1[tx] || offset > seq->cds_bytes_len ||
        length > seq->cds_bytes_len - offset) r = DUCKVEP_REASON_TRANSCRIPT_NOT_CODING;
    else if (seq->codon_table && !duckvep_codon_table_supported((duckvep_codon_table_t)seq->codon_table[tx]))
        r = DUCKVEP_REASON_NON_STANDARD_CODON_TABLE;
    else if ((flags & curated) || (seq->peptide_edit_offset &&
             seq->peptide_edit_offset[tx + 1u] != seq->peptide_edit_offset[tx]))
        r = DUCKVEP_REASON_CURATED_TRANSCRIPT;
    else if (length < 6u || (length % 3u && !(flags & DUCKVEP_TX_CDS_END_NF)) ||
             !phase_padding_ok(seq->cds_bytes + (size_t)offset, m->cds_phase_offset ? m->cds_phase_offset[tx] : 0u,
                               (flags & DUCKVEP_TX_CDS_START_NF) != 0u))
        r = DUCKVEP_REASON_INCOMPLETE_CDS; /* too short, truncated without the flag, or an unpadded later phase */
    else {
        const uint8_t *c = seq->cds_bytes + (size_t)offset;
        duckvep_codon_table_t table = seq->codon_table
            ? (duckvep_codon_table_t)seq->codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
        size_t codons = length / 3u, open_end = (flags & DUCKVEP_TX_CDS_END_NF) != 0u;
        if (!(flags & DUCKVEP_TX_CDS_START_NF) && !start_codon_of(c, table)) r = DUCKVEP_REASON_NONCANONICAL_START;
        else if (!open_end && !stop_codon_of(c + length - 3u, table)) r = DUCKVEP_REASON_NONCANONICAL_STOP;
        else for (size_t i = 1u; i + (open_end ? 0u : 1u) < codons; i++)
            if (stop_codon_of(c + 3u * i, table)) { r = DUCKVEP_REASON_INTERNAL_STOP; break; }
        if (r == DUCKVEP_REASON_SUPPORTED_DOMAIN && open_end && stop_codon_of(c, table))
            r = DUCKVEP_REASON_INTERNAL_STOP;
    }
    s->have_domain = 1u;
    s->domain_transcript = tx;
    s->domain_reason = r;
    return r;
}

/* The reference lane has no source event or carrier evidence. It reuses the
 * model-owned CDS and the stream's cached reference translation. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_reference(
    duckvep_haplotype_stream_t *s, uint32_t tx, duckvep_haplotype_leaf_t *leaf) {
    const duckvep_sequence_pool_t *seq;
    duckvep_prediction_reason_t reason;
    duckvep_haplotype_stream_status_t status;
    size_t length;
    uint64_t offset;
    if (!s || !s->initialized || !leaf || !s->sequences || tx >= s->sequences->transcript_count)
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    seq = s->sequences;
    length = seq->cds_length[tx];
    offset = seq->cds_offset[tx];
    if (offset > seq->cds_bytes_len || length > seq->cds_bytes_len - offset)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    status = prepare_reference_protein(s, tx);
    if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    reason = transcript_domain(s, tx);
    memset(leaf, 0, sizeof(*leaf));
    leaf->carriers.transcript_index = tx;
    leaf->reference_cds = seq->cds_bytes + (size_t)offset;
    leaf->cds = leaf->reference_cds;
    leaf->cds_length = length;
    leaf->reference_protein = s->reference_protein_known ? s->buffers.reference_protein : NULL;
    leaf->reference_protein_length = s->reference_protein_length;
    leaf->reference_coding_protein = s->buffers.reference_coding_protein;
    leaf->reference_coding_translation = s->reference_coding_translation;
    leaf->protein = leaf->reference_protein;
    leaf->protein_length = s->reference_protein_length;
    leaf->translation = s->reference_coding_translation;
    leaf->projection_status = DUCKVEP_CDS_EDIT_OK;
    leaf->sequence_status = s->reference_protein_known ? DUCKVEP_HAPLOTYPE_OK : DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE;
    leaf->prediction_reason = reason;
    leaf->path_reason = reason;
    leaf->prediction_status = (reason == DUCKVEP_REASON_SUPPORTED_DOMAIN && s->reference_protein_known)
        ? DUCKVEP_PREDICTION_ELIGIBLE : DUCKVEP_PREDICTION_INCOMPLETE_INPUT;
    leaf->path_status = leaf->prediction_status;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

static int literal_acgt(const uint8_t *b, size_t n) {
    for (size_t i = 0u; i < n; i++) {
        uint8_t c = b[i] & 0xDFu;
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return 0;
    }
    return 1;
}

/* First offending adjacent pair of ascending-CDS edits (descending array in the buffer). */
static duckvep_prediction_reason_t edit_relation(const duckvep_haplotype_edit_t *edits, size_t n) {
    duckvep_prediction_reason_t overlap = DUCKVEP_REASON_SUPPORTED_DOMAIN;
    for (size_t i = 1u; i < n; i++) {
        const duckvep_haplotype_edit_t *lo = &edits[i], *hi = &edits[i - 1u];
        if (lo->cds_start == hi->cds_start && !lo->ref_len && !hi->ref_len) {
            if (!overlap) overlap = DUCKVEP_REASON_SAME_GAP_INSERTIONS;
        } else if (lo->cds_start == hi->cds_start && lo->ref_len && lo->ref_len == hi->ref_len) {
            if (lo->variant_strand == hi->variant_strand && lo->alt_len == hi->alt_len &&
                (!lo->alt_len || !memcmp(lo->alt, hi->alt, lo->alt_len))) {
                if (!overlap) overlap = DUCKVEP_REASON_DUPLICATE_EDITS;
            } else return DUCKVEP_REASON_CONTRADICTORY_EDITS;
        } else if (lo->ref_len > hi->cds_start - lo->cds_start) {
            if (!overlap) overlap = DUCKVEP_REASON_OVERLAPPING_EDITS;
        }
    }
    return overlap;
}

/* Whole-haplotype classifier. Runs only on a path that is inside the supported
 * domain. It classifies the whole edited peptide against the uncurated reference peptide, never an edit
 * alone. The decision follows the translated sequence of the edited CDS, not the nominal net frame offset
 * of its edits, in this order:
 *   edited CDS does not begin with a start codon (ATG in the
 *     standard code, any start codon of another code)        -> start_lost (alone: initiation is unknown, so
 *                                                               start loss suppresses every other prediction)
 *   no stop in the edited CDS                                -> stop_lost, plus frameshift_variant when the
 *                                                               frame is still displaced when the CDS runs out;
 *                                                               no downstream extension is ever invented
 *   first stop starting before the homologous reference terminator, unless the peptide before it is the
 *     unchanged reference peptide (an inserted stop codon next to the terminator removes no residue)
 *                                                            -> stop_gained, plus frameshift_variant when the
 *                                                               stop codon intersects a displaced-frame interval
 *   first stop at or inside the terminator's window (or after it, which is stop_lost as above):
 *     stop in a displaced-frame interval                     -> frameshift_variant
 *     identical peptide, terminal codon unchanged            -> synonymous_variant
 *     identical peptide, terminal codon changed or moved     -> stop_retained_variant (LOW, precedes synonymous)
 *     substitutions only, changed peptide                    -> missense_variant
 *     every edit a pure insertion / pure deletion            -> inframe_insertion / inframe_deletion
 *     any other change, including a frame that was displaced and restored before termination
 *                                                            -> protein_altering_variant
 * The terminator window is the interval of the edited CDS that holds the reference terminator (the last three
 * reference bases) after the ascending edit islands are applied: it starts at the last unedited boundary
 * before the terminator, or earlier when an island starts before it and overlaps it, and ends at the CDS end.
 * Displaced-frame intervals are geometric facts of the ascending edit islands (leaf.stop_in_displaced_frame,
 * from duckvep_haplotype_block_frame_intersects): they start at the first frame-changing edit and end
 * after the ALT bases of the restoring edit, so an early stop inside an open interval keeps the frame
 * term even when a later edit would have restored the frame, while a restoration before termination
 * removes it. Edits after the first stop stay contributors (role post_stop) and never an expressed effect.
 * Purity is a property of the normalized edit path (differing islands), as the contract pins for frame SO.
 * No eligible path is left pending: the reference is a complete table-1 CDS with its first stop at the
 * terminator (domain check), so the remaining guards are defensive and report unsupported_context. */
/* The edited peptide before its first stop equals the reference peptide before the terminator. The first
 * residue is the initiator: the caller has established that both first codons are start codons, so it is
 * the same residue even where a genetic code translates two of its start codons differently internally. */
static int same_peptide(const duckvep_haplotype_stream_t *s, const duckvep_haplotype_leaf_t *leaf,
                        size_t first_stop) {
    uint32_t tx = leaf->carriers.transcript_index;
    uint64_t flags = s->carriers.model->flags ? s->carriers.model->flags[tx] : 0u;
    /* Without an annotated start the first residue is an ordinary one, and without an annotated end the
     * reference peptide is every complete codon and the edited one must have as many, with no stop. */
    size_t skip = (flags & DUCKVEP_TX_CDS_START_NF) ? 0u : 1u;
    size_t ref_n = s->sequences->cds_length[tx] / 3u, alt_n;
    if (flags & DUCKVEP_TX_CDS_END_NF) {
        if (first_stop) return 0;
        alt_n = leaf->translation.length;
    } else {
        ref_n -= 1u;
        alt_n = first_stop ? first_stop - 1u : SIZE_MAX;
    }
    return ref_n == alt_n && ref_n >= skip &&
        !memcmp(leaf->reference_coding_protein + skip, s->buffers.protein + skip, ref_n - skip);
}

/* ejc50, decided on the edited spliced transcript of the shared path and never per
 * allele. The 5' UTR is not edited (only CDS edits are applied), so the edited stop sits at the CDS cDNA origin
 * plus its edited CDS offset, and the penultimate exon's last base moves by the length change of every edit
 * that starts at or before it: indels upstream and inside that exon shift J, indels after it (the last exon)
 * do not, and an insertion between the two exons belongs to the last exon. Edits after the stop are part of
 * the edited transcript, so a post-stop indel before J shifts J (S never moves); one at or after J leaves the
 * prediction unchanged. Unresolvable exon topology is unknown, never a guess. */
static void nmd_ejc50(const duckvep_haplotype_stream_t *s, duckvep_haplotype_leaf_t *leaf, uint64_t mask,
                      size_t first_stop) {
    const duckvep_transcript_model_t *m = s->carriers.model;
    const duckvep_exon_model_t *x = s->exons;
    leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_UNKNOWN;
    leaf->nmd_stop_valid = leaf->nmd_junction_valid = 0u;
    if (mask & DUCKVEP_SO(DUCKVEP_SO_START_LOST)) return;          /* lost initiation */
    const int open_end = m && m->flags && (m->flags[leaf->carriers.transcript_index] & DUCKVEP_TX_CDS_END_NF);
    if (open_end && !(mask & (DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED) | DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)))) {
        leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE;          /* termination was not annotated and is not changed */
        return;
    }
    if (!first_stop) return;                                        /* no termination: run-off, stop_lost */
    if (!(mask & DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED))) {             /* known termination, nothing premature */
        leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE;
        return;
    }
    uint32_t tx = leaf->carriers.transcript_index;
    if (open_end) return;                                           /* the last junction is not annotated */
    if (!m || !x || !m->exon_offset || !m->exon_count || !x->cdna_start1 || !x->cdna_end1 || !x->start1 ||
        !x->end1 || !m->strand || !m->cds_start1 || !m->cds_end1) return;
    size_t ref_length = s->sequences->cds_length[tx];
    size_t n = m->exon_count[tx], first = m->exon_offset[tx];
    if (!n || first > x->exon_count || n > x->exon_count - first) return;
    /* cDNA position of the first CDS base: the model's cache when present, else from the exon holding it. */
    uint64_t origin = m->cds_cdna_start1 ? m->cds_cdna_start1[tx] : 0u;
    if (!origin) {
        uint32_t coding_first = m->strand[tx] < 0 ? m->cds_end1[tx] : m->cds_start1[tx];
        for (size_t i = first; i < first + n && !origin; i++)
            if (coding_first >= x->start1[i] && coding_first <= x->end1[i])
                origin = (uint64_t)x->cdna_start1[i] + (m->strand[tx] < 0 ? x->end1[i] - coding_first
                                                                            : coding_first - x->start1[i]);
    }
    if (!origin) return;
    /* Phase padding (leading N) is in the stored CDS but not in the transcript. */
    const uint64_t padding = m->cds_phase_offset ? m->cds_phase_offset[tx] : 0u;
    if ((uint64_t)first_stop * 3u <= padding || ref_length <= padding) return;
    ref_length -= (size_t)padding;
    uint64_t stop = origin + (uint64_t)first_stop * 3u - 1u - padding;
    uint64_t total = 0u, last_start = 0u, pen_start = 0u, last_end = 0u, pen_end = 0u;
    for (size_t i = first; i < first + n; i++) {
        uint64_t a = x->cdna_start1[i], z = x->cdna_end1[i];
        if (!a || z < a) return;
        total += z - a + 1u;
        if (a > last_start) { pen_start = last_start; pen_end = last_end; last_start = a; last_end = z; }
        else if (a > pen_start) { pen_start = a; pen_end = z; }
    }
    /* Exons tile the spliced transcript from cDNA 1 without gaps or overlaps, and hold the CDS. */
    if (total != last_end || (n > 1u && (!pen_start || pen_end + 1u != last_start)) ||
        origin + ref_length - 1u > last_end) { leaf->nmd_stop_valid = 0u; return; }
    leaf->nmd_stop_position1 = stop;
    leaf->nmd_stop_valid = 1u;
    const duckvep_haplotype_stream_buffers_t *eb = &s->buffers;
    if ((uint64_t)first_stop * 3u - padding <= DUCKVEP_HAPLOTYPE_NMD_START_PROXIMAL_BASES)
        leaf->nmd_exceptions |= DUCKVEP_HAPLOTYPE_NMD_EXCEPTION_START_PROXIMAL;
    /* The exon of the edited transcript that holds the stop. Exons are stored in cDNA order. */
    uint64_t edited_end = 0u, previous_end = 0u;
    for (size_t i = first; i < first + n; i++) {
        uint64_t a = x->cdna_start1[i], z = x->cdna_end1[i];
        if (a != previous_end + 1u) break;                           /* not in cDNA order: no statement */
        previous_end = z;
        int64_t size = (int64_t)(z - a + 1u);
        for (size_t j = 0u; j < leaf->edit_count; j++) {
            uint64_t q0 = origin + (uint64_t)eb->edits[j].cds_start - 1u - padding;
            if (q0 >= a && q0 <= z) size += (int64_t)eb->edits[j].alt_len - (int64_t)eb->edits[j].ref_len;
        }
        if (size <= 0) break;
        if (stop <= edited_end + (uint64_t)size) {
            if (size > DUCKVEP_HAPLOTYPE_NMD_LONG_EXON_BASES)
                leaf->nmd_exceptions |= DUCKVEP_HAPLOTYPE_NMD_EXCEPTION_LONG_EXON;
            break;
        }
        edited_end += (uint64_t)size;
    }
    if (n == 1u) {                                                  /* intronless: no junction, always escapes */
        leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_ESCAPE;
        return;
    }
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    int64_t junction = (int64_t)pen_end, last_size = (int64_t)(last_end - last_start + 1u),
            pen_size = (int64_t)(pen_end - pen_start + 1u);
    for (size_t i = 0u; i < leaf->edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &b->edits[i];
        uint64_t q0 = origin + (uint64_t)e->cds_start - 1u - padding;
        int64_t change = (int64_t)e->alt_len - (int64_t)e->ref_len;
        if (q0 <= pen_end) {
            if (e->ref_len && q0 + e->ref_len - 1u > pen_end) { leaf->nmd_stop_valid = 0u; return; } /* an edit may not span the junction */
            junction += change;
            if (change)
                for (size_t j = 0u; j < leaf->contributor_count; j++)
                    if (b->contributors[j].source.event_id == b->edit_event_ids[i]) b->contributors[j].nmd_moved_junction = 1u;
            if (q0 >= pen_start) pen_size += change;
        } else if (q0 >= last_start) {
            if (e->ref_len && q0 + e->ref_len - 1u > last_end) { leaf->nmd_stop_valid = 0u; return; }
            last_size += change;
        }
    }
    if (pen_size <= 0 || last_size <= 0 || junction < 1) { leaf->nmd_stop_valid = 0u; return; } /* an exon was deleted whole */
    leaf->nmd_junction_position1 = (uint64_t)junction;
    leaf->nmd_junction_valid = 1u;
    leaf->nmd = junction - (int64_t)stop > DUCKVEP_HAPLOTYPE_NMD_THRESHOLD ? DUCKVEP_HAPLOTYPE_NMD_TRIGGER
                                                                            : DUCKVEP_HAPLOTYPE_NMD_ESCAPE;
}

/* Reading through a lost stop: translation continues from the last bases of the edited CDS into the
 * transcript's stored 3' flank, to the first stop or the end of the flank. The residues are appended to
 * the path's protein. Nothing is appended when the model has no flank for the transcript or the protein
 * buffer cannot hold the whole possible extension, so a protein is never cut short silently. */
static void extend_past_cds(duckvep_haplotype_stream_t *s, duckvep_haplotype_leaf_t *leaf,
                            duckvep_codon_table_t table) {
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    const duckvep_sequence_pool_t *seq = s->sequences;
    uint32_t tx = leaf->carriers.transcript_index;
    if (leaf->protein != b->protein || !seq->flank_bytes || !seq->post_cds_offset || !seq->post_cds_length) return;
    uint64_t offset = seq->post_cds_offset[tx];
    size_t flank = seq->post_cds_length[tx], tail = leaf->cds_length % 3u, at = leaf->translation.length;
    if (!flank || offset > seq->flank_bytes_len || flank > seq->flank_bytes_len - offset ||
        at >= b->protein_capacity || (tail + flank) / 3u + 1u > b->protein_capacity - at) return;
    const uint8_t *bases = seq->flank_bytes + (size_t)offset;
    char codon[3];
    size_t filled = 0u;
    for (size_t i = 0u; i < tail; i++) codon[filled++] = (char)leaf->cds[leaf->cds_length - tail + i];
    for (size_t i = 0u; i < flank; i++) {
        codon[filled++] = (char)bases[i];
        if (filled < 3u) continue;
        filled = 0u;
        char residue = duckvep_translate_codon(codon, table);
        b->protein[at++] = (uint8_t)residue;
        if (residue == '*') break;
    }
    b->protein[at] = 0u;
    leaf->protein_length = at;
    leaf->flags |= DUCKVEP_HAPLOTYPE_FLAG_EXTENDED;
}

static void classify_haplotype(duckvep_haplotype_stream_t *s, duckvep_haplotype_leaf_t *leaf) {
    if (leaf->path_status != DUCKVEP_PREDICTION_ELIGIBLE || !leaf->cds || leaf->ordered_replacements) return;
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    uint32_t tx = leaf->carriers.transcript_index;
    size_t ref_length = s->sequences->cds_length[tx], alt_length = leaf->cds_length;
    if (!leaf->edit_count) {
        if (alt_length != ref_length || memcmp(leaf->cds, leaf->reference_cds, ref_length)) return;
        leaf->path_status = DUCKVEP_PREDICTION_PREDICTED;
        leaf->haplotype_so_mask = 0u;
        leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE; /* a reference lane has no termination change */
        return;
    }
    const uint64_t tx_flags = s->carriers.model->flags ? s->carriers.model->flags[tx] : 0u;
    const int open_start = (tx_flags & DUCKVEP_TX_CDS_START_NF) != 0u, open_end = (tx_flags & DUCKVEP_TX_CDS_END_NF) != 0u;
    if (!leaf->reference_coding_protein || !leaf->reference_cds || ref_length < 6u ||
        (!open_end && (ref_length % 3u || leaf->reference_coding_translation.first_stop_position1 != ref_length / 3u)) ||
        (open_end && leaf->reference_coding_translation.first_stop_position1 != 0u)) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_INVALID_SEQUENCE;
        return;
    }
    /* The phase padding is the one place an N is expected; every other base must be literal. */
    const size_t padding = s->carriers.model->cds_phase_offset ? s->carriers.model->cds_phase_offset[tx] : 0u;
    if ((!leaf->translation.unambiguous || !leaf->reference_coding_translation.unambiguous) &&
        (!padding || alt_length <= padding || !literal_acgt(leaf->cds + padding, alt_length - padding) ||
         !literal_acgt(leaf->reference_cds + padding, ref_length - padding))) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_INVALID_BASE;
        return;
    }
    /* Edit-path facts and the terminator window [window_begin, window_end) of the edited CDS. */
    const size_t terminator0 = ref_length - 3u;
    int frame = 0, substitutions = 1, insertions = 1, deletions = 1, touched = 0;
    int64_t shift_before = 0, after_change = 0;
    size_t window_begin = 0u;
    for (size_t i = 0u; i < leaf->edit_count; i++) {
        const duckvep_haplotype_edit_t *e = &b->edits[i];
        size_t r0 = (size_t)e->cds_start - 1u, r1 = r0 + e->ref_len;
        int64_t change = (int64_t)e->alt_len - (int64_t)e->ref_len;
        if ((e->alt_len % 3u) != (e->ref_len % 3u)) frame = 1;
        substitutions &= e->ref_len == e->alt_len;
        insertions &= e->ref_len == 0u && e->alt_len != 0u;
        deletions &= e->alt_len == 0u && e->ref_len != 0u;
        if (r1 <= terminator0) shift_before += change;
        else if (r0 >= ref_length) after_change += change;
        else if (!touched) {
            touched = 1;
            window_begin = (size_t)((int64_t)(r0 < terminator0 ? r0 : terminator0) + shift_before);
        }
    }
    if (!touched) window_begin = (size_t)((int64_t)terminator0 + shift_before);
    /* An edit inside the padded first codon changes a residue that is unknown on both sides: an unchanged
     * peptide is then a coding_sequence_variant, as VEP reports it, not a synonymous one. */
    int unknown_first = 0;
    for (size_t i = 0u; padding && i < leaf->edit_count; i++)
        unknown_first |= b->edits[i].cds_start <= 3u;
    const size_t window_end = (size_t)((int64_t)alt_length - after_change);
    const size_t first_stop = leaf->translation.first_stop_position1;
    uint64_t mask;
    const duckvep_codon_table_t table = s->sequences->codon_table
        ? (duckvep_codon_table_t)s->sequences->codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
    if (!open_start && (alt_length < 3u || !start_codon_of(leaf->cds, table))) {
        mask = DUCKVEP_SO(DUCKVEP_SO_START_LOST);
    } else if (open_end) {
        /* The reference has no annotated termination, so nothing can be lost or retained: a stop is new, a
         * frame still displaced where the annotation ends is a frameshift, and otherwise the peptides of
         * the complete codons are compared. */
        if (first_stop) {
            mask = DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
            if (leaf->stop_in_displaced_frame) mask |= DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
        } else if ((alt_length % 3u) != (ref_length % 3u)) mask = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
        else if (same_peptide(s, leaf, 0u)) {
            /* An unchanged peptide whose only edits lie in the trailing partial codon says nothing about
             * that residue: the term VEP uses for it, not synonymous. */
            int partial_only = (ref_length % 3u) != 0u;
            for (size_t i = 0u; i < leaf->edit_count && partial_only; i++)
                partial_only = (size_t)b->edits[i].cds_start - 1u >= ref_length - ref_length % 3u;
            mask = DUCKVEP_SO(partial_only ? DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON :
                              unknown_first ? DUCKVEP_SO_CODING_SEQUENCE : DUCKVEP_SO_SYNONYMOUS);
        }
        else if (frame) mask = DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
        else if (substitutions) mask = DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
        else if (insertions) mask = DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION);
        else if (deletions) mask = DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION);
        else mask = DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
    } else if (!first_stop) {
        mask = DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
        if (alt_length % 3u) mask |= DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
    } else if ((first_stop - 1u) * 3u < window_begin && !same_peptide(s, leaf, first_stop)) {
        /* A new first stop before the homologous reference terminator that truncates the reference
         * peptide. A stop that follows an unchanged peptide (an inserted stop codon next to the terminator)
         * removes no residue and is judged as the terminal codon below. */
        mask = DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
        if (leaf->stop_in_displaced_frame) mask |= DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
    } else if ((first_stop - 1u) * 3u >= window_end) {
        /* The first stop lies beyond the terminator's window (only after-edit bases can be there). */
        mask = DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
        if (leaf->stop_in_displaced_frame) mask |= DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
    } else {
        /* Termination is present at the terminator: read in frame, or a displaced-frame stop inside its window. */
        size_t stop0 = (first_stop - 1u) * 3u;
        int same = same_peptide(s, leaf, first_stop);
        int terminal_changed = stop0 != window_begin;
        for (size_t k = 0u; k < 3u && !terminal_changed; k++)
            terminal_changed = (leaf->cds[stop0 + k] & 0xDFu) != (leaf->reference_cds[terminator0 + k] & 0xDFu);
        if (leaf->stop_in_displaced_frame) mask = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
        else if (same) mask = DUCKVEP_SO(terminal_changed ? DUCKVEP_SO_STOP_RETAINED :
                                         unknown_first ? DUCKVEP_SO_CODING_SEQUENCE : DUCKVEP_SO_SYNONYMOUS);
        else if (frame) mask = DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
        else if (substitutions) mask = DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
        else if (insertions) mask = DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION);
        else if (deletions) mask = DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION);
        else mask = DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
    }
    leaf->path_status = DUCKVEP_PREDICTION_PREDICTED;
    leaf->haplotype_so_mask = mask;
    if ((mask & DUCKVEP_SO(DUCKVEP_SO_STOP_LOST)) && !first_stop) extend_past_cds(s, leaf, table);
    nmd_ejc50(s, leaf, mask, first_stop);
}

static uint8_t complement_base(uint8_t base) {
    switch (base & 0xDFu) {
    case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    default: return base;
    }
}

/* Build a final, joint cDNA allele only for literal, length-preserving substitutions
 * confined to one exon. The caller-owned workspace isolates RNA replay from the
 * raw CDS and protein output buffers. */
static int noncoding_haplotype_mask(duckvep_haplotype_stream_t *s,
    const duckvep_haplotype_leaf_t *leaf, uint64_t *mask) {
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    const uint32_t tx = leaf->carriers.transcript_index;
    const duckvep_sequence_pool_t *seq = s->sequences;
    const duckvep_transcript_model_t *m = s->carriers.model;
    size_t pre, cds, post, total;
    uint32_t coding_begin = 0u, coding_end = 0u;
    int transcript_noncoding;

    if (!leaf->contributor_count || tx >= seq->transcript_count || !m->cds_start1 ||
        !m->cds_end1)
        return 0;
    for (size_t i = 0u; i < leaf->contributor_count; i++)
        if (leaf->contributors[i].projection_status != DUCKVEP_CDS_EDIT_OUT_OF_CDS)
            return 0;
    transcript_noncoding = m->cds_start1[tx] == 0u && m->cds_end1[tx] == 0u;
    if (!b->noncoding)
        return -1;
    if (transcript_noncoding) {
        if (!seq->cdna_provided || !seq->cdna_bytes || !seq->cdna_offset ||
            !seq->cdna_length || seq->cdna_offset[tx] > seq->cdna_bytes_len ||
            seq->cdna_length[tx] == 0u || seq->cdna_length[tx] >
            seq->cdna_bytes_len - seq->cdna_offset[tx])
            return -1;
        pre = cds = post = 0u;
        total = seq->cdna_length[tx];
        if (total > b->noncoding_capacity / 2u)
            return -1;
        memcpy(b->noncoding, seq->cdna_bytes + seq->cdna_offset[tx], total);
    } else {
        duckvep_transcript_coordinate_t coding_first, coding_last;
        size_t transcript_length = 0u;
        size_t exon_end;

        if (!seq->flanks_complete || !seq->flank_bytes || !seq->pre_cds_offset ||
            !seq->pre_cds_length || !seq->post_cds_offset || !seq->post_cds_length ||
            !seq->cds_length)
            return -1;
        pre = seq->pre_cds_length[tx];
        cds = seq->cds_length[tx];
        post = seq->post_cds_length[tx];
        if (pre > SIZE_MAX - cds || pre + cds > SIZE_MAX - post ||
            seq->pre_cds_offset[tx] > seq->flank_bytes_len ||
            pre > seq->flank_bytes_len - seq->pre_cds_offset[tx] ||
            seq->post_cds_offset[tx] > seq->flank_bytes_len ||
            post > seq->flank_bytes_len - seq->post_cds_offset[tx])
            return -1;
        total = pre + cds + post;
        if (total > b->noncoding_capacity / 2u || pre >= UINT32_MAX ||
            cds > UINT32_MAX - pre || !leaf->reference_cds ||
            leaf->cds_length != cds)
            return -1;
        if (duckvep_project_transcript_coordinate(m, s->exons, tx,
            m->cds_start1[tx], &coding_first) ||
            duckvep_project_transcript_coordinate(m, s->exons, tx,
            m->cds_end1[tx], &coding_last) || !coding_first.exonic ||
            !coding_last.exonic)
            return -1;
        coding_begin = coding_first.cdna_anchor1 < coding_last.cdna_anchor1
            ? coding_first.cdna_anchor1 : coding_last.cdna_anchor1;
        coding_end = coding_first.cdna_anchor1 > coding_last.cdna_anchor1
            ? coding_first.cdna_anchor1 : coding_last.cdna_anchor1;
        exon_end = (size_t)m->exon_offset[tx] + m->exon_count[tx];
        for (size_t e = m->exon_offset[tx]; e < exon_end; e++)
            if (s->exons->cdna_end1[e] > transcript_length)
                transcript_length = s->exons->cdna_end1[e];
        if (coding_begin != pre + 1u || coding_end != pre + cds ||
            transcript_length != total)
            return -1;
        memcpy(b->noncoding, seq->flank_bytes + seq->pre_cds_offset[tx], pre);
        memcpy(b->noncoding + pre, leaf->reference_cds, cds);
        memcpy(b->noncoding + pre + cds,
            seq->flank_bytes + seq->post_cds_offset[tx], post);
    }
    memcpy(b->noncoding + total, b->noncoding, total);
    int changed = 0;
    for (size_t i = 0u; i < leaf->contributor_count; i++) {
        const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
        const duckvep_event_t *event = c->prepared;
        duckvep_transcript_coordinate_t first, last;
        size_t n, start0;
        if (!event || !event->ref_diff_length ||
            event->ref_diff_length != event->alt_diff_length || event->interbase ||
            event->start1 > UINT32_MAX - event->ref_diff_length + 1u ||
            duckvep_project_transcript_coordinate(s->carriers.model, s->exons, tx, event->start1, &first) ||
            duckvep_project_transcript_coordinate(s->carriers.model, s->exons, tx,
                event->start1 + event->ref_diff_length - 1u, &last) ||
            !first.exonic || !last.exonic || first.exon_idx != last.exon_idx || !first.cdna_anchor1 ||
            !last.cdna_anchor1)
            return 0;
        n = event->ref_diff_length;
        start0 = first.cdna_anchor1 < last.cdna_anchor1 ? first.cdna_anchor1 - 1u : last.cdna_anchor1 - 1u;
        if (start0 > total || n > total - start0 ||
            (!transcript_noncoding &&
            !(start0 + n <= pre || start0 >= pre + cds)) ||
            !literal_acgt(c->source.ref + event->ref_diff_offset, n) ||
            !literal_acgt(c->source.alt + event->alt_diff_offset, n))
            return 0;
        size_t exon_begin = s->exons->cdna_start1[first.exon_idx];
        size_t exon_limit = s->exons->cdna_end1[first.exon_idx];
        if ((exon_begin > 1u && start0 - (exon_begin - 1u) < 3u) ||
            (exon_limit < total && exon_limit - (start0 + n) < 3u)) return 0;
        for (size_t j = 0u; j < n; j++) {
            size_t source_j = s->carriers.model->strand[tx] > 0 ? j : n - 1u - j;
            uint8_t ref = c->source.ref[event->ref_diff_offset + source_j];
            uint8_t alt = c->source.alt[event->alt_diff_offset + source_j];
            if (s->carriers.model->strand[tx] < 0) { ref = complement_base(ref); alt = complement_base(alt); }
            if ((b->noncoding[total + start0 + j] & 0xDFu) != (ref & 0xDFu)) return -1;
            if ((b->noncoding[total + start0 + j] & 0xDFu) != (alt & 0xDFu)) changed = 1;
            b->noncoding[total + start0 + j] = alt;
        }
    }
    if (!changed) { *mask = 0u; return 1; }
    return duckvep_effect_eval_haplotype_noncoding(b->noncoding, total,
        b->noncoding + total, total, transcript_noncoding ? 0u :
        (uint32_t)pre + 1u, transcript_noncoding ? 0u :
        (uint32_t)(pre + cds), mask) == DUCKVEP_HAPLOTYPE_NONCODING_OK ? 1 : -1;
}

/* Runs on every leaf after sequence construction. It changes no existing field except
 * ordering the (already exposed) edit buffers of failed decoded-call leaves. */
static int conditional_reference_is_closed(const duckvep_haplotype_leaf_t *leaf)
{
    size_t i;

    if (!leaf->reference_protein || leaf->reference_protein_length == 0u ||
        leaf->reference_protein[leaf->reference_protein_length - 1u] != (uint8_t)'*')
        return 0;
    for (i = 0u; i + 1u < leaf->reference_protein_length; i++)
        if (leaf->reference_protein[i] == (uint8_t)'*')
            return 0;
    return 1;
}

static int conditional_contributors_are_cds_substitutions(
    const duckvep_haplotype_leaf_t *leaf)
{
    size_t i;

    if (!leaf->contributor_count || leaf->ordered_replacements || !leaf->edit_count)
        return 0;
    for (i = 0u; i < leaf->contributor_count; i++) {
        const duckvep_haplotype_contributor_t *contributor = &leaf->contributors[i];
        const duckvep_haplotype_edit_t *edit = contributor->projected;

        if (contributor->projection_status != DUCKVEP_CDS_EDIT_OK || !edit ||
            !contributor->edit_count || edit->ref_len == 0u ||
            edit->ref_len != edit->alt_len)
            return 0;
    }
    for (i = 0u; i < leaf->block_count; i++)
        if (leaf->blocks[i].ref_len == 0u ||
            leaf->blocks[i].ref_len != leaf->blocks[i].alt_len)
            return 0;
    return 1;
}

static void conditional_prediction(duckvep_haplotype_stream_t *s,
    duckvep_haplotype_leaf_t *leaf)
{
    const duckvep_sequence_pool_t *seq = s->sequences;
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    uint32_t tx = leaf->carriers.transcript_index;
    uint64_t offset;
    size_t begin, end, prediction_length, first_stop;
    uint64_t flags;
    duckvep_codon_table_t table;
    duckvep_haplotype_status_t status;

    if (leaf->path_reason != DUCKVEP_REASON_CURATED_TRANSCRIPT)
        return;
    if (!seq->peptide_edit_code) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_UNTYPED_CURATED_METADATA;
        return;
    }
    flags = s->carriers.model->flags ? s->carriers.model->flags[tx] : 0u;
    if ((flags & DUCKVEP_TX_RNA_EDIT) ||
        !leaf->cds || !leaf->reference_cds || !s->reference_protein_known ||
        leaf->cds_length != seq->cds_length[tx] || leaf->cds_length % 3u ||
        !conditional_reference_is_closed(leaf) ||
        !conditional_contributors_are_cds_substitutions(leaf)) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_UNSUPPORTED_CURATED_EDIT;
        return;
    }
    offset = seq->cds_offset[tx];
    begin = seq->peptide_edit_offset[tx];
    end = seq->peptide_edit_offset[tx + 1u];
    table = seq->codon_table ? (duckvep_codon_table_t)seq->codon_table[tx] :
        DUCKVEP_CODON_TABLE_STANDARD;
    if (offset > seq->cds_bytes_len || leaf->cds_length > seq->cds_bytes_len - offset ||
        begin >= end || end > seq->peptide_edit_count) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_UNSUPPORTED_CURATED_EDIT;
        return;
    }
    status = duckvep_haplotype_conditional_peptide(leaf->reference_cds, leaf->cds,
        leaf->cds_length, table, seq->peptide_edit_position1 + begin,
        seq->peptide_edit_alt + begin, seq->peptide_edit_code + begin, end - begin,
        (flags & (DUCKVEP_TX_CDS_START_NF | DUCKVEP_TX_CDS_END_NF)) == 0u,
        b->prediction_protein, b->prediction_protein_capacity, &prediction_length,
        &first_stop);
    if (status != DUCKVEP_HAPLOTYPE_OK) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_UNSUPPORTED_CURATED_EDIT;
        return;
    }
    leaf->prediction_reference_protein = leaf->reference_protein;
    leaf->prediction_reference_protein_length = leaf->reference_protein_length;
    leaf->prediction_protein = b->prediction_protein;
    leaf->prediction_protein_length = prediction_length;
    leaf->path_status = DUCKVEP_PREDICTION_CONDITIONAL_RECODING;
    leaf->path_reason = DUCKVEP_REASON_ASSUMED_RECODING_PROGRAMME;
    leaf->haplotype_so_mask = 0u;
    leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_UNKNOWN;
    leaf->nmd_stop_valid = 0u;
    leaf->nmd_junction_valid = 0u;
    leaf->nmd_exceptions = 0u;
    leaf->nmd_stop_position1 = first_stop;
    leaf->nmd_junction_position1 = 0u;
}

static void finish_prediction(duckvep_haplotype_stream_t *s, duckvep_haplotype_leaf_t *leaf) {
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    uint32_t tx = leaf->carriers.transcript_index;
    int raw = leaf->contributor_count && leaf->contributors[0].source.source_record;
    duckvep_prediction_reason_t relation = DUCKVEP_REASON_SUPPORTED_DOMAIN;
    /* Failed decoded-call leaves keep their edit islands in ascending CDS order;
     * conflicts are classified on the descending array before it is reversed. */
    if (!raw && !leaf->cds && leaf->edit_count) {
        sort_edits_descending(b->edits, b->edit_event_ids, leaf->edit_count, NULL);
        relation = edit_relation(b->edits, leaf->edit_count);
        for (size_t i = 0u; i < leaf->edit_count / 2u; i++) {
            duckvep_haplotype_edit_t edit = b->edits[i];
            uint64_t id = b->edit_event_ids[i];
            size_t j = leaf->edit_count - 1u - i;
            b->edits[i] = b->edits[j]; b->edit_event_ids[i] = b->edit_event_ids[j];
            b->edits[j] = edit; b->edit_event_ids[j] = id;
        }
    }
    leaf->listed_edit_count = raw && !leaf->cds ? 0u : leaf->edit_count;

    /* Contributor roles. An edit island is post-stop when it starts after the first stop codon of the
     * rebuilt protein, decided per edit from its own position in the edited CDS (the interaction block it
     * shares with earlier edits does not matter). A source is post_stop only when every island it
     * contributes is. */
    const size_t first_stop = leaf->translation.first_stop_position1;
    if (leaf->cds && !raw && first_stop) {
        for (size_t k = 0u; k < leaf->block_count; k++) {
            const duckvep_haplotype_block_t *block = &leaf->blocks[k];
            int64_t shift = (int64_t)block->alt_start0 - ((int64_t)block->cds_start - 1);
            for (size_t e = block->edit_begin; e < block->edit_begin + block->edit_count; e++) {
                const duckvep_haplotype_edit_t *edit = &b->edits[e];
                int64_t alt_start0 = (int64_t)edit->cds_start - 1 + shift;
                shift += (int64_t)edit->alt_len - (int64_t)edit->ref_len;
                if (alt_start0 < (int64_t)(first_stop * 3u)) continue;
                for (size_t i = 0u; i < leaf->contributor_count; i++)
                    if (b->contributors[i].source.event_id == b->edit_event_ids[e])
                        b->contributors[i].post_stop_edits++;
            }
        }
    }
    for (size_t i = 0u; i < leaf->contributor_count; i++) {
        duckvep_haplotype_contributor_t *c = &b->contributors[i];
        c->role = c->projection_status == DUCKVEP_CDS_EDIT_SOURCE_SHADOWED ? DUCKVEP_ROLE_SHADOWED :
            (c->edit_count || c->source_replaced) ?
                (!leaf->cds ? DUCKVEP_ROLE_UNAPPLIED :
                 c->edit_count && c->post_stop_edits == c->edit_count ? DUCKVEP_ROLE_POST_STOP :
                 DUCKVEP_ROLE_APPLIED) : DUCKVEP_ROLE_OMITTED;
    }

    /* Path eligibility, in contract order: incomplete evidence, policy, projection, domain,
     * conflict. Carrier-specific phase-domain and ploidy checks follow per key. */
    uint64_t noncoding_mask = 0u;
    int noncoding = noncoding_haplotype_mask(s, leaf, &noncoding_mask);
    if (noncoding == 1) {
        for (size_t i = 0u; i < leaf->contributor_count; i++)
            b->contributors[i].role = DUCKVEP_ROLE_APPLIED;
    }
    duckvep_prediction_status_t status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
    duckvep_prediction_reason_t reason = DUCKVEP_REASON_SUPPORTED_DOMAIN;
    if (leaf->sequence_status == DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE) {
        status = DUCKVEP_PREDICTION_INCOMPLETE_INPUT;
        reason = (leaf->evidence_flags & DUCKVEP_CARRIER_MISSING) || !(leaf->evidence_flags & DUCKVEP_CARRIER_UNPHASED)
            ? DUCKVEP_REASON_MISSING_CALL : DUCKVEP_REASON_UNPHASED_HETEROZYGOUS;
    } else if (!s->have_phase_policy || s->phase_policy != DUCKVEP_PHASE_STRICT || raw) {
        reason = DUCKVEP_REASON_NON_STRICT_PHASE_POLICY;
    } else {
        duckvep_prediction_reason_t domain = noncoding == 1 ?
            DUCKVEP_REASON_SUPPORTED_DOMAIN : transcript_domain(s, tx);
        for (size_t i = 0u; i < leaf->contributor_count && reason == DUCKVEP_REASON_SUPPORTED_DOMAIN; i++) {
            const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
            if (c->projection_status != DUCKVEP_CDS_EDIT_OK && noncoding != 1) {
                reason = DUCKVEP_REASON_PROJECTION;
                leaf->prediction_projection = c->projection_status;
            }
        }
        if (reason == DUCKVEP_REASON_SUPPORTED_DOMAIN &&
            leaf->projection_status != DUCKVEP_CDS_EDIT_OK && noncoding != 1) {
            reason = DUCKVEP_REASON_PROJECTION;
            leaf->prediction_projection = leaf->projection_status;
        }
        if (reason != DUCKVEP_REASON_SUPPORTED_DOMAIN) {
            /* Projection reasons are preserved verbatim. */
        } else if (domain != DUCKVEP_REASON_SUPPORTED_DOMAIN) reason = domain;
        else if (relation == DUCKVEP_REASON_CONTRADICTORY_EDITS) {
            status = DUCKVEP_PREDICTION_EDIT_CONFLICT; reason = relation;
        } else if (relation != DUCKVEP_REASON_SUPPORTED_DOMAIN) {
            status = DUCKVEP_PREDICTION_UNSUPPORTED_OVERLAP; reason = relation;
        } else if (leaf->sequence_status == DUCKVEP_HAPLOTYPE_REF_MISMATCH) {
            reason = DUCKVEP_REASON_REFERENCE_MISMATCH;
        } else if (leaf->sequence_status == DUCKVEP_HAPLOTYPE_EDIT_ORDER) {
            status = DUCKVEP_PREDICTION_UNSUPPORTED_OVERLAP; reason = DUCKVEP_REASON_OVERLAPPING_EDITS;
        } else if (leaf->sequence_status == DUCKVEP_HAPLOTYPE_INVALID_BASE) {
            reason = DUCKVEP_REASON_INVALID_BASE;
        } else if (leaf->sequence_status != DUCKVEP_HAPLOTYPE_OK ||
            (!leaf->cds && noncoding != 1)) {
            reason = DUCKVEP_REASON_INVALID_SEQUENCE;
        } else {
            for (size_t i = 0u; i < leaf->contributor_count && reason == DUCKVEP_REASON_SUPPORTED_DOMAIN; i++) {
                const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
                if (!literal_acgt(c->source.ref, c->source.ref_len) ||
                    !literal_acgt(c->source.alt, c->source.alt_len))
                    reason = DUCKVEP_REASON_NON_LITERAL_ALLELE;
            }
            if (reason == DUCKVEP_REASON_SUPPORTED_DOMAIN) status = DUCKVEP_PREDICTION_ELIGIBLE;
        }
    }
    leaf->path_status = status;
    leaf->path_reason = reason;
    classify_haplotype(s, leaf);
    conditional_prediction(s, leaf);
    if (noncoding == 1 && (leaf->path_status == DUCKVEP_PREDICTION_ELIGIBLE ||
                           leaf->path_status == DUCKVEP_PREDICTION_PREDICTED)) {
        leaf->path_status = DUCKVEP_PREDICTION_PREDICTED;
        leaf->haplotype_so_mask = noncoding_mask;
        leaf->nmd = DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE;
    }
    leaf->prediction_status = leaf->path_status;
    leaf->prediction_reason = leaf->path_reason;
    uint32_t id = leaf->carriers.first_call;
    for (uint32_t i = 0u; i < leaf->carriers.call_count; i++) {
        const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, id);
        if (!call) break;
        duckvep_prediction_status_t cs;
        duckvep_prediction_reason_t cr;
        duckvep_haplotype_carrier_prediction(leaf, call, &cs, &cr);
        if (cs != DUCKVEP_PREDICTION_ELIGIBLE && cs != DUCKVEP_PREDICTION_PREDICTED) {
            leaf->prediction_status = cs;
            leaf->prediction_reason = cr;
            break;
        }
        id = call->next_leaf;
    }
}

void duckvep_haplotype_carrier_prediction(const duckvep_haplotype_leaf_t *leaf,
    const duckvep_carrier_call_t *call, duckvep_prediction_status_t *status,
    duckvep_prediction_reason_t *reason) {
    *status = leaf->path_status;
    *reason = leaf->path_reason;
    if (leaf->path_status == DUCKVEP_PREDICTION_INCOMPLETE_INPUT) return;
    if (call->key.domain_split) {
        *status = DUCKVEP_PREDICTION_INCOMPLETE_INPUT;
        *reason = DUCKVEP_REASON_CROSS_PS_UNRESOLVED;
    }
    /* A lane of a complete, phased call is one haplotype whatever the call's ploidy (a haploid call has one
     * lane and nothing to phase), so its prediction is the path's. */
}

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_next(
    duckvep_haplotype_stream_t *s, duckvep_haplotype_leaf_t *out) {
    if (out) memset(out, 0, sizeof(*out));
    if (!s || !s->initialized) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    if (s->error) return s->error;
    if (!out || !s->closing) return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
    duckvep_haplotype_leaf_t leaf = {0};
    duckvep_carriers_status_t status = duckvep_carriers_next_leaf(&s->carriers, &leaf.carriers);
    if (status == DUCKVEP_CARRIERS_DONE) {
        status = duckvep_carriers_release(&s->carriers);
        if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
        s->closing = 0u;
        return DUCKVEP_HAPLOTYPE_STREAM_DONE;
    }
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    const duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    size_t count;
    status = duckvep_carriers_leaf_events(&s->carriers, leaf.carriers.id,
        b->leaf_events, b->leaf_capacity, &count);
    if (status == DUCKVEP_CARRIERS_OUTPUT_FULL)
        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL);
    if (status != DUCKVEP_CARRIERS_OK) return carrier_fail(s, status);
    uint32_t tx = leaf.carriers.transcript_index;
    int raw_records = 0, retained_call = 0;
    for (size_t i = 0u; i < count; i++) {
        const duckvep_haplotype_stored_event_t *e = find_event(s, b->leaf_events[i].event_id);
        const duckvep_haplotype_projection_t *p = e ? find_projection(s, e, tx) : NULL;
        if (!p) return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
        uint8_t evidence = b->leaf_events[i].evidence_flags;
        if (!i) raw_records = e->source.source_record;
        if (raw_records != e->source.source_record)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG);
        /* Exon-admitted genotype retention selects the upstream mutation
         * route even for shadowed, unmapped, skipped-allele or UTR sources. */
        retained_call |= p->source_exonic && (e->source.allele_index != 0u ||
            (evidence & (DUCKVEP_CARRIER_CALLED | DUCKVEP_CARRIER_REFERENCE_REPLAY)) != 0u);
        b->contributors[i] = (duckvep_haplotype_contributor_t){
            .source = e->source, .projection_status = p->status, .evidence_flags = evidence,
            .prepared = &e->prepared,
            .projected = p->status == DUCKVEP_CDS_EDIT_OK ? &p->edit : NULL};
        if (raw_records && !p->source_selected) {
            b->contributors[i].projection_status = DUCKVEP_CDS_EDIT_SOURCE_SHADOWED;
            continue;
        }
        if (raw_records && (p->status == DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED ||
                            p->status == DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED)) {
            /* A checked REF slot excluded from mutation is a nonmutating
             * reference observation. Its source ordinal, not byte equality,
             * distinguishes it from a selected ALT with the same spelling. */
            if (p->status == DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED && !e->source.allele_index) {
                b->contributors[i].projection_status = DUCKVEP_CDS_EDIT_OK;
                continue;
            }
            b->contributors[i].evidence_flags |= DUCKVEP_CARRIER_CONDITIONAL;
            continue;
        }
        if (!raw_records) leaf.evidence_flags |= evidence;
        if ((evidence & (DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_UNPHASED)) &&
            !(evidence & DUCKVEP_CARRIER_CONDITIONAL))
            leaf.sequence_status = DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE;
        if (leaf.projection_status == DUCKVEP_CDS_EDIT_OK && !p->cds_unaffected)
            leaf.projection_status = p->status;
        if (p->status == DUCKVEP_CDS_EDIT_OK &&
            (evidence & (DUCKVEP_CARRIER_CALLED | DUCKVEP_CARRIER_CONDITIONAL |
                         DUCKVEP_CARRIER_REFERENCE_REPLAY)) &&
            (!raw_records || e->source.allele_index || (evidence & DUCKVEP_CARRIER_REFERENCE_REPLAY))) {
            if (raw_records) {
                if (leaf.edit_count == b->edit_capacity)
                    return fail(s, DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL);
                b->edits[leaf.edit_count] = p->edit;
                b->edit_event_ids[leaf.edit_count++] = i;
            } else {
                size_t before = leaf.edit_count;
                duckvep_haplotype_stream_status_t added = append_differing_edits(s, p,
                    e->source.event_id, &leaf);
                if (added != DUCKVEP_HAPLOTYPE_STREAM_OK) return added;
                b->contributors[i].edit_count = (uint32_t)(leaf.edit_count - before);
            }
        }
    }
    if (raw_records && leaf.projection_status == DUCKVEP_CDS_EDIT_OK &&
        leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
        sort_edits_descending(b->edits, b->edit_event_ids, leaf.edit_count, b->contributors);
        for (size_t i = 1u; i < leaf.edit_count; i++) {
            if (b->edits[i].ref_len > b->edits[i - 1u].cds_start - b->edits[i].cds_start) {
                leaf.ordered_replacements = 1u;
                break;
            }
        }
        if (!leaf.ordered_replacements) {
            /* Disjoint full records have the same literal replay as their
             * differing islands, which also retain local codon/frame facts. */
            leaf.edit_count = 0u;
            for (size_t i = 0u; i < count; i++) {
                const duckvep_haplotype_contributor_t *c = &b->contributors[i];
                if (c->projection_status != DUCKVEP_CDS_EDIT_OK || !c->source.allele_index) continue;
                const duckvep_haplotype_stored_event_t *e = find_event(s, b->leaf_events[i].event_id);
                const duckvep_haplotype_projection_t *p = find_projection(s, e, tx);
                size_t before = leaf.edit_count;
                duckvep_haplotype_stream_status_t added = append_differing_edits(s, p,
                    c->source.event_id, &leaf);
                if (added != DUCKVEP_HAPLOTYPE_STREAM_OK) return added;
                b->contributors[i].edit_count = (uint32_t)(leaf.edit_count - before);
            }
        }
    }
    leaf.contributors = b->contributors;
    leaf.contributor_count = count;
    if (leaf.projection_status == DUCKVEP_CDS_EDIT_OK && leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
        const duckvep_sequence_pool_t *seq = s->sequences;
        uint64_t offset = seq->cds_offset[tx];
        size_t length = seq->cds_length[tx];
        if (offset > seq->cds_bytes_len || length > seq->cds_bytes_len - offset)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
        duckvep_haplotype_result_t applied;
        if (leaf.ordered_replacements) {
            leaf.sequence_status = duckvep_haplotype_compose_replacements(
                seq->cds_bytes + (size_t)offset, length, b->edits, leaf.edit_count,
                s->carriers.model->strand[tx], b->edit_event_ids, b->cds, b->cds_capacity,
                b->blocks, b->edit_capacity, &leaf.block_count, &applied);
            leaf.cds_length = applied.cds_len;
            if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
                leaf.edit_count = applied.applied_edits;
                for (size_t i = 0u; i < leaf.edit_count; i++) {
                    size_t source = (size_t)b->edit_event_ids[i];
                    if (source >= count) return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                    b->contributors[source].source_replaced = 1u;
                    b->edit_event_ids[i] = b->contributors[source].source.event_id;
                }
            }
        } else {
            sort_edits_descending(b->edits, b->edit_event_ids, leaf.edit_count, NULL);
            leaf.sequence_status = duckvep_haplotype_apply_cds_edits(seq->cds_bytes + (size_t)offset,
                length, b->edits, leaf.edit_count, s->carriers.model->strand[tx], b->cds, b->cds_capacity,
                &leaf.cds_length, &applied);
        }
        if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
            duckvep_codon_table_t table = seq->codon_table
                ? (duckvep_codon_table_t)seq->codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
            duckvep_translation_status_t translation = duckvep_translate_cds(b->cds,
                leaf.cds_length, table,
                b->protein, b->protein_capacity, &leaf.translation);
            switch (translation) {
            case DUCKVEP_TRANSLATION_OK: break;
            case DUCKVEP_TRANSLATION_BUFFER_TOO_SMALL:
                return fail(s, DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL);
            case DUCKVEP_TRANSLATION_INVALID_BASE:
                leaf.sequence_status = DUCKVEP_HAPLOTYPE_INVALID_BASE; break;
            default: leaf.sequence_status = DUCKVEP_HAPLOTYPE_INVALID_ARG; break;
            }
            if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK) {
                duckvep_haplotype_stream_status_t reference_status = prepare_reference_protein(s, tx);
                if (reference_status != DUCKVEP_HAPLOTYPE_STREAM_OK) return reference_status;
                /* Replay consumes descending coordinates; interaction discovery
                 * consumes ascending coordinates. Reverse descriptors, not bases,
                 * and borrow both sequences without rebuilding each block. */
                for (size_t i = 0u; i < leaf.edit_count / 2u; i++) {
                    duckvep_haplotype_edit_t edit = b->edits[i];
                    uint64_t id = b->edit_event_ids[i];
                    b->edits[i] = b->edits[leaf.edit_count - 1u - i];
                    b->edits[leaf.edit_count - 1u - i] = edit;
                    b->edit_event_ids[i] = b->edit_event_ids[leaf.edit_count - 1u - i];
                    b->edit_event_ids[leaf.edit_count - 1u - i] = id;
                }
                if (leaf.ordered_replacements) {
                    for (size_t i = 0u; i < leaf.block_count / 2u; i++) {
                        duckvep_haplotype_block_t block = b->blocks[i];
                        b->blocks[i] = b->blocks[leaf.block_count - 1u - i];
                        b->blocks[leaf.block_count - 1u - i] = block;
                    }
                    for (size_t i = 0u; i < leaf.block_count; i++)
                        b->blocks[i].edit_begin = leaf.edit_count -
                            (b->blocks[i].edit_begin + b->blocks[i].edit_count);
                } else if (duckvep_haplotype_partition(b->edits, leaf.edit_count,
                        b->blocks, b->edit_capacity, &leaf.block_count) != DUCKVEP_HAPLOTYPE_OK)
                    return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                size_t first_stop = leaf.translation.first_stop_position1;
                if (first_stop > leaf.cds_length / 3u)
                    return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                for (size_t i = 0u; i < leaf.block_count; i++) {
                    const duckvep_haplotype_block_t *block = &b->blocks[i];
                    size_t start0 = (size_t)block->cds_start - 1u;
                    if (start0 > length || block->ref_len > length - start0 ||
                        block->alt_start0 > leaf.cds_length ||
                        block->alt_len > leaf.cds_length - block->alt_start0)
                        return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                    if (first_stop && !leaf.ordered_replacements) {
                        int intersects;
                        if (duckvep_haplotype_block_frame_intersects(b->edits,
                                leaf.edit_count, block, (first_stop - 1u) * 3u, 3u,
                                &intersects) != DUCKVEP_HAPLOTYPE_OK)
                            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                        leaf.stop_in_displaced_frame |= (uint8_t)intersects;
                    }
                }
                leaf.blocks = b->blocks;
                leaf.edit_event_ids = b->edit_event_ids;
                leaf.reference_cds = seq->cds_bytes + (size_t)offset;
                leaf.cds = b->cds;
                leaf.protein = b->protein;
                leaf.protein_length = leaf.translation.first_stop_position1
                    ? leaf.translation.first_stop_position1 : leaf.translation.length;
                leaf.nominal_length_diff = applied.length_diff;
                leaf.flags = applied.flags;
                if (leaf.protein_length < leaf.translation.length)
                    leaf.flags |= DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED;
                leaf.reference_protein = s->reference_protein_known ? b->reference_protein : NULL;
                leaf.reference_protein_length = s->reference_protein_length;
                leaf.reference_coding_protein = b->reference_coding_protein;
                leaf.reference_coding_translation = s->reference_coding_translation;
                if (raw_records && !retained_call && leaf.reference_protein) {
                    leaf.protein = leaf.reference_protein;
                    leaf.protein_length = leaf.reference_protein_length;
                    leaf.flags &= ~(uint32_t)DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED;
                }
                if (leaf.cds_length > UINT64_MAX - s->translated_bases)
                    return fail(s, DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR);
                s->translated_bases += leaf.cds_length;
            }
        }
        if (leaf.sequence_status == DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL)
            return fail(s, DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL);
        if (leaf.sequence_status != DUCKVEP_HAPLOTYPE_OK) {
            leaf.cds_length = leaf.protein_length = 0u;
            leaf.block_count = 0u;
        }
    }
    if (raw_records) {
        size_t kept = 0u;
        for (size_t i = 0u; i < count; i++) {
            duckvep_haplotype_contributor_t c = b->contributors[i];
            if (c.evidence_flags & DUCKVEP_CARRIER_REFERENCE_REPLAY) {
                c.evidence_flags &= (uint8_t)~DUCKVEP_CARRIER_REFERENCE_REPLAY;
                if (!c.evidence_flags && !c.source_replaced && leaf.cds) continue;
                if (!c.evidence_flags) c.evidence_flags = DUCKVEP_CARRIER_CALLED;
            }
            leaf.evidence_flags |= c.evidence_flags;
            b->contributors[kept++] = c;
        }
        leaf.contributor_count = kept;
    }
    if (leaf.cds && leaf.sequence_status == DUCKVEP_HAPLOTYPE_OK &&
        (leaf.evidence_flags & DUCKVEP_CARRIER_CONDITIONAL))
        leaf.sequence_status = DUCKVEP_HAPLOTYPE_CONDITIONAL;
    finish_prediction(s, &leaf);
    s->completed_leaves++;
    *out = leaf;
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}
