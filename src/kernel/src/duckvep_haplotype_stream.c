#include "duckvep_haplotype_stream.h"
#include "duckvep_classify.h"
#include "duckvep_so.h"

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
        !valid_array(b->reference_coding_protein, b->reference_protein_capacity, 1u))
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
        ? duckvep_compat_vep116_source_cds_edit_build(model, s->exons, s->sequences,
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
        call->alt_index > INT32_MAX || call->phase_set.present > 1u ||
        s->buffers.events[s->current_event].source.source_record ||
        (call->policy != DUCKVEP_PHASE_STRICT && call->policy != DUCKVEP_PHASE_VEP116_COMPAT) ||
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
    if (call->policy == DUCKVEP_PHASE_VEP116_COMPAT &&
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
        if (call->policy == DUCKVEP_PHASE_VEP116_COMPAT) {
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
        size_t first = broadcast || call->policy == DUCKVEP_PHASE_VEP116_COMPAT ? 0u : declared_set;
        size_t end = broadcast ? set_count : first + 1u;
        for (size_t i = first; i < end; i++) {
            duckvep_haplotype_stream_status_t status = push_call_lane(s, tx, call, sets[i],
                assignment.lane, evidence, split);
            if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        }
    }
    /* Compaction loses the original missing slot positions, not their evidence.
     * Include the remaining lanes so an all-missing GT is never implicit REF. */
    if (call->policy == DUCKVEP_PHASE_VEP116_COMPAT && missing) {
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
        (s->have_phase_policy && s->phase_policy != DUCKVEP_PHASE_VEP116_RAW))
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
    s->phase_policy = DUCKVEP_PHASE_VEP116_RAW;
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

/* ---- coding-v1 eligibility (slice 2): status/reason only, no classifier ---- */

static int stop_codon(const uint8_t *c) {
    uint8_t a = c[0] & 0xDFu, b = c[1] & 0xDFu, d = c[2] & 0xDFu;
    return a == 'T' && ((b == 'A' && (d == 'A' || d == 'G')) || (b == 'G' && d == 'A'));
}

/* Complete table-1 CDS with phase-zero start, canonical ATG start and terminal stop,
 * no internal stop, and no curated RNA/peptide edit or recoding. Cached per transcript. */
static duckvep_prediction_reason_t transcript_domain(duckvep_haplotype_stream_t *s, uint32_t tx) {
    if (s->have_domain && s->domain_transcript == tx) return s->domain_reason;
    const duckvep_transcript_model_t *m = s->carriers.model;
    const duckvep_sequence_pool_t *seq = s->sequences;
    duckvep_prediction_reason_t r = DUCKVEP_REASON_SUPPORTED_DOMAIN;
    size_t length = seq->cds_length ? seq->cds_length[tx] : 0u;
    uint64_t offset = seq->cds_offset ? seq->cds_offset[tx] : 0u;
    const uint64_t curated = DUCKVEP_TX_SELENOCYSTEINE | DUCKVEP_TX_STOP_CODON_READTHROUGH |
        DUCKVEP_TX_RNA_EDIT | DUCKVEP_TX_AMINO_ACID_SUB;
    const uint64_t incomplete = DUCKVEP_TX_CDS_START_NF | DUCKVEP_TX_CDS_END_NF;
    uint64_t flags = m->flags ? m->flags[tx] : 0u;
    if (!length || !m->cds_start1 || !m->cds_start1[tx] || offset > seq->cds_bytes_len ||
        length > seq->cds_bytes_len - offset) r = DUCKVEP_REASON_TRANSCRIPT_NOT_CODING;
    else if (seq->codon_table && seq->codon_table[tx] != DUCKVEP_CODON_TABLE_STANDARD)
        r = DUCKVEP_REASON_NON_STANDARD_CODON_TABLE;
    else if ((flags & curated) || (seq->peptide_edit_offset &&
             seq->peptide_edit_offset[tx + 1u] != seq->peptide_edit_offset[tx]))
        r = DUCKVEP_REASON_CURATED_TRANSCRIPT;
    else if ((flags & incomplete) || length % 3u || length < 6u ||
             (m->cds_phase_offset && m->cds_phase_offset[tx]))
        r = DUCKVEP_REASON_INCOMPLETE_CDS;
    else {
        const uint8_t *c = seq->cds_bytes + (size_t)offset;
        if ((c[0] & 0xDFu) != 'A' || (c[1] & 0xDFu) != 'T' || (c[2] & 0xDFu) != 'G')
            r = DUCKVEP_REASON_NONCANONICAL_START;
        else if (!stop_codon(c + length - 3u)) r = DUCKVEP_REASON_NONCANONICAL_STOP;
        else for (size_t i = 3u; i + 3u < length; i += 3u)
            if (stop_codon(c + i)) { r = DUCKVEP_REASON_INTERNAL_STOP; break; }
    }
    s->have_domain = 1u;
    s->domain_transcript = tx;
    s->domain_reason = r;
    return r;
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

/* Whole-haplotype classifier (coding-v1 slices 3 to 5). Runs only on a path that is inside the supported
 * domain. It classifies the whole edited peptide against the uncurated reference peptide, never an edit
 * alone. The decision follows the translated sequence of the edited CDS, not the nominal net frame offset
 * of its edits, in this order:
 *   edited CDS does not begin with ATG                       -> start_lost (alone: initiation is unknown, so
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
static int base_matches(uint8_t c, char upper) { return (c & 0xDFu) == (uint8_t)upper; }

/* The edited peptide before its first stop equals the reference peptide before the terminator. */
static int same_peptide(const duckvep_haplotype_stream_t *s, const duckvep_haplotype_leaf_t *leaf,
                        size_t first_stop) {
    size_t ref_n = s->sequences->cds_length[leaf->carriers.transcript_index] / 3u - 1u;
    return ref_n == first_stop - 1u && !memcmp(leaf->reference_coding_protein, s->buffers.protein, ref_n);
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
        return;
    }
    if (!leaf->reference_coding_protein || !leaf->reference_cds || ref_length % 3u || ref_length < 6u ||
        leaf->reference_coding_translation.first_stop_position1 != ref_length / 3u) {
        leaf->path_status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        leaf->path_reason = DUCKVEP_REASON_INVALID_SEQUENCE;
        return;
    }
    if (!leaf->translation.unambiguous || !leaf->reference_coding_translation.unambiguous) {
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
    const size_t window_end = (size_t)((int64_t)alt_length - after_change);
    const size_t first_stop = leaf->translation.first_stop_position1;
    uint64_t mask;
    if (alt_length < 3u || !base_matches(leaf->cds[0], 'A') || !base_matches(leaf->cds[1], 'T') ||
        !base_matches(leaf->cds[2], 'G')) {
        mask = DUCKVEP_SO(DUCKVEP_SO_START_LOST);
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
        else if (same) mask = DUCKVEP_SO(terminal_changed ? DUCKVEP_SO_STOP_RETAINED : DUCKVEP_SO_SYNONYMOUS);
        else if (frame) mask = DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
        else if (substitutions) mask = DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
        else if (insertions) mask = DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION);
        else if (deletions) mask = DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION);
        else mask = DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
    }
    leaf->path_status = DUCKVEP_PREDICTION_PREDICTED;
    leaf->haplotype_so_mask = mask;
}

/* Runs on every leaf after sequence construction. It changes no existing field except
 * ordering the (already exposed) edit buffers of failed decoded-call leaves. */
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
    duckvep_prediction_status_t status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
    duckvep_prediction_reason_t reason = DUCKVEP_REASON_SUPPORTED_DOMAIN;
    if (leaf->sequence_status == DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE) {
        status = DUCKVEP_PREDICTION_INCOMPLETE_INPUT;
        reason = (leaf->evidence_flags & DUCKVEP_CARRIER_MISSING) || !(leaf->evidence_flags & DUCKVEP_CARRIER_UNPHASED)
            ? DUCKVEP_REASON_MISSING_CALL : DUCKVEP_REASON_UNPHASED_HETEROZYGOUS;
    } else if (!s->have_phase_policy || s->phase_policy != DUCKVEP_PHASE_STRICT || raw) {
        reason = DUCKVEP_REASON_NON_STRICT_PHASE_POLICY;
    } else {
        duckvep_prediction_reason_t domain = transcript_domain(s, tx);
        for (size_t i = 0u; i < leaf->contributor_count && reason == DUCKVEP_REASON_SUPPORTED_DOMAIN; i++) {
            const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
            if (c->projection_status != DUCKVEP_CDS_EDIT_OK) {
                reason = DUCKVEP_REASON_PROJECTION;
                leaf->prediction_projection = c->projection_status;
            }
        }
        if (reason == DUCKVEP_REASON_SUPPORTED_DOMAIN && leaf->projection_status != DUCKVEP_CDS_EDIT_OK) {
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
        } else if (leaf->sequence_status != DUCKVEP_HAPLOTYPE_OK || !leaf->cds) {
            reason = DUCKVEP_REASON_INVALID_SEQUENCE;
        } else {
            for (size_t i = 0u; i < leaf->contributor_count && reason == DUCKVEP_REASON_SUPPORTED_DOMAIN; i++) {
                const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
                if (!literal_acgt(c->source.ref, c->source.ref_len) ||
                    !literal_acgt(c->source.alt, c->source.alt_len))
                    reason = DUCKVEP_REASON_NON_LITERAL_ALLELE;
                else if (c->prepared && (c->prepared->ref_diff_length > 50u || c->prepared->alt_diff_length > 50u))
                    reason = DUCKVEP_REASON_ALLELE_OVER_50;
            }
            if (reason == DUCKVEP_REASON_SUPPORTED_DOMAIN) status = DUCKVEP_PREDICTION_ELIGIBLE;
        }
    }
    leaf->path_status = status;
    leaf->path_reason = reason;
    classify_haplotype(s, leaf);
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
    } else if ((leaf->path_status == DUCKVEP_PREDICTION_ELIGIBLE ||
                leaf->path_status == DUCKVEP_PREDICTION_PREDICTED) && call->key.ploidy != 2u) {
        *status = DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT;
        *reason = DUCKVEP_REASON_NON_DIPLOID_CALL;
    }
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
