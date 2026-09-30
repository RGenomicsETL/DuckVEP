#include "core/duckvep_core_discovery.h"
#include "kernel/src/duckvep_budget.h"
#include "kernel/src/duckvep_event.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Transcript discovery for phased calls: the transcripts of a loaded model whose coding sequence the event overlaps, in
 * ascending model ordinal, from the resident transcript interval index and the model's exon arrays.
 *
 * The event is a VCF record's REF and one ALT, and its span is normalized exactly as the annotation builder normalizes it
 * (duckvep_event_prepare_small): the shared prefix and suffix are trimmed, so the anchor base of an indel never creates an
 * overlap by itself. A substitution or deletion overlaps a CDS exon when a base of its differing span lies in a CDS base of an
 * exon; an insertion (an interbase event) is placed on the base the annotation builder places it on (the left flank, or the
 * right flank at an internal coding exon's entrance, at a transcript end or after the last CDS base) and is not a CDS event
 * at the outer CDS edges, which VEP assigns to the UTR. The list is the set of (event, transcript) pairs the annotation
 * builder reports with the CDS region bit; records that overlap no such exon, the large majority of a genome, give an empty
 * list. Transcripts without a CDS are never returned. Alleles that are not literal replacements (symbolic, breakend,
 * missing or identical to REF) give an empty list. Lifted circular models are refused, as by duckvep_haplotypes. */
/* Defined in the vendored cgranges.c; the public header only exposes the by-name form. */
int64_t cr_overlap_int(const cgranges_t *cr, int32_t ctg_id, int32_t st, int32_t en, int64_t **b_, int64_t *m_b_);

void duckvep_discovery_init(duckvep_discovery_t *d) {
    d->hits = NULL; d->hit_capacity = 0; d->cached_region = UINT32_MAX; d->cached_contig = -1;
}

void duckvep_discovery_release(duckvep_discovery_t *d) {
    duckvep_budget_free(d->hits);
    duckvep_discovery_init(d);
}

void duckvep_discovery_reset_model(duckvep_discovery_t *d) {
    d->cached_region = UINT32_MAX; d->cached_contig = -1;
}

void duckvep_u32_list_release(duckvep_u32_list_t *list) {
    duckvep_budget_free(list->items);
    list->items = NULL; list->count = list->capacity = 0u;
}

static int exon_u32_compare(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

/* Is genomic base p a coding base of the exon layout of transcript t? */
static int cds_exon_base(const duckvep_owned_model_t *m, uint32_t t, int64_t p) {
    if (p < (int64_t)m->cds_starts[t] || p > (int64_t)m->cds_ends[t]) return 0;
    uint32_t offset = m->exon_offsets[t], exons = m->exon_counts[t];
    for (uint32_t e = 0u; e < exons; e++)
        if ((int64_t)m->exon_starts[offset + e] <= p && (int64_t)m->exon_ends[offset + e] >= p) return 1;
    return 0;
}

/* Does [start, end] overlap (or, with inside, lie within) a gap of at most 13 bases between consecutive exons inside the CDS interval? */
static int short_intron_overlap(const duckvep_owned_model_t *m, uint32_t t, int64_t start, int64_t end, int inside) {
    uint32_t offset = m->exon_offsets[t], exons = m->exon_counts[t];
    int64_t cs = m->cds_starts[t], ce = m->cds_ends[t];
    for (uint32_t e = 0u; e + 1u < exons; e++) {
        int64_t gap_start, gap_end;
        if (m->strands[t] >= 0) {
            gap_start = (int64_t)m->exon_ends[offset + e] + 1; gap_end = (int64_t)m->exon_starts[offset + e + 1] - 1;
        } else {
            gap_start = (int64_t)m->exon_ends[offset + e + 1] + 1; gap_end = (int64_t)m->exon_starts[offset + e] - 1;
        }
        if (gap_end < gap_start || gap_end - gap_start > 12) continue;
        int64_t lo = gap_start > cs ? gap_start : cs, hi = gap_end < ce ? gap_end : ce;
        if (inside ? (gap_start <= start && end <= gap_end && lo <= hi && lo <= start && end <= hi)
                   : (lo <= hi && lo <= end && hi >= start)) return 1;
    }
    return 0;
}

/* One normalized event against one transcript that has a CDS. */
static int event_overlaps_cds(const duckvep_owned_model_t *m, uint32_t t, const duckvep_event_t *event) {
    if (!event->interbase) {
        int64_t low = event->start1 > m->cds_starts[t] ? (int64_t)event->start1 : (int64_t)m->cds_starts[t];
        int64_t high = event->end1 < m->cds_ends[t] ? (int64_t)event->end1 : (int64_t)m->cds_ends[t];
        uint32_t offset = m->exon_offsets[t], exons = m->exon_counts[t];
        if (low <= high)
            for (uint32_t e = 0u; e < exons; e++)
                if ((int64_t)m->exon_starts[offset + e] <= high && (int64_t)m->exon_ends[offset + e] >= low) return 1;
        /* VEP treats an intron of at most 13 bases (a frameshift intron) inside the CDS as coding sequence: an event
         * that overlaps such a gap, a deletion of the whole intron included, is a CDS event. */
        return short_intron_overlap(m, t, (int64_t)event->start1, (int64_t)event->end1, 0);
    }
    int64_t left = (int64_t)event->insertion_boundary0, right = left + 1, cs = m->cds_starts[t], ce = m->cds_ends[t];
    if (short_intron_overlap(m, t, left, right, 1)) return 1;   /* an insertion inside a frameshift intron */
    if (left == ce || right == cs || left == (int64_t)m->transcript_ends[t]) return 0;   /* UTR: after_coding / before_coding */
    int64_t point = left;
    if (event->anchor_side == (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT && !cds_exon_base(m, t, left) && cds_exon_base(m, t, right))
        point = right;   /* an internal coding exon's entrance */
    return cds_exon_base(m, t, point);
}

duckvep_discovery_status_t duckvep_discover_coding(const duckvep_owned_model_t *m, duckvep_discovery_t *d,
    int64_t region, int64_t start, const uint8_t *ref, size_t ref_length, const uint8_t *alt, size_t alt_length,
    duckvep_u32_list_t *out, size_t *appended) {
    *appended = 0u;
    if (start < 1 || start > (int64_t)INT32_MAX || ref_length > UINT16_MAX || alt_length > UINT16_MAX)
        return DUCKVEP_DISCOVERY_BAD_SPAN;
    duckvep_event_t event;
    if (!ref_length || !alt_length ||
        !duckvep_event_prepare_small((uint32_t)start, ref, (uint16_t)ref_length, alt, (uint16_t)alt_length, &event))
        return DUCKVEP_DISCOVERY_OK;
    {   /* only literal replacements are events: symbolic, breakend and unspecified alleles are not */
        int literal = 1;
        for (size_t i = 0u; i < alt_length && literal; i++) literal = (alt[i] | 0x20) >= 'a' && (alt[i] | 0x20) <= 'z';
        for (size_t i = 0u; i < ref_length && literal; i++) literal = (ref[i] | 0x20) >= 'a' && (ref[i] | 0x20) <= 'z';
        if (!literal) return DUCKVEP_DISCOVERY_OK;
    }
    if (!m->interval_index_complete || region < 0 || region > UINT16_MAX) return DUCKVEP_DISCOVERY_OK;
    if ((uint32_t)region != d->cached_region) {
        char region_name[16];
        snprintf(region_name, sizeof region_name, "%u", (unsigned)region);
        d->cached_contig = cr_get_ctg(m->interval_index, region_name);
        d->cached_region = (uint32_t)region;
    }
    int64_t query_start = event.interbase ? (int64_t)event.insertion_boundary0 : (int64_t)event.start1;
    int64_t query_end = event.interbase ? (int64_t)event.insertion_boundary0 + 1 : (int64_t)event.end1;
    if (query_start < 1) query_start = 1;
    int64_t count = cr_overlap_int(m->interval_index, d->cached_contig, (int32_t)(query_start - 1), (int32_t)query_end,
        &d->hits, &d->hit_capacity);
    if (count < 0) return DUCKVEP_DISCOVERY_NOMEM;
    size_t first = out->count;
    for (int64_t i = 0; i < count; i++) {
        uint32_t transcript = (uint32_t)cr_label(m->interval_index, d->hits[i]);
        if (!m->cds_sequence_lengths[transcript] || !m->cds_starts[transcript]) continue;
        if (!event_overlaps_cds(m, transcript, &event)) continue;
        if (out->count == out->capacity) {
            size_t next = out->capacity ? out->capacity * 2u : 256u;
            uint32_t *grown = duckvep_budget_realloc(DUCKVEP_OWNER_SCRATCH, out->items, next * sizeof *out->items);
            if (!grown) return DUCKVEP_DISCOVERY_NOMEM;
            out->items = grown; out->capacity = next;
        }
        out->items[out->count++] = transcript;
    }
    if (out->count - first > 1u) qsort(out->items + first, out->count - first, sizeof *out->items, exon_u32_compare);
    *appended = out->count - first;
    return DUCKVEP_DISCOVERY_OK;
}
