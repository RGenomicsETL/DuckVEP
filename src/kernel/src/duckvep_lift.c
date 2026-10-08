/* duckvep_lift.c — see duckvep_lift.h. Pure C over borrowed views. */
#include "duckvep_lift.h"
#include "duckvep_budget.h"
#include "duckvep_hgvs.h"

#include <stdlib.h>
#include <string.h>

#define LIFT_WINDOW \
    ((uint64_t)DUCKVEP_HGVS_SHIFT_LIMIT + 2u * (uint64_t)UINT16_MAX + 1024u)
#define LIFT_MAX_VIRTUAL ((uint64_t)INT32_MAX - 1u)

enum {
    DVW_LIFT_ARGS = 300u,
    DVW_LIFT_REGION = 301u,
    DVW_LIFT_OOM = 302u,
    DVW_LIFT_NONE = 303u,
    DVW_LIFT_ROW = 304u
};

typedef struct lift_entry {
    uint16_t region;
    uint32_t start;
    uint32_t source;
    int8_t   image;
} lift_entry_t;

static duckvep_status_t lift_fail(duckvep_error_t *error, duckvep_status_t status,
                                  uint32_t where, const char *message) {
    if (error != NULL) {
        error->status = status;
        error->where_code = where;
        strncpy(error->message, message, sizeof error->message - 1u);
        error->message[sizeof error->message - 1u] = '\0';
    }
    return status;
}

static int lift_entry_compare(const void *left, const void *right) {
    const lift_entry_t *a = left;
    const lift_entry_t *b = right;

    if (a->region != b->region) return a->region < b->region ? -1 : 1;
    if (a->start != b->start) return a->start < b->start ? -1 : 1;
    if (a->source != b->source) return a->source < b->source ? -1 : 1;
    return (int)a->image - (int)b->image;
}

static int lift_region_index(const duckvep_lift_t *lift, uint32_t chrom_id,
                             size_t *index) {
    size_t lo = 0u;
    size_t hi = lift->region_count;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;

        if (lift->region_chrom_id[mid] < chrom_id) lo = mid + 1u;
        else hi = mid;
    }
    if (lo == lift->region_count || lift->region_chrom_id[lo] != chrom_id) {
        return 0;
    }
    *index = lo;
    return 1;
}

int duckvep_lift_region(const duckvep_lift_t *lift, uint16_t chrom_id,
                        uint32_t *length, uint32_t *base,
                        uint32_t *virtual_length) {
    size_t index;

    if (lift == NULL || !lift_region_index(lift, chrom_id, &index) ||
        lift->region_length[index] == 0u) {
        return 0;
    }
    if (length != NULL) *length = lift->region_length[index];
    if (base != NULL) *base = lift->region_base[index];
    if (virtual_length != NULL) {
        *virtual_length = lift->region_virtual_length[index];
    }
    return 1;
}

static uint32_t lift_point(uint32_t value, int wrapped, uint32_t wrap_start,
                           uint32_t length, uint64_t shift) {
    uint64_t lifted = value;

    if (wrapped && value < wrap_start) lifted += length;
    return (uint32_t)(lifted + shift);
}

/* One block holds every lifted array so close is a single free. */
typedef struct lift_arena {
    void  *block;
    size_t used;
    size_t capacity;
} lift_arena_t;

static void *lift_take(lift_arena_t *arena, size_t count, size_t size) {
    size_t bytes;
    void *pointer;

    if (count == 0u) count = 1u;
    if (size != 0u && count > (SIZE_MAX - 16u) / size) return NULL;
    bytes = (count * size + 15u) & ~(size_t)15u;
    if (arena->block == NULL) {
        /* Sizing pass: only accumulate. */
        arena->used += bytes;
        return (void *)1;
    }
    if (arena->used + bytes > arena->capacity) return NULL;
    pointer = (unsigned char *)arena->block + arena->used;
    arena->used += bytes;
    return pointer;
}

void duckvep_lift_close(duckvep_lift_t *lift) {
    if (lift == NULL) return;
    duckvep_budget_free(lift->region_chrom_id);
    duckvep_budget_free(lift->region_length);
    duckvep_budget_free(lift->region_base);
    duckvep_budget_free(lift->region_virtual_length);
    duckvep_budget_free(lift->storage);
    duckvep_budget_free(lift);
}

typedef struct lift_plan {
    lift_entry_t *entries;
    lift_entry_t *feature_entries;
    size_t        transcripts;
    size_t        features;
    size_t        exons;
    size_t        mirna;
    size_t        edits;
} lift_plan_t;

/* Allocate (or, when the arena has no block yet, size) every lifted array. */
static int lift_layout(lift_arena_t *arena, duckvep_lift_t *lift,
                       const lift_plan_t *plan, int have_mirna, int have_edits,
                       int have_cdna, int have_phase, int have_full_cdna) {
    size_t t = plan->transcripts;
    size_t e = plan->exons;
    size_t f = plan->features;
    duckvep_transcript_model_t *tx = &lift->transcripts;
    duckvep_exon_model_t *ex = &lift->exons;
    duckvep_sequence_pool_t *sq = &lift->sequences;
    duckvep_interval_feature_model_t *ft = &lift->interval_features;
    int ok = 1;

#define TAKE(target, count) \
    do { \
        void *p_ = lift_take(arena, (count), sizeof(*(target))); \
        if (p_ == NULL) ok = 0; \
        else if (arena->block != NULL) (target) = p_; \
    } while (0)
    TAKE(tx->chrom_id, t);
    TAKE(tx->start1, t);
    TAKE(tx->end1, t);
    TAKE(tx->strand, t);
    TAKE(tx->flags, t);
    TAKE(tx->exon_offset, t);
    TAKE(tx->exon_count, t);
    TAKE(tx->cds_start1, t);
    TAKE(tx->cds_end1, t);
    TAKE(ex->start1, e);
    TAKE(ex->end1, e);
    if (have_cdna) {
        TAKE(ex->cdna_start1, e);
        TAKE(ex->cdna_end1, e);
    }
    if (have_phase) {
        TAKE(ex->phase, e);
        TAKE(ex->end_phase, e);
    }
    if (have_mirna) {
        TAKE(tx->mature_mirna_offset, t + 1u);
        TAKE(tx->mature_mirna_start1, plan->mirna);
        TAKE(tx->mature_mirna_end1, plan->mirna);
    }
    TAKE(sq->cds_offset, t);
    TAKE(sq->cds_length, t);
    if (have_full_cdna) {
        TAKE(sq->cdna_offset, t);
        TAKE(sq->cdna_length, t);
    }
    TAKE(sq->codon_table, t);
    TAKE(sq->pre_cds_offset, t);
    TAKE(sq->pre_cds_length, t);
    TAKE(sq->post_cds_offset, t);
    TAKE(sq->post_cds_length, t);
    if (have_edits) {
        TAKE(sq->peptide_edit_offset, t + 1u);
        TAKE(sq->peptide_edit_position1, plan->edits);
        TAKE(sq->peptide_edit_alt, plan->edits);
    }
    TAKE(ft->chrom_id, f);
    TAKE(ft->start1, f);
    TAKE(ft->end1, f);
    TAKE(ft->kind, f);
    TAKE(lift->transcript_source, t);
    TAKE(lift->transcript_image, t);
    TAKE(lift->feature_source, f);
    TAKE(lift->feature_image, f);
#undef TAKE
    return ok;
}

duckvep_status_t duckvep_lift_open(
    const duckvep_lift_regions_t           *regions,
    const duckvep_transcript_model_t       *transcripts,
    const duckvep_exon_model_t             *exons,
    const duckvep_sequence_pool_t          *seq,
    const duckvep_interval_feature_model_t *interval_features,
    duckvep_lift_t                        **out_lift,
    duckvep_error_t                        *error) {

    duckvep_lift_t *lift = NULL;
    lift_plan_t plan;
    lift_arena_t arena;
    uint8_t *wrapped_region = NULL;
    size_t region_count, transcript_count, feature_count;
    size_t index, total, feature_total, out_exon, out_mirna, out_edit;
    int lifted_any, have_mirna, have_edits, have_cdna, have_phase;
    int have_full_cdna;
    const duckvep_lift_t *view;
    uint64_t *bases = NULL;
    duckvep_status_t status = DUCKVEP_ERR_INTERNAL;

    if (out_lift == NULL) {
        return lift_fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_LIFT_ARGS,
                         "out_lift is NULL");
    }
    *out_lift = NULL;
    if (regions == NULL || transcripts == NULL || exons == NULL) {
        return lift_fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_LIFT_ARGS,
                         "lift requires regions, transcripts and exons");
    }
    region_count = regions->count;
    transcript_count = transcripts->transcript_count;
    feature_count = interval_features != NULL
        ? interval_features->feature_count : 0u;
    if (region_count > UINT16_MAX + 1u ||
        (region_count != 0u && (regions->chrom_id == NULL ||
                                regions->length == NULL ||
                                regions->circular == NULL))) {
        return lift_fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_LIFT_ARGS,
                         "lift region view is incomplete");
    }
    for (index = 1u; index < region_count; index++) {
        if (regions->chrom_id[index] <= regions->chrom_id[index - 1u]) {
            return lift_fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_LIFT_REGION,
                             "lift regions are not sorted and unique");
        }
    }
    lift = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, 1u, sizeof *lift);
    wrapped_region = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, region_count + 1u, 1u);
    bases = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, region_count + 1u, sizeof *bases);
    if (lift == NULL || wrapped_region == NULL || bases == NULL) {
        status = lift_fail(error, DUCKVEP_ERR_INTERNAL, DVW_LIFT_OOM,
                           "out of memory lifting circular regions");
        goto done;
    }
    lift->region_count = region_count;
    lift->region_chrom_id = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, region_count + 1u, sizeof(uint16_t));
    lift->region_length = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, region_count + 1u, sizeof(uint32_t));
    lift->region_base = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, region_count + 1u, sizeof(uint32_t));
    lift->region_virtual_length = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, region_count + 1u, sizeof(uint32_t));
    if (lift->region_chrom_id == NULL || lift->region_length == NULL ||
        lift->region_base == NULL || lift->region_virtual_length == NULL) {
        status = lift_fail(error, DUCKVEP_ERR_INTERNAL, DVW_LIFT_OOM,
                           "out of memory lifting circular regions");
        goto done;
    }
    /* The arena is one allocation, so the region arrays above are allocated
     * separately and freed with the arena owner below. */
    for (index = 0u; index < region_count; index++) {
        lift->region_chrom_id[index] = regions->chrom_id[index];
    }
    for (index = 0u; index < transcript_count; index++) {
        size_t region, exon;

        if (!lift_region_index(lift, transcripts->chrom_id[index], &region)) {
            continue;
        }
        if (transcripts->start1[index] > transcripts->end1[index]) {
            wrapped_region[region] = 1u;
        }
        for (exon = transcripts->exon_offset[index];
             exon < (size_t)transcripts->exon_offset[index] +
                    transcripts->exon_count[index] &&
             exon < exons->exon_count; exon++) {
            if (exons->start1[exon] > exons->end1[exon]) {
                wrapped_region[region] = 1u;
            }
        }
    }
    for (index = 0u; index < feature_count; index++) {
        size_t region;

        if (lift_region_index(lift, interval_features->chrom_id[index], &region) &&
            interval_features->start1[index] > interval_features->end1[index]) {
            wrapped_region[region] = 1u;
        }
    }
    lifted_any = 0;
    for (index = 0u; index < region_count; index++) {
        uint64_t length, multiple, base, virtual_length;

        if (!wrapped_region[index]) continue;
        if (!regions->circular[index] || regions->length[index] == 0u) {
            status = lift_fail(error, DUCKVEP_ERR_MODEL_INVALID, DVW_LIFT_REGION,
                               "wrapped coordinates require a circular region with a sequence length");
            goto done;
        }
        length = regions->length[index];
        multiple = (LIFT_WINDOW + length - 1u) / length;
        base = multiple * length;
        virtual_length = 2u * base + 3u * length;
        if (virtual_length >= LIFT_MAX_VIRTUAL) {
            status = lift_fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_LIFT_REGION,
                               "circular region is too long for lifted interval execution");
            goto done;
        }
        lift->region_length[index] = (uint32_t)length;
        lift->region_base[index] = (uint32_t)base;
        lift->region_virtual_length[index] = (uint32_t)virtual_length;
        bases[index] = base;
        lifted_any = 1;
    }
    if (!lifted_any) {
        status = lift_fail(error, DUCKVEP_ERR_UNSUPPORTED, DVW_LIFT_NONE,
                           "no wrapped circular region to lift");
        goto done;
    }
    view = lift;

    /* Admission: one entry per object image, ordered by (region, lifted start,
     * source, image) so the lifted model satisfies the kernel sort contract and
     * an unlifted region keeps its source order. */
    memset(&plan, 0, sizeof plan);
    total = 0u;
    for (index = 0u; index < transcript_count; index++) {
        size_t region;

        total += (lift_region_index(view, transcripts->chrom_id[index], &region) &&
                  lift->region_length[region] != 0u) ? DUCKVEP_LIFT_IMAGES : 1u;
    }
    feature_total = 0u;
    for (index = 0u; index < feature_count; index++) {
        size_t region;

        feature_total += (lift_region_index(view, interval_features->chrom_id[index],
                                            &region) &&
                          lift->region_length[region] != 0u)
            ? DUCKVEP_LIFT_IMAGES : 1u;
    }
    plan.entries = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, total + 1u, sizeof *plan.entries);
    plan.feature_entries = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, feature_total + 1u, sizeof *plan.feature_entries);
    if (plan.entries == NULL || plan.feature_entries == NULL) {
        status = lift_fail(error, DUCKVEP_ERR_INTERNAL, DVW_LIFT_OOM,
                           "out of memory lifting circular regions");
        goto plan_done;
    }
    total = 0u;
    for (index = 0u; index < transcript_count; index++) {
        size_t region = 0u;
        int image, first, images;

        (void)lift_region_index(view, transcripts->chrom_id[index], &region);
        images = lift->region_length[region] != 0u ? (int)DUCKVEP_LIFT_IMAGES : 1;
        first = images == (int)DUCKVEP_LIFT_IMAGES ? -1 : 0;
        for (image = first; image < first + images; image++) {
            uint64_t shift = lift->region_length[region] == 0u ? 0u
                : (uint64_t)((int64_t)bases[region] +
                             (int64_t)image * (int64_t)lift->region_length[region]);

            plan.entries[total].region = transcripts->chrom_id[index];
            plan.entries[total].source = (uint32_t)index;
            plan.entries[total].image = (int8_t)image;
            plan.entries[total].start = (uint32_t)(transcripts->start1[index] + shift);
            total++;
        }
    }
    qsort(plan.entries, total, sizeof *plan.entries, lift_entry_compare);
    feature_total = 0u;
    for (index = 0u; index < feature_count; index++) {
        size_t region = 0u;
        int image, first, images;

        (void)lift_region_index(view, interval_features->chrom_id[index], &region);
        images = lift->region_length[region] != 0u ? (int)DUCKVEP_LIFT_IMAGES : 1;
        first = images == (int)DUCKVEP_LIFT_IMAGES ? -1 : 0;
        for (image = first; image < first + images; image++) {
            uint64_t shift = lift->region_length[region] == 0u ? 0u
                : (uint64_t)((int64_t)bases[region] +
                             (int64_t)image * (int64_t)lift->region_length[region]);

            plan.feature_entries[feature_total].region =
                interval_features->chrom_id[index];
            plan.feature_entries[feature_total].source = (uint32_t)index;
            plan.feature_entries[feature_total].image = (int8_t)image;
            plan.feature_entries[feature_total].start =
                (uint32_t)(interval_features->start1[index] + shift);
            feature_total++;
        }
    }
    qsort(plan.feature_entries, feature_total, sizeof *plan.feature_entries,
          lift_entry_compare);
    plan.transcripts = total;
    plan.features = feature_total;
    have_mirna = transcripts->mature_mirna_offset != NULL;
    have_edits = seq != NULL && seq->peptide_edit_offset != NULL;
    have_cdna = exons->cdna_start1 != NULL;
    have_phase = exons->phase != NULL;
    have_full_cdna = seq != NULL && seq->cdna_provided != 0 &&
        seq->cdna_offset != NULL && seq->cdna_length != NULL;
    for (index = 0u; index < total; index++) {
        uint32_t t = plan.entries[index].source;

        plan.exons += transcripts->exon_count[t];
        if (have_mirna) {
            plan.mirna += transcripts->mature_mirna_offset[t + 1u] -
                          transcripts->mature_mirna_offset[t];
        }
        if (have_edits) {
            plan.edits += seq->peptide_edit_offset[t + 1u] -
                          seq->peptide_edit_offset[t];
        }
    }
    if (seq == NULL) {
        status = lift_fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_LIFT_ARGS,
                           "lift requires the sequence pool");
        goto plan_done;
    }

    memset(&arena, 0, sizeof arena);
    (void)lift_layout(&arena, lift, &plan, have_mirna, have_edits, have_cdna,
                      have_phase, have_full_cdna);
    arena.capacity = arena.used;
    arena.used = 0u;
    arena.block = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, 1u, arena.capacity != 0u ? arena.capacity : 16u);
    if (arena.block == NULL) {
        status = lift_fail(error, DUCKVEP_ERR_INTERNAL, DVW_LIFT_OOM,
                           "out of memory lifting circular regions");
        goto plan_done;
    }
    lift->storage = arena.block;
    if (!lift_layout(&arena, lift, &plan, have_mirna, have_edits, have_cdna,
                     have_phase, have_full_cdna)) {
        status = lift_fail(error, DUCKVEP_ERR_INTERNAL, DVW_LIFT_OOM,
                           "lifted arena layout failed");
        goto plan_done;
    }

    {
        /* The arena arrays were published as const views; fill through
         * mutable aliases of the same storage. */
        uint16_t *tx_chrom = (uint16_t *)lift->transcripts.chrom_id;
        uint32_t *tx_start = (uint32_t *)lift->transcripts.start1;
        uint32_t *tx_end = (uint32_t *)lift->transcripts.end1;
        int8_t *tx_strand = (int8_t *)lift->transcripts.strand;
        uint64_t *tx_flags = (uint64_t *)lift->transcripts.flags;
        uint32_t *tx_exon_offset = (uint32_t *)lift->transcripts.exon_offset;
        uint16_t *tx_exon_count = (uint16_t *)lift->transcripts.exon_count;
        uint32_t *tx_cds_start = (uint32_t *)lift->transcripts.cds_start1;
        uint32_t *tx_cds_end = (uint32_t *)lift->transcripts.cds_end1;
        uint32_t *ex_start = (uint32_t *)lift->exons.start1;
        uint32_t *ex_end = (uint32_t *)lift->exons.end1;
        uint32_t *ex_cdna_start = (uint32_t *)lift->exons.cdna_start1;
        uint32_t *ex_cdna_end = (uint32_t *)lift->exons.cdna_end1;
        int8_t *ex_phase = (int8_t *)lift->exons.phase;
        int8_t *ex_end_phase = (int8_t *)lift->exons.end_phase;
        uint32_t *mirna_offset = (uint32_t *)lift->transcripts.mature_mirna_offset;
        uint32_t *mirna_start = (uint32_t *)lift->transcripts.mature_mirna_start1;
        uint32_t *mirna_end = (uint32_t *)lift->transcripts.mature_mirna_end1;
        uint64_t *cds_offset = (uint64_t *)lift->sequences.cds_offset;
        uint32_t *cds_length = (uint32_t *)lift->sequences.cds_length;
        uint64_t *full_cdna_offset = (uint64_t *)lift->sequences.cdna_offset;
        uint32_t *full_cdna_length = (uint32_t *)lift->sequences.cdna_length;
        uint8_t *codon_table = (uint8_t *)lift->sequences.codon_table;
        uint64_t *pre_offset = (uint64_t *)lift->sequences.pre_cds_offset;
        uint32_t *pre_length = (uint32_t *)lift->sequences.pre_cds_length;
        uint64_t *post_offset = (uint64_t *)lift->sequences.post_cds_offset;
        uint32_t *post_length = (uint32_t *)lift->sequences.post_cds_length;
        uint32_t *edit_offset = (uint32_t *)lift->sequences.peptide_edit_offset;
        uint32_t *edit_position = (uint32_t *)lift->sequences.peptide_edit_position1;
        uint8_t *edit_alt = (uint8_t *)lift->sequences.peptide_edit_alt;
        uint16_t *ft_chrom = (uint16_t *)lift->interval_features.chrom_id;
        uint32_t *ft_start = (uint32_t *)lift->interval_features.start1;
        uint32_t *ft_end = (uint32_t *)lift->interval_features.end1;
        uint8_t *ft_kind = (uint8_t *)lift->interval_features.kind;
        int have_flanks = seq->pre_cds_offset != NULL && seq->pre_cds_length != NULL &&
                          seq->post_cds_offset != NULL && seq->post_cds_length != NULL;

        out_exon = out_mirna = out_edit = 0u;
        if (have_mirna) mirna_offset[0] = 0u;
        if (have_edits) edit_offset[0] = 0u;
        for (index = 0u; index < total; index++) {
            uint32_t t = plan.entries[index].source;
            size_t region = 0u, exon, segment, edit;
            uint32_t length, wrap_start;
            uint64_t shift;
            int wrapped;

            (void)lift_region_index(view, transcripts->chrom_id[t], &region);
            length = lift->region_length[region];
            shift = length == 0u ? 0u
                : (uint64_t)((int64_t)bases[region] +
                             (int64_t)plan.entries[index].image * (int64_t)length);
            wrapped = length != 0u && transcripts->start1[t] > transcripts->end1[t];
            wrap_start = transcripts->start1[t];
            tx_chrom[index] = transcripts->chrom_id[t];
            tx_start[index] = plan.entries[index].start;
            tx_end[index] = (uint32_t)((uint64_t)transcripts->end1[t] +
                                       (wrapped ? length : 0u) + shift);
            tx_strand[index] = transcripts->strand[t];
            tx_flags[index] = transcripts->flags[t];
            if (transcripts->cds_start1[t] != 0u) {
                tx_cds_start[index] = lift_point(transcripts->cds_start1[t],
                                                 wrapped, wrap_start, length, shift);
                tx_cds_end[index] = lift_point(transcripts->cds_end1[t], wrapped,
                                               wrap_start, length, shift);
            }
            tx_exon_offset[index] = (uint32_t)out_exon;
            tx_exon_count[index] = transcripts->exon_count[t];
            for (exon = transcripts->exon_offset[t];
                 exon < (size_t)transcripts->exon_offset[t] +
                        transcripts->exon_count[t]; exon++, out_exon++) {
                ex_start[out_exon] = lift_point(exons->start1[exon], wrapped,
                                                wrap_start, length, shift);
                ex_end[out_exon] = lift_point(exons->end1[exon], wrapped,
                                              wrap_start, length, shift);
                if (have_cdna) {
                    ex_cdna_start[out_exon] = exons->cdna_start1[exon];
                    ex_cdna_end[out_exon] = exons->cdna_end1[exon];
                }
                if (have_phase) {
                    ex_phase[out_exon] = exons->phase[exon];
                    ex_end_phase[out_exon] = exons->end_phase[exon];
                }
            }
            if (have_mirna) {
                for (segment = transcripts->mature_mirna_offset[t];
                     segment < transcripts->mature_mirna_offset[t + 1u];
                     segment++, out_mirna++) {
                    mirna_start[out_mirna] = lift_point(
                        transcripts->mature_mirna_start1[segment], wrapped,
                        wrap_start, length, shift);
                    mirna_end[out_mirna] = lift_point(
                        transcripts->mature_mirna_end1[segment], wrapped,
                        wrap_start, length, shift);
                }
                mirna_offset[index + 1u] = (uint32_t)out_mirna;
            }
            cds_offset[index] = seq->cds_offset[t];
            cds_length[index] = seq->cds_length[t];
            if (have_full_cdna) {
                full_cdna_offset[index] = seq->cdna_offset[t];
                full_cdna_length[index] = seq->cdna_length[t];
            }
            codon_table[index] = seq->codon_table[t];
            if (have_flanks) {
                pre_offset[index] = seq->pre_cds_offset[t];
                pre_length[index] = seq->pre_cds_length[t];
                post_offset[index] = seq->post_cds_offset[t];
                post_length[index] = seq->post_cds_length[t];
            }
            if (have_edits) {
                for (edit = seq->peptide_edit_offset[t];
                     edit < seq->peptide_edit_offset[t + 1u]; edit++, out_edit++) {
                    edit_position[out_edit] = seq->peptide_edit_position1[edit];
                    edit_alt[out_edit] = seq->peptide_edit_alt[edit];
                }
                edit_offset[index + 1u] = (uint32_t)out_edit;
            }
            lift->transcript_source[index] = t;
            lift->transcript_image[index] = plan.entries[index].image;
        }
        for (index = 0u; index < feature_total; index++) {
            uint32_t f = plan.feature_entries[index].source;
            size_t region = 0u;
            uint32_t length;
            uint64_t shift;
            int wrapped;

            (void)lift_region_index(view, interval_features->chrom_id[f], &region);
            length = lift->region_length[region];
            shift = length == 0u ? 0u
                : (uint64_t)((int64_t)bases[region] +
                             (int64_t)plan.feature_entries[index].image * (int64_t)length);
            wrapped = length != 0u &&
                      interval_features->start1[f] > interval_features->end1[f];
            ft_chrom[index] = interval_features->chrom_id[f];
            ft_start[index] = plan.feature_entries[index].start;
            ft_end[index] = (uint32_t)((uint64_t)interval_features->end1[f] +
                                       (wrapped ? length : 0u) + shift);
            ft_kind[index] = interval_features->kind[f];
            lift->feature_source[index] = f;
            lift->feature_image[index] = plan.feature_entries[index].image;
        }
    }
    lift->transcripts.transcript_count = total;
    lift->transcripts.mature_mirna_count = plan.mirna;
    lift->exons.exon_count = plan.exons;
    lift->sequences.transcript_count = total;
    lift->sequences.cds_bytes = seq->cds_bytes;
    lift->sequences.cds_bytes_len = seq->cds_bytes_len;
    lift->sequences.cdna_bytes = seq->cdna_bytes;
    lift->sequences.cdna_bytes_len = seq->cdna_bytes_len;
    lift->sequences.cdna_provided = (uint8_t)have_full_cdna;
    lift->sequences.peptide_edit_count = plan.edits;
    lift->sequences.flank_bytes = seq->flank_bytes;
    lift->sequences.flank_bytes_len = seq->flank_bytes_len;
    lift->sequences.flanks_complete = seq->flanks_complete;
    if (!(seq->pre_cds_offset != NULL && seq->pre_cds_length != NULL &&
          seq->post_cds_offset != NULL && seq->post_cds_length != NULL)) {
        lift->sequences.pre_cds_offset = NULL;
        lift->sequences.pre_cds_length = NULL;
        lift->sequences.post_cds_offset = NULL;
        lift->sequences.post_cds_length = NULL;
    }
    lift->interval_features.feature_count = feature_total;
    *out_lift = lift;
    lift = NULL;
    status = DUCKVEP_OK;
plan_done:
    duckvep_budget_free(plan.entries);
    duckvep_budget_free(plan.feature_entries);
done:
    duckvep_budget_free(wrapped_region);
    duckvep_budget_free(bases);
    if (lift != NULL) {
        duckvep_budget_free(lift->region_chrom_id);
        duckvep_budget_free(lift->region_length);
        duckvep_budget_free(lift->region_base);
        duckvep_budget_free(lift->region_virtual_length);
        duckvep_budget_free(lift->storage);
        duckvep_budget_free(lift);
    }
    return status;
}

typedef struct lift_row {
    uint32_t variant;
    uint32_t object;
    uint64_t gap;
    uint64_t overlap_rank;
    uint8_t  kind;
    uint8_t  after;
    uint8_t  image_rank;
    size_t   row;
} lift_row_t;

static int lift_row_compare(const void *left, const void *right) {
    const lift_row_t *a = left;
    const lift_row_t *b = right;

    if (a->variant != b->variant) return a->variant < b->variant ? -1 : 1;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    if (a->object != b->object) return a->object < b->object ? -1 : 1;
    if (a->gap != b->gap) return a->gap < b->gap ? -1 : 1;
    if (a->overlap_rank != b->overlap_rank) {
        return a->overlap_rank < b->overlap_rank ? -1 : 1;
    }
    if (a->after != b->after) return a->after < b->after ? -1 : 1;
    if (a->image_rank != b->image_rank) {
        return a->image_rank < b->image_rank ? -1 : 1;
    }
    return a->row < b->row ? -1 : a->row > b->row;
}

duckvep_status_t duckvep_lift_resolve(
    const duckvep_lift_t *lift, duckvep_consequence_t *rows, size_t count,
    const uint32_t *event_start1, const uint32_t *event_end1,
    size_t *order, size_t *kept_out, duckvep_error_t *error) {

    lift_row_t *keys;
    duckvep_consequence_t *copy;
    size_t index, kept;

    if (kept_out != NULL) *kept_out = 0u;
    if (lift == NULL || (count != 0u && (rows == NULL || event_start1 == NULL ||
                                         event_end1 == NULL || order == NULL)) ||
        kept_out == NULL) {
        return lift_fail(error, DUCKVEP_ERR_INVALID_ARG, DVW_LIFT_ROW,
                         "lift resolve requires rows, spans and an order buffer");
    }
    if (count == 0u) return DUCKVEP_OK;
    keys = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, count * sizeof *keys);
    copy = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, count * sizeof *copy);
    if (keys == NULL || copy == NULL) {
        duckvep_budget_free(keys);
        duckvep_budget_free(copy);
        return lift_fail(error, DUCKVEP_ERR_INTERNAL, DVW_LIFT_OOM,
                         "out of memory resolving lifted rows");
    }
    for (index = 0u; index < count; index++) {
        const duckvep_consequence_t *row = &rows[index];
        lift_row_t *key = &keys[index];
        uint64_t a = event_start1[row->variant_idx];
        uint64_t b = event_end1[row->variant_idx];
        uint64_t start, end, low, high;
        int8_t image;
        int is_transcript = row->overlap_object_kind ==
                            (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT;
        uint32_t lifted = is_transcript ? row->tx_idx : row->interval_feature_idx;

        if (lifted >= (is_transcript ? lift->transcripts.transcript_count
                                     : lift->interval_features.feature_count)) {
            duckvep_budget_free(keys);
            duckvep_budget_free(copy);
            return lift_fail(error, DUCKVEP_ERR_OUT_OF_RANGE, DVW_LIFT_ROW,
                             "lifted row references an unknown object");
        }
        key->variant = row->variant_idx;
        key->kind = row->overlap_object_kind;
        key->row = index;
        if (is_transcript) {
            start = lift->transcripts.start1[lifted];
            end = lift->transcripts.end1[lifted];
            image = lift->transcript_image[lifted];
            key->object = lift->transcript_source[lifted];
        } else {
            start = lift->interval_features.start1[lifted];
            end = lift->interval_features.end1[lifted];
            image = lift->feature_image[lifted];
            key->object = lift->feature_source[lifted];
        }
        key->gap = start > b ? start - b : (a > end ? a - end : 0u);
        low = a > start ? a : start;
        high = b < end ? b : end;
        key->overlap_rank = (key->gap == 0u && high >= low)
            ? UINT64_MAX - (high - low + 1u) : UINT64_MAX;
        key->after = start >= a ? 0u : 1u;
        key->image_rank = image == 0 ? 0u : (image < 0 ? 1u : 2u);
    }
    qsort(keys, count, sizeof *keys, lift_row_compare);
    memcpy(copy, rows, count * sizeof *rows);
    kept = 0u;
    for (index = 0u; index < count; index++) {
        if (index != 0u && keys[index].variant == keys[index - 1u].variant &&
            keys[index].kind == keys[index - 1u].kind &&
            keys[index].object == keys[index - 1u].object) {
            continue;
        }
        keys[kept++] = keys[index];
    }
    for (index = 0u; index < kept; index++) {
        duckvep_consequence_t *row = &rows[index];

        order[index] = keys[index].row;
        *row = copy[keys[index].row];
        if (row->overlap_object_kind ==
            (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT) {
            row->tx_idx = lift->transcript_source[row->tx_idx];
        } else {
            row->interval_feature_idx = lift->feature_source[row->interval_feature_idx];
        }
    }
    duckvep_budget_free(keys);
    duckvep_budget_free(copy);
    *kept_out = kept;
    return DUCKVEP_OK;
}
