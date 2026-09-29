/*
 * duckvep_lift.h — lifted-interval execution of circular sequence regions.
 *
 * The consequence kernel is linear: coordinates are unsigned 32-bit positions
 * and a span with start > end is an insertion, never a wrap. Circular topology
 * is therefore resolved before the kernel, once, into a second linear model.
 *
 * A circular region of length L that carries a wrapped object (start > end) is
 * shifted by a base B, a multiple of L at least as wide as any reference
 * window, so an event at source position p executes at p + B. Reference
 * windows and 5'/3' flanks around an event never leave the lifted interval.
 * Every object of the region is admitted at the three images p + B + k*L,
 * k in {-1, 0, +1}, which cover every relative placement of a span shorter
 * than L on the circle; a wrapped span [s, e] becomes [s, e + L]. The lifted
 * model is opened by the ordinary kernel, so candidate discovery, splice and
 * CDS projection, sequence edits, NMD, HGVS and regulation rows all consume
 * the same lifted coordinates. Regions without wrapped objects keep one image
 * at base zero, so a model without wrapped objects is unchanged.
 *
 * Kernel rows refer to lifted object indices. duckvep_lift_resolve maps them
 * to the source index and keeps one image per (event, object) pair: the image
 * nearest the event, ties broken by relative geometry only, so the choice is
 * invariant under rotation of the reference.
 */
#ifndef DUCKVEP_LIFT_H
#define DUCKVEP_LIFT_H

#include "duckvep_kernel.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DUCKVEP_LIFT_IMAGES 3u

/* Region topology, sorted ascending by chrom_id. length may be zero only for a
 * non-circular region. */
typedef struct duckvep_lift_regions {
    const uint16_t *chrom_id;
    const uint32_t *length;
    const uint8_t  *circular;
    size_t          count;
} duckvep_lift_regions_t;

typedef struct duckvep_lift {
    /* Lifted borrowed views: pass to duckvep_model_open. The byte pools of the
     * sequence view are the source pools; source arrays outlive the lift. */
    duckvep_transcript_model_t       transcripts;
    duckvep_exon_model_t             exons;
    duckvep_sequence_pool_t          sequences;
    duckvep_interval_feature_model_t interval_features;
    /* Source object and image (-1, 0, +1) of each lifted object. */
    uint32_t *transcript_source;
    int8_t   *transcript_image;
    uint32_t *feature_source;
    int8_t   *feature_image;
    /* Per region, parallel to the input regions. region_length is zero for a
     * region executed without lifting. */
    uint16_t *region_chrom_id;
    uint32_t *region_length;
    uint32_t *region_base;
    uint32_t *region_virtual_length;
    size_t    region_count;
    void     *storage; /* private arena of lifted arrays */
} duckvep_lift_t;

/* Lift wrapped circular regions. Returns DUCKVEP_ERR_UNSUPPORTED when no
 * region needs lifting; the caller then uses the source model directly. */
duckvep_status_t duckvep_lift_open(
    const duckvep_lift_regions_t           *regions,
    const duckvep_transcript_model_t       *transcripts,
    const duckvep_exon_model_t             *exons,
    const duckvep_sequence_pool_t          *seq,
    const duckvep_interval_feature_model_t *interval_features,
    duckvep_lift_t                        **out_lift,
    duckvep_error_t                        *error);
void duckvep_lift_close(duckvep_lift_t *lift);

/* Lift parameters of a region; returns 0 when the region runs unlifted. Any of
 * the outputs may be NULL. */
int duckvep_lift_region(const duckvep_lift_t *lift, uint16_t chrom_id,
                        uint32_t *length, uint32_t *base,
                        uint32_t *virtual_length);

/* Collapse the kernel rows of one lifted run to one row per (event, source
 * object). Rows are rewritten in place to source object indices and grouped by
 * variant_idx, then object kind, then source index. `order[i]` receives the
 * original position of kept row i so the caller can gather parallel streams;
 * `order` holds `count` entries. event_start1/end1 are the lifted event spans
 * indexed by row.variant_idx. */
duckvep_status_t duckvep_lift_resolve(
    const duckvep_lift_t *lift, duckvep_consequence_t *rows, size_t count,
    const uint32_t *event_start1, const uint32_t *event_end1,
    size_t *order, size_t *kept, duckvep_error_t *error);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_LIFT_H */
