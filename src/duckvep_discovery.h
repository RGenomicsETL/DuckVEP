/* Coding-transcript discovery for one VCF record allele, shared by the scalar duckvep_coding_transcripts and the
 * fused reader duckvep_coding_calls so that both apply exactly the same normalization and overlap rules. */
#ifndef DUCKVEP_DISCOVERY_H
#define DUCKVEP_DISCOVERY_H

#include <stddef.h>
#include <stdint.h>

#include "duckdb_extension.h"
#include "duckvep_model.h"

typedef struct {
    int64_t *hits;          /* cgranges hit buffer (budget-owned, reused across calls) */
    int64_t hit_capacity;
    uint32_t cached_region;
    int32_t cached_contig;
} duckvep_discovery_t;

typedef struct {
    uint32_t *items;        /* budget-owned, grown by duckvep_discovery_append */
    size_t count, capacity;
} duckvep_u32_list_t;

typedef enum {
    DUCKVEP_DISCOVERY_OK = 0,
    DUCKVEP_DISCOVERY_BAD_SPAN,   /* position outside 1..2147483647 or an allele over 65535 bases */
    DUCKVEP_DISCOVERY_NOMEM
} duckvep_discovery_status_t;

void duckvep_discovery_init(duckvep_discovery_t *);
void duckvep_discovery_release(duckvep_discovery_t *);
/* Forget the cached region lookup (a different model was pinned). */
void duckvep_discovery_reset_model(duckvep_discovery_t *);
void duckvep_u32_list_release(duckvep_u32_list_t *);

/* Appends to `out`, in ascending transcript ordinal, the transcripts whose CDS the normalized event of (position, ref, alt)
 * overlaps; `*appended` receives how many were added (0 for alleles that are not literal replacements, unknown regions
 * and models without a complete interval index). */
duckvep_discovery_status_t duckvep_discover_coding(const duckvep_owned_model_t *model, duckvep_discovery_t *scratch,
    int64_t region, int64_t position, const uint8_t *ref, size_t ref_length, const uint8_t *alt, size_t alt_length,
    duckvep_u32_list_t *out, size_t *appended);

#endif
