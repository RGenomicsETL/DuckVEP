/* Host-neutral row logic for duckvep_so_terms, duckvep_allele_geometry and
 * duckvep_breakend_geometry. No DuckDB symbol may appear in this file or its
 * .c; the v2 host (host_v2/) adapts DuckDB vectors to these functions.
 * The v1 host keeps its own copy of this logic in src/duckvep_annotate.c and
 * src/duckvep_sql.c until a later slice moves it here; the v1-vs-v2 equality
 * test is the parity guard until then. */
#ifndef DUCKVEP_CORE_GEOMETRY_H
#define DUCKVEP_CORE_GEOMETRY_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "kernel/src/duckvep_sv.h"

/* duckvep_so_terms: one row per consequence bit. */
typedef struct {
    uint8_t bit_index;
    uint64_t consequence_mask;
    const char *consequence;
    uint8_t impact_code;
    const char *impact;
    uint8_t severity_rank;
    uint8_t evaluator_tier;
} duckvep_core_so_term_t;

size_t duckvep_core_so_term_count(void);
bool duckvep_core_so_term(size_t index, duckvep_core_so_term_t *out);

/* duckvep_allele_geometry: the 14-field result. */
typedef struct {
    uint8_t kind_code;
    bool interbase;
    uint8_t anchor_side_code;
    uint64_t raw_start0, raw_end0, feature_start0, feature_end0;
    uint64_t edit_start0, edit_end0;
    bool has_insertion_boundary0;
    uint64_t insertion_boundary0;
    uint16_t reference_difference_offset, reference_difference_length;
    uint16_t alternate_difference_offset, alternate_difference_length;
} duckvep_core_allele_geometry_t;

#define DUCKVEP_CORE_GEOMETRY_FIELD_COUNT 14
extern const char *const duckvep_core_allele_geometry_fields[DUCKVEP_CORE_GEOMETRY_FIELD_COUNT];
extern const char duckvep_core_allele_geometry_error[];

/* Returns false for an invalid input; the caller reports the message above. */
bool duckvep_core_allele_geometry(uint64_t position, const char *reference,
    size_t reference_length, const char *alternate, size_t alternate_length,
    duckvep_core_allele_geometry_t *out);

/* duckvep_breakend_geometry: a parse plus the two error messages. */
const char *duckvep_core_breakend_error(duckvep_breakend_status_t status);

#endif
