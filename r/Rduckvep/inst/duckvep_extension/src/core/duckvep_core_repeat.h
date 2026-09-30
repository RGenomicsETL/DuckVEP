/* Host-neutral row logic of duckvep_repeat_alleles. The host flattens each
 * axis's list row into elements and reads the count of each into a cell. */
#ifndef DUCKVEP_CORE_REPEAT_H
#define DUCKVEP_CORE_REPEAT_H

#include "duckvep_core_cells.h"

typedef struct {
    bool present;            /* record, unit and count all non-NULL */
    const char *unit;
    uint32_t unit_length;
    duckvep_cell_t count;
} duckvep_repeat_element_t;

typedef struct {
    bool usable;             /* the list is non-NULL and its elements are {unit, count} records */
    const duckvep_repeat_element_t *elements;
    size_t count;
} duckvep_repeat_axis_t;

typedef struct {
    const char *error;       /* non-NULL: the row fails with this message */
    const char *status;      /* ok, summary_only, incomplete_input, nonintegral_count */
    long double required[2]; /* bases of the reference and alternate alleles */
} duckvep_repeat_plan_t;

extern const char duckvep_core_repeat_expected_lists[];
extern const char duckvep_core_repeat_exact_required[];
extern const char duckvep_core_repeat_cap_invalid[];

/* The max_allele_bases option: NULL cell means the option is absent (5000). */
bool duckvep_core_repeat_cap(const duckvep_cell_t *cell, long double *cap);

/* Validates and sizes one row; `error_buffer` backs a formatted capacity error. */
void duckvep_core_repeat_plan(const duckvep_repeat_axis_t axes[2], bool exact, long double cap,
    duckvep_repeat_plan_t *plan, char *error_buffer, size_t error_buffer_size);

/* Writes the axis's text (plan->required[axis] bytes plus a NUL) into out. */
size_t duckvep_core_repeat_render(const duckvep_repeat_axis_t *axis, char *out);

const char *duckvep_core_repeat_direction(const duckvep_repeat_plan_t *plan);

#endif
