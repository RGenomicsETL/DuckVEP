#include "duckvep_core_repeat.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

const char duckvep_core_repeat_expected_lists[] =
    "duckvep_repeat_alleles: expected lists of {unit VARCHAR, count numeric}";
const char duckvep_core_repeat_exact_required[] =
    "duckvep_repeat_alleles: sequence_exact is required";
const char duckvep_core_repeat_cap_invalid[] =
    "duckvep_repeat_alleles: max_allele_bases must be an integer from 0 through 2147483647";

bool duckvep_core_repeat_cap(const duckvep_cell_t *cell, long double *cap)
{
    bool fractional = false;
    *cap = 5000;
    if (!cell) return true;
    return duckvep_core_cell_number(cell, cap, &fractional) &&
        isfinite(*cap) && *cap >= 0 && *cap <= INT32_MAX && !fractional;
}

static bool dna(const char *data, uint32_t length)
{
    if (!length) return false;
    for (uint32_t i = 0; i < length; i++) {
        char c = data[i];
        if (!strchr("ACGTRYSWKMBDHVNacgtryswkmbdhvn", c) || !c) return false;
    }
    return true;
}

void duckvep_core_repeat_plan(const duckvep_repeat_axis_t axes[2], bool exact, long double cap,
    duckvep_repeat_plan_t *plan, char *error_buffer, size_t error_buffer_size)
{
    bool incomplete = false, fractional = false, invalid_unit = false, invalid_count = false;
    plan->error = NULL;
    plan->status = NULL;
    plan->required[0] = plan->required[1] = 0;
    for (size_t axis = 0; axis < 2; axis++) {
        const duckvep_repeat_axis_t *part = &axes[axis];
        if (!part->usable) { incomplete = true; continue; }
        for (size_t i = 0; i < part->count; i++) {
            const duckvep_repeat_element_t *element = &part->elements[i];
            long double n = 0;
            bool frac = false;
            if (!element->present) { incomplete = true; continue; }
            if (!dna(element->unit, element->unit_length)) invalid_unit = true;
            if (!duckvep_core_cell_number(&element->count, &n, &frac)) {
                incomplete = true;
                continue;
            }
            if (!isfinite(n) || n < 0) invalid_count = true;
            if (frac) fractional = true;
            plan->required[axis] += element->unit_length * n;
        }
    }
    if (invalid_unit || invalid_count) {
        plan->error = invalid_unit ? "duckvep_repeat_alleles: repeat units must contain non-empty IUPAC DNA" :
            "duckvep_repeat_alleles: repeat counts must be finite and nonnegative";
        return;
    }
    plan->status = !exact ? "summary_only" : incomplete ? "incomplete_input" :
        fractional ? "nonintegral_count" : "ok";
    if (strcmp(plan->status, "ok") == 0) {
        size_t bad = plan->required[0] > cap ? 0 : plan->required[1] > cap ? 1 : 2;
        if (bad != 2) {
            /* Print as double: MinGW's 80-bit long double does not match the
             * Windows C runtime's printf, which reads long double as double. */
            snprintf(error_buffer, error_buffer_size,
                "duckvep_repeat_alleles: %s requires %.6e bases which exceeds max_allele_bases=%.0f",
                bad ? "alternate" : "reference", (double)plan->required[bad], (double)cap);
            plan->error = error_buffer;
        }
    }
}

size_t duckvep_core_repeat_render(const duckvep_repeat_axis_t *axis, char *out)
{
    size_t written = 0;
    for (size_t i = 0; i < axis->count; i++) {
        const duckvep_repeat_element_t *element = &axis->elements[i];
        long double n = 0;
        bool frac = false;
        (void)duckvep_core_cell_number(&element->count, &n, &frac);
        for (size_t j = 0; j < (size_t)n; j++) {
            memcpy(out + written, element->unit, element->unit_length);
            written += element->unit_length;
        }
    }
    out[written] = '\0';
    return written;
}

const char *duckvep_core_repeat_direction(const duckvep_repeat_plan_t *plan)
{
    return plan->required[1] > plan->required[0] ? "GAIN" :
        plan->required[1] < plan->required[0] ? "LOSS" : "NEUTRAL";
}
