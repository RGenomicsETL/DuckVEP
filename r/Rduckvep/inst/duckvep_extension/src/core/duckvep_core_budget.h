/* The native allocation budget as a SQL surface, host-neutral: the rows of duckvep_native_budget and the checked
 * setters behind duckvep_native_budget_set, duckvep_worker_limits_set and duckvep_native_budget_reset_high_water.
 * The budget itself is src/kernel/src/duckvep_budget.c. Both hosts format these messages identically. */
#ifndef DUCKVEP_CORE_BUDGET_H
#define DUCKVEP_CORE_BUDGET_H

#include <stddef.h>
#include <stdint.h>

#include "kernel/src/duckvep_budget.h"

typedef struct {
    const char *owner;
    uint64_t current, high_water, limit, charges, refusals;
} duckvep_core_budget_row_t;

/* One row per owner, then "total". Charges and refusals are process-wide: reported on the total row, zero elsewhere. */
#define DUCKVEP_CORE_BUDGET_ROWS ((unsigned)DUCKVEP_OWNER_COUNT + 1u)
void duckvep_core_budget_row(unsigned index, const duckvep_budget_stats_t *stats, duckvep_core_budget_row_t *row);

/* `present` is 0 for a NULL argument. Return 1 on success, 0 with `error` set. */
int duckvep_core_budget_set(int present, int64_t bytes, char *error, size_t error_size);
int duckvep_core_worker_limits(const int present[4], const int64_t value[4], char *error, size_t error_size);

#endif
