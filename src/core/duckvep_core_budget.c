#include "core/duckvep_core_budget.h"

#include <stdio.h>

void duckvep_core_budget_row(unsigned index, const duckvep_budget_stats_t *stats, duckvep_core_budget_row_t *row) {
    int total = index == (unsigned)DUCKVEP_OWNER_COUNT;
    row->owner = total ? "total" : duckvep_budget_owner_name((duckvep_budget_owner_t)index);
    row->current = total ? stats->total_current : stats->current[index];
    row->high_water = total ? stats->total_high_water : stats->high_water[index];
    row->limit = stats->limit;
    row->charges = total ? stats->charges : 0;
    row->refusals = total ? stats->refusals : 0;
}

int duckvep_core_budget_set(int present, int64_t bytes, char *error, size_t error_size) {
    if (!present || bytes <= 0) {
        snprintf(error, error_size, "duckvep_native_budget_set: the budget must be a positive byte count");
        return 0;
    }
    if (!duckvep_budget_set_limit((uint64_t)bytes)) {
        duckvep_budget_stats_t stats;
        duckvep_budget_stats(&stats);
        snprintf(error, error_size,
            "duckvep_native_budget_set: capacity error: %lld bytes is below the %llu bytes already in use",
            (long long)bytes, (unsigned long long)stats.total_current);
        return 0;
    }
    return 1;
}

int duckvep_core_worker_limits(const int present[4], const int64_t value[4], char *error, size_t error_size) {
    duckvep_budget_worker_limits_t limits;
    for (unsigned i = 0u; i < 4u; i++) {
        if (!present[i]) goto invalid;
    }
    if (value[0] <= 0 || value[0] > 1024 || value[1] < 0 || value[2] < 0 || value[3] < 0) goto invalid;
    limits.max_workers = (uint32_t)value[0];
    limits.scratch_bytes = (uint64_t)value[1];
    limits.emit_bytes = (uint64_t)value[2];
    limits.idle_bytes = (uint64_t)value[3];
    (void)duckvep_budget_set_worker_limits(&limits);
    return 1;
invalid:
    snprintf(error, error_size, "duckvep_worker_limits_set: workers must be 1..1024 and byte limits non-negative");
    return 0;
}
