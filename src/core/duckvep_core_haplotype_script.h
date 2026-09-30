/* The statements a caller runs to feed duckvep_haplotype_scan on a host whose input sink is COPY (the v2 host): the
 * caller-side normalization that the v1 host runs on a private connection (phase-domain aggregation, duplicate and
 * version checks, record ordering), wrapped around the caller's own calls query, and the COPY into the staging format.
 * The text is host-neutral; only the v2 host serves it. */
#ifndef DUCKVEP_CORE_HAPLOTYPE_SCRIPT_H
#define DUCKVEP_CORE_HAPLOTYPE_SCRIPT_H

#include "duckvep_core_haplotypes.h"
#include "duckvep_core_sql.h"

#define DUCKVEP_HAPLOTYPE_SCRIPT_MAX 6

typedef struct {
    const char *phase_policy;   /* "strict" (default) or "vep116_compat" */
    const char *input_mode;     /* "alt_events" (default) or "source_records" */
    int hgvs;                   /* -1 absent, 0 false, 1 true */
    int has_limit[DUCKVEP_HAP_LIMIT_COUNT];
    uint64_t limit[DUCKVEP_HAP_LIMIT_COUNT];
} duckvep_haplotype_script_options_t;

/* Fills out[0..*count). Returns 0 on success, 1 on allocation failure, 2 for an option combination the function
 * rejects (`error` set). */
int duckvep_core_haplotype_script(const char *query, const char *model, const char *job,
    const duckvep_haplotype_script_options_t *options, duckvep_sql_text out[DUCKVEP_HAPLOTYPE_SCRIPT_MAX],
    size_t *count, char *error, size_t error_size);

/* The SQL type text of staged input column `index` for the mode (NULL past the last column). */
const char *duckvep_hap_input_type(int source_records, unsigned index);
unsigned duckvep_hap_input_columns(int source_records);
/* The column types of the plan input (event_index, seq_region, position, reference, versions). */
const char *duckvep_hap_plan_input_type(unsigned index);

#endif
