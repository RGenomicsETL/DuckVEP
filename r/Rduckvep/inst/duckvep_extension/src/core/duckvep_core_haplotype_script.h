/* The statements a caller runs to feed duckvep_haplotype_scan on a host whose input sink is COPY (the v2 host): the
 * caller-side normalization that the v1 host runs on a private connection (phase-domain aggregation, duplicate and
 * version checks, record ordering), wrapped around the caller's own calls query, and the COPY into the staging format.
 * The text is host-neutral; only the v2 host serves it. */
#ifndef DUCKVEP_CORE_HAPLOTYPE_SCRIPT_H
#define DUCKVEP_CORE_HAPLOTYPE_SCRIPT_H

#include "duckvep_core_haplotypes.h"
#include "duckvep_core_sql.h"

/* The alt_events wrapper between the caller's query and the phase domains, shared by both hosts.
 *
 * Under the strict policy a diploid heterozygous call that is the only heterozygous or missing call of its
 * sample on a transcript needs no phase: the two haplotypes are the one that carries the allele and the one
 * that does not, whichever is written first. Such a call is read as phased in slot order. Any second
 * heterozygous or missing call of the sample on the transcript, phased or not, leaves it unresolved. */
#define DUCKVEP_HAP_ALT_CALLS(source) \
    "calls AS MATERIALIZED (SELECT *, " \
    "list_contains(list_transform(duckvep_phase_call(alleles,phase_before,{phase_set: phase_set}), " \
    "lambda a: a.phase_scope), 'phase_set') scoped FROM " source "), domains AS (SELECT transcript_index, sample_index, "
#define DUCKVEP_HAP_ALT_MIDDLE_COMPAT ") source), " DUCKVEP_HAP_ALT_CALLS("raw")
#define DUCKVEP_HAP_ALT_MIDDLE_STRICT \
    ") source), resolved AS (SELECT * REPLACE (CASE WHEN len(alleles)=2 AND alleles[1] IS NOT NULL AND " \
    "alleles[2] IS NOT NULL AND alleles[1]<>alleles[2] AND NOT coalesce(phase_before[2],false) AND " \
    "count(*) FILTER(WHERE len(list_filter(alleles,lambda a: a IS NULL))>0 OR len(list_distinct(alleles))>1) " \
    "OVER(PARTITION BY transcript_index,sample_index)=1 THEN [false,true] ELSE phase_before END AS phase_before) " \
    "FROM raw), " DUCKVEP_HAP_ALT_CALLS("resolved")

#define DUCKVEP_HAPLOTYPE_SCRIPT_MAX 6

typedef struct {
    const char *phase_policy;   /* "strict" (default) or "vep_compat" */
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
