#include "duckvep_core_haplotype_script.h"
#include "duckvep_core_phase.h"

#include <stdio.h>
#include <string.h>

/* The alt_events wrapper of the v1 host, verbatim in meaning: type the caller's columns, derive each call's phase
 * scopes with duckvep_phase_call, aggregate the phase domains per (transcript, sample) and attach the identity facts
 * (duplicate copies, event versions, ploidies) that the native scan validates. */
static const char alt_prefix[] =
    "WITH raw AS MATERIALIZED (SELECT event_index::UBIGINT event_index, seq_region::UINTEGER seq_region, "
    "position::UBIGINT AS position, reference::VARCHAR AS reference, alternate::VARCHAR AS alternate, "
    "alt_index::UINTEGER alt_index, transcript_index::UINTEGER transcript_index, sample_index::UINTEGER sample_index, "
    "alleles::INTEGER[] alleles, phase_before::BOOLEAN[] phase_before, phase_set::BIGINT phase_set FROM (";
static const char alt_middle[] =
    ") source), calls AS MATERIALIZED (SELECT *, "
    "list_contains(list_transform(duckvep_phase_call(alleles,phase_before,{phase_set: phase_set}), "
    "lambda a: a.phase_scope), 'phase_set') scoped FROM raw), domains AS (SELECT transcript_index, sample_index, ";
static const char alt_domain_strict[] =
    "coalesce(list(DISTINCT phase_set ORDER BY phase_set NULLS FIRST) FILTER(WHERE scoped), [NULL]::BIGINT[]) AS domain_sets ";
static const char alt_domain_compat[] = "[NULL]::BIGINT[] AS domain_sets ";
static const char alt_suffix[] =
    ",count(DISTINCT len(alleles)) ploidies FROM calls GROUP BY transcript_index,sample_index), "
    "event_versions AS (SELECT event_index, "
    "count(DISTINCT (seq_region,position,reference,alternate,alt_index)) versions FROM raw GROUP BY event_index) "
    "SELECT c.* EXCLUDE(scoped), d.domain_sets, "
    "count(*) OVER(PARTITION BY c.event_index,transcript_index,sample_index) copies, v.versions, d.ploidies "
    "FROM calls c LEFT JOIN domains d USING(transcript_index,sample_index) "
    "LEFT JOIN event_versions v USING(event_index) "
    "ORDER BY seq_region,position,event_index,transcript_index,sample_index";

static const char raw_columns[] =
    "SELECT event_index::UBIGINT event_index, "
    "seq_region::UINTEGER seq_region, position::UBIGINT AS position, reference::VARCHAR AS reference, "
    "alternates::VARCHAR[] alternates, transcript_index::UINTEGER transcript_index, "
    "sample_index::UINTEGER sample_index, gt::VARCHAR gt FROM (";

static const char plan_query_head[] =
    "SELECT event_index, first(seq_region), first(position), first(reference), "
    "count(DISTINCT (seq_region,position,reference,alternates)) versions FROM ";
static const char plan_query_tail[] = " GROUP BY event_index ORDER BY 2,3,1";

/* The ordering, raw-GT parsing and selection of the source_records mode of the v1 host, reading the raw relation and
 * the record plan (_duckvep_haplotype_plan) instead of two private TEMP tables. */
static const char final_head[] =
    "WITH ordered AS (SELECT event_index, CASE WHEN buffer_id=0 THEN source_ordinal ELSE "
    "source_ordinal-ordinal+_duckvep_record_order(count(*) OVER(PARTITION BY buffer_id)::UBIGINT,ordinal) "
    "END replay_order FROM _duckvep_haplotype_plan(";
static const char final_mid1[] =
    ")), genotypes AS (SELECT event_index,sample_index,count(DISTINCT gt) gt_versions FROM ";
static const char final_mid2[] =
    " GROUP BY event_index,sample_index), "
    "calls AS MATERIALIZED (SELECT *, count(*) OVER(PARTITION BY event_index,transcript_index,sample_index) copies, "
    "CASE WHEN alternates IS NULL OR len(alternates)>2147483647 OR "
    "len(list_filter(alternates,lambda a: a IS NULL OR len(a)=0 OR len(a)>65535))>0 "
    "THEN error('duckvep_haplotypes: invalid source ALT list') ELSE len(alternates) END alt_count FROM ";
static const char final_tail[] =
    " JOIN ordered USING(event_index)), "
    "parsed AS MATERIALIZED (SELECT *, _duckvep_raw_gt(gt,alt_count::UINTEGER) raw_gt FROM calls), "
    "selected AS (SELECT *, max(replay_order) FILTER(WHERE raw_gt.disposition=3) OVER "
    "(PARTITION BY seq_region,position,reference,alternates,transcript_index) selected_order FROM parsed) "
    "SELECT c.event_index,seq_region,position,reference, "
    "CASE WHEN a.i=0 THEN reference WHEN a.i>alt_count THEN '' ELSE alternates[a.i] END alternate, "
    "(CASE WHEN a.i>alt_count THEN 4294967295 ELSE a.i END)::UINTEGER alt_index, "
    "transcript_index,c.sample_index,raw_gt,NULL::BOOLEAN[] phase_before,NULL::BIGINT phase_set, "
    "alt_count,copies,1::BIGINT versions,gt_versions,replay_order, "
    "(selected_order IS NULL OR replay_order=selected_order) source_selected "
    "FROM selected c LEFT JOIN genotypes g USING(event_index,sample_index),range(0,alt_count+2) a(i) "
    "ORDER BY seq_region,position,event_index,alt_index,transcript_index,sample_index";

static bool append_copy_options(duckvep_sql_text *sql, const char *model, const char *job, const char *stage,
    const duckvep_haplotype_script_options_t *o) {
    char number[48];
    bool ok = duckvep_sql_append(sql, "\n) TO 'duckvep_stage' (FORMAT duckvep_stage, JOB ") && duckvep_sql_literal(sql, job) &&
        duckvep_sql_append(sql, ", MODEL ") && duckvep_sql_literal(sql, model);
    if (ok && stage) ok = duckvep_sql_append(sql, ", STAGE ") && duckvep_sql_literal(sql, stage);
    if (ok && !stage) {
        if (o->phase_policy) ok = duckvep_sql_append(sql, ", PHASE_POLICY ") && duckvep_sql_literal(sql, o->phase_policy);
        if (ok && o->input_mode) ok = duckvep_sql_append(sql, ", INPUT_MODE ") && duckvep_sql_literal(sql, o->input_mode);
        if (ok && o->hgvs >= 0) ok = duckvep_sql_append(sql, o->hgvs ? ", HGVS true" : ", HGVS false");
        for (unsigned i = 0u; ok && i < DUCKVEP_HAP_LIMIT_COUNT; i++) {
            if (!o->has_limit[i]) continue;
            snprintf(number, sizeof number, " %llu", (unsigned long long)o->limit[i]);
            ok = duckvep_sql_append(sql, ", ") && duckvep_sql_append(sql, duckvep_hap_limit_names[i]) &&
                duckvep_sql_append(sql, number);
        }
    }
    return ok && duckvep_sql_append(sql, ", USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)");
}

int duckvep_core_haplotype_script(const char *query, const char *model, const char *job,
    const duckvep_haplotype_script_options_t *o, duckvep_sql_text out[DUCKVEP_HAPLOTYPE_SCRIPT_MAX],
    size_t *count, char *error, size_t error_size) {
    int source_records = o->input_mode && !strcmp(o->input_mode, "source_records");
    duckvep_phase_policy_t policy = DUCKVEP_PHASE_STRICT;
    int policy_valid = !o->phase_policy || duckvep_core_phase_policy(o->phase_policy, strlen(o->phase_policy), &policy);
    int compat = policy == DUCKVEP_PHASE_VEP_COMPAT;
    bool ok = true;
    size_t n = 0u;
    *count = 0u;
    if (!policy_valid ||
        (o->input_mode && strcmp(o->input_mode, "alt_events") && !source_records)) {
        snprintf(error, error_size, "%s", !policy_valid
            ? "duckvep_haplotype_load_sql: phase_policy must be 'strict' or 'vep_compat'"
            : "duckvep_haplotype_load_sql: input_mode must be 'alt_events' or 'source_records'; source_records requires phase_policy='vep_compat'");
        return 2;
    }
    if (source_records && !compat) {
        snprintf(error, error_size,
            "duckvep_haplotype_load_sql: input_mode must be 'alt_events' or 'source_records'; source_records requires phase_policy='vep_compat'");
        return 2;
    }
    for (unsigned i = 0u; i < DUCKVEP_HAP_LIMIT_COUNT; i++) {
        if (o->has_limit[i] && !duckvep_hap_limit_valid(i, o->limit[i])) {
            snprintf(error, error_size, "duckvep_haplotype_load_sql: invalid %s", duckvep_hap_limit_names[i]);
            return 2;
        }
    }
    if (!source_records) {
        duckvep_sql_text *sql = &out[n++];
        /* The newlines keep a trailing line comment in a query from swallowing the parenthesis. */
        ok = duckvep_sql_append(sql, "COPY (\n") && duckvep_sql_append(sql, alt_prefix) && duckvep_sql_append(sql, "\n") &&
            duckvep_sql_append(sql, query) && duckvep_sql_append(sql, "\n") && duckvep_sql_append(sql, alt_middle) &&
            duckvep_sql_append(sql, compat ? alt_domain_compat : alt_domain_strict) && duckvep_sql_append(sql, alt_suffix) &&
            append_copy_options(sql, model, job, NULL, o);
    } else {
        duckvep_sql_text raw = {0}, plain = {0};
        /* The raw relation is a TEMP table: the caller's own transaction sees its uncommitted rows, and the later
         * statements read it twice (the record plan and the call expansion). */
        ok = duckvep_sql_append(&plain, "__duckvep_haplotype_raw_") && duckvep_sql_append(&plain, job) &&
            duckvep_sql_identifier(&raw, plain.data);
        duckvep_sql_free(&plain);
        if (!ok) { duckvep_sql_free(&raw); return 1; }
        duckvep_sql_text *sql = &out[n++];
        ok = duckvep_sql_append(sql, "CREATE TEMP TABLE ") && duckvep_sql_append(sql, raw.data) &&
            duckvep_sql_append(sql, " AS\n") && duckvep_sql_append(sql, raw_columns) && duckvep_sql_append(sql, "\n") &&
            duckvep_sql_append(sql, query) && duckvep_sql_append(sql, "\n) source");
        if (ok) {
            sql = &out[n++];
            ok = duckvep_sql_append(sql, "COPY (\n") && duckvep_sql_append(sql, plan_query_head) &&
                duckvep_sql_append(sql, raw.data) && duckvep_sql_append(sql, plan_query_tail) &&
                append_copy_options(sql, model, job, "plan_input", o);
        }
        if (ok) {
            sql = &out[n++];
            ok = duckvep_sql_append(sql, "COPY (\n") && duckvep_sql_append(sql, final_head) && duckvep_sql_literal(sql, job) &&
                duckvep_sql_append(sql, final_mid1) && duckvep_sql_append(sql, raw.data) &&
                duckvep_sql_append(sql, final_mid2) && duckvep_sql_append(sql, raw.data) &&
                duckvep_sql_append(sql, final_tail) && append_copy_options(sql, model, job, NULL, o);
        }
        if (ok) {
            sql = &out[n++];
            ok = duckvep_sql_append(sql, "DROP TABLE IF EXISTS ") && duckvep_sql_append(sql, raw.data);
        }
        duckvep_sql_free(&raw);
    }
    *count = n;
    return ok ? 0 : 1;
}

static const char *const alt_types[15] = {"UBIGINT", "UINTEGER", "UBIGINT", "VARCHAR", "VARCHAR", "UINTEGER",
    "UINTEGER", "UINTEGER", "INTEGER[]", "BOOLEAN[]", "BIGINT", "BIGINT[]", "BIGINT", "BIGINT", "BIGINT"};
/* raw_gt (index 8) is a STRUCT checked by shape, not by text. */
static const char *const source_types[17] = {"UBIGINT", "UINTEGER", "UBIGINT", "VARCHAR", "VARCHAR", "UINTEGER",
    "UINTEGER", "UINTEGER", NULL, "BOOLEAN[]", "BIGINT", "BIGINT", "BIGINT", "BIGINT", "BIGINT", "UBIGINT", "BOOLEAN"};
static const char *const plan_types[5] = {"UBIGINT", "UINTEGER", "UBIGINT", "VARCHAR", "BIGINT"};

unsigned duckvep_hap_input_columns(int source_records) { return source_records ? 17u : 15u; }

const char *duckvep_hap_input_type(int source_records, unsigned index) {
    if (index >= duckvep_hap_input_columns(source_records)) return NULL;
    return source_records ? source_types[index] : alt_types[index];
}

const char *duckvep_hap_plan_input_type(unsigned index) { return index < 5u ? plan_types[index] : NULL; }
