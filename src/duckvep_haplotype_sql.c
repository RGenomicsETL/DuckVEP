/* Query-scoped DuckDB (v1 C API) adapter for the native phased replay stream. DuckDB owns input materialization,
 * phase-domain aggregation, sorting and output vectors; the replay itself, its workspace and the result writers are the
 * host-neutral src/core/duckvep_core_haplotypes.c. */
#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
extern duckdb_ext_api_v1 duckdb_ext_api;
#include "duckvep_list.h"

#include "duckvep_model.h"
#include "core/duckvep_core_haplotypes.h"
#include "core/duckvep_core_arrangements.h"
#include "core/duckvep_core_haplotype_script.h"
#include "core/duckvep_core_phase.h"
#include "core/duckvep_core_discovery.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIMIT_COUNT DUCKVEP_HAP_LIMIT_COUNT
#define LIMIT_WORKSPACE DUCKVEP_HAP_LIMIT_WORKSPACE
#define LIMIT_PROJECTIONS DUCKVEP_HAP_LIMIT_PROJECTIONS
#define LIMIT_PLOIDY DUCKVEP_HAP_LIMIT_PLOIDY
#define limit_names duckvep_hap_limit_names
#define limit_defaults duckvep_hap_limit_defaults
enum { HAPLOTYPE_LIST_COLUMN = 9, HAPLOTYPE_STOP_COLUMN = 14,
    HAPLOTYPE_PROVENANCE_COLUMN = 22, HAPLOTYPE_CARRIER_PREDICTION_COLUMN = 24,
    HAPLOTYPE_CONSEQUENCES_COLUMN = 25, HAPLOTYPE_NMD_CONTRIBUTORS_COLUMN = 31,
    HAPLOTYPE_OUTPUT_COLUMNS = DUCKVEP_HAP_OUTPUT_COLUMNS };
enum { HAPLOTYPE_BLOCK_EVENT_FIELD = 9, HAPLOTYPE_BLOCK_FIELDS = 10, HAPLOTYPE_PROVENANCE_FIELDS = 10,
    HAPLOTYPE_EDIT_FIELDS = 7, HAPLOTYPE_CARRIER_PREDICTION_FIELDS = 10 };

typedef struct {
    duckvep_registry_t *registry;
    duckvep_model_entry_t *entry;
    char *query;
    duckvep_hap_config_t cfg;
} haplotype_bind_t;

typedef struct {
    duckdb_result input;
    duckdb_data_chunk chunk;
    idx_t row;
    int have_result, eof, source_records;
    duckvep_hap_state_t *core;
} haplotype_state_t;


static void haplotype_bind_destroy(void *pointer) {
    haplotype_bind_t *b = pointer;
    if (!b) return;
    duckvep_registry_unpin(b->registry, b->entry);
    duckvep_registry_release(b->registry);
    duckdb_free(b->query);
    duckvep_budget_free(b);
}

static duckdb_logical_type record_type(const char *const *names, const duckdb_type *ids,
    idx_t count, int list_field) {
    duckdb_logical_type types[HAPLOTYPE_BLOCK_FIELDS];
    const char *field_names[HAPLOTYPE_BLOCK_FIELDS];
    for (idx_t i = 0u; i < count; i++) {
        types[i] = duckdb_create_logical_type(ids[i]); field_names[i] = names[i];
        if (list_field >= 0 && i == (idx_t)list_field) {
            duckdb_logical_type element = types[i];
            types[i] = duckdb_create_list_type(element);
            duckdb_destroy_logical_type(&element);
        }
    }
    duckdb_logical_type type = duckdb_create_struct_type(types, field_names, count);
    for (idx_t i = 0u; i < count; i++) duckdb_destroy_logical_type(&types[i]);
    return type;
}

static void bind_record_list(duckdb_bind_info info, const char *name,
    const char *const *fields, const duckdb_type *ids, idx_t count, int list_field) {
    duckdb_logical_type record = record_type(fields, ids, count, list_field);
    duckdb_logical_type list = duckdb_create_list_type(record);
    duckdb_bind_add_result_column(info, name, list);
    duckdb_destroy_logical_type(&list);
    duckdb_destroy_logical_type(&record);
}

static void haplotype_bind(duckdb_bind_info info) {
    haplotype_bind_t *b = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1u, sizeof(*b));
    if (!b) { duckdb_bind_set_error(info, "duckvep_haplotypes: bind allocation failed"); return; }
    b->registry = duckdb_bind_get_extra_info(info);
    duckvep_registry_retain(b->registry);
    duckdb_value value = duckdb_bind_get_parameter(info, 0u);
    if (value && !duckdb_is_null_value(value)) b->query = duckdb_get_varchar(value);
    duckdb_destroy_value(&value);
    value = duckdb_bind_get_parameter(info, 1u);
    char *name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    duckdb_destroy_value(&value);
    if (name) b->entry = duckvep_registry_pin(b->registry, name);
    duckdb_free(name);
    if (b->entry && b->entry->model.lifted) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: phased edit sets are not supported for models with wrapped circular objects");
        haplotype_bind_destroy(b); return;
    }
    if (!b->entry || !b->query || !b->query[0]) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: require a nonempty calls query and loaded model name");
        haplotype_bind_destroy(b); return;
    }
    value = duckdb_bind_get_named_parameter(info, "phase_policy");
    name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    int valid;
    b->cfg.policy = DUCKVEP_PHASE_STRICT;
    valid = !value || (name && duckvep_core_phase_policy(name, strlen(name), &b->cfg.policy));
    duckdb_free(name); duckdb_destroy_value(&value);
    if (!valid) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: phase_policy must be 'strict' or 'vep_compat'");
        haplotype_bind_destroy(b); return;
    }
    value = duckdb_bind_get_named_parameter(info, "input_mode");
    name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    valid = !value || (name && (!strcmp(name, "alt_events") || !strcmp(name, "source_records")));
    b->cfg.source_records = name && !strcmp(name, "source_records");
    duckdb_free(name); duckdb_destroy_value(&value);
    if (!valid || (b->cfg.source_records && b->cfg.policy != DUCKVEP_PHASE_VEP_COMPAT)) {
        duckdb_bind_set_error(info, "duckvep_haplotypes: input_mode must be 'alt_events' or 'source_records'; source_records requires phase_policy='vep_compat'");
        haplotype_bind_destroy(b); return;
    }
    value = duckdb_bind_get_named_parameter(info, "hgvs");
    b->cfg.hgvs = value && !duckdb_is_null_value(value) && duckdb_get_bool(value);
    duckdb_destroy_value(&value);
    for (unsigned i = 0u; i < LIMIT_COUNT; i++) {
        value = duckdb_bind_get_named_parameter(info, limit_names[i]);
        uint64_t n = value && !duckdb_is_null_value(value) ? duckdb_get_uint64(value) : limit_defaults[i];
        valid = (!value || !duckdb_is_null_value(value)) && duckvep_hap_limit_valid(i, n);
        duckdb_destroy_value(&value);
        if (!valid) {
            char error[128]; snprintf(error, sizeof(error), "duckvep_haplotypes: invalid %s", limit_names[i]);
            duckdb_bind_set_error(info, error); haplotype_bind_destroy(b); return;
        }
        b->cfg.limits[i] = (size_t)n;
    }
    const char *const names[] = {"transcript_index", "cds", "protein", "sequence_flags",
        "evidence_flags", "projection_status", "sequence_status", "edit_count", "carrier_count"};
    const duckdb_type ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UINTEGER};
    for (unsigned i = 0u; i < 9u; i++) {
        duckdb_logical_type type = duckdb_create_logical_type(ids[i]);
        duckdb_bind_add_result_column(info, names[i], type); duckdb_destroy_logical_type(&type);
    }
    const char *const carrier_names[] = {"sample_index", "phase_set", "haplotype_lane", "ploidy"};
    const duckdb_type carrier_ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_BIGINT,
        DUCKDB_TYPE_USMALLINT, DUCKDB_TYPE_USMALLINT};
    const char *const event_names[] = {"event_index", "seq_region", "position", "reference", "alternate",
        "evidence_flags", "projection_status", "alt_index"};
    const duckdb_type event_ids[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT,
        DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UINTEGER};
    const char *const block_names[] = {"cds_start", "reference", "alternate", "alt_start0",
        "length_change", "sequence_flags", "coding_status", "local_consequence_mask",
        "after_first_stop", "event_indices"};
    const duckdb_type block_ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_BIGINT, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_UBIGINT};
    bind_record_list(info, "carriers", carrier_names, carrier_ids, 4u, -1);
    bind_record_list(info, "contributors", event_names, event_ids, b->cfg.source_records ? 8u : 7u, -1);
    bind_record_list(info, "coding_blocks", block_names, block_ids, HAPLOTYPE_BLOCK_FIELDS, HAPLOTYPE_BLOCK_EVENT_FIELD);
    const char *const difference_names[] = {"ref_start0", "alt_start0", "reference", "alternate",
        "alignment_start0"};
    const duckdb_type difference_ids[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
        DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UBIGINT};
    bind_record_list(info, "cds_differences", difference_names, difference_ids, 5u, -1);
    bind_record_list(info, "protein_differences", difference_names, difference_ids, 5u, -1);
    duckdb_logical_type stop_in_frame_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_bind_add_result_column(info, "stop_in_displaced_frame", stop_in_frame_type);
    duckdb_destroy_logical_type(&stop_in_frame_type);
    duckdb_logical_type string_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "hgvsc", string_type);
    duckdb_bind_add_result_column(info, "hgvsc_status", string_type);
    duckdb_bind_add_result_column(info, "hgvsp", string_type);
    duckdb_bind_add_result_column(info, "hgvsp_status", string_type);
    duckdb_destroy_logical_type(&string_type);
    /* eligibility and provenance only. */
    string_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "prediction_policy", string_type);
    duckdb_bind_add_result_column(info, "prediction_status", string_type);
    duckdb_bind_add_result_column(info, "prediction_reason", string_type);
    duckdb_destroy_logical_type(&string_type);
    const char *const provenance_names[] = {"event_index", "alt_index", "seq_region", "position",
        "reference", "alternate", "evidence_flags", "projection_status", "role", "edit_count"};
    const duckdb_type provenance_ids[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UINTEGER,
        DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UTINYINT, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UINTEGER};
    const char *const edit_names[] = {"edit_index", "event_index", "block_index", "cds_start",
        "reference", "alternate", "variant_strand"};
    const duckdb_type edit_ids[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT,
        DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_TINYINT};
    bind_record_list(info, "contributor_provenance", provenance_names, provenance_ids,
        HAPLOTYPE_PROVENANCE_FIELDS, -1);
    bind_record_list(info, "normalized_edits", edit_names, edit_ids, HAPLOTYPE_EDIT_FIELDS, -1);
    const char *const carrier_prediction_names[] = {"sample_index", "phase_set", "haplotype_lane",
        "prediction_status", "prediction_reason", "haplotype_impact", "haplotype_consequences",
        "nmd_prediction", "nmd_stop_position", "nmd_junction_position"};
    const duckdb_type carrier_prediction_ids[] = {DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_BIGINT,
        DUCKDB_TYPE_USMALLINT, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UBIGINT};
    /* Slice 4: the consequence set and IMPACT of the shared edited sequence, decided per carrier so an
     * ineligible carrier of the row (for example a triploid call) never hides an eligible one. */
    bind_record_list(info, "carrier_predictions", carrier_prediction_names, carrier_prediction_ids,
        HAPLOTYPE_CARRIER_PREDICTION_FIELDS, 6);
    /* Slice 3: whole-haplotype reduced SO set and IMPACT of the same-codon classifier. */
    duckdb_logical_type so_element = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type so_list = duckdb_create_list_type(so_element);
    duckdb_bind_add_result_column(info, "haplotype_consequences", so_list);
    duckdb_destroy_logical_type(&so_list);
    duckdb_destroy_logical_type(&so_element);
    string_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "haplotype_impact", string_type);
    duckdb_destroy_logical_type(&string_type);
    /* Slice 6: ejc50 NMD of the shared edited transcript (row summary; per carrier in carrier_predictions). */
    string_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "nmd_rule", string_type);
    duckdb_bind_add_result_column(info, "nmd_prediction", string_type);
    duckdb_destroy_logical_type(&string_type);
    duckdb_logical_type length_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_bind_add_result_column(info, "nmd_stop_position", length_type);
    duckdb_bind_add_result_column(info, "nmd_junction_position", length_type);
    duckdb_logical_type index_list = duckdb_create_list_type(length_type);
    duckdb_bind_add_result_column(info, "nmd_contributors", index_list);
    duckdb_destroy_logical_type(&index_list);
    duckdb_destroy_logical_type(&length_type);
    string_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "nmd_exceptions", string_type);
    duckdb_bind_add_result_column(info, "prediction_reference_protein", string_type);
    duckdb_bind_add_result_column(info, "prediction_protein", string_type);
    duckdb_destroy_logical_type(&string_type);
    length_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "nominal_length_diff", length_type);
    duckdb_destroy_logical_type(&length_type);
    duckdb_bind_set_bind_data(info, b, haplotype_bind_destroy);
}

static void haplotype_state_destroy(void *pointer) {
    haplotype_state_t *s = pointer;
    if (!s) return;
    if (s->chunk) duckdb_destroy_data_chunk(&s->chunk);
    if (s->have_result) duckdb_destroy_result(&s->input);
    duckvep_hap_close(s->core);
    duckvep_budget_free(s);
}


static int raw_prepare_query(duckdb_connection connection, const char *sql,
    duckdb_result *result, char *error, size_t error_size) {
    if (duckdb_query(connection, sql, result) == DuckDBSuccess) return 1;
    duckvep_sql_set_error(error, error_size, duckdb_result_error(result));
    return 0;
}

/* Called under the registry preparation lock. The SELECT is evaluated once;
 * DuckDB owns both temporary relations and their spill. No native record array
 * grows with the source stream. Drop only tables created by this invocation. */
static int raw_prepare(haplotype_state_t *s, const haplotype_bind_t *b,
    char *error, size_t error_size) {
    duckdb_connection connection = b->registry->query_connection;
    duckdb_result result = {0};
    duckdb_appender appender = NULL;
    duckdb_data_chunk chunk = NULL;
    int have_raw = 0, have_order = 0, ok = 0;
    if (!raw_prepare_query(connection,
        "CREATE TEMP TABLE __duckvep_haplotype_raw(event_index UBIGINT, seq_region UINTEGER, "
        "position UBIGINT, reference VARCHAR, alternates VARCHAR[], transcript_index UINTEGER, "
        "sample_index UINTEGER, gt VARCHAR)", &result, error, error_size)) goto cleanup;
    have_raw = 1; duckdb_destroy_result(&result);
    if (duckdb_appender_create(connection, NULL, "__duckvep_haplotype_raw", &appender) != DuckDBSuccess)
        goto append_failed;
    while ((chunk = duckdb_fetch_chunk(s->input))) {
        if (duckdb_append_data_chunk(appender, chunk) != DuckDBSuccess) goto append_failed;
        duckdb_destroy_data_chunk(&chunk);
    }
    if (duckdb_result_error(&s->input)) {
        duckvep_sql_set_error(error, error_size, duckdb_result_error(&s->input)); goto cleanup;
    }
    if (duckdb_appender_close(appender) != DuckDBSuccess) goto append_failed;
    duckdb_appender_destroy(&appender);
    duckdb_destroy_result(&s->input); s->have_result = 0;
    if (!raw_prepare_query(connection,
        "CREATE TEMP TABLE __duckvep_haplotype_order(event_index UBIGINT, buffer_id UBIGINT, "
        "ordinal UBIGINT, source_ordinal UBIGINT)", &result, error, error_size)) goto cleanup;
    have_order = 1; duckdb_destroy_result(&result);
    if (!raw_prepare_query(connection,
        "SELECT event_index, first(seq_region), first(position), first(reference), "
        "count(DISTINCT (seq_region,position,reference,alternates)) versions "
        "FROM __duckvep_haplotype_raw GROUP BY event_index ORDER BY 2,3,1",
        &result, error, error_size)) goto cleanup;
    if (duckdb_appender_create(connection, NULL, "__duckvep_haplotype_order", &appender) != DuckDBSuccess)
        goto append_failed;
    duckvep_haplotype_record_plan_t plan;
    if (!duckvep_haplotype_record_plan_init(&plan, &b->entry->model.transcripts)) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid source-planning model");
        goto cleanup;
    }
    uint64_t source_ordinal = 0u;
    while ((chunk = duckdb_fetch_chunk(result))) {
        duckdb_vector v[5];
        for (unsigned i = 0u; i < 5u; i++) v[i] = duckdb_data_chunk_get_vector(chunk, i);
        for (idx_t row = 0u; row < duckdb_data_chunk_get_size(chunk); row++) {
            for (unsigned i = 0u; i < 4u; i++) if (duckvep_row_is_null(v[i], row)) {
                duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: required input source column is NULL");
                goto cleanup;
            }
            if (((int64_t *)duckdb_vector_get_data(v[4]))[row] != 1) {
                duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: inconsistent source record identity");
                goto cleanup;
            }
            uint32_t chrom = ((uint32_t *)duckdb_vector_get_data(v[1]))[row];
            uint64_t pos = ((uint64_t *)duckdb_vector_get_data(v[2]))[row];
            duckdb_string_t ref = ((duckdb_string_t *)duckdb_vector_get_data(v[3]))[row];
            uint32_t length = duckdb_string_t_length(ref);
            uint64_t buffer, ordinal;
            if (chrom > UINT16_MAX || !pos || pos > UINT32_MAX || !length || length > UINT16_MAX ||
                length - 1u > UINT32_MAX - pos || source_ordinal == UINT64_MAX ||
                !duckvep_haplotype_record_plan_next(&plan, (uint16_t)chrom, (uint32_t)pos,
                    (uint32_t)(pos + length - 1u), &buffer, &ordinal)) {
                duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid source span or record count");
                goto cleanup;
            }
            source_ordinal++;
            if (duckdb_append_uint64(appender, ((uint64_t *)duckdb_vector_get_data(v[0]))[row]) != DuckDBSuccess ||
                duckdb_append_uint64(appender, buffer) != DuckDBSuccess ||
                duckdb_append_uint64(appender, ordinal) != DuckDBSuccess ||
                duckdb_append_uint64(appender, source_ordinal) != DuckDBSuccess ||
                duckdb_appender_end_row(appender) != DuckDBSuccess) goto append_failed;
        }
        duckdb_destroy_data_chunk(&chunk);
    }
    if (duckdb_result_error(&result)) {
        duckvep_sql_set_error(error, error_size, duckdb_result_error(&result)); goto cleanup;
    }
    if (duckdb_appender_close(appender) != DuckDBSuccess) goto append_failed;
    duckdb_appender_destroy(&appender); duckdb_destroy_result(&result);
    s->have_result = 1;
    ok = raw_prepare_query(connection,
        "WITH ordered AS (SELECT event_index, CASE WHEN buffer_id=0 THEN source_ordinal ELSE "
        "source_ordinal-ordinal+_duckvep_record_order(count(*) OVER(PARTITION BY buffer_id)::UBIGINT,ordinal) "
        "END replay_order FROM __duckvep_haplotype_order), "
        "genotypes AS (SELECT event_index,sample_index,count(DISTINCT gt) gt_versions "
        "FROM __duckvep_haplotype_raw GROUP BY event_index,sample_index), "
        "calls AS MATERIALIZED (SELECT *, count(*) OVER(PARTITION BY event_index,transcript_index,sample_index) copies, "
        "CASE WHEN alternates IS NULL OR len(alternates)>2147483647 OR "
        "len(list_filter(alternates,lambda a: a IS NULL OR len(a)=0 OR len(a)>65535))>0 "
        "THEN error('duckvep_haplotypes: invalid source ALT list') ELSE len(alternates) END alt_count "
        "FROM __duckvep_haplotype_raw JOIN ordered USING(event_index)), "
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
        "ORDER BY seq_region,position,event_index,alt_index,transcript_index,sample_index",
        &s->input, error, error_size);
    goto cleanup;
append_failed:
    duckvep_sql_set_error(error, error_size, appender ? duckdb_appender_error(appender)
        : "duckvep_haplotypes: could not stage source records");
cleanup:
    if (chunk) duckdb_destroy_data_chunk(&chunk);
    if (appender) duckdb_appender_destroy(&appender);
    duckdb_destroy_result(&result);
    const char *drops[] = {"DROP TABLE temp.main.__duckvep_haplotype_order", "DROP TABLE temp.main.__duckvep_haplotype_raw"};
    const int created[] = {have_order, have_raw};
    for (unsigned i = 0u; i < 2u; i++) if (created[i]) {
        if (duckdb_query(connection, drops[i], &result) != DuckDBSuccess && ok) {
            duckvep_sql_set_error(error, error_size, duckdb_result_error(&result)); ok = 0;
        }
        duckdb_destroy_result(&result);
    }
    return ok;
}

static int input_open(haplotype_state_t *s, const haplotype_bind_t *b, char *error, size_t error_size) {
    /* One SELECT statement and snapshot on the registry's retained connection.
     * Only query preparation/materialization is serialized. The returned result
     * is owned by this scan and fetched without the registry mutex. TEMP objects
     * and uncommitted writes on the caller's connection are not visible here. */
    const char *prefix =
        "WITH raw AS MATERIALIZED (SELECT event_index::UBIGINT event_index, seq_region::UINTEGER seq_region, "
        "position::UBIGINT AS position, reference::VARCHAR AS reference, alternate::VARCHAR AS alternate, "
        "alt_index::UINTEGER alt_index, transcript_index::UINTEGER transcript_index, sample_index::UINTEGER sample_index, "
        "alleles::INTEGER[] alleles, phase_before::BOOLEAN[] phase_before, phase_set::BIGINT phase_set FROM (";
    const char *middle = b->cfg.policy == DUCKVEP_PHASE_STRICT ? DUCKVEP_HAP_ALT_MIDDLE_STRICT
                                                                 : DUCKVEP_HAP_ALT_MIDDLE_COMPAT;
    const char *domain = b->cfg.policy == DUCKVEP_PHASE_STRICT ?
        "coalesce(list(DISTINCT phase_set ORDER BY phase_set NULLS FIRST) FILTER(WHERE scoped), [NULL]::BIGINT[]) AS domain_sets " :
        "[NULL]::BIGINT[] AS domain_sets ";
    /* Validate each identity/domain once, then attach its facts to every call.
     * Duplicate calls remain visible, including calls carrying only REF. */
    const char *suffix =
        ",count(DISTINCT len(alleles)) ploidies FROM calls GROUP BY transcript_index,sample_index), "
        "event_versions AS (SELECT event_index, "
        "count(DISTINCT (seq_region,position,reference,alternate,alt_index)) versions FROM raw GROUP BY event_index) "
        "SELECT c.* EXCLUDE(scoped), d.domain_sets, "
        "count(*) OVER(PARTITION BY c.event_index,transcript_index,sample_index) copies, v.versions, d.ploidies "
        "FROM calls c LEFT JOIN domains d USING(transcript_index,sample_index) "
        "LEFT JOIN event_versions v USING(event_index) "
        "ORDER BY seq_region,position,event_index,transcript_index,sample_index";
    if (b->cfg.source_records) {
        prefix = "SELECT event_index::UBIGINT event_index, "
            "seq_region::UINTEGER seq_region, position::UBIGINT AS position, reference::VARCHAR AS reference, "
            "alternates::VARCHAR[] alternates, transcript_index::UINTEGER transcript_index, "
            "sample_index::UINTEGER sample_index, gt::VARCHAR gt FROM (";
        middle = ") source";
        domain = ""; suffix = "";
    }
    size_t qlen = strlen(b->query), overhead = strlen(prefix) + strlen(middle) + strlen(domain) + strlen(suffix) + 1u;
    if (qlen > SIZE_MAX - overhead || qlen + overhead > b->cfg.limits[LIMIT_WORKSPACE] - duckvep_hap_workspace_bytes(s->core)) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: workspace_limit exceeded by calls query text");
        return 0;
    }
    char *sql = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, qlen + overhead);
    if (!sql) return 0;
    snprintf(sql, qlen + overhead, "%s%s%s%s%s", prefix, b->query, middle, domain, suffix);
    duckdb_prepared_statement statement = NULL;
    if (!duckvep_registry_query_acquire(b->registry, error, error_size)) { duckvep_budget_free(sql); return 0; }
    duckdb_extracted_statements extracted = NULL;
    idx_t statements = duckdb_extract_statements(b->registry->query_connection, sql, &extracted);
    int ok = statements == 1u && duckdb_prepare_extracted_statement(
        b->registry->query_connection, extracted, 0u, &statement) == DuckDBSuccess;
    duckvep_budget_free(sql);
    if (!ok) {
        const char *message = statement ? duckdb_prepare_error(statement) :
            extracted ? duckdb_extract_statements_error(extracted) : NULL;
        duckvep_sql_set_error(error, error_size, message ? message : "calls query must be one SELECT statement");
    } else if (duckdb_prepared_statement_type(statement) != DUCKDB_STATEMENT_TYPE_SELECT) {
        duckvep_sql_set_error(error, error_size, "calls query must be one SELECT statement"); ok = 0;
    } else {
        s->have_result = 1;
        ok = duckdb_execute_prepared(statement, &s->input) == DuckDBSuccess;
        if (!ok) duckvep_sql_set_error(error, error_size, duckdb_result_error(&s->input));
    }
    if (statement) duckdb_destroy_prepare(&statement);
    if (extracted) duckdb_destroy_extracted(&extracted);
    if (ok && b->cfg.source_records) ok = raw_prepare(s, b, error, error_size);
    pthread_mutex_unlock(&b->registry->query_mutex);
    return ok;
}

static void haplotype_init(duckdb_init_info info) {
    const haplotype_bind_t *bind = duckdb_init_get_bind_data(info);
    haplotype_state_t *s = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1u, sizeof(*s));
    char error[DUCKVEP_SQL_ERROR_SIZE] = "duckvep_haplotypes: workspace allocation or configured limit exceeded";
    duckdb_init_set_max_threads(info, 1u);
    duckvep_budget_clear_failure();
    if (!s) goto failed;
    s->source_records = bind->cfg.source_records;
    duckvep_hap_config_t cfg = bind->cfg;
    cfg.model = &bind->entry->model;
    s->core = duckvep_hap_open(&cfg, error, sizeof(error));
    if (!s->core) goto failed;
    if (!input_open(s, bind, error, sizeof(error))) goto failed;
    duckdb_init_set_init_data(info, s, haplotype_state_destroy);
    return;
failed:
    { char final_error[DUCKVEP_SQL_ERROR_SIZE + 256];
      duckdb_init_set_error(info, duckvep_sql_final_error(final_error, sizeof final_error, error, error)); }
    haplotype_state_destroy(s);
}


/* Fills the neutral input row from the current row of the fetched chunk. */
static int input_next(void *context, duckvep_hap_row_t *out, char *error, size_t error_size) {
    (void)error; (void)error_size;
    haplotype_state_t *s = context;
    while (!s->chunk || s->row == duckdb_data_chunk_get_size(s->chunk)) {
        if (s->eof) return 0;
        if (s->chunk) duckdb_destroy_data_chunk(&s->chunk);
        s->chunk = duckdb_fetch_chunk(s->input); s->row = 0u;
        if (!s->chunk) { s->eof = 1; return 0; }
    }
    const idx_t row = s->row++;
    const int source_records = s->source_records;
    const unsigned columns = source_records ? 17u : 15u;
    duckdb_vector v[17];
    memset(out, 0, sizeof *out);
    for (unsigned i = 0u; i < columns; i++) {
        v[i] = duckdb_data_chunk_get_vector(s->chunk, i);
        if (duckvep_row_is_null(v[i], row)) out->null_mask |= UINT32_C(1) << i;
    }
    const uint32_t required = ((UINT32_C(1) << columns) - 1u) & ~((UINT32_C(1) << 9) | (UINT32_C(1) << 10));
    if (out->null_mask & required) return 1;   /* the core reports the first NULL column */
    out->event_id = ((uint64_t *)duckdb_vector_get_data(v[0]))[row];
    out->chrom = ((uint32_t *)duckdb_vector_get_data(v[1]))[row];
    out->pos = ((uint64_t *)duckdb_vector_get_data(v[2]))[row];
    /* Point into the vector: an inlined string's bytes live in its cell, so a local copy would dangle. */
    duckdb_string_t *ref = &((duckdb_string_t *)duckdb_vector_get_data(v[3]))[row];
    duckdb_string_t *alt = &((duckdb_string_t *)duckdb_vector_get_data(v[4]))[row];
    out->ref_len = duckdb_string_t_length(*ref); out->alt_len = duckdb_string_t_length(*alt);
    out->ref = (const uint8_t *)duckdb_string_t_data(ref); out->alt = (const uint8_t *)duckdb_string_t_data(alt);
    out->allele_index = ((uint32_t *)duckdb_vector_get_data(v[5]))[row];
    out->transcript = ((uint32_t *)duckdb_vector_get_data(v[6]))[row];
    out->sample = ((uint32_t *)duckdb_vector_get_data(v[7]))[row];
    out->copies = ((int64_t *)duckdb_vector_get_data(v[12]))[row];
    out->versions = ((int64_t *)duckdb_vector_get_data(v[13]))[row];
    out->ploidies = ((int64_t *)duckdb_vector_get_data(v[14]))[row];
    if (source_records) {
        for (unsigned i = 0u; i < 7u; i++)
            out->raw[i] = ((uint32_t *)duckdb_vector_get_data(duckdb_struct_vector_get_child(v[8], i)))[row];
        out->replay_order = ((uint64_t *)duckdb_vector_get_data(v[15]))[row];
        out->source_selected = ((bool *)duckdb_vector_get_data(v[16]))[row];
        return 1;
    }
    duckdb_list_entry gt = ((duckdb_list_entry *)duckdb_vector_get_data(v[8]))[row];
    duckdb_vector g = duckdb_list_vector_get_child(v[8]);
    out->gt = duckdb_vector_get_data(g); out->gt_validity = duckdb_vector_get_validity(g);
    out->gt_offset = gt.offset; out->gt_length = gt.length;
    out->have_phase = !(out->null_mask >> 9 & 1u);
    if (out->have_phase) {
        duckdb_list_entry phases = ((duckdb_list_entry *)duckdb_vector_get_data(v[9]))[row];
        duckdb_vector p = duckdb_list_vector_get_child(v[9]);
        out->phase = duckdb_vector_get_data(p); out->phase_validity = duckdb_vector_get_validity(p);
        out->phase_offset = phases.offset; out->phase_length = phases.length;
    }
    out->phase_set_present = !(out->null_mask >> 10 & 1u);
    if (out->phase_set_present) out->phase_set = ((int64_t *)duckdb_vector_get_data(v[10]))[row];
    duckdb_list_entry sets = ((duckdb_list_entry *)duckdb_vector_get_data(v[11]))[row];
    duckdb_vector domains = duckdb_list_vector_get_child(v[11]);
    out->sets = duckdb_vector_get_data(domains); out->sets_validity = duckdb_vector_get_validity(domains);
    out->sets_offset = sets.offset; out->sets_length = sets.length;
    return 1;
}

static void haplotype_scan(duckdb_function_info info, duckdb_data_chunk output) {
    haplotype_state_t *s = duckdb_function_get_init_data(info);
    idx_t capacity = duckdb_vector_size();
    for (unsigned i = 0u; i < HAPLOTYPE_OUTPUT_COLUMNS; i++)
        duckdb_vector_ensure_validity_writable(duckdb_data_chunk_get_vector(output, i));
    for (unsigned i = HAPLOTYPE_LIST_COLUMN; i < HAPLOTYPE_STOP_COLUMN; i++)
        if (duckdb_list_vector_set_size(duckdb_data_chunk_get_vector(output, i), 0u) != DuckDBSuccess) {
            duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset output list"); return;
        }
    for (unsigned i = HAPLOTYPE_PROVENANCE_COLUMN; i <= HAPLOTYPE_CONSEQUENCES_COLUMN; i++)
        if (duckdb_list_vector_set_size(duckdb_data_chunk_get_vector(output, i), 0u) != DuckDBSuccess) {
            duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset output list"); return;
        }
    duckdb_vector carrier_records = duckdb_list_vector_get_child(
        duckdb_data_chunk_get_vector(output, HAPLOTYPE_CARRIER_PREDICTION_COLUMN));
    if (duckdb_list_vector_set_size(duckdb_struct_vector_get_child(carrier_records, 6u), 0u) != DuckDBSuccess) {
        duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset carrier consequence list"); return;
    }
    if (duckdb_list_vector_set_size(duckdb_data_chunk_get_vector(output, HAPLOTYPE_NMD_CONTRIBUTORS_COLUMN), 0u) != DuckDBSuccess) {
        duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset NMD contributor list"); return;
    }
    duckdb_vector blocks = duckdb_list_vector_get_child(duckdb_data_chunk_get_vector(output, 11u));
    if (duckdb_list_vector_set_size(duckdb_struct_vector_get_child(blocks, HAPLOTYPE_BLOCK_EVENT_FIELD), 0u) != DuckDBSuccess) {
        duckdb_function_set_error(info, "duckvep_haplotypes: cannot reset block event list"); return;
    }
    char error[DUCKVEP_SQL_ERROR_SIZE] = {0};
    duckvep_hap_input_t input = {s, input_next};
    size_t rows = 0u;
    if (!duckvep_hap_scan(s->core, &input, output, capacity, &rows, error, sizeof(error))) {
        duckdb_function_set_error(info, error); return;
    }
    duckdb_data_chunk_set_size(output, rows);
}


/* Transcript discovery for phased calls: the transcripts of a loaded model whose coding sequence the event (a VCF record's
 * REF and one ALT) overlaps, in ascending model ordinal. The rules live in duckvep_discovery.c, shared with the fused
 * reader duckvep_coding_calls. Lifted circular models are refused, as by duckvep_haplotypes. */

/* The three numeric parameters are HUGEINT so that every integer type, and sums such as UBIGINT + BIGINT that DuckDB
 * widens, bind without a cast. Values outside int64 are rejected as malformed spans (regions: no such region). */
static int hugeint_to_i64(duckdb_hugeint value, int64_t *out) {
    if (value.upper == 0 && value.lower <= (uint64_t)INT64_MAX) { *out = (int64_t)value.lower; return 1; }
    if (value.upper == -1 && value.lower >= (uint64_t)INT64_MAX + 1u) { *out = (int64_t)value.lower; return 1; }
    return 0;
}

typedef struct {
    duckvep_registry_t *registry;
    duckvep_model_entry_t *entry;
    duckvep_discovery_t scratch;
} exon_discovery_t;

static void exon_discovery_release(exon_discovery_t *d) {
    if (d->entry) duckvep_registry_unpin(d->registry, d->entry);
    d->entry = NULL;
    duckvep_discovery_release(&d->scratch);
}

static void coding_transcripts_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    exon_discovery_t d = {duckdb_scalar_function_get_extra_info(info), NULL, {NULL, 0, UINT32_MAX, -1}};
    idx_t rows = duckdb_data_chunk_get_size(input);
    duckdb_vector name_vector = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector region_vector = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector start_vector = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector ref_vector = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector alt_vector = duckdb_data_chunk_get_vector(input, 4);
    const duckdb_hugeint *regions = duckdb_vector_get_data(region_vector);
    const duckdb_hugeint *starts = duckdb_vector_get_data(start_vector);
    const duckdb_string_t *refs = duckdb_vector_get_data(ref_vector), *alts = duckdb_vector_get_data(alt_vector);
    duckdb_list_entry *entries = duckdb_vector_get_data(output);
    duckvep_u32_list_t out = {NULL, 0u, 0u};
    char *current_name = NULL;
    const char *error = NULL;
    duckdb_vector_ensure_validity_writable(output);
    uint64_t *validity = duckdb_vector_get_validity(output);
    const duckdb_string_t *names = duckdb_vector_get_data(name_vector);
    uint64_t *input_validity[5] = {duckdb_vector_get_validity(name_vector), duckdb_vector_get_validity(region_vector),
        duckdb_vector_get_validity(start_vector), duckdb_vector_get_validity(ref_vector), duckdb_vector_get_validity(alt_vector)};
    size_t current_length = 0u;
    for (idx_t row = 0; row < rows && !error; row++) {
        int missing = 0;
        for (unsigned k = 0u; k < 5u; k++)
            missing |= input_validity[k] && !duckdb_validity_row_is_valid(input_validity[k], row);
        if (missing) {
            duckdb_validity_set_row_invalid(validity, row);
            entries[row].offset = out.count; entries[row].length = 0u;
            continue;
        }
        /* The model name is almost always one constant: compare in place and copy it only when it changes. */
        duckdb_string_t name_cell = names[row];
        size_t name_length = duckdb_string_t_length(name_cell);
        const char *name_data = duckdb_string_t_data(&name_cell);
        if (!current_name || name_length != current_length || memcmp(current_name, name_data, name_length)) {
            char *name = duckvep_vector_string(name_vector, row);
            if (!name) { error = "duckvep_coding_transcripts: out of memory copying the model name"; break; }
            if (d.entry) duckvep_registry_unpin(d.registry, d.entry);
            d.entry = duckvep_registry_pin(d.registry, name);
            duckvep_budget_free(current_name);
            current_name = name; current_length = strlen(name);
            duckvep_discovery_reset_model(&d.scratch);
            if (!d.entry) { error = "duckvep_coding_transcripts: unknown model name"; break; }
            if (d.entry->model.lifted) {
                error = "duckvep_coding_transcripts: transcript discovery is not supported for models with wrapped circular objects";
                break;
            }
        }
        const duckvep_owned_model_t *m = &d.entry->model;
        int64_t start = 0, region = -1;
        duckdb_string_t ref_cell = refs[row], alt_cell = alts[row];
        entries[row].offset = out.count; entries[row].length = 0u;
        if (!hugeint_to_i64(starts[row], &start)) start = 0;
        if (!hugeint_to_i64(regions[row], &region)) region = -1;
        size_t added = 0u;
        duckvep_discovery_status_t status = duckvep_discover_coding(m, &d.scratch, region, start,
            (const uint8_t *)duckdb_string_t_data(&ref_cell), duckdb_string_t_length(ref_cell),
            (const uint8_t *)duckdb_string_t_data(&alt_cell), duckdb_string_t_length(alt_cell), &out, &added);
        if (status == DUCKVEP_DISCOVERY_BAD_SPAN)
            error = "duckvep_coding_transcripts: position must be from 1 through 2147483647 and alleles at most 65535 bases";
        else if (status == DUCKVEP_DISCOVERY_NOMEM)
            error = "duckvep_coding_transcripts: out of memory collecting transcripts";
        else entries[row].length = added;
    }
    if (!error && out.count) {
        if (duckdb_list_vector_reserve(output, out.count) != DuckDBSuccess) error = "duckvep_coding_transcripts: out of memory";
        else {
            memcpy(duckdb_vector_get_data(duckdb_list_vector_get_child(output)), out.items, out.count * sizeof *out.items);
            duckdb_list_vector_set_size(output, out.count);
        }
    } else if (!error) duckdb_list_vector_set_size(output, 0);
    duckvep_u32_list_release(&out);
    duckvep_budget_free(current_name);
    exon_discovery_release(&d);
    if (error) duckdb_scalar_function_set_error(info, error);
}

static void register_coding_transcripts(duckdb_connection connection, duckvep_registry_t *registry) {
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_logical_type string = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type region = duckdb_create_logical_type(DUCKDB_TYPE_HUGEINT);
    duckdb_logical_type position = duckdb_create_logical_type(DUCKDB_TYPE_HUGEINT);
    duckdb_logical_type element = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type result = duckdb_create_list_type(element);
    duckdb_scalar_function_set_name(function, "duckvep_coding_transcripts");
    duckdb_scalar_function_add_parameter(function, string);
    duckdb_scalar_function_add_parameter(function, region);
    duckdb_scalar_function_add_parameter(function, position);
    duckdb_scalar_function_add_parameter(function, string);
    duckdb_scalar_function_add_parameter(function, string);
    duckdb_scalar_function_set_return_type(function, result);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_scalar_function_set_volatile(function);
    duckvep_registry_retain(registry);
    duckdb_scalar_function_set_extra_info(function, registry, duckvep_registry_release);
    duckdb_scalar_function_set_function(function, coding_transcripts_scalar);
    (void)duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&result); duckdb_destroy_logical_type(&element);
    duckdb_destroy_logical_type(&position); duckdb_destroy_logical_type(&region);
    duckdb_destroy_logical_type(&string);
}

typedef struct {
    duckvep_registry_t *registry;
    duckvep_model_entry_t *entry;
    char *query;
    duckvep_arrangement_config_t config;
} arrangement_bind_t;

typedef struct {
    haplotype_state_t input;
    duckvep_arrangement_state_t *core;
} arrangement_state_t;

static void arrangement_bind_destroy(void *pointer) {
    arrangement_bind_t *bind = pointer;
    if (!bind) return;
    duckvep_registry_unpin(bind->registry, bind->entry);
    duckvep_registry_release(bind->registry);
    duckdb_free(bind->query);
    duckvep_budget_free(bind);
}

static void arrangement_state_destroy(void *pointer) {
    arrangement_state_t *state = pointer;
    if (!state) return;
    if (state->input.chunk) duckdb_destroy_data_chunk(&state->input.chunk);
    if (state->input.have_result) duckdb_destroy_result(&state->input.input);
    duckvep_arrangement_close(state->core);
    duckvep_budget_free(state);
}

static int arrangement_limit(duckdb_bind_info info, const char *name, size_t fallback, size_t *value) {
    duckdb_value raw = duckdb_bind_get_named_parameter(info, name);
    uint64_t number = raw && !duckdb_is_null_value(raw) ? duckdb_get_uint64(raw) : fallback;
    int ok = (!raw || !duckdb_is_null_value(raw)) && number && number <= SIZE_MAX;
    duckdb_destroy_value(&raw);
    if (ok) *value = (size_t)number;
    return ok;
}

static void arrangement_bind(duckdb_bind_info info) {
    arrangement_bind_t *bind = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1u, sizeof(*bind));
    duckdb_value value;
    char *name;
    if (!bind) { duckdb_bind_set_error(info, "duckvep_haplotype_arrangements: bind allocation failed"); return; }
    bind->registry = duckdb_bind_get_extra_info(info);
    duckvep_registry_retain(bind->registry);
    value = duckdb_bind_get_parameter(info, 0u);
    if (value && !duckdb_is_null_value(value)) bind->query = duckdb_get_varchar(value);
    duckdb_destroy_value(&value);
    value = duckdb_bind_get_parameter(info, 1u);
    name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    duckdb_destroy_value(&value);
    if (name) bind->entry = duckvep_registry_pin(bind->registry, name);
    duckdb_free(name);
    if (!bind->entry || !bind->query || !bind->query[0]) {
        duckdb_bind_set_error(info, "duckvep_haplotype_arrangements: require a nonempty alt_events query and loaded model name");
        arrangement_bind_destroy(bind); return;
    }
    if (bind->entry->model.lifted) {
        duckdb_bind_set_error(info, "duckvep_haplotype_arrangements: replay is not supported for models with wrapped circular objects");
        arrangement_bind_destroy(bind); return;
    }
    duckvep_arrangement_config_defaults(&bind->config, &bind->entry->model);
    if (!arrangement_limit(info, "max_sites", bind->config.max_sites, &bind->config.max_sites) ||
        !arrangement_limit(info, "max_calls", bind->config.max_calls, &bind->config.max_calls) ||
        !arrangement_limit(info, "max_arrangements", bind->config.max_arrangements, &bind->config.max_arrangements) ||
        !arrangement_limit(info, "max_replays", bind->config.max_replays, &bind->config.max_replays) ||
        bind->config.max_sites > DUCKVEP_PHASE_ARRANGEMENT_MAX_SITES) {
        duckdb_bind_set_error(info, "duckvep_haplotype_arrangements: invalid bounded limit");
        arrangement_bind_destroy(bind); return;
    }
    const char *const names[] = {"hypothesis_id", "hypothesis_lane", "hypothesis_reference_lane",
        "event_index", "seq_region", "position", "reference", "alternate", "alt_index", "transcript_index",
        "sample_index", "original_allele0", "original_allele1", "original_phase_before0",
        "original_phase_before1", "original_phase_before_present", "original_phase_set_present",
        "original_phase_set", "assigned_allele", "contributes", "cds", "protein", "nominal_length_diff",
        "prediction_status", "consequence_mask", "prediction_semantics"};
    const duckdb_type types[] = {DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_USMALLINT, DUCKDB_TYPE_BOOLEAN,
        DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_UINTEGER, DUCKDB_TYPE_INTEGER, DUCKDB_TYPE_INTEGER,
        DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_BIGINT,
        DUCKDB_TYPE_INTEGER, DUCKDB_TYPE_BOOLEAN, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_BIGINT,
        DUCKDB_TYPE_VARCHAR, DUCKDB_TYPE_UBIGINT, DUCKDB_TYPE_VARCHAR};
    if (sizeof(types) / sizeof(types[0]) != DUCKVEP_ARRANGEMENT_COLUMN_COUNT) {
        duckdb_bind_set_error(info, "duckvep_haplotype_arrangements: shared output schema mismatch");
        arrangement_bind_destroy(bind); return;
    }
    for (unsigned i = 0u; i < DUCKVEP_ARRANGEMENT_COLUMN_COUNT; i++) {
        duckdb_logical_type type = duckdb_create_logical_type(types[i]);
        duckdb_bind_add_result_column(info, names[i], type);
        duckdb_destroy_logical_type(&type);
    }
    duckdb_bind_set_bind_data(info, bind, arrangement_bind_destroy);
}

static int arrangement_input_open(arrangement_state_t *state, const arrangement_bind_t *bind,
    char *error, size_t error_size) {
    duckvep_sql_text sql = {0};
    duckdb_prepared_statement statement = NULL;
    duckdb_extracted_statements extracted = NULL;
    char *query_sql;
    int ok;
    if (duckvep_core_arrangement_input_sql(bind->query, &sql, error, error_size) != 0) {
        return 0;
    }
    query_sql = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, sql.length + 1u);
    if (!query_sql) { duckvep_sql_free(&sql); return 0; }
    memcpy(query_sql, sql.data, sql.length + 1u);
    duckvep_sql_free(&sql);
    if (!duckvep_registry_query_acquire(bind->registry, error, error_size)) { duckvep_budget_free(query_sql); return 0; }
    idx_t statements = duckdb_extract_statements(bind->registry->query_connection, query_sql, &extracted);
    duckvep_budget_free(query_sql);
    ok = statements == 1u && duckdb_prepare_extracted_statement(bind->registry->query_connection, extracted, 0u, &statement) == DuckDBSuccess;
    if (!ok) {
        const char *message = statement ? duckdb_prepare_error(statement) :
            extracted ? duckdb_extract_statements_error(extracted) : NULL;
        duckvep_sql_set_error(error, error_size, message ? message : "calls query must be one SELECT statement");
    } else if (duckdb_prepared_statement_type(statement) != DUCKDB_STATEMENT_TYPE_SELECT) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotype_arrangements: calls query must be one SELECT statement");
        ok = 0;
    } else {
        state->input.have_result = 1;
        state->input.source_records = 0;
        ok = duckdb_execute_prepared(statement, &state->input.input) == DuckDBSuccess;
        if (!ok) duckvep_sql_set_error(error, error_size, duckdb_result_error(&state->input.input));
    }
    if (statement) duckdb_destroy_prepare(&statement);
    if (extracted) duckdb_destroy_extracted(&extracted);
    pthread_mutex_unlock(&bind->registry->query_mutex);
    return ok;
}

static void arrangement_init(duckdb_init_info info) {
    const arrangement_bind_t *bind = duckdb_init_get_bind_data(info);
    arrangement_state_t *state = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1u, sizeof(*state));
    char error[DUCKVEP_SQL_ERROR_SIZE] = "duckvep_haplotype_arrangements: bounded workspace allocation failed";
    duckvep_hap_input_t input;
    duckdb_init_set_max_threads(info, 1u);
    duckvep_budget_clear_failure();
    if (!state) goto failed;
    state->core = duckvep_arrangement_open(&bind->config, error, sizeof(error));
    if (!state->core || !arrangement_input_open(state, bind, error, sizeof(error))) goto failed;
    input.context = &state->input; input.next = input_next;
    if (!duckvep_arrangement_load(state->core, &input, error, sizeof(error))) goto failed;
    if (state->input.chunk) duckdb_destroy_data_chunk(&state->input.chunk);
    if (state->input.have_result) { duckdb_destroy_result(&state->input.input); state->input.have_result = 0; }
    duckdb_init_set_init_data(info, state, arrangement_state_destroy);
    return;
failed:
    { char final_error[DUCKVEP_SQL_ERROR_SIZE + 256];
      duckdb_init_set_error(info, duckvep_sql_final_error(final_error, sizeof final_error, error, error)); }
    arrangement_state_destroy(state);
}

static void arrangement_scan(duckdb_function_info info, duckdb_data_chunk output) {
    arrangement_state_t *state = duckdb_function_get_init_data(info);
    duckdb_vector vectors[DUCKVEP_ARRANGEMENT_COLUMN_COUNT];
    idx_t rows = 0u;
    for (unsigned i = 0u; i < DUCKVEP_ARRANGEMENT_COLUMN_COUNT; i++) {
        vectors[i] = duckdb_data_chunk_get_vector(output, i);
        duckdb_vector_ensure_validity_writable(vectors[i]);
    }
    while (rows < duckdb_vector_size()) {
        duckvep_arrangement_row_t row;
        if (!duckvep_arrangement_next(state->core, &row)) break;
        ((uint64_t *)duckdb_vector_get_data(vectors[0]))[rows] = row.hypothesis_id;
        ((uint16_t *)duckdb_vector_get_data(vectors[1]))[rows] = row.lane;
        ((bool *)duckdb_vector_get_data(vectors[2]))[rows] = row.reference_lane != 0;
        ((uint64_t *)duckdb_vector_get_data(vectors[3]))[rows] = row.event_id;
        ((uint32_t *)duckdb_vector_get_data(vectors[4]))[rows] = row.chrom;
        ((uint64_t *)duckdb_vector_get_data(vectors[5]))[rows] = row.position;
        duckdb_vector_assign_string_element_len(vectors[6], rows, (const char *)row.reference, row.reference_length);
        duckdb_vector_assign_string_element_len(vectors[7], rows, (const char *)row.alternate, row.alternate_length);
        ((uint32_t *)duckdb_vector_get_data(vectors[8]))[rows] = row.allele_index;
        ((uint32_t *)duckdb_vector_get_data(vectors[9]))[rows] = row.transcript;
        ((uint32_t *)duckdb_vector_get_data(vectors[10]))[rows] = row.sample;
        ((int32_t *)duckdb_vector_get_data(vectors[11]))[rows] = row.allele0;
        ((int32_t *)duckdb_vector_get_data(vectors[12]))[rows] = row.allele1;
        ((bool *)duckdb_vector_get_data(vectors[13]))[rows] = row.phase0 != 0;
        ((bool *)duckdb_vector_get_data(vectors[14]))[rows] = row.phase1 != 0;
        ((bool *)duckdb_vector_get_data(vectors[15]))[rows] = row.phase_present != 0;
        ((bool *)duckdb_vector_get_data(vectors[16]))[rows] = row.phase_set_present != 0;
        if (row.phase_set_present) ((int64_t *)duckdb_vector_get_data(vectors[17]))[rows] = row.phase_set;
        else duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vectors[17]), rows);
        ((int32_t *)duckdb_vector_get_data(vectors[18]))[rows] = row.assigned_allele;
        ((bool *)duckdb_vector_get_data(vectors[19]))[rows] = row.contributes != 0;
        duckdb_vector_assign_string_element_len(vectors[20], rows, (const char *)row.cds, row.cds_length);
        duckdb_vector_assign_string_element_len(vectors[21], rows, (const char *)row.protein, row.protein_length);
        ((int64_t *)duckdb_vector_get_data(vectors[22]))[rows] = row.nominal_length_diff;
        duckdb_vector_assign_string_element(vectors[23], rows, duckvep_arrangement_prediction_status(row.prediction_status));
        ((uint64_t *)duckdb_vector_get_data(vectors[24]))[rows] = row.consequence_mask;
        duckdb_vector_assign_string_element(vectors[25], rows,
            row.hypothetical ? "hypothetical_assignment" : "reference_replay");
        rows++;
    }
    duckdb_data_chunk_set_size(output, rows);
}

static void register_arrangements(duckdb_connection connection, duckvep_registry_t *registry) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type string = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_table_function_set_name(function, "duckvep_haplotype_arrangements");
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_named_parameter(function, "max_sites", integer);
    duckdb_table_function_add_named_parameter(function, "max_calls", integer);
    duckdb_table_function_add_named_parameter(function, "max_arrangements", integer);
    duckdb_table_function_add_named_parameter(function, "max_replays", integer);
    duckvep_registry_retain(registry);
    duckdb_table_function_set_extra_info(function, registry, duckvep_registry_release);
    duckdb_table_function_set_bind(function, arrangement_bind);
    duckdb_table_function_set_init(function, arrangement_init);
    duckdb_table_function_set_function(function, arrangement_scan);
    (void)duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&integer); duckdb_destroy_logical_type(&string);
}

void duckvep_register_haplotypes(duckdb_connection connection, duckvep_registry_t *registry) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type string = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_set_name(function, "duckvep_haplotypes");
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_named_parameter(function, "phase_policy", string);
    duckdb_table_function_add_named_parameter(function, "input_mode", string);
    duckdb_table_function_add_named_parameter(function, "hgvs", boolean);
    for (unsigned i = 0u; i < LIMIT_COUNT; i++)
        duckdb_table_function_add_named_parameter(function, limit_names[i], integer);
    duckvep_registry_retain(registry);
    duckdb_table_function_set_extra_info(function, registry, duckvep_registry_release);
    duckdb_table_function_set_bind(function, haplotype_bind);
    duckdb_table_function_set_init(function, haplotype_init);
    duckdb_table_function_set_function(function, haplotype_scan);
    (void)duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&string); duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&boolean);
    register_arrangements(connection, registry);
    register_coding_transcripts(connection, registry);
    duckvep_register_coding_calls(connection, registry);
}
