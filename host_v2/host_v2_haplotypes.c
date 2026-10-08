/* The haplotype functions of the v2 host: duckvep_haplotype_load_sql (the statements to run), the job options of the
 * COPY format duckvep_stage, duckvep_haplotype_scan (the native replay over the captured input), the record-plan scan
 * _duckvep_haplotype_plan, duckvep_haplotype_drop and duckvep_coding_transcripts.
 *
 * v1 runs the caller's calls query on a private connection and normalizes it with temporary tables. v2 cannot, so the
 * normalization is SQL the caller runs (duckvep_haplotype_load_sql wraps the caller's query), and the normalized rows are
 * captured by COPY into spillable column collections. The scan replays them through src/core/duckvep_core_haplotypes.c,
 * the same code the v1 host runs, a full chunk of rows per exec call. A job is scanned once and is released when the
 * scan ends or fails. */
#include "duckvep_host.h"

#include "host_v2_columns.h"
#include "host_v2_stage.h"

#include "core/duckvep_core_arrangements.h"
#include "core/duckvep_core_discovery.h"
#include "core/duckvep_core_haplotype_script.h"
#include "core/duckvep_core_haplotypes.h"
#include "core/duckvep_core_phase.h"
#include "kernel/src/duckvep_budget.h"
#include "kernel/src/duckvep_haplotype_stream.h"

/* ---------------------------------------------------------------------------
 * The COPY options of a job: JOB, STAGE, PHASE_POLICY, INPUT_MODE, HGVS and the limits
 * ------------------------------------------------------------------------- */

struct hap_bind {
    char *job;
    char *stage;
    char *policy;
    char *mode;
    bool has_hgvs, hgvs;
    bool has_limit[DUCKVEP_HAP_LIMIT_COUNT];
    uint64_t limit[DUCKVEP_HAP_LIMIT_COUNT];
    bool has_arrangement_limit[DUCKVEP_ARRANGEMENT_LIMIT_COUNT];
    uint64_t arrangement_limit[DUCKVEP_ARRANGEMENT_LIMIT_COUNT];
};

static bool same_name(const duckdb_v2_identifier_t *name, const char *text) {
    size_t length = strlen(text);
    if (name->len != length) {
        return false;
    }
    for (size_t i = 0; i < length; ++i) {
        char c = name->ptr[i];
        if (c >= 'A' && c <= 'Z') {
            c = (char)(c + ('a' - 'A'));
        }
        if (c != text[i]) {
            return false;
        }
    }
    return true;
}

static char *text_copy(const char *text, size_t length) {
    char *copy = malloc(length + 1);
    if (copy) {
        memcpy(copy, text, length);
        copy[length] = '\0';
    }
    return copy;
}

hap_bind *host_v2_hap_bind_create(void) {
    return calloc(1, sizeof(hap_bind));
}

void host_v2_hap_bind_destroy(hap_bind *bind) {
    if (bind) {
        free(bind->job);
        free(bind->stage);
        free(bind->policy);
        free(bind->mode);
        free(bind);
    }
}

bool host_v2_hap_bind_is_job(const hap_bind *bind) {
    return bind && bind->job != NULL;
}

/* A non-negative integer option value of any integer type. */
static bool value_to_u64(duckdb_v2_value_handle value, uint64_t *out, duckdb_v2_error_info_handle *detail) {
    duckdb_v2_logical_type_handle type = NULL;
    DUCKDB_V2_LOGICAL_TYPE_ID id = DUCKDB_V2_LOGICAL_TYPE_ID_INVALID;
    bool ok = false;
    int64_t signed_value = 0;
    if (duckdb_v2_value_get_logical_type(value, &type, detail) || duckdb_v2_logical_type_get_id(type, &id, detail)) {
        goto done;
    }
    switch (id) {
    case DUCKDB_V2_LOGICAL_TYPE_ID_TINYINT: {
        int8_t v;
        ok = !duckdb_v2_value_get_tinyint(value, &v, detail);
        signed_value = v;
        break;
    }
    case DUCKDB_V2_LOGICAL_TYPE_ID_SMALLINT: {
        int16_t v;
        ok = !duckdb_v2_value_get_smallint(value, &v, detail);
        signed_value = v;
        break;
    }
    case DUCKDB_V2_LOGICAL_TYPE_ID_INTEGER: {
        int32_t v;
        ok = !duckdb_v2_value_get_int(value, &v, detail);
        signed_value = v;
        break;
    }
    case DUCKDB_V2_LOGICAL_TYPE_ID_BIGINT:
        ok = !duckdb_v2_value_get_bigint(value, &signed_value, detail);
        break;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UTINYINT: {
        uint8_t v;
        ok = !duckdb_v2_value_get_utinyint(value, &v, detail);
        *out = v;
        goto done;
    }
    case DUCKDB_V2_LOGICAL_TYPE_ID_USMALLINT: {
        uint16_t v;
        ok = !duckdb_v2_value_get_usmallint(value, &v, detail);
        *out = v;
        goto done;
    }
    case DUCKDB_V2_LOGICAL_TYPE_ID_UINTEGER: {
        uint32_t v;
        ok = !duckdb_v2_value_get_uint(value, &v, detail);
        *out = v;
        goto done;
    }
    case DUCKDB_V2_LOGICAL_TYPE_ID_UBIGINT:
        ok = !duckdb_v2_value_get_ubigint(value, out, detail);
        goto done;
    default:
        break;
    }
    if (ok && signed_value >= 0) {
        *out = (uint64_t)signed_value;
    } else {
        ok = false;
    }
done:
    (void)duckdb_v2_logical_type_destroy(&type);
    return ok;
}

static bool value_to_text(duckdb_v2_value_handle value, char **out, duckdb_v2_error_info_handle *detail) {
    duckdb_v2_logical_type_handle type = NULL;
    DUCKDB_V2_LOGICAL_TYPE_ID id = DUCKDB_V2_LOGICAL_TYPE_ID_INVALID;
    duckdb_v2_str text;
    bool ok = false;
    DUCKDB_V2_ERROR a = duckdb_v2_value_get_logical_type(value, &type, detail);
    DUCKDB_V2_ERROR b = a ? a : duckdb_v2_logical_type_get_id(type, &id, detail);
    if (a || b || id != DUCKDB_V2_LOGICAL_TYPE_ID_VARCHAR || duckdb_v2_value_get_varchar(value, &text, detail)) {
        goto done;
    }
    free(*out);
    *out = text_copy(text.ptr, text.len);
    ok = *out != NULL;
done:
    (void)duckdb_v2_logical_type_destroy(&type);
    return ok;
}

bool host_v2_hap_bind_option(hap_bind *bind, const duckdb_v2_identifier_t *name, duckdb_v2_value_handle value,
                             bool *consumed, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    char message[160];
    bool is_null = false;
    bool ok = false;
    char **target = NULL;
    const char *label = NULL;
    *consumed = false;
    if (same_name(name, "job")) { target = &bind->job; label = "JOB"; }
    else if (same_name(name, "stage")) { target = &bind->stage; label = "STAGE"; }
    else if (same_name(name, "phase_policy")) { target = &bind->policy; label = "PHASE_POLICY"; }
    else if (same_name(name, "input_mode")) { target = &bind->mode; label = "INPUT_MODE"; }
    if (target) {
        *consumed = true;
        if (duckdb_v2_value_is_null(value, &is_null, &detail) || is_null || !value_to_text(value, target, &detail)) {
            (void)snprintf(message, sizeof message, "duckvep_haplotypes: option %s must be a non-NULL string", label);
            INPUT_ERROR(message);
        }
        ok = true;
        goto cleanup;
    }
    if (same_name(name, "hgvs")) {
        bool flag = false;
        *consumed = true;
        if (duckdb_v2_value_is_null(value, &is_null, &detail) || is_null ||
            duckdb_v2_value_get_bool(value, &flag, &detail)) {
            INPUT_ERROR("duckvep_haplotypes: option HGVS must be a non-NULL boolean");
        }
        bind->has_hgvs = true;
        bind->hgvs = flag;
        ok = true;
        goto cleanup;
    }
    for (unsigned i = 0; i < DUCKVEP_HAP_LIMIT_COUNT; ++i) {
        if (same_name(name, duckvep_hap_limit_names[i])) {
            uint64_t number = 0;
            *consumed = true;
            if (duckdb_v2_value_is_null(value, &is_null, &detail) || is_null ||
                !value_to_u64(value, &number, &detail) || !duckvep_hap_limit_valid(i, number)) {
                (void)snprintf(message, sizeof message, "duckvep_haplotypes: invalid %s", duckvep_hap_limit_names[i]);
                INPUT_ERROR(message);
            }
            bind->has_limit[i] = true;
            bind->limit[i] = number;
            ok = true;
            goto cleanup;
        }
    }
    for (unsigned i = 0; i < DUCKVEP_ARRANGEMENT_LIMIT_COUNT; ++i) {
        if (same_name(name, duckvep_arrangement_limit_names[i])) {
            uint64_t number = 0;
            *consumed = true;
            if (duckdb_v2_value_is_null(value, &is_null, &detail) || is_null || !value_to_u64(value, &number, &detail)) {
                (void)snprintf(message, sizeof message, "duckvep_haplotype_arrangements: invalid %s",
                               duckvep_arrangement_limit_names[i]);
                INPUT_ERROR(message);
            }
            bind->has_arrangement_limit[i] = true;
            bind->arrangement_limit[i] = number;
            ok = true;
            goto cleanup;
        }
    }
    ok = true; /* not a job option */
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return ok;
}

/* Whether the staged column's type is `expected` (SQL text); a NULL expectation accepts a STRUCT. */
static bool type_matches(duckdb_v2_logical_type_handle type, const char *expected) {
    char text[512];
    idx_t length = 0;
    duckdb_v2_error_info_handle detail = NULL;
    bool ok = false;
    if (duckdb_v2_logical_type_to_text(type, text, sizeof text, &length, &detail) || length >= sizeof text) {
        goto done;
    }
    text[length] = '\0';
    ok = expected ? strcmp(text, expected) == 0 : strncmp(text, "STRUCT(", 7) == 0;
done:
    (void)duckdb_v2_error_info_destroy(&detail);
    return ok;
}

bool host_v2_hap_bind_finish(hap_bind *bind, model_state *state, duckdb_v2_context_handle context, stage *layout,
                             const char *model, duckdb_v2_error_info_handle *error) {
    (void)context;
    char message[DUCKVEP_SQL_ERROR_SIZE + 128];
    duckvep_model_entry_t *entry = NULL;
    duckvep_hap_config_t config;
    duckvep_arrangement_config_t arrangement_config;
    bool plan = false, arrangements = false, ok = false, policy_valid = true;
    if (!*bind->job || !model || !*model) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                  "duckvep_haplotypes: the JOB and MODEL options are required and must be non-empty");
        return false;
    }
    if (bind->stage && strcmp(bind->stage, "plan_input") != 0 && strcmp(bind->stage, "calls") != 0 &&
        strcmp(bind->stage, "arrangements") != 0) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                  "duckvep_haplotypes: STAGE must be 'calls', 'plan_input' or 'arrangements'");
        return false;
    }
    plan = bind->stage && strcmp(bind->stage, "plan_input") == 0;
    arrangements = bind->stage && strcmp(bind->stage, "arrangements") == 0;
    entry = duckvep_registry_pin(state->registry, model);
    if (!entry) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, "duckvep_haplotypes: unknown model name");
        return false;
    }
    if (arrangements) {
        size_t *limits[DUCKVEP_ARRANGEMENT_LIMIT_COUNT] = {
            &arrangement_config.max_sites, &arrangement_config.max_calls,
            &arrangement_config.max_arrangements, &arrangement_config.max_replays
        };
        if (bind->policy || bind->mode || bind->has_hgvs) {
            set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                      "duckvep_haplotype_arrangements: phase_policy, input_mode and hgvs are not supported");
            goto cleanup;
        }
        for (unsigned i = 0; i < DUCKVEP_HAP_LIMIT_COUNT; ++i) {
            if (bind->has_limit[i]) {
                set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                          "duckvep_haplotype_arrangements: haplotype replay limits are not supported");
                goto cleanup;
            }
        }
        duckvep_arrangement_config_defaults(&arrangement_config, &entry->model);
        for (unsigned i = 0; i < DUCKVEP_ARRANGEMENT_LIMIT_COUNT; ++i) {
            if (bind->has_arrangement_limit[i]) {
                if (bind->arrangement_limit[i] > SIZE_MAX) {
                    (void)snprintf(message, sizeof message, "duckvep_haplotype_arrangements: invalid %s",
                                   duckvep_arrangement_limit_names[i]);
                    set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
                    goto cleanup;
                }
                *limits[i] = (size_t)bind->arrangement_limit[i];
            }
        }
        if (!duckvep_arrangement_config_check(&arrangement_config, message, sizeof message)) {
            set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
            goto cleanup;
        }
        if (layout->columns != duckvep_hap_input_columns(0)) {
            set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                      "duckvep_haplotype_arrangements: the staged query must return 15 columns (run the statements of "
                      "duckvep_haplotype_arrangements_load_sql)");
            goto cleanup;
        }
        for (unsigned i = 0; i < layout->columns; ++i) {
            const char *expected = duckvep_hap_input_type(0, i);
            if (!type_matches(layout->types[i], expected)) {
                (void)snprintf(message, sizeof message,
                               "duckvep_haplotype_arrangements: staged column %u must be %s (run the statements of "
                               "duckvep_haplotype_arrangements_load_sql)", i + 1, expected);
                set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
                goto cleanup;
            }
        }
        layout->model = text_copy(bind->job, strlen(bind->job));
        layout->relation = text_copy(HAP_RELATION_ARRANGEMENTS, strlen(HAP_RELATION_ARRANGEMENTS));
        layout->hap_model = text_copy(model, strlen(model));
        if (!layout->model || !layout->relation || !layout->hap_model) {
            set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
            goto cleanup;
        }
        arrangement_config.model = NULL;
        layout->arrangements = arrangement_config;
        ok = true;
        goto cleanup;
    }
    duckvep_hap_config_defaults(&config, &entry->model);
    config.policy = DUCKVEP_PHASE_STRICT;
    policy_valid = !bind->policy || duckvep_core_phase_policy(bind->policy, strlen(bind->policy), &config.policy);
    config.source_records = bind->mode && strcmp(bind->mode, "source_records") == 0;
    config.hgvs = bind->has_hgvs && bind->hgvs;
    if (!policy_valid ||
        (bind->mode && strcmp(bind->mode, "alt_events") && strcmp(bind->mode, "source_records"))) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                  !policy_valid ? "duckvep_haplotypes: phase_policy must be 'strict' or 'vep_compat'"
                      : "duckvep_haplotypes: input_mode must be 'alt_events' or 'source_records'; source_records "
                        "requires phase_policy='vep_compat'");
        goto cleanup;
    }
    for (unsigned i = 0; i < DUCKVEP_HAP_LIMIT_COUNT; ++i) {
        if (bind->has_limit[i]) {
            config.limits[i] = (size_t)bind->limit[i];
        }
    }
    if (!duckvep_hap_config_check(&config, message, sizeof message)) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
        goto cleanup;
    }
    /* The staged rows must be what duckvep_haplotype_load_sql produces. */
    if (plan) {
        if (layout->columns != 5) {
            set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID,
                      "duckvep_haplotypes: the record plan input must have 5 columns");
            goto cleanup;
        }
        for (unsigned i = 0; i < 5; ++i) {
            if (!type_matches(layout->types[i], duckvep_hap_plan_input_type(i))) {
                (void)snprintf(message, sizeof message,
                               "duckvep_haplotypes: record plan input column %u must be %s", i + 1,
                               duckvep_hap_plan_input_type(i));
                set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
                goto cleanup;
            }
        }
    } else {
        unsigned columns = duckvep_hap_input_columns(config.source_records);
        if (layout->columns != columns) {
            (void)snprintf(message, sizeof message,
                           "duckvep_haplotypes: the staged query must return %u columns (run the statements of "
                           "duckvep_haplotype_load_sql)", columns);
            set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
            goto cleanup;
        }
        for (unsigned i = 0; i < columns; ++i) {
            const char *expected = duckvep_hap_input_type(config.source_records, i);
            if (!type_matches(layout->types[i], expected)) {
                (void)snprintf(message, sizeof message,
                               "duckvep_haplotypes: staged column %u must be %s (run the statements of "
                               "duckvep_haplotype_load_sql)", i + 1, expected ? expected : "a STRUCT");
                set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
                goto cleanup;
            }
        }
    }
    layout->model = text_copy(bind->job, strlen(bind->job));
    layout->relation = text_copy(plan ? HAP_RELATION_PLAN : HAP_RELATION_CALLS,
                                 strlen(plan ? HAP_RELATION_PLAN : HAP_RELATION_CALLS));
    layout->hap_model = text_copy(model, strlen(model));
    if (!layout->model || !layout->relation || !layout->hap_model) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        goto cleanup;
    }
    config.model = NULL;
    layout->hap = config;
    ok = true;
cleanup:
    duckvep_registry_unpin(state->registry, entry);
    return ok;
}

/* ---------------------------------------------------------------------------
 * A scan over a staged job
 * ------------------------------------------------------------------------- */

typedef struct {
    stage *staged;
    duckdb_v2_column_data_collection_shared_scan_state_handle shared;
    duckdb_v2_column_data_collection_worker_scan_state_handle worker;
    duckdb_v2_data_chunk_handle chunk;
    idx_t rows, cursor;
    bool done;
} job_scan;

static void job_scan_close(job_scan *scan) {
    (void)duckdb_v2_column_data_collection_worker_scan_state_destroy(&scan->worker);
    (void)duckdb_v2_column_data_collection_shared_scan_state_destroy(&scan->shared);
    (void)duckdb_v2_data_chunk_destroy(&scan->chunk);
}

static bool job_scan_open(job_scan *scan, duckdb_v2_context_handle context, char *message, size_t size) {
    duckdb_v2_error_info_handle detail = NULL;
    stage *s = scan->staged;
    bool ok = false;
    if (!s->collection) {
        scan->done = true; /* no rows were staged */
        return true;
    }
    if (duckdb_v2_column_data_collection_shared_scan_state_create(s->collection, &scan->shared, &detail) ||
        duckdb_v2_column_data_collection_worker_scan_state_create(s->collection, &scan->worker, &detail) ||
        duckdb_v2_data_chunk_create_with_context(context, s->types, s->columns, &scan->chunk, &detail)) {
        duckdb_v2_str text = {"could not scan the staged rows", 30};
        if (detail) {
            (void)duckdb_v2_error_info_get_text(detail, &text);
        }
        (void)snprintf(message, size, "%.*s", (int)text.len, text.ptr);
        job_scan_close(scan);
        goto done;
    }
    ok = true;
done:
    (void)duckdb_v2_error_info_destroy(&detail);
    return ok;
}

/* Fetches the next captured chunk; false at the end (or on failure, `message` set). */
static bool job_scan_chunk(job_scan *scan, char *message, size_t size, bool *failed) {
    duckdb_v2_error_info_handle detail = NULL;
    bool produced = false;
    idx_t rows = 0;
    *failed = false;
    for (;;) {
        if (scan->done || !scan->chunk) {
            goto finished;
        }
        if (duckdb_v2_column_data_collection_scan(scan->staged->collection, scan->shared, scan->worker, scan->chunk,
                                                  &produced, &detail) ||
            (produced && duckdb_v2_data_chunk_get_size(scan->chunk, &rows, &detail))) {
            duckdb_v2_str text = {"could not scan the staged rows", 30};
            if (detail) {
                (void)duckdb_v2_error_info_get_text(detail, &text);
            }
            (void)snprintf(message, size, "%.*s", (int)text.len, text.ptr);
            *failed = true;
            goto finished;
        }
        if (!produced) {
            scan->done = true;
            goto finished;
        }
        if (rows == 0) {
            continue;
        }
        scan->rows = rows;
        scan->cursor = 0;
        (void)duckdb_v2_error_info_destroy(&detail);
        return true;
    }
finished:
    (void)duckdb_v2_error_info_destroy(&detail);
    return false;
}

/* ---------------------------------------------------------------------------
 * duckvep_haplotype_scan(job)
 * ------------------------------------------------------------------------- */

typedef struct {
    char *job;
    int source_records;
} scan_bind;

typedef struct {
    model_state *state;
    duckvep_model_entry_t *entry;
    duckvep_hap_state_t *core;
    job_scan scan;
    int source_records;
    column cols[17];
    column gt, phase, sets, raw[7];
} scan_state;

static void scan_bind_destroy(void *pointer) {
    scan_bind *bind = pointer;
    if (bind) {
        free(bind->job);
        free(bind);
    }
}

static void scan_state_destroy(void *pointer) {
    scan_state *s = pointer;
    if (!s) {
        return;
    }
    duckvep_hap_close(s->core);
    job_scan_close(&s->scan);
    host_v2_stage_destroy(s->scan.staged);
    if (s->entry) {
        duckvep_registry_unpin(s->state->registry, s->entry);
    }
    host_v2_state_release(s->state);
    free(s);
}

static char *bind_argument(duckdb_v2_table_function_bind_info_handle info, idx_t index,
                           duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_value_handle value = NULL;
    duckdb_v2_str text;
    char *copy = NULL;
    bool is_null = false;
    DUCKDB_CALL(duckdb_v2_table_function_bind_get_arg_value(info, index, &value, &detail));
    DUCKDB_CALL(duckdb_v2_value_is_null(value, &is_null, &detail));
    if (!is_null) {
        DUCKDB_CALL(duckdb_v2_value_get_varchar(value, &text, &detail));
        copy = text_copy(text.ptr, text.len);
    }
cleanup:
    (void)duckdb_v2_value_destroy(&value);
    (void)duckdb_v2_error_info_destroy(&detail);
    return copy;
}

static void scan_bind_exec(duckdb_v2_table_function_bind_info_handle info, duckdb_v2_context_handle context,
                           duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    model_state *state = NULL;
    scan_bind *bind = calloc(1, sizeof(*bind));
    duckvep_hap_config_t config;
    char *model = NULL;
    char message[256];
    bool owned = false;
    if (!bind) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_scan: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_bind_get_user_data(info, (void **)&state, &detail));
    bind->job = bind_argument(info, 0, error);
    if (!bind->job || !*bind->job) {
        INPUT_ERROR("duckvep_haplotype_scan: job must be a non-empty string");
    }
    if (!host_v2_stage_peek_job(state, bind->job, HAP_RELATION_CALLS, &model, &config)) {
        (void)snprintf(message, sizeof message,
                       "duckvep_haplotype_scan: job '%s' is not staged (run the statements of "
                       "duckvep_haplotype_load_sql first; a job is scanned once)", bind->job);
        INPUT_ERROR(message);
    }
    bind->source_records = config.source_records;
    for (unsigned i = 0; i < DUCKVEP_HAP_OUTPUT_COLUMNS; ++i) {
        DUCKDB_CALL(make_type(context, duckvep_hap_column_type(i, bind->source_records), &type, &detail));
        DUCKDB_CALL(duckdb_v2_table_function_bind_add_result_column(
            info, string_view(duckvep_hap_column_name(i)), type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    {
        duckdb_v2_opaque data = {bind, scan_bind_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_bind_set_bind_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        scan_bind_destroy(bind);
    }
    free(model);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void scan_init_exec(duckdb_v2_table_function_init_global_info_handle info, duckdb_v2_context_handle context,
                           duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    scan_bind *bind = NULL;
    scan_state *s = calloc(1, sizeof(*s));
    duckvep_hap_config_t config;
    char message[DUCKVEP_SQL_ERROR_SIZE + 256];
    bool owned = false;
    if (!s) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_scan: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_init_global_get_bind_data(info, (void **)&bind, &detail));
    s->state = state;
    host_v2_state_retain(state);
    s->scan.staged = host_v2_stage_take(state, bind->job, HAP_RELATION_CALLS);
    if (!s->scan.staged) {
        (void)snprintf(message, sizeof message,
                       "duckvep_haplotype_scan: job '%s' is not staged (a job is scanned once)", bind->job);
        INPUT_ERROR(message);
    }
    s->source_records = s->scan.staged->hap.source_records;
    s->entry = duckvep_registry_pin(state->registry, s->scan.staged->hap_model);
    if (!s->entry) {
        INPUT_ERROR("duckvep_haplotypes: unknown model name");
    }
    config = s->scan.staged->hap;
    config.model = &s->entry->model;
    if (!duckvep_hap_config_check(&config, message, sizeof message)) {
        INPUT_ERROR(message);
    }
    s->core = duckvep_hap_open(&config, message, sizeof message);
    if (!s->core) {
        char final_message[DUCKVEP_SQL_ERROR_SIZE + 256];
        report(*error, duckvep_sql_final_error(final_message, sizeof final_message, message, message));
        goto cleanup;
    }
    if (!job_scan_open(&s->scan, context, message, sizeof message)) {
        INPUT_ERROR(message);
    }
    /* The scan position is unsynchronized state: one thread only. */
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_max_threads(info, 1, &detail));
    {
        duckdb_v2_opaque data = {s, scan_state_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_init_global_set_global_state(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        scan_state_destroy(s);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

#define AT(type, c, row) (((const type *)(c).view.data)[physical_row(&(c).view, (row))])

/* Refreshes the typed views of the chunk just fetched. */
static bool scan_load_chunk(scan_state *s, char *message, size_t size) {
    duckdb_v2_error_info_handle error_holder = NULL;
    duckdb_v2_error_info_handle *error = &error_holder;
    unsigned columns = duckvep_hap_input_columns(s->source_records);
    (void)size;
    error_holder = NULL;
    for (unsigned i = 0; i < columns; ++i) {
        duckdb_v2_vector_handle vector = NULL;
        duckdb_v2_error_info_handle detail = NULL;
        if (duckdb_v2_data_chunk_get_vector(s->scan.chunk, i, &vector, &detail) ||
            !column_open(vector, true, &s->cols[i], error)) {
            (void)snprintf(message, size, "duckvep_haplotype_scan: could not read staged column %u", i + 1);
            (void)duckdb_v2_error_info_destroy(&detail);
            return false;
        }
        (void)duckdb_v2_error_info_destroy(&detail);
    }
    if (s->source_records) {
        for (unsigned i = 0; i < 7; ++i) {
            if (!column_child(&s->cols[8], i, &s->raw[i], error)) {
                (void)snprintf(message, size, "duckvep_haplotype_scan: could not read the raw GT struct");
                return false;
            }
        }
    } else if (!column_child(&s->cols[8], 0, &s->gt, error) || !column_child(&s->cols[9], 0, &s->phase, error) ||
               !column_child(&s->cols[11], 0, &s->sets, error)) {
        (void)snprintf(message, size, "duckvep_haplotype_scan: could not read the staged lists");
        return false;
    }
    return true;
}

static int scan_input_next(void *context, duckvep_hap_row_t *out, char *error, size_t error_size) {
    scan_state *s = context;
    job_scan *scan = &s->scan;
    bool failed = false;
    idx_t row;
    unsigned columns = duckvep_hap_input_columns(s->source_records);
    if (scan->cursor >= scan->rows || !scan->rows) {
        if (!job_scan_chunk(scan, error, error_size, &failed)) {
            return failed ? -1 : 0;
        }
        if (!scan_load_chunk(s, error, error_size)) {
            return -1;
        }
    }
    row = scan->cursor++;
    memset(out, 0, sizeof *out);
    for (unsigned i = 0; i < columns; ++i) {
        if (!column_valid(&s->cols[i], row)) {
            out->null_mask |= UINT32_C(1) << i;
        }
    }
    if (out->null_mask & (((UINT32_C(1) << columns) - 1u) & ~((UINT32_C(1) << 9) | (UINT32_C(1) << 10)))) {
        return 1; /* the core names the first NULL column */
    }
    out->event_id = AT(uint64_t, s->cols[0], row);
    out->chrom = AT(uint32_t, s->cols[1], row);
    out->pos = AT(uint64_t, s->cols[2], row);
    {
        duckdb_v2_str ref = string_at(&s->cols[3].view, row), alt = string_at(&s->cols[4].view, row);
        out->ref = (const uint8_t *)ref.ptr;
        out->ref_len = (uint32_t)ref.len;
        out->alt = (const uint8_t *)alt.ptr;
        out->alt_len = (uint32_t)alt.len;
    }
    out->allele_index = AT(uint32_t, s->cols[5], row);
    out->transcript = AT(uint32_t, s->cols[6], row);
    out->sample = AT(uint32_t, s->cols[7], row);
    out->copies = AT(int64_t, s->cols[12], row);
    out->versions = AT(int64_t, s->cols[13], row);
    out->ploidies = AT(int64_t, s->cols[14], row);
    if (s->source_records) {
        for (unsigned i = 0; i < 7; ++i) {
            out->raw[i] = AT(uint32_t, s->raw[i], row);
        }
        out->replay_order = AT(uint64_t, s->cols[15], row);
        out->source_selected = AT(uint8_t, s->cols[16], row) != 0;
        return 1;
    }
    {
        duckdb_v2_list_entry gt = list_at(&s->cols[8], row), sets = list_at(&s->cols[11], row);
        out->gt = (const int32_t *)s->gt.view.data;
        out->gt_validity = s->gt.view.validity;
        out->gt_offset = gt.offset;
        out->gt_length = gt.length;
        out->have_phase = !(out->null_mask >> 9 & 1u);
        if (out->have_phase) {
            duckdb_v2_list_entry phases = list_at(&s->cols[9], row);
            out->phase = (const uint8_t *)s->phase.view.data;
            out->phase_validity = s->phase.view.validity;
            out->phase_offset = phases.offset;
            out->phase_length = phases.length;
        }
        out->phase_set_present = !(out->null_mask >> 10 & 1u);
        if (out->phase_set_present) {
            out->phase_set = AT(int64_t, s->cols[10], row);
        }
        out->sets = (const int64_t *)s->sets.view.data;
        out->sets_validity = s->sets.view.validity;
        out->sets_offset = sets.offset;
        out->sets_length = sets.length;
    }
    return 1;
}

static void scan_exec(duckdb_v2_table_function_exec_info_handle info, duckdb_v2_context_handle context,
                      duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    scan_state *s = NULL;
    duckdb_v2_data_chunk_handle chunk = NULL;
    v2_call *call = calloc(1, sizeof(*call));
    duckvep_hap_input_t input;
    char message[DUCKVEP_SQL_ERROR_SIZE + 256];
    const size_t capacity = duckvep_h_vector_size();
    size_t rows = 0;
    if (!call) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_scan: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_global_state(info, (void **)&s, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_output_chunk(info, &chunk, &detail));
    call->error = error;
    call->rows = capacity;
    call->argc = DUCKVEP_HAP_OUTPUT_COLUMNS;
    for (unsigned i = 0; i < DUCKVEP_HAP_OUTPUT_COLUMNS; ++i) {
        duckdb_v2_vector_handle vector = NULL;
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, i, &vector, &detail));
        if (!v2_open_writable(call, &call->inputs[i], vector, capacity, true)) {
            goto cleanup;
        }
    }
    input.context = s;
    input.next = scan_input_next;
    /* A full chunk of rows, as on v1: list children grow as rows are appended (duckvep_host.h). */
    if (!duckvep_hap_scan(s->core, &input, call, capacity, &rows, message, sizeof message)) {
        if (!call->failed) {
            report(*error, message);
        }
        goto cleanup;
    }
    v2_finish(call);
    if (call->failed) {
        goto cleanup;
    }
    if (rows != capacity) {
        /* The rows written; an empty batch ends the scan. */
        for (unsigned i = 0; i < DUCKVEP_HAP_OUTPUT_COLUMNS; ++i) {
            DUCKDB_CALL(duckdb_v2_vector_set_size(call->inputs[i].handle, rows, &detail));
        }
    }
cleanup:
    free(call);
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_haplotype_arrangements(job): replay a strict staged diploid input
 * ------------------------------------------------------------------------- */

typedef struct {
    scan_state input;
    duckvep_arrangement_state_t *core;
} arrangement_scan_state;

static void arrangement_scan_state_destroy(void *pointer) {
    arrangement_scan_state *s = pointer;
    if (!s) return;
    duckvep_arrangement_close(s->core);
    job_scan_close(&s->input.scan);
    host_v2_stage_destroy(s->input.scan.staged);
    if (s->input.entry) duckvep_registry_unpin(s->input.state->registry, s->input.entry);
    host_v2_state_release(s->input.state);
    free(s);
}

static void arrangement_bind_exec(duckdb_v2_table_function_bind_info_handle info, duckdb_v2_context_handle context,
                                  duckdb_v2_error_info_handle *error) {
    static const char *const names[DUCKVEP_ARRANGEMENT_COLUMN_COUNT] = {
        "hypothesis_id", "hypothesis_lane", "hypothesis_reference_lane", "event_index", "seq_region", "position",
        "reference", "alternate", "alt_index", "transcript_index", "sample_index", "original_allele0",
        "original_allele1", "original_phase_before0", "original_phase_before1", "original_phase_before_present",
        "original_phase_set_present", "original_phase_set", "assigned_allele", "contributes", "cds", "protein",
        "nominal_length_diff", "prediction_status", "consequence_mask", "prediction_semantics"
    };
    static const char *const types[DUCKVEP_ARRANGEMENT_COLUMN_COUNT] = {
        "UBIGINT", "USMALLINT", "BOOLEAN", "UBIGINT", "UINTEGER", "UBIGINT", "VARCHAR", "VARCHAR",
        "UINTEGER", "UINTEGER", "UINTEGER", "INTEGER", "INTEGER", "BOOLEAN", "BOOLEAN", "BOOLEAN",
        "BOOLEAN", "BIGINT", "INTEGER", "BOOLEAN", "VARCHAR", "VARCHAR", "BIGINT", "VARCHAR", "UBIGINT",
        "VARCHAR"
    };
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    model_state *state = NULL;
    scan_bind *bind = calloc(1, sizeof(*bind));
    stage *staged = NULL;
    char message[256];
    bool owned = false;
    if (!bind) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_arrangements: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_bind_get_user_data(info, (void **)&state, &detail));
    bind->job = bind_argument(info, 0, error);
    if (!bind->job || !*bind->job) INPUT_ERROR("duckvep_haplotype_arrangements: job must be a non-empty string");
    staged = host_v2_stage_acquire(state, bind->job, HAP_RELATION_ARRANGEMENTS);
    if (!staged) {
        (void)snprintf(message, sizeof message,
                       "duckvep_haplotype_arrangements: job '%s' is not staged (run the statements of "
                       "duckvep_haplotype_arrangements_load_sql first; a job is scanned once)", bind->job);
        INPUT_ERROR(message);
    }
    host_v2_stage_destroy(staged);
    staged = NULL;
    for (unsigned i = 0; i < DUCKVEP_ARRANGEMENT_COLUMN_COUNT; ++i) {
        DUCKDB_CALL(make_type(context, types[i], &type, &detail));
        DUCKDB_CALL(duckdb_v2_table_function_bind_add_result_column(info, string_view(names[i]), type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    {
        duckdb_v2_opaque data = {bind, scan_bind_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_bind_set_bind_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (staged) host_v2_stage_destroy(staged);
    if (!owned) scan_bind_destroy(bind);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void arrangement_init_exec(duckdb_v2_table_function_init_global_info_handle info, duckdb_v2_context_handle context,
                                  duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    scan_bind *bind = NULL;
    arrangement_scan_state *s = calloc(1, sizeof(*s));
    duckvep_hap_input_t input;
    duckvep_arrangement_config_t config;
    char message[DUCKVEP_SQL_ERROR_SIZE + 256];
    bool owned = false;
    if (!s) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_arrangements: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_init_global_get_bind_data(info, (void **)&bind, &detail));
    s->input.state = state;
    host_v2_state_retain(state);
    s->input.scan.staged = host_v2_stage_take(state, bind->job, HAP_RELATION_ARRANGEMENTS);
    if (!s->input.scan.staged) {
        (void)snprintf(message, sizeof message,
                       "duckvep_haplotype_arrangements: job '%s' is not staged (a job is scanned once)", bind->job);
        INPUT_ERROR(message);
    }
    s->input.entry = duckvep_registry_pin(state->registry, s->input.scan.staged->hap_model);
    if (!s->input.entry) INPUT_ERROR("duckvep_haplotype_arrangements: unknown model name");
    if (s->input.entry->model.lifted) {
        INPUT_ERROR("duckvep_haplotype_arrangements: replay is not supported for models with wrapped circular objects");
    }
    config = s->input.scan.staged->arrangements;
    config.model = &s->input.entry->model;
    s->core = duckvep_arrangement_open(&config, message, sizeof message);
    if (!s->core) {
        char final_message[DUCKVEP_SQL_ERROR_SIZE + 256];
        report(*error, duckvep_sql_final_error(final_message, sizeof final_message, message, message));
        goto cleanup;
    }
    s->input.source_records = false;
    if (!job_scan_open(&s->input.scan, context, message, sizeof message)) INPUT_ERROR(message);
    input.context = &s->input;
    input.next = scan_input_next;
    if (!duckvep_arrangement_load(s->core, &input, message, sizeof message)) {
        char final_message[DUCKVEP_SQL_ERROR_SIZE + 256];
        report(*error, duckvep_sql_final_error(final_message, sizeof final_message, message, message));
        goto cleanup;
    }
    job_scan_close(&s->input.scan);
    host_v2_stage_destroy(s->input.scan.staged);
    s->input.scan.staged = NULL;
    /* The arrangement engine has fully validated its charged input before it yields a row. */
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_max_threads(info, 1, &detail));
    {
        duckdb_v2_opaque data = {s, arrangement_scan_state_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_init_global_set_global_state(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) arrangement_scan_state_destroy(s);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void arrangement_exec(duckdb_v2_table_function_exec_info_handle info, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    arrangement_scan_state *s = NULL;
    duckdb_v2_data_chunk_handle chunk = NULL;
    v2_call *call = calloc(1, sizeof(*call));
    const size_t capacity = duckvep_h_vector_size();
    size_t rows = 0;
    if (!call) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_arrangements: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_global_state(info, (void **)&s, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_output_chunk(info, &chunk, &detail));
    call->error = error;
    call->rows = capacity;
    call->argc = DUCKVEP_ARRANGEMENT_COLUMN_COUNT;
    for (unsigned i = 0; i < DUCKVEP_ARRANGEMENT_COLUMN_COUNT; ++i) {
        duckdb_v2_vector_handle vector = NULL;
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, i, &vector, &detail));
        if (!v2_open_writable(call, &call->inputs[i], vector, capacity, true)) goto cleanup;
    }
    while (rows < capacity) {
        duckvep_arrangement_row_t row;
        v2_vec *v = call->inputs;
        if (!duckvep_arrangement_next(s->core, &row)) break;
        ((uint64_t *)v[0].data)[rows] = row.hypothesis_id;
        ((uint16_t *)v[1].data)[rows] = row.lane;
        ((bool *)v[2].data)[rows] = row.reference_lane != 0;
        ((uint64_t *)v[3].data)[rows] = row.event_id;
        ((uint32_t *)v[4].data)[rows] = row.chrom;
        ((uint64_t *)v[5].data)[rows] = row.position;
        duckvep_h_assign_string(&v[6], rows, (const char *)row.reference, row.reference_length);
        duckvep_h_assign_string(&v[7], rows, (const char *)row.alternate, row.alternate_length);
        ((uint32_t *)v[8].data)[rows] = row.allele_index;
        ((uint32_t *)v[9].data)[rows] = row.transcript;
        ((uint32_t *)v[10].data)[rows] = row.sample;
        ((int32_t *)v[11].data)[rows] = row.allele0;
        ((int32_t *)v[12].data)[rows] = row.allele1;
        ((bool *)v[13].data)[rows] = row.phase0 != 0;
        ((bool *)v[14].data)[rows] = row.phase1 != 0;
        ((bool *)v[15].data)[rows] = row.phase_present != 0;
        ((bool *)v[16].data)[rows] = row.phase_set_present != 0;
        if (row.phase_set_present) ((int64_t *)v[17].data)[rows] = row.phase_set;
        else mark_null(v[17].validity, rows);
        ((int32_t *)v[18].data)[rows] = row.assigned_allele;
        ((bool *)v[19].data)[rows] = row.contributes != 0;
        duckvep_h_assign_string(&v[20], rows, (const char *)row.cds, row.cds_length);
        duckvep_h_assign_string(&v[21], rows, (const char *)row.protein, row.protein_length);
        ((int64_t *)v[22].data)[rows] = row.nominal_length_diff;
        {
            const char *status = duckvep_arrangement_prediction_status(row.prediction_status);
            duckvep_h_assign_string(&v[23], rows, status, strlen(status));
        }
        ((uint64_t *)v[24].data)[rows] = row.consequence_mask;
        {
            const char *semantics = row.hypothetical ? "hypothetical_assignment" : "reference_replay";
            duckvep_h_assign_string(&v[25], rows, semantics, strlen(semantics));
        }
        if (call->failed) goto cleanup;
        ++rows;
    }
    if (rows != capacity) {
        for (unsigned i = 0; i < DUCKVEP_ARRANGEMENT_COLUMN_COUNT; ++i) {
            DUCKDB_CALL(duckdb_v2_vector_set_size(call->inputs[i].handle, rows, &detail));
        }
    }
cleanup:
    free(call);
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * _duckvep_haplotype_plan(job): the record plan of a source_records job
 * ------------------------------------------------------------------------- */

typedef struct {
    model_state *state;
    duckvep_model_entry_t *entry;
    duckvep_haplotype_record_plan_t plan;
    job_scan scan;
    uint64_t ordinal;
    column cols[5];
} plan_state;

static void plan_state_destroy(void *pointer) {
    plan_state *s = pointer;
    if (!s) {
        return;
    }
    job_scan_close(&s->scan);
    host_v2_stage_destroy(s->scan.staged);
    if (s->entry) {
        duckvep_registry_unpin(s->state->registry, s->entry);
    }
    host_v2_state_release(s->state);
    free(s);
}

static void plan_bind_exec(duckdb_v2_table_function_bind_info_handle info, duckdb_v2_context_handle context,
                           duckdb_v2_error_info_handle *error) {
    static const char *const names[4] = {"event_index", "buffer_id", "ordinal", "source_ordinal"};
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    model_state *state = NULL;
    duckvep_hap_config_t config;
    char *job = NULL, *model = NULL;
    char message[256];
    scan_bind *bind = calloc(1, sizeof(*bind));
    bool owned = false;
    if (!bind) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "_duckvep_haplotype_plan: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_bind_get_user_data(info, (void **)&state, &detail));
    job = bind_argument(info, 0, error);
    if (!job || !*job) {
        INPUT_ERROR("_duckvep_haplotype_plan: job must be a non-empty string");
    }
    bind->job = text_copy(job, strlen(job));
    if (!bind->job) {
        INPUT_ERROR("_duckvep_haplotype_plan: out of memory");
    }
    if (!host_v2_stage_peek_job(state, job, HAP_RELATION_PLAN, &model, &config)) {
        (void)snprintf(message, sizeof message,
                       "duckvep_haplotypes: the record plan of job '%s' is not staged (run the statements of "
                       "duckvep_haplotype_load_sql in order)", job);
        INPUT_ERROR(message);
    }
    for (unsigned i = 0; i < 4; ++i) {
        DUCKDB_CALL(make_type(context, "UBIGINT", &type, &detail));
        DUCKDB_CALL(duckdb_v2_table_function_bind_add_result_column(info, string_view(names[i]), type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    {
        duckdb_v2_opaque data = {bind, scan_bind_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_bind_set_bind_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        scan_bind_destroy(bind);
    }
    free(job);
    free(model);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void plan_init_exec(duckdb_v2_table_function_init_global_info_handle info, duckdb_v2_context_handle context,
                           duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    plan_state *s = calloc(1, sizeof(*s));
    char *job = NULL;
    char message[DUCKVEP_SQL_ERROR_SIZE + 256];
    bool owned = false;
    if (!s) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "_duckvep_haplotype_plan: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_get_user_data(info, (void **)&state, &detail));
    s->state = state;
    host_v2_state_retain(state);
    {
        void *bind_data = NULL;
        DUCKDB_CALL(duckdb_v2_table_function_init_global_get_bind_data(info, &bind_data, &detail));
        job = bind_data ? text_copy(((scan_bind *)bind_data)->job, strlen(((scan_bind *)bind_data)->job)) : NULL;
    }
    s->scan.staged = job ? host_v2_stage_acquire(state, job, HAP_RELATION_PLAN) : NULL;
    if (!s->scan.staged) {
        INPUT_ERROR("duckvep_haplotypes: the record plan is not staged (a job is scanned once)");
    }
    s->entry = duckvep_registry_pin(state->registry, s->scan.staged->hap_model);
    if (!s->entry) {
        INPUT_ERROR("duckvep_haplotypes: unknown model name");
    }
    if (s->entry->model.lifted) {
        INPUT_ERROR("duckvep_haplotypes: phased edit sets are not supported for models with wrapped circular objects");
    }
    if (!duckvep_haplotype_record_plan_init(&s->plan, &s->entry->model.transcripts)) {
        INPUT_ERROR("duckvep_haplotypes: invalid source-planning model");
    }
    if (!job_scan_open(&s->scan, context, message, sizeof message)) {
        INPUT_ERROR(message);
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_max_threads(info, 1, &detail));
    {
        duckdb_v2_opaque data = {s, plan_state_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_init_global_set_global_state(info, &data, &detail));
        owned = true;
    }
cleanup:
    free(job);
    if (!owned) {
        plan_state_destroy(s);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void plan_exec(duckdb_v2_table_function_exec_info_handle info, duckdb_v2_context_handle context,
                      duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    plan_state *s = NULL;
    duckdb_v2_data_chunk_handle chunk = NULL;
    duckdb_v2_vector_handle outputs[4] = {0};
    uint64_t *data[4] = {0};
    char message[256];
    bool failed = false;
    idx_t produced = 0;
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_global_state(info, (void **)&s, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_output_chunk(info, &chunk, &detail));
    for (unsigned i = 0; i < 4; ++i) {
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, i, &outputs[i], &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(outputs[i], (void **)&data[i], &detail));
    }
    if (job_scan_chunk(&s->scan, message, sizeof message, &failed)) {
        for (unsigned i = 0; i < 5; ++i) {
            duckdb_v2_vector_handle vector = NULL;
            DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(s->scan.chunk, i, &vector, &detail));
            if (!column_open(vector, true, &s->cols[i], error)) {
                goto cleanup;
            }
        }
        for (idx_t row = 0; row < s->scan.rows; ++row) {
            uint64_t buffer, ordinal;
            uint32_t length;
            uint64_t chrom, pos;
            duckdb_v2_str ref;
            for (unsigned i = 0; i < 4; ++i) {
                if (!column_valid(&s->cols[i], row)) {
                    INPUT_ERROR("duckvep_haplotypes: required input source column is NULL");
                }
            }
            if (!column_valid(&s->cols[4], row) || AT(int64_t, s->cols[4], row) != 1) {
                INPUT_ERROR("duckvep_haplotypes: inconsistent source record identity");
            }
            chrom = AT(uint32_t, s->cols[1], row);
            pos = AT(uint64_t, s->cols[2], row);
            ref = string_at(&s->cols[3].view, row);
            length = (uint32_t)ref.len;
            if (chrom > UINT16_MAX || !pos || pos > UINT32_MAX || !length || length > UINT16_MAX ||
                length - 1u > UINT32_MAX - pos || s->ordinal == UINT64_MAX ||
                !duckvep_haplotype_record_plan_next(&s->plan, (uint16_t)chrom, (uint32_t)pos,
                                                    (uint32_t)(pos + length - 1u), &buffer, &ordinal)) {
                INPUT_ERROR("duckvep_haplotypes: invalid source span or record count");
            }
            s->ordinal++;
            data[0][row] = AT(uint64_t, s->cols[0], row);
            data[1][row] = buffer;
            data[2][row] = ordinal;
            data[3][row] = s->ordinal;
        }
        produced = s->scan.rows;
    } else if (failed) {
        INPUT_ERROR(message);
    }
    for (unsigned i = 0; i < 4; ++i) {
        uint64_t *validity = NULL;
        DUCKDB_CALL(duckdb_v2_vector_set_size(outputs[i], produced, &detail));
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(outputs[i], &validity, &detail));
        for (idx_t row = 0; row < produced; ++row) {
            mark_valid(validity, row);
        }
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_haplotype_load_sql(calls_query, model, job [, options]) -> VARCHAR[]
 * ------------------------------------------------------------------------- */

#define SCRIPT_OPTIONS (2 + 1 + DUCKVEP_HAP_LIMIT_COUNT)

typedef struct {
    duckvep_sql_text statements[DUCKVEP_HAPLOTYPE_SCRIPT_MAX];
    size_t count;
} hap_script;

static char *cell_text(const duckvep_cell_t *cell) {
    return duckvep_core_string_copy(cell->text, cell->text_length);
}

static void load_sql_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                          duckdb_v2_error_info_handle *error) {
    (void)context;
    static const char *keys[SCRIPT_OPTIONS];
    static duckvep_core_option_kind_t kinds[SCRIPT_OPTIONS];
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle arguments[4] = {0};
    column strings[3];
    options option;
    duckdb_v2_vector_handle output = NULL, child = NULL;
    void *list_data = NULL, *child_data = NULL;
    uint64_t *list_validity = NULL;
    duckdb_v2_arena_handle arena = NULL;
    idx_t rows = 0, at = 0;
    uint32_t argc = 0;
    hap_script *scripts = NULL;
    size_t total = 0;
    keys[0] = "phase_policy";
    keys[1] = "input_mode";
    keys[2] = "hgvs";
    kinds[0] = kinds[1] = DUCKVEP_CORE_OPTION_TEXT;
    kinds[2] = DUCKVEP_CORE_OPTION_BOOLEAN;
    for (unsigned i = 0; i < DUCKVEP_HAP_LIMIT_COUNT; ++i) {
        keys[3 + i] = duckvep_hap_limit_names[i];
        kinds[3 + i] = DUCKVEP_CORE_OPTION_INTEGER;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    for (uint32_t i = 0; i < argc; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
    }
    for (uint32_t i = 0; i < 3; ++i) {
        if (!column_open(arguments[i], false, &strings[i], error)) {
            goto cleanup;
        }
    }
    if (!options_open(argc == 4 ? arguments[3] : NULL, keys, kinds, SCRIPT_OPTIONS, &option, error)) {
        goto cleanup;
    }
    scripts = calloc(rows ? rows : 1, sizeof(*scripts));
    if (!scripts) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_load_sql: out of memory");
        goto cleanup;
    }
    for (idx_t row = 0; row < rows; ++row) {
        char *values[3] = {0};
        char *policy = NULL, *mode = NULL;
        duckvep_haplotype_script_options_t settings;
        duckvep_cell_t cell;
        char message[256];
        const char *failure = NULL;
        int built = 0;
        memset(&settings, 0, sizeof settings);
        settings.hgvs = -1;
        if (!options_check(&option, row, error)) {
            goto cleanup;
        }
        for (uint32_t i = 0; i < 3; ++i) {
            if (column_valid(&strings[i], row)) {
                duckdb_v2_str text = string_at(&strings[i].view, row);
                values[i] = duckvep_core_string_copy(text.ptr, text.len);
            }
            if (!values[i] || !*values[i]) {
                failure = "duckvep_haplotype_load_sql: the calls query, model and job must be non-empty strings";
            }
        }
        if (!failure && option_cell(&option, 0, row, &cell) && cell.valid) {
            policy = cell_text(&cell);
            settings.phase_policy = policy;
        }
        if (!failure && option_cell(&option, 1, row, &cell) && cell.valid) {
            mode = cell_text(&cell);
            settings.input_mode = mode;
        }
        if (!failure && option_cell(&option, 2, row, &cell)) {
            if (!cell.valid) {
                failure = "duckvep_haplotype_load_sql: hgvs cannot be NULL";
            } else {
                settings.hgvs = cell.boolean ? 1 : 0;
            }
        }
        for (unsigned i = 0; !failure && i < DUCKVEP_HAP_LIMIT_COUNT; ++i) {
            if (option_cell(&option, 3 + i, row, &cell)) {
                bool negative = cell.kind == DUCKVEP_CELL_TINYINT || cell.kind == DUCKVEP_CELL_SMALLINT ||
                                cell.kind == DUCKVEP_CELL_INTEGER || cell.kind == DUCKVEP_CELL_BIGINT
                                    ? cell.i < 0 : false;
                if (!cell.valid || negative) {
                    (void)snprintf(message, sizeof message, "duckvep_haplotype_load_sql: invalid %s",
                                   duckvep_hap_limit_names[i]);
                    failure = message;
                } else {
                    settings.has_limit[i] = 1;
                    settings.limit[i] = cell.kind == DUCKVEP_CELL_TINYINT || cell.kind == DUCKVEP_CELL_SMALLINT ||
                                        cell.kind == DUCKVEP_CELL_INTEGER || cell.kind == DUCKVEP_CELL_BIGINT
                                            ? (uint64_t)cell.i : cell.u;
                }
            }
        }
        if (!failure) {
            built = duckvep_core_haplotype_script(values[0], values[1], values[2], &settings, scripts[row].statements,
                                                  &scripts[row].count, message, sizeof message);
            if (built == 1) {
                failure = "duckvep_haplotype_load_sql: out of memory";
            } else if (built == 2) {
                failure = message;
            }
        }
        if (failure) {
            report(*error, failure);
        }
        for (size_t i = 0; i < 3; ++i) {
            duckvep_budget_free(values[i]);
        }
        duckvep_budget_free(policy);
        duckvep_budget_free(mode);
        if (failure) {
            goto cleanup;
        }
        total += scripts[row].count;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &list_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &list_validity, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_child(output, 0, &child, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(child, total, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(child, &child_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(child, &arena, &detail));
    {
        uint64_t *child_validity = NULL;
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(child, &child_validity, &detail));
        for (size_t i = 0; i < total; ++i) {
            mark_valid(child_validity, i);
        }
    }
    for (idx_t row = 0; row < rows; ++row) {
        ((duckdb_v2_list_entry *)list_data)[row] = (duckdb_v2_list_entry){at, scripts[row].count};
        mark_valid(list_validity, row);
        for (size_t i = 0; i < scripts[row].count; ++i, ++at) {
            DUCKDB_CALL(write_string(arena, &((duckdb_v2_bytes *)child_data)[at], scripts[row].statements[i].data,
                                     scripts[row].statements[i].length, &detail));
        }
    }
cleanup:
    if (scripts) {
        for (idx_t row = 0; row < rows; ++row) {
            for (size_t i = 0; i < scripts[row].count; ++i) {
                duckvep_sql_free(&scripts[row].statements[i]);
            }
        }
        free(scripts);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_haplotype_arrangements_load_sql(calls_query, model, job [, options])
 * ------------------------------------------------------------------------- */

#define ARRANGEMENT_SCRIPT_OPTIONS DUCKVEP_ARRANGEMENT_LIMIT_COUNT

typedef struct {
    bool has_limit[DUCKVEP_ARRANGEMENT_LIMIT_COUNT];
    uint64_t limit[DUCKVEP_ARRANGEMENT_LIMIT_COUNT];
} arrangement_script_options;

static int arrangement_script(const char *query, const char *model, const char *job,
                              const arrangement_script_options *settings, duckvep_sql_text *out,
                              char *message, size_t message_size) {
    duckvep_sql_text input = {0};
    char number[48];
    bool ok;
    int code = duckvep_core_arrangement_input_sql(query, &input, message, message_size);
    if (code) return code;
    ok = duckvep_sql_append(out, "COPY (\n") && duckvep_sql_append(out, input.data) &&
        duckvep_sql_append(out, "\n) TO 'duckvep_stage' (FORMAT duckvep_stage, JOB ") &&
        duckvep_sql_literal(out, job) && duckvep_sql_append(out, ", MODEL ") && duckvep_sql_literal(out, model) &&
        duckvep_sql_append(out, ", STAGE 'arrangements'");
    for (unsigned i = 0; ok && i < DUCKVEP_ARRANGEMENT_LIMIT_COUNT; ++i) {
        if (!settings->has_limit[i]) continue;
        (void)snprintf(number, sizeof number, " %llu", (unsigned long long)settings->limit[i]);
        ok = duckvep_sql_append(out, ", ") && duckvep_sql_append(out, duckvep_arrangement_limit_names[i]) &&
            duckvep_sql_append(out, number);
    }
    ok = ok && duckvep_sql_append(out, ", USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)");
    duckvep_sql_free(&input);
    if (!ok) {
        duckvep_sql_free(out);
        (void)snprintf(message, message_size, "%s", "duckvep_haplotype_arrangements_load_sql: out of memory");
        return 1;
    }
    return 0;
}

static void arrangement_load_sql_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                                      duckdb_v2_error_info_handle *error) {
    const char *keys[ARRANGEMENT_SCRIPT_OPTIONS];
    duckvep_core_option_kind_t kinds[ARRANGEMENT_SCRIPT_OPTIONS];
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle arguments[4] = {0};
    column strings[3];
    options option;
    duckdb_v2_vector_handle output = NULL, child = NULL;
    void *list_data = NULL, *child_data = NULL;
    uint64_t *list_validity = NULL;
    duckdb_v2_arena_handle arena = NULL;
    idx_t rows = 0;
    uint32_t argc = 0;
    hap_script *scripts = NULL;
    (void)context;
    for (unsigned i = 0; i < ARRANGEMENT_SCRIPT_OPTIONS; ++i) {
        keys[i] = duckvep_arrangement_limit_names[i];
        kinds[i] = DUCKVEP_CORE_OPTION_INTEGER;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    for (uint32_t i = 0; i < argc; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
    }
    for (uint32_t i = 0; i < 3; ++i) {
        if (!column_open(arguments[i], false, &strings[i], error)) goto cleanup;
    }
    if (!options_open(argc == 4 ? arguments[3] : NULL, keys, kinds, ARRANGEMENT_SCRIPT_OPTIONS, &option, error)) {
        goto cleanup;
    }
    scripts = calloc(rows ? rows : 1, sizeof(*scripts));
    if (!scripts) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_haplotype_arrangements_load_sql: out of memory");
        goto cleanup;
    }
    for (idx_t row = 0; row < rows; ++row) {
        char *values[3] = {0};
        arrangement_script_options settings = {0};
        duckvep_cell_t cell;
        char message[256];
        const char *failure = NULL;
        int built = 0;
        if (!options_check(&option, row, error)) goto cleanup;
        for (uint32_t i = 0; i < 3; ++i) {
            if (column_valid(&strings[i], row)) {
                duckdb_v2_str text = string_at(&strings[i].view, row);
                values[i] = duckvep_core_string_copy(text.ptr, text.len);
            }
            if (!values[i] || !*values[i]) {
                failure = "duckvep_haplotype_arrangements_load_sql: the calls query, model and job must be non-empty strings";
            }
        }
        for (unsigned i = 0; !failure && i < DUCKVEP_ARRANGEMENT_LIMIT_COUNT; ++i) {
            if (option_cell(&option, i, row, &cell)) {
                bool negative = cell.kind == DUCKVEP_CELL_TINYINT || cell.kind == DUCKVEP_CELL_SMALLINT ||
                                cell.kind == DUCKVEP_CELL_INTEGER || cell.kind == DUCKVEP_CELL_BIGINT
                                    ? cell.i < 0 : false;
                if (!cell.valid || negative ||
                    ((cell.kind == DUCKVEP_CELL_TINYINT || cell.kind == DUCKVEP_CELL_SMALLINT ||
                      cell.kind == DUCKVEP_CELL_INTEGER || cell.kind == DUCKVEP_CELL_BIGINT) && cell.i == 0) ||
                    (!(cell.kind == DUCKVEP_CELL_TINYINT || cell.kind == DUCKVEP_CELL_SMALLINT ||
                       cell.kind == DUCKVEP_CELL_INTEGER || cell.kind == DUCKVEP_CELL_BIGINT) && cell.u == 0)) {
                    (void)snprintf(message, sizeof message, "duckvep_haplotype_arrangements_load_sql: invalid %s",
                                   duckvep_arrangement_limit_names[i]);
                    failure = message;
                } else {
                    settings.has_limit[i] = true;
                    settings.limit[i] = cell.kind == DUCKVEP_CELL_TINYINT || cell.kind == DUCKVEP_CELL_SMALLINT ||
                                        cell.kind == DUCKVEP_CELL_INTEGER || cell.kind == DUCKVEP_CELL_BIGINT
                                            ? (uint64_t)cell.i : cell.u;
                }
            }
        }
        if (!failure) {
            built = arrangement_script(values[0], values[1], values[2], &settings, &scripts[row].statements[0],
                                       message, sizeof message);
            scripts[row].count = built ? 0u : 1u;
            if (built == 1) failure = "duckvep_haplotype_arrangements_load_sql: out of memory";
            else if (built == 2) failure = message;
        }
        if (failure) report(*error, failure);
        for (size_t i = 0; i < 3; ++i) duckvep_budget_free(values[i]);
        if (failure) goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &list_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &list_validity, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_child(output, 0, &child, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(child, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(child, &child_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(child, &arena, &detail));
    {
        uint64_t *child_validity = NULL;
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(child, &child_validity, &detail));
        for (idx_t row = 0; row < rows; ++row) mark_valid(child_validity, row);
    }
    for (idx_t row = 0; row < rows; ++row) {
        ((duckdb_v2_list_entry *)list_data)[row] = (duckdb_v2_list_entry){row, 1};
        mark_valid(list_validity, row);
        DUCKDB_CALL(write_string(arena, &((duckdb_v2_bytes *)child_data)[row], scripts[row].statements[0].data,
                                 scripts[row].statements[0].length, &detail));
    }
cleanup:
    if (scripts) {
        for (idx_t row = 0; row < rows; ++row) duckvep_sql_free(&scripts[row].statements[0]);
        free(scripts);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_haplotype_drop(job): releases a staged job that will not be scanned
 * ------------------------------------------------------------------------- */

static void drop_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                      duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    duckdb_v2_vector_view views[1];
    duckdb_v2_vector_handle output = NULL;
    void *output_data = NULL;
    uint64_t *output_validity = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 1, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &output_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &output_validity, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        char *name = NULL;
        bool found = false;
        stage *s;
        if (row_is_valid(&views[0], row)) {
            duckdb_v2_str text = string_at(&views[0], row);
            name = duckvep_core_string_copy(text.ptr, text.len);
        }
        if (!name || !*name) {
            duckvep_budget_free(name);
            INPUT_ERROR("duckvep_haplotype_drop: job must be a non-empty string");
        }
        while ((s = host_v2_stage_take(state, name, HAP_RELATION_CALLS)) != NULL ||
               (s = host_v2_stage_take(state, name, HAP_RELATION_PLAN)) != NULL ||
               (s = host_v2_stage_take(state, name, HAP_RELATION_ARRANGEMENTS)) != NULL) {
            host_v2_stage_destroy(s);
            found = true;
        }
        duckvep_budget_free(name);
        ((bool *)output_data)[row] = found;
        mark_valid(output_validity, row);
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_coding_transcripts(model, region, position, reference, alternate) -> UINTEGER[]
 * ------------------------------------------------------------------------- */

static bool hugeint_to_i64(const duckdb_v2_hugeint_t *value, int64_t *out) {
    if (value->upper == 0 && value->lower <= (uint64_t)INT64_MAX) {
        *out = (int64_t)value->lower;
        return true;
    }
    if (value->upper == -1 && value->lower >= (uint64_t)INT64_MAX + 1u) {
        *out = (int64_t)value->lower;
        return true;
    }
    return false;
}

static void coding_transcripts_exec(duckdb_v2_scalar_function_exec_info_handle info,
                                    duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    duckdb_v2_vector_handle arguments[5] = {0};
    column inputs[5];
    duckdb_v2_vector_handle output = NULL, child = NULL;
    duckdb_v2_list_entry *entries = NULL;
    uint64_t *validity = NULL;
    uint32_t *child_data = NULL;
    idx_t rows = 0;
    duckvep_model_entry_t *entry = NULL;
    duckvep_discovery_t scratch;
    duckvep_u32_list_t out = {NULL, 0u, 0u};
    char *current_name = NULL;
    const char *failure = NULL;
    duckvep_discovery_init(&scratch);
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    for (uint32_t i = 0; i < 5; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
        if (!column_open(arguments[i], false, &inputs[i], error)) {
            goto cleanup;
        }
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, (void **)&entries, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &validity, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        mark_valid(validity, row);
    }
    for (idx_t row = 0; row < rows && !failure; ++row) {
        bool missing = false;
        duckdb_v2_str name_text, ref, alt;
        int64_t start = 0, region = -1;
        size_t added = 0;
        duckvep_discovery_status_t status;
        for (unsigned k = 0; k < 5; ++k) {
            missing = missing || !column_valid(&inputs[k], row);
        }
        entries[row].offset = out.count;
        entries[row].length = 0;
        if (missing) {
            mark_null(validity, row);
            continue;
        }
        /* The model name is almost always one constant: compare in place and copy it only when it changes. */
        name_text = string_at(&inputs[0].view, row);
        if (!current_name || strlen(current_name) != name_text.len || memcmp(current_name, name_text.ptr, name_text.len)) {
            char *name = duckvep_core_string_copy(name_text.ptr, name_text.len);
            if (!name) {
                failure = "duckvep_coding_transcripts: out of memory copying the model name";
                break;
            }
            if (entry) {
                duckvep_registry_unpin(state->registry, entry);
            }
            entry = duckvep_registry_pin(state->registry, name);
            duckvep_budget_free(current_name);
            current_name = name;
            duckvep_discovery_reset_model(&scratch);
            if (!entry) {
                failure = "duckvep_coding_transcripts: unknown model name";
                break;
            }
            if (entry->model.lifted) {
                failure = "duckvep_coding_transcripts: transcript discovery is not supported for models with wrapped "
                          "circular objects";
                break;
            }
        }
        if (!hugeint_to_i64(&((const duckdb_v2_hugeint_t *)inputs[2].view.data)[physical_row(&inputs[2].view, row)],
                            &start)) {
            start = 0;
        }
        if (!hugeint_to_i64(&((const duckdb_v2_hugeint_t *)inputs[1].view.data)[physical_row(&inputs[1].view, row)],
                            &region)) {
            region = -1;
        }
        ref = string_at(&inputs[3].view, row);
        alt = string_at(&inputs[4].view, row);
        status = duckvep_discover_coding(&entry->model, &scratch, region, start, (const uint8_t *)ref.ptr, ref.len,
                                         (const uint8_t *)alt.ptr, alt.len, &out, &added);
        if (status == DUCKVEP_DISCOVERY_BAD_SPAN) {
            failure = "duckvep_coding_transcripts: position must be from 1 through 2147483647 and alleles at most "
                      "65535 bases";
        } else if (status == DUCKVEP_DISCOVERY_NOMEM) {
            failure = "duckvep_coding_transcripts: out of memory collecting transcripts";
        } else {
            entries[row].length = added;
        }
    }
    if (failure) {
        report(*error, failure);
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_vector_get_child(output, 0, &child, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(child, out.count, &detail));
    if (out.count) {
        uint64_t *child_validity = NULL;
        DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(child, (void **)&child_data, &detail));
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(child, &child_validity, &detail));
        memcpy(child_data, out.items, out.count * sizeof *out.items);
        for (size_t i = 0; i < out.count; ++i) {
            mark_valid(child_validity, i);
        }
    }
cleanup:
    if (entry) {
        duckvep_registry_unpin(state->registry, entry);
    }
    duckvep_u32_list_release(&out);
    duckvep_budget_free(current_name);
    duckvep_discovery_release(&scratch);
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * Registration
 * ------------------------------------------------------------------------- */

static bool register_scalar_with_state(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                                       model_state *state, const char *name, const char *const *parameter_types,
                                       const char *const *parameter_names, idx_t parameter_count,
                                       const char *return_type,
                                       duckdb_v2_scalar_function_exec_callback_fn callback,
                                       duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_scalar_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_str text = string_view(name);
    bool success = false, retained = false;
    DUCKDB_CALL(duckdb_v2_scalar_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_name(function, &text, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_get_signature(function, &signature, &detail));
    for (idx_t i = 0; i < parameter_count; ++i) {
        DUCKDB_CALL(make_type(context, parameter_types[i], &type, &detail));
        DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view(parameter_names[i]), type,
                                                               NULL, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    DUCKDB_CALL(make_type(context, return_type, &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_set_return_type(signature, type, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(
        function, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY_VOLATILE, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(
        function, DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING, DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING_SPECIAL,
        &detail));
    if (state) {
        duckdb_v2_opaque data = {state, host_v2_state_release, NULL};
        host_v2_state_retain(state);
        retained = true;
        DUCKDB_CALL(duckdb_v2_scalar_function_set_user_data(function, &data, &detail));
        retained = false;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_set_exec_callback(function, callback, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_register(function, &detail));
    success = true;
cleanup:
    if (retained) {
        host_v2_state_release(state);
    }
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_scalar_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

static bool register_job_table(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                               model_state *state, const char *name,
                               duckdb_v2_table_function_bind_callback_fn bind,
                               duckdb_v2_table_function_init_global_callback_fn init,
                               duckdb_v2_table_function_exec_callback_fn exec, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_table_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_str text = string_view(name);
    bool success = false, retained = false;
    DUCKDB_CALL(duckdb_v2_table_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_name(function, &text, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_get_signature(function, &signature, &detail));
    DUCKDB_CALL(make_type(context, "VARCHAR", &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view("job"), type, NULL, &detail));
    {
        duckdb_v2_opaque data = {state, host_v2_state_release, NULL};
        host_v2_state_retain(state);
        retained = true;
        DUCKDB_CALL(duckdb_v2_table_function_set_user_data(function, &data, &detail));
        retained = false;
    }
    DUCKDB_CALL(duckdb_v2_table_function_set_bind_callback(function, bind, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_init_global_callback(function, init, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_exec_callback(function, exec, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_register(function, &detail));
    success = true;
cleanup:
    if (retained) {
        host_v2_state_release(state);
    }
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_table_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

bool host_v2_register_coding_calls(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                                   model_state *state, duckdb_v2_error_info_handle *error);

bool host_v2_register_haplotypes(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                                 model_state *state, duckdb_v2_error_info_handle *error) {
    static const char *const script_types[] = {"VARCHAR", "VARCHAR", "VARCHAR", "ANY"};
    static const char *const script_names[] = {"calls_query", "model", "job", "options"};
    static const char *const one_type[] = {"VARCHAR"};
    static const char *const one_name[] = {"job"};
    static const char *const discovery_types[] = {"VARCHAR", "HUGEINT", "HUGEINT", "VARCHAR", "VARCHAR"};
    static const char *const discovery_names[] = {"model", "region", "position", "reference", "alternate"};
    return register_scalar_with_state(extension, context, NULL, "duckvep_haplotype_load_sql", script_types,
                                      script_names, 3, "VARCHAR[]", load_sql_exec, error) &&
           register_scalar_with_state(extension, context, NULL, "duckvep_haplotype_load_sql", script_types,
                                      script_names, 4, "VARCHAR[]", load_sql_exec, error) &&
           register_scalar_with_state(extension, context, NULL, "duckvep_haplotype_arrangements_load_sql",
                                      script_types, script_names, 3, "VARCHAR[]", arrangement_load_sql_exec, error) &&
           register_scalar_with_state(extension, context, NULL, "duckvep_haplotype_arrangements_load_sql",
                                      script_types, script_names, 4, "VARCHAR[]", arrangement_load_sql_exec, error) &&
           register_scalar_with_state(extension, context, state, "duckvep_haplotype_drop", one_type, one_name, 1,
                                      "BOOLEAN", drop_exec, error) &&
           register_scalar_with_state(extension, context, state, "duckvep_coding_transcripts", discovery_types,
                                      discovery_names, 5, "UINTEGER[]", coding_transcripts_exec, error) &&
           register_job_table(extension, context, state, "duckvep_haplotype_scan", scan_bind_exec, scan_init_exec,
                              scan_exec, error) &&
           register_job_table(extension, context, state, "duckvep_haplotype_arrangements", arrangement_bind_exec,
                              arrangement_init_exec, arrangement_exec, error) &&
           register_job_table(extension, context, state, "_duckvep_haplotype_plan", plan_bind_exec, plan_init_exec,
                              plan_exec, error) &&
           host_v2_register_coding_calls(extension, context, state, error);
}
