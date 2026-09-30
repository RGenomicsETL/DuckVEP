/* The resource-control functions of the v2 host: duckvep_native_budget (table), duckvep_native_budget_set,
 * duckvep_native_budget_reset_high_water and duckvep_worker_limits_set. The rows, the checks and the messages are
 * src/core/duckvep_core_budget.c, shared with v1; the budget is process-wide state of the extension binary, so the
 * annotation workers, the model publish, the haplotype scan and the coding_calls reader that run on v2 are charged to
 * the same budget these functions report. */
#include "host_v2_columns.h"
#include "host_v2_stage.h"

#include "core/duckvep_core_budget.h"

bool host_v2_register_budget(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error);

static const char *const table_names[6] = {"owner", "current_bytes", "high_water_bytes", "limit_bytes", "charges",
                                           "refusals"};

typedef struct {
    bool done;
} budget_scan;

static void scan_destroy(void *pointer) {
    free(pointer);
}

static void budget_bind(duckdb_v2_table_function_bind_info_handle info, duckdb_v2_context_handle context,
                        duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    for (unsigned i = 0; i < 6; ++i) {
        DUCKDB_CALL(make_type(context, i == 0 ? "VARCHAR" : "UBIGINT", &type, &detail));
        DUCKDB_CALL(duckdb_v2_table_function_bind_add_result_column(info, string_view(table_names[i]), type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
cleanup:
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void budget_init(duckdb_v2_table_function_init_global_info_handle info, duckdb_v2_context_handle context,
                        duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    budget_scan *scan = calloc(1, sizeof(*scan));
    bool owned = false;
    if (!scan) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY,
                  "duckvep_native_budget: capacity error: cannot allocate scan state");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_max_threads(info, 1, &detail));
    {
        duckdb_v2_opaque data = {scan, scan_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_init_global_set_global_state(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        scan_destroy(scan);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void budget_exec(duckdb_v2_table_function_exec_info_handle info, duckdb_v2_context_handle context,
                        duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    budget_scan *scan = NULL;
    duckdb_v2_data_chunk_handle chunk = NULL;
    duckdb_v2_vector_handle columns[6] = {0};
    void *data[6] = {0};
    uint64_t *validity[6] = {0};
    duckdb_v2_arena_handle arena = NULL;
    duckvep_budget_stats_t stats;
    idx_t count = 0;
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_global_state(info, (void **)&scan, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_output_chunk(info, &chunk, &detail));
    if (scan->done) {
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, 0, &columns[0], &detail));
        DUCKDB_CALL(duckdb_v2_vector_set_size(columns[0], 0, &detail));
        goto cleanup;
    }
    duckvep_budget_stats(&stats);
    for (unsigned i = 0; i < 6; ++i) {
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, i, &columns[i], &detail));
        DUCKDB_CALL(duckdb_v2_vector_set_size(columns[i], DUCKVEP_CORE_BUDGET_ROWS, &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(columns[i], &data[i], &detail));
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(columns[i], &validity[i], &detail));
        for (idx_t row = 0; row < DUCKVEP_CORE_BUDGET_ROWS; ++row) {
            mark_valid(validity[i], row);
        }
    }
    DUCKDB_CALL(duckdb_v2_vector_get_arena(columns[0], &arena, &detail));
    for (; count < DUCKVEP_CORE_BUDGET_ROWS; ++count) {
        duckvep_core_budget_row_t row;
        duckvep_core_budget_row((unsigned)count, &stats, &row);
        DUCKDB_CALL(write_string(arena, &((duckdb_v2_bytes *)data[0])[count], row.owner, strlen(row.owner), &detail));
        ((uint64_t *)data[1])[count] = row.current;
        ((uint64_t *)data[2])[count] = row.high_water;
        ((uint64_t *)data[3])[count] = row.limit;
        ((uint64_t *)data[4])[count] = row.charges;
        ((uint64_t *)data[5])[count] = row.refusals;
    }
    scan->done = true;
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* Reads up to four BIGINT arguments; `present` tells which are non-NULL for the row. */
typedef enum { KIND_SET, KIND_RESET, KIND_LIMITS } budget_kind;

static void scalar_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_error_info_handle *error,
                        budget_kind kind, idx_t arguments) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[4];
    duckdb_v2_vector_handle output = NULL;
    void *data = NULL;
    uint64_t *validity = NULL;
    idx_t rows = 0;
    char message[256];
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (arguments && !load_views(info, arguments, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &validity, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        int present[4] = {0, 0, 0, 0};
        int64_t value[4] = {0, 0, 0, 0};
        for (idx_t k = 0; k < arguments; ++k) {
            present[k] = row_is_valid(&views[k], row);
            if (present[k]) {
                value[k] = (int64_t)integer_at(&views[k], row);
            }
        }
        if (kind == KIND_SET) {
            if (!duckvep_core_budget_set(present[0], value[0], message, sizeof message)) {
                INPUT_ERROR(message);
            }
            ((int64_t *)data)[row] = value[0];
        } else if (kind == KIND_LIMITS) {
            if (!duckvep_core_worker_limits(present, value, message, sizeof message)) {
                INPUT_ERROR(message);
            }
            ((bool *)data)[row] = true;
        } else {
            duckvep_budget_reset_high_water();
            ((bool *)data)[row] = true;
        }
        mark_valid(validity, row);
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void set_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                     duckdb_v2_error_info_handle *error) {
    (void)context;
    scalar_exec(info, error, KIND_SET, 1);
}

static void limits_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                        duckdb_v2_error_info_handle *error) {
    (void)context;
    scalar_exec(info, error, KIND_LIMITS, 4);
}

static void reset_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                       duckdb_v2_error_info_handle *error) {
    (void)context;
    scalar_exec(info, error, KIND_RESET, 0);
}

static bool register_bigint_scalar(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                                   const char *name, idx_t parameters, const char *return_type,
                                   duckdb_v2_scalar_function_exec_callback_fn callback,
                                   duckdb_v2_error_info_handle *error) {
    static const char *const types[4] = {"BIGINT", "BIGINT", "BIGINT", "BIGINT"};
    static const char *const names[4] = {"a", "b", "c", "d"};
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_scalar_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_str text = string_view(name);
    bool success = false;
    DUCKDB_CALL(duckdb_v2_scalar_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_name(function, &text, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_get_signature(function, &signature, &detail));
    for (idx_t i = 0; i < parameters; ++i) {
        DUCKDB_CALL(make_type(context, types[i], &type, &detail));
        DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view(names[i]), type, NULL, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    DUCKDB_CALL(make_type(context, return_type, &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_set_return_type(signature, type, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(
        function, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY_VOLATILE, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_exec_callback(function, callback, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_register(function, &detail));
    success = true;
cleanup:
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_scalar_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

bool host_v2_register_budget(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_table_function_handle function = NULL;
    duckdb_v2_str name = string_view("duckvep_native_budget");
    bool success = false;
    DUCKDB_CALL(duckdb_v2_table_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_name(function, &name, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_bind_callback(function, budget_bind, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_init_global_callback(function, budget_init, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_exec_callback(function, budget_exec, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_register(function, &detail));
    success = register_bigint_scalar(extension, context, "duckvep_native_budget_set", 1, "BIGINT", set_exec, error) &&
              register_bigint_scalar(extension, context, "duckvep_worker_limits_set", 4, "BOOLEAN", limits_exec,
                                     error) &&
              register_bigint_scalar(extension, context, "duckvep_native_budget_reset_high_water", 0, "BOOLEAN",
                                     reset_exec, error);
cleanup:
    (void)duckdb_v2_table_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}
