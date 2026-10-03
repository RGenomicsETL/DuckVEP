/* duckvep_coding_calls(model, path) on the v2 host: a table function over the host-neutral fused VCF/BCF reader of
 * src/core/duckvep_core_coding_calls.c. It needs no caller relation, so it binds and scans directly. Each exec call
 * emits a full chunk of rows; the list children grow as rows are appended (duckvep_host.h). */
#include "duckvep_host.h"

#include "host_v2_columns.h"
#include "host_v2_stage.h"

#include "core/duckvep_core_coding_calls.h"
#include "kernel/src/duckvep_budget.h"

#include <htslib/hts.h>

#define CC "duckvep_coding_calls: "

typedef struct {
    model_state *state;
    duckvep_model_entry_t *entry;
    duckvep_cc_plan_t *plan;
} cc_bind;

typedef struct {
    duckvep_cc_reader_t *reader;
} cc_scan;

static void cc_bind_destroy(void *pointer) {
    cc_bind *b = pointer;
    if (!b) {
        return;
    }
    duckvep_cc_plan_destroy(b->plan);
    if (b->entry) {
        duckvep_registry_unpin(b->state->registry, b->entry);
    }
    host_v2_state_release(b->state);
    free(b);
}

static char *argument_text(duckdb_v2_table_function_bind_info_handle info, idx_t index,
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
        copy = malloc(text.len + 1);
        if (copy) {
            memcpy(copy, text.ptr, text.len);
            copy[text.len] = '\0';
        }
    }
cleanup:
    (void)duckdb_v2_value_destroy(&value);
    (void)duckdb_v2_error_info_destroy(&detail);
    return copy;
}

static void cc_bind_exec(duckdb_v2_table_function_bind_info_handle info, duckdb_v2_context_handle context,
                         duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    model_state *state = NULL;
    cc_bind *b = calloc(1, sizeof(*b));
    char *name = NULL, *path = NULL;
    char message[DUCKVEP_SQL_ERROR_SIZE];
    bool owned = false;
    if (!b) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, CC "out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_bind_get_user_data(info, (void **)&state, &detail));
    b->state = state;
    host_v2_state_retain(state);
    name = argument_text(info, 0, error);
    path = argument_text(info, 1, error);
    if (name) {
        b->entry = duckvep_registry_pin(state->registry, name);
    }
    if (!b->entry || !path || !path[0]) {
        INPUT_ERROR(CC "require a loaded model name and a nonempty VCF or BCF path");
    }
    b->plan = duckvep_cc_plan_create(&b->entry->model, path, message, sizeof message);
    if (!b->plan) {
        INPUT_ERROR(message);
    }
    for (unsigned i = 0; i < DUCKVEP_CC_COLUMNS; ++i) {
        DUCKDB_CALL(make_type(context, duckvep_cc_column_type(i), &type, &detail));
        DUCKDB_CALL(duckdb_v2_table_function_bind_add_result_column(
            info, string_view(duckvep_cc_column_name(i)), type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    {
        duckdb_v2_opaque data = {b, cc_bind_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_bind_set_bind_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        cc_bind_destroy(b);
    }
    free(name);
    free(path);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void cc_scan_destroy(void *pointer) {
    cc_scan *s = pointer;
    if (s) {
        duckvep_cc_close(s->reader);
        free(s);
    }
}

static void cc_init_exec(duckdb_v2_table_function_init_global_info_handle info, duckdb_v2_context_handle context,
                         duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    cc_bind *b = NULL;
    cc_scan *s = calloc(1, sizeof(*s));
    char message[DUCKVEP_SQL_ERROR_SIZE];
    bool owned = false;
    if (!s) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, CC "out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_get_bind_data(info, (void **)&b, &detail));
    s->reader = duckvep_cc_open(b->plan, message, sizeof message);
    if (!s->reader) {
        report(*error, message);
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_max_threads(info, 1, &detail));
    {
        duckdb_v2_opaque data = {s, cc_scan_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_table_function_init_global_set_global_state(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        cc_scan_destroy(s);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void cc_exec(duckdb_v2_table_function_exec_info_handle info, duckdb_v2_context_handle context,
                    duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    cc_scan *s = NULL;
    duckdb_v2_data_chunk_handle chunk = NULL;
    v2_call *call = calloc(1, sizeof(*call));
    duckvep_cc_row_t row;
    char message[DUCKVEP_SQL_ERROR_SIZE];
    const size_t capacity = duckvep_h_vector_size();
    size_t rows = 0;
    v2_vec *v;
    if (!call) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, CC "out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_global_state(info, (void **)&s, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_output_chunk(info, &chunk, &detail));
    call->error = error;
    call->rows = capacity;
    for (unsigned i = 0; i < DUCKVEP_CC_COLUMNS; ++i) {
        duckdb_v2_vector_handle vector = NULL;
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, i, &vector, &detail));
        if (!v2_open_writable(call, &call->inputs[i], vector, capacity, true)) {
            goto cleanup;
        }
    }
    v = call->inputs;
    /* A full chunk of rows, as on v1: the two list children grow as rows are appended (duckvep_host.h). */
    while (rows < capacity) {
        duckdb_v2_list_entry alleles_entry, phase_entry;
        int fetched;
        message[0] = '\0';
        fetched = duckvep_cc_next(s->reader, &row, message, sizeof message);
        if (fetched < 0) {
            report(*error, message);
            goto cleanup;
        }
        if (fetched == 0) {
            break;
        }
        ((int64_t *)v[0].data)[rows] = row.event_index;
        ((int32_t *)v[1].data)[rows] = row.region;
        ((int64_t *)v[2].data)[rows] = row.position;
        duckvep_h_assign_string(&v[3], rows, row.ref, row.ref_length);
        duckvep_h_assign_string(&v[4], rows, row.alt, row.alt_length);
        ((int32_t *)v[5].data)[rows] = row.alt_index;
        ((int32_t *)v[6].data)[rows] = row.transcript;
        ((int32_t *)v[7].data)[rows] = row.sample;
        if (!duckvep_h_list_extend(&v[8], row.lanes, &alleles_entry) ||
            !duckvep_h_list_extend(&v[9], row.lanes, &phase_entry)) {
            INPUT_ERROR(CC "output list allocation failed");
        }
        ((duckdb_v2_list_entry *)v[8].data)[rows] = alleles_entry;
        ((duckdb_v2_list_entry *)v[9].data)[rows] = phase_entry;
        if (row.lanes) {
            v2_vec *allele_child = duckvep_h_list_values(&v[8]), *phase_child = duckvep_h_list_values(&v[9]);
            if (call->failed) {
                goto cleanup;
            }
            for (uint32_t k = 0; k < row.lanes; ++k) {
                int32_t a = row.alleles[k];
                ((int32_t *)allele_child->data)[alleles_entry.offset + k] = a < 0 ? 0 : a;
                if (a < 0) {
                    mark_null(allele_child->validity, alleles_entry.offset + k);
                }
                ((bool *)phase_child->data)[phase_entry.offset + k] = row.phase[k] != 0u;
            }
        }
        if (row.phase_set_present) {
            ((int64_t *)v[10].data)[rows] = row.phase_set;
        } else {
            mark_null(v[10].validity, rows);
        }
        ++rows;
    }
    v2_finish(call);
    if (!call->failed && rows != capacity) {
        /* The rows written; an empty batch ends the scan. */
        for (unsigned i = 0; i < DUCKVEP_CC_COLUMNS; ++i) {
            DUCKDB_CALL(duckdb_v2_vector_set_size(call->inputs[i].handle, rows, &detail));
        }
    }
cleanup:
    free(call);
    (void)duckdb_v2_error_info_destroy(&detail);
}

bool host_v2_register_coding_calls(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                                   model_state *state, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_table_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_str name = string_view("duckvep_coding_calls");
    bool success = false, retained = false;
    DUCKDB_CALL(duckdb_v2_table_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_name(function, &name, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_get_signature(function, &signature, &detail));
    DUCKDB_CALL(make_type(context, "VARCHAR", &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view("model"), type, NULL, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view("path"), type, NULL, &detail));
    {
        duckdb_v2_opaque data = {state, host_v2_state_release, NULL};
        host_v2_state_retain(state);
        retained = true;
        DUCKDB_CALL(duckdb_v2_table_function_set_user_data(function, &data, &detail));
        retained = false;
    }
    DUCKDB_CALL(duckdb_v2_table_function_set_bind_callback(function, cc_bind_exec, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_init_global_callback(function, cc_init_exec, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_exec_callback(function, cc_exec, &detail));
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
