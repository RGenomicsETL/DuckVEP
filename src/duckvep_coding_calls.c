/* duckvep_coding_calls(model, path): the v1 (stable C API) table function over the host-neutral fused VCF/BCF reader of
 * src/core/duckvep_core_coding_calls.c. The reader yields one calls row at a time; this file binds the function, pins the
 * model and writes the rows into the output vectors. */
#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
DUCKDB_EXTENSION_EXTERN
#include "duckvep_list.h"

#include "core/duckvep_core_coding_calls.h"
#include "duckvep_model.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define CC "duckvep_coding_calls: "

typedef struct {
    duckvep_registry_t *registry;
    duckvep_model_entry_t *entry;
    duckvep_cc_plan_t *plan;
} cc_bind_t;

typedef struct {
    duckvep_cc_reader_t *reader;
} cc_state_t;

static void cc_bind_destroy(void *pointer) {
    cc_bind_t *b = pointer;
    if (!b) return;
    duckvep_cc_plan_destroy(b->plan);
    duckvep_registry_unpin(b->registry, b->entry);
    duckvep_registry_release(b->registry);
    duckvep_budget_free(b);
}

static void cc_add_column(duckdb_bind_info info, const char *name, duckdb_type type, int list) {
    duckdb_logical_type element = duckdb_create_logical_type(type), result = element;
    if (list) result = duckdb_create_list_type(element);
    duckdb_bind_add_result_column(info, name, result);
    if (list) duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&element);
}

static void cc_bind(duckdb_bind_info info) {
    cc_bind_t *b = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1u, sizeof(*b));
    if (!b) { duckdb_bind_set_error(info, CC "bind allocation failed"); return; }
    b->registry = duckdb_bind_get_extra_info(info);
    duckvep_registry_retain(b->registry);
    duckdb_value value = duckdb_bind_get_parameter(info, 0u);
    char *name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    duckdb_destroy_value(&value);
    value = duckdb_bind_get_parameter(info, 1u);
    char *path = NULL;
    if (value && !duckdb_is_null_value(value)) path = duckdb_get_varchar(value);
    duckdb_destroy_value(&value);
    if (name) b->entry = duckvep_registry_pin(b->registry, name);
    duckdb_free(name);
    if (!b->entry || !path || !path[0]) {
        duckdb_bind_set_error(info, CC "require a loaded model name and a nonempty VCF or BCF path");
        duckdb_free(path); cc_bind_destroy(b); return;
    }
    char error[DUCKVEP_SQL_ERROR_SIZE] = {0};
    b->plan = duckvep_cc_plan_create(&b->entry->model, path, error, sizeof error);
    duckdb_free(path);
    if (!b->plan) { duckdb_bind_set_error(info, error); cc_bind_destroy(b); return; }
    cc_add_column(info, "event_index", DUCKDB_TYPE_BIGINT, 0);
    cc_add_column(info, "seq_region", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "position", DUCKDB_TYPE_BIGINT, 0);
    cc_add_column(info, "reference", DUCKDB_TYPE_VARCHAR, 0);
    cc_add_column(info, "alternate", DUCKDB_TYPE_VARCHAR, 0);
    cc_add_column(info, "alt_index", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "transcript_index", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "sample_index", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "alleles", DUCKDB_TYPE_INTEGER, 1);
    cc_add_column(info, "phase_before", DUCKDB_TYPE_BOOLEAN, 1);
    cc_add_column(info, "phase_set", DUCKDB_TYPE_BIGINT, 0);
    duckdb_bind_set_bind_data(info, b, cc_bind_destroy);
}

static void cc_state_destroy(void *pointer) {
    cc_state_t *s = pointer;
    if (!s) return;
    duckvep_cc_close(s->reader);
    duckvep_budget_free(s);
}

static void cc_init(duckdb_init_info info) {
    const cc_bind_t *bind = duckdb_init_get_bind_data(info);
    cc_state_t *s = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1u, sizeof(*s));
    char error[DUCKVEP_SQL_ERROR_SIZE] = CC "allocation failed or native budget exceeded";
    duckdb_init_set_max_threads(info, 1u);
    duckvep_budget_clear_failure();
    if (!s) goto failed;
    s->reader = duckvep_cc_open(bind->plan, error, sizeof error);
    if (!s->reader) goto failed;
    duckdb_init_set_init_data(info, s, cc_state_destroy);
    return;
failed:
    duckdb_init_set_error(info, error);
    cc_state_destroy(s);
}

static void cc_scan(duckdb_function_info info, duckdb_data_chunk output) {
    cc_state_t *s = duckdb_function_get_init_data(info);
    char error[DUCKVEP_SQL_ERROR_SIZE] = {0};
    idx_t capacity = duckdb_vector_size(), rows = 0u;
    duckdb_vector v[11];
    for (unsigned i = 0u; i < 11u; i++) {
        v[i] = duckdb_data_chunk_get_vector(output, i);
        duckdb_vector_ensure_validity_writable(v[i]);
    }
    if (duckdb_list_vector_set_size(v[8], 0u) != DuckDBSuccess || duckdb_list_vector_set_size(v[9], 0u) != DuckDBSuccess) {
        duckdb_function_set_error(info, CC "cannot reset output list"); return;
    }
    while (rows < capacity) {
        duckvep_cc_row_t row;
        int r = duckvep_cc_next(s->reader, &row, error, sizeof error);
        if (r < 0) { duckdb_function_set_error(info, error); return; }
        if (r == 0) break;
        ((int64_t *)duckdb_vector_get_data(v[0]))[rows] = row.event_index;
        ((int32_t *)duckdb_vector_get_data(v[1]))[rows] = row.region;
        ((int64_t *)duckdb_vector_get_data(v[2]))[rows] = row.position;
        duckdb_vector_assign_string_element_len(v[3], rows, row.ref, row.ref_length);
        duckdb_vector_assign_string_element_len(v[4], rows, row.alt, row.alt_length);
        ((int32_t *)duckdb_vector_get_data(v[5]))[rows] = row.alt_index;
        ((int32_t *)duckdb_vector_get_data(v[6]))[rows] = row.transcript;
        ((int32_t *)duckdb_vector_get_data(v[7]))[rows] = row.sample;
        duckdb_list_entry *alleles_entry = (duckdb_list_entry *)duckdb_vector_get_data(v[8]) + rows;
        duckdb_list_entry *phase_entry = (duckdb_list_entry *)duckdb_vector_get_data(v[9]) + rows;
        if (!duckvep_list_extend(v[8], row.lanes, alleles_entry) || !duckvep_list_extend(v[9], row.lanes, phase_entry)) {
            duckdb_function_set_error(info, CC "output list allocation failed"); return;
        }
        duckdb_vector allele_child = duckdb_list_vector_get_child(v[8]), phase_child = duckdb_list_vector_get_child(v[9]);
        duckdb_vector_ensure_validity_writable(allele_child);
        uint64_t *allele_validity = duckdb_vector_get_validity(allele_child);
        int32_t *allele_data = duckdb_vector_get_data(allele_child);
        bool *phase_data = duckdb_vector_get_data(phase_child);
        for (uint32_t k = 0u; k < row.lanes; k++) {
            int32_t a = row.alleles[k];
            if (a < 0) { allele_data[alleles_entry->offset + k] = 0; duckdb_validity_set_row_invalid(allele_validity, alleles_entry->offset + k); }
            else { allele_data[alleles_entry->offset + k] = a; duckdb_validity_set_row_valid(allele_validity, alleles_entry->offset + k); }
            phase_data[phase_entry->offset + k] = row.phase[k] != 0u;
        }
        if (row.phase_set_present) {
            ((int64_t *)duckdb_vector_get_data(v[10]))[rows] = row.phase_set;
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(v[10]), rows);
        } else duckdb_validity_set_row_invalid(duckdb_vector_get_validity(v[10]), rows);
        rows++;
    }
    duckdb_data_chunk_set_size(output, rows);
}

void duckvep_register_coding_calls(duckdb_connection connection, duckvep_registry_t *registry) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type string = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_set_name(function, "duckvep_coding_calls");
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_parameter(function, string);
    duckvep_registry_retain(registry);
    duckdb_table_function_set_extra_info(function, registry, duckvep_registry_release);
    duckdb_table_function_set_bind(function, cc_bind);
    duckdb_table_function_set_init(function, cc_init);
    duckdb_table_function_set_function(function, cc_scan);
    (void)duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&string);
}
