#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_builder.h"
#include "core/duckvep_core_repeat.h"
DUCKDB_EXTENSION_EXTERN
#include "duckvep_v1_cells.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    duckdb_vector list, units, counts, records;
    duckdb_list_entry *entries;
    duckdb_string_t *strings;
    duckdb_logical_type count_type;
} repeat_axis;

static bool valid(duckdb_vector vector, idx_t row) {
    uint64_t *mask = duckdb_vector_get_validity(vector);
    return !mask || duckdb_validity_row_is_valid(mask, row);
}

static bool axis_init(duckdb_vector vector, repeat_axis *axis) {
    axis->list = vector;
    duckdb_logical_type type = duckdb_vector_get_column_type(vector);
    bool is_list = duckdb_get_type_id(type) == DUCKDB_TYPE_LIST;
    duckdb_destroy_logical_type(&type);
    if (!is_list) return !valid(vector, 0);
    axis->entries = duckdb_vector_get_data(vector);
    axis->records = duckdb_list_vector_get_child(vector);
    type = duckdb_vector_get_column_type(axis->records);
    duckdb_type record_id = duckdb_get_type_id(type);
    if (record_id == DUCKDB_TYPE_SQLNULL) {
        duckdb_destroy_logical_type(&type);
        return true;
    }
    bool is_struct = record_id == DUCKDB_TYPE_STRUCT && duckdb_struct_type_child_count(type) == 2;
    idx_t unit_index = 0;
    if (is_struct) {
        char *first = duckdb_struct_type_child_name(type, 0);
        char *second = duckdb_struct_type_child_name(type, 1);
        unit_index = strcmp(first, "unit") == 0 ? 0 : 1;
        is_struct = strcmp(first, unit_index ? "count" : "unit") == 0 &&
            strcmp(second, unit_index ? "unit" : "count") == 0;
        duckdb_free(first);
        duckdb_free(second);
    }
    duckdb_destroy_logical_type(&type);
    if (!is_struct) return false;
    axis->units = duckdb_struct_vector_get_child(axis->records, unit_index);
    axis->counts = duckdb_struct_vector_get_child(axis->records, 1 - unit_index);
    type = duckdb_vector_get_column_type(axis->units);
    bool is_text = duckdb_get_type_id(type) == DUCKDB_TYPE_VARCHAR;
    duckdb_destroy_logical_type(&type);
    if (!is_text) return false;
    axis->strings = duckdb_vector_get_data(axis->units);
    axis->count_type = duckdb_vector_get_column_type(axis->counts);
    return true;
}

static void axis_destroy(repeat_axis *axis) {
    if (axis->count_type) duckdb_destroy_logical_type(&axis->count_type);
}

static void set_error(duckdb_function_info info, const char *message) {
    duckvep_builder_set_error(info, message);
}

/* Reads one list row of an axis into neutral elements (grown on demand). */
static bool axis_row(repeat_axis *axis, idx_t row, duckdb_type count_id, uint8_t count_scale,
                     duckvep_repeat_element_t **elements, size_t *capacity,
                     duckvep_repeat_axis_t *out) {
    *out = (duckvep_repeat_axis_t){0};
    if (!axis->entries || !valid(axis->list, row) || !axis->units) return true;
    duckdb_list_entry list = axis->entries[row];
    if (list.length > *capacity) {
        duckvep_repeat_element_t *grown = duckvep_budget_realloc(DUCKVEP_OWNER_CONTROL, *elements,
            list.length * sizeof(**elements));
        if (!grown) return false;
        *elements = grown;
        *capacity = list.length;
    }
    for (idx_t i = 0; i < list.length; i++) {
        idx_t at = list.offset + i;
        duckvep_repeat_element_t *element = &(*elements)[i];
        *element = (duckvep_repeat_element_t){0};
        element->present = valid(axis->records, at) && valid(axis->units, at) && valid(axis->counts, at);
        element->unit = duckdb_string_t_data(&axis->strings[at]);
        element->unit_length = duckdb_string_t_length(axis->strings[at]);
        duckvep_v1_fill_cell(axis->counts, count_id, count_scale, at, &element->count);
    }
    out->usable = true;
    out->elements = *elements;
    out->count = list.length;
    return true;
}

static void repeat_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector args[4] = {0};
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    repeat_axis axes[2] = {{0}};
    duckvep_repeat_element_t *buffers[2] = {NULL, NULL};
    size_t capacities[2] = {0, 0};
    duckdb_type count_ids[2] = {DUCKDB_TYPE_INVALID, DUCKDB_TYPE_INVALID};
    uint8_t count_scales[2] = {0, 0};
    if (!axis_init(args[0], &axes[0]) || !axis_init(args[1], &axes[1])) {
        set_error(info, duckvep_core_repeat_expected_lists);
        axis_destroy(&axes[0]); axis_destroy(&axes[1]);
        return;
    }
    for (idx_t axis = 0; axis < 2; axis++) {
        if (!axes[axis].count_type) continue;
        count_ids[axis] = duckdb_get_type_id(axes[axis].count_type);
        if (count_ids[axis] == DUCKDB_TYPE_DECIMAL) {
            count_scales[axis] = duckdb_decimal_scale(axes[axis].count_type);
            count_ids[axis] = duckdb_decimal_internal_type(axes[axis].count_type);
        }
    }
    duckdb_vector fields[7];
    for (idx_t i = 0; i < 7; i++) {
        fields[i] = duckdb_struct_vector_get_child(output, i);
        duckdb_vector_ensure_validity_writable(fields[i]);
    }
    bool *exact = duckdb_vector_get_data(args[2]);
    const char *const keys[] = {"max_allele_bases"};
    const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_NUMERIC};
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        duckdb_vector cap_vector = NULL;
        if (argc == 4 && !duckvep_builder_option_vectors(info, args[3], row, keys, kinds, 1, &cap_vector)) break;
        if (!valid(args[2], row)) {
            set_error(info, duckvep_core_repeat_exact_required);
            break;
        }
        long double cap = 5000;
        duckvep_cell_t cap_cell;
        if (cap_vector) {
            duckdb_logical_type cap_type = duckdb_vector_get_column_type(cap_vector);
            duckdb_type cap_id = duckdb_get_type_id(cap_type);
            uint8_t cap_scale = 0;
            if (cap_id == DUCKDB_TYPE_DECIMAL) {
                cap_scale = duckdb_decimal_scale(cap_type);
                cap_id = duckdb_decimal_internal_type(cap_type);
            }
            duckdb_destroy_logical_type(&cap_type);
            duckvep_v1_fill_cell(cap_vector, cap_id, cap_scale, row, &cap_cell);
        }
        if (!duckvep_core_repeat_cap(cap_vector ? &cap_cell : NULL, &cap)) {
            set_error(info, duckvep_core_repeat_cap_invalid);
            break;
        }
        duckvep_repeat_axis_t rows[2];
        if (!axis_row(&axes[0], row, count_ids[0], count_scales[0], &buffers[0], &capacities[0], &rows[0]) ||
            !axis_row(&axes[1], row, count_ids[1], count_scales[1], &buffers[1], &capacities[1], &rows[1])) {
            set_error(info, "duckvep_repeat_alleles: allocation failed");
            break;
        }
        duckvep_repeat_plan_t plan;
        char error[240];
        duckvep_core_repeat_plan(rows, exact[row], cap, &plan, error, sizeof(error));
        if (plan.error) {
            set_error(info, plan.error);
            break;
        }
        for (idx_t i = 0; i < 7; i++)
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[i]), row);
        duckdb_vector_assign_string_element(fields[6], row, plan.status);
        if (strcmp(plan.status, "ok") != 0) {
            for (idx_t i = 0; i < 6; i++)
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[i]), row);
            continue;
        }
        bool failed = false;
        for (idx_t axis = 0; axis < 2; axis++) {
            size_t length = (size_t)plan.required[axis];
            char *text = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
            if (!text) {
                set_error(info, "duckvep_repeat_alleles: allocation failed");
                failed = true;
                break;
            }
            size_t written = duckvep_core_repeat_render(&rows[axis], text);
            duckdb_vector_assign_string_element_len(fields[axis], row, text, written);
            duckvep_budget_free(text);
        }
        if (failed) break;
        ((uint64_t *)duckdb_vector_get_data(fields[2]))[row] = (uint64_t)plan.required[0];
        ((uint64_t *)duckdb_vector_get_data(fields[3]))[row] = (uint64_t)plan.required[1];
        ((int64_t *)duckdb_vector_get_data(fields[4]))[row] = (int64_t)plan.required[1] - (int64_t)plan.required[0];
        duckdb_vector_assign_string_element(fields[5], row, duckvep_core_repeat_direction(&plan));
    }
    duckvep_budget_free(buffers[0]); duckvep_budget_free(buffers[1]);
    axis_destroy(&axes[0]); axis_destroy(&axes[1]);
}

bool duckvep_register_repeat_alleles(duckdb_connection connection) {
    duckdb_logical_type any = duckdb_create_logical_type(DUCKDB_TYPE_ANY);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type text = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type u64 = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type i64 = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type types[] = {text, text, u64, u64, i64, text, text};
    const char *names[] = {"reference", "alternate", "reference_length", "alternate_length",
        "length_change", "length_direction", "status"};
    duckdb_logical_type result = duckdb_create_struct_type(types, names, 7);
    duckdb_scalar_function_set overloads = duckdb_create_scalar_function_set("duckvep_repeat_alleles");
    for (idx_t arity = 3; arity <= 4; arity++) {
        duckdb_scalar_function function = duckdb_create_scalar_function();
        duckdb_scalar_function_set_name(function, "duckvep_repeat_alleles");
        duckdb_scalar_function_add_parameter(function, any);
        duckdb_scalar_function_add_parameter(function, any);
        duckdb_scalar_function_add_parameter(function, boolean);
        if (arity == 4) duckdb_scalar_function_add_parameter(function, any);
        duckdb_scalar_function_set_return_type(function, result);
        duckdb_scalar_function_set_special_handling(function);
        duckdb_scalar_function_set_function(function, repeat_scalar);
        duckdb_add_scalar_function_to_set(overloads, function);
        duckdb_destroy_scalar_function(&function);
    }
    duckdb_state state = duckdb_register_scalar_function_set(connection, overloads);
    duckdb_destroy_scalar_function_set(&overloads);
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&i64);
    duckdb_destroy_logical_type(&u64);
    duckdb_destroy_logical_type(&text);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&any);
    return state == DuckDBSuccess;
}
