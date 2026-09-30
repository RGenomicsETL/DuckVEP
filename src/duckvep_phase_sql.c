/* GT phase assignments use DuckDB list storage and the native lane reducer. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "duckvep_v1_cells.h"

#include "duckvep_phase.h"
#include "duckvep_sql.h"
#include "duckvep_builder.h"
#include "core/duckvep_core_phase.h"
#include "kernel/src/duckvep_haplotype_stream.h"

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

static bool phase_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return !validity || duckdb_validity_row_is_valid(validity, row);
}

typedef struct {
    duckdb_function_info info;
    duckdb_vector allele_values, phase_values;
    duckdb_type allele_id, phase_id;
    uint8_t allele_scale;
    bool null_alleles, null_phases;
    idx_t allele_base, phase_base, at;
    duckdb_vector fields[7];
    bool has_phase_set;
    int64_t phase_set;
} phase_row_context;

static void phase_allele_cell(void *pointer, size_t slot, duckvep_cell_t *cell) {
    phase_row_context *context = pointer;
    if (context->null_alleles) {
        *cell = (duckvep_cell_t){0};
        return;
    }
    duckvep_v1_fill_cell(context->allele_values, context->allele_id, context->allele_scale,
                         context->allele_base + slot, cell);
}

static void phase_flag_cell(void *pointer, size_t slot, duckvep_cell_t *cell) {
    phase_row_context *context = pointer;
    if (context->null_phases) {
        *cell = (duckvep_cell_t){0};
        return;
    }
    duckvep_v1_fill_cell(context->phase_values, context->phase_id, 0,
                         context->phase_base + slot, cell);
}

static bool phase_emit(void *pointer, size_t slot, const duckvep_core_phase_slot_t *out) {
    phase_row_context *context = pointer;
    idx_t at = context->at + slot;
    duckdb_vector *fields = context->fields;
    for (idx_t i = 0u; i < 7u; i++)
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[i]), at);
    ((uint16_t *)duckdb_vector_get_data(fields[0]))[at] = out->input_slot;
    ((int32_t *)duckdb_vector_get_data(fields[1]))[at] = out->allele_index;
    if (!out->allele_called) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[1]), at);
    ((uint16_t *)duckdb_vector_get_data(fields[2]))[at] = out->lane;
    if (!out->lane) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[2]), at);
    ((uint16_t *)duckdb_vector_get_data(fields[3]))[at] = out->ploidy;
    if (out->phase_set_applies && context->has_phase_set) {
        ((int64_t *)duckdb_vector_get_data(fields[4]))[at] = context->phase_set;
    } else {
        duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[4]), at);
    }
    duckdb_vector_assign_string_element(fields[5], at, out->scope);
    duckdb_vector_assign_string_element(fields[6], at, out->status);
    return true;
}

static void phase_scalar(duckdb_function_info info, duckdb_data_chunk input,
                         duckdb_vector output) {
    duckdb_vector alleles = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector phases = duckdb_data_chunk_get_vector(input, 1);
    idx_t arity = duckdb_data_chunk_get_column_count(input);
    duckdb_vector options = arity == 3 ? duckdb_data_chunk_get_vector(input, 2) : NULL;
    const char *const keys[] = {"phase_set", "phase_policy"};
    const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_INTEGER, DUCKVEP_OPTION_TEXT};
    duckdb_logical_type allele_type = duckdb_vector_get_column_type(alleles);
    duckdb_logical_type phase_type = duckdb_vector_get_column_type(phases);
    bool allele_list = duckdb_get_type_id(allele_type) == DUCKDB_TYPE_LIST;
    bool phase_list = duckdb_get_type_id(phase_type) == DUCKDB_TYPE_LIST;
    duckdb_vector allele_values = allele_list ? duckdb_list_vector_get_child(alleles) : NULL;
    duckdb_vector phase_values = phase_list ? duckdb_list_vector_get_child(phases) : NULL;
    duckdb_list_entry *allele_lists = allele_list ? duckdb_vector_get_data(alleles) : NULL;
    duckdb_list_entry *phase_lists = phase_list ? duckdb_vector_get_data(phases) : NULL;
    duckdb_logical_type allele_child_type = allele_list ? duckdb_vector_get_column_type(allele_values) : NULL;
    duckdb_logical_type phase_child_type = phase_list ? duckdb_vector_get_column_type(phase_values) : NULL;
    duckdb_type allele_id = allele_child_type ? duckdb_get_type_id(allele_child_type) : DUCKDB_TYPE_SQLNULL;
    uint8_t allele_scale = allele_id == DUCKDB_TYPE_DECIMAL ? duckdb_decimal_scale(allele_child_type) : 0;
    if (allele_id == DUCKDB_TYPE_DECIMAL) allele_id = duckdb_decimal_internal_type(allele_child_type);
    duckdb_type phase_id = phase_child_type ? duckdb_get_type_id(phase_child_type) : DUCKDB_TYPE_SQLNULL;
    bool null_alleles = allele_id == DUCKDB_TYPE_SQLNULL;
    bool null_phases = phase_id == DUCKDB_TYPE_SQLNULL;
    if (allele_child_type) duckdb_destroy_logical_type(&allele_child_type);
    if (phase_child_type) duckdb_destroy_logical_type(&phase_child_type);
    duckdb_destroy_logical_type(&allele_type);
    duckdb_destroy_logical_type(&phase_type);
    idx_t rows = duckdb_data_chunk_get_size(input);
    size_t total = 0u;

    for (idx_t row = 0u; row < rows; row++) {
        bool have_gt = allele_list && phase_valid(alleles, row), have_phase = phase_list && phase_valid(phases, row);
        const char *message = duckvep_core_phase_check_row(have_gt, have_phase,
            have_gt ? allele_lists[row].length : 0, have_phase ? phase_lists[row].length : 0, &total);
        if (message) {
            duckdb_scalar_function_set_error(info, message);
            return;
        }
    }
    if (duckdb_list_vector_reserve(output, total) != DuckDBSuccess ||
        duckdb_list_vector_set_size(output, total) != DuckDBSuccess) {
        duckdb_scalar_function_set_error(info, "duckvep_phase_call: could not reserve output allele slots");
        return;
    }
    duckdb_vector_ensure_validity_writable(output);
    duckdb_list_entry *lists = duckdb_vector_get_data(output);
    duckdb_vector records = duckdb_list_vector_get_child(output);
    phase_row_context context = {info, allele_values, phase_values, allele_id, phase_id, allele_scale,
        null_alleles, null_phases, 0, 0, 0, {0}, false, 0};
    for (idx_t i = 0u; i < 7u; i++) {
        context.fields[i] = duckdb_struct_vector_get_child(records, i);
        duckdb_vector_ensure_validity_writable(context.fields[i]);
    }
    idx_t at = 0u;
    for (idx_t row = 0u; row < rows; row++) {
        duckdb_vector fields_option[2] = {NULL, NULL};
        if (options && !duckvep_builder_option_vectors(info, options, row, keys, kinds, 2, fields_option)) return;
        duckdb_vector ps = fields_option[0] && phase_valid(fields_option[0], row) ? fields_option[0] : NULL;
        duckdb_vector policies = fields_option[1] && phase_valid(fields_option[1], row) ? fields_option[1] : NULL;
        if (ps) {
            duckdb_logical_type type = duckdb_vector_get_column_type(ps);
            if (duckdb_get_type_id(type) == DUCKDB_TYPE_SQLNULL) ps = NULL;
            duckdb_destroy_logical_type(&type);
        }
        if (policies) {
            duckdb_logical_type type = duckdb_vector_get_column_type(policies);
            if (duckdb_get_type_id(type) == DUCKDB_TYPE_SQLNULL) policies = NULL;
            duckdb_destroy_logical_type(&type);
        }
        int64_t phase_set = 0;
        if (ps) {
            duckdb_logical_type set_type = duckdb_vector_get_column_type(ps);
            duckvep_cell_t cell;
            duckvep_v1_fill_cell(ps, duckdb_get_type_id(set_type), 0, row, &cell);
            duckdb_destroy_logical_type(&set_type);
            if (!duckvep_core_phase_set(&cell, &phase_set)) {
                duckdb_scalar_function_set_error(info, duckvep_core_phase_set_error);
                return;
            }
        }
        duckdb_string_t *policy_names = policies ? duckdb_vector_get_data(policies) : NULL;
        duckvep_phase_policy_t policy = DUCKVEP_PHASE_STRICT;
        if (policies &&
            !duckvep_core_phase_policy(duckdb_string_t_data(&policy_names[row]),
                                       duckdb_string_t_length(policy_names[row]), &policy)) {
            duckdb_scalar_function_set_error(info, duckvep_core_phase_policy_error);
            return;
        }
        lists[row] = (duckdb_list_entry){at, 0u};
        if (!allele_list || !phase_valid(alleles, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        idx_t count = allele_lists[row].length;
        bool have_phase = phase_list && phase_valid(phases, row);
        context.allele_base = allele_lists[row].offset;
        context.phase_base = have_phase ? phase_lists[row].offset : 0u;
        context.at = at;
        context.has_phase_set = ps != NULL;
        context.phase_set = phase_set;
        const char *message = NULL;
        duckvep_core_phase_result_t result = duckvep_core_phase_row(&(duckvep_core_phase_reader_t){
            &context, phase_allele_cell, phase_flag_cell, phase_emit}, count, have_phase, policy, &message);
        if (result != DUCKVEP_CORE_PHASE_OK) {
            duckdb_scalar_function_set_error(info, message);
            return;
        }
        lists[row].length = count;
        at += count;
    }
}

static void raw_gt_scalar(duckdb_function_info info, duckdb_data_chunk input,
    duckdb_vector output) {
    (void)info;
    duckdb_vector text = duckdb_data_chunk_get_vector(input, 0u);
    duckdb_vector count_vector = duckdb_data_chunk_get_vector(input, 1u);
    duckdb_string_t *strings = duckdb_vector_get_data(text);
    uint32_t *counts = duckdb_vector_get_data(count_vector);
    uint32_t *fields[7];
    duckdb_vector children[7];
    for (idx_t i = 0u; i < 7u; i++) {
        children[i] = duckdb_struct_vector_get_child(output, i);
        fields[i] = duckdb_vector_get_data(children[i]);
        duckdb_vector_ensure_validity_writable(children[i]);
    }
    duckdb_vector_ensure_validity_writable(output);
    for (idx_t row = 0u; row < duckdb_data_chunk_get_size(input); row++) {
        if (!phase_valid(text, row) || !phase_valid(count_vector, row)) {
            /* A NULL STRUCT row must null its children too, or field reads
             * and vector copies see uninitialized values. */
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            for (idx_t i = 0u; i < 7u; i++)
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(children[i]), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        uint32_t values[7];
        duckvep_core_raw_gt(duckdb_string_t_data(&strings[row]),
            duckdb_string_t_length(strings[row]), counts[row], values);
        for (idx_t i = 0u; i < 7u; i++) fields[i][row] = values[i];
    }
}

static void record_order_scalar(duckdb_function_info info, duckdb_data_chunk input,
    duckdb_vector output) {
    duckdb_vector count_vector = duckdb_data_chunk_get_vector(input, 0u);
    duckdb_vector ordinal_vector = duckdb_data_chunk_get_vector(input, 1u);
    uint64_t *counts = duckdb_vector_get_data(count_vector);
    uint64_t *ordinals = duckdb_vector_get_data(ordinal_vector);
    uint64_t *ranks = duckdb_vector_get_data(output);
    duckdb_vector_ensure_validity_writable(output);
    for (idx_t row = 0u; row < duckdb_data_chunk_get_size(input); row++) {
        if (!phase_valid(count_vector, row) || !phase_valid(ordinal_vector, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        ranks[row] = duckvep_core_record_order(counts[row], ordinals[row]);
        if (!ranks[row]) {
            duckdb_scalar_function_set_error(info, duckvep_core_record_order_error);
            return;
        }
    }
}

static bool register_raw_preparation(duckdb_connection connection) {
    duckdb_logical_type uinteger = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
    duckdb_logical_type ubigint = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type fields[7];
    for (idx_t i = 0u; i < 7u; i++) fields[i] = uinteger;
    const char *names[] = {"status", "allele0", "allele1", "parsed_slots", "source_ploidy",
        "source_has_missing", "disposition"};
    duckdb_logical_type result = duckdb_create_struct_type(fields, names, 7u);
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(function, "_duckvep_raw_gt");
    duckdb_scalar_function_add_parameter(function, varchar);
    duckdb_scalar_function_add_parameter(function, uinteger);
    duckdb_scalar_function_set_return_type(function, result);
    duckdb_scalar_function_set_function(function, raw_gt_scalar);
    duckdb_scalar_function_set_special_handling(function);
    duckdb_state state = duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    if (state == DuckDBSuccess) {
        function = duckdb_create_scalar_function();
        duckdb_scalar_function_set_name(function, "_duckvep_record_order");
        duckdb_scalar_function_add_parameter(function, ubigint);
        duckdb_scalar_function_add_parameter(function, ubigint);
        duckdb_scalar_function_set_return_type(function, ubigint);
        duckdb_scalar_function_set_function(function, record_order_scalar);
        duckdb_scalar_function_set_special_handling(function);
        state = duckdb_register_scalar_function(connection, function);
        duckdb_destroy_scalar_function(&function);
    }
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&ubigint);
    duckdb_destroy_logical_type(&uinteger);
    return state == DuckDBSuccess;
}

bool duckvep_register_phase_kernels(duckdb_connection connection) {
    if (!register_raw_preparation(connection)) {
        return false;
    }
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type ushort = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type types[] = {ushort, integer, ushort, ushort, bigint, varchar, varchar};
    const char *names[] = {"input_slot", "allele_index", "haplotype_lane", "ploidy",
        "phase_set", "phase_scope", "status"};
    duckdb_logical_type record = duckdb_create_struct_type(types, names, 7u);
    duckdb_logical_type result = duckdb_create_list_type(record);
    duckdb_logical_type any = duckdb_create_logical_type(DUCKDB_TYPE_ANY);
    duckdb_scalar_function_set overloads = duckdb_create_scalar_function_set("duckvep_phase_call");
    for (idx_t arity = 2; arity <= 3; arity++) {
            duckdb_scalar_function function = duckdb_create_scalar_function();
            duckdb_scalar_function_set_name(function, "duckvep_phase_call");
            duckdb_scalar_function_add_parameter(function, any);
            duckdb_scalar_function_add_parameter(function, any);
            if (arity == 3) duckdb_scalar_function_add_parameter(function, any);
            duckdb_scalar_function_set_return_type(function, result);
            duckdb_scalar_function_set_special_handling(function);
            duckdb_scalar_function_set_function(function, phase_scalar);
            duckdb_add_scalar_function_to_set(overloads, function);
            duckdb_destroy_scalar_function(&function);
    }
    duckdb_state state = duckdb_register_scalar_function_set(connection, overloads);
    duckdb_destroy_scalar_function_set(&overloads);
    duckdb_destroy_logical_type(&any);
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&record);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&ushort);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&integer);
    return state == DuckDBSuccess;
}
