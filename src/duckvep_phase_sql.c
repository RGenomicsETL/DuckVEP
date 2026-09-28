/* GT phase assignments use DuckDB list storage and the native lane reducer. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckvep_phase.h"
#include "duckvep_sql.h"
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

static bool phase_allele(duckdb_vector vector, duckdb_type type, uint8_t scale, idx_t at, int32_t *out) {
    void *data = duckdb_vector_get_data(vector);
    long double number;
    switch (type) {
    case DUCKDB_TYPE_TINYINT: number = ((int8_t *)data)[at]; break;
    case DUCKDB_TYPE_SMALLINT: number = ((int16_t *)data)[at]; break;
    case DUCKDB_TYPE_INTEGER: number = ((int32_t *)data)[at]; break;
    case DUCKDB_TYPE_BIGINT: number = ((int64_t *)data)[at]; break;
    case DUCKDB_TYPE_UTINYINT: number = ((uint8_t *)data)[at]; break;
    case DUCKDB_TYPE_USMALLINT: number = ((uint16_t *)data)[at]; break;
    case DUCKDB_TYPE_UINTEGER: number = ((uint32_t *)data)[at]; break;
    case DUCKDB_TYPE_UBIGINT: number = ((uint64_t *)data)[at]; break;
    case DUCKDB_TYPE_FLOAT: number = ((float *)data)[at]; break;
    case DUCKDB_TYPE_DOUBLE: number = ((double *)data)[at]; break;
    case DUCKDB_TYPE_HUGEINT: {
        duckdb_hugeint value = ((duckdb_hugeint *)data)[at];
        number = (long double)value.upper * 18446744073709551616.0L + value.lower;
        break;
    }
    case DUCKDB_TYPE_VARCHAR: {
        duckdb_string_t string = ((duckdb_string_t *)data)[at];
        size_t length = duckdb_string_t_length(string);
        if (!length || length >= 32) return false;
        char text[32];
        memcpy(text, duckdb_string_t_data(&string), length);
        text[length] = '\0';
        char *end;
        long value = strtol(text, &end, 10);
        if (end != text + length || value < INT32_MIN || value > INT32_MAX) return false;
        *out = (int32_t)value;
        return true;
    }
    default: return false;
    }
    for (uint8_t i = 0; i < scale; i++) number /= 10;
    if (!isfinite(number) || number < INT32_MIN - 0.5L || number >= INT32_MAX + 0.5L)
        return false;
    *out = (int32_t)roundl(number);
    return true;
}

static bool phase_flag(duckdb_vector vector, duckdb_type type, idx_t at, bool *result) {
    if (type == DUCKDB_TYPE_BOOLEAN) {
        *result = ((bool *)duckdb_vector_get_data(vector))[at];
        return true;
    }
    if (type == DUCKDB_TYPE_VARCHAR) {
        duckdb_string_t string = ((duckdb_string_t *)duckdb_vector_get_data(vector))[at];
        const char *data = duckdb_string_t_data(&string);
        uint32_t length = duckdb_string_t_length(string);
        if (length == 4 || length == 5) {
            const char *expected = length == 4 ? "true" : "false";
            bool matches = true;
            for (uint32_t i = 0; i < length; i++) {
                char c = data[i];
                if (c >= 'A' && c <= 'Z') c = (char)(c + ('a' - 'A'));
                if (c != expected[i]) { matches = false; break; }
            }
            if (matches) { *result = length == 4; return true; }
        }
        if (length == 1 && (*data == '1' || *data == '0')) {
            *result = *data == '1'; return true;
        }
        return false;
    }
    int32_t value = 0;
    if (!phase_allele(vector, type, 0, at, &value)) return false;
    *result = value != 0;
    return true;
}

static void phase_scalar(duckdb_function_info info, duckdb_data_chunk input,
                         duckdb_vector output) {
    duckdb_vector alleles = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector phases = duckdb_data_chunk_get_vector(input, 1);
    idx_t arity = duckdb_data_chunk_get_column_count(input);
    duckdb_vector third = arity >= 3 ? duckdb_data_chunk_get_vector(input, 2) : NULL;
    duckdb_vector ps = third;
    duckdb_vector policies = arity == 4 ? duckdb_data_chunk_get_vector(input, 3) : NULL;
    if (arity == 3) {
        duckdb_logical_type type = duckdb_vector_get_column_type(third);
        if (duckdb_get_type_id(type) == DUCKDB_TYPE_VARCHAR) {
            policies = third;
            ps = NULL;
        }
        duckdb_destroy_logical_type(&type);
    }
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
    int64_t *sets = ps ? duckdb_vector_get_data(ps) : NULL;
    duckdb_string_t *policy_names = policies ? duckdb_vector_get_data(policies) : NULL;
    idx_t rows = duckdb_data_chunk_get_size(input), total = 0u;

    for (idx_t row = 0u; row < rows; row++) {
        bool have_gt = allele_list && phase_valid(alleles, row), have_phase = phase_list && phase_valid(phases, row);
        if ((!have_gt && have_phase) ||
            (have_gt && have_phase && allele_lists[row].length != phase_lists[row].length)) {
            duckdb_scalar_function_set_error(info, "duckvep_phase_call: allele and phase lists must have equal length");
            return;
        }
        if (!have_gt) continue;
        idx_t count = allele_lists[row].length;
        if (!count || count > UINT16_MAX || count > UINT64_MAX - total) {
            duckdb_scalar_function_set_error(info, "duckvep_phase_call: ploidy must be between 1 and 65535");
            return;
        }
        total += count;
    }
    if (duckdb_list_vector_reserve(output, total) != DuckDBSuccess ||
        duckdb_list_vector_set_size(output, total) != DuckDBSuccess) {
        duckdb_scalar_function_set_error(info, "duckvep_phase_call: could not reserve output allele slots");
        return;
    }
    duckdb_vector_ensure_validity_writable(output);
    duckdb_list_entry *lists = duckdb_vector_get_data(output);
    duckdb_vector records = duckdb_list_vector_get_child(output), fields[7];
    for (idx_t i = 0u; i < 7u; i++) {
        fields[i] = duckdb_struct_vector_get_child(records, i);
        duckdb_vector_ensure_validity_writable(fields[i]);
    }
    idx_t at = 0u;
    for (idx_t row = 0u; row < rows; row++) {
        duckvep_phase_policy_t policy = DUCKVEP_PHASE_STRICT;
        if (policies && phase_valid(policies, row)) {
            const char *name = duckdb_string_t_data(&policy_names[row]);
            uint32_t length = duckdb_string_t_length(policy_names[row]);
            if (length == 13u && !memcmp(name, "vep116_compat", 13u)) {
                policy = DUCKVEP_PHASE_VEP116_COMPAT;
            } else if (length != 6u || memcmp(name, "strict", 6u)) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: phase_policy must be 'strict' or 'vep116_compat'");
                return;
            }
        }
        lists[row] = (duckdb_list_entry){at, 0u};
        if (!allele_list || !phase_valid(alleles, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        duckvep_phase_summary_t summary = {0};
        idx_t count = allele_lists[row].length, base = allele_lists[row].offset;
        bool have_phase = phase_list && phase_valid(phases, row);
        idx_t phase_base = have_phase ? phase_lists[row].offset : 0u;
        for (idx_t slot = 0u; slot < count; slot++) {
            bool called = !null_alleles && phase_valid(allele_values, base + slot);
            int32_t allele = -1;
            if (called && !phase_allele(allele_values, allele_id, allele_scale, base + slot, &allele)) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: called allele indices must be non-negative INTEGER values");
                return;
            }
            bool phased = false;
            if (have_phase && !null_phases && phase_valid(phase_values, phase_base + slot) &&
                !phase_flag(phase_values, phase_id, phase_base + slot, &phased)) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: phase_before elements must be BOOLEAN values");
                return;
            }
            if ((called && allele < 0) ||
                duckvep_phase_observe(&summary, allele, phased) != DUCKVEP_PHASE_OK) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: called allele indices must be non-negative INTEGER values");
                return;
            }
        }
        lists[row].length = count;
        uint16_t called_before = 0u;
        for (idx_t slot = 0u; slot < count; slot++, at++) {
            bool called = !null_alleles && phase_valid(allele_values, base + slot);
            int32_t allele = -1;
            if (called && !phase_allele(allele_values, allele_id, allele_scale, base + slot, &allele)) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: called allele indices must be non-negative INTEGER values");
                return;
            }
            bool phased = false;
            if (have_phase && !null_phases && phase_valid(phase_values, phase_base + slot) &&
                !phase_flag(phase_values, phase_id, phase_base + slot, &phased)) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: phase_before elements must be BOOLEAN values");
                return;
            }
            duckvep_phase_assignment_t assignment;
            if (duckvep_phase_assign(&summary, (uint16_t)(slot + 1u), called_before, allele, phased,
                                     policy, &assignment) != DUCKVEP_PHASE_OK) {
                duckdb_scalar_function_set_error(info, "duckvep_phase_call: invalid decoded phase state");
                return;
            }
            if (called) called_before++;
            for (idx_t i = 0u; i < 7u; i++)
                duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[i]), at);
            ((uint16_t *)duckdb_vector_get_data(fields[0]))[at] = (uint16_t)(slot + 1u);
            ((int32_t *)duckdb_vector_get_data(fields[1]))[at] = allele;
            if (!called) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[1]), at);
            ((uint16_t *)duckdb_vector_get_data(fields[2]))[at] = assignment.lane;
            if (!assignment.lane) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[2]), at);
            ((uint16_t *)duckdb_vector_get_data(fields[3]))[at] = summary.ploidy;
            if (assignment.scope == DUCKVEP_PHASE_SET && ps && phase_valid(ps, row)) {
                ((int64_t *)duckdb_vector_get_data(fields[4]))[at] = sets[row];
            } else {
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[4]), at);
            }
            const char *scope = assignment.scope == DUCKVEP_PHASE_SET ? "phase_set" :
                assignment.scope == DUCKVEP_PHASE_ALL_SETS ? "all_phase_sets" :
                assignment.scope == DUCKVEP_PHASE_ALLELE_SLOT ? "allele_slot" : "unresolved";
            const char *status = assignment.status == DUCKVEP_PHASE_CALLED ? "called" :
                assignment.status == DUCKVEP_PHASE_MISSING ? "missing" : "unphased";
            duckdb_vector_assign_string_element(fields[5], at, scope);
            duckdb_vector_assign_string_element(fields[6], at, status);
        }
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
    for (idx_t i = 0u; i < 7u; i++)
        fields[i] = duckdb_vector_get_data(duckdb_struct_vector_get_child(output, i));
    duckdb_vector_ensure_validity_writable(output);
    for (idx_t row = 0u; row < duckdb_data_chunk_get_size(input); row++) {
        if (!phase_valid(text, row) || !phase_valid(count_vector, row)) {
            duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
            continue;
        }
        duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
        duckvep_raw_gt_t call = {0};
        duckvep_raw_gt_status_t status = duckvep_phase_parse_vep116_raw(
            (const uint8_t *)duckdb_string_t_data(&strings[row]),
            duckdb_string_t_length(strings[row]), counts[row], &call);
        fields[0][row] = (uint32_t)status;
        fields[1][row] = call.allele_index[0]; fields[2][row] = call.allele_index[1];
        fields[3][row] = call.parsed_slots; fields[4][row] = call.source_ploidy;
        fields[5][row] = call.source_has_missing; fields[6][row] = (uint32_t)call.disposition;
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
        ranks[row] = duckvep_haplotype_record_order(counts[row], ordinals[row]);
        if (!ranks[row]) {
            duckdb_scalar_function_set_error(info, "duckvep_haplotypes: invalid source-buffer ordinal");
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
    for (idx_t arity = 2; arity <= 4; arity++) {
        idx_t alternatives = arity == 3 ? 2 : 1;
        for (idx_t alternative = 0; alternative < alternatives; alternative++) {
            duckdb_scalar_function function = duckdb_create_scalar_function();
            duckdb_scalar_function_set_name(function, "duckvep_phase_call");
            duckdb_scalar_function_add_parameter(function, any);
            duckdb_scalar_function_add_parameter(function, any);
            if (arity >= 3) duckdb_scalar_function_add_parameter(function,
                alternative ? varchar : bigint);
            if (arity == 4) duckdb_scalar_function_add_parameter(function, varchar);
            duckdb_scalar_function_set_return_type(function, result);
            duckdb_scalar_function_set_special_handling(function);
            duckdb_scalar_function_set_function(function, phase_scalar);
            duckdb_add_scalar_function_to_set(overloads, function);
            duckdb_destroy_scalar_function(&function);
        }
    }
    duckdb_state state = duckdb_register_scalar_function_set(connection, overloads);
    duckdb_destroy_scalar_function_set(&overloads);
    duckdb_destroy_logical_type(&any);
    duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&record);
    duckdb_destroy_logical_type(&varchar);
    duckdb_destroy_logical_type(&bigint);
    duckdb_destroy_logical_type(&ushort);
    duckdb_destroy_logical_type(&integer);
    return state == DuckDBSuccess;
}
