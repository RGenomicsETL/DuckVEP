#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_builder.h"
DUCKDB_EXTENSION_EXTERN

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

static bool huge_fractional(duckdb_hugeint n, uint8_t scale) {
    uint64_t high = (uint64_t)n.upper, low = n.lower;
    if (n.upper < 0) {
        low = ~low + 1;
        high = ~high + (low == 0);
    }
    uint32_t limbs[4] = {(uint32_t)(high >> 32), (uint32_t)high,
        (uint32_t)(low >> 32), (uint32_t)low};
    for (uint8_t digit = 0; digit < scale; digit++) {
        uint64_t remainder = 0;
        for (idx_t i = 0; i < 4; i++) {
            uint64_t value = (remainder << 32) | limbs[i];
            limbs[i] = (uint32_t)(value / 10);
            remainder = value % 10;
        }
        if (remainder) return true;
    }
    return false;
}

static bool numeric(duckdb_vector vector, duckdb_logical_type type, idx_t at,
                    long double *number, bool *fractional) {
    if (!vector || !valid(vector, at)) return false;
    duckdb_type id = duckdb_get_type_id(type);
    uint8_t scale = 0;
    if (id == DUCKDB_TYPE_DECIMAL) {
        scale = duckdb_decimal_scale(type);
        id = duckdb_decimal_internal_type(type);
    }
    void *data = duckdb_vector_get_data(vector);
    switch (id) {
    case DUCKDB_TYPE_TINYINT: *number = ((int8_t *)data)[at]; break;
    case DUCKDB_TYPE_SMALLINT: *number = ((int16_t *)data)[at]; break;
    case DUCKDB_TYPE_INTEGER: *number = ((int32_t *)data)[at]; break;
    case DUCKDB_TYPE_BIGINT: *number = ((int64_t *)data)[at]; break;
    case DUCKDB_TYPE_UTINYINT: *number = ((uint8_t *)data)[at]; break;
    case DUCKDB_TYPE_USMALLINT: *number = ((uint16_t *)data)[at]; break;
    case DUCKDB_TYPE_UINTEGER: *number = ((uint32_t *)data)[at]; break;
    case DUCKDB_TYPE_UBIGINT: *number = ((uint64_t *)data)[at]; break;
    case DUCKDB_TYPE_FLOAT: *number = ((float *)data)[at]; break;
    case DUCKDB_TYPE_DOUBLE: *number = ((double *)data)[at]; break;
    case DUCKDB_TYPE_HUGEINT: {
        duckdb_hugeint n = ((duckdb_hugeint *)data)[at];
        *number = (long double)n.upper * 18446744073709551616.0L + n.lower;
        if (scale) *fractional = huge_fractional(n, scale);
        break;
    }
    default: return false;
    }
    if (scale) {
        long double base = 1;
        for (uint8_t i = 0; i < scale; i++) base *= 10;
        if (id != DUCKDB_TYPE_HUGEINT)
            *fractional = fmodl(*number, base) != 0;
        *number /= base;
    } else *fractional = isfinite(*number) && truncl(*number) != *number;
    return true;
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

static bool dna(duckdb_string_t unit) {
    const char *data = duckdb_string_t_data(&unit);
    uint32_t length = duckdb_string_t_length(unit);
    if (!length) return false;
    for (uint32_t i = 0; i < length; i++) {
        char c = data[i];
        if (!strchr("ACGTRYSWKMBDHVNacgtryswkmbdhvn", c) || !c) return false;
    }
    return true;
}

static void set_error(duckdb_function_info info, const char *message) {
    duckvep_builder_set_error(info, message);
}

static void repeat_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector args[4] = {0};
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    repeat_axis axes[2] = {{0}};
    if (!axis_init(args[0], &axes[0]) || !axis_init(args[1], &axes[1])) {
        set_error(info, "duckvep_repeat_alleles: expected lists of {unit VARCHAR, count numeric}");
        axis_destroy(&axes[0]); axis_destroy(&axes[1]);
        return;
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
        duckdb_logical_type cap_type = cap_vector ? duckdb_vector_get_column_type(cap_vector) : NULL;
        if (!valid(args[2], row)) {
            set_error(info, "duckvep_repeat_alleles: sequence_exact is required");
            break;
        }
        long double cap = 5000;
        bool cap_frac = false;
        bool cap_ok = !cap_vector || (numeric(cap_vector, cap_type, row, &cap, &cap_frac) &&
            isfinite(cap) && cap >= 0 && cap <= INT32_MAX && !cap_frac);
        if (cap_type) duckdb_destroy_logical_type(&cap_type);
        if (!cap_ok) {
            set_error(info, "duckvep_repeat_alleles: max_allele_bases must be an integer from 0 through 2147483647");
            break;
        }
        bool incomplete = false, fractional = false, invalid_unit = false, invalid_count = false;
        long double required[2] = {0, 0};
        for (idx_t axis = 0; axis < 2; axis++) {
            repeat_axis *part = &axes[axis];
            if (!part->entries || !valid(part->list, row) || !part->units) { incomplete = true; continue; }
            duckdb_list_entry list = part->entries[row];
            for (idx_t i = 0; i < list.length; i++) {
                idx_t at = list.offset + i;
                if (!valid(part->records, at) || !valid(part->units, at) ||
                    !valid(part->counts, at)) { incomplete = true; continue; }
                duckdb_string_t unit = part->strings[at];
                if (!dna(unit)) invalid_unit = true;
                long double n = 0;
                bool frac = false;
                if (!numeric(part->counts, part->count_type, at, &n, &frac)) {
                    incomplete = true;
                    continue;
                }
                if (!isfinite(n) || n < 0) invalid_count = true;
                if (frac) fractional = true;
                required[axis] += duckdb_string_t_length(unit) * n;
            }
        }
        if (invalid_unit || invalid_count) {
            set_error(info, invalid_unit ? "duckvep_repeat_alleles: repeat units must contain non-empty IUPAC DNA" :
                "duckvep_repeat_alleles: repeat counts must be finite and nonnegative");
            break;
        }
        const char *status = !exact[row] ? "summary_only" : incomplete ? "incomplete_input" :
            fractional ? "nonintegral_count" : "ok";
        if (strcmp(status, "ok") == 0) {
            idx_t bad = required[0] > cap ? 0 : required[1] > cap ? 1 : 2;
            if (bad != 2) {
                char error[240];
                /* Print as double: MinGW's 80-bit long double does not match the
                 * Windows C runtime's printf, which reads long double as double. */
                snprintf(error, sizeof(error), "duckvep_repeat_alleles: %s requires %.6e bases which exceeds max_allele_bases=%.0f",
                    bad ? "alternate" : "reference", (double)required[bad], (double)cap);
                set_error(info, error);
                break;
            }
        }
        for (idx_t i = 0; i < 7; i++)
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[i]), row);
        duckdb_vector_assign_string_element(fields[6], row, status);
        if (strcmp(status, "ok") != 0) {
            for (idx_t i = 0; i < 6; i++)
                duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[i]), row);
            continue;
        }
        for (idx_t axis = 0; axis < 2; axis++) {
            repeat_axis *part = &axes[axis];
            size_t length = (size_t)required[axis];
            char *text = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
            if (!text) {
                set_error(info, "duckvep_repeat_alleles: allocation failed");
                if (cap_type) duckdb_destroy_logical_type(&cap_type);
                axis_destroy(&axes[0]); axis_destroy(&axes[1]);
                return;
            }
            size_t written = 0;
            duckdb_list_entry list = part->entries[row];
            for (idx_t i = 0; i < list.length; i++) {
                idx_t at = list.offset + i;
                duckdb_string_t unit = part->strings[at];
                size_t width = duckdb_string_t_length(unit);
                long double n = 0;
                bool frac = false;
                numeric(part->counts, part->count_type, at, &n, &frac);
                for (size_t j = 0; j < (size_t)n; j++) {
                    memcpy(text + written, duckdb_string_t_data(&unit), width);
                    written += width;
                }
            }
            text[written] = '\0';
            duckdb_vector_assign_string_element_len(fields[axis], row, text, written);
            duckvep_budget_free(text);
        }
        ((uint64_t *)duckdb_vector_get_data(fields[2]))[row] = (uint64_t)required[0];
        ((uint64_t *)duckdb_vector_get_data(fields[3]))[row] = (uint64_t)required[1];
        ((int64_t *)duckdb_vector_get_data(fields[4]))[row] = (int64_t)required[1] - (int64_t)required[0];
        duckdb_vector_assign_string_element(fields[5], row, required[1] > required[0] ? "GAIN" :
            required[1] < required[0] ? "LOSS" : "NEUTRAL");
    }
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
