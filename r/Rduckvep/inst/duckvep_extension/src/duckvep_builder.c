#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "core/duckvep_core_cells.h"
DUCKDB_EXTENSION_EXTERN

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Reports a builder failure; if the native budget refused an allocation on this
 * thread the message becomes an explicit capacity error. */
void
duckvep_builder_set_error(duckdb_function_info info, const char *message)
{
    char capacity[192], text[768];

    if (duckvep_budget_take_failure(capacity, sizeof capacity)) {
        (void)snprintf(text, sizeof text, "%s (%s)", capacity, message);
        message = text;
    }
    duckdb_scalar_function_set_error(info, message);
}

char *duckvep_builder_string(duckdb_string_t string) {
    return duckvep_core_string_copy(duckdb_string_t_data(&string), duckdb_string_t_length(string));
}

bool duckvep_builder_option_vectors(duckdb_function_info info, duckdb_vector vector,
                                    idx_t row, const char *const *names,
                                    const duckvep_option_kind *kinds, size_t count,
                                    duckdb_vector *values) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    if (validity && !duckdb_validity_row_is_valid(validity, row)) return true;
    duckdb_logical_type type = duckdb_vector_get_column_type(vector);
    if (duckdb_get_type_id(type) != DUCKDB_TYPE_STRUCT) {
        duckdb_destroy_logical_type(&type);
        duckdb_scalar_function_set_error(info, "DuckVEP builder: options must be a STRUCT");
        return false;
    }
    idx_t fields = duckdb_struct_type_child_count(type);
    for (idx_t i = 0; i < fields; i++) {
        char *key = duckdb_struct_type_child_name(type, i);
        duckdb_logical_type field_type = duckdb_struct_type_child_type(type, i);
        size_t at = duckvep_core_option_index(key, names, count);
        duckdb_type id = duckdb_get_type_id(field_type);
        bool integer = id >= DUCKDB_TYPE_TINYINT && id <= DUCKDB_TYPE_UBIGINT;
        bool numeric = integer || id == DUCKDB_TYPE_HUGEINT || id == DUCKDB_TYPE_DECIMAL ||
            id == DUCKDB_TYPE_FLOAT || id == DUCKDB_TYPE_DOUBLE;
        duckvep_core_field_type_t field = id == DUCKDB_TYPE_SQLNULL ? DUCKVEP_CORE_FIELD_SQLNULL :
            id == DUCKDB_TYPE_VARCHAR ? DUCKVEP_CORE_FIELD_VARCHAR :
            id == DUCKDB_TYPE_BOOLEAN ? DUCKVEP_CORE_FIELD_BOOLEAN :
            integer ? DUCKVEP_CORE_FIELD_INTEGER : numeric ? DUCKVEP_CORE_FIELD_NUMERIC :
            DUCKVEP_CORE_FIELD_OTHER;
        duckvep_core_option_kind_t core_kinds[8];
        for (size_t k = 0; k < count && k < 8; k++)
            core_kinds[k] = kinds[k] == DUCKVEP_OPTION_TEXT ? DUCKVEP_CORE_OPTION_TEXT :
                kinds[k] == DUCKVEP_OPTION_INTEGER ? DUCKVEP_CORE_OPTION_INTEGER :
                kinds[k] == DUCKVEP_OPTION_BOOLEAN ? DUCKVEP_CORE_OPTION_BOOLEAN :
                DUCKVEP_CORE_OPTION_NUMERIC;
        char message[256];
        if (!duckvep_core_option_permitted(key, at, count, core_kinds, field, message, sizeof(message))) {
            duckdb_scalar_function_set_error(info, message);
            duckdb_free(key);
            duckdb_destroy_logical_type(&field_type);
            duckdb_destroy_logical_type(&type);
            return false;
        }
        duckdb_vector child = duckdb_struct_vector_get_child(vector, i);
        values[at] = child;
        duckdb_free(key);
        duckdb_destroy_logical_type(&field_type);
    }
    duckdb_destroy_logical_type(&type);
    return true;
}

bool duckvep_builder_options(duckdb_function_info info, duckdb_vector vector,
                             idx_t row, const char *const *names, size_t count,
                             char **values) {
    duckvep_option_kind *kinds = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, count, sizeof(*kinds));
    duckdb_vector *fields = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, count, sizeof(*fields));
    if (!kinds || !fields) {
        duckvep_budget_free(kinds); duckvep_budget_free(fields);
        duckvep_builder_set_error(info, "DuckVEP builder: allocation failed");
        return false;
    }
    for (size_t i = 0; i < count; i++) kinds[i] = DUCKVEP_OPTION_TEXT;
    bool ok = duckvep_builder_option_vectors(info, vector, row, names, kinds, count, fields);
    for (size_t i = 0; ok && i < count; i++) {
        if (!fields[i]) continue;
        uint64_t *validity = duckdb_vector_get_validity(fields[i]);
        if (validity && !duckdb_validity_row_is_valid(validity, row)) continue;
        duckdb_logical_type type = duckdb_vector_get_column_type(fields[i]);
        bool string_field = duckdb_get_type_id(type) == DUCKDB_TYPE_VARCHAR;
        duckdb_destroy_logical_type(&type);
        if (!string_field) continue;
        values[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(fields[i]))[row]);
        if (!values[i]) {
            duckvep_builder_set_error(info, "DuckVEP builder: invalid option string or allocation failure");
            ok = false;
        }
    }
    duckvep_budget_free(kinds); duckvep_budget_free(fields);
    return ok;
}

bool duckvep_register_builder(duckdb_connection connection, const char *name,
                              idx_t required, duckdb_scalar_function_t callback) {
    duckdb_logical_type text = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type options = duckdb_create_logical_type(DUCKDB_TYPE_ANY);
    duckdb_scalar_function_set set = duckdb_create_scalar_function_set(name);
    for (idx_t arity = required; arity <= required + 1; arity++) {
        duckdb_scalar_function function = duckdb_create_scalar_function();
        duckdb_scalar_function_set_name(function, name);
        for (idx_t i = 0; i < required; i++) duckdb_scalar_function_add_parameter(function, text);
        if (arity > required) duckdb_scalar_function_add_parameter(function, options);
        duckdb_scalar_function_set_return_type(function, text);
        duckdb_scalar_function_set_special_handling(function);
        duckdb_scalar_function_set_function(function, callback);
        duckdb_add_scalar_function_to_set(set, function);
        duckdb_destroy_scalar_function(&function);
    }
    duckdb_state status = duckdb_register_scalar_function_set(connection, set);
    duckdb_destroy_scalar_function_set(&set);
    duckdb_destroy_logical_type(&options);
    duckdb_destroy_logical_type(&text);
    return status == DuckDBSuccess;
}
