#include "duckvep_builder.h"
DUCKDB_EXTENSION_EXTERN

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static bool reserve(duckvep_sql_text *text, size_t extra) {
    if (extra > SIZE_MAX - text->length - 1) return false;
    size_t needed = text->length + extra + 1;
    if (needed <= text->capacity) return true;
    size_t capacity = text->capacity ? text->capacity : 128;
    while (capacity < needed) {
        if (capacity > SIZE_MAX / 2) { capacity = needed; break; }
        capacity *= 2;
    }
    char *data = realloc(text->data, capacity);
    if (!data) return false;
    text->data = data;
    text->capacity = capacity;
    return true;
}

bool duckvep_sql_append(duckvep_sql_text *text, const char *part) {
    size_t size = strlen(part);
    if (!reserve(text, size)) return false;
    memcpy(text->data + text->length, part, size + 1);
    text->length += size;
    return true;
}

static bool quoted(duckvep_sql_text *text, const char *part, char mark) {
    if (!part) return duckvep_sql_append(text, "NULL");
    size_t size = strlen(part), duplicates = 0;
    for (size_t i = 0; i < size; i++) if (part[i] == mark) duplicates++;
    if (size > SIZE_MAX - duplicates - 2 || !reserve(text, size + duplicates + 2)) return false;
    text->data[text->length++] = mark;
    for (size_t i = 0; i < size; i++) {
        text->data[text->length++] = part[i];
        if (part[i] == mark) text->data[text->length++] = mark;
    }
    text->data[text->length++] = mark;
    text->data[text->length] = '\0';
    return true;
}

bool duckvep_sql_identifier(duckvep_sql_text *text, const char *name) {
    return name && quoted(text, name, '"');
}
bool duckvep_sql_literal(duckvep_sql_text *text, const char *value) {
    return quoted(text, value, '\'');
}
void duckvep_sql_free(duckvep_sql_text *text) {
    free(text->data);
    *text = (duckvep_sql_text){0};
}

char *duckvep_builder_string(duckdb_string_t string) {
    size_t length = duckdb_string_t_length(string);
    const char *data = duckdb_string_t_data(&string);
    if (memchr(data, 0, length)) return NULL;
    char *copy = malloc(length + 1);
    if (copy) {
        memcpy(copy, data, length);
        copy[length] = '\0';
    }
    return copy;
}

bool duckvep_builder_options(duckdb_function_info info, duckdb_vector vector,
                             idx_t row, const char *const *names, size_t count,
                             char **values) {
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
        size_t at = 0;
        while (at < count && strcmp(key, names[at]) != 0) at++;
        duckdb_type id = duckdb_get_type_id(field_type);
        if (at == count || (id != DUCKDB_TYPE_VARCHAR && id != DUCKDB_TYPE_SQLNULL)) {
            char message[256];
            snprintf(message, sizeof(message), "DuckVEP builder: %s option '%s'",
                     at == count ? "unknown" : "expected VARCHAR for", key);
            duckdb_scalar_function_set_error(info, message);
            duckdb_free(key);
            duckdb_destroy_logical_type(&field_type);
            duckdb_destroy_logical_type(&type);
            return false;
        }
        duckdb_vector child = duckdb_struct_vector_get_child(vector, i);
        uint64_t *child_validity = duckdb_vector_get_validity(child);
        if (id == DUCKDB_TYPE_VARCHAR && (!child_validity || duckdb_validity_row_is_valid(child_validity, row))) {
            duckdb_string_t *strings = duckdb_vector_get_data(child);
            values[at] = duckvep_builder_string(strings[row]);
            if (!values[at]) {
                duckdb_scalar_function_set_error(info, "DuckVEP builder: invalid option string or allocation failure");
                duckdb_free(key);
                duckdb_destroy_logical_type(&field_type);
                duckdb_destroy_logical_type(&type);
                return false;
            }
        }
        duckdb_free(key);
        duckdb_destroy_logical_type(&field_type);
    }
    duckdb_destroy_logical_type(&type);
    return true;
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
