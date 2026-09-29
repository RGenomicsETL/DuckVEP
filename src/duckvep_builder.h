#ifndef DUCKVEP_BUILDER_H
#define DUCKVEP_BUILDER_H

#include "duckdb_extension.h"
#include <stdbool.h>
#include <stddef.h>

/* SQL builders own their text until the output vector has copied it. */
typedef struct {
    char *data;
    size_t length;
    size_t capacity;
} duckvep_sql_text;

bool duckvep_sql_append(duckvep_sql_text *text, const char *part);
bool duckvep_sql_identifier(duckvep_sql_text *text, const char *name);
bool duckvep_sql_literal(duckvep_sql_text *text, const char *value);
bool duckvep_sql_relation(duckvep_sql_text *text, const char *name);
void duckvep_sql_free(duckvep_sql_text *text);

/* NULL options use defaults; every supplied field must have its declared type. */
typedef enum { DUCKVEP_OPTION_TEXT, DUCKVEP_OPTION_INTEGER, DUCKVEP_OPTION_NUMERIC, DUCKVEP_OPTION_BOOLEAN } duckvep_option_kind;
bool duckvep_builder_option_vectors(duckdb_function_info info, duckdb_vector vector,
                                    idx_t row, const char *const *names,
                                    const duckvep_option_kind *kinds, size_t count,
                                    duckdb_vector *values);

bool duckvep_builder_options(duckdb_function_info info, duckdb_vector vector,
                             idx_t row, const char *const *names, size_t count,
                             char **values);
char *duckvep_builder_string(duckdb_string_t string);
void duckvep_builder_set_error(duckdb_function_info info, const char *message);
bool duckvep_register_builder(duckdb_connection connection, const char *name,
                              idx_t required, duckdb_scalar_function_t callback);

#endif
