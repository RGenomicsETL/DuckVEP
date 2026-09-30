/* Host-neutral SQL text assembly for the builders (duckvep_*_sql). The text is
 * the product: both hosts must return the same bytes for the same arguments. */
#ifndef DUCKVEP_CORE_SQL_H
#define DUCKVEP_CORE_SQL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

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

/* A NUL-terminated budget-owned copy of a VARCHAR; NULL when it contains a NUL
 * byte or the allocation fails. Release with duckvep_budget_free. */
char *duckvep_core_string_copy(const char *data, size_t length);

#endif
