#include "duckvep_core_sql.h"
#include "kernel/src/duckvep_budget.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static bool reserve(duckvep_sql_text *text, size_t extra) {
    if (extra > SIZE_MAX - text->length - 1) return false;
    size_t needed = text->length + extra + 1;
    if (needed <= text->capacity) return true;
    size_t capacity = text->capacity ? text->capacity : 128;
    while (capacity < needed) {
        if (capacity > SIZE_MAX / 2) { capacity = needed; break; }
        capacity *= 2;
    }
    char *data = duckvep_budget_realloc(DUCKVEP_OWNER_CONTROL, text->data, capacity);
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
    duckvep_budget_free(text->data);
    *text = (duckvep_sql_text){0};
}

/* A qualified name has exactly one schema separator; each component is quoted. */
bool duckvep_sql_relation(duckvep_sql_text *sql, const char *name) {
    const char *dot = strchr(name, '.');
    if (!*name || (dot && (!dot[1] || dot == name || strchr(dot + 1, '.')))) return false;
    if (!dot) return duckvep_sql_identifier(sql, name);
    size_t size = (size_t)(dot - name);
    char *schema = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, size + 1);
    if (!schema) return false;
    memcpy(schema, name, size); schema[size] = 0;
    bool ok = duckvep_sql_identifier(sql, schema) && duckvep_sql_append(sql, ".") &&
        duckvep_sql_identifier(sql, dot + 1);
    duckvep_budget_free(schema);
    return ok;
}

char *duckvep_core_string_copy(const char *data, size_t length) {
    if (memchr(data, 0, length)) return NULL;
    char *copy = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
    if (copy) {
        memcpy(copy, data, length);
        copy[length] = '\0';
    }
    return copy;
}
