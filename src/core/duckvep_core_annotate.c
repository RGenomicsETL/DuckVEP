#include "duckvep_core_annotate.h"

#include "kernel/src/duckvep_budget.h"

#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#include "duckvep_annotate_template.h"
#include "duckvep_projected_template.h"
#include "duckvep_projection_template.h"

const char duckvep_core_projected_model_empty[] =
    "duckvep_annotate_projected: model_name must be non-empty";
const char duckvep_core_annotate_failed[] =
    "duckvep_annotate_sql: invalid input or allocation failure";
const char duckvep_core_projection_failed[] =
    "duckvep_transcript_projection_sql: invalid table name, options, or allocation failure";

/* A relation name: an optional schema before the first dot, each part quoted. */
static bool projection_table(duckvep_sql_text *sql, const char *name)
{
    const char *dot = strchr(name, '.');
    if (dot) {
        size_t length = (size_t)(dot - name);
        char *schema = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
        if (!schema) return false;
        memcpy(schema, name, length);
        schema[length] = 0;
        bool ok = duckvep_sql_identifier(sql, schema) && duckvep_sql_append(sql, ".") &&
            duckvep_sql_identifier(sql, dot + 1);
        duckvep_budget_free(schema);
        return ok;
    }
    return duckvep_sql_identifier(sql, name);
}

static bool option_number(const duckvep_cell_t *cell, duckvep_sql_text *sql)
{
    if (!cell->valid || cell->kind == DUCKVEP_CELL_UNSUPPORTED)
        return duckvep_sql_append(sql, "NULL");
    char buffer[32];
    switch (cell->kind) {
    case DUCKVEP_CELL_TINYINT: case DUCKVEP_CELL_SMALLINT:
    case DUCKVEP_CELL_INTEGER: case DUCKVEP_CELL_BIGINT:
        snprintf(buffer, sizeof(buffer), "%" PRId64, cell->i);
        break;
    case DUCKVEP_CELL_UTINYINT: case DUCKVEP_CELL_USMALLINT:
    case DUCKVEP_CELL_UINTEGER: case DUCKVEP_CELL_UBIGINT:
        snprintf(buffer, sizeof(buffer), "%" PRIu64, cell->u);
        break;
    default: return false;
    }
    return duckvep_sql_append(sql, buffer);
}

static bool option_boolean(const duckvep_cell_t *cell, duckvep_sql_text *sql)
{
    if (!cell->valid || cell->kind == DUCKVEP_CELL_UNSUPPORTED)
        return duckvep_sql_append(sql, "NULL");
    return duckvep_sql_append(sql, cell->boolean ? "true" : "false");
}

bool duckvep_core_annotate_sql(bool projected, const char *events, const char *model,
    const duckvep_cell_t *const *options, duckvep_sql_text *out)
{
    static const char *const tokens[] = {"__DUCKVEP_EVENTS__", "__DUCKVEP_MODEL__",
        "__DUCKVEP_HGVS__", "__DUCKVEP_UPSTREAM__", "__DUCKVEP_DOWNSTREAM__", "__DUCKVEP_RICH__"};
    duckvep_sql_text values[6] = {{0}};
    bool ok = projection_table(&values[0], events) &&
        duckvep_sql_literal(&values[1], model) &&
        duckvep_sql_append(&values[2], "false") &&
        duckvep_sql_append(&values[3], "5000") &&
        duckvep_sql_append(&values[4], "5000") &&
        duckvep_sql_append(&values[5], "false");
    for (size_t i = 0; ok && i < (projected ? 2 : 4); i++) {
        if (!options || !options[i]) continue;
        size_t dest = projected ? i + 3 : (i == 3 ? 5 : i + 2);
        duckvep_sql_free(&values[dest]);
        ok = !projected && (i == 0 || i == 3) ? option_boolean(options[i], &values[dest])
                                              : option_number(options[i], &values[dest]);
    }
    size_t part_count = projected ? sizeof(duckvep_projected_parts) / sizeof(*duckvep_projected_parts)
                                  : sizeof(duckvep_annotate_parts) / sizeof(*duckvep_annotate_parts);
    for (size_t i = 0; ok && i < part_count; i++) {
        const char *part = projected ? duckvep_projected_parts[i] : duckvep_annotate_parts[i];
        while (ok && *part) {
            const char *mark = strstr(part, "__DUCKVEP_");
            if (!mark) { ok = duckvep_sql_append(out, part); break; }
            size_t length = (size_t)(mark - part);
            char *prefix = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
            if (!prefix) { ok = false; break; }
            memcpy(prefix, part, length); prefix[length] = 0;
            ok = duckvep_sql_append(out, prefix); duckvep_budget_free(prefix);
            size_t index = 0;
            while (index < 6 && strncmp(mark, tokens[index], strlen(tokens[index])) != 0) index++;
            if (index == 6) { ok = false; break; }
            if (ok) ok = duckvep_sql_append(out, values[index].data);
            part = mark + strlen(tokens[index]);
        }
    }
    for (size_t i = 0; i < 6; i++) duckvep_sql_free(&values[i]);
    return ok;
}

bool duckvep_core_projection_sql(const char *const names[3], duckvep_sql_text *out)
{
    bool ok = true;
    for (size_t i = 0; ok && i < sizeof(duckvep_projection_parts) / sizeof(*duckvep_projection_parts); i++) {
        const char *part = duckvep_projection_parts[i];
        while (ok && *part) {
            const char *mark = strstr(part, "__DUCKVEP_");
            if (!mark) { ok = duckvep_sql_append(out, part); break; }
            size_t length = (size_t)(mark - part);
            char *prefix = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
            if (!prefix) { ok = false; break; }
            memcpy(prefix, part, length); prefix[length] = 0;
            ok = duckvep_sql_append(out, prefix); duckvep_budget_free(prefix);
            const char *token = NULL; size_t which = 0;
            if (!strncmp(mark, "__DUCKVEP_EVENTS_TABLE__", strlen("__DUCKVEP_EVENTS_TABLE__"))) token = "__DUCKVEP_EVENTS_TABLE__";
            else if (!strncmp(mark, "__DUCKVEP_ANNOTATIONS_TABLE__", strlen("__DUCKVEP_ANNOTATIONS_TABLE__"))) { token = "__DUCKVEP_ANNOTATIONS_TABLE__"; which = 1; }
            else if (!strncmp(mark, "__DUCKVEP_TRANSCRIPTS_TABLE__", strlen("__DUCKVEP_TRANSCRIPTS_TABLE__"))) { token = "__DUCKVEP_TRANSCRIPTS_TABLE__"; which = 2; }
            else ok = false;
            if (ok) ok = projection_table(out, names[which]);
            if (token) part = mark + strlen(token);
        }
    }
    return ok;
}
