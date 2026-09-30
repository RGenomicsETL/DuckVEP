/* The SQL-builder family of the v2 host: native scalars that return SQL text
 * (duckvep_*_sql). The text is assembled by src/core; this file reads the
 * string and options-STRUCT arguments and writes the text back. The order of
 * the argument and option checks per row is the v1 host's, so that both hosts
 * fail identically. */
#include "host_v2_columns.h"

#include "core/duckvep_core_annotate.h"
#include "core/duckvep_core_ensembl.h"
#include "core/duckvep_core_lof.h"
#include "core/duckvep_core_prepare.h"

#define MAX_ARGS 10

typedef struct {
    duckdb_v2_error_info_handle *error;
    duckdb_v2_error_info_handle detail;
    idx_t rows;
    uint32_t argc;
    column args[MAX_ARGS];
    duckdb_v2_vector_handle output;
    duckdb_v2_arena_handle arena;
    void *data;
} call;

/* Opens the string arguments [0, strings) and leaves the options argument (if
 * any) for options_open; prepares the VARCHAR result. */
static bool call_open(duckdb_v2_scalar_function_exec_info_handle info, idx_t strings, call *c,
                      duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    bool success = false;
    uint64_t *validity = NULL;
    memset(c, 0, sizeof(*c));
    c->error = error;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &c->rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &c->argc, &detail));
    for (idx_t i = 0; i < strings && i < c->argc; ++i) {
        duckdb_v2_vector_handle vector = NULL;
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &vector, &detail));
        if (!column_open(vector, false, &c->args[i], error)) {
            goto cleanup;
        }
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &c->output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(c->output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(c->output, c->rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(c->output, &c->data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(c->output, &validity, &detail));
    for (idx_t row = 0; row < c->rows; ++row) {
        mark_valid(validity, row);
    }
    DUCKDB_CALL(duckdb_v2_vector_get_arena(c->output, &c->arena, &detail));
    success = true;
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

static duckdb_v2_vector_handle call_argument(duckdb_v2_scalar_function_exec_info_handle info,
                                             idx_t index, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle vector = NULL;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, index, &vector, &detail));
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return vector;
}

/* A budget-owned copy of a string argument; NULL for a NULL row or an embedded NUL. */
static char *copy_argument(const call *c, idx_t index, idx_t row) {
    if (!column_valid(&c->args[index], row)) {
        return NULL;
    }
    duckdb_v2_str text = string_at(&c->args[index].view, row);
    return duckvep_core_string_copy(text.ptr, text.len);
}

static bool write_sql(call *c, idx_t row, const duckvep_sql_text *sql) {
    duckdb_v2_error_info_handle detail = NULL;
    DUCKDB_V2_ERROR status = write_string(c->arena, &((duckdb_v2_bytes *)c->data)[row], sql->data,
                                          sql->length, &detail);
    if (status != DUCKDB_V2_ERROR_NONE) {
        copy_duckdb_error(*c->error, status, detail);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
    return status == DUCKDB_V2_ERROR_NONE;
}

static void plain_error(const call *c, const char *message) {
    set_error(*c->error, DUCKDB_V2_ERROR_INPUT_INVALID, message);
}

static void free_all(char **values, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        duckvep_budget_free(values[i]);
    }
}

/* The cell of option `at` for a row, from an opened options argument. */
static const duckvep_cell_t *option_or_null(const options *o, size_t at, idx_t row,
                                            duckvep_cell_t *storage) {
    return option_cell(o, at, row, storage) ? storage : NULL;
}

/* ---------------------------------------------------------------------------
 * duckvep_annotate_sql(events, model [, options]) and the projected variant
 * ------------------------------------------------------------------------- */

static void annotate_common(duckdb_v2_scalar_function_exec_info_handle info,
                            duckdb_v2_error_info_handle *error, bool projected) {
    static const char *const keys[] = {"hgvs", "upstream_distance", "downstream_distance", "rich"};
    static const duckvep_core_option_kind_t kinds[] = {
        DUCKVEP_CORE_OPTION_BOOLEAN, DUCKVEP_CORE_OPTION_INTEGER, DUCKVEP_CORE_OPTION_INTEGER,
        DUCKVEP_CORE_OPTION_BOOLEAN};
    static const char *const projected_keys[] = {"upstream_distance", "downstream_distance"};
    static const duckvep_core_option_kind_t projected_kinds[] = {DUCKVEP_CORE_OPTION_INTEGER,
                                                                 DUCKVEP_CORE_OPTION_INTEGER};
    call c;
    options option;
    if (!call_open(info, 2, &c, error)) {
        return;
    }
    if (!options_open(c.argc == 3 ? call_argument(info, 2, error) : NULL, projected ? projected_keys : keys,
                      projected ? projected_kinds : kinds, projected ? 2 : 4, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *names[2] = {0};
        bool ok = true;
        duckvep_cell_t cells[4];
        const duckvep_cell_t *values[4] = {0};
        duckvep_sql_text sql = {0};
        for (idx_t i = 0; i < 2; ++i) {
            names[i] = copy_argument(&c, i, row);
            if (!names[i]) {
                ok = false;
                break;
            }
        }
        if (ok && projected && names[1][0] == '\0') {
            plain_error(&c, duckvep_core_projected_model_empty);
            free_all(names, 2);
            return;
        }
        if (ok && c.argc == 3 && !options_check(&option, row, error)) {
            free_all(names, 2);
            return;
        }
        for (size_t i = 0; ok && i < (projected ? 2u : 4u); ++i) {
            values[i] = option_or_null(&option, i, row, &cells[i]);
        }
        if (ok) {
            ok = duckvep_core_annotate_sql(projected, names[0], names[1], values, &sql);
        }
        if (ok) {
            ok = write_sql(&c, row, &sql);
            if (!ok) {
                duckvep_sql_free(&sql);
                free_all(names, 2);
                return;
            }
        } else {
            report(*error, duckvep_core_annotate_failed);
        }
        duckvep_sql_free(&sql);
        free_all(names, 2);
        if (!ok) {
            return;
        }
    }
}

static void annotate_exec(duckdb_v2_scalar_function_exec_info_handle info,
                          duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    annotate_common(info, error, false);
}

static void projected_exec(duckdb_v2_scalar_function_exec_info_handle info,
                           duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    annotate_common(info, error, true);
}

/* ---------------------------------------------------------------------------
 * duckvep_transcript_projection_sql(events, annotations, transcripts [, options])
 * ------------------------------------------------------------------------- */

static void projection_exec(duckdb_v2_scalar_function_exec_info_handle info,
                            duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    static const char *const keys[] = {"unused"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_TEXT};
    call c;
    options option;
    if (!call_open(info, 3, &c, error)) {
        return;
    }
    if (!options_open(c.argc == 4 ? call_argument(info, 3, error) : NULL, keys, kinds, 0, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *names[3] = {0};
        bool ok = true;
        duckvep_sql_text sql = {0};
        for (idx_t i = 0; i < 3; ++i) {
            names[i] = copy_argument(&c, i, row);
            if (!names[i] || !*names[i]) {
                ok = false;
                break;
            }
        }
        if (ok && c.argc == 4 && !options_check(&option, row, error)) {
            free_all(names, 3);
            return;
        }
        if (ok) {
            ok = duckvep_core_projection_sql((const char *const *)names, &sql);
        }
        if (ok) {
            ok = write_sql(&c, row, &sql);
            if (!ok) {
                duckvep_sql_free(&sql);
                free_all(names, 3);
                return;
            }
        } else {
            report(*error, duckvep_core_projection_failed);
        }
        duckvep_sql_free(&sql);
        free_all(names, 3);
        if (!ok) {
            return;
        }
    }
}

/* ---------------------------------------------------------------------------
 * duckvep_ensembl_{regions,transcripts,regulation_features}_sql
 * ------------------------------------------------------------------------- */

static void ensembl_common(duckdb_v2_scalar_function_exec_info_handle info,
                           duckdb_v2_error_info_handle *error, duckvep_ensembl_kind_t kind) {
    static const char *const keys[] = {"species_id"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_INTEGER};
    idx_t required = duckvep_core_ensembl_required(kind);
    call c;
    options option;
    bool with_options;
    if (!call_open(info, required, &c, error)) {
        return;
    }
    with_options = c.argc > required;
    if (!options_open(with_options && required == 3 ? call_argument(info, required, error) : NULL, keys,
                      kinds, 1, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *values[3] = {0};
        bool ok = true;
        char species[40] = "1";
        duckvep_sql_text sql = {0};
        bool building;
        for (idx_t i = 0; i < required; ++i) {
            values[i] = copy_argument(&c, i, row);
            if (!values[i]) {
                break;
            }
        }
        for (idx_t i = 0; i < required; ++i) {
            if (!values[i] || !*values[i]) {
                ok = false;
            }
        }
        if (!ok) {
            plain_error(&c, duckvep_core_ensembl_names_required);
        }
        if (ok && with_options) {
            if (required == 3) {
                ok = options_check(&option, row, error);
            } else {
                plain_error(&c, duckvep_core_ensembl_no_options);
                ok = false;
            }
        }
        if (ok && required == 3 && with_options) {
            duckvep_cell_t cell;
            if (option_cell(&option, 0, row, &cell)) {
                ok = duckvep_core_ensembl_species(&cell, species);
                if (!ok) {
                    plain_error(&c, duckvep_core_ensembl_species_range);
                }
            }
        }
        building = ok;
        if (ok) {
            ok = duckvep_core_ensembl_sql(kind, (const char *const *)values, species, &sql);
        }
        if (ok) {
            ok = write_sql(&c, row, &sql);
        } else if (building) {
            plain_error(&c, duckvep_core_ensembl_failed(kind));
        }
        duckvep_sql_free(&sql);
        free_all(values, required);
        if (!ok) {
            return;
        }
    }
}

static void regions_exec(duckdb_v2_scalar_function_exec_info_handle info,
                         duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    ensembl_common(info, error, DUCKVEP_ENSEMBL_REGIONS);
}
static void transcripts_exec(duckdb_v2_scalar_function_exec_info_handle info,
                             duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    ensembl_common(info, error, DUCKVEP_ENSEMBL_TRANSCRIPTS);
}
static void regulation_exec(duckdb_v2_scalar_function_exec_info_handle info,
                            duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    ensembl_common(info, error, DUCKVEP_ENSEMBL_REGULATION);
}

/* ---------------------------------------------------------------------------
 * duckvep_model_receipt_sql(regions, model, 6 parameters [, options])
 * ------------------------------------------------------------------------- */

static void receipt_exec(duckdb_v2_scalar_function_exec_info_handle info,
                         duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    static const char *const keys[] = {"regulation_features_table"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_TEXT};
    call c;
    options option;
    if (!call_open(info, 8, &c, error)) {
        return;
    }
    if (!options_open(c.argc == 9 ? call_argument(info, 8, error) : NULL, keys, kinds, 1, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *values[8] = {0};
        char *option_table = NULL;
        duckvep_sql_text sql = {0};
        bool ok;
        for (idx_t i = 0; i < 8; ++i) {
            if (!column_valid(&c.args[i], row)) {
                continue;
            }
            values[i] = copy_argument(&c, i, row);
            if (!values[i]) {
                report(*error, "duckvep_model_receipt_sql: invalid argument string or allocation failure");
                free_all(values, 8);
                return;
            }
        }
        ok = c.argc != 9 || options_check(&option, row, error);
        if (ok && c.argc == 9) {
            duckvep_cell_t cell;
            if (option_cell(&option, 0, row, &cell) && cell.valid) {
                option_table = duckvep_core_string_copy(cell.text, cell.text_length);
                if (!option_table) {
                    report(*error, "DuckVEP builder: invalid option string or allocation failure");
                    ok = false;
                }
            }
        }
        if (ok && (!values[0] || !values[1])) {
            plain_error(&c, duckvep_core_receipt_tables_required);
            ok = false;
        }
        if (!ok) {
            free_all(values, 8);
            duckvep_budget_free(option_table);
            return;
        }
        ok = duckvep_core_receipt_sql((const char *const *)values, option_table, &sql);
        if (!ok) {
            report(*error, "duckvep_model_receipt_sql: allocation failed");
        } else {
            ok = write_sql(&c, row, &sql);
        }
        duckvep_sql_free(&sql);
        free_all(values, 8);
        duckvep_budget_free(option_table);
        if (!ok) {
            return;
        }
    }
}

/* ---------------------------------------------------------------------------
 * duckvep_prepare_sv_geometry_sql(events [, options]),
 * duckvep_prepare_expansionhunter_sql(events, reference [, options])
 * ------------------------------------------------------------------------- */

static void prepare_common(duckdb_v2_scalar_function_exec_info_handle info,
                           duckdb_v2_error_info_handle *error, idx_t relations) {
    static const char *const keys[] = {"unused"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_TEXT};
    call c;
    options option;
    if (!call_open(info, relations, &c, error)) {
        return;
    }
    if (!options_open(c.argc == relations + 1 ? call_argument(info, relations, error) : NULL, keys, kinds,
                      0, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *names[2] = {0};
        bool ok = true;
        duckvep_sql_text sql = {0};
        /* The options are checked first, then the relation names. */
        if (c.argc == relations + 1 && !options_check(&option, row, error)) {
            return;
        }
        for (idx_t i = 0; i < relations; ++i) {
            names[i] = copy_argument(&c, i, row);
            if (!names[i]) {
                ok = false;
                break;
            }
        }
        if (ok) {
            ok = relations == 1 ? duckvep_core_prepare_sv_sql(names[0], &sql)
                                : duckvep_core_prepare_expansionhunter_sql(names[0], names[1], &sql);
        }
        if (ok) {
            ok = write_sql(&c, row, &sql);
            if (!ok) {
                duckvep_sql_free(&sql);
                free_all(names, 2);
                return;
            }
        } else {
            report(*error, relations == 1 ? duckvep_core_prepare_sv_failed
                                          : duckvep_core_prepare_expansionhunter_failed);
        }
        duckvep_sql_free(&sql);
        free_all(names, 2);
        if (!ok) {
            return;
        }
    }
}

static void prepare_sv_exec(duckdb_v2_scalar_function_exec_info_handle info,
                            duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    prepare_common(info, error, 1);
}
static void prepare_expansionhunter_exec(duckdb_v2_scalar_function_exec_info_handle info,
                                         duckdb_v2_context_handle context,
                                         duckdb_v2_error_info_handle *error) {
    (void)context;
    prepare_common(info, error, 2);
}

/* ---------------------------------------------------------------------------
 * duckvep_prepare_breakend_pairs_sql, _breakend_fusion_sql, _structural_hgvs_sql
 * ------------------------------------------------------------------------- */

static void structural_common(duckdb_v2_scalar_function_exec_info_handle info,
                              duckdb_v2_error_info_handle *error, duckvep_structural_kind_t kind) {
    static const char *const keys[] = {"max_span"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_INTEGER};
    size_t relations = duckvep_core_structural_relations(kind);
    bool has_max_span = duckvep_core_structural_has_max_span(kind);
    call c;
    options option;
    if (!call_open(info, relations, &c, error)) {
        return;
    }
    if (!options_open(c.argc == relations + 1 ? call_argument(info, relations, error) : NULL, keys, kinds,
                      has_max_span ? 1 : 0, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *names[2] = {0};
        int64_t max_span = 5000;
        bool ok = true;
        duckvep_sql_text sql = {0};
        if (c.argc == relations + 1) {
            if (!options_check(&option, row, error)) {
                return;
            }
            if (has_max_span) {
                duckvep_cell_t cell;
                char message[160];
                bool have = option_cell(&option, 0, row, &cell);
                if (!duckvep_core_structural_max_span(have ? &cell : NULL, &max_span, message, sizeof message)) {
                    plain_error(&c, message);
                    return;
                }
            }
        }
        for (size_t i = 0; ok && i < relations; ++i) {
            names[i] = copy_argument(&c, i, row);
            if (!names[i]) {
                ok = false;
            }
        }
        if (ok) {
            ok = duckvep_core_structural_sql(kind, names, max_span, &sql);
        }
        if (ok) {
            ok = write_sql(&c, row, &sql);
            if (!ok) {
                duckvep_sql_free(&sql);
                free_all(names, 2);
                return;
            }
        } else {
            char message[160];
            (void)snprintf(message, sizeof message, "%s: invalid relation name or allocation failure",
                           duckvep_core_structural_name(kind));
            plain_error(&c, message);
        }
        duckvep_sql_free(&sql);
        free_all(names, 2);
        if (!ok) {
            return;
        }
    }
}

static void pairs_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                       duckdb_v2_error_info_handle *error) {
    (void)context;
    structural_common(info, error, DUCKVEP_STRUCTURAL_PAIRS);
}
static void fusion_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                        duckdb_v2_error_info_handle *error) {
    (void)context;
    structural_common(info, error, DUCKVEP_STRUCTURAL_FUSION);
}
static void hgvs_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                      duckdb_v2_error_info_handle *error) {
    (void)context;
    structural_common(info, error, DUCKVEP_STRUCTURAL_HGVS);
}

/* ---------------------------------------------------------------------------
 * duckvep_lof_sql(annotations, transcripts, reference [, options])
 * ------------------------------------------------------------------------- */

static void lof_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                     duckdb_v2_error_info_handle *error) {
    (void)context;
    call c;
    options option;
    if (!call_open(info, 3, &c, error)) {
        return;
    }
    if (!options_open(c.argc == 4 ? call_argument(info, 3, error) : NULL, duckvep_core_lof_option_names,
                      duckvep_core_lof_option_kinds, DUCKVEP_LOF_OPTIONS, &option, error)) {
        return;
    }
    for (idx_t row = 0; row < c.rows; ++row) {
        char *names[3] = {0};
        duckvep_lof_options_t parsed;
        bool ok = true;
        duckvep_sql_text sql = {0};
        memset(&parsed, 0, sizeof parsed);
        for (idx_t i = 0; ok && i < 3; ++i) {
            names[i] = copy_argument(&c, i, row);
            ok = names[i] != NULL;
        }
        if (!ok) {
            report(*error, duckvep_core_lof_names_failed);
        }
        if (ok && c.argc == 4) {
            duckvep_cell_t cells[DUCKVEP_LOF_OPTIONS];
            const duckvep_cell_t *pointers[DUCKVEP_LOF_OPTIONS] = {0};
            const char *message;
            if (!options_check(&option, row, error)) {
                free_all(names, 3);
                return;
            }
            for (size_t i = 0; i < DUCKVEP_LOF_OPTIONS; ++i) {
                pointers[i] = option_or_null(&option, i, row, &cells[i]);
            }
            message = duckvep_core_lof_read_options(pointers, &parsed);
            if (message) {
                report(*error, message);
                ok = false;
            }
        } else if (ok) {
            const duckvep_cell_t *none[DUCKVEP_LOF_OPTIONS] = {0};
            (void)duckvep_core_lof_read_options(none, &parsed);
        }
        if (ok && !duckvep_core_lof_sql(names, &parsed, &sql)) {
            report(*error, duckvep_core_lof_names_failed);
            ok = false;
        }
        if (ok) {
            ok = write_sql(&c, row, &sql);
        }
        duckvep_sql_free(&sql);
        duckvep_core_lof_free_options(&parsed);
        free_all(names, 3);
        if (!ok) {
            return;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Registration: `required` VARCHAR parameters, plus an optional options STRUCT
 * (an ANY parameter), as two overloads.
 * ------------------------------------------------------------------------- */

static bool register_builder(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                             const char *name, idx_t required,
                             duckdb_v2_scalar_function_exec_callback_fn callback,
                             duckdb_v2_error_info_handle *error) {
    static const char *const names[MAX_ARGS] = {"arg1", "arg2", "arg3", "arg4", "arg5",
                                                "arg6", "arg7", "arg8", "options", "options"};
    char const *parameter_types[MAX_ARGS];
    for (idx_t arity = required; arity <= required + 1; ++arity) {
        for (idx_t i = 0; i < arity; ++i) {
            parameter_types[i] = i < required ? "VARCHAR" : "ANY";
        }
        if (!register_scalar(extension, context, name, parameter_types, names, arity, "VARCHAR", callback,
                             error)) {
            return false;
        }
    }
    return true;
}

bool host_v2_register_builders(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                               duckdb_v2_error_info_handle *error) {
    return register_builder(extension, context, "duckvep_ensembl_regions_sql", 3, regions_exec, error) &&
           register_builder(extension, context, "duckvep_ensembl_transcripts_sql", 3, transcripts_exec, error) &&
           register_builder(extension, context, "duckvep_ensembl_regulation_features_sql", 2,
                            regulation_exec, error) &&
           register_builder(extension, context, "duckvep_model_receipt_sql", 8, receipt_exec, error) &&
           register_builder(extension, context, "duckvep_annotate_sql", 2, annotate_exec, error) &&
           register_builder(extension, context, "duckvep_annotate_projected_sql", 2, projected_exec, error) &&
           register_builder(extension, context, "duckvep_transcript_projection_sql", 3, projection_exec, error) &&
           register_builder(extension, context, "duckvep_prepare_sv_geometry_sql", 1, prepare_sv_exec, error) &&
           register_builder(extension, context, "duckvep_prepare_expansionhunter_sql", 2,
                            prepare_expansionhunter_exec, error) &&
           register_builder(extension, context, "duckvep_prepare_breakend_pairs_sql", 1, pairs_exec, error) &&
           register_builder(extension, context, "duckvep_prepare_breakend_fusion_sql", 2, fusion_exec, error) &&
           register_builder(extension, context, "duckvep_prepare_structural_hgvs_sql", 2, hgvs_exec, error) &&
           register_builder(extension, context, "duckvep_lof_sql", 3, lof_exec, error);
}
