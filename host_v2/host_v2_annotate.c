/* The annotation natives of the v2 host: _duckvep_annotate_{small,structural,breakend}_*,
 * _duckvep_annotate_small_{hgvs,projected,...} and __duckvep_projection_code. The work is
 * src/core/duckvep_core_annotate_run.c, compiled unchanged for both hosts over the vector layer
 * of duckvep_host.h; this file opens the call's vectors, runs the core and registers the functions. */
#include "duckvep_host.h"

#include "host_v2_model.h"
#include "core/duckvep_core_annotate_run.h"
#include "core/duckvep_core_annotate_types.h"
#include "kernel/src/duckvep_codon.h"

typedef void (*native_fn)(v2_call *info, v2_call *input, v2_vec *output);

static void run_native(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_error_info_handle *error,
                       native_fn function) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle result = NULL;
    v2_call *call = calloc(1, sizeof(*call));
    uint32_t argc = 0;
    idx_t rows = 0;
    if (!call) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_annotate: out of memory");
        return;
    }
    call->info = info;
    call->error = error;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, &call->user_data, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    call->rows = rows;
    call->argc = argc;
    if (argc > V2_CALL_INPUTS) {
        INPUT_ERROR("duckvep_annotate: too many arguments");
    }
    for (uint32_t i = 0; i < argc; ++i) {
        duckdb_v2_vector_handle vector = NULL;
        duckdb_v2_vector_view view;
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &vector, &detail));
        DUCKDB_CALL(duckdb_v2_vector_flatten(vector, &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_view(vector, &view, &detail));
        call->inputs[i].handle = vector;
        call->inputs[i].data = (void *)view.data;
        call->inputs[i].validity = (uint64_t *)view.validity;
        call->inputs[i].size = rows;
        call->inputs[i].call = call;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &result, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(result, &detail));
    if (!v2_open_writable(call, &call->output, result, rows, true)) {
        goto cleanup;
    }
    function(call, call, &call->output);
    v2_finish(call);
cleanup:
    free(call);
    (void)duckdb_v2_error_info_destroy(&detail);
}

#define NATIVE(name, function)                                                                     \
    static void name(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context, \
                     duckdb_v2_error_info_handle *error) {                                          \
        (void)context;                                                                              \
        run_native(info, error, function);                                                          \
    }

NATIVE(exec_small_rich, duckvep_annotate_scalar)
NATIVE(exec_small_compact, duckvep_annotate_compact_scalar)
NATIVE(exec_small_hgvs, duckvep_annotate_hgvs_scalar)
NATIVE(exec_small_rich_hgvs, duckvep_annotate_rich_hgvs_scalar)
NATIVE(exec_structural_rich, duckvep_annotate_sv_scalar)
NATIVE(exec_structural_compact, duckvep_annotate_sv_compact_scalar)
NATIVE(exec_breakend_rich, duckvep_annotate_breakend_scalar)
NATIVE(exec_breakend_compact, duckvep_annotate_breakend_compact_scalar)
NATIVE(exec_small_projected, duckvep_annotate_projected_scalar)
NATIVE(exec_small_projected_hgvs, duckvep_annotate_projected_hgvs_scalar)

/* __duckvep_projection_code(UTINYINT) -> VARCHAR: the 64-character amino-acid string of a codon table. */
static void projection_code_exec(duckdb_v2_scalar_function_exec_info_handle info,
                                 duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[1];
    duckdb_v2_vector_handle output = NULL;
    duckdb_v2_arena_handle arena = NULL;
    void *data = NULL;
    uint64_t *validity = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 1, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &validity, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(output, &arena, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        const char *amino_acids;
        mark_valid(validity, row);
        if (!row_is_valid(&views[0], row)) {
            mark_null(validity, row);
            continue;
        }
        amino_acids = duckvep_codon_table_amino_acids(
            (duckvep_codon_table_t)((const uint8_t *)views[0].data)[physical_row(&views[0], row)]);
        if (amino_acids == NULL) {
            INPUT_ERROR("duckvep_transcript_projection: unsupported genetic code");
        }
        DUCKDB_CALL(write_string(arena, &((duckdb_v2_bytes *)data)[row], amino_acids, 64, &detail));
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * Registration
 * ------------------------------------------------------------------------- */

static bool result_type_text(char *out, size_t size, duckvep_result_kind_t kind, int hgvs, int projection) {
    duckvep_result_column_t columns[DUCKVEP_RESULT_COLUMNS_MAX];
    size_t count = duckvep_core_result_columns(kind, hgvs, projection, columns);
    size_t at = (size_t)snprintf(out, size, "STRUCT(");
    for (size_t i = 0; i < count; ++i) {
        int written = snprintf(out + at, size - at, "%s%s %s", i ? ", " : "", columns[i].name,
                               duckvep_core_column_type_name(columns[i].type));
        if (written < 0 || (size_t)written >= size - at) {
            return false;
        }
        at += (size_t)written;
    }
    return snprintf(out + at, size - at, ")[]") > 0;
}

typedef enum { FAMILY_SMALL, FAMILY_STRUCTURAL, FAMILY_BREAKEND } family_t;

static bool register_native(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                            model_state *state, const char *name, family_t family, const char *result,
                            duckdb_v2_scalar_function_exec_callback_fn callback,
                            duckdb_v2_error_info_handle *error) {
    static const char *const small[] = {"VARCHAR", "UINTEGER", "UBIGINT", "VARCHAR", "VARCHAR", "UBIGINT", "UBIGINT"};
    static const char *const structural[] = {"VARCHAR", "UINTEGER", "UBIGINT", "UBIGINT", "VARCHAR",
                                             "VARCHAR", "UBIGINT", "UBIGINT"};
    static const char *const breakend[] = {"VARCHAR", "UINTEGER", "UBIGINT", "UINTEGER", "UBIGINT",
                                           "UBIGINT", "UBIGINT"};
    static const char *const names[] = {"model", "seq_region", "position", "p4", "p5", "p6", "p7", "p8"};
    const char *const *types = family == FAMILY_STRUCTURAL ? structural : family == FAMILY_BREAKEND ? breakend
                                                                                                  : small;
    idx_t base = family == FAMILY_STRUCTURAL ? 6 : 5;
    for (idx_t distance = 0; distance <= 2; ++distance) {
        duckdb_v2_error_info_handle detail = NULL;
        duckdb_v2_scalar_function_handle function = NULL;
        duckdb_v2_function_signature_handle signature = NULL;
        duckdb_v2_logical_type_handle type = NULL;
        duckdb_v2_str text = string_view(name);
        bool retained = false, ok = false;
        DUCKDB_CALL(duckdb_v2_scalar_function_create_with_extension(extension, &function, &detail));
        DUCKDB_CALL(duckdb_v2_scalar_function_set_name(function, &text, &detail));
        DUCKDB_CALL(duckdb_v2_scalar_function_get_signature(function, &signature, &detail));
        for (idx_t i = 0; i < base + distance; ++i) {
            DUCKDB_CALL(make_type(context, types[i], &type, &detail));
            DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view(names[i]), type, NULL,
                                                                   &detail));
            DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
        }
        DUCKDB_CALL(make_type(context, result, &type, &detail));
        DUCKDB_CALL(duckdb_v2_function_signature_set_return_type(signature, type, &detail));
        DUCKDB_CALL(duckdb_v2_scalar_function_set_property(function, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY,
                                                           DUCKDB_V2_FUNCTION_PROPERTY_STABILITY_VOLATILE,
                                                           &detail));
        {
            duckdb_v2_opaque data = {state, host_v2_state_release, NULL};
            host_v2_state_retain(state);
            retained = true;
            DUCKDB_CALL(duckdb_v2_scalar_function_set_user_data(function, &data, &detail));
            retained = false;
        }
        DUCKDB_CALL(duckdb_v2_scalar_function_set_exec_callback(function, callback, &detail));
        DUCKDB_CALL(duckdb_v2_scalar_function_register(function, &detail));
        ok = true;
cleanup:
        if (retained) {
            host_v2_state_release(state);
        }
        (void)duckdb_v2_logical_type_destroy(&type);
        (void)duckdb_v2_scalar_function_destroy(&function);
        (void)duckdb_v2_error_info_destroy(&detail);
        if (!ok) {
            return false;
        }
    }
    return true;
}

bool host_v2_register_annotate(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                               model_state *state, duckdb_v2_error_info_handle *error) {
    char rich[2048], rich_hgvs[2048], projected[2048], projected_hgvs[2048], compact[1024], compact_hgvs[1024];
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_scalar_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_str name = string_view("__duckvep_projection_code");
    bool ok = false;
    if (!result_type_text(rich, sizeof rich, DUCKVEP_RESULT_RICH, 0, 0) ||
        !result_type_text(rich_hgvs, sizeof rich_hgvs, DUCKVEP_RESULT_RICH, 1, 0) ||
        !result_type_text(projected, sizeof projected, DUCKVEP_RESULT_RICH, 0, 1) ||
        !result_type_text(projected_hgvs, sizeof projected_hgvs, DUCKVEP_RESULT_RICH, 1, 1) ||
        !result_type_text(compact, sizeof compact, DUCKVEP_RESULT_COMPACT, 0, 0) ||
        !result_type_text(compact_hgvs, sizeof compact_hgvs, DUCKVEP_RESULT_COMPACT_HGVS, 0, 0)) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, "duckvep: result type text too long");
        return false;
    }
    if (!(register_native(extension, context, state, "_duckvep_annotate_small_rich", FAMILY_SMALL, rich,
                          exec_small_rich, error) &&
          register_native(extension, context, state, "_duckvep_annotate_small_compact", FAMILY_SMALL, compact,
                          exec_small_compact, error) &&
          register_native(extension, context, state, "_duckvep_annotate_small_hgvs", FAMILY_SMALL, compact_hgvs,
                          exec_small_hgvs, error) &&
          register_native(extension, context, state, "_duckvep_annotate_small_rich_hgvs", FAMILY_SMALL,
                          rich_hgvs, exec_small_rich_hgvs, error) &&
          register_native(extension, context, state, "_duckvep_annotate_small_projected", FAMILY_SMALL,
                          projected, exec_small_projected, error) &&
          register_native(extension, context, state, "_duckvep_annotate_small_projected_hgvs", FAMILY_SMALL,
                          projected_hgvs, exec_small_projected_hgvs, error) &&
          register_native(extension, context, state, "_duckvep_annotate_structural_rich", FAMILY_STRUCTURAL,
                          rich, exec_structural_rich, error) &&
          register_native(extension, context, state, "_duckvep_annotate_structural_compact", FAMILY_STRUCTURAL,
                          compact, exec_structural_compact, error) &&
          register_native(extension, context, state, "_duckvep_annotate_breakend_rich", FAMILY_BREAKEND, rich,
                          exec_breakend_rich, error) &&
          register_native(extension, context, state, "_duckvep_annotate_breakend_compact", FAMILY_BREAKEND,
                          compact, exec_breakend_compact, error))) {
        return false;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_name(function, &name, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_get_signature(function, &signature, &detail));
    DUCKDB_CALL(make_type(context, "UTINYINT", &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view("code"), type, NULL, &detail));
    DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    DUCKDB_CALL(make_type(context, "VARCHAR", &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_set_return_type(signature, type, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(function, DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING,
                                                       DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING_SPECIAL,
                                                       &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_exec_callback(function, projection_code_exec, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_register(function, &detail));
    ok = true;
cleanup:
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_scalar_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return ok;
}
