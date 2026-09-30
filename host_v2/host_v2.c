/* DuckDB C API v2 host for DuckVEP (preview; slice 1 of issue #8).
 *
 * This file is a v2-only implementation: it includes duckdb_extension_v2.h,
 * calls only duckdb_v2_* functions, compiles with the unstable and deprecated
 * surfaces disabled, and runs no SQL at LOAD. The row logic lives in
 * src/core/ and src/kernel/, which do not mention DuckDB.
 *
 * Ported so far: duckvep_so_terms (table), duckvep_allele_geometry and
 * duckvep_breakend_geometry (scalars); host_v2_nested.c adds the nested family. See docs/v2-host.md. */
#include "host_v2_common.h"

#include "core/duckvep_core_geometry.h"

bool host_v2_register_builders(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                               duckdb_v2_error_info_handle *error);
bool host_v2_register_model(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                            duckdb_v2_error_info_handle *error);
bool host_v2_register_nested(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error);

/* ---------------------------------------------------------------------------
 * duckvep_so_terms(): table function, resumable scan state.
 * ------------------------------------------------------------------------- */

/* Deliberately below the engine's vector size so the 41 rows take several
 * exec calls and exercise the scan state; the tests pin this. */
#define SO_CHUNK_ROWS 16

enum { SO_BIT, SO_MASK, SO_CONSEQUENCE, SO_IMPACT_CODE, SO_IMPACT, SO_RANK, SO_TIER, SO_COLUMNS };

static const char *const so_names[SO_COLUMNS] = {
    "bit_index", "consequence_mask", "consequence", "impact_code", "impact",
    "severity_rank", "evaluator_tier"};
static const char *const so_types[SO_COLUMNS] = {
    "UTINYINT", "UBIGINT", "VARCHAR", "UTINYINT", "VARCHAR", "UTINYINT", "UTINYINT"};

typedef struct {
    size_t offset;
} so_scan;

static void so_bind(duckdb_v2_table_function_bind_info_handle info,
                    duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    for (idx_t column = 0; column < SO_COLUMNS; ++column) {
        DUCKDB_CALL(make_type(context, so_types[column], &type, &detail));
        DUCKDB_CALL(duckdb_v2_table_function_bind_add_result_column(
            info, string_view(so_names[column]), type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
cleanup:
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void so_scan_destroy(void *pointer) {
    free(pointer);
}

static void so_init(duckdb_v2_table_function_init_global_info_handle info,
                    duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    so_scan *scan = (so_scan *)calloc(1, sizeof(*scan));
    bool owned = false;
    if (!scan) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY,
                  "duckvep_so_terms: failed to allocate scan state");
        return;
    }
    duckdb_v2_opaque state = {scan, so_scan_destroy, NULL};
    /* The scan offset is unsynchronized global state: one thread only. */
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_max_threads(info, 1, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_init_global_set_global_state(info, &state, &detail));
    owned = true;
cleanup:
    if (!owned) {
        so_scan_destroy(scan);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void so_exec(duckdb_v2_table_function_exec_info_handle info,
                    duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    so_scan *scan = NULL;
    duckdb_v2_data_chunk_handle chunk = NULL;
    duckdb_v2_vector_handle columns[SO_COLUMNS] = {0};
    void *data[SO_COLUMNS] = {0};
    duckdb_v2_arena_handle consequence_arena = NULL, impact_arena = NULL;
    idx_t count = 0;
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_global_state(info, (void **)&scan, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_exec_get_output_chunk(info, &chunk, &detail));
    if (!scan) {
        INPUT_ERROR("duckvep_so_terms: missing scan state");
    }
    for (idx_t column = 0; column < SO_COLUMNS; ++column) {
        DUCKDB_CALL(duckdb_v2_data_chunk_get_vector(chunk, column, &columns[column], &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(columns[column], &data[column], &detail));
    }
    DUCKDB_CALL(duckdb_v2_vector_get_arena(columns[SO_CONSEQUENCE], &consequence_arena, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(columns[SO_IMPACT], &impact_arena, &detail));
    while (count < SO_CHUNK_ROWS && scan->offset < duckvep_core_so_term_count()) {
        duckvep_core_so_term_t term;
        if (!duckvep_core_so_term(scan->offset, &term)) {
            INPUT_ERROR("duckvep_so_terms: term index out of range");
        }
        ((uint8_t *)data[SO_BIT])[count] = term.bit_index;
        ((uint64_t *)data[SO_MASK])[count] = term.consequence_mask;
        ((uint8_t *)data[SO_IMPACT_CODE])[count] = term.impact_code;
        ((uint8_t *)data[SO_RANK])[count] = term.severity_rank;
        ((uint8_t *)data[SO_TIER])[count] = term.evaluator_tier;
        DUCKDB_CALL(write_string(consequence_arena, &((duckdb_v2_bytes *)data[SO_CONSEQUENCE])[count],
                                 term.consequence, strlen(term.consequence), &detail));
        DUCKDB_CALL(write_string(impact_arena, &((duckdb_v2_bytes *)data[SO_IMPACT])[count],
                                 term.impact, strlen(term.impact), &detail));
        count++;
        scan->offset++;
    }
    /* An empty batch ends the scan. */
    DUCKDB_CALL(duckdb_v2_vector_set_size(columns[0], count, &detail));
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

static bool register_so_terms(duckdb_v2_extension_handle extension,
                              duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_table_function_handle function = NULL;
    bool success = false;
    duckdb_v2_str name = string_view("duckvep_so_terms");
    DUCKDB_CALL(duckdb_v2_table_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_name(function, &name, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_bind_callback(function, so_bind, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_init_global_callback(function, so_init, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_set_exec_callback(function, so_exec, &detail));
    DUCKDB_CALL(duckdb_v2_table_function_register(function, &detail));
    success = true;
cleanup:
    (void)duckdb_v2_table_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

/* ---------------------------------------------------------------------------
 * duckvep_allele_geometry(position UBIGINT, ref VARCHAR, alt VARCHAR) -> STRUCT
 * ------------------------------------------------------------------------- */

static void allele_geometry_exec(duckdb_v2_scalar_function_exec_info_handle info,
                                 duckdb_v2_context_handle context,
                                 duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[3];
    duckdb_v2_vector_handle result = NULL;
    duckdb_v2_vector_handle fields[DUCKVEP_CORE_GEOMETRY_FIELD_COUNT] = {0};
    void *data[DUCKVEP_CORE_GEOMETRY_FIELD_COUNT] = {0};
    uint64_t *validity[DUCKVEP_CORE_GEOMETRY_FIELD_COUNT] = {0};
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 3, views, error)) {
        goto cleanup;
    }
    if (!open_struct_result(info, rows, DUCKVEP_CORE_GEOMETRY_FIELD_COUNT, &result, fields, data,
                            validity, error)) {
        goto cleanup;
    }
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_core_allele_geometry_t g;
        duckdb_v2_str reference, alternate;
        if (!row_is_valid(&views[0], row) || !row_is_valid(&views[1], row) ||
            !row_is_valid(&views[2], row)) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(result, row, &detail));
            continue;
        }
        reference = string_at(&views[1], row);
        alternate = string_at(&views[2], row);
        if (!duckvep_core_allele_geometry(integer_at(&views[0], row), reference.ptr, reference.len,
                                          alternate.ptr, alternate.len, &g)) {
            INPUT_ERROR(duckvep_core_allele_geometry_error);
        }
        ((uint8_t *)data[0])[row] = g.kind_code;
        ((bool *)data[1])[row] = g.interbase;
        ((uint8_t *)data[2])[row] = g.anchor_side_code;
        ((uint64_t *)data[3])[row] = g.raw_start0;
        ((uint64_t *)data[4])[row] = g.raw_end0;
        ((uint64_t *)data[5])[row] = g.feature_start0;
        ((uint64_t *)data[6])[row] = g.feature_end0;
        ((uint64_t *)data[7])[row] = g.edit_start0;
        ((uint64_t *)data[8])[row] = g.edit_end0;
        if (g.has_insertion_boundary0) {
            ((uint64_t *)data[9])[row] = g.insertion_boundary0;
        } else {
            mark_null(validity[9], row);
        }
        ((uint16_t *)data[10])[row] = g.reference_difference_offset;
        ((uint16_t *)data[11])[row] = g.reference_difference_length;
        ((uint16_t *)data[12])[row] = g.alternate_difference_offset;
        ((uint16_t *)data[13])[row] = g.alternate_difference_length;
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_breakend_geometry(alt VARCHAR) -> STRUCT
 * ------------------------------------------------------------------------- */

enum { BND_MATE_CHROM, BND_MATE_POSITION, BND_LOCAL_JOIN_AFTER, BND_MATE_EXTENDS_RIGHT,
       BND_REPLACEMENT, BND_FIELDS };

static void breakend_geometry_exec(duckdb_v2_scalar_function_exec_info_handle info,
                                   duckdb_v2_context_handle context,
                                   duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[1];
    duckdb_v2_vector_handle result = NULL;
    duckdb_v2_vector_handle fields[BND_FIELDS] = {0};
    void *data[BND_FIELDS] = {0};
    uint64_t *validity[BND_FIELDS] = {0};
    duckdb_v2_arena_handle chrom_arena = NULL, replacement_arena = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 1, views, error)) {
        goto cleanup;
    }
    if (!open_struct_result(info, rows, BND_FIELDS, &result, fields, data, validity, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_vector_get_arena(fields[BND_MATE_CHROM], &chrom_arena, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(fields[BND_REPLACEMENT], &replacement_arena, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_breakend_t parsed;
        duckvep_breakend_status_t status;
        duckdb_v2_str alternate;
        if (!row_is_valid(&views[0], row)) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(result, row, &detail));
            continue;
        }
        alternate = string_at(&views[0], row);
        status = duckvep_breakend_parse((const uint8_t *)alternate.ptr, alternate.len, &parsed);
        if (status == DUCKVEP_BREAKEND_NOT_BREAKEND) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(result, row, &detail));
            continue;
        }
        if (status != DUCKVEP_BREAKEND_OK) {
            INPUT_ERROR(duckvep_core_breakend_error(status));
        }
        if (parsed.has_mate) {
            DUCKDB_CALL(write_string(chrom_arena, &((duckdb_v2_bytes *)data[BND_MATE_CHROM])[row],
                                     (const char *)parsed.mate_chrom, parsed.mate_chrom_length,
                                     &detail));
            ((uint64_t *)data[BND_MATE_POSITION])[row] = parsed.mate_position;
            ((bool *)data[BND_MATE_EXTENDS_RIGHT])[row] = parsed.mate_extends_right != 0u;
        } else {
            mark_null(validity[BND_MATE_CHROM], row);
            mark_null(validity[BND_MATE_POSITION], row);
            mark_null(validity[BND_MATE_EXTENDS_RIGHT], row);
        }
        ((bool *)data[BND_LOCAL_JOIN_AFTER])[row] = parsed.local_join_after != 0u;
        DUCKDB_CALL(write_string(replacement_arena,
                                 &((duckdb_v2_bytes *)data[BND_REPLACEMENT])[row],
                                 (const char *)parsed.replacement, parsed.replacement_length,
                                 &detail));
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * Scalar registration
 * ------------------------------------------------------------------------- */

static bool register_allele_geometry(duckdb_v2_extension_handle extension,
                                     duckdb_v2_context_handle context,
                                     duckdb_v2_error_info_handle *error) {
    static const char *const parameter_types[] = {"UBIGINT", "VARCHAR", "VARCHAR"};
    static const char *const parameter_names[] = {"position", "reference", "alternate"};
    return register_scalar(
        extension, context, "duckvep_allele_geometry", parameter_types, parameter_names, 3,
        "STRUCT(kind_code UTINYINT, interbase BOOLEAN, anchor_side_code UTINYINT, "
        "raw_start0 UBIGINT, raw_end0 UBIGINT, feature_start0 UBIGINT, feature_end0 UBIGINT, "
        "edit_start0 UBIGINT, edit_end0 UBIGINT, insertion_boundary0 UBIGINT, "
        "reference_difference_offset USMALLINT, reference_difference_length USMALLINT, "
        "alternate_difference_offset USMALLINT, alternate_difference_length USMALLINT)",
        allele_geometry_exec, error);
}

static bool register_breakend_geometry(duckdb_v2_extension_handle extension,
                                       duckdb_v2_context_handle context,
                                       duckdb_v2_error_info_handle *error) {
    static const char *const parameter_types[] = {"VARCHAR"};
    static const char *const parameter_names[] = {"alternate"};
    return register_scalar(
        extension, context, "duckvep_breakend_geometry", parameter_types, parameter_names, 1,
        "STRUCT(mate_chrom VARCHAR, mate_position UBIGINT, local_join_after BOOLEAN, "
        "mate_extends_right BOOLEAN, replacement_sequence VARCHAR)",
        breakend_geometry_exec, error);
}

/* LOAD registers native functions in the running database and runs no SQL:
 * nothing is written to the catalog, so a read-only primary database works. */
DUCKDB_EXTENSION_ENTRYPOINT(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                            duckdb_v2_error_info_handle *error) {
    if (!register_so_terms(extension, error)) {
        return;
    }
    if (!register_allele_geometry(extension, context, error)) {
        return;
    }
    if (!register_breakend_geometry(extension, context, error)) {
        return;
    }
    if (!host_v2_register_nested(extension, context, error)) {
        return;
    }
    if (!host_v2_register_builders(extension, context, error)) {
        return;
    }
    if (!host_v2_register_model(extension, context, error)) {
        return;
    }
}
