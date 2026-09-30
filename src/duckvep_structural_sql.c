#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_registration.h"
#include "core/duckvep_core_prepare.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
DUCKDB_EXTENSION_EXTERN
#include "duckvep_v1_cells.h"

static void run_builder(duckvep_structural_kind_t kind, duckdb_function_info info,
                        duckdb_data_chunk input, duckdb_vector output) {
    size_t relation_count = duckvep_core_structural_relations(kind);
    bool has_max_span = duckvep_core_structural_has_max_span(kind);
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[3];
    for (idx_t i = 0; i < argc && i < 3; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        const char *const names[] = {"max_span"};
        const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_INTEGER};
        duckdb_vector option_fields[1] = {NULL};
        int64_t max_span = 5000;
        char *relations[2] = {0};
        bool ok = true;
        if (argc == relation_count + 1) {
            ok = duckvep_builder_option_vectors(info, args[relation_count], row, names, kinds,
                                                has_max_span ? 1u : 0u, option_fields);
            if (ok && has_max_span) {
                duckvep_cell_t cell;
                char message[160];
                if (option_fields[0]) {
                    duckdb_logical_type type = duckdb_vector_get_column_type(option_fields[0]);
                    duckvep_v1_fill_cell(option_fields[0], duckdb_get_type_id(type), 0, row, &cell);
                    duckdb_destroy_logical_type(&type);
                }
                ok = duckvep_core_structural_max_span(option_fields[0] ? &cell : NULL, &max_span,
                                                      message, sizeof(message));
                if (!ok) duckdb_scalar_function_set_error(info, message);
            }
            if (!ok) return;
        }
        for (size_t i = 0; ok && i < relation_count; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) { ok = false; break; }
            relations[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!relations[i]) ok = false;
        }
        duckvep_sql_text sql = {0};
        if (ok) ok = duckvep_core_structural_sql(kind, relations, max_span, &sql);
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else {
            char message[160];
            snprintf(message, sizeof(message),
                     "%s: invalid relation name or allocation failure", duckvep_core_structural_name(kind));
            duckdb_scalar_function_set_error(info, message);
        }
        duckvep_sql_free(&sql);
        for (size_t i = 0; i < 2; i++) duckvep_budget_free(relations[i]);
        if (!ok) return;
    }
}

static void prepare_pairs(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    run_builder(DUCKVEP_STRUCTURAL_PAIRS, info, input, output);
}
static void prepare_fusion(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    run_builder(DUCKVEP_STRUCTURAL_FUSION, info, input, output);
}
static void prepare_hgvs(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    run_builder(DUCKVEP_STRUCTURAL_HGVS, info, input, output);
}

bool register_duckvep_structural_sql(duckdb_connection connection) {
    return duckvep_register_builder(connection, "duckvep_prepare_breakend_pairs_sql", 1, prepare_pairs) &&
        duckvep_register_builder(connection, "duckvep_prepare_breakend_fusion_sql", 2, prepare_fusion) &&
        duckvep_register_builder(connection, "duckvep_prepare_structural_hgvs_sql", 2, prepare_hgvs);
}
