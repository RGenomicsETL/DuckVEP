#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_registration.h"
#include "core/duckvep_core_prepare.h"
#include <stdlib.h>
#include <string.h>
DUCKDB_EXTENSION_EXTERN

static void prepare_sv(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector name_vector = duckdb_data_chunk_get_vector(input, 0);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        duckdb_vector option = NULL;
        const char *const names[] = {"unused"};
        const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_TEXT};
        if (argc == 2 && !duckvep_builder_option_vectors(info,
            duckdb_data_chunk_get_vector(input, 1), row, names, kinds, 0, &option)) return;
        uint64_t *valid = duckdb_vector_get_validity(name_vector);
        char *name = valid && !duckdb_validity_row_is_valid(valid, row) ? NULL :
            duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(name_vector))[row]);
        duckvep_sql_text sql = {0};
        bool ok = duckvep_core_prepare_sv_sql(name, &sql);
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckvep_builder_set_error(info, duckvep_core_prepare_sv_failed);
        duckvep_sql_free(&sql);
        duckvep_budget_free(name);
        if (!ok) return;
    }
}

static void prepare_str(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    const char *const names[] = {"unused"};
    const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_TEXT};
    duckdb_vector args[3];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        duckdb_vector option = NULL;
        if (argc == 3 && !duckvep_builder_option_vectors(info, args[2], row, names, kinds, 0, &option)) return;
        char *tables[2] = {0};
        bool ok = true;
        for (idx_t i = 0; i < 2; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) { ok = false; break; }
            tables[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!tables[i]) { ok = false; break; }
        }
        duckvep_sql_text sql = {0};
        if (ok) ok = duckvep_core_prepare_expansionhunter_sql(tables[0], tables[1], &sql);
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckvep_builder_set_error(info, duckvep_core_prepare_expansionhunter_failed);
        duckvep_sql_free(&sql);
        for (size_t i = 0; i < 2; i++) duckvep_budget_free(tables[i]);
        if (!ok) return;
    }
}

bool register_duckvep_prepare_sql(duckdb_connection connection) {
    return duckvep_register_builder(connection, "duckvep_prepare_sv_geometry_sql", 1, prepare_sv) &&
        duckvep_register_builder(connection, "duckvep_prepare_expansionhunter_sql", 2, prepare_str);
}
