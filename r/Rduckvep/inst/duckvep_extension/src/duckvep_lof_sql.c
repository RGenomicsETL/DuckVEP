#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_registration.h"
#include "core/duckvep_core_lof.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
DUCKDB_EXTENSION_EXTERN
#include "duckvep_v1_cells.h"

/*
 * duckvep_lof_sql: the LOFTEE loss-of-function relation as SQL over rows DuckVEP
 * already produces. The builder only names the relations and the constants; the
 * rules are ordinary joins and CASE expressions (see docs/functions.md and
 * benchmarks/duckvep_lof.md). The rules follow konradjk/loftee at a46b502.
 */

static void lof_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[4];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *names[3] = {0};
        duckvep_lof_options_t options;
        memset(&options, 0, sizeof(options));
        bool ok = true;
        for (size_t i = 0; ok && i < 3; i++) {
            uint64_t *validity = duckdb_vector_get_validity(args[i]);
            ok = !validity || duckdb_validity_row_is_valid(validity, row);
            if (ok) names[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            ok = ok && names[i] != NULL;
        }
        if (!ok) duckvep_builder_set_error(info, duckvep_core_lof_names_failed);
        if (ok && argc == 4) {
            duckdb_vector fields[DUCKVEP_LOF_OPTIONS] = {0};
            duckvep_cell_t cells[DUCKVEP_LOF_OPTIONS];
            const duckvep_cell_t *pointers[DUCKVEP_LOF_OPTIONS] = {0};
            const char *message;
            ok = duckvep_builder_option_vectors(info, args[3], row, duckvep_core_lof_option_names,
                (const duckvep_option_kind[]){DUCKVEP_OPTION_TEXT, DUCKVEP_OPTION_TEXT, DUCKVEP_OPTION_TEXT,
                    DUCKVEP_OPTION_INTEGER, DUCKVEP_OPTION_NUMERIC, DUCKVEP_OPTION_BOOLEAN},
                DUCKVEP_LOF_OPTIONS, fields);
            for (size_t i = 0; ok && i < DUCKVEP_LOF_OPTIONS; i++) {
                duckdb_logical_type type;
                if (!fields[i]) continue;
                type = duckdb_vector_get_column_type(fields[i]);
                duckvep_v1_fill_cell(fields[i], duckdb_get_type_id(type), 0, row, &cells[i]);
                duckdb_destroy_logical_type(&type);
                pointers[i] = &cells[i];
            }
            if (ok) {
                message = duckvep_core_lof_read_options(pointers, &options);
                if (message) {
                    duckvep_builder_set_error(info, message);
                    ok = false;
                }
            }
        } else if (ok) {
            (void)duckvep_core_lof_read_options((const duckvep_cell_t *const[DUCKVEP_LOF_OPTIONS]){0}, &options);
        }
        duckvep_sql_text sql = {0};
        if (ok && !duckvep_core_lof_sql(names, &options, &sql)) {
            duckvep_builder_set_error(info, duckvep_core_lof_names_failed);
            ok = false;
        }
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        duckvep_sql_free(&sql);
        duckvep_core_lof_free_options(&options);
        for (size_t i = 0; i < 3; i++) duckvep_budget_free(names[i]);
        if (!ok) return;
    }
}

bool register_duckvep_lof_sql(duckdb_connection connection) {
    return duckvep_register_builder(connection, "duckvep_lof_sql", 3, lof_builder);
}
