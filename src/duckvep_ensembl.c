#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
DUCKDB_EXTENSION_EXTERN

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckvep_sql.h"
#include "duckvep_builder.h"
#include "duckvep_v1_cells.h"
#include "core/duckvep_core_ensembl.h"

static void ensembl_builder(duckdb_function_info info, duckdb_data_chunk input,
                            duckdb_vector output, duckvep_ensembl_kind_t kind) {
    static const char *const keys[] = {"species_id"};
    static const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_INTEGER};
    idx_t required = duckvep_core_ensembl_required(kind);
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[4];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *values[3] = {0};
        for (idx_t i = 0; i < required; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) break;
            values[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!values[i]) break;
        }
        bool ok = true;
        for (idx_t i = 0; i < required; i++) if (!values[i] || !*values[i]) ok = false;
        if (!ok) duckdb_scalar_function_set_error(info, duckvep_core_ensembl_names_required);
        duckdb_vector species_field = NULL;
        if (ok && argc > required) {
            if (required == 3) ok = duckvep_builder_option_vectors(info, args[required], row, keys, kinds, 1, &species_field);
            else {
                duckdb_scalar_function_set_error(info, duckvep_core_ensembl_no_options);
                ok = false;
            }
        }
        char species[40] = "1";
        if (ok && species_field) {
            duckdb_logical_type type = duckdb_vector_get_column_type(species_field);
            duckvep_cell_t cell;
            duckvep_v1_fill_cell(species_field, duckdb_get_type_id(type), 0, row, &cell);
            duckdb_destroy_logical_type(&type);
            ok = duckvep_core_ensembl_species(&cell, species);
            if (!ok) duckdb_scalar_function_set_error(info, duckvep_core_ensembl_species_range);
        }
        duckvep_sql_text sql = {0};
        bool building = ok;
        if (ok) ok = duckvep_core_ensembl_sql(kind, (const char *const *)values, species, &sql);
        if (ok) {
            duckdb_vector_ensure_validity_writable(output);
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
            duckdb_vector_assign_string_element(output, row, sql.data);
        } else if (building) duckdb_scalar_function_set_error(info, duckvep_core_ensembl_failed(kind));
        duckvep_sql_free(&sql);
        for (idx_t i = 0; i < required; i++) duckvep_budget_free(values[i]);
        if (!ok) return;
    }
}

static void regions_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    ensembl_builder(info, input, output, DUCKVEP_ENSEMBL_REGIONS);
}
static void transcripts_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    ensembl_builder(info, input, output, DUCKVEP_ENSEMBL_TRANSCRIPTS);
}
static void regulation_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    ensembl_builder(info, input, output, DUCKVEP_ENSEMBL_REGULATION);
}

static void duckvep_model_receipt_sql(duckdb_function_info info,
                                     duckdb_data_chunk input, duckdb_vector output) {
    static const char *const options[] = {"regulation_features_table"};
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[9];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *values[8] = {0};
        for (idx_t i = 0; i < 8; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) continue;
            duckdb_string_t *strings = duckdb_vector_get_data(args[i]);
            values[i] = duckvep_builder_string(strings[row]);
            if (!values[i]) {
                duckvep_builder_set_error(info, "duckvep_model_receipt_sql: invalid argument string or allocation failure");
                for (idx_t j = 0; j <= i; j++) duckvep_budget_free(values[j]);
                return;
            }
        }
        char *option_values[1] = {0};
        bool options_ok = argc != 9 || duckvep_builder_options(info, args[8], row, options, 1, option_values);
        if (options_ok && (!values[0] || !values[1]))
            duckdb_scalar_function_set_error(info, duckvep_core_receipt_tables_required);
        if (!options_ok || !values[0] || !values[1]) {
            for (idx_t i = 0; i < 8; i++) duckvep_budget_free(values[i]);
            duckvep_budget_free(option_values[0]);
            return;
        }
        duckvep_sql_text sql = {0};
        bool ok = duckvep_core_receipt_sql((const char *const *)values, option_values[0], &sql);
        if (!ok) duckvep_builder_set_error(info, "duckvep_model_receipt_sql: allocation failed");
        else {
            duckdb_vector_ensure_validity_writable(output);
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
            duckdb_vector_assign_string_element(output, row, sql.data);
        }
        duckvep_sql_free(&sql);
        for (idx_t i = 0; i < 8; i++) duckvep_budget_free(values[i]);
        duckvep_budget_free(option_values[0]);
        if (!ok) return;
    }
}

bool
register_duckvep_ensembl_functions(duckhts_registration_t *registration)
{
	return duckvep_register_builder(registration->connection, "duckvep_ensembl_regions_sql", 3, regions_builder) &&
	    duckvep_register_builder(registration->connection, "duckvep_ensembl_transcripts_sql", 3, transcripts_builder) &&
	    duckvep_register_builder(registration->connection, "duckvep_ensembl_regulation_features_sql", 2, regulation_builder) &&
	    duckvep_register_builder(registration->connection, "duckvep_model_receipt_sql", 8,
	        duckvep_model_receipt_sql);
}
