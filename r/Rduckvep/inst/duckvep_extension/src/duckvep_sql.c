/* Relation-oriented DuckVEP SQL surface and static annotation metadata. */
#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
DUCKDB_EXTENSION_EXTERN

#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckvep_so.h"
#include "duckvep_codon.h"
#include "duckvep_sql.h"
#include "duckvep_builder.h"
#include "core/duckvep_core_geometry.h"

typedef struct {
	idx_t offset;
} duckvep_so_scan_t;

static void
duckvep_so_terms_bind(duckdb_bind_info info)
{
	duckdb_logical_type utinyint_type, ubigint_type, varchar_type;

	utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_bind_add_result_column(info, "bit_index", utinyint_type);
	duckdb_bind_add_result_column(info, "consequence_mask", ubigint_type);
	duckdb_bind_add_result_column(info, "consequence", varchar_type);
	duckdb_bind_add_result_column(info, "impact_code", utinyint_type);
	duckdb_bind_add_result_column(info, "impact", varchar_type);
	duckdb_bind_add_result_column(info, "severity_rank", utinyint_type);
	duckdb_bind_add_result_column(info, "evaluator_tier", utinyint_type);
	duckdb_destroy_logical_type(&utinyint_type);
	duckdb_destroy_logical_type(&ubigint_type);
	duckdb_destroy_logical_type(&varchar_type);
}

static void
duckvep_so_scan_destroy(void *data)
{
	duckdb_free(data);
}

static void
duckvep_so_terms_init(duckdb_init_info info)
{
	duckvep_so_scan_t *scan;

	scan = duckdb_malloc(sizeof(*scan));
	if (scan == NULL) {
		duckdb_init_set_error(info,
		    "duckvep_so_terms: failed to allocate scan state");
		return;
	}
	scan->offset = 0;
	duckdb_init_set_init_data(info, scan, duckvep_so_scan_destroy);
}

static void
duckvep_so_terms_scan(duckdb_function_info info, duckdb_data_chunk output)
{
	duckvep_so_scan_t *scan;
	duckdb_vector bit_vector, mask_vector, consequence_vector;
	duckdb_vector impact_code_vector, impact_vector, rank_vector, tier_vector;
	uint8_t *bits, *impact_codes, *ranks, *tiers;
	uint64_t *masks;
	idx_t count, vector_size;

	scan = duckdb_function_get_init_data(info);
	if (scan == NULL) {
		duckdb_data_chunk_set_size(output, 0);
		return;
	}
	bit_vector = duckdb_data_chunk_get_vector(output, 0);
	mask_vector = duckdb_data_chunk_get_vector(output, 1);
	consequence_vector = duckdb_data_chunk_get_vector(output, 2);
	impact_code_vector = duckdb_data_chunk_get_vector(output, 3);
	impact_vector = duckdb_data_chunk_get_vector(output, 4);
	rank_vector = duckdb_data_chunk_get_vector(output, 5);
	tier_vector = duckdb_data_chunk_get_vector(output, 6);
	bits = duckdb_vector_get_data(bit_vector);
	masks = duckdb_vector_get_data(mask_vector);
	impact_codes = duckdb_vector_get_data(impact_code_vector);
	ranks = duckdb_vector_get_data(rank_vector);
	tiers = duckdb_vector_get_data(tier_vector);
	vector_size = duckdb_vector_size();
	count = 0;
	while (count < vector_size && scan->offset < duckvep_core_so_term_count()) {
		duckvep_core_so_term_t term;

		(void)duckvep_core_so_term((size_t)scan->offset, &term);
		bits[count] = term.bit_index;
		masks[count] = term.consequence_mask;
		impact_codes[count] = term.impact_code;
		ranks[count] = term.severity_rank;
		tiers[count] = term.evaluator_tier;
		duckdb_vector_assign_string_element(consequence_vector, count,
		    term.consequence);
		duckdb_vector_assign_string_element(impact_vector, count,
		    term.impact);
		count++;
		scan->offset++;
	}
	duckdb_data_chunk_set_size(output, count);
}

static bool
duckvep_register_so_terms(duckdb_connection connection)
{
	duckdb_table_function function;
	duckdb_state state;

	function = duckdb_create_table_function();
	duckdb_table_function_set_name(function, "duckvep_so_terms");
	duckdb_table_function_set_bind(function, duckvep_so_terms_bind);
	duckdb_table_function_set_init(function, duckvep_so_terms_init);
	duckdb_table_function_set_function(function, duckvep_so_terms_scan);
	state = duckdb_register_table_function(connection, function);
	duckdb_destroy_table_function(&function);
	return state == DuckDBSuccess;
}

static bool duckvep_projection_table(duckvep_sql_text *sql, const char *name);

/* Annotation SQL keeps the event relation in the caller's transaction. */
#include "duckvep_annotate_template.h"
#include "duckvep_projected_template.h"

static bool
duckvep_annotate_number(duckdb_vector vector, idx_t row, duckvep_sql_text *sql)
{
    uint64_t *validity = duckdb_vector_get_validity(vector);
    duckdb_logical_type type = duckdb_vector_get_column_type(vector);
    duckdb_type id = duckdb_get_type_id(type);
    duckdb_destroy_logical_type(&type);
    if (id == DUCKDB_TYPE_SQLNULL || (validity && !duckdb_validity_row_is_valid(validity, row)))
        return duckvep_sql_append(sql, "NULL");
    const void *data = duckdb_vector_get_data(vector);
    int64_t signed_value = 0;
    uint64_t unsigned_value = 0;
    bool is_unsigned = false;
    switch (id) {
    case DUCKDB_TYPE_TINYINT: signed_value = ((const int8_t *)data)[row]; break;
    case DUCKDB_TYPE_SMALLINT: signed_value = ((const int16_t *)data)[row]; break;
    case DUCKDB_TYPE_INTEGER: signed_value = ((const int32_t *)data)[row]; break;
    case DUCKDB_TYPE_BIGINT: signed_value = ((const int64_t *)data)[row]; break;
    case DUCKDB_TYPE_UTINYINT: unsigned_value = ((const uint8_t *)data)[row]; is_unsigned = true; break;
    case DUCKDB_TYPE_USMALLINT: unsigned_value = ((const uint16_t *)data)[row]; is_unsigned = true; break;
    case DUCKDB_TYPE_UINTEGER: unsigned_value = ((const uint32_t *)data)[row]; is_unsigned = true; break;
    case DUCKDB_TYPE_UBIGINT: unsigned_value = ((const uint64_t *)data)[row]; is_unsigned = true; break;
    default: return false;
    }
    char buffer[32];
    if (is_unsigned) snprintf(buffer, sizeof(buffer), "%" PRIu64, unsigned_value);
    else snprintf(buffer, sizeof(buffer), "%" PRId64, signed_value);
    return duckvep_sql_append(sql, buffer);
}

static bool
duckvep_annotate_boolean(duckdb_vector vector, idx_t row, duckvep_sql_text *sql)
{
    uint64_t *validity = duckdb_vector_get_validity(vector);
    duckdb_logical_type type = duckdb_vector_get_column_type(vector);
    duckdb_type id = duckdb_get_type_id(type);
    duckdb_destroy_logical_type(&type);
    if (id == DUCKDB_TYPE_SQLNULL || (validity && !duckdb_validity_row_is_valid(validity, row)))
        return duckvep_sql_append(sql, "NULL");
    return duckvep_sql_append(sql, ((const uint8_t *)duckdb_vector_get_data(vector))[row] ? "true" : "false");
}

static void
duckvep_annotate_builder_impl(duckdb_function_info info, duckdb_data_chunk input,
    duckdb_vector output, bool projected)
{
    const char *const keys[] = {"hgvs", "upstream_distance", "downstream_distance", "rich"};
    const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_BOOLEAN, DUCKVEP_OPTION_INTEGER,
        DUCKVEP_OPTION_INTEGER, DUCKVEP_OPTION_BOOLEAN};
    const char *const projected_keys[] = {"upstream_distance", "downstream_distance"};
    const duckvep_option_kind projected_kinds[] = {DUCKVEP_OPTION_INTEGER, DUCKVEP_OPTION_INTEGER};
    const char *const tokens[] = {"__DUCKVEP_EVENTS__", "__DUCKVEP_MODEL__", "__DUCKVEP_HGVS__",
        "__DUCKVEP_UPSTREAM__", "__DUCKVEP_DOWNSTREAM__", "__DUCKVEP_RICH__"};
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[3];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *names[2] = {0};
        bool ok = true;
        for (idx_t i = 0; i < 2; i++) {
            uint64_t *validity = duckdb_vector_get_validity(args[i]);
            if (validity && !duckdb_validity_row_is_valid(validity, row)) { ok = false; break; }
            names[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!names[i]) { ok = false; break; }
        }
        if (ok && projected && names[1][0] == '\0') {
            duckdb_scalar_function_set_error(info, "duckvep_annotate_projected: model_name must be non-empty");
            for (idx_t i = 0; i < 2; i++) duckvep_budget_free(names[i]);
            return;
        }
        duckdb_vector fields[4] = {0};
        if (ok && argc == 3) {
            ok = projected ? duckvep_builder_option_vectors(info, args[2], row,
                projected_keys, projected_kinds, 2, fields) :
                duckvep_builder_option_vectors(info, args[2], row, keys, kinds, 4, fields);
            if (!ok) { for (idx_t i = 0; i < 2; i++) duckvep_budget_free(names[i]); return; }
        }
        duckvep_sql_text values[6] = {{0}};
        if (ok) ok = duckvep_projection_table(&values[0], names[0]) &&
            duckvep_sql_literal(&values[1], names[1]) &&
            duckvep_sql_append(&values[2], "false") &&
            duckvep_sql_append(&values[3], "5000") &&
            duckvep_sql_append(&values[4], "5000") &&
            duckvep_sql_append(&values[5], "false");
        for (size_t i = 0; ok && i < (projected ? 2 : 4); i++) {
            if (!fields[i]) continue;
            size_t dest = projected ? i + 3 : (i == 3 ? 5 : i + 2);
            duckvep_sql_free(&values[dest]);
            ok = !projected && (i == 0 || i == 3) ?
                duckvep_annotate_boolean(fields[i], row, &values[dest]) :
                duckvep_annotate_number(fields[i], row, &values[dest]);
        }
        duckvep_sql_text sql = {0};
        size_t part_count = projected ? sizeof(duckvep_projected_parts) / sizeof(*duckvep_projected_parts) :
            sizeof(duckvep_annotate_parts) / sizeof(*duckvep_annotate_parts);
        for (size_t i = 0; ok && i < part_count; i++) {
            const char *part = projected ? duckvep_projected_parts[i] : duckvep_annotate_parts[i];
            while (ok && *part) {
                const char *mark = strstr(part, "__DUCKVEP_");
                if (!mark) { ok = duckvep_sql_append(&sql, part); break; }
                size_t length = (size_t)(mark - part);
                char *prefix = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
                if (!prefix) { ok = false; break; }
                memcpy(prefix, part, length); prefix[length] = 0;
                ok = duckvep_sql_append(&sql, prefix); duckvep_budget_free(prefix);
                size_t index = 0;
                while (index < 6 && strncmp(mark, tokens[index], strlen(tokens[index])) != 0) index++;
                if (index == 6) { ok = false; break; }
                if (ok) ok = duckvep_sql_append(&sql, values[index].data);
                part = mark + strlen(tokens[index]);
            }
        }
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckvep_builder_set_error(info, "duckvep_annotate_sql: invalid input or allocation failure");
        duckvep_sql_free(&sql);
        for (size_t i = 0; i < 6; i++) duckvep_sql_free(&values[i]);
        for (idx_t i = 0; i < 2; i++) duckvep_budget_free(names[i]);
        if (!ok) return;
    }
}

static void
duckvep_annotate_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output)
{
    duckvep_annotate_builder_impl(info, input, output, false);
}

static void
duckvep_projected_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output)
{
    duckvep_annotate_builder_impl(info, input, output, true);
}

static bool
duckvep_register_annotate_relation(duckhts_registration_t *registration)
{
    return duckvep_register_builder(registration->connection, "duckvep_annotate_sql", 2,
        duckvep_annotate_builder) &&
        duckvep_register_builder(registration->connection, "duckvep_annotate_projected_sql", 2,
        duckvep_projected_builder);
}

/* SQL presentation uses the kernel's immutable genetic-code authority. The
 * string is indexed by T/C/A/G at each of the three codon positions. */
static void
duckvep_projection_code(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	duckdb_vector vector = duckdb_data_chunk_get_vector(input, 0);
	const uint8_t *codes = duckdb_vector_get_data(vector);
	uint64_t *validity = duckdb_vector_get_validity(vector);
	idx_t count = duckdb_data_chunk_get_size(input);

	duckdb_vector_ensure_validity_writable(output);
	for (idx_t row = 0; row < count; row++) {
		const char *amino_acids;

		if (validity != NULL && !duckdb_validity_row_is_valid(validity, row)) {
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
			continue;
		}
		amino_acids = duckvep_codon_table_amino_acids((duckvep_codon_table_t)codes[row]);
		if (amino_acids == NULL) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_transcript_projection: unsupported genetic code");
			return;
		}
		duckdb_vector_assign_string_element_len(output, row, amino_acids, 64);
	}
}

static bool
duckvep_register_projection_code(duckdb_connection connection)
{
	duckdb_scalar_function function = duckdb_create_scalar_function();
	duckdb_logical_type code_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	duckdb_logical_type text_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_state state;

	duckdb_scalar_function_set_name(function, "__duckvep_projection_code");
	duckdb_scalar_function_add_parameter(function, code_type);
	duckdb_scalar_function_set_return_type(function, text_type);
	duckdb_scalar_function_set_function(function, duckvep_projection_code);
	state = duckdb_register_scalar_function(connection, function);
	duckdb_destroy_scalar_function(&function);
	duckdb_destroy_logical_type(&code_type);
	duckdb_destroy_logical_type(&text_type);
	return state == DuckDBSuccess;
}

bool
register_duckvep_sql_kernels(duckhts_registration_t *registration)
{
	if (!duckvep_register_so_terms(registration->connection)) {
		return duckhts_registration_error(registration,
		    "DuckHTS could not register duckvep_so_terms");
	}
	if (!duckvep_register_projection_code(registration->connection)) {
		return duckhts_registration_error(registration,
		    "DuckHTS could not register __duckvep_projection_code");
	}
	if (!duckvep_register_phase_kernels(registration->connection)) {
		return duckhts_registration_error(registration,
		    "DuckHTS could not register DuckVEP phase kernels");
	}
	return true;
}

static bool
duckvep_projection_table(duckvep_sql_text *sql, const char *name)
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

#include "duckvep_projection_template.h"

static void
duckvep_projection_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output)
{
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[4];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *names[3] = {0};
        bool ok = true;
        for (idx_t i = 0; i < 3; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) { ok = false; break; }
            names[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!names[i] || !*names[i]) { ok = false; break; }
        }
        if (ok && argc == 4) {
            duckdb_vector field = NULL;
            const char *const keys[] = {"unused"};
            const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_TEXT};
            ok = duckvep_builder_option_vectors(info, args[3], row, keys, kinds, 0, &field);
            if (!ok) {
                for (idx_t i = 0; i < 3; i++) duckvep_budget_free(names[i]);
                return;
            }
        }
        duckvep_sql_text sql = {0};
        for (size_t i = 0; ok && i < sizeof(duckvep_projection_parts) / sizeof(*duckvep_projection_parts); i++) {
            const char *part = duckvep_projection_parts[i];
            while (ok && *part) {
                const char *mark = strstr(part, "__DUCKVEP_");
                if (!mark) { ok = duckvep_sql_append(&sql, part); break; }
                size_t length = (size_t)(mark - part);
                char *prefix = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
                if (!prefix) { ok = false; break; }
                memcpy(prefix, part, length); prefix[length] = 0;
                ok = duckvep_sql_append(&sql, prefix); duckvep_budget_free(prefix);
                const char *token = NULL; idx_t which = 0;
                if (!strncmp(mark, "__DUCKVEP_EVENTS_TABLE__", strlen("__DUCKVEP_EVENTS_TABLE__"))) token = "__DUCKVEP_EVENTS_TABLE__";
                else if (!strncmp(mark, "__DUCKVEP_ANNOTATIONS_TABLE__", strlen("__DUCKVEP_ANNOTATIONS_TABLE__"))) { token = "__DUCKVEP_ANNOTATIONS_TABLE__"; which = 1; }
                else if (!strncmp(mark, "__DUCKVEP_TRANSCRIPTS_TABLE__", strlen("__DUCKVEP_TRANSCRIPTS_TABLE__"))) { token = "__DUCKVEP_TRANSCRIPTS_TABLE__"; which = 2; }
                else ok = false;
                if (ok) ok = duckvep_projection_table(&sql, names[which]);
                if (token) part = mark + strlen(token);
            }
        }
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckvep_builder_set_error(info, "duckvep_transcript_projection_sql: invalid table name, options, or allocation failure");
        duckvep_sql_free(&sql);
        for (idx_t i = 0; i < 3; i++) duckvep_budget_free(names[i]);
        if (!ok) return;
    }
}

static bool
duckvep_register_projection_relation(duckhts_registration_t *registration)
{
    return duckvep_register_builder(registration->connection,
        "duckvep_transcript_projection_sql", 3, duckvep_projection_builder);
}

bool
register_duckvep_sql_functions(duckhts_registration_t *registration)
{
	return duckvep_register_repeat_alleles(registration->connection) &&
	    duckvep_register_annotate_relation(registration) &&
	    duckvep_register_projection_relation(registration);
}
