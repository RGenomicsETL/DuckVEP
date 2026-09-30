/* Stable-ABI DuckDB vector adapter for the resident DuckVEP kernel: the geometry scalars and the
 * registration of the annotation natives, whose implementation is src/core/duckvep_core_annotate_run.c. */
#include "duckvep_host.h"

#include "duckvep_model.h"
#include "kernel/src/duckvep_budget.h"
#include "kernel/src/duckvep_event.h"
#include "kernel/src/duckvep_sv.h"
#include "core/duckvep_core_geometry.h"
#include "core/duckvep_core_annotate_run.h"
#include "core/duckvep_core_annotate_types.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool duckvep_register_budget(duckdb_connection connection);

static int
duckvep_validity_is_null(const uint64_t *validity, idx_t row)
{
	return validity != NULL &&
	    ((validity[row / 64] >> (row % 64)) & UINT64_C(1)) == 0;
}

enum duckvep_allele_geometry_field {
	DUCKVEP_GEOMETRY_KIND_CODE = 0,
	DUCKVEP_GEOMETRY_INTERBASE,
	DUCKVEP_GEOMETRY_ANCHOR_SIDE_CODE,
	DUCKVEP_GEOMETRY_RAW_START0,
	DUCKVEP_GEOMETRY_RAW_END0,
	DUCKVEP_GEOMETRY_FEATURE_START0,
	DUCKVEP_GEOMETRY_FEATURE_END0,
	DUCKVEP_GEOMETRY_EDIT_START0,
	DUCKVEP_GEOMETRY_EDIT_END0,
	DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0,
	DUCKVEP_GEOMETRY_REF_DIFF_OFFSET,
	DUCKVEP_GEOMETRY_REF_DIFF_LENGTH,
	DUCKVEP_GEOMETRY_ALT_DIFF_OFFSET,
	DUCKVEP_GEOMETRY_ALT_DIFF_LENGTH,
	DUCKVEP_GEOMETRY_FIELD_COUNT
};

#define duckvep_allele_geometry_field_names ((const char **)duckvep_core_allele_geometry_fields)

static void
duckvep_allele_geometry_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckdb_vector position_vector, reference_vector, alternate_vector;
	duckdb_vector fields[DUCKVEP_GEOMETRY_FIELD_COUNT];
	uint64_t *positions;
	duckdb_string_t *references, *alternates;
	idx_t rows, row;
	size_t field;

	position_vector = duckdb_data_chunk_get_vector(input, 0);
	reference_vector = duckdb_data_chunk_get_vector(input, 1);
	alternate_vector = duckdb_data_chunk_get_vector(input, 2);
	positions = duckdb_vector_get_data(position_vector);
	references = duckdb_vector_get_data(reference_vector);
	alternates = duckdb_vector_get_data(alternate_vector);
	for (field = 0; field < DUCKVEP_GEOMETRY_FIELD_COUNT; field++)
		fields[field] = duckdb_struct_vector_get_child(output, (idx_t)field);
	duckdb_vector_ensure_validity_writable(output);
	duckdb_vector_ensure_validity_writable(
	    fields[DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0]);
	rows = duckdb_data_chunk_get_size(input);
	for (row = 0; row < rows; row++) {
		duckvep_core_allele_geometry_t g;
		const char *reference, *alternate;
		size_t reference_length, alternate_length;

		if (duckvep_validity_is_null(
		        duckdb_vector_get_validity(position_vector), row) ||
		    duckvep_validity_is_null(
		        duckdb_vector_get_validity(reference_vector), row) ||
		    duckvep_validity_is_null(
		        duckdb_vector_get_validity(alternate_vector), row)) {
			duckdb_validity_set_row_invalid(
			    duckdb_vector_get_validity(output), row);
			continue;
		}
		reference = duckdb_string_t_data(&references[row]);
		alternate = duckdb_string_t_data(&alternates[row]);
		reference_length = (size_t)duckdb_string_t_length(references[row]);
		alternate_length = (size_t)duckdb_string_t_length(alternates[row]);
		if (!duckvep_core_allele_geometry(positions[row], reference,
		    reference_length, alternate, alternate_length, &g)) {
			duckdb_scalar_function_set_error(info,
			    duckvep_core_allele_geometry_error);
			return;
		}
		((uint8_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_KIND_CODE]))[row] = g.kind_code;
		((bool *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_INTERBASE]))[row] = g.interbase;
		((uint8_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_ANCHOR_SIDE_CODE]))[row] = g.anchor_side_code;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_RAW_START0]))[row] = g.raw_start0;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_RAW_END0]))[row] = g.raw_end0;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_FEATURE_START0]))[row] = g.feature_start0;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_FEATURE_END0]))[row] = g.feature_end0;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_EDIT_START0]))[row] = g.edit_start0;
		((uint64_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_EDIT_END0]))[row] = g.edit_end0;
		if (g.has_insertion_boundary0) {
			((uint64_t *)duckdb_vector_get_data(
			    fields[DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0]))[row] =
			    g.insertion_boundary0;
		} else {
			duckdb_validity_set_row_invalid(
			    duckdb_vector_get_validity(
			        fields[DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0]), row);
		}
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_REF_DIFF_OFFSET]))[row] =
		    g.reference_difference_offset;
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_REF_DIFF_LENGTH]))[row] =
		    g.reference_difference_length;
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_ALT_DIFF_OFFSET]))[row] =
		    g.alternate_difference_offset;
		((uint16_t *)duckdb_vector_get_data(
		    fields[DUCKVEP_GEOMETRY_ALT_DIFF_LENGTH]))[row] =
		    g.alternate_difference_length;
	}
}

/* A NULL STRUCT row must also invalidate its children, or copying the vector
 * (TRY, projections) reads uninitialized child strings. */
static void
duckvep_breakend_null_row(duckdb_vector output, duckdb_vector *fields, idx_t row)
{
	size_t field;

	duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
	for (field = 0; field < 5; field++)
		duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[field]), row);
}

static void
duckvep_breakend_geometry_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckdb_vector source = duckdb_data_chunk_get_vector(input, 0);
	duckdb_string_t *alternates = duckdb_vector_get_data(source);
	duckdb_vector fields[5];
	idx_t rows = duckdb_data_chunk_get_size(input), row;
	size_t field;

	duckdb_vector_ensure_validity_writable(output);
	for (field = 0; field < 5; field++) {
		fields[field] = duckdb_struct_vector_get_child(output, field);
		duckdb_vector_ensure_validity_writable(fields[field]);
	}
	for (row = 0; row < rows; row++) {
		duckvep_breakend_t parsed;
		duckvep_breakend_status_t status;

		if (duckvep_validity_is_null(duckdb_vector_get_validity(source), row)) {
			duckvep_breakend_null_row(output, fields, row);
			continue;
		}
		status = duckvep_breakend_parse(
		    (const uint8_t *)duckdb_string_t_data(&alternates[row]),
		    (size_t)duckdb_string_t_length(alternates[row]), &parsed);
		if (status == DUCKVEP_BREAKEND_NOT_BREAKEND) {
			duckvep_breakend_null_row(output, fields, row);
			continue;
		}
		if (status != DUCKVEP_BREAKEND_OK) {
			duckdb_scalar_function_set_error(info,
			    duckvep_core_breakend_error(status));
			return;
		}
		duckdb_validity_set_row_valid(duckdb_vector_get_validity(output), row);
		for (field = 0; field < 5; field++)
			duckdb_validity_set_row_valid(duckdb_vector_get_validity(fields[field]), row);
		if (parsed.has_mate) {
			duckdb_vector_assign_string_element_len(fields[0], row,
			    (const char *)parsed.mate_chrom, parsed.mate_chrom_length);
			((uint64_t *)duckdb_vector_get_data(fields[1]))[row] = parsed.mate_position;
			((bool *)duckdb_vector_get_data(fields[3]))[row] = parsed.mate_extends_right != 0u;
		} else {
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[0]), row);
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[1]), row);
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(fields[3]), row);
		}
		((bool *)duckdb_vector_get_data(fields[2]))[row] = parsed.local_join_after != 0u;
		duckdb_vector_assign_string_element_len(fields[4], row,
		    (const char *)parsed.replacement, parsed.replacement_length);
	}
}

static void
duckvep_register_breakend_geometry_scalar(duckdb_connection connection)
{
	const char *names[] = {"mate_chrom", "mate_position", "local_join_after",
	    "mate_extends_right", "replacement_sequence"};
	duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_logical_type position_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	duckdb_logical_type fields[] = {varchar_type, position_type, bool_type, bool_type, varchar_type};
	duckdb_logical_type result = duckdb_create_struct_type(fields, names, 5);
	duckdb_scalar_function scalar = duckdb_create_scalar_function();

	duckdb_scalar_function_set_name(scalar, "duckvep_breakend_geometry");
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_set_return_type(scalar, result);
	duckdb_scalar_function_set_special_handling(scalar);
	duckdb_scalar_function_set_function(scalar, duckvep_breakend_geometry_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);
	duckdb_destroy_logical_type(&result);
	duckdb_destroy_logical_type(&bool_type);
	duckdb_destroy_logical_type(&position_type);
	duckdb_destroy_logical_type(&varchar_type);
}

static void
duckvep_register_allele_geometry_scalar(duckdb_connection connection)
{
	duckdb_scalar_function scalar;
	duckdb_logical_type position_type, varchar_type, result_type;
	duckdb_logical_type field_types[DUCKVEP_GEOMETRY_FIELD_COUNT];
	size_t field;

	scalar = duckdb_create_scalar_function();
	position_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	field_types[DUCKVEP_GEOMETRY_KIND_CODE] =
	    duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	field_types[DUCKVEP_GEOMETRY_INTERBASE] =
	    duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	field_types[DUCKVEP_GEOMETRY_ANCHOR_SIDE_CODE] =
	    duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	for (field = DUCKVEP_GEOMETRY_RAW_START0;
	    field <= DUCKVEP_GEOMETRY_INSERTION_BOUNDARY0; field++)
		field_types[field] =
		    duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	for (field = DUCKVEP_GEOMETRY_REF_DIFF_OFFSET;
	    field < DUCKVEP_GEOMETRY_FIELD_COUNT; field++)
		field_types[field] =
		    duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);
	result_type = duckdb_create_struct_type(field_types,
	    duckvep_allele_geometry_field_names, DUCKVEP_GEOMETRY_FIELD_COUNT);
	duckdb_scalar_function_set_name(scalar, "duckvep_allele_geometry");
	duckdb_scalar_function_add_parameter(scalar, position_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_set_return_type(scalar, result_type);
	duckdb_scalar_function_set_special_handling(scalar);
	duckdb_scalar_function_set_function(scalar,
	    duckvep_allele_geometry_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);
	duckdb_destroy_logical_type(&position_type);
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&result_type);
	for (field = 0; field < DUCKVEP_GEOMETRY_FIELD_COUNT; field++)
		duckdb_destroy_logical_type(&field_types[field]);
}

static duckdb_logical_type
duckvep_list_type(duckvep_result_kind_t kind, int with_hgvs, int with_projection)
{
	duckvep_result_column_t columns[DUCKVEP_RESULT_COLUMNS_MAX];
	duckdb_logical_type types[DUCKVEP_RESULT_COLUMNS_MAX], structure, list;
	const char *names[DUCKVEP_RESULT_COLUMNS_MAX];
	size_t index, count;

	count = duckvep_core_result_columns(kind, with_hgvs, with_projection, columns);
	for (index = 0; index < count; index++) {
		names[index] = columns[index].name;
		types[index] = duckdb_create_logical_type(
		    columns[index].type == DUCKVEP_COLUMN_UINTEGER ? DUCKDB_TYPE_UINTEGER :
		    columns[index].type == DUCKVEP_COLUMN_VARCHAR ? DUCKDB_TYPE_VARCHAR :
		    columns[index].type == DUCKVEP_COLUMN_BOOLEAN ? DUCKDB_TYPE_BOOLEAN :
		    columns[index].type == DUCKVEP_COLUMN_UBIGINT ? DUCKDB_TYPE_UBIGINT :
		    DUCKDB_TYPE_UTINYINT);
	}
	structure = duckdb_create_struct_type(types, names, count);
	list = duckdb_create_list_type(structure);
	for (index = 0; index < count; index++)
		duckdb_destroy_logical_type(&types[index]);
	duckdb_destroy_logical_type(&structure);
	return list;
}

static duckdb_logical_type
duckvep_annotation_list_type(int with_hgvs, int with_projection)
{
	return duckvep_list_type(DUCKVEP_RESULT_RICH, with_hgvs, with_projection);
}

static duckdb_logical_type
duckvep_compact_annotation_list_type(void)
{
	return duckvep_list_type(DUCKVEP_RESULT_COMPACT, 0, 0);
}

static duckdb_logical_type
duckvep_hgvs_annotation_list_type(void)
{
	return duckvep_list_type(DUCKVEP_RESULT_COMPACT_HGVS, 0, 0);
}

static bool
duckvep_register_annotate_scalar(duckdb_connection connection,
	duckvep_registry_t *registry, duckdb_logical_type varchar_type,
	duckdb_logical_type uinteger_type, duckdb_logical_type ubigint_type,
	int compact, duckvep_scalar_event_family_t event_family)
{
	duckdb_scalar_function_set functions;
	duckdb_logical_type result_type;
	bool ok = true;
	const char *name;

	result_type = compact ? duckvep_compact_annotation_list_type() :
	    duckvep_annotation_list_type(0, 0);
	if (event_family == DUCKVEP_SCALAR_STRUCTURAL)
		name = compact ? "_duckvep_annotate_structural_compact" :
		    "_duckvep_annotate_structural_rich";
	else if (event_family == DUCKVEP_SCALAR_BREAKEND)
		name = compact ? "_duckvep_annotate_breakend_compact" :
		    "_duckvep_annotate_breakend_rich";
	else
		name = compact ? "_duckvep_annotate_small_compact" :
		    "_duckvep_annotate_small_rich";
	functions = duckdb_create_scalar_function_set(name);
	for (int distance_parameters = 0; distance_parameters <= 2; distance_parameters++) {
		duckdb_scalar_function scalar = duckdb_create_scalar_function();
		duckdb_scalar_function_set_name(scalar, name);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, uinteger_type);
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		if (event_family == DUCKVEP_SCALAR_STRUCTURAL) {
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
			duckdb_scalar_function_add_parameter(scalar, varchar_type);
			duckdb_scalar_function_add_parameter(scalar, varchar_type);
		} else if (event_family == DUCKVEP_SCALAR_BREAKEND) {
			duckdb_scalar_function_add_parameter(scalar, uinteger_type);
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		} else {
			duckdb_scalar_function_add_parameter(scalar, varchar_type);
			duckdb_scalar_function_add_parameter(scalar, varchar_type);
		}
		if (distance_parameters >= 1)
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		if (distance_parameters >= 2)
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		duckdb_scalar_function_set_return_type(scalar, result_type);
		duckdb_scalar_function_set_volatile(scalar);
		duckvep_registry_retain(registry);
		duckdb_scalar_function_set_extra_info(scalar, registry,
		    duckvep_registry_release);
		if (event_family == DUCKVEP_SCALAR_STRUCTURAL)
			duckdb_scalar_function_set_function(scalar, compact ?
			    duckvep_annotate_sv_compact_scalar : duckvep_annotate_sv_scalar);
		else if (event_family == DUCKVEP_SCALAR_BREAKEND)
			duckdb_scalar_function_set_function(scalar, compact ?
			    duckvep_annotate_breakend_compact_scalar :
			    duckvep_annotate_breakend_scalar);
		else
			duckdb_scalar_function_set_function(scalar, compact ?
			    duckvep_annotate_compact_scalar : duckvep_annotate_scalar);
		ok = duckdb_add_scalar_function_to_set(functions, scalar) == DuckDBSuccess && ok;
		duckdb_destroy_scalar_function(&scalar);
	}
	ok = ok && duckdb_register_scalar_function_set(connection, functions) == DuckDBSuccess;
	duckdb_destroy_scalar_function_set(&functions);
	duckdb_destroy_logical_type(&result_type);
	return ok;
}

static bool
duckvep_register_hgvs_scalar(duckdb_connection connection,
	duckvep_registry_t *registry, duckdb_logical_type varchar_type,
	duckdb_logical_type uinteger_type, duckdb_logical_type ubigint_type,
	int rich)
{
	bool ok = true;
	const char *name = rich ? "_duckvep_annotate_small_rich_hgvs" :
	    "_duckvep_annotate_small_hgvs";
	duckdb_scalar_function_set functions = duckdb_create_scalar_function_set(name);
	duckdb_logical_type result_type = rich ? duckvep_annotation_list_type(1, 0) :
	    duckvep_hgvs_annotation_list_type();

	for (int distance_parameters = 0; distance_parameters <= 2; distance_parameters++) {
		duckdb_scalar_function scalar = duckdb_create_scalar_function();
		duckdb_scalar_function_set_name(scalar, name);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, uinteger_type);
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		if (distance_parameters >= 1)
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		if (distance_parameters >= 2)
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		duckdb_scalar_function_set_return_type(scalar, result_type);
		duckdb_scalar_function_set_volatile(scalar);
		duckvep_registry_retain(registry);
		duckdb_scalar_function_set_extra_info(scalar, registry,
		    duckvep_registry_release);
		duckdb_scalar_function_set_function(scalar, rich ?
		    duckvep_annotate_rich_hgvs_scalar : duckvep_annotate_hgvs_scalar);
		ok = duckdb_add_scalar_function_to_set(functions, scalar) == DuckDBSuccess && ok;
		duckdb_destroy_scalar_function(&scalar);
	}
	ok = ok && duckdb_register_scalar_function_set(connection, functions) == DuckDBSuccess;
	duckdb_destroy_scalar_function_set(&functions);
	duckdb_destroy_logical_type(&result_type);
	return ok;
}

static bool
duckvep_register_projected_scalar(duckdb_connection connection,
	duckvep_registry_t *registry, duckdb_logical_type varchar_type,
	duckdb_logical_type uinteger_type, duckdb_logical_type ubigint_type,
	int with_hgvs)
{
	bool ok = true;
	const char *name = with_hgvs ? "_duckvep_annotate_small_projected_hgvs" :
	    "_duckvep_annotate_small_projected";
	duckdb_scalar_function_set functions = duckdb_create_scalar_function_set(name);
	duckdb_logical_type result_type = duckvep_annotation_list_type(with_hgvs, 1);

	for (int distance_parameters = 0; distance_parameters <= 2; distance_parameters++) {
		duckdb_scalar_function scalar = duckdb_create_scalar_function();
		duckdb_scalar_function_set_name(scalar, name);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, uinteger_type);
		duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		if (distance_parameters >= 1)
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		if (distance_parameters >= 2)
			duckdb_scalar_function_add_parameter(scalar, ubigint_type);
		duckdb_scalar_function_set_return_type(scalar, result_type);
		duckdb_scalar_function_set_volatile(scalar);
		duckvep_registry_retain(registry);
		duckdb_scalar_function_set_extra_info(scalar, registry,
		    duckvep_registry_release);
		duckdb_scalar_function_set_function(scalar, with_hgvs ?
		    duckvep_annotate_projected_hgvs_scalar :
		    duckvep_annotate_projected_scalar);
		ok = duckdb_add_scalar_function_to_set(functions, scalar) == DuckDBSuccess && ok;
		duckdb_destroy_scalar_function(&scalar);
	}
	ok = ok && duckdb_register_scalar_function_set(connection, functions) == DuckDBSuccess;
	duckdb_destroy_scalar_function_set(&functions);
	duckdb_destroy_logical_type(&result_type);
	return ok;
}

bool
register_duckvep_functions(duckdb_connection connection,
	duckdb_database database)
{
	duckvep_registry_t *registry;
	duckdb_logical_type varchar_type, uinteger_type, ubigint_type;

	registry = duckvep_registry_create(database);
	if (registry == NULL)
		return false;
	if (!duckvep_register_budget(connection))
		return false;
	duckvep_register_model_functions(connection, registry);
	duckvep_register_haplotypes(connection, registry);
	duckvep_register_allele_geometry_scalar(connection);
	duckvep_register_breakend_geometry_scalar(connection);

	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	uinteger_type = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
	ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	bool ok = duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, DUCKVEP_SCALAR_SMALL) &&
	    duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, DUCKVEP_SCALAR_SMALL) &&
	    duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0) &&
	    duckvep_register_hgvs_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1) &&
	    duckvep_register_projected_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0) &&
	    duckvep_register_projected_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1) &&
	    duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, DUCKVEP_SCALAR_STRUCTURAL) &&
	    duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, DUCKVEP_SCALAR_STRUCTURAL) &&
	    duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 0, DUCKVEP_SCALAR_BREAKEND) &&
	    duckvep_register_annotate_scalar(connection, registry, varchar_type,
	    uinteger_type, ubigint_type, 1, DUCKVEP_SCALAR_BREAKEND);
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&uinteger_type);
	duckdb_destroy_logical_type(&ubigint_type);

	duckvep_registry_release(registry);
	return ok;
}
