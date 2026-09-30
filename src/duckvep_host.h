/* The v1 host layer for src/core/duckvep_core_annotate_run.c: the few DuckDB vector operations
 * the annotation natives use, as the v1 C API calls (so the v1 build compiles to what it did
 * when this code lived in duckvep_annotate.c). The v2 host has its own duckvep_host.h. */
#ifndef DUCKVEP_HOST_H
#define DUCKVEP_HOST_H

#include "duckdb_extension.h"
#include "core/duckvep_core_model.h"
DUCKDB_EXTENSION_EXTERN
#include "duckvep_list.h"

typedef duckdb_vector duckvep_h_vector;
typedef duckdb_data_chunk duckvep_h_chunk;
typedef duckdb_function_info duckvep_h_info;
typedef duckdb_list_entry duckvep_h_list_entry;

#define duckvep_h_data(v) duckdb_vector_get_data(v)
#define duckvep_h_validity(v) duckdb_vector_get_validity(v)
#define duckvep_h_ensure_validity(v) duckdb_vector_ensure_validity_writable(v)
#define duckvep_h_chunk_vector(c, i) duckdb_data_chunk_get_vector((c), (idx_t)(i))
#define duckvep_h_chunk_rows(c) ((size_t)duckdb_data_chunk_get_size(c))
#define duckvep_h_chunk_columns(c) ((size_t)duckdb_data_chunk_get_column_count(c))
#define duckvep_h_struct_child(v, i) duckdb_struct_vector_get_child((v), (idx_t)(i))
#define duckvep_h_list_child(v) duckdb_list_vector_get_child(v)
#define duckvep_h_list_values(v) duckdb_list_vector_get_child(v)
#define duckvep_h_list_reserve(v, n) (duckdb_list_vector_reserve((v), (idx_t)(n)) != DuckDBError)
#define duckvep_h_list_set_size(v, n) (duckdb_list_vector_set_size((v), (idx_t)(n)) != DuckDBError)
#define duckvep_h_assign_string(v, row, text, len) \
	duckdb_vector_assign_string_element_len((v), (idx_t)(row), (text), (idx_t)(len))
#define duckvep_h_set_valid(v, row) duckdb_validity_set_row_valid(duckdb_vector_get_validity(v), (idx_t)(row))
#define duckvep_h_set_null(v, row) duckdb_validity_set_row_invalid(duckdb_vector_get_validity(v), (idx_t)(row))
/* Appends `count` elements to a result list and returns where they start; growth within a chunk is incremental. */
#define duckvep_h_list_extend(v, count, entry) duckvep_list_extend((v), (idx_t)(count), (entry))
#define duckvep_h_set_error(info, message) duckdb_scalar_function_set_error((info), (message))
#define duckvep_h_extra_info(info) duckdb_scalar_function_get_extra_info(info)
#define duckvep_h_vector_size() duckdb_vector_size()

static inline char *
duckvep_h_string(duckdb_vector vector, size_t row)
{
	duckvep_col_t column;

	column.data = duckdb_vector_get_data(vector);
	column.validity = duckdb_vector_get_validity(vector);
	return duckvep_col_string(&column, row);
}

static inline int
duckvep_h_string_wellformed(duckdb_vector vector, size_t row)
{
	duckvep_col_t column;

	column.data = duckdb_vector_get_data(vector);
	column.validity = duckdb_vector_get_validity(vector);
	return duckvep_col_string_wellformed(&column, row);
}

#endif
