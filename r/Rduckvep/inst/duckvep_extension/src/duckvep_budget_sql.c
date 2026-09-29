/* SQL surface for the native allocation budget and bounded-execution limits. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "kernel/src/duckvep_budget.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

typedef struct {
	idx_t offset;
	duckvep_budget_stats_t stats;
} duckvep_budget_scan_t;

static void
duckvep_budget_bind(duckdb_bind_info info)
{
	duckdb_logical_type varchar_type, ubigint_type;

	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	duckdb_bind_add_result_column(info, "owner", varchar_type);
	duckdb_bind_add_result_column(info, "current_bytes", ubigint_type);
	duckdb_bind_add_result_column(info, "high_water_bytes", ubigint_type);
	duckdb_bind_add_result_column(info, "limit_bytes", ubigint_type);
	duckdb_bind_add_result_column(info, "charges", ubigint_type);
	duckdb_bind_add_result_column(info, "refusals", ubigint_type);
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&ubigint_type);
}

static void
duckvep_budget_scan_destroy(void *data)
{
	duckvep_budget_free(data);
}

static void
duckvep_budget_init(duckdb_init_info info)
{
	duckvep_budget_scan_t *scan;

	duckdb_init_set_max_threads(info, 1);
	scan = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*scan));
	if (scan == NULL) {
		duckdb_init_set_error(info,
		    "duckvep_native_budget: capacity error: cannot allocate scan state");
		return;
	}
	duckvep_budget_stats(&scan->stats);
	duckdb_init_set_init_data(info, scan, duckvep_budget_scan_destroy);
}

static void
duckvep_budget_scan(duckdb_function_info info, duckdb_data_chunk output)
{
	duckvep_budget_scan_t *scan;
	uint64_t *current, *high, *limit, *charges, *refusals;
	idx_t count;

	scan = duckdb_function_get_init_data(info);
	current = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 1));
	high = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2));
	limit = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 3));
	charges = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 4));
	refusals = duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 5));
	count = 0;
	while (scan->offset <= DUCKVEP_OWNER_COUNT) {
		int total = scan->offset == DUCKVEP_OWNER_COUNT;

		duckdb_vector_assign_string_element(
		    duckdb_data_chunk_get_vector(output, 0), count,
		    total ? "total" :
		    duckvep_budget_owner_name((duckvep_budget_owner_t)scan->offset));
		current[count] = total ? scan->stats.total_current :
		    scan->stats.current[scan->offset];
		high[count] = total ? scan->stats.total_high_water :
		    scan->stats.high_water[scan->offset];
		limit[count] = scan->stats.limit;
		/* Process-wide counters: reported on the total row, zero elsewhere. */
		charges[count] = total ? scan->stats.charges : 0;
		refusals[count] = total ? scan->stats.refusals : 0;
		count++;
		scan->offset++;
	}
	duckdb_data_chunk_set_size(output, count);
}

static bool
duckvep_budget_int_arg(duckdb_data_chunk input, idx_t column, idx_t row,
	int64_t *value)
{
	duckdb_vector vector;
	uint64_t *validity;

	vector = duckdb_data_chunk_get_vector(input, column);
	validity = duckdb_vector_get_validity(vector);
	if (validity != NULL && !duckdb_validity_row_is_valid(validity, row))
		return false;
	*value = ((int64_t *)duckdb_vector_get_data(vector))[row];
	return true;
}

static void
duckvep_budget_set_scalar(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	idx_t row, rows;
	int64_t *result;

	rows = duckdb_data_chunk_get_size(input);
	result = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		int64_t value;
		char message[256];

		if (!duckvep_budget_int_arg(input, 0, row, &value) || value <= 0) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_native_budget_set: the budget must be a positive byte count");
			return;
		}
		if (!duckvep_budget_set_limit((uint64_t)value)) {
			duckvep_budget_stats_t stats;

			duckvep_budget_stats(&stats);
			(void)snprintf(message, sizeof(message),
			    "duckvep_native_budget_set: capacity error: %lld bytes is below the %llu bytes already in use",
			    (long long)value, (unsigned long long)stats.total_current);
			duckdb_scalar_function_set_error(info, message);
			return;
		}
		result[row] = value;
	}
}

static void
duckvep_worker_limits_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	idx_t row, rows;
	bool *result;

	rows = duckdb_data_chunk_get_size(input);
	result = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		int64_t workers, scratch, emit, idle;
		duckvep_budget_worker_limits_t limits;

		if (!duckvep_budget_int_arg(input, 0, row, &workers) ||
		    !duckvep_budget_int_arg(input, 1, row, &scratch) ||
		    !duckvep_budget_int_arg(input, 2, row, &emit) ||
		    !duckvep_budget_int_arg(input, 3, row, &idle) ||
		    workers <= 0 || workers > 1024 || scratch < 0 || emit < 0 ||
		    idle < 0) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_worker_limits_set: workers must be 1..1024 and byte limits non-negative");
			return;
		}
		limits.max_workers = (uint32_t)workers;
		limits.scratch_bytes = (uint64_t)scratch;
		limits.emit_bytes = (uint64_t)emit;
		limits.idle_bytes = (uint64_t)idle;
		(void)duckvep_budget_set_worker_limits(&limits);
		result[row] = true;
	}
}

static void
duckvep_budget_reset_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	idx_t row, rows;
	bool *result;

	(void)info;
	rows = duckdb_data_chunk_get_size(input);
	result = duckdb_vector_get_data(output);
	duckvep_budget_reset_high_water();
	for (row = 0; row < rows; row++)
		result[row] = true;
}

#ifdef DUCKVEP_FAULT_INJECTION
/* Test builds only: arm(n) fails the nth budget allocation from now on and
 * returns how many allocations were counted since the previous arm. */
static void
duckvep_fault_arm_scalar(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	idx_t row, rows;
	int64_t *result;

	rows = duckdb_data_chunk_get_size(input);
	result = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		int64_t nth;

		if (!duckvep_budget_int_arg(input, 0, row, &nth) || nth < 0) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_fault_arm: expected a non-negative index");
			return;
		}
		result[row] = (int64_t)duckvep_fault_allocations();
		duckvep_fault_arm((uint64_t)nth);
	}
}

static void
duckvep_fault_count_scalar(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	idx_t row, rows;
	int64_t *result;

	(void)info;
	rows = duckdb_data_chunk_get_size(input);
	result = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++)
		result[row] = (int64_t)duckvep_fault_allocations();
}

static void
duckvep_fault_fired_scalar(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	idx_t row, rows;
	int64_t *result;

	(void)info;
	rows = duckdb_data_chunk_get_size(input);
	result = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++)
		result[row] = (int64_t)duckvep_fault_fired();
}
#endif

static bool
duckvep_register_bigint_scalar(duckdb_connection connection, const char *name,
	idx_t parameters, duckdb_type return_type,
	duckdb_scalar_function_t callback)
{
	duckdb_scalar_function function;
	duckdb_logical_type bigint_type, result_type;
	duckdb_state state;
	idx_t index;

	function = duckdb_create_scalar_function();
	bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
	result_type = duckdb_create_logical_type(return_type);
	duckdb_scalar_function_set_name(function, name);
	for (index = 0; index < parameters; index++)
		duckdb_scalar_function_add_parameter(function, bigint_type);
	duckdb_scalar_function_set_return_type(function, result_type);
	duckdb_scalar_function_set_volatile(function);
	duckdb_scalar_function_set_function(function, callback);
	state = duckdb_register_scalar_function(connection, function);
	duckdb_destroy_scalar_function(&function);
	duckdb_destroy_logical_type(&bigint_type);
	duckdb_destroy_logical_type(&result_type);
	return state == DuckDBSuccess;
}

bool
duckvep_register_budget(duckdb_connection connection)
{
	duckdb_table_function function;
	duckdb_state state;
	bool ok;

	function = duckdb_create_table_function();
	duckdb_table_function_set_name(function, "duckvep_native_budget");
	duckdb_table_function_set_bind(function, duckvep_budget_bind);
	duckdb_table_function_set_init(function, duckvep_budget_init);
	duckdb_table_function_set_function(function, duckvep_budget_scan);
	state = duckdb_register_table_function(connection, function);
	duckdb_destroy_table_function(&function);
	ok = state == DuckDBSuccess &&
	    duckvep_register_bigint_scalar(connection, "duckvep_native_budget_set",
	    1, DUCKDB_TYPE_BIGINT, duckvep_budget_set_scalar) &&
	    duckvep_register_bigint_scalar(connection, "duckvep_worker_limits_set",
	    4, DUCKDB_TYPE_BOOLEAN, duckvep_worker_limits_scalar) &&
	    duckvep_register_bigint_scalar(connection,
	    "duckvep_native_budget_reset_high_water", 0, DUCKDB_TYPE_BOOLEAN,
	    duckvep_budget_reset_scalar);
#ifdef DUCKVEP_FAULT_INJECTION
	ok = ok &&
	    duckvep_register_bigint_scalar(connection, "duckvep_fault_arm", 1,
	    DUCKDB_TYPE_BIGINT, duckvep_fault_arm_scalar) &&
	    duckvep_register_bigint_scalar(connection, "duckvep_fault_allocations",
	    0, DUCKDB_TYPE_BIGINT, duckvep_fault_count_scalar) &&
	    duckvep_register_bigint_scalar(connection, "duckvep_fault_fired", 0,
	    DUCKDB_TYPE_BIGINT, duckvep_fault_fired_scalar);
#endif
	return ok;
}
