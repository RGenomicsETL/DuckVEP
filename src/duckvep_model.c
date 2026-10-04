#include "duckvep_model.h"
#include "core/duckvep_core_snapshot.h"
#include "kernel/src/duckvep_model_internal.h"
#include "kernel/src/duckvep_budget.h"

#include <htslib/faidx.h>

DUCKDB_EXTENSION_EXTERN

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#if defined(_WIN32)
#include <io.h>
#include <share.h>
#include <windows.h>
#else
#include <unistd.h>
#endif


typedef struct duckvep_load_bind {
	duckvep_registry_t *registry;
	char *arguments[4];
	char *mature_mirna_query;
	char *peptide_edit_query;
	char *interval_feature_query;
	char *reference_fasta;
	int transcript_coverage_complete;
} duckvep_load_bind_t;

typedef struct duckvep_load_state {
	int emitted;
} duckvep_load_state_t;

static char *duckvep_string_copy(const char *);

/* ---- v1 row source: a query on the registry's private connection ------------- */

typedef struct duckvep_query_source {
	duckdb_connection connection;
	const char *query;
	duckdb_prepared_statement statement;
	duckdb_result result;
	int have_result;
	duckdb_data_chunk chunk;
	duckvep_batch_t batch;
	duckvep_col_t columns[16];
} duckvep_query_source_t;

static void
duckvep_query_result_close(duckvep_query_source_t *source)
{
	if (source == NULL)
		return;
	if (source->have_result)
		duckdb_destroy_result(&source->result);
	if (source->statement != NULL)
		duckdb_destroy_prepare(&source->statement);
	source->have_result = 0;
	source->statement = NULL;
}

static int
duckvep_query_result_open(duckdb_connection connection, const char *query,
	duckvep_query_source_t *source, char *error, size_t error_size)
{
	duckdb_state state;
	const char *message;

	state = duckdb_prepare(connection, query, &source->statement);
	if (state != DuckDBSuccess) {
		message = duckdb_prepare_error(source->statement);
		duckvep_sql_set_error(error, error_size, message);
		duckvep_query_result_close(source);
		return 0;
	}
	state = duckdb_execute_prepared(source->statement, &source->result);
	source->have_result = 1;
	if (state != DuckDBSuccess) {
		duckvep_sql_set_error(error, error_size,
		    duckdb_result_error(&source->result));
		duckvep_query_result_close(source);
		return 0;
	}
	return 1;
}

int
duckvep_row_is_null(duckdb_vector vector, idx_t row)
{
	uint64_t *validity;

	validity = duckdb_vector_get_validity(vector);
	return validity != NULL &&
	    ((validity[row / 64] >> (row % 64)) & UINT64_C(1)) == 0;
}


static int
duckvep_query_count(duckdb_connection connection, const char *query,
    size_t *count, char *error, size_t error_size)
{
	duckdb_result result;
	duckdb_data_chunk chunk;
	char *sql;
	size_t length;
	duckdb_state state;

	length = strlen(query);
	if (length > SIZE_MAX - sizeof("SELECT CAST(count(*) AS UBIGINT) FROM () q")) {
		duckvep_sql_set_error(error, error_size, "model query is too long");
		return 0;
	}
	sql = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + sizeof("SELECT CAST(count(*) AS UBIGINT) FROM () q"));
	if (sql == NULL) {
		duckvep_sql_set_error(error, error_size, "out of memory counting model rows");
		return 0;
	}
	(void)sprintf(sql, "SELECT CAST(count(*) AS UBIGINT) FROM (%s) q", query);
	memset(&result, 0, sizeof(result));
	state = duckdb_query(connection, sql, &result);
	duckvep_budget_free(sql);
	if (state != DuckDBSuccess) {
		duckvep_sql_set_error(error, error_size, duckdb_result_error(&result));
		duckdb_destroy_result(&result);
		return 0;
	}
	chunk = duckdb_fetch_chunk(result);
	if (chunk == NULL || duckdb_data_chunk_get_size(chunk) != 1 ||
	    duckdb_column_type(&result, 0) != DUCKDB_TYPE_UBIGINT ||
	    duckvep_row_is_null(duckdb_data_chunk_get_vector(chunk, 0), 0) ||
	    *(uint64_t *)duckdb_vector_get_data(
	    duckdb_data_chunk_get_vector(chunk, 0)) > SIZE_MAX) {
		if (chunk != NULL)
			duckdb_destroy_data_chunk(&chunk);
		duckvep_sql_set_error(error, error_size, "invalid model row count");
		duckdb_destroy_result(&result);
		return 0;
	}
	*count = (size_t)*(uint64_t *)duckdb_vector_get_data(
	    duckdb_data_chunk_get_vector(chunk, 0));
	duckdb_destroy_data_chunk(&chunk);
	duckdb_destroy_result(&result);
	return 1;
}
static int
duckvep_qs_count(void *context, size_t *count, char *error, size_t error_size)
{
	duckvep_query_source_t *source = context;

	return duckvep_query_count(source->connection, source->query, count, error, error_size);
}

static int
duckvep_qs_open(void *context, char *error, size_t error_size)
{
	duckvep_query_source_t *source = context;

	return duckvep_query_result_open(source->connection, source->query, source, error, error_size);
}

static size_t
duckvep_qs_columns(void *context)
{
	duckvep_query_source_t *source = context;

	return (size_t)duckdb_column_count(&source->result);
}

static const char *
duckvep_qs_name(void *context, size_t column)
{
	duckvep_query_source_t *source = context;

	return duckdb_column_name(&source->result, (idx_t)column);
}

static duckvep_ctype_t
duckvep_qs_type(void *context, size_t column)
{
	duckvep_query_source_t *source = context;

	switch (duckdb_column_type(&source->result, (idx_t)column)) {
	case DUCKDB_TYPE_BOOLEAN: return DUCKVEP_CT_BOOLEAN;
	case DUCKDB_TYPE_TINYINT: return DUCKVEP_CT_TINYINT;
	case DUCKDB_TYPE_UTINYINT: return DUCKVEP_CT_UTINYINT;
	case DUCKDB_TYPE_UINTEGER: return DUCKVEP_CT_UINTEGER;
	case DUCKDB_TYPE_UBIGINT: return DUCKVEP_CT_UBIGINT;
	case DUCKDB_TYPE_VARCHAR: return DUCKVEP_CT_VARCHAR;
	case DUCKDB_TYPE_BLOB: return DUCKVEP_CT_BLOB;
	default: return DUCKVEP_CT_OTHER;
	}
}

static void
duckvep_qs_release(duckvep_batch_t *batch)
{
	duckvep_query_source_t *source = batch->context;

	duckdb_destroy_data_chunk(&source->chunk);
}

static duckvep_batch_t *
duckvep_qs_next(void *context)
{
	duckvep_query_source_t *source = context;
	idx_t column, columns;

	source->chunk = duckdb_fetch_chunk(source->result);
	if (source->chunk == NULL)
		return NULL;
	columns = duckdb_column_count(&source->result);
	if (columns > 16)
		columns = 16;
	for (column = 0; column < columns; column++) {
		duckdb_vector vector = duckdb_data_chunk_get_vector(source->chunk, column);

		source->columns[column].data = duckdb_vector_get_data(vector);
		source->columns[column].validity = duckdb_vector_get_validity(vector);
	}
	source->batch.rows = (size_t)duckdb_data_chunk_get_size(source->chunk);
	source->batch.columns = source->columns;
	source->batch.context = source;
	source->batch.release = duckvep_qs_release;
	return &source->batch;
}

static void
duckvep_qs_close(void *context)
{
	duckvep_query_result_close(context);
}

static void
duckvep_query_source_init(duckvep_query_source_t *query, duckvep_source_t *source,
	duckdb_connection connection, const char *text)
{
	memset(query, 0, sizeof(*query));
	query->connection = connection;
	query->query = text;
	source->context = query;
	source->count = duckvep_qs_count;
	source->open = duckvep_qs_open;
	source->column_count = duckvep_qs_columns;
	source->column_name = duckvep_qs_name;
	source->column_type = duckvep_qs_type;
	source->next = duckvep_qs_next;
	source->close = duckvep_qs_close;
}


static int
duckvep_query_command(duckdb_connection connection, const char *sql,
	char *error, size_t error_size)
{
	duckdb_result result;
	duckdb_state state;

	memset(&result, 0, sizeof(result));
	state = duckdb_query(connection, sql, &result);
	if (state != DuckDBSuccess)
		duckvep_sql_set_error(error, error_size,
		    duckdb_result_error(&result));
	duckdb_destroy_result(&result);
	return state == DuckDBSuccess;
}

int
duckvep_registry_query_acquire(duckvep_registry_t *registry, char *error, size_t error_size)
{
	/* Query execution can call another DuckVEP table function on a different
	 * DuckDB worker. Waiting for this same connection would deadlock; a busy
	 * preparation slot is a named capacity failure, never recursive execution. */
	if (pthread_mutex_trylock(&registry->query_mutex) == 0)
		return 1;
	duckvep_sql_set_error(error, error_size,
	    "DuckVEP query preparation is busy; nested preparation is unsupported, concurrent callers may retry");
	return 0;
}
static int
duckvep_model_load_queries(duckdb_connection connection,
	const char *region_query, const char *transcript_query,
	const char *exon_query, const char *mature_mirna_query,
	const char *peptide_edit_query, const char *interval_feature_query,
	const char *reference_fasta, int transcript_coverage_complete,
	duckvep_owned_model_t *model,
	char *error, size_t error_size)
{
	duckvep_query_source_t queries[6];
	duckvep_source_t sources[6];
	duckvep_model_sources_t model_sources;
	int ok;

	memset(model, 0, sizeof(*model));
	if (connection == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "model-loading connection is unavailable");
		return 0;
	}
	duckvep_query_source_init(&queries[0], &sources[0], connection, region_query);
	duckvep_query_source_init(&queries[1], &sources[1], connection, transcript_query);
	duckvep_query_source_init(&queries[2], &sources[2], connection, exon_query);
	duckvep_query_source_init(&queries[3], &sources[3], connection, mature_mirna_query);
	duckvep_query_source_init(&queries[4], &sources[4], connection, peptide_edit_query);
	duckvep_query_source_init(&queries[5], &sources[5], connection, interval_feature_query);
	memset(&model_sources, 0, sizeof(model_sources));
	model_sources.regions = &sources[0];
	model_sources.transcripts = &sources[1];
	model_sources.exons = &sources[2];
	model_sources.mature_mirna = mature_mirna_query != NULL ? &sources[3] : NULL;
	model_sources.peptide_edits = peptide_edit_query != NULL ? &sources[4] : NULL;
	model_sources.interval_features = interval_feature_query != NULL ? &sources[5] : NULL;
	model_sources.reference_fasta = reference_fasta;
	model_sources.transcript_coverage_complete = transcript_coverage_complete;
	if (!duckvep_query_command(connection, "BEGIN TRANSACTION", error,
	    error_size))
		return 0;
	ok = duckvep_core_model_load_relations(model, &model_sources, error, error_size);
	if (ok)
		ok = duckvep_query_command(connection, "COMMIT", error, error_size);
	if (!ok)
		(void)duckvep_query_command(connection, "ROLLBACK", NULL, 0);
	if (!ok) {
		duckvep_owned_model_destroy(model);
		return 0;
	}
	return duckvep_core_model_finish(model, error, error_size);
}


/* True for a non-NULL, non-empty string without embedded NUL bytes, i.e. one
 * duckvep_vector_string can only fail to copy for lack of memory. */
int
duckvep_vector_string_wellformed(duckdb_vector vector, idx_t row)
{
	duckdb_string_t *strings;
	uint32_t length;

	if (duckvep_row_is_null(vector, row))
		return 0;
	strings = duckdb_vector_get_data(vector);
	length = duckdb_string_t_length(strings[row]);
	return length != 0 &&
	    memchr(duckdb_string_t_data(&strings[row]), '\0', length) == NULL;
}

char *
duckvep_vector_string(duckdb_vector vector, idx_t row)
{
	duckdb_string_t *strings;
	const char *data;
	uint32_t length;
	char *copy;

	if (duckvep_row_is_null(vector, row))
		return NULL;
	strings = duckdb_vector_get_data(vector);
	length = duckdb_string_t_length(strings[row]);
	data = duckdb_string_t_data(&strings[row]);
	if (memchr(data, '\0', length) != NULL)
		return NULL;
	copy = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, (size_t)length + 1);
	if (copy == NULL)
		return NULL;
	memcpy(copy, data, length);
	copy[length] = '\0';
	return copy;
}

static char *
duckvep_bind_string(duckdb_bind_info info, idx_t parameter)
{
	duckdb_value value;
	char *string;

	value = duckdb_bind_get_parameter(info, parameter);
	if (value == NULL || duckdb_is_null_value(value)) {
		if (value != NULL)
			duckdb_destroy_value(&value);
		return NULL;
	}
	string = duckdb_get_varchar(value);
	duckdb_destroy_value(&value);
	return string;
}

static char *
duckvep_string_copy(const char *source)
{
	size_t length;
	char *copy;

	if (source == NULL)
		return NULL;
	length = strlen(source);
	copy = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
	if (copy != NULL)
		memcpy(copy, source, length + 1);
	return copy;
}

static void
duckvep_model_load_bind_destroy(void *pointer)
{
	duckvep_load_bind_t *bind;
	size_t index;

	bind = pointer;
	if (bind == NULL)
		return;
	for (index = 0; index < 4; index++) {
		if (bind->arguments[index] != NULL)
			duckdb_free(bind->arguments[index]);
	}
	if (bind->mature_mirna_query != NULL)
		duckdb_free(bind->mature_mirna_query);
	if (bind->peptide_edit_query != NULL)
		duckdb_free(bind->peptide_edit_query);
	if (bind->interval_feature_query != NULL)
		duckdb_free(bind->interval_feature_query);
	if (bind->reference_fasta != NULL)
		duckdb_free(bind->reference_fasta);
	duckvep_budget_free(bind);
}

static void
duckvep_model_load_bind(duckdb_bind_info info)
{
	duckvep_load_bind_t *bind;
	duckdb_value complete_value;
	duckdb_value mature_mirna_value;
	duckdb_value peptide_edit_value;
	duckdb_value interval_feature_value;
	duckdb_value reference_fasta_value;
	duckdb_logical_type bool_type;
	idx_t parameter_count;
	size_t index;

	parameter_count = duckdb_bind_get_parameter_count(info);
	if (parameter_count != 4) {
		duckdb_bind_set_error(info,
		    "duckvep_model_load: expected four positional arguments");
		return;
	}
	bind = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*bind));
	if (bind == NULL) {
		duckdb_bind_set_error(info, "duckvep_model_load: out of memory");
		return;
	}
	bind->registry = duckdb_bind_get_extra_info(info);
	for (index = 0; index < 4; index++) {
		bind->arguments[index] = duckvep_bind_string(info, (idx_t)index);
		if (bind->arguments[index] == NULL ||
		    bind->arguments[index][0] == '\0') {
			duckdb_bind_set_error(info,
			    "duckvep_model_load: arguments must be non-empty strings");
			duckvep_model_load_bind_destroy(bind);
			return;
		}
	}
	mature_mirna_value = duckdb_bind_get_named_parameter(
	    info, "mature_mirna_query");
	if (mature_mirna_value != NULL) {
		if (duckdb_is_null_value(mature_mirna_value)) {
			duckdb_destroy_value(&mature_mirna_value);
		} else {
			bind->mature_mirna_query = duckdb_get_varchar(
			    mature_mirna_value);
			duckdb_destroy_value(&mature_mirna_value);
			if (bind->mature_mirna_query == NULL ||
			    bind->mature_mirna_query[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: mature_mirna_query must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	peptide_edit_value = duckdb_bind_get_named_parameter(
	    info, "peptide_edit_query");
	if (peptide_edit_value != NULL) {
		if (duckdb_is_null_value(peptide_edit_value)) {
			duckdb_destroy_value(&peptide_edit_value);
		} else {
			bind->peptide_edit_query = duckdb_get_varchar(
			    peptide_edit_value);
			duckdb_destroy_value(&peptide_edit_value);
			if (bind->peptide_edit_query == NULL ||
			    bind->peptide_edit_query[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: peptide_edit_query must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	interval_feature_value = duckdb_bind_get_named_parameter(
	    info, "interval_feature_query");
	if (interval_feature_value != NULL) {
		if (duckdb_is_null_value(interval_feature_value)) {
			duckdb_destroy_value(&interval_feature_value);
		} else {
			bind->interval_feature_query = duckdb_get_varchar(
			    interval_feature_value);
			duckdb_destroy_value(&interval_feature_value);
			if (bind->interval_feature_query == NULL ||
			    bind->interval_feature_query[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: interval_feature_query must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	reference_fasta_value = duckdb_bind_get_named_parameter(
	    info, "reference_fasta");
	if (reference_fasta_value != NULL) {
		if (duckdb_is_null_value(reference_fasta_value)) {
			duckdb_destroy_value(&reference_fasta_value);
		} else {
			bind->reference_fasta = duckdb_get_varchar(
			    reference_fasta_value);
			duckdb_destroy_value(&reference_fasta_value);
			if (bind->reference_fasta == NULL ||
			    bind->reference_fasta[0] == '\0') {
				duckdb_bind_set_error(info,
				    "duckvep_model_load: reference_fasta must be a non-empty string");
				duckvep_model_load_bind_destroy(bind);
				return;
			}
		}
	}
	complete_value = duckdb_bind_get_named_parameter(
	    info, "transcript_coverage_complete");
	if (complete_value != NULL) {
		if (duckdb_is_null_value(complete_value)) {
			duckdb_destroy_value(&complete_value);
			duckdb_bind_set_error(info,
			    "duckvep_model_load: transcript_coverage_complete cannot be NULL");
			duckvep_model_load_bind_destroy(bind);
			return;
		}
		bind->transcript_coverage_complete = duckdb_get_bool(complete_value);
		duckdb_destroy_value(&complete_value);
	}
	bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	duckdb_bind_add_result_column(info, "loaded", bool_type);
	duckdb_destroy_logical_type(&bool_type);
	duckdb_bind_set_bind_data(info, bind,
	    duckvep_model_load_bind_destroy);
}

static void
duckvep_model_load_state_destroy(void *pointer)
{
	duckvep_budget_free(pointer);
}

static void
duckvep_model_load_init(duckdb_init_info info)
{
	duckvep_load_bind_t *bind;
	duckvep_load_state_t *state;
	duckvep_registry_t *registry;
	duckvep_model_entry_t *entry;
	char error[DUCKVEP_SQL_ERROR_SIZE];
	char final_error[DUCKVEP_SQL_ERROR_SIZE + 256];
	int loaded;

	duckdb_init_set_max_threads(info, 1);
	duckvep_budget_clear_failure();
	bind = duckdb_init_get_bind_data(info);
	if (bind == NULL || bind->registry == NULL) {
		duckdb_init_set_error(info,
		    "duckvep_model_load: missing bind state");
		return;
	}
	state = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*state));
	if (state == NULL) {
		duckdb_init_set_error(info, duckvep_sql_final_error(final_error,
		    sizeof(final_error), NULL, "duckvep_model_load: out of memory"));
		return;
	}
	registry = bind->registry;
	pthread_mutex_lock(&registry->mutex);
	entry = duckvep_registry_find_locked(registry, bind->arguments[0]);
	pthread_mutex_unlock(&registry->mutex);
	if (entry != NULL) {
		duckvep_budget_free(state);
		duckdb_init_set_error(info,
		    "duckvep_model_load: model name already exists");
		return;
	}
	entry = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*entry));
	if (entry == NULL ||
	    (entry->name = duckvep_string_copy(bind->arguments[0])) == NULL) {
		duckvep_model_entry_destroy(entry);
		duckvep_budget_free(state);
		duckdb_init_set_error(info, duckvep_sql_final_error(final_error,
		    sizeof(final_error), NULL, "duckvep_model_load: out of memory"));
		return;
	}
	memset(error, 0, sizeof(error));
	if (!duckvep_registry_query_acquire(registry, error, sizeof(error))) {
		duckvep_model_entry_destroy(entry);
		duckvep_budget_free(state);
		duckdb_init_set_error(info, error);
		return;
	}
	loaded = duckvep_model_load_queries(registry->query_connection,
	    bind->arguments[1], bind->arguments[2], bind->arguments[3],
	    bind->mature_mirna_query,
	    bind->peptide_edit_query,
	    bind->interval_feature_query,
	    bind->reference_fasta, bind->transcript_coverage_complete,
	    &entry->model, error,
	    sizeof(error));
	pthread_mutex_unlock(&registry->query_mutex);
	if (!loaded) {
		duckvep_model_entry_destroy(entry);
		duckvep_budget_free(state);
		duckdb_init_set_error(info, duckvep_sql_final_error(final_error,
		    sizeof(final_error), error,
		    "duckvep_model_load: model load failed"));
		return;
	}
	pthread_mutex_lock(&registry->mutex);
	if (duckvep_registry_find_locked(registry, entry->name) != NULL) {
		pthread_mutex_unlock(&registry->mutex);
		duckvep_model_entry_destroy(entry);
		duckvep_budget_free(state);
		duckdb_init_set_error(info,
		    "duckvep_model_load: model name was created concurrently");
		return;
	}
	entry->next = registry->models;
	registry->models = entry;
	pthread_mutex_unlock(&registry->mutex);
	duckdb_init_set_init_data(info, state, duckvep_model_load_state_destroy);
}

static void
duckvep_model_load_scan(duckdb_function_info info, duckdb_data_chunk output)
{
	duckvep_load_state_t *state;
	duckdb_vector vector;
	bool *values;

	state = duckdb_function_get_init_data(info);
	if (state->emitted) {
		duckdb_data_chunk_set_size(output, 0);
		return;
	}
	vector = duckdb_data_chunk_get_vector(output, 0);
	values = duckdb_vector_get_data(vector);
	values[0] = true;
	state->emitted = 1;
	duckdb_data_chunk_set_size(output, 1);
}

static void
duckvep_model_drop_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_registry_t *registry;
	duckdb_vector name_vector;
	bool *values;
	idx_t row, rows;

	registry = duckdb_scalar_function_get_extra_info(info);
	rows = duckdb_data_chunk_get_size(input);
	name_vector = duckdb_data_chunk_get_vector(input, 0);
	values = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		duckvep_model_entry_t *entry, *previous;
		char *name;

		name = duckvep_vector_string(name_vector, row);
		if (name == NULL || *name == '\0') {
			duckvep_budget_free(name);
			duckdb_scalar_function_set_error(info,
			    "duckvep_model_drop: name must be a non-empty string");
			return;
		}
		entry = NULL;
		previous = NULL;
		pthread_mutex_lock(&registry->mutex);
		for (entry = registry->models; entry != NULL;
		    previous = entry, entry = entry->next) {
			if (strcmp(entry->name, name) == 0)
				break;
		}
		if (entry != NULL && entry->pins == 0) {
			if (previous != NULL)
				previous->next = entry->next;
			else
				registry->models = entry->next;
			entry->next = NULL;
		} else {
			entry = NULL;
		}
		pthread_mutex_unlock(&registry->mutex);
		duckvep_budget_free(name);
		values[row] = entry != NULL;
		duckvep_model_entry_destroy(entry);
	}
}



/* duckvep_model_save(name, path): writes the loaded model's arrays to a snapshot file. */
static void
duckvep_model_save_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_registry_t *registry;
	duckdb_vector name_vector, path_vector;
	bool *values;
	idx_t row, rows;

	registry = duckdb_scalar_function_get_extra_info(info);
	rows = duckdb_data_chunk_get_size(input);
	name_vector = duckdb_data_chunk_get_vector(input, 0);
	path_vector = duckdb_data_chunk_get_vector(input, 1);
	values = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		duckvep_model_entry_t *entry;
		char error[DUCKVEP_SQL_ERROR_SIZE];
		char *name, *path;
		int saved;

		name = duckvep_vector_string(name_vector, row);
		path = duckvep_vector_string(path_vector, row);
		if (name == NULL || *name == '\0' || path == NULL || *path == '\0') {
			duckvep_budget_free(name);
			duckvep_budget_free(path);
			duckdb_scalar_function_set_error(info,
			    "duckvep_model_save: name and path must be non-empty strings");
			return;
		}
		entry = duckvep_registry_pin(registry, name);
		duckvep_budget_free(name);
		if (entry == NULL) {
			duckvep_budget_free(path);
			duckdb_scalar_function_set_error(info, "duckvep_model_save: unknown model name");
			return;
		}
		memset(error, 0, sizeof(error));
		saved = duckvep_core_model_snapshot_save(&entry->model, path, error, sizeof(error));
		duckvep_registry_unpin(registry, entry);
		duckvep_budget_free(path);
		if (!saved) {
			duckdb_scalar_function_set_error(info, error);
			return;
		}
		values[row] = true;
	}
}

/* duckvep_model_restore(name, path): maps a snapshot and installs it as a model. */
static void
duckvep_model_restore_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_registry_t *registry;
	duckdb_vector name_vector, path_vector;
	bool *values;
	idx_t row, rows;

	registry = duckdb_scalar_function_get_extra_info(info);
	rows = duckdb_data_chunk_get_size(input);
	name_vector = duckdb_data_chunk_get_vector(input, 0);
	path_vector = duckdb_data_chunk_get_vector(input, 1);
	values = duckdb_vector_get_data(output);
	for (row = 0; row < rows; row++) {
		char error[DUCKVEP_SQL_ERROR_SIZE + 256];
		char *name, *path;
		int restored;

		name = duckvep_vector_string(name_vector, row);
		path = duckvep_vector_string(path_vector, row);
		if (name == NULL || *name == '\0' || path == NULL || *path == '\0') {
			duckvep_budget_free(name);
			duckvep_budget_free(path);
			duckdb_scalar_function_set_error(info,
			    "duckvep_model_restore: name and path must be non-empty strings");
			return;
		}
		memset(error, 0, sizeof(error));
		restored = duckvep_core_model_snapshot_install(registry, name, path, error, sizeof(error));
		duckvep_budget_free(name);
		duckvep_budget_free(path);
		if (!restored) {
			duckdb_scalar_function_set_error(info, error);
			return;
		}
		values[row] = true;
	}
}

/* Internal diagnostic: FNV-1a of the loaded arrays, for comparing loads across hosts. */
static void
duckvep_model_fingerprint_scalar(duckdb_function_info info,
	duckdb_data_chunk input, duckdb_vector output)
{
	duckvep_registry_t *registry;
	duckdb_vector name_vector;
	uint64_t *values;
	idx_t row, rows;

	registry = duckdb_scalar_function_get_extra_info(info);
	rows = duckdb_data_chunk_get_size(input);
	name_vector = duckdb_data_chunk_get_vector(input, 0);
	values = duckdb_vector_get_data(output);
	duckdb_vector_ensure_validity_writable(output);
	for (row = 0; row < rows; row++) {
		duckvep_model_entry_t *entry;
		char *name;

		name = duckvep_vector_string(name_vector, row);
		entry = name != NULL ? duckvep_registry_pin(registry, name) : NULL;
		duckvep_budget_free(name);
		if (entry == NULL) {
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
			continue;
		}
		values[row] = duckvep_core_model_fingerprint(&entry->model);
		duckvep_registry_unpin(registry, entry);
	}
}

static void
duckvep_v1_registry_release(duckvep_registry_t *registry)
{
	duckdb_connection connection = registry->query_connection;

	if (connection != NULL)
		duckdb_disconnect(&connection);
	registry->query_connection = NULL;
}

duckvep_registry_t *
duckvep_registry_create(duckdb_database database)
{
	duckvep_registry_t *registry;
	duckdb_connection connection = NULL;

	registry = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*registry));
	if (registry == NULL)
		return NULL;
	(void)pthread_mutex_init(&registry->mutex, NULL);
	(void)pthread_mutex_init(&registry->query_mutex, NULL);
	(void)pthread_cond_init(&registry->admission, NULL);
	registry->references = 1;
	registry->host_release = duckvep_v1_registry_release;
	if (duckdb_connect(database, &connection) !=
	    DuckDBSuccess || connection == NULL) {
		pthread_cond_destroy(&registry->admission);
		pthread_mutex_destroy(&registry->query_mutex);
		pthread_mutex_destroy(&registry->mutex);
		duckvep_budget_free(registry);
		return NULL;
	}
	registry->query_connection = connection;
	return registry;
}

void
duckvep_register_model_functions(duckdb_connection connection,
	duckvep_registry_t *registry)
{
	duckdb_table_function table;
	duckdb_scalar_function scalar;
	duckdb_logical_type varchar_type, bool_type;

	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
	table = duckdb_create_table_function();
	duckdb_table_function_set_name(table, "duckvep_model_load");
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_parameter(table, varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "mature_mirna_query", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "peptide_edit_query", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "interval_feature_query", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "reference_fasta", varchar_type);
	duckdb_table_function_add_named_parameter(table,
	    "transcript_coverage_complete", bool_type);
	duckvep_registry_retain(registry);
	duckdb_table_function_set_extra_info(table, registry,
	    duckvep_registry_release);
	duckdb_table_function_set_bind(table, duckvep_model_load_bind);
	duckdb_table_function_set_init(table, duckvep_model_load_init);
	duckdb_table_function_set_function(table, duckvep_model_load_scan);
	(void)duckdb_register_table_function(connection, table);
	duckdb_destroy_table_function(&table);

	scalar = duckdb_create_scalar_function();
	duckdb_scalar_function_set_name(scalar, "duckvep_model_drop");
	duckdb_scalar_function_add_parameter(scalar, varchar_type);
	duckdb_scalar_function_set_return_type(scalar, bool_type);
	duckdb_scalar_function_set_volatile(scalar);
	duckvep_registry_retain(registry);
	duckdb_scalar_function_set_extra_info(scalar, registry,
	    duckvep_registry_release);
	duckdb_scalar_function_set_function(scalar, duckvep_model_drop_scalar);
	(void)duckdb_register_scalar_function(connection, scalar);
	duckdb_destroy_scalar_function(&scalar);

	{
		static const struct { const char *name; duckdb_scalar_function_t function; } snapshots[2] = {
			{"duckvep_model_save", duckvep_model_save_scalar},
			{"duckvep_model_restore", duckvep_model_restore_scalar}};
		size_t i;

		for (i = 0; i < 2; i++) {
			scalar = duckdb_create_scalar_function();
			duckdb_scalar_function_set_name(scalar, snapshots[i].name);
			duckdb_scalar_function_add_parameter(scalar, varchar_type);
			duckdb_scalar_function_add_parameter(scalar, varchar_type);
			duckdb_scalar_function_set_return_type(scalar, bool_type);
			duckdb_scalar_function_set_volatile(scalar);
			duckvep_registry_retain(registry);
			duckdb_scalar_function_set_extra_info(scalar, registry,
			    duckvep_registry_release);
			duckdb_scalar_function_set_function(scalar, snapshots[i].function);
			(void)duckdb_register_scalar_function(connection, scalar);
			duckdb_destroy_scalar_function(&scalar);
		}
	}

	{
		duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);

		scalar = duckdb_create_scalar_function();
		duckdb_scalar_function_set_name(scalar, "_duckvep_model_fingerprint");
		duckdb_scalar_function_add_parameter(scalar, varchar_type);
		duckdb_scalar_function_set_return_type(scalar, ubigint_type);
		duckdb_scalar_function_set_volatile(scalar);
		duckdb_scalar_function_set_special_handling(scalar);
		duckvep_registry_retain(registry);
		duckdb_scalar_function_set_extra_info(scalar, registry,
		    duckvep_registry_release);
		duckdb_scalar_function_set_function(scalar, duckvep_model_fingerprint_scalar);
		(void)duckdb_register_scalar_function(connection, scalar);
		duckdb_destroy_scalar_function(&scalar);
		duckdb_destroy_logical_type(&ubigint_type);
	}
	duckdb_destroy_logical_type(&varchar_type);
	duckdb_destroy_logical_type(&bool_type);
}
