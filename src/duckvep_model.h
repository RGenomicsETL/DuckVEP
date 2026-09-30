/* Private resident-model registry shared by the DuckDB adapter functions. */
#ifndef DUCKVEP_MODEL_H
#define DUCKVEP_MODEL_H

#include "duckdb_extension.h"
#include "core/duckvep_core_model.h"

int duckvep_row_is_null(duckdb_vector, idx_t);
char *duckvep_vector_string(duckdb_vector, idx_t);
int duckvep_vector_string_wellformed(duckdb_vector, idx_t);

duckvep_registry_t *duckvep_registry_create(duckdb_database);
int duckvep_registry_query_acquire(duckvep_registry_t *, char *, size_t);
void duckvep_register_model_functions(duckdb_connection,
	duckvep_registry_t *);
void duckvep_register_haplotypes(duckdb_connection, duckvep_registry_t *);
void duckvep_register_coding_calls(duckdb_connection, duckvep_registry_t *);

#endif /* DUCKVEP_MODEL_H */
