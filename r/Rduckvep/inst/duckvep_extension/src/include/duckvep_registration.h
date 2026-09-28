#ifndef DUCKVEP_REGISTRATION_H
#define DUCKVEP_REGISTRATION_H

#include "duckdb_extension.h"
#include <stdbool.h>
#include <stddef.h>

typedef struct {
    duckdb_connection connection;
    duckdb_extension_info info;
    struct duckdb_extension_access *access;
} duckvep_registration_t;

bool duckvep_registration_error(duckvep_registration_t *registration, const char *message);
bool duckvep_register_sql(duckvep_registration_t *registration, const char *sql);
bool duckvep_register_sql_parts(duckvep_registration_t *registration,
                               const char *const *parts, size_t count);

#define duckhts_registration_t duckvep_registration_t
#define duckhts_registration_error duckvep_registration_error
#define duckhts_register_sql duckvep_register_sql
#define duckhts_register_sql_parts duckvep_register_sql_parts

#endif
