#define DUCKDB_EXTENSION_NAME duckvep
#include "duckdb_extension.h"
#include "duckvep_registration.h"
#include "duckvep_sql.h"
#include <stdio.h>

DUCKDB_EXTENSION_EXTERN

extern bool register_duckvep_functions(duckdb_connection, duckdb_database);
extern bool register_duckvep_sql_kernels(duckvep_registration_t *);
extern bool register_duckvep_ensembl_functions(duckvep_registration_t *);
extern bool register_duckvep_sql_functions(duckvep_registration_t *);

DUCKDB_EXTENSION_ENTRYPOINT(duckdb_connection connection,
                            duckdb_extension_info info,
                            struct duckdb_extension_access *access) {
    duckvep_registration_t registration = {connection, info, access};
    const char *version = duckdb_library_version();
    unsigned major = 0, minor = 0;
    if (sscanf(version, "v%u.%u", &major, &minor) != 2 ||
        major < 1 || (major == 1 && minor < 4)) {
        char message[160];
        snprintf(message, sizeof(message),
                 "DuckVEP requires DuckDB 1.4.0 or newer; loaded runtime is %s", version);
        return duckvep_registration_error(&registration, message);
    }
    return register_duckvep_functions(connection, *access->get_database(info)) &&
           register_duckvep_sql_kernels(&registration) &&
           register_duckvep_sql_functions(&registration) &&
           register_duckvep_ensembl_functions(&registration);
}
