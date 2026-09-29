#define DUCKDB_EXTENSION_NAME duckvep
#include "duckdb_extension.h"
#include "duckvep_registration.h"
#include "duckvep_sql.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

DUCKDB_EXTENSION_EXTERN

static void duckvep_revcomp(duckdb_function_info info, duckdb_data_chunk input,
                            duckdb_vector output) {
    duckdb_vector sequences = duckdb_data_chunk_get_vector(input, 0);
    duckdb_string_t *values = duckdb_vector_get_data(sequences);
    const char *source = "ACGTRYSWKMBDHVNacgtryswkmbdhvn";
    const char *target = "TGCAYRSWMKVHDBNtgcayrswmkvhdbn";
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        idx_t length = duckdb_string_t_length(values[row]);
        const char *sequence = duckdb_string_t_data(&values[row]);
        char *reversed = malloc(length ? length : 1);
        if (!reversed) {
            duckdb_scalar_function_set_error(info, "_duckvep_revcomp: allocation failed");
            return;
        }
        idx_t end = length, at = 0;
        while (end) {
            idx_t first = end - 1;
            while (first && ((unsigned char)sequence[first] & 0xc0) == 0x80)
                first--;
            idx_t width = end - first;
            if (width == 1) {
                char base = sequence[first];
                const char *match = base ? strchr(source, base) : NULL;
                reversed[at] = match ? target[match - source] : base;
            } else {
                memcpy(reversed + at, sequence + first, width);
            }
            at += width;
            end = first;
        }
        duckdb_vector_assign_string_element_len(output, row, reversed, length);
        free(reversed);
    }
}

static bool duckvep_register_revcomp(duckdb_connection connection) {
    duckdb_scalar_function function = duckdb_create_scalar_function();
    duckdb_logical_type text = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_scalar_function_set_name(function, "_duckvep_revcomp");
    duckdb_scalar_function_add_parameter(function, text);
    duckdb_scalar_function_set_return_type(function, text);
    duckdb_scalar_function_set_function(function, duckvep_revcomp);
    duckdb_state state = duckdb_register_scalar_function(connection, function);
    duckdb_destroy_scalar_function(&function);
    duckdb_destroy_logical_type(&text);
    return state == DuckDBSuccess;
}
extern bool register_duckvep_functions(duckdb_connection, duckdb_database);
extern bool register_duckvep_sql_kernels(duckvep_registration_t *);
extern bool register_duckvep_ensembl_functions(duckvep_registration_t *);
extern bool register_duckvep_sql_functions(duckvep_registration_t *);
extern bool register_duckvep_prepare_sql(duckdb_connection);
extern bool register_duckvep_structural_sql(duckdb_connection);

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
    return duckvep_register_revcomp(connection) &&
           register_duckvep_functions(connection, *access->get_database(info)) &&
           register_duckvep_sql_kernels(&registration) &&
           register_duckvep_sql_functions(&registration) &&
           register_duckvep_prepare_sql(connection) &&
           register_duckvep_structural_sql(connection) &&
           register_duckvep_ensembl_functions(&registration);
}
