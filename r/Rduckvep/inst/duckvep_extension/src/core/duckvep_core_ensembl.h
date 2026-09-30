#ifndef DUCKVEP_CORE_ENSEMBL_H
#define DUCKVEP_CORE_ENSEMBL_H

#include "duckvep_core_cells.h"
#include "duckvep_core_sql.h"

/* duckvep_ensembl_regions_sql and duckvep_ensembl_transcripts_sql take
 * (core_schema, reference_chunks_table, assembly [, {species_id}]);
 * duckvep_ensembl_regulation_features_sql takes (funcgen_schema, regions_table). */
typedef enum {
    DUCKVEP_ENSEMBL_REGIONS, DUCKVEP_ENSEMBL_TRANSCRIPTS, DUCKVEP_ENSEMBL_REGULATION
} duckvep_ensembl_kind_t;

/* The required string arguments: 3, 3, 2. */
size_t duckvep_core_ensembl_required(duckvep_ensembl_kind_t kind);
/* The message when SQL construction fails. */
const char *duckvep_core_ensembl_failed(duckvep_ensembl_kind_t kind);
extern const char duckvep_core_ensembl_names_required[];
extern const char duckvep_core_ensembl_no_options[];
extern const char duckvep_core_ensembl_species_range[];

/* The species_id option as SQL text ("1" when absent, "NULL" for a NULL value).
 * False when it does not fit a BIGINT. `species` holds 40 bytes. */
bool duckvep_core_ensembl_species(const duckvep_cell_t *cell, char *species);

bool duckvep_core_ensembl_sql(duckvep_ensembl_kind_t kind, const char *const values[3],
    const char *species, duckvep_sql_text *out);

/* duckvep_model_receipt_sql: values[0..7] are regions table, model table, then the six
 * literal parameters (NULL means SQL NULL); option_table is regulation_features_table. */
extern const char duckvep_core_receipt_tables_required[];
bool duckvep_core_receipt_sql(const char *const values[8], const char *option_table,
    duckvep_sql_text *out);

#endif
