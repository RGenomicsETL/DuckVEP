/* Text of the preparation builders: duckvep_prepare_sv_geometry_sql,
 * duckvep_prepare_expansionhunter_sql and the BND identity, breakend gene and
 * structural HGVS builders. */
#ifndef DUCKVEP_CORE_PREPARE_H
#define DUCKVEP_CORE_PREPARE_H

#include "duckvep_core_cells.h"
#include "duckvep_core_sql.h"

bool duckvep_core_prepare_sv_sql(const char *relation, duckvep_sql_text *out);
bool duckvep_core_prepare_expansionhunter_sql(const char *events, const char *reference,
    duckvep_sql_text *out);
extern const char duckvep_core_prepare_sv_failed[];
extern const char duckvep_core_prepare_expansionhunter_failed[];

typedef enum {
    DUCKVEP_STRUCTURAL_PAIRS, DUCKVEP_STRUCTURAL_FUSION, DUCKVEP_STRUCTURAL_HGVS
} duckvep_structural_kind_t;

/* Relation arguments: 1, 2, 2. Only the HGVS builder has the max_span option. */
size_t duckvep_core_structural_relations(duckvep_structural_kind_t kind);
bool duckvep_core_structural_has_max_span(duckvep_structural_kind_t kind);
const char *duckvep_core_structural_name(duckvep_structural_kind_t kind);

/* The max_span option: absent or NULL keeps 5000; otherwise 1..60000. False with
 * `message` set when out of range. */
bool duckvep_core_structural_max_span(const duckvep_cell_t *cell, int64_t *max_span,
    char *message, size_t message_size);

bool duckvep_core_structural_sql(duckvep_structural_kind_t kind, char *const relations[2],
    int64_t max_span, duckvep_sql_text *out);

#endif
