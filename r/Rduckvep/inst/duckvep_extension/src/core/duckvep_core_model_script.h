/* The statements a caller runs to load a model on a host whose model sink is
 * COPY (the v2 host): one COPY per relation into the staging format, then a
 * publish call. The text is host-neutral; only the v2 host serves it. */
#ifndef DUCKVEP_CORE_MODEL_SCRIPT_H
#define DUCKVEP_CORE_MODEL_SCRIPT_H

#include "duckvep_core_sql.h"

#define DUCKVEP_MODEL_SCRIPT_MAX 8

extern const char *const duckvep_model_relations[6];

/* queries: regions, transcripts, exons (required) then mature_mirna, peptide_edits,
 * interval_features (NULL when absent). coverage: -1 absent, 0 false, 1 true.
 * Fills out[0..*count). Returns false on allocation failure. */
bool duckvep_core_model_script(const char *name, const char *const queries[6],
    const char *reference_fasta, int coverage, duckvep_sql_text out[DUCKVEP_MODEL_SCRIPT_MAX],
    size_t *count);

#endif
