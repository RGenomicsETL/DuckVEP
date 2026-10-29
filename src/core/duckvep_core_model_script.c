#include "duckvep_core_model_script.h"

#include <inttypes.h>
#include <stdio.h>

const char *const duckvep_model_relations[7] = {
    "regions", "transcripts", "exons", "mature_mirna", "peptide_edits", "interval_features", "rna_edits"};

bool duckvep_core_model_script(const char *name, const char *const queries[7],
    const char *reference_fasta, int coverage, duckvep_sql_text out[DUCKVEP_MODEL_SCRIPT_MAX],
    size_t *count)
{
    static uint64_t scripts_built;
    bool ok = true;
    size_t n = 0;
    char token[24];
    duckvep_sql_text key = {0};

    /* Staged relations are shared by every connection to the database. Each script stages
     * under its own key, name<US>sequence, and publishes that key, so two loads of one
     * name cannot publish each other's rows. */
    (void)snprintf(token, sizeof token, "\x1f%" PRIu64,
                   __atomic_add_fetch(&scripts_built, 1, __ATOMIC_RELAXED));
    ok = duckvep_sql_append(&key, name) && duckvep_sql_append(&key, token);
    for (size_t i = 0; ok && i < 7; i++) {
        if (!queries[i]) continue;
        duckvep_sql_text *sql = &out[n++];
        /* The newlines keep a trailing line comment in a query from swallowing the parenthesis. */
        ok = duckvep_sql_append(sql, "COPY (\n") && duckvep_sql_append(sql, queries[i]) &&
            duckvep_sql_append(sql, "\n) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL ") &&
            duckvep_sql_literal(sql, key.data) && duckvep_sql_append(sql, ", RELATION ") &&
            duckvep_sql_literal(sql, duckvep_model_relations[i]) &&
            duckvep_sql_append(sql, ", USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)");
    }
    if (ok) {
        duckvep_sql_text *sql = &out[n++];
        ok = duckvep_sql_append(sql, "SELECT duckvep_model_publish(") && duckvep_sql_literal(sql, key.data);
        if (ok && (reference_fasta || coverage >= 0)) {
            ok = duckvep_sql_append(sql, ", {");
            if (ok && reference_fasta)
                ok = duckvep_sql_append(sql, "'reference_fasta': ") && duckvep_sql_literal(sql, reference_fasta);
            if (ok && coverage >= 0)
                ok = duckvep_sql_append(sql, reference_fasta ? ", 'transcript_coverage_complete': " :
                                                               "'transcript_coverage_complete': ") &&
                    duckvep_sql_append(sql, coverage ? "true" : "false");
            if (ok) ok = duckvep_sql_append(sql, "}");
        }
        if (ok) ok = duckvep_sql_append(sql, ")");
    }
    duckvep_sql_free(&key);
    *count = n;
    return ok;
}
