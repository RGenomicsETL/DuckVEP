#include "duckvep_core_model_script.h"

const char *const duckvep_model_relations[6] = {
    "regions", "transcripts", "exons", "mature_mirna", "peptide_edits", "interval_features"};

bool duckvep_core_model_script(const char *name, const char *const queries[6],
    const char *reference_fasta, int coverage, duckvep_sql_text out[DUCKVEP_MODEL_SCRIPT_MAX],
    size_t *count)
{
    bool ok = true;
    size_t n = 0;

    for (size_t i = 0; ok && i < 6; i++) {
        if (!queries[i]) continue;
        duckvep_sql_text *sql = &out[n++];
        /* The newlines keep a trailing line comment in a query from swallowing the parenthesis. */
        ok = duckvep_sql_append(sql, "COPY (\n") && duckvep_sql_append(sql, queries[i]) &&
            duckvep_sql_append(sql, "\n) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL ") &&
            duckvep_sql_literal(sql, name) && duckvep_sql_append(sql, ", RELATION ") &&
            duckvep_sql_literal(sql, duckvep_model_relations[i]) &&
            duckvep_sql_append(sql, ", USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)");
    }
    if (ok) {
        duckvep_sql_text *sql = &out[n++];
        ok = duckvep_sql_append(sql, "SELECT duckvep_model_publish(") && duckvep_sql_literal(sql, name);
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
    *count = n;
    return ok;
}
