/* duckvep_lof_sql: the LOFTEE loss-of-function relation as SQL over rows DuckVEP already
 * produces. The text and the option rules are host-neutral; each host reads the arguments. */
#ifndef DUCKVEP_CORE_LOF_H
#define DUCKVEP_CORE_LOF_H

#include "duckvep_core_cells.h"
#include "duckvep_core_sql.h"

enum { DUCKVEP_LOF_GERP, DUCKVEP_LOF_ANCESTOR, DUCKVEP_LOF_PHYLOCSF, DUCKVEP_LOF_MIN_INTRON,
       DUCKVEP_LOF_CUTOFF, DUCKVEP_LOF_CHECK_CDS, DUCKVEP_LOF_OPTIONS };

typedef struct {
    char *relation[3];       /* gerp, ancestor, phylocsf; NULL when omitted */
    int64_t min_intron_size;
    double gerp_cutoff;
    bool check_complete_cds;
} duckvep_lof_options_t;

extern const char *const duckvep_core_lof_option_names[DUCKVEP_LOF_OPTIONS];
extern const duckvep_core_option_kind_t duckvep_core_lof_option_kinds[DUCKVEP_LOF_OPTIONS];
extern const char duckvep_core_lof_names_failed[];

/* Reads the options from cells (NULL: absent). Returns an error message or NULL. */
const char *duckvep_core_lof_read_options(const duckvep_cell_t *const cells[DUCKVEP_LOF_OPTIONS],
    duckvep_lof_options_t *out);
void duckvep_core_lof_free_options(duckvep_lof_options_t *options);

/* names: annotations, transcripts, reference relation names. */
bool duckvep_core_lof_sql(char *const names[3], const duckvep_lof_options_t *options,
    duckvep_sql_text *out);

#endif
