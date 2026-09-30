/* Text of duckvep_annotate_sql, duckvep_annotate_projected_sql and
 * duckvep_transcript_projection_sql. */
#ifndef DUCKVEP_CORE_ANNOTATE_H
#define DUCKVEP_CORE_ANNOTATE_H

#include "duckvep_core_cells.h"
#include "duckvep_core_sql.h"

extern const char duckvep_core_projected_model_empty[];
extern const char duckvep_core_annotate_failed[];
extern const char duckvep_core_projection_failed[];

/* The options of duckvep_annotate_sql are hgvs, upstream_distance,
 * downstream_distance, rich (index 0..3); the projected builder has only
 * upstream_distance and downstream_distance (index 0..1). A NULL cell pointer
 * means the option is absent (the default); an invalid cell renders as NULL.
 * Returns false on an allocation failure or an unrenderable option. */
bool duckvep_core_annotate_sql(bool projected, const char *events, const char *model,
    const duckvep_cell_t *const *options, duckvep_sql_text *out);

/* events, annotations and transcripts relation names. */
bool duckvep_core_projection_sql(const char *const names[3], duckvep_sql_text *out);

#endif
