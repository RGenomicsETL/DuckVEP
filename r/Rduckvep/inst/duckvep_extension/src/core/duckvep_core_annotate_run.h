/* The annotation natives' entry points. Each takes the host layer's info, input chunk and
 * output vector (see duckvep_host.h) and is registered by the host under its own name. */
#ifndef DUCKVEP_CORE_ANNOTATE_RUN_H
#define DUCKVEP_CORE_ANNOTATE_RUN_H

#include "duckvep_host.h"

typedef enum duckvep_scalar_event_family {
	DUCKVEP_SCALAR_SMALL = 0,
	DUCKVEP_SCALAR_STRUCTURAL,
	DUCKVEP_SCALAR_BREAKEND
} duckvep_scalar_event_family_t;

void duckvep_annotate_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_compact_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_hgvs_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_rich_hgvs_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_sv_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_sv_compact_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_breakend_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_breakend_compact_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_projected_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);
void duckvep_annotate_projected_hgvs_scalar(duckvep_h_info info, duckvep_h_chunk input, duckvep_h_vector output);

#endif
