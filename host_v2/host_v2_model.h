/* The per-database model state shared by the model functions and the annotation natives. */
#ifndef DUCKVEP_HOST_V2_MODEL_H
#define DUCKVEP_HOST_V2_MODEL_H

#include "host_v2_common.h"
#include "core/duckvep_core_model.h"

typedef struct model_state model_state;

void host_v2_state_retain(model_state *state);
void host_v2_state_release(void *state);
duckvep_registry_t *host_v2_registry_of(void *user_data);
bool host_v2_register_annotate(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                               model_state *state, duckdb_v2_error_info_handle *error);

#endif
