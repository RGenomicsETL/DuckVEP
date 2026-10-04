/* Model snapshots: save the arrays of a loaded model to a file, and restore a model by mapping it. */
#ifndef DUCKVEP_CORE_SNAPSHOT_H
#define DUCKVEP_CORE_SNAPSHOT_H

#include "duckvep_core_model.h"

/* Writes `model` to `path` (through `path`.partial and a rename). Returns 1, or 0 with `error` set. */
int duckvep_core_model_snapshot_save(const duckvep_owned_model_t *model, const char *path,
	char *error, size_t error_size);

/* Maps the snapshot at `path`, validates it and installs it under `name`. The error is final
 * (budget refusals named). */
int duckvep_core_model_snapshot_install(duckvep_registry_t *registry, const char *name, const char *path,
	char *final_error, size_t final_error_size);

/* Unmaps (or frees) the storage of a restored model and returns its bytes to the budget. */
void duckvep_core_snapshot_release(void *base, size_t bytes, int mapped);

#endif
