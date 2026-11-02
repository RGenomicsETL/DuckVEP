/* The staging area of the v2 host, shared by the model sink (host_v2_model.c) and the haplotype job sink and scans
 * (host_v2_haplotypes.c): typed caller rows captured by COPY into DuckDB-managed column collections (they spill under the
 * memory limit), keyed by (literal model name, optional load key, relation). Haplotype jobs remain unkeyed. */
#ifndef DUCKVEP_HOST_V2_STAGE_H
#define DUCKVEP_HOST_V2_STAGE_H

#include "host_v2_model.h"

#include "core/duckvep_core_cells.h"
#include "core/duckvep_core_arrangements.h"
#include "core/duckvep_core_haplotypes.h"

#include <pthread.h>

#define MAX_STAGE_COLUMNS 20

/* The relations of a haplotype job. */
#define HAP_RELATION_CALLS "haplotype_calls"
#define HAP_RELATION_PLAN "haplotype_plan_input"
#define HAP_RELATION_ARRANGEMENTS "haplotype_arrangements"

typedef struct stage {
    char *model;    /* literal model name, or an unkeyed job name */
    char *load_key; /* NULL for unkeyed manual stages and haplotype jobs */
    char *relation;
    duckdb_v2_column_data_collection_handle collection; /* NULL: no rows */
    size_t columns;
    duckvep_ctype_t ctypes[MAX_STAGE_COLUMNS];
    char *names[MAX_STAGE_COLUMNS];
    duckdb_v2_logical_type_handle types[MAX_STAGE_COLUMNS];
    uint64_t rows;
    /* A haplotype job's settings (the options of its COPY). */
    char *hap_model;
    duckvep_hap_config_t hap;
    duckvep_arrangement_config_t arrangements;
    int refs; /* the registry, plus each scan that shares the stage */
    struct stage *next;
} stage;

struct model_state {
    pthread_mutex_t lock;
    stage *stages;
    duckvep_registry_t *registry;
    size_t references;
};

/* Drops one reference; the last one frees the stage and its rows. */
void host_v2_stage_destroy(stage *s);
/* A new reference to the unkeyed stage of (model, relation) that leaves it staged, or NULL. */
stage *host_v2_stage_acquire(model_state *state, const char *model, const char *relation);
/* Removes and returns the unkeyed stage of (model, relation), or NULL. */
stage *host_v2_stage_take(model_state *state, const char *model, const char *relation);
/* Removes and returns the stage of (model, load_key, relation), or NULL. */
stage *host_v2_stage_take_key(model_state *state, const char *model, const char *load_key,
                              const char *relation);
/* Remove all non-job stages belonging to exactly (model, load_key). */
void host_v2_stage_discard_key(model_state *state, const char *model, const char *load_key);
/* Remove only unkeyed, non-job stages belonging to model. */
void host_v2_stage_discard_unkeyed_model(model_state *state, const char *model);
/* Copies the settings of a staged job without removing it; false when the relation is not staged. The model name is a
 * heap copy (free). */
bool host_v2_stage_peek_job(model_state *state, const char *job, const char *relation, char **model,
                            duckvep_hap_config_t *config);

/* COPY options of a haplotype job (implemented in host_v2_haplotypes.c). */
typedef struct hap_bind hap_bind;
hap_bind *host_v2_hap_bind_create(void);
void host_v2_hap_bind_destroy(hap_bind *bind);
/* Consumes the option `name` when it is a job option; `*consumed` tells. */
bool host_v2_hap_bind_option(hap_bind *bind, const duckdb_v2_identifier_t *name, duckdb_v2_value_handle value,
                             bool *consumed, duckdb_v2_error_info_handle *error);
/* True when the COPY is a haplotype job (a JOB option was given). */
bool host_v2_hap_bind_is_job(const hap_bind *bind);
/* Validates the job against the staged layout and fills the stage's key, relation and settings. */
bool host_v2_hap_bind_finish(hap_bind *bind, model_state *state, duckdb_v2_context_handle context, stage *layout,
                             const char *model, duckdb_v2_error_info_handle *error);

#endif
