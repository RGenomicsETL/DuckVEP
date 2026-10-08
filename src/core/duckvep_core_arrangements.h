/* Bounded hypothetical diploid-arrangement replay shared by both hosts. */
#ifndef DUCKVEP_CORE_ARRANGEMENTS_H
#define DUCKVEP_CORE_ARRANGEMENTS_H

#include "core/duckvep_core_haplotypes.h"

enum { DUCKVEP_ARRANGEMENT_DEFAULT_SITES = 16u, DUCKVEP_ARRANGEMENT_DEFAULT_CALLS = 64u,
    DUCKVEP_ARRANGEMENT_DEFAULT_ALTERNATIVES = 128u, DUCKVEP_ARRANGEMENT_DEFAULT_REPLAYS = 256u,
    DUCKVEP_ARRANGEMENT_LIMIT_COUNT = 4u, DUCKVEP_ARRANGEMENT_COLUMN_COUNT = 26u };

extern const char *const duckvep_arrangement_limit_names[DUCKVEP_ARRANGEMENT_LIMIT_COUNT];

typedef struct {
    const duckvep_owned_model_t *model;
    size_t max_sites, max_calls, max_arrangements, max_replays;
} duckvep_arrangement_config_t;

typedef struct {
    uint64_t hypothesis_id, event_id;
    uint32_t transcript, sample, chrom, allele_index;
    uint64_t position;
    const uint8_t *reference, *alternate, *cds, *protein;
    size_t reference_length, alternate_length, cds_length, protein_length;
    int32_t allele0, allele1, assigned_allele;
    uint8_t phase0, phase1, phase_present, phase_set_present, contributes, reference_lane;
    int64_t phase_set, nominal_length_diff;
    uint64_t consequence_mask;
    uint16_t lane;
    uint8_t prediction_status, prediction_reason, hypothetical;
} duckvep_arrangement_row_t;

typedef struct duckvep_arrangement_state duckvep_arrangement_state_t;

void duckvep_arrangement_config_defaults(duckvep_arrangement_config_t *config,
    const duckvep_owned_model_t *model);
int duckvep_arrangement_config_check(const duckvep_arrangement_config_t *config,
    char *error, size_t error_size);
duckvep_arrangement_state_t *duckvep_arrangement_open(const duckvep_arrangement_config_t *config,
    char *error, size_t error_size);
void duckvep_arrangement_close(duckvep_arrangement_state_t *state);
/* Loads, validates, enumerates and replays all hypotheses before returning. */
int duckvep_arrangement_load(duckvep_arrangement_state_t *state, const duckvep_hap_input_t *input,
    char *error, size_t error_size);
/* Rows are hypothesis-major, then transcript, lane and original source call. */
int duckvep_arrangement_next(duckvep_arrangement_state_t *state, duckvep_arrangement_row_t *row);
size_t duckvep_arrangement_count(const duckvep_arrangement_state_t *state);
size_t duckvep_arrangement_workspace_bytes(const duckvep_arrangement_state_t *state);
const char *duckvep_arrangement_prediction_status(uint8_t status);

#endif
