/* The phased-haplotype replay of duckvep_haplotypes, host-neutral: the bounded workspace, the stream kernel driver, the
 * per-leaf result writers and the output schema. A host supplies the model, the options, a row-at-a-time input
 * (duckvep_hap_input_t over its own query result or staged rows) and an output chunk written through the duckvep_h_*
 * host layer of "duckvep_host.h". The core never mentions DuckDB. */
#ifndef DUCKVEP_CORE_HAPLOTYPES_H
#define DUCKVEP_CORE_HAPLOTYPES_H

#include <stddef.h>
#include <stdint.h>

#include "core/duckvep_core_model.h"
#include "kernel/src/duckvep_haplotype_stream.h"

enum { DUCKVEP_HAP_LIMIT_EVENTS, DUCKVEP_HAP_LIMIT_TRANSCRIPTS, DUCKVEP_HAP_LIMIT_CARRIERS,
    DUCKVEP_HAP_LIMIT_PREFIXES, DUCKVEP_HAP_LIMIT_PROJECTIONS, DUCKVEP_HAP_LIMIT_ALLELES,
    DUCKVEP_HAP_LIMIT_LEAF_EVENTS, DUCKVEP_HAP_LIMIT_LEAF_EDITS, DUCKVEP_HAP_LIMIT_SEQUENCE,
    DUCKVEP_HAP_LIMIT_PLOIDY, DUCKVEP_HAP_LIMIT_PHASE_SETS, DUCKVEP_HAP_LIMIT_ALIGNMENT,
    DUCKVEP_HAP_LIMIT_DIFFERENCES, DUCKVEP_HAP_LIMIT_HGVS_OPERATIONS, DUCKVEP_HAP_LIMIT_HGVS_BYTES,
    DUCKVEP_HAP_LIMIT_HGVS_REFERENCE, DUCKVEP_HAP_LIMIT_WORKSPACE, DUCKVEP_HAP_LIMIT_COUNT };

extern const char *const duckvep_hap_limit_names[DUCKVEP_HAP_LIMIT_COUNT];
extern const uint64_t duckvep_hap_limit_defaults[DUCKVEP_HAP_LIMIT_COUNT];
/* Whether `value` is an acceptable setting of limit `index` (a named option of the function). */
int duckvep_hap_limit_valid(unsigned index, uint64_t value);

#define DUCKVEP_HAP_OUTPUT_COLUMNS 32u
/* The output column `index` (0..30): its name and its type as SQL text; `source_records` selects the eight-field
 * contributors record. */
const char *duckvep_hap_column_name(unsigned index);
const char *duckvep_hap_column_type(unsigned index, int source_records);

typedef struct {
    const duckvep_owned_model_t *model;
    duckvep_phase_policy_t policy;
    int source_records, hgvs;
    size_t limits[DUCKVEP_HAP_LIMIT_COUNT];
} duckvep_hap_config_t;

/* A config with every option at its default. */
void duckvep_hap_config_defaults(duckvep_hap_config_t *config, const duckvep_owned_model_t *model);
/* Checks the option combination; false with `error` set. */
int duckvep_hap_config_check(const duckvep_hap_config_t *config, char *error, size_t error_size);

/* One normalized input row. Its views are borrowed from the host and stay valid until the host's next fetch.
 * Columns, in order: alt_events 0 event_index, 1 seq_region, 2 position, 3 reference, 4 alternate, 5 alt_index,
 * 6 transcript_index, 7 sample_index, 8 alleles, 9 phase_before, 10 phase_set, 11 domain_sets, 12 copies,
 * 13 versions, 14 ploidies; source_records 0..7 as above, 8 raw_gt, 9/10 NULL, 11 alt_count, 12 copies, 13 versions,
 * 14 gt_versions, 15 replay_order, 16 source_selected. `null_mask` has bit i set when column i is NULL. */
typedef struct {
    uint32_t null_mask;
    uint64_t event_id, pos, replay_order;
    uint32_t chrom, allele_index, transcript, sample;
    const uint8_t *ref, *alt;
    uint32_t ref_len, alt_len;
    int64_t copies, versions, ploidies;
    /* alt_events: the genotype lists; a view is (values, validity, first element, length). */
    const int32_t *gt; const uint64_t *gt_validity; size_t gt_offset, gt_length;
    int have_phase;
    const uint8_t *phase; const uint64_t *phase_validity; size_t phase_offset, phase_length;
    int phase_set_present; int64_t phase_set;
    const int64_t *sets; const uint64_t *sets_validity; size_t sets_offset, sets_length;
    /* source_records: the parsed raw GT struct (seven UINTEGER fields) and the selection flag. */
    uint32_t raw[7];
    int source_selected;
} duckvep_hap_row_t;

typedef struct {
    void *context;
    /* Fills `row` with the next input row: 1 a row, 0 end of input, -1 failure (`error` set). */
    int (*next)(void *context, duckvep_hap_row_t *row, char *error, size_t error_size);
} duckvep_hap_input_t;

typedef struct duckvep_hap_state duckvep_hap_state_t;

/* Allocates the bounded workspace and initializes the stream; NULL with `error` set (config limits exceeded, invalid
 * model). The state keeps a copy of the config; the model must outlive it. */
duckvep_hap_state_t *duckvep_hap_open(const duckvep_hap_config_t *config, char *error, size_t error_size);
void duckvep_hap_close(duckvep_hap_state_t *state);
/* Bytes of the bounded workspace in use (the calls query text is charged against the same limit). */
size_t duckvep_hap_workspace_bytes(const duckvep_hap_state_t *state);

/* Writes up to `capacity` result rows to the host chunk `output` (opaque here: duckvep_h_chunk), pulling input as it
 * needs it. Returns 1 with *rows the number written (0 at the end of the stream), 0 with `error` set. */
int duckvep_hap_scan(duckvep_hap_state_t *state, const duckvep_hap_input_t *input, void *output, size_t capacity,
    size_t *rows, char *error, size_t error_size);

#endif
