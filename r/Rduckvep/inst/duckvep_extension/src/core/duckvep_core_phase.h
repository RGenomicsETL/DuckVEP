/* Host-neutral row logic of duckvep_phase_call, plus the small text kernels
 * (_duckvep_revcomp, _duckvep_raw_gt, _duckvep_record_order). */
#ifndef DUCKVEP_CORE_PHASE_H
#define DUCKVEP_CORE_PHASE_H

#include "duckvep_core_cells.h"
#include "kernel/src/duckvep_phase.h"

/* One output slot of a phased call. */
typedef struct {
    uint16_t input_slot;
    bool allele_called;
    int32_t allele_index;
    uint16_t lane;           /* 0: no lane */
    uint16_t ploidy;
    bool phase_set_applies;  /* the slot's scope is a phase set */
    const char *scope;
    const char *status;
} duckvep_core_phase_slot_t;

/* How a host supplies one row. `allele` and `phase` fill the cell of a slot
 * (valid=false for a NULL element, and for every slot of an untyped-NULL list);
 * `phase` is only called when the phase list is present. `emit` receives the
 * finished slots in order and returns false when the host failed (and has
 * already reported why). */
typedef struct {
    void *context;
    void (*allele)(void *context, size_t slot, duckvep_cell_t *cell);
    void (*phase)(void *context, size_t slot, duckvep_cell_t *cell);
    bool (*emit)(void *context, size_t slot, const duckvep_core_phase_slot_t *out);
} duckvep_core_phase_reader_t;

typedef enum { DUCKVEP_CORE_PHASE_OK, DUCKVEP_CORE_PHASE_ERROR, DUCKVEP_CORE_PHASE_HOST_FAILED }
    duckvep_core_phase_result_t;

extern const char duckvep_core_phase_policy_error[];
extern const char duckvep_core_phase_set_error[];

/* Row-level list checks: equal lengths, ploidy 1..65535. Adds the row's slots to
 * *total. Returns an error message or NULL. */
const char *duckvep_core_phase_check_row(bool have_gt, bool have_phase, size_t gt_length,
    size_t phase_length, size_t *total);

bool duckvep_core_phase_policy(const char *name, size_t length, duckvep_phase_policy_t *policy);

/* The phase_set option as a BIGINT; false when an unsigned value does not fit. */
bool duckvep_core_phase_set(const duckvep_cell_t *cell, int64_t *phase_set);

duckvep_core_phase_result_t duckvep_core_phase_row(const duckvep_core_phase_reader_t *reader,
    size_t count, bool have_phase, duckvep_phase_policy_t policy, const char **error);

/* _duckvep_revcomp: writes length bytes (reverse complement of IUPAC bases;
 * UTF-8 sequences and unknown bytes are kept). */
void duckvep_core_revcomp(const char *sequence, size_t length, char *out);

/* _duckvep_raw_gt: the seven UINTEGER fields. */
void duckvep_core_raw_gt(const char *gt, size_t length, uint32_t source_alt_count, uint32_t fields[7]);

/* _duckvep_record_order: 0 means an invalid ordinal (see the message). */
uint64_t duckvep_core_record_order(uint64_t count, uint64_t ordinal);
extern const char duckvep_core_record_order_error[];

#endif
