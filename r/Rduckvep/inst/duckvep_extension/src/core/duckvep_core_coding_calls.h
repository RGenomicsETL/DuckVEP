/* duckvep_coding_calls, host-neutral: the fused VCF/BCF reader behind the table function of both hosts. A plan binds a
 * pinned model to a path (the bind-time checks); a reader yields one calls row at a time. The views of a row stay valid
 * until the next duckvep_cc_next. See duckvep_core_coding_calls.c. */
#ifndef DUCKVEP_CORE_CODING_CALLS_H
#define DUCKVEP_CORE_CODING_CALLS_H

#include <stddef.h>
#include <stdint.h>

#include "core/duckvep_core_discovery.h"
#include "core/duckvep_core_model.h"

typedef struct duckvep_cc_plan duckvep_cc_plan_t;
typedef struct duckvep_cc_reader duckvep_cc_reader_t;

typedef struct {
    int64_t event_index, position, phase_set;
    int32_t region, alt_index, transcript, sample;
    const char *ref, *alt;
    size_t ref_length, alt_length;
    uint32_t lanes;
    const int32_t *alleles;   /* lanes values, -1 for a missing allele */
    const uint8_t *phase;     /* lanes flags, the first always 0 */
    int phase_set_present;
} duckvep_cc_row_t;

/* The output schema: names and SQL types of the eleven columns. */
#define DUCKVEP_CC_COLUMNS 11u
const char *duckvep_cc_column_name(unsigned index);
const char *duckvep_cc_column_type(unsigned index);

/* NULL with `error` set: a lifted model, a model with no (or duplicate) seq_region_name values, allocation. */
duckvep_cc_plan_t *duckvep_cc_plan_create(const duckvep_owned_model_t *model, const char *path, char *error, size_t error_size);
void duckvep_cc_plan_destroy(duckvep_cc_plan_t *plan);

/* Opens the file: NULL with `error` set. The plan must outlive the reader. */
duckvep_cc_reader_t *duckvep_cc_open(const duckvep_cc_plan_t *plan, char *error, size_t error_size);
void duckvep_cc_close(duckvep_cc_reader_t *reader);
/* 1 a row, 0 end of file, -1 failure (`error` set). */
int duckvep_cc_next(duckvep_cc_reader_t *reader, duckvep_cc_row_t *row, char *error, size_t error_size);

#endif
