#include "duckvep_core_phase.h"

#include "kernel/src/duckvep_haplotype_stream.h"

#include <stdlib.h>
#include <string.h>

const char duckvep_core_phase_policy_error[] =
    "duckvep_phase_call: phase_policy must be 'strict' or 'vep_compat'";
const char duckvep_core_phase_set_error[] = "duckvep_phase_call: phase_set exceeds BIGINT range";
const char duckvep_core_record_order_error[] = "duckvep_haplotypes: invalid source-buffer ordinal";

static const char allele_error[] =
    "duckvep_phase_call: called allele indices must be non-negative INTEGER values";
static const char flag_error[] = "duckvep_phase_call: phase_before elements must be BOOLEAN values";
static const char state_error[] = "duckvep_phase_call: invalid decoded phase state";

const char *duckvep_core_phase_check_row(bool have_gt, bool have_phase, size_t gt_length,
    size_t phase_length, size_t *total)
{
    if ((!have_gt && have_phase) || (have_gt && have_phase && gt_length != phase_length))
        return "duckvep_phase_call: allele and phase lists must have equal length";
    if (!have_gt) return NULL;
    if (!gt_length || gt_length > UINT16_MAX || gt_length > UINT64_MAX - *total)
        return "duckvep_phase_call: ploidy must be between 1 and 65535";
    *total += gt_length;
    return NULL;
}

bool duckvep_core_phase_policy(const char *name, size_t length, duckvep_phase_policy_t *policy)
{
    *policy = DUCKVEP_PHASE_STRICT;
    /* "vep_compat" names the pinned executable VEP release; "vep116_compat" is its older spelling. */
    if ((length == 10u && !memcmp(name, "vep_compat", 10u)) ||
        (length == 13u && !memcmp(name, "vep116_compat", 13u))) {
        *policy = DUCKVEP_PHASE_VEP116_COMPAT;
        return true;
    }
    return length == 6u && !memcmp(name, "strict", 6u);
}

bool duckvep_core_phase_set(const duckvep_cell_t *cell, int64_t *phase_set)
{
    *phase_set = 0;
    switch (cell->kind) {
    case DUCKVEP_CELL_TINYINT: case DUCKVEP_CELL_SMALLINT:
    case DUCKVEP_CELL_INTEGER: case DUCKVEP_CELL_BIGINT:
        *phase_set = cell->i; break;
    case DUCKVEP_CELL_UTINYINT: case DUCKVEP_CELL_USMALLINT: case DUCKVEP_CELL_UINTEGER:
        *phase_set = (int64_t)cell->u; break;
    case DUCKVEP_CELL_UBIGINT:
        if (cell->u > INT64_MAX) return false;
        *phase_set = (int64_t)cell->u;
        break;
    default: break;
    }
    return true;
}

/* Reads one slot; returns an error message or NULL. */
static const char *read_slot(const duckvep_core_phase_reader_t *reader, size_t slot,
    bool have_phase, bool *called, int32_t *allele, bool *phased)
{
    duckvep_cell_t cell = {0};
    *called = false;
    *allele = -1;
    *phased = false;
    reader->allele(reader->context, slot, &cell);
    *called = cell.valid;
    if (*called && !duckvep_core_phase_allele(&cell, allele)) return allele_error;
    if (have_phase) {
        duckvep_cell_t flag = {0};
        reader->phase(reader->context, slot, &flag);
        if (flag.valid && !duckvep_core_phase_flag(&flag, phased)) return flag_error;
    }
    return NULL;
}

duckvep_core_phase_result_t duckvep_core_phase_row(const duckvep_core_phase_reader_t *reader,
    size_t count, bool have_phase, duckvep_phase_policy_t policy, const char **error)
{
    duckvep_phase_summary_t summary = {0};
    uint16_t called_before = 0u;
    bool called, phased;
    int32_t allele;

    for (size_t slot = 0u; slot < count; slot++) {
        *error = read_slot(reader, slot, have_phase, &called, &allele, &phased);
        if (*error) return DUCKVEP_CORE_PHASE_ERROR;
        if ((called && allele < 0) ||
            duckvep_phase_observe(&summary, allele, phased) != DUCKVEP_PHASE_OK) {
            *error = allele_error;
            return DUCKVEP_CORE_PHASE_ERROR;
        }
    }
    for (size_t slot = 0u; slot < count; slot++) {
        duckvep_phase_assignment_t assignment;
        duckvep_core_phase_slot_t out;
        *error = read_slot(reader, slot, have_phase, &called, &allele, &phased);
        if (*error) return DUCKVEP_CORE_PHASE_ERROR;
        if (duckvep_phase_assign(&summary, (uint16_t)(slot + 1u), called_before, allele, phased,
                                 policy, &assignment) != DUCKVEP_PHASE_OK) {
            *error = state_error;
            return DUCKVEP_CORE_PHASE_ERROR;
        }
        if (called) called_before++;
        out.input_slot = (uint16_t)(slot + 1u);
        out.allele_called = called;
        out.allele_index = allele;
        out.lane = assignment.lane;
        out.ploidy = summary.ploidy;
        out.phase_set_applies = assignment.scope == DUCKVEP_PHASE_SET;
        out.scope = assignment.scope == DUCKVEP_PHASE_SET ? "phase_set" :
            assignment.scope == DUCKVEP_PHASE_ALL_SETS ? "all_phase_sets" :
            assignment.scope == DUCKVEP_PHASE_ALLELE_SLOT ? "allele_slot" : "unresolved";
        out.status = assignment.status == DUCKVEP_PHASE_CALLED ? "called" :
            assignment.status == DUCKVEP_PHASE_MISSING ? "missing" : "unphased";
        if (!reader->emit(reader->context, slot, &out)) return DUCKVEP_CORE_PHASE_HOST_FAILED;
    }
    return DUCKVEP_CORE_PHASE_OK;
}

void duckvep_core_revcomp(const char *sequence, size_t length, char *reversed)
{
    const char *source = "ACGTRYSWKMBDHVNacgtryswkmbdhvn";
    const char *target = "TGCAYRSWMKVHDBNtgcayrswmkvhdbn";
    size_t end = length, at = 0;
    while (end) {
        size_t first = end - 1;
        while (first && ((unsigned char)sequence[first] & 0xc0) == 0x80)
            first--;
        size_t width = end - first;
        if (width == 1) {
            char base = sequence[first];
            const char *match = base ? strchr(source, base) : NULL;
            reversed[at] = match ? target[match - source] : base;
        } else {
            memcpy(reversed + at, sequence + first, width);
        }
        at += width;
        end = first;
    }
}

void duckvep_core_raw_gt(const char *gt, size_t length, uint32_t source_alt_count, uint32_t fields[7])
{
    duckvep_raw_gt_t call = {0};
    duckvep_raw_gt_status_t status = duckvep_phase_parse_vep116_raw(
        (const uint8_t *)gt, length, source_alt_count, &call);
    fields[0] = (uint32_t)status;
    fields[1] = call.allele_index[0]; fields[2] = call.allele_index[1];
    fields[3] = call.parsed_slots; fields[4] = call.source_ploidy;
    fields[5] = call.source_has_missing; fields[6] = (uint32_t)call.disposition;
}

uint64_t duckvep_core_record_order(uint64_t count, uint64_t ordinal)
{
    return duckvep_haplotype_record_order(count, ordinal);
}
