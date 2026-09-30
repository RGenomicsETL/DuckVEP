#include "duckvep_core_geometry.h"

#include "kernel/include/duckvep_so.h"
#include "kernel/src/duckvep_event.h"

const char *const duckvep_core_allele_geometry_fields[DUCKVEP_CORE_GEOMETRY_FIELD_COUNT] = {
    "kind_code", "interbase", "anchor_side_code", "raw_start0",
    "raw_end0", "feature_start0", "feature_end0", "edit_start0",
    "edit_end0", "insertion_boundary0", "reference_difference_offset",
    "reference_difference_length", "alternate_difference_offset",
    "alternate_difference_length"
};

const char duckvep_core_allele_geometry_error[] =
    "duckvep_allele_geometry: position must fit UINTEGER and REF/ALT must be distinct non-empty A/C/G/T/N alleles of at most 65,535 bases";

size_t duckvep_core_so_term_count(void)
{
    return (size_t)DUCKVEP_SO_BIT_COUNT;
}

bool duckvep_core_so_term(size_t index, duckvep_core_so_term_t *out)
{
    duckvep_so_bit_t bit;
    duckvep_impact_t impact;
    const char *name;

    if (index >= (size_t)DUCKVEP_SO_BIT_COUNT || out == NULL)
        return false;
    bit = (duckvep_so_bit_t)index;
    impact = duckvep_so_bit_impact(bit);
    name = duckvep_so_name(bit);
    out->bit_index = (uint8_t)bit;
    out->consequence_mask = DUCKVEP_SO(bit);
    out->consequence = name != NULL ? name : "";
    out->impact_code = (uint8_t)impact;
    out->impact = duckvep_impact_name(impact);
    out->severity_rank = duckvep_so_rank(bit);
    out->evaluator_tier = duckvep_so_tier(bit);
    return true;
}

static bool dna_valid(const char *sequence, size_t length)
{
    size_t index;

    if (sequence == NULL || length == 0 || length > UINT16_MAX)
        return false;
    for (index = 0; index < length; index++) {
        unsigned char base = (unsigned char)sequence[index];
        if (base >= 'a' && base <= 'z')
            base = (unsigned char)(base - ('a' - 'A'));
        if (base != 'A' && base != 'C' && base != 'G' &&
            base != 'T' && base != 'N')
            return false;
    }
    return true;
}

bool duckvep_core_allele_geometry(uint64_t position, const char *reference,
    size_t reference_length, const char *alternate, size_t alternate_length,
    duckvep_core_allele_geometry_t *out)
{
    duckvep_event_t event;

    if (position == 0 || position > UINT32_MAX ||
        !dna_valid(reference, reference_length) ||
        !dna_valid(alternate, alternate_length) ||
        !duckvep_event_prepare_small((uint32_t)position,
            (const uint8_t *)reference, (uint16_t)reference_length,
            (const uint8_t *)alternate, (uint16_t)alternate_length, &event))
        return false;
    out->kind_code = event.kind;
    out->interbase = event.interbase != 0u;
    out->anchor_side_code = event.anchor_side;
    out->raw_start0 = (uint64_t)event.raw_start1 - 1u;
    out->raw_end0 = event.raw_end1;
    out->feature_start0 = (uint64_t)event.feature_start1 - 1u;
    out->feature_end0 = event.feature_end1;
    out->edit_start0 = event.interbase ? event.insertion_boundary0
        : (uint64_t)event.start1 - 1u;
    out->edit_end0 = event.interbase ? event.insertion_boundary0
        : (uint64_t)event.end1;
    out->has_insertion_boundary0 = event.interbase != 0u;
    out->insertion_boundary0 = event.interbase ? event.insertion_boundary0 : 0u;
    out->reference_difference_offset = event.ref_diff_offset;
    out->reference_difference_length = event.ref_diff_length;
    out->alternate_difference_offset = event.alt_diff_offset;
    out->alternate_difference_length = event.alt_diff_length;
    return true;
}

const char *duckvep_core_breakend_error(duckvep_breakend_status_t status)
{
    return status == DUCKVEP_BREAKEND_POSITION_OVERFLOW
        ? "duckvep_breakend_geometry: mate position exceeds UBIGINT"
        : "duckvep_breakend_geometry: malformed BND ALT; expected two matching brackets around chrom:position, or a single leading/trailing dot, with non-empty A/C/G/T/N replacement sequence";
}
