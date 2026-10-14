#include "duckvep_product_profile.h"
#include <limits.h>
#include <string.h>

static int dna_base(uint8_t b) {
    return b == 'A' || b == 'C' || b == 'G' || b == 'T';
}

static int inside(uint32_t pos, uint32_t lo, uint32_t hi) {
    return lo <= pos && pos <= hi;
}

static uint32_t position_flags(const duckvep_profile_v1_input_t *in, uint32_t pos) {
    uint32_t flags = 0u;
    for (size_t j = 1u; j < in->exon_count; j++) {
        uint32_t left = in->exons[j - 1u].end1;
        uint32_t right = in->exons[j].start1;
        uint32_t gap_start = left + 1u, gap_end = right - 1u;
        if (inside(pos, gap_start, gap_end)) {
            flags |= DUCKVEP_PROFILE_INTRON;
        }
    }
    return flags;
}

duckvep_profile_status_t duckvep_product_position_v1(
    const duckvep_profile_v1_input_t *in, duckvep_profile_v1_fact_t *out,
    size_t capacity, size_t *count) {
    if (!count) return DUCKVEP_PROFILE_INVALID;
    *count = 0u;
    if (!in || (capacity && !out) || !in->reference_start1 ||
        !in->dna_length || in->dna_length > UINT32_MAX - in->reference_start1 + (uint64_t)1u ||
        !in->reference_dna || !in->reference_witness || !in->final_dna ||
        (in->exon_count && !in->exons) || (in->feature_count && !in->features) ||
        (in->strand != 1 && in->strand != -1) ||
        in->transcript.start1 < in->reference_start1 ||
        in->transcript.end1 < in->transcript.start1 ||
        (uint64_t)in->transcript.end1 >= (uint64_t)in->reference_start1 + in->dna_length)
        return DUCKVEP_PROFILE_INVALID;
    if (in->joint_status != DUCKVEP_PROFILE_OK)
        return in->joint_status == DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION
            ? DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION : DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE;
    if (!in->model_complete || !in->exon_count) return DUCKVEP_PROFILE_PARTIAL_MODEL;
    if (in->witness_length != in->dna_length || in->final_length != in->dna_length)
        return DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE;
    for (size_t j = 0u; j < in->exon_count; j++) {
        duckvep_profile_interval_t e = in->exons[j];
        if (e.start1 < in->transcript.start1 || e.end1 > in->transcript.end1 ||
            e.start1 > e.end1 ||
            (j && in->exons[j - 1u].end1 >= e.start1 - 1u))
            return DUCKVEP_PROFILE_PARTIAL_MODEL;
    }
    for (size_t j = 0u; j < in->feature_count; j++) {
        duckvep_profile_feature_t f = in->features[j];
        if (!f.feature_id || f.interval.start1 < in->reference_start1 ||
            f.interval.end1 < f.interval.start1 ||
            (uint64_t)f.interval.end1 >= (uint64_t)in->reference_start1 + in->dna_length ||
            (j && (f.interval.start1 < in->features[j - 1u].interval.start1 ||
              (f.interval.start1 == in->features[j - 1u].interval.start1 &&
               (f.interval.end1 < in->features[j - 1u].interval.end1 ||
                (f.interval.end1 == in->features[j - 1u].interval.end1 &&
                 f.feature_id <= in->features[j - 1u].feature_id))))))
            return DUCKVEP_PROFILE_INVALID;
    }
    if (memcmp(in->reference_witness, in->reference_dna, in->dna_length))
        return DUCKVEP_PROFILE_REFERENCE_MISMATCH;
    for (size_t i = 0u; i < in->dna_length; i++)
        if (!dna_base(in->reference_dna[i]) || !dna_base(in->final_dna[i]))
            return DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE;
    size_t needed = 0u;
    for (size_t i = 0u; i < in->dna_length; i++) {
        if (in->reference_dna[i] == in->final_dna[i]) continue;
        if (needed == SIZE_MAX) return DUCKVEP_PROFILE_SHORT_BUFFER;
        needed++;
        uint32_t pos = in->reference_start1 + (uint32_t)i;
        for (size_t j = 0u; j < in->feature_count; j++)
            if (inside(pos, in->features[j].interval.start1, in->features[j].interval.end1)) {
                if (needed == SIZE_MAX) return DUCKVEP_PROFILE_SHORT_BUFFER;
                needed++;
            }
    }
    if (needed > capacity) return DUCKVEP_PROFILE_SHORT_BUFFER;
    size_t k = 0u;
    for (size_t i = 0u; i < in->dna_length; i++) {
        if (in->reference_dna[i] == in->final_dna[i]) continue;
        uint32_t pos = in->reference_start1 + (uint32_t)i;
        out[k++] = (duckvep_profile_v1_fact_t){pos, position_flags(in, pos), 0u};
        for (size_t j = 0u; j < in->feature_count; j++)
            if (inside(pos, in->features[j].interval.start1, in->features[j].interval.end1))
                out[k++] = (duckvep_profile_v1_fact_t){pos,
                    DUCKVEP_PROFILE_REGULATORY_INTERVAL, in->features[j].feature_id};
    }
    *count = needed;
    return DUCKVEP_PROFILE_OK;
}
