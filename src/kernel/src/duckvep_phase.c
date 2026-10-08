#include "duckvep_phase.h"

#include <limits.h>
#include <stddef.h>
#include <string.h>

duckvep_raw_gt_status_t duckvep_phase_parse_vep_raw(
    const uint8_t *gt, size_t length, uint32_t source_alt_count, duckvep_raw_gt_t *out) {
    if (!out) return DUCKVEP_RAW_GT_INVALID_ARG;
    memset(out, 0, sizeof(*out));
    if (!gt || !length || source_alt_count > INT32_MAX) return DUCKVEP_RAW_GT_INVALID_ARG;
    duckvep_raw_gt_t result = {0};
    result.allele_index[0] = result.allele_index[1] = UINT32_MAX;
    size_t at = 0u;
    int has_pipe = 0, hom_ref = 1;
    if (gt[0] == '|' || gt[0] == '/') {
        has_pipe = gt[0] == '|';
        hom_ref = 0; /* The upstream non-reference filter does not accept a prefix. */
        at++;
    }
    /* Validate the source grammar and its complete ordinal domain before
     * reproducing the parser's single-separator split and numeric coercion. */
    for (;;) {
        if (at == length) return DUCKVEP_RAW_GT_INVALID_SYNTAX;
        if (result.source_ploidy == UINT16_MAX) return DUCKVEP_RAW_GT_PLOIDY_LIMIT;
        result.source_ploidy++;
        if (gt[at] == '.') {
            result.source_has_missing = 1u;
            hom_ref = 0;
            at++;
        } else {
            if (gt[at] < '0' || gt[at] > '9') return DUCKVEP_RAW_GT_INVALID_SYNTAX;
            uint32_t allele = 0u;
            while (at < length && gt[at] >= '0' && gt[at] <= '9') {
                uint32_t digit = gt[at++] - '0';
                if (allele > (UINT32_MAX - digit) / 10u)
                    return DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE;
                allele = allele * 10u + digit;
                if (allele > source_alt_count) return DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE;
            }
            if (allele) hom_ref = 0;
        }
        if (at == length) break;
        if (gt[at] != '|' && gt[at] != '/') return DUCKVEP_RAW_GT_INVALID_SYNTAX;
        has_pipe |= gt[at++] == '|';
    }
    if (hom_ref) {
        result.disposition = DUCKVEP_RAW_GT_OMITTED_REFERENCE;
        *out = result;
        return DUCKVEP_RAW_GT_OK;
    }
    uint8_t separator = has_pipe ? '|' : '/';
    at = 0u;
    while (at < length) {
        size_t end = at;
        while (end < length && gt[end] != separator) end++;
        /* Only an exact dot is omitted. Empty or non-numeric-leading tokens
         * coerce to REF; e.g. "|0|1", "/1|0", or "./1|2". All original
         * numeric atoms were checked above, so this prefix conversion cannot
         * overflow or address an undeclared source allele. */
        if (!(end - at == 1u && gt[at] == '.')) {
            uint32_t allele = 0u;
            for (size_t i = at; i < end && gt[i] >= '0' && gt[i] <= '9'; i++)
                allele = allele * 10u + (uint32_t)(gt[i] - '0');
            if (result.parsed_slots < 2u) result.allele_index[result.parsed_slots] = allele;
            result.parsed_slots++;
        }
        if (end == length) break;
        at = end + 1u;
    }
    result.disposition = result.parsed_slots ? DUCKVEP_RAW_GT_RETAINED : DUCKVEP_RAW_GT_OMITTED_EMPTY;
    *out = result;
    return DUCKVEP_RAW_GT_OK;
}

duckvep_phase_status_t duckvep_phase_observe(
    duckvep_phase_summary_t *s, int32_t allele, uint8_t phase_before) {
    if (!s || allele < -1 || phase_before > 1u) return DUCKVEP_PHASE_INVALID_ARG;
    if (s->ploidy == UINT16_MAX) return DUCKVEP_PHASE_PLOIDY_LIMIT;
    if (!s->ploidy) {
        s->first_allele = allele;
        s->homozygous = allele >= 0;
    } else if (allele < 0 || allele != s->first_allele) {
        s->homozygous = 0u;
    }
    if (!phase_before) {
        if (!s->unphased_count) {
            s->unphased_allele = allele;
            s->unphased_equal = allele >= 0;
        } else if (allele < 0 || allele != s->unphased_allele) {
            s->unphased_equal = 0u;
        }
        s->unphased_count++;
    }
    s->ploidy++;
    return DUCKVEP_PHASE_OK;
}

duckvep_phase_status_t duckvep_phase_assign(
    const duckvep_phase_summary_t *s, uint16_t slot1, uint16_t called_before,
    int32_t allele, uint8_t phase_before, duckvep_phase_policy_t policy,
    duckvep_phase_assignment_t *out) {
    if (!out) return DUCKVEP_PHASE_INVALID_ARG;
    memset(out, 0, sizeof(*out));
    if (!s || !s->ploidy || !slot1 || slot1 > s->ploidy || called_before >= slot1 ||
        allele < -1 || phase_before > 1u ||
        (policy != DUCKVEP_PHASE_STRICT && policy != DUCKVEP_PHASE_VEP_COMPAT))
        return DUCKVEP_PHASE_INVALID_ARG;
    out->lane = slot1;
    out->status = allele < 0 ? DUCKVEP_PHASE_MISSING : DUCKVEP_PHASE_CALLED;
    if (policy == DUCKVEP_PHASE_VEP_COMPAT) {
        /* BaseVCF4::get_samples_genotypes omits missing entries before
         * Transcript::get_genotypes assigns allele slots and phased=1. */
        out->lane = allele < 0 ? 0u : (uint16_t)(called_before + 1u);
        out->scope = DUCKVEP_PHASE_ALLELE_SLOT;
    } else if (s->ploidy == 1u || s->homozygous) {
        out->scope = DUCKVEP_PHASE_ALL_SETS;
    } else if (phase_before || s->unphased_count == 1u || s->unphased_equal) {
        /* Permuting the remaining unphased slots cannot change their allele
         * assignment if there is only one, or all contain the same known allele. */
        out->scope = DUCKVEP_PHASE_SET;
    } else {
        out->lane = 0u;
        out->scope = DUCKVEP_PHASE_UNRESOLVED;
        if (allele >= 0) out->status = DUCKVEP_PHASE_UNPHASED;
    }
    return DUCKVEP_PHASE_OK;
}

static int arrangement_same_phase_block(const duckvep_phase_arrangement_site_t *left,
                                        const duckvep_phase_arrangement_site_t *right) {
    return left->phase_set_present == right->phase_set_present &&
        (!left->phase_set_present || left->phase_set == right->phase_set);
}

static int arrangement_is_unit_start(const duckvep_phase_arrangement_site_t *sites,
                                     size_t index) {
    if (!sites[index].phase_before[1]) return 1;
    for (size_t prior = 0u; prior < index; prior++) {
        if (sites[prior].phase_before[1] &&
            arrangement_same_phase_block(&sites[prior], &sites[index])) return 0;
    }
    return 1;
}

static size_t arrangement_unit_index(const duckvep_phase_arrangement_site_t *sites,
                                     size_t index) {
    size_t unit = 0u;
    for (size_t site = 0u; site <= index; site++) {
        if (!arrangement_is_unit_start(sites, site)) continue;
        if (site == index || (sites[index].phase_before[1] && sites[site].phase_before[1] &&
            arrangement_same_phase_block(&sites[site], &sites[index]))) return unit;
        unit++;
    }
    return SIZE_MAX;
}

duckvep_phase_arrangement_status_t duckvep_phase_arrange_diploid(
    const duckvep_phase_arrangement_site_t *sites, size_t site_count,
    size_t alternative_limit, uint16_t *alt_lanes, size_t out_capacity,
    size_t *required_alternatives) {
    size_t unit_count = 0u, alternative_count, required_entries;

    if (required_alternatives) *required_alternatives = 0u;
    if (!sites || !site_count || site_count > DUCKVEP_PHASE_ARRANGEMENT_MAX_SITES ||
        !alternative_limit || !required_alternatives || (!alt_lanes && out_capacity))
        return DUCKVEP_PHASE_ARRANGEMENT_INVALID_ARG;
    for (size_t site = 0u; site < site_count; site++) {
        const duckvep_phase_arrangement_site_t *current = &sites[site];
        if (site && current->sample_index != sites[0].sample_index)
            return DUCKVEP_PHASE_ARRANGEMENT_MIXED_SAMPLE;
        if (!current->source_id || current->phase_before[0] > 1u ||
            current->phase_before[1] > 1u || current->phase_set_present > 1u ||
            current->phase_before[0])
            return DUCKVEP_PHASE_ARRANGEMENT_INVALID_ARG;
        if (current->allele[0] < 0 || current->allele[1] < 0)
            return DUCKVEP_PHASE_ARRANGEMENT_MISSING;
        if (current->source_alt_count != 1u || current->allele[0] > 1 ||
            current->allele[1] > 1)
            return DUCKVEP_PHASE_ARRANGEMENT_MULTIALLELIC;
        if (current->allele[0] == current->allele[1])
            return DUCKVEP_PHASE_ARRANGEMENT_HOMOZYGOUS;
        for (size_t prior = 0u; prior < site; prior++) {
            if (sites[prior].source_id == current->source_id)
                return DUCKVEP_PHASE_ARRANGEMENT_DUPLICATE_SOURCE_ID;
        }
        if (arrangement_is_unit_start(sites, site)) unit_count++;
    }
    if (unit_count > sizeof(size_t) * CHAR_BIT)
        return DUCKVEP_PHASE_ARRANGEMENT_EXPLOSION;
    alternative_count = (size_t)1u << (unit_count - 1u);
    if (alternative_count > alternative_limit)
        return DUCKVEP_PHASE_ARRANGEMENT_EXPLOSION;
    if (alternative_count > SIZE_MAX / site_count)
        return DUCKVEP_PHASE_ARRANGEMENT_OVERFLOW;
    required_entries = alternative_count * site_count;
    *required_alternatives = alternative_count;
    if (!alt_lanes && !out_capacity) return DUCKVEP_PHASE_ARRANGEMENT_OK;
    if (out_capacity < required_entries) return DUCKVEP_PHASE_ARRANGEMENT_CAPACITY;

    for (size_t alternative = 0u; alternative < alternative_count; alternative++) {
        for (size_t site = 0u; site < site_count; site++) {
            size_t unit = arrangement_unit_index(sites, site);
            uint16_t lane = sites[site].phase_before[1] && sites[site].allele[0] == 1 ? 1u : 2u;
            if (unit && ((alternative >> (unit - 1u)) & 1u)) lane = (uint16_t)(3u - lane);
            alt_lanes[alternative * site_count + site] = lane;
        }
    }
    return DUCKVEP_PHASE_ARRANGEMENT_OK;
}
