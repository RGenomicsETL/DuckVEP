#include <stdint.h>
#include <stdio.h>

#include "duckvep_phase.h"

static duckvep_phase_arrangement_site_t site(uint64_t source_id, int32_t first,
    int32_t second, uint8_t phased, uint8_t phase_set_present, int64_t phase_set) {
    duckvep_phase_arrangement_site_t out = {0};
    out.source_id = source_id;
    out.allele[0] = first;
    out.allele[1] = second;
    out.source_alt_count = 1u;
    out.phase_before[1] = phased;
    out.phase_set_present = phase_set_present;
    out.phase_set = phase_set;
    return out;
}

static int emit(const char *name, const duckvep_phase_arrangement_site_t *sites, size_t count) {
    uint16_t lanes[16];
    size_t required = 0u;
    if (duckvep_phase_arrange_diploid(sites, count, 8u, lanes, sizeof(lanes) / sizeof(lanes[0]),
            &required) != DUCKVEP_PHASE_ARRANGEMENT_OK) return 1;
    for (size_t arrangement = 0u; arrangement < required; arrangement++) {
        printf("%s\t%zu\t", name, arrangement);
        for (size_t i = 0u; i < count; i++) {
            printf("%s%llu:%u", i ? "," : "", (unsigned long long)sites[i].source_id,
                (unsigned)lanes[arrangement * count + i]);
        }
        printf("\n");
    }
    return 0;
}

int main(void) {
    duckvep_phase_arrangement_site_t unphased[] = {
        site(101u, 0, 1, 0u, 0u, 0), site(102u, 1, 0, 0u, 0u, 0)
    };
    duckvep_phase_arrangement_site_t three[] = {
        site(111u, 0, 1, 0u, 0u, 0), site(112u, 0, 1, 0u, 0u, 0),
        site(113u, 1, 0, 0u, 0u, 0)
    };
    duckvep_phase_arrangement_site_t block[] = {
        site(121u, 0, 1, 1u, 1u, 8), site(122u, 1, 0, 1u, 1u, 8),
        site(123u, 1, 0, 0u, 0u, 0)
    };
    duckvep_phase_arrangement_site_t separate_blocks[] = {
        site(131u, 0, 1, 1u, 1u, 8), site(132u, 1, 0, 1u, 1u, 9)
    };
    duckvep_phase_arrangement_site_t mixed_default_block[] = {
        site(141u, 0, 1, 0u, 0u, 0), site(142u, 0, 1, 1u, 0u, 0),
        site(143u, 1, 0, 1u, 0u, 0)
    };
    return emit("unphased", unphased, 2u) || emit("three", three, 3u) ||
        emit("block", block, 3u) || emit("separate", separate_blocks, 2u) ||
        emit("mixed_default", mixed_default_block, 3u);
}
