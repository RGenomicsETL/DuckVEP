#include "duckvep_property.h"

static duckvep_phase_arrangement_site_t arrangement_site(uint64_t source_id,
    int32_t first, int32_t second, uint8_t phased, uint8_t ps_present, int64_t ps) {
    duckvep_phase_arrangement_site_t site = {0};
    site.source_id = source_id;
    site.allele[0] = first;
    site.allele[1] = second;
    site.source_alt_count = 1u;
    site.phase_before[1] = phased;
    site.phase_set_present = ps_present;
    site.phase_set = ps;
    return site;
}

TEST diploid_arrangements_enumerate_complete_canonical_phase_sets(void) {
    duckvep_phase_arrangement_site_t unphased[] = {
        arrangement_site(11u, 0, 1, 0u, 0u, 0),
        arrangement_site(12u, 1, 0, 0u, 0u, 0),
        arrangement_site(13u, 0, 1, 0u, 0u, 0)
    };
    static const uint16_t expected_unphased[] = {
        2u, 2u, 2u,
        2u, 1u, 2u,
        2u, 2u, 1u,
        2u, 1u, 1u
    };
    uint16_t lanes[12], untouched[12], before[12];
    size_t required = 0u;
    for (size_t i = 0u; i < 12u; i++) untouched[i] = UINT16_MAX;
    memcpy(before, untouched, sizeof(before));
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_OK, duckvep_phase_arrange_diploid(
        unphased, 3u, 4u, NULL, 0u, &required));
    ASSERT_EQ(4u, required);
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_CAPACITY, duckvep_phase_arrange_diploid(
        unphased, 3u, 4u, untouched, 11u, &required));
    ASSERT_EQ(4u, required);
    ASSERT_MEM_EQ(before, untouched, sizeof(before));
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_OK, duckvep_phase_arrange_diploid(
        unphased, 3u, 4u, lanes, 12u, &required));
    ASSERT_EQ(4u, required);
    ASSERT_MEM_EQ(expected_unphased, lanes, sizeof(lanes));

    duckvep_phase_arrangement_site_t block[] = {
        arrangement_site(21u, 0, 1, 1u, 1u, 8),
        arrangement_site(22u, 1, 0, 1u, 1u, 8),
        arrangement_site(23u, 1, 0, 0u, 0u, 0)
    };
    static const uint16_t expected_block[] = {2u, 1u, 2u, 2u, 1u, 1u};
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_OK, duckvep_phase_arrange_diploid(
        block, 3u, 2u, lanes, 6u, &required));
    ASSERT_EQ(2u, required);
    ASSERT_MEM_EQ(expected_block, lanes, sizeof(expected_block));

    duckvep_phase_arrangement_site_t mixed_default_block[] = {
        arrangement_site(31u, 0, 1, 0u, 0u, 0),
        arrangement_site(32u, 0, 1, 1u, 0u, 0),
        arrangement_site(33u, 1, 0, 1u, 0u, 0)
    };
    static const uint16_t expected_mixed[] = {2u, 2u, 1u, 2u, 1u, 2u};
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_OK, duckvep_phase_arrange_diploid(
        mixed_default_block, 3u, 2u, lanes, 6u, &required));
    ASSERT_EQ(2u, required);
    ASSERT_MEM_EQ(expected_mixed, lanes, sizeof(expected_mixed));
    PASS();
}

TEST diploid_arrangements_reject_unsupported_calls_without_output(void) {
    duckvep_phase_arrangement_site_t sites[] = {
        arrangement_site(31u, 0, 1, 0u, 0u, 0),
        arrangement_site(32u, 0, 1, 0u, 0u, 0)
    };
    uint16_t lanes[4] = {99u, 99u, 99u, 99u}, before[4];
    size_t required = 77u;
    memcpy(before, lanes, sizeof(lanes));
    sites[1].allele[1] = -1;
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_MISSING, duckvep_phase_arrange_diploid(
        sites, 2u, 2u, lanes, 4u, &required));
    ASSERT_EQ(0u, required);
    ASSERT_MEM_EQ(before, lanes, sizeof(lanes));
    sites[1].allele[1] = 0;
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_HOMOZYGOUS, duckvep_phase_arrange_diploid(
        sites, 2u, 2u, lanes, 4u, &required));
    sites[1].allele[1] = 2;
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_MULTIALLELIC, duckvep_phase_arrange_diploid(
        sites, 2u, 2u, lanes, 4u, &required));
    sites[1].allele[1] = 1;
    sites[1].sample_index = 1u;
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_MIXED_SAMPLE, duckvep_phase_arrange_diploid(
        sites, 2u, 2u, lanes, 4u, &required));
    sites[1].sample_index = 0u;
    sites[1] = sites[0];
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_DUPLICATE_SOURCE_ID, duckvep_phase_arrange_diploid(
        sites, 2u, 2u, lanes, 4u, &required));
    sites[1].source_id = 32u;
    sites[1].allele[1] = 1;
    ASSERT_EQ(DUCKVEP_PHASE_ARRANGEMENT_EXPLOSION, duckvep_phase_arrange_diploid(
        sites, 2u, 1u, lanes, 4u, &required));
    ASSERT_EQ(0u, required);
    ASSERT_MEM_EQ(before, lanes, sizeof(lanes));
    PASS();
}
