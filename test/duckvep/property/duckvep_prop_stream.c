#include "duckvep_property.h"

/* Exhaust actual slot permutations, independently of the reducer's equality
 * shortcut. A missing allele permits every allele in this finite alphabet. */
static int phase_next_permutation(uint8_t *order, uint16_t count) {
    int left = (int)count - 2;
    while (left >= 0 && order[left] >= order[left + 1]) left--;
    if (left < 0) return 0;
    int right = (int)count - 1;
    while (order[right] <= order[left]) right--;
    uint8_t tmp = order[left]; order[left] = order[right]; order[right] = tmp;
    for (int a = left + 1, b = (int)count - 1; a < b; a++, b--) {
        tmp = order[a]; order[a] = order[b]; order[b] = tmp;
    }
    return 1;
}

TEST phase_assignments_match_all_permitted_slot_permutations(void) {
    uint32_t cases = 0u;
    for (uint16_t ploidy = 1u; ploidy <= 5u; ploidy++) {
        uint32_t genotypes = UINT32_C(1) << (2u * ploidy);
        for (uint32_t genotype = 0u; genotype < genotypes; genotype++) {
            for (uint32_t bits = 0u; bits < (UINT32_C(1) << ploidy); bits++) {
                cases++;
                int32_t alleles[5];
                uint8_t phase[5], order[5], possible[5] = {0};
                uint8_t global = 0u;
                uint16_t unphased = 0u;
                duckvep_phase_summary_t summary = {0};
                for (uint16_t i = 0u; i < ploidy; i++) {
                    alleles[i] = (int32_t)((genotype >> (2u * i)) & 3u) - 1;
                    phase[i] = (uint8_t)((bits >> i) & 1u);
                    order[i] = (uint8_t)i;
                    if (!phase[i]) unphased++;
                    global |= alleles[i] < 0 ? 7u : (uint8_t)(1u << alleles[i]);
                    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_observe(&summary, alleles[i], phase[i]));
                }
                do {
                    int allowed = 1;
                    for (uint16_t i = 0u; i < ploidy; i++)
                        if (phase[i] && order[i] != i) allowed = 0;
                    if (!allowed) continue;
                    for (uint16_t i = 0u; i < ploidy; i++) {
                        int32_t value = alleles[order[i]];
                        possible[i] |= value < 0 ? 7u : (uint8_t)(1u << value);
                    }
                } while (phase_next_permutation(order, ploidy));
                int invariant = ploidy == 1u || (global & (global - 1u)) == 0u;
                uint16_t called_before = 0u;
                for (uint16_t i = 0u; i < ploidy; i++) {
                    duckvep_phase_assignment_t out;
                    int resolved = invariant || phase[i] || unphased == 1u ||
                        (possible[i] & (possible[i] - 1u)) == 0u;
                    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_assign(&summary, i + 1u, called_before,
                        alleles[i], phase[i], DUCKVEP_PHASE_STRICT, &out));
                    ASSERT_EQ(resolved ? i + 1u : 0u, out.lane);
                    ASSERT_EQ(invariant ? DUCKVEP_PHASE_ALL_SETS :
                        resolved ? DUCKVEP_PHASE_SET : DUCKVEP_PHASE_UNRESOLVED, out.scope);
                    ASSERT_EQ(alleles[i] < 0 ? DUCKVEP_PHASE_MISSING :
                        resolved ? DUCKVEP_PHASE_CALLED : DUCKVEP_PHASE_UNPHASED, out.status);
                    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_assign(&summary, i + 1u, called_before,
                        alleles[i], phase[i], DUCKVEP_PHASE_VEP_COMPAT, &out));
                    ASSERT_EQ(alleles[i] < 0 ? 0u : called_before + 1u, out.lane);
                    ASSERT_EQ(DUCKVEP_PHASE_ALLELE_SLOT, out.scope);
                    duckvep_phase_call_status_t expected_status = alleles[i] < 0
                        ? DUCKVEP_PHASE_MISSING : DUCKVEP_PHASE_CALLED;
                    ASSERT_EQ(expected_status, out.status);
                    if (alleles[i] >= 0) called_before++;
                }
            }
        }
    }
    ASSERT_EQ(37448u, cases);
    PASS();
}

TEST raw_gt_source_spelling_and_file_slots_are_distinct(void) {
    static const struct {
        const char *gt;
        uint16_t ploidy;
        uint8_t missing;
        uint32_t slots, first, second;
        duckvep_raw_gt_disposition_t disposition;
    } cases[] = {
        {"0|1", 2, 0, 2, 0, 1, DUCKVEP_RAW_GT_RETAINED},
        {"|0|1", 2, 0, 3, 0, 0, DUCKVEP_RAW_GT_RETAINED},
        {"0/1|2", 3, 0, 2, 0, 2, DUCKVEP_RAW_GT_RETAINED},
        {"0|1/2", 3, 0, 2, 0, 1, DUCKVEP_RAW_GT_RETAINED},
        {".|1", 2, 1, 1, 1, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"1|.", 2, 1, 1, 1, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"./.", 2, 1, 0, UINT32_MAX, UINT32_MAX, DUCKVEP_RAW_GT_OMITTED_EMPTY},
        {"|./.", 2, 1, 2, 0, 0, DUCKVEP_RAW_GT_RETAINED},
        {"./.|1", 3, 1, 2, 0, 1, DUCKVEP_RAW_GT_RETAINED},
        {"|.", 1, 1, 1, 0, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"1", 1, 0, 1, 1, UINT32_MAX, DUCKVEP_RAW_GT_RETAINED},
        {"0", 1, 0, 0, UINT32_MAX, UINT32_MAX, DUCKVEP_RAW_GT_OMITTED_REFERENCE},
        {"0|0/0", 3, 0, 0, UINT32_MAX, UINT32_MAX, DUCKVEP_RAW_GT_OMITTED_REFERENCE},
        {"|0|0", 2, 0, 3, 0, 0, DUCKVEP_RAW_GT_RETAINED},
        {"01/2", 2, 0, 2, 1, 2, DUCKVEP_RAW_GT_RETAINED},
        {"/1|0", 2, 0, 2, 0, 0, DUCKVEP_RAW_GT_RETAINED}
    };
    for (size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); i++) {
        duckvep_raw_gt_t parsed;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
            (const uint8_t *)cases[i].gt, strlen(cases[i].gt), 2u, &parsed));
        ASSERT_EQ(cases[i].ploidy, parsed.source_ploidy);
        ASSERT_EQ(cases[i].missing, parsed.source_has_missing);
        ASSERT_EQ(cases[i].slots, parsed.parsed_slots);
        ASSERT_EQ(cases[i].first, parsed.allele_index[0]);
        ASSERT_EQ(cases[i].second, parsed.allele_index[1]);
        ASSERT_EQ(cases[i].disposition, parsed.disposition);
    }
    static const char *invalid[] = {"|", "/", "||1", "1|", "1.", "-1", "0//1", "x", "1\\2"};
    duckvep_raw_gt_t parsed, zero = {0};
    memset(&zero, 0, sizeof(zero));
    for (size_t i = 0u; i < sizeof(invalid) / sizeof(invalid[0]); i++) {
        memset(&parsed, 0xff, sizeof(parsed));
        ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_SYNTAX, duckvep_phase_parse_vep_raw(
            (const uint8_t *)invalid[i], strlen(invalid[i]), 2u, &parsed));
        ASSERT_MEM_EQ(&zero, &parsed, sizeof(zero));
    }
    ASSERT_EQ(DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE, duckvep_phase_parse_vep_raw(
        (const uint8_t *)"0|3", 3u, 2u, &parsed));
    ASSERT_MEM_EQ(&zero, &parsed, sizeof(zero));
    ASSERT_EQ(DUCKVEP_RAW_GT_ALLELE_OUT_OF_RANGE, duckvep_phase_parse_vep_raw(
        (const uint8_t *)"2147483648", 10u, INT32_MAX, &parsed));
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_ARG, duckvep_phase_parse_vep_raw(NULL, 1u, 2u, &parsed));
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_ARG, duckvep_phase_parse_vep_raw(
        (const uint8_t *)"0", 1u, UINT32_MAX, &parsed));
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_ARG, duckvep_phase_parse_vep_raw(
        (const uint8_t *)"0", 1u, 2u, NULL));
    const uint8_t embedded_nul[] = {'0', '|', 0};
    ASSERT_EQ(DUCKVEP_RAW_GT_INVALID_SYNTAX, duckvep_phase_parse_vep_raw(
        embedded_nul, sizeof(embedded_nul), 2u, &parsed));
    uint8_t *many = malloc(2u * (size_t)UINT16_MAX + 1u);
    ASSERT(many != NULL);
    for (size_t i = 0u; i < 2u * (size_t)UINT16_MAX + 1u; i++) many[i] = i % 2u ? '|' : '0';
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
        many, 2u * (size_t)UINT16_MAX - 1u, 0u, &parsed));
    ASSERT_EQ(UINT16_MAX, parsed.source_ploidy);
    ASSERT_EQ(DUCKVEP_RAW_GT_PLOIDY_LIMIT, duckvep_phase_parse_vep_raw(
        many, 2u * (size_t)UINT16_MAX + 1u, 0u, &parsed));
    ASSERT_MEM_EQ(&zero, &parsed, sizeof(zero));
    free(many);
    PASS();
}

TEST phase_ploidy_and_invalid_observations_are_checked(void) {
    duckvep_phase_summary_t summary = {0}, before;
    duckvep_phase_assignment_t out;
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG,
        duckvep_phase_assign(&summary, 1u, 0u, 0, 0u, DUCKVEP_PHASE_STRICT, &out));
    for (uint32_t i = 0u; i < UINT16_MAX; i++)
        ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_observe(&summary, INT32_MAX, 0u));
    before = summary;
    ASSERT_EQ(DUCKVEP_PHASE_PLOIDY_LIMIT, duckvep_phase_observe(&summary, 1, 1u));
    ASSERT_EQ(0, memcmp(&before, &summary, sizeof summary));
    ASSERT_EQ(DUCKVEP_PHASE_OK, duckvep_phase_assign(&summary, UINT16_MAX, UINT16_MAX - 1u,
        INT32_MAX, 0u, DUCKVEP_PHASE_STRICT, &out));
    ASSERT_EQ(UINT16_MAX, out.lane);
    ASSERT_EQ(DUCKVEP_PHASE_ALL_SETS, out.scope);
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_observe(&summary, -2, 0u));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_observe(&summary, 0, 2u));
    ASSERT_EQ(0, memcmp(&before, &summary, sizeof summary));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_observe(NULL, 0, 0u));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_assign(&summary, 0u, 0u,
        INT32_MAX, 0u, DUCKVEP_PHASE_STRICT, &out));
    ASSERT_EQ(DUCKVEP_PHASE_INVALID_ARG, duckvep_phase_assign(&summary, 1u, 0u,
        INT32_MAX, 0u, (duckvep_phase_policy_t)99, &out));
    PASS();
}

struct carrier_test_pool {
    duckvep_carrier_transcript_t transcripts[3];
    duckvep_carrier_call_t calls[16];
    duckvep_carrier_prefix_t prefixes[64];
    uint32_t active[3];
    duckvep_carrier_bucket_t tx_index[8], call_index[32], prefix_index[128];
};

static duckvep_carrier_buffers_t carrier_test_buffers(struct carrier_test_pool *pool) {
    duckvep_carrier_buffers_t b = {pool->transcripts, pool->calls, pool->prefixes,
        pool->active, pool->tx_index, pool->call_index, pool->prefix_index,
        3u, 16u, 64u, 8u, 32u, 128u};
    return b;
}

static duckvep_carrier_key_t carrier_test_key(uint32_t sample) {
    duckvep_carrier_key_t key = {sample, sample % 3u == 2u ? -5 : 0,
        (uint16_t)(sample % 4u + 1u), 4u, (uint8_t)(sample % 3u != 0u), 0u};
    return key;
}

struct haplotype_stream_scene {
    duckvep_transcript_model_t model;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t sequences;
    uint16_t chrom[64], exon_count[64];
    uint32_t starts[64], ends[64], offsets[64], cdna_starts[64], cdna_ends[64], lengths[64];
    uint64_t sequence_offsets[64], pre_flank_offsets[64], post_flank_offsets[64], cdna_offsets[64];
    uint32_t pre_flank_lengths[64], post_flank_lengths[64], cdna_lengths[64];
    int8_t strands[64];
    uint8_t reference[12], flank[128], cdna[128], noncoding[256];
    struct carrier_test_pool carrier_pool;
    duckvep_haplotype_stored_event_t events[8];
    duckvep_haplotype_projection_t projections[16];
    duckvep_haplotype_contributor_t contributors[8];
    duckvep_haplotype_edit_t edits[8];
    uint64_t edit_event_ids[9]; /* Last slot is a provenance-write canary. */
    size_t edit_source_index[9];
    uint64_t edit_source_cdna[9];
    uint32_t edit_source_exon[9];
    duckvep_haplotype_block_t blocks[8];
    duckvep_carrier_event_t leaf_events[8];
    uint8_t alleles[128], cds[128], protein[128], reference_protein[128];
    uint8_t reference_coding_protein[128];
    duckvep_haplotype_stream_buffers_t buffers;
    duckvep_haplotype_stream_t stream;
};

static void haplotype_stream_scene_prepare(struct haplotype_stream_scene *f, uint32_t count) {
    memset(f, 0, sizeof(*f));
    memset(f->reference, 'A', sizeof(f->reference));
    for (uint32_t i = 0u; i < count; i++) {
        f->starts[i] = 100u + (i / 2u) * 12u + (i % 2u) * 4u;
        f->ends[i] = f->starts[i] + 11u;
        f->strands[i] = 1;
        f->offsets[i] = i;
        f->exon_count[i] = 1u;
        f->cdna_starts[i] = 1u;
        f->cdna_ends[i] = f->lengths[i] = 12u;
    }
    f->model.transcript_count = count; f->model.chrom_id = f->chrom;
    f->model.start1 = f->model.cds_start1 = f->starts;
    f->model.end1 = f->model.cds_end1 = f->ends;
    f->model.strand = f->strands; f->model.exon_offset = f->offsets;
    f->model.exon_count = f->exon_count;
    f->exons.exon_count = count; f->exons.start1 = f->starts; f->exons.end1 = f->ends;
    f->exons.cdna_start1 = f->cdna_starts; f->exons.cdna_end1 = f->cdna_ends;
    f->sequences.transcript_count = count; f->sequences.cds_bytes = f->reference;
    f->sequences.cds_bytes_len = sizeof(f->reference);
    f->sequences.cds_offset = f->sequence_offsets; f->sequences.cds_length = f->lengths;
    f->buffers = (duckvep_haplotype_stream_buffers_t){
        .carriers = carrier_test_buffers(&f->carrier_pool),
        .events = f->events, .projections = f->projections, .alleles = f->alleles,
        .event_capacity = 8u, .projection_capacity = 16u, .allele_capacity = sizeof(f->alleles),
        .leaf_events = f->leaf_events, .contributors = f->contributors, .edits = f->edits,
        .blocks = f->blocks, .edit_event_ids = f->edit_event_ids,
        .edit_source_index = f->edit_source_index, .edit_source_cdna = f->edit_source_cdna,
        .edit_source_exon = f->edit_source_exon,
        .leaf_capacity = 8u, .edit_capacity = 8u, .cds = f->cds, .protein = f->protein,
        .cds_capacity = sizeof(f->cds), .protein_capacity = sizeof(f->protein),
        .reference_protein = f->reference_protein,
        .reference_coding_protein = f->reference_coding_protein,
        .reference_protein_capacity = sizeof(f->reference_protein)};
}

static duckvep_haplotype_stream_status_t haplotype_test_begin_candidates(
    duckvep_haplotype_stream_t *s, const duckvep_haplotype_source_t *source,
    const uint32_t *candidates, size_t count) {
    duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_begin(s, source);
    for (size_t i = 0u; status == DUCKVEP_HAPLOTYPE_STREAM_OK && i < count; i++)
        status = duckvep_haplotype_stream_project(s, candidates[i]);
    return status;
}

/* A complete, fully phased test GT has one ALT lane and REF elsewhere.
 * Exercise the same call interpretation used by the SQL adapter. */
static duckvep_haplotype_stream_status_t haplotype_test_push_called(
    duckvep_haplotype_stream_t *s, const duckvep_carrier_key_t *key) {
    if (s->error) return s->error;
    if (!s->have_current) return DUCKVEP_HAPLOTYPE_STREAM_OK;
    int32_t gt[8] = {0};
    uint8_t phase[8] = {1u,1u,1u,1u,1u,1u,1u,1u};
    if (!key->lane || key->lane > key->ploidy || key->ploidy > 8u)
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    gt[key->lane - 1u] = 1;
    duckvep_haplotype_phase_set_t set = {key->phase_set, key->phase_set_present};
    duckvep_haplotype_call_t call = {gt, phase, key->sample_index, 1u,
        key->ploidy, set, DUCKVEP_PHASE_STRICT};
    const duckvep_haplotype_stored_event_t *event = &s->buffers.events[s->current_event];
    for (uint32_t i = 0u; i < event->projection_count; i++) {
        size_t at = (event->projection_begin + i) % s->buffers.projection_capacity;
        duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_push_call(s,
            s->buffers.projections[at].transcript_index, &call, &set, 1u);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

TEST haplotype_record_order_matches_red_black_insertions(void) {
    /* Independent pointer-tree insertion and preorder traversal. Sorted keys
     * append on the right spine; fix-up uses rotations and uncle recolouring. */
    struct node { unsigned parent, left, right, red; } nodes[2049] = {{0}};
    unsigned root = 0u, stack[2048];
    ASSERT_EQ(0u, duckvep_haplotype_record_order(0u, 1u));
    ASSERT_EQ(0u, duckvep_haplotype_record_order(1u, 0u));
    ASSERT_EQ(0u, duckvep_haplotype_record_order(1u, 2u));
    for (unsigned n = 1u; n <= 2048u; n++) {
        nodes[n].parent = n - 1u; nodes[n].red = 1u;
        if (n == 1u) root = n;
        else nodes[n - 1u].right = n;
        unsigned current = n;
        while (nodes[nodes[current].parent].red) {
            unsigned parent = nodes[current].parent, grand = nodes[parent].parent;
            ASSERT_EQ(parent, nodes[grand].right);
            unsigned uncle = nodes[grand].left;
            if (nodes[uncle].red) {
                nodes[parent].red = nodes[uncle].red = 0u;
                nodes[grand].red = 1u; current = grand;
            } else {
                nodes[parent].red = 0u; nodes[grand].red = 1u;
                nodes[grand].right = nodes[parent].left;
                if (nodes[parent].left) nodes[nodes[parent].left].parent = grand;
                unsigned above = nodes[grand].parent;
                nodes[parent].parent = above;
                if (!above) root = parent;
                else if (nodes[above].left == grand) nodes[above].left = parent;
                else nodes[above].right = parent;
                nodes[parent].left = grand; nodes[grand].parent = parent;
            }
        }
        nodes[root].red = 0u;
        size_t pending = 1u; stack[0] = root;
        uint64_t rank = 0u;
        while (pending) {
            unsigned at = stack[--pending];
            ASSERT_EQ(++rank, duckvep_haplotype_record_order(n, at));
            if (nodes[at].right) stack[pending++] = nodes[at].right;
            if (nodes[at].left) stack[pending++] = nodes[at].left;
        }
        ASSERT_EQ(n, rank);
    }
    for (unsigned bit = 1u; bit < 64u; bit++) {
        uint64_t count = UINT64_MAX >> (63u - bit);
        uint64_t ranks[64];
        ASSERT_EQ(count, duckvep_haplotype_record_order(count, count));
        for (unsigned i = 0u; i <= bit; i++) {
            ranks[i] = duckvep_haplotype_record_order(count, UINT64_C(1) << i);
            ASSERT(ranks[i] && ranks[i] <= count);
            for (unsigned j = 0u; j < i; j++) ASSERT(ranks[j] != ranks[i]);
        }
    }
    PASS();
}

TEST haplotype_record_plan_closes_transcripts_not_later_record_ends(void) {
    uint16_t chrom[] = {0u, 0u, 0u, 1u};
    uint32_t starts[] = {10u, 18u, 40u, 5u}, ends[] = {20u, 30u, 50u, 9u};
    duckvep_transcript_model_t model = {.chrom_id = chrom, .start1 = starts,
        .end1 = ends, .transcript_count = 4u};
    duckvep_haplotype_record_plan_t plan;
    uint64_t buffer, ordinal;
    ASSERT(duckvep_haplotype_record_plan_init(&plan, &model));
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 1u, 9u, &buffer, &ordinal));
    ASSERT_EQ(0u, buffer); ASSERT_EQ(0u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 9u, 10u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(1u, ordinal); ASSERT_EQ(30u, plan.end1);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 29u, 100u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(2u, ordinal); ASSERT_EQ(30u, plan.end1);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 31u, 39u, &buffer, &ordinal));
    ASSERT_EQ(0u, buffer); ASSERT_EQ(0u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 40u, 40u, &buffer, &ordinal));
    ASSERT_EQ(2u, buffer); ASSERT_EQ(1u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 1u, 5u, 5u, &buffer, &ordinal));
    ASSERT_EQ(3u, buffer); ASSERT_EQ(1u, ordinal);
    ASSERT_FALSE(duckvep_haplotype_record_plan_next(&plan, 1u, 4u, 5u, &buffer, &ordinal));
    ASSERT_EQ(0u, buffer); ASSERT_EQ(0u, ordinal);
    ASSERT(duckvep_haplotype_record_plan_init(&plan, &model));
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 9u, 41u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(50u, plan.end1);
    ASSERT(duckvep_haplotype_record_plan_next(&plan, 0u, 35u, 35u, &buffer, &ordinal));
    ASSERT_EQ(1u, buffer); ASSERT_EQ(2u, ordinal);
    starts[1] = 1u;
    ASSERT_FALSE(duckvep_haplotype_record_plan_init(&plan, &model));
    PASS();
}

TEST haplotype_record_plan_matches_dense_closure(void) {
    uint64_t seeds[] = {173u, 20260906u};
    for (unsigned seed = 0u; seed < 2u; seed++) {
        uint64_t state = seeds[seed];
        for (unsigned trial = 0u; trial < 4096u; trial++) {
            uint16_t chrom[32];
            uint32_t starts[32], ends[32];
            for (unsigned i = 0u; i < 32u; i++) {
                state = state * UINT64_C(6364136223846793005) + 1u;
                chrom[i] = (uint16_t)(i / 16u);
                starts[i] = 1u + (i % 16u) * 16u;
                ends[i] = starts[i] + (uint32_t)((state >> 32u) % 32u);
            }
            duckvep_transcript_model_t model = {.chrom_id = chrom, .start1 = starts,
                .end1 = ends, .transcript_count = 32u};
            duckvep_haplotype_record_plan_t plan;
            ASSERT(duckvep_haplotype_record_plan_init(&plan, &model));
            uint64_t expected_buffer = 0u, expected_ordinal = 0u;
            for (uint16_t region = 0u; region < 2u; region++) {
                uint32_t start = 1u, expected_end = 0u;
                for (unsigned record = 0u; record < 64u; record++) {
                    state = state * UINT64_C(6364136223846793005) + 1u;
                    start += (uint32_t)((state >> 32u) % 8u);
                    uint32_t end = start + (uint32_t)((state >> 48u) % 64u);
                    if (record % 13u == 0u) end += 200u;
                    if (start > expected_end) {
                        uint32_t reachable = end, previous;
                        int overlaps = 0;
                        do {
                            previous = reachable;
                            for (unsigned i = 0u; i < 32u; i++) {
                                if (chrom[i] != region || starts[i] > reachable || ends[i] < start) continue;
                                overlaps = 1;
                                if (ends[i] > reachable) reachable = ends[i];
                            }
                        } while (reachable != previous);
                        expected_end = overlaps ? reachable : 0u;
                        if (overlaps) { expected_buffer++; expected_ordinal = 0u; }
                    }
                    uint64_t buffer, ordinal;
                    ASSERT(duckvep_haplotype_record_plan_next(&plan, region, start, end, &buffer, &ordinal));
                    if (expected_end) expected_ordinal++;
                    ASSERT_EQ(expected_end ? expected_buffer : 0u, buffer);
                    ASSERT_EQ(expected_end ? expected_ordinal : 0u, ordinal);
                    ASSERT_EQ(expected_end, plan.end1);
                }
            }
        }
    }
    PASS();
}

TEST haplotype_stream_raw_record_selection_and_tie_order(void) {
    for (unsigned planned = 0u; planned < 2u; planned++) {
        for (unsigned shadowed = 0u; shadowed < 2u; shadowed++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            duckvep_raw_gt_t call;
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                (const uint8_t *)"1|1", 3u, 1u, &call));
            for (unsigned i = 0u; i < 2u; i++) {
                duckvep_haplotype_source_t source = {i + 1u, (const uint8_t *)"A",
                    (const uint8_t *)(i ? "G" : "C"), 100u, 0u, 1u, 1u, 1u, 1u,
                    planned ? 2u - i : 0u, 0u};
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                    &f.stream, 0u, 0u, &call, (uint8_t)!(i && shadowed)));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(planned || shadowed ? 'C' : 'G', leaf.cds[0]);
            ASSERT_EQ(shadowed ? 1u : 2u, leaf.edit_count);
            ASSERT_EQ(2u, leaf.contributor_count);
            ASSERT_EQ(shadowed ? DUCKVEP_CDS_EDIT_SOURCE_SHADOWED : DUCKVEP_CDS_EDIT_OK,
                leaf.contributors[1].projection_status);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
        }
    }
    PASS();
}

TEST haplotype_stream_raw_record_slots_keep_conditional_deletion_provenance(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    duckvep_raw_gt_t calls[2];
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
        (const uint8_t *)".|1", 3u, 1u, &calls[0]));
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
        (const uint8_t *)"1|1", 3u, 1u, &calls[1]));
    duckvep_haplotype_source_t sources[] = {
        {11u, (const uint8_t *)"A", (const uint8_t *)"C", 101u, 0u, 1u, 1u, 1u, 1u, 0u, 0u},
        {11u, (const uint8_t *)"A", (const uint8_t *)"", 101u, 0u, 1u, 0u, UINT32_MAX, 1u, 0u, 0u},
        {12u, (const uint8_t *)"A", (const uint8_t *)"G", 104u, 0u, 1u, 1u, 1u, 1u, 0u, 0u}
    };
    for (size_t i = 0; i < 3u; i++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &sources[i]));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
            s, 0u, 7u, &calls[i == 2u], 1u));
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    unsigned lanes = 0u;
    while (duckvep_haplotype_stream_next(s, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(2u, leaf.contributor_count);
        ASSERT_EQ(2u, leaf.edit_count);
        ASSERT_EQ(11u, leaf.edit_event_ids[0]); ASSERT_EQ(12u, leaf.edit_event_ids[1]);
        const duckvep_carrier_call_t *carrier = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
        ASSERT(carrier != NULL); ASSERT_EQ(7u, carrier->key.sample_index);
        ASSERT_EQ(2u, carrier->key.ploidy); ASSERT_EQ(0u, carrier->key.phase_set_present);
        ASSERT_EQ(1u, leaf.carriers.call_count);
        unsigned lane = carrier->key.lane;
        ASSERT(lane == 1u || lane == 2u); lanes |= 1u << lane;
        ASSERT_EQ(lane == 1u ? 0 : -1, leaf.nominal_length_diff);
        const char *expected = lane == 1u ? "ACAAGAAAAAAA" : "AAAGAAAAAAA";
        ASSERT_EQ(strlen(expected), leaf.cds_length);
        ASSERT_MEM_EQ(expected, leaf.cds, leaf.cds_length);
        ASSERT_EQ(lane == 1u ? 1u : UINT32_MAX, leaf.contributors[0].source.allele_index);
        ASSERT_EQ(lane == 1u ? 1u : 0u, leaf.contributors[0].source.alt_len);
        ASSERT_EQ(lane == 1u ? DUCKVEP_CARRIER_CALLED : 0u,
            leaf.contributors[0].evidence_flags & DUCKVEP_CARRIER_CALLED);
        ASSERT(leaf.evidence_flags & DUCKVEP_CARRIER_CONDITIONAL);
        ASSERT(leaf.evidence_flags & DUCKVEP_CARRIER_MISSING);
    }
    ASSERT_EQ(6u, lanes);
    ASSERT_EQ(2u, s->completed_leaves); ASSERT_EQ(3u, s->projected_events);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    PASS();
}

TEST haplotype_stream_raw_record_reference_observations_require_matching_cds(void) {
    for (unsigned mismatch = 0u; mismatch < 2u; mismatch++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        const uint8_t *ref = (const uint8_t *)(mismatch ? "CAA" : "AAA");
        duckvep_haplotype_source_t source = {1u, ref, ref, 100u, 0u, 3u, 3u, 0u, 1u, 0u, 0u};
        duckvep_raw_gt_t call;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw((const uint8_t *)".", 1u, 1u, &call));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(2u, leaf.carriers.call_count); ASSERT_EQ(1u, leaf.contributor_count);
        ASSERT_EQ(0u, leaf.edit_count); ASSERT_EQ(0u, leaf.block_count);
        ASSERT_EQ(0, leaf.nominal_length_diff);
        ASSERT_EQ(DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_CONDITIONAL, leaf.evidence_flags);
        ASSERT_EQ(0u, leaf.contributors[0].source.allele_index);
        if (mismatch) {
            ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH, leaf.projection_status);
            ASSERT(leaf.cds == NULL); ASSERT(leaf.protein == NULL);
        } else {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
            ASSERT_EQ(12u, leaf.cds_length); ASSERT_MEM_EQ(f.reference, leaf.cds, 12u);
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    }
    /* Haplosaurus omits unprefixed homozygous REF calls. Its retained REF slots
     * still have no mutation when their literal bases fail /^[ACGT]*$/. A
     * positive ALT ordinal with identical bytes remains a selected ALT. */
    static const struct {
        uint8_t source, model, alternate;
        int skipped;
    } bases[] = {{'A', 'A', 'C', 0}, {'N', 'N', 'N', 1}, {'n', 'N', 'n', 1},
                 {'a', 'A', 'a', 1}, {'U', 'T', 'U', 1}};
    static const struct {
        const char *gt;
        uint32_t allele_index;
        int missing, omitted;
    } calls[] = {{"0|0", 0u, 0, 1}, {"|0|0", 0u, 0, 0}, {"1|1", 1u, 0, 0},
                 {".", 0u, 1, 0}, {"|./.", 0u, 1, 0}};
    for (size_t b = 0u; b < sizeof(bases) / sizeof(bases[0]); b++) {
        for (size_t c = 0u; c < sizeof(calls) / sizeof(calls[0]); c++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.reference[0] = bases[b].model;
            duckvep_haplotype_stream_t *s = &f.stream;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                s, &f.model, &f.exons, &f.sequences, &f.buffers));
            duckvep_haplotype_source_t source = {
                .event_id = 1u, .ref = &bases[b].source,
                .alt = calls[c].allele_index ? &bases[b].alternate : &bases[b].source,
                .pos1 = 100u, .ref_len = 1u, .alt_len = 1u,
                .allele_index = calls[c].allele_index, .source_record = 1u};
            duckvep_raw_gt_t call;
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                (const uint8_t *)calls[c].gt, strlen(calls[c].gt), 1u, &call));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_raw_call(s, 0u, 7u, &call, 1u));
            if (calls[c].omitted) {
                ASSERT_EQ(DUCKVEP_RAW_GT_OMITTED_REFERENCE, call.disposition);
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
                continue;
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
            int skipped_alt = calls[c].allele_index && bases[b].skipped;
            int conditional = calls[c].missing || skipped_alt;
            if (leaf.projection_status != DUCKVEP_CDS_EDIT_OK ||
                leaf.sequence_status != (conditional ? DUCKVEP_HAPLOTYPE_CONDITIONAL : DUCKVEP_HAPLOTYPE_OK))
                fprintf(stderr, "reference slot base=%c GT=%s ordinal=%u projection=%d sequence=%d\n",
                    bases[b].source, calls[c].gt, calls[c].allele_index,
                    leaf.projection_status, leaf.sequence_status);
            ASSERT_EQ(conditional ? DUCKVEP_HAPLOTYPE_CONDITIONAL : DUCKVEP_HAPLOTYPE_OK,
                leaf.sequence_status);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(12u, leaf.cds_length);
            int applied = calls[c].allele_index && !bases[b].skipped;
            uint8_t expected[12];
            memcpy(expected, f.reference, sizeof(expected));
            if (applied) expected[0] = bases[b].alternate;
            ASSERT_MEM_EQ(expected, leaf.cds, 12u);
            ASSERT_EQ(applied ? 1u : 0u, leaf.edit_count);
            ASSERT_EQ(applied ? 1u : 0u, leaf.block_count);
            ASSERT_EQ(2u, leaf.carriers.call_count);
            ASSERT_EQ(calls[c].missing || calls[c].allele_index ? 1u : 0u, leaf.contributor_count);
            uint8_t evidence = calls[c].missing
                ? DUCKVEP_CARRIER_MISSING | DUCKVEP_CARRIER_CONDITIONAL
                : calls[c].allele_index ? DUCKVEP_CARRIER_CALLED : 0u;
            if (skipped_alt) evidence |= DUCKVEP_CARRIER_CONDITIONAL;
            ASSERT_EQ(evidence, leaf.evidence_flags);
            if (leaf.contributor_count) {
                ASSERT_EQ(calls[c].allele_index, leaf.contributors[0].source.allele_index);
                ASSERT_EQ(evidence, leaf.contributors[0].evidence_flags);
                ASSERT_EQ(skipped_alt ? DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED : DUCKVEP_CDS_EDIT_OK,
                    leaf.contributors[0].projection_status);
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
        }
    }
    PASS();
}

TEST haplotype_lowercase_source_preserves_bases_and_raw_identity(void) {
    ASSERT(duckvep_codon_is_start((const uint8_t *)"ATG", DUCKVEP_CODON_TABLE_VERT_MITO));
    ASSERT(duckvep_codon_is_start((const uint8_t *)"ATA", DUCKVEP_CODON_TABLE_VERT_MITO));
    ASSERT_EQ('M', duckvep_translate_codon("ATA", DUCKVEP_CODON_TABLE_VERT_MITO));
    ASSERT(!duckvep_codon_is_start((const uint8_t *)"ATA", DUCKVEP_CODON_TABLE_STANDARD));
    static const struct {
        uint32_t position;
        const char *ref, *alt, *cds, *protein;
        uint64_t mask;
    } cases[] = {
        {105u, "t", "c", "ATGGCCGAATAA", "MAE*", DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)},
        {104u, "c", "t", "ATGGTTGAATAA", "MVE*", DUCKVEP_SO(DUCKVEP_SO_MISSENSE)}
    };
    for (size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); i++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "atggctgaataa", 12u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_haplotype_source_t source = {i + 1u, (const uint8_t *)cases[i].ref,
            (const uint8_t *)cases[i].alt, cases[i].position, 0u, 1u, 1u, 0u, 0u, 0u, 1u};
        duckvep_carrier_key_t key = {0u, 0, 1u, 2u, 1u, 0u};
        uint32_t tx = 0u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
            haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
            duckvep_haplotype_stream_finish(&f.stream));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.prediction_status);
        ASSERT_MEM_EQ("atggctgaataa", leaf.reference_cds, 12u);
        ASSERT_MEM_EQ(cases[i].cds, leaf.cds, 12u);
        ASSERT_MEM_EQ(cases[i].protein, leaf.protein, 4u);
        ASSERT_EQ(cases[i].mask, leaf.haplotype_so_mask);
        ASSERT_EQ(duckvep_so_impact(cases[i].mask), duckvep_so_impact(leaf.haplotype_so_mask));
        ASSERT_EQ(1u, leaf.contributor_count);
        ASSERT_EQ(source.event_id, leaf.contributors[0].source.event_id);
        ASSERT_MEM_EQ(cases[i].ref, leaf.contributors[0].source.ref, 1u);
        ASSERT_MEM_EQ(cases[i].alt, leaf.contributors[0].source.alt, 1u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }
    PASS();
}

TEST haplotype_stream_skipped_sources_match_raw_haplosaurus_slots(void) {
    /* Actual 116 Runner/container observations: raw_indel_oracle.Vla8xO,
     * final receipt 96fef08540301180ea41b70e4c0342889286563cbc3f6ed46648f7cab0368944;
     * GT=1 companion e3efba82cbb9d53623753b37f8747be580e48180e0700c558c60aec5f6a2f51a.
     * GT=1 still has two upstream slots: its undefined second ALT deletes REF.
     * These are raw replay expectations, not independent HGVSp expectations. */
    static const struct {
        uint64_t id;
        uint32_t pos1;
        const char *reference, *ref, *alt, *reference_protein;
        const char *cds, *protein, *undefined_cds, *undefined_protein;
        int skipped;
    } cases[] = {
        {110256u, 16u, "ATGGCTGCCTAA", "T", "TGCC", "MAA*",
            "ATGGCTGCCGCCTAA", "MAAA*", "ATGGCGCCTAA", "MAP", 0},
        {143856u, 16u, "ATGGCNGCCTAA", "N", "NGCC", "MAA*",
            "ATGGCNGCCTAA", "MAA*", "ATGGCGCCTAA", "MAP", 1},
        {146516u, 14u, "ATGNCNGCCTAA", "N", "NGCC", "MXA*",
            "ATGNCNGCCTAA", "MXA*", "ATGCNGCCTAA", "MXP", 1},
        {5400u, 14u, "ATGNAAGCCTAA", "NAAG", "N", "MXA*",
            "ATGNAAGCCTAA", "MXA*", "ATGCCTAA", "MP", 1},
        {143844u, 15u, "ATGGCNGCCTAA", "CN", "C", "MAA*",
            "ATGGCGCCTAA", "MAP", "ATGGGCCTAA", "MGL", 0}
    };
    unsigned leaves = 0u, carriers = 0u;
    for (size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); i++) {
        for (unsigned haploid_gt = 0u; haploid_gt < 2u; haploid_gt++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.starts[0] = 11u;
            f.ends[0] = 22u;
            memcpy(f.reference, cases[i].reference, 12u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            duckvep_raw_gt_t call;
            const char *gt = haploid_gt ? "1" : "1|1";
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                (const uint8_t *)gt, strlen(gt), 1u, &call));
            ASSERT_EQ(haploid_gt ? 1u : 2u, call.source_ploidy);
            for (unsigned undefined = 0u; undefined <= haploid_gt; undefined++) {
                duckvep_haplotype_source_t source = {
                    .event_id = cases[i].id, .pos1 = cases[i].pos1,
                    .ref = (const uint8_t *)cases[i].ref,
                    .alt = (const uint8_t *)(undefined ? "" : cases[i].alt),
                    .ref_len = (uint16_t)strlen(cases[i].ref),
                    .alt_len = undefined ? 0u : (uint16_t)strlen(cases[i].alt),
                    .allele_index = undefined ? UINT32_MAX : 1u, .source_record = 1u};
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                    duckvep_haplotype_stream_begin(&f.stream, &source));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                    duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                    duckvep_haplotype_stream_push_raw_call(&f.stream, 0u, 7u, &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
                duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            duckvep_haplotype_stream_status_t status;
            unsigned lanes = 0u, groups = 0u;
            while ((status = duckvep_haplotype_stream_next(&f.stream, &leaf)) ==
                   DUCKVEP_HAPLOTYPE_STREAM_OK) {
                const duckvep_carrier_call_t *first = duckvep_carriers_call(
                    &f.stream.carriers, leaf.carriers.first_call);
                ASSERT(first != NULL);
                int undefined = haploid_gt && first->key.lane == 2u;
                int skipped = !undefined && cases[i].skipped;
                const char *cds = undefined ? cases[i].undefined_cds : cases[i].cds;
                const char *protein = undefined ? cases[i].undefined_protein : cases[i].protein;
                int64_t length_diff = (int64_t)strlen(cds) - 12;
                uint32_t flags = length_diff ? DUCKVEP_HAPLOTYPE_FLAG_INDEL : 0u;
                if (length_diff % 3) flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
                ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                ASSERT_EQ(skipped || undefined ? DUCKVEP_HAPLOTYPE_CONDITIONAL :
                    DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
                ASSERT_EQ(strlen(cds), leaf.cds_length);
                ASSERT_MEM_EQ(cds, leaf.cds, leaf.cds_length);
                ASSERT_EQ(strlen(protein), leaf.protein_length);
                ASSERT_MEM_EQ(protein, leaf.protein, leaf.protein_length);
                ASSERT_MEM_EQ(cases[i].reference, leaf.reference_cds, 12u);
                ASSERT_EQ(strlen(cases[i].reference_protein), leaf.reference_protein_length);
                ASSERT_MEM_EQ(cases[i].reference_protein, leaf.reference_protein,
                    leaf.reference_protein_length);
                ASSERT_EQ(length_diff, leaf.nominal_length_diff);
                ASSERT_EQ(flags, leaf.flags);
                ASSERT_EQ(1u, leaf.contributor_count);
                ASSERT_EQ(skipped ? 0u : 1u, leaf.edit_count);
                ASSERT_EQ(skipped ? 0u : 1u, leaf.block_count);
                const duckvep_haplotype_contributor_t *contributor = leaf.contributors;
                ASSERT_EQ(cases[i].id, contributor->source.event_id);
                ASSERT_EQ(cases[i].pos1, contributor->source.pos1);
                ASSERT_EQ(1u, contributor->source.source_record);
                ASSERT_EQ(undefined ? UINT32_MAX : 1u, contributor->source.allele_index);
                ASSERT_EQ(strlen(cases[i].ref), contributor->source.ref_len);
                ASSERT_MEM_EQ(cases[i].ref, contributor->source.ref, contributor->source.ref_len);
                ASSERT_EQ(undefined ? 0u : strlen(cases[i].alt), contributor->source.alt_len);
                if (!undefined)
                    ASSERT_MEM_EQ(cases[i].alt, contributor->source.alt, contributor->source.alt_len);
                ASSERT(contributor->prepared != NULL);
                ASSERT_EQ(skipped ? DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED :
                    DUCKVEP_CDS_EDIT_OK, contributor->projection_status);
                /* Disjoint replay reports applied sources through edit IDs;
                 * source_replaced is reserved for ordered overlapping replay. */
                ASSERT_EQ(0u, leaf.ordered_replacements);
                ASSERT_EQ(0u, contributor->source_replaced);
                ASSERT_EQ(skipped, contributor->projected == NULL);
                uint8_t evidence = undefined ? DUCKVEP_CARRIER_CONDITIONAL : DUCKVEP_CARRIER_CALLED;
                if (skipped) evidence |= DUCKVEP_CARRIER_CONDITIONAL;
                ASSERT_EQ(evidence, contributor->evidence_flags);
                ASSERT_EQ(evidence, leaf.evidence_flags);
                if (!skipped) ASSERT_EQ(cases[i].id, leaf.edit_event_ids[0]);
                ASSERT_EQ(haploid_gt ? 1u : 2u, leaf.carriers.call_count);
                unsigned group_carriers = 0u;
                for (uint32_t id = leaf.carriers.first_call; id;) {
                    const duckvep_carrier_call_t *carrier =
                        duckvep_carriers_call(&f.stream.carriers, id);
                    ASSERT(carrier != NULL);
                    ASSERT_EQ(7u, carrier->key.sample_index);
                    ASSERT_EQ(2u, carrier->key.ploidy);
                    ASSERT_EQ(0u, carrier->key.phase_set_present);
                    ASSERT(carrier->key.lane == 1u || carrier->key.lane == 2u);
                    ASSERT_EQ(0u, lanes & (1u << carrier->key.lane));
                    lanes |= 1u << carrier->key.lane;
                    id = carrier->next_leaf;
                    group_carriers++;
                }
                ASSERT_EQ(leaf.carriers.call_count, group_carriers);
                carriers += group_carriers;
                groups++;
                leaves++;
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
            ASSERT_EQ(haploid_gt ? 2u : 1u, groups);
            ASSERT_EQ(6u, lanes);
            ASSERT_MEM_EQ(cases[i].reference, f.reference, 12u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(&f.stream));
        }
    }
    ASSERT_EQ(15u, leaves);
    ASSERT_EQ(20u, carriers);
    PASS();
}

TEST haplotype_stream_skipped_source_does_not_hide_valid_replay(void) {
    /* Actual two-source GT=1|1 observation, receipt
     * 0eb53b96f6fab741bed65faf42fcaebed6eb4a27125e17fa136ff16d84cc58e2:
     * N>NGCC is observed but skipped; C>A alone produces MAD*. The minus
     * model is the exact reverse-complement coordinate transform. */
    for (int strand = -1; strand <= 1; strand += 2) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.starts[0] = 11u;
        f.ends[0] = 22u;
        f.strands[0] = (int8_t)strand;
        memcpy(f.reference, "ATGGCNGCCTAA", 12u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_raw_gt_t call;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
            (const uint8_t *)"1|1", 3u, 1u, &call));
        for (unsigned i = 0u; i < 2u; i++) {
            int skipped = strand > 0 ? i == 0u : i == 1u;
            duckvep_haplotype_source_t source = {
                .event_id = skipped ? 1u : 2u, .source_record = 1u, .allele_index = 1u,
                .pos1 = skipped ? (strand > 0 ? 16u : 17u) : (strand > 0 ? 18u : 15u),
                .ref = (const uint8_t *)(skipped ? "N" : strand > 0 ? "C" : "G"),
                .alt = (const uint8_t *)(skipped ? (strand > 0 ? "NGCC" : "GGCN") :
                    strand > 0 ? "A" : "T"),
                .ref_len = 1u, .alt_len = skipped ? 4u : 1u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_begin(&f.stream, &source));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_project(&f.stream, 0u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_raw_call(&f.stream, 0u, 7u, &call, 1u));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
            duckvep_haplotype_stream_finish(&f.stream));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
        ASSERT_EQ(12u, leaf.cds_length);
        ASSERT_MEM_EQ("ATGGCNGACTAA", leaf.cds, 12u);
        ASSERT_EQ(4u, leaf.protein_length);
        ASSERT_MEM_EQ("MAD*", leaf.protein, 4u);
        ASSERT_EQ(2u, leaf.carriers.call_count);
        ASSERT_EQ(2u, leaf.contributor_count);
        ASSERT_EQ(1u, leaf.edit_count);
        ASSERT_EQ(1u, leaf.block_count);
        ASSERT_EQ(2u, leaf.edit_event_ids[0]);
        ASSERT_EQ(0, leaf.nominal_length_diff);
        ASSERT_EQ(0u, leaf.flags);
        for (size_t i = 0u; i < leaf.contributor_count; i++) {
            const duckvep_haplotype_contributor_t *contributor = &leaf.contributors[i];
            int skipped = contributor->source.event_id == 1u;
            ASSERT_EQ(skipped ? DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED :
                DUCKVEP_CDS_EDIT_OK, contributor->projection_status);
            ASSERT_EQ(0u, leaf.ordered_replacements);
            ASSERT_EQ(0u, contributor->source_replaced);
            ASSERT_EQ(skipped, contributor->projected == NULL);
            ASSERT_EQ(DUCKVEP_CARRIER_CALLED | (skipped ? DUCKVEP_CARRIER_CONDITIONAL : 0u),
                contributor->evidence_flags);
        }
        ASSERT_MEM_EQ("ATGGCNGCCTAA", f.reference, 12u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(&f.stream));
    }
    PASS();
}

TEST haplotype_stream_retained_reference_replaces_but_omitted_missing_does_not(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        for (unsigned missing = 0u; missing < 2u; missing++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.strands[0] = reverse ? -1 : 1;
            duckvep_haplotype_stream_t *s = &f.stream;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                s, &f.model, &f.exons, &f.sequences, &f.buffers));
            const uint8_t *ref = (const uint8_t *)(reverse ? "TTT" : "AAA");
            const uint8_t *alt = (const uint8_t *)(reverse ? "G" : "C");
            duckvep_haplotype_source_t sources[] = {
                {1u, ref, ref, 100u, 0u, 3u, 3u, 0u, 1u, 0u, 0u},
                {2u, ref, alt, 101u, 0u, 1u, 1u, 1u, 1u, 0u, 0u}
            };
            for (size_t i = 0u; i < 2u; i++) {
                duckvep_raw_gt_t call;
                const char *gt = i ? "1|1" : missing ? ".|." : "0|1";
                ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                    (const uint8_t *)gt, strlen(gt), 1u, &call));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &sources[i]));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
            duckvep_haplotype_leaf_t leaf;
            unsigned checked = 0u;
            while (duckvep_haplotype_stream_next(s, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                const duckvep_carrier_call_t *carrier = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
                ASSERT(carrier != NULL);
                if (!missing && carrier->key.lane == 2u) continue;
                ASSERT_EQ(12u, leaf.cds_length);
                ASSERT_EQ(0, leaf.nominal_length_diff);
                ASSERT_EQ(2u, leaf.contributor_count);
                ASSERT_EQ(missing ? 1u : 2u, leaf.edit_count);
                ASSERT_EQ(missing ? 0u : 1u, leaf.ordered_replacements);
                ASSERT_EQ(missing ? DUCKVEP_HAPLOTYPE_CONDITIONAL : DUCKVEP_HAPLOTYPE_OK,
                    leaf.sequence_status);
                uint8_t expected[12]; memcpy(expected, f.reference, 12u);
                if (missing) expected[reverse ? 10u : 1u] = 'C';
                ASSERT_MEM_EQ(expected, leaf.cds, 12u);
                ASSERT_EQ(0u, leaf.evidence_flags & DUCKVEP_CARRIER_REFERENCE_REPLAY);
                if (!missing) {
                    ASSERT_EQ(1u, leaf.block_count);
                    ASSERT_EQ(2u, leaf.blocks[0].edit_count);
                    ASSERT_EQ(0, leaf.blocks[0].length_diff);
                    ASSERT_EQ(1u, leaf.contributors[0].source_replaced);
                }
                checked++;
            }
            ASSERT_EQ(1u, checked);
        }
    }
    PASS();
}

TEST haplotype_stream_nominal_lengths_survive_clipped_and_disjoint_replay(void) {
    static const struct {
        uint32_t start[2];
        uint16_t ref_len[2];
        const char *alt[2];
        int64_t nominal;
        size_t length, edits, sources;
        uint32_t flags;
    } cases[] = {
        /* The undefined slot deletes the last base; the following substitution
         * inserts at that emptied position while retaining nominal REF length. */
        {{12u, 12u}, {1u, 1u}, {"", "C"}, -1, 12u, 2u, 2u, 3u},
        /* The six-base source sees only three bases after the other deletion.
         * Its ALT already matches, but its nominal -3 remains in the sum. */
        {{7u, 9u}, {6u, 4u}, {"AAA", "A"}, -6, 9u, 1u, 2u, 1u},
        {{9u, 0u}, {4u, 0u}, {"A", ""}, -3, 9u, 1u, 1u, 1u},
        {{4u, 9u}, {1u, 2u}, {"AAA", "A"}, 1, 13u, 2u, 2u, 3u}
    };
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        for (size_t scenario = 0u; scenario < sizeof(cases) / sizeof(cases[0]); scenario++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.strands[0] = reverse ? -1 : 1;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            uint8_t refs[2][6], alts[2][3];
            duckvep_haplotype_source_t sources[2];
            for (size_t i = 0u; i < cases[scenario].sources; i++) {
                size_t alt_length = strlen(cases[scenario].alt[i]);
                memset(refs[i], reverse ? 'T' : 'A', cases[scenario].ref_len[i]);
                for (size_t j = 0u; j < alt_length; j++) alts[i][j] = (uint8_t)
                    haplo_test_variant_from_tx_base(cases[scenario].alt[i][reverse ? alt_length - 1u - j : j],
                        (int)reverse);
                sources[i] = (duckvep_haplotype_source_t){0u, refs[i], alts[i],
                    reverse ? 113u - cases[scenario].start[i] - cases[scenario].ref_len[i]
                            : 99u + cases[scenario].start[i],
                    0u, cases[scenario].ref_len[i], (uint16_t)alt_length,
                    scenario == 0u && i == 0u ? UINT32_MAX : 1u, 1u, 0u, 0u};
            }
            if (cases[scenario].sources == 2u && sources[0].pos1 > sources[1].pos1) {
                duckvep_haplotype_source_t tmp = sources[0]; sources[0] = sources[1]; sources[1] = tmp;
            }
            for (size_t i = 0u; i < cases[scenario].sources; i++) {
                sources[i].event_id = i + 1u;
                const char *gt = scenario ? "1|1" : sources[i].allele_index == UINT32_MAX ? ".|1" : "0|1";
                duckvep_raw_gt_t call;
                ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                    (const uint8_t *)gt, strlen(gt), 1u, &call));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &sources[i]));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                    &f.stream, 0u, 0u, &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
            ASSERT_EQ(scenario ? DUCKVEP_HAPLOTYPE_OK : DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
            ASSERT_EQ(cases[scenario].nominal, leaf.nominal_length_diff);
            ASSERT_EQ(cases[scenario].length, leaf.cds_length);
            ASSERT_EQ(cases[scenario].edits, leaf.edit_count);
            ASSERT_EQ(cases[scenario].sources, leaf.contributor_count);
            ASSERT_EQ(cases[scenario].flags, leaf.flags);
            ASSERT_EQ(scenario < 2u, leaf.ordered_replacements);
            ASSERT_EQ(scenario ? 2u : 1u, leaf.carriers.call_count);
            for (size_t i = 0u; i < leaf.cds_length; i++)
                ASSERT_EQ(!scenario && i == 11u ? 'C' : 'A', leaf.cds[i]);
            if (scenario < 2u) ASSERT(leaf.nominal_length_diff != (int64_t)leaf.cds_length - 12);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
        }
    }
    PASS();
}

TEST haplotype_stream_raw_record_identity_and_evidence_are_checked(void) {
    for (unsigned scenario = 0u; scenario < 6u; scenario++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)"A", (const uint8_t *)"C",
            101u, 0u, 1u, 1u, 1u, 1u, 0u, 0u};
        duckvep_raw_gt_t call;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw((const uint8_t *)"1", 1u, 1u, &call));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
        duckvep_haplotype_stream_status_t expected = DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        if (scenario == 0u) {
            source.allele_index = 2u; source.ref = (const uint8_t *)"G";
            ASSERT_EQ(expected, duckvep_haplotype_stream_begin(s, &source));
        } else if (scenario == 1u) {
            expected = DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER;
            ASSERT_EQ(expected, duckvep_haplotype_stream_begin(s, &source));
        } else if (scenario == 2u) {
            source.allele_index = UINT32_MAX; source.alt = (const uint8_t *)""; source.alt_len = 0u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
            duckvep_carrier_key_t key = {0u, 0, 1u, 2u, 0u, 0u};
            ASSERT_EQ(expected, haplotype_test_push_called(s, &key));
        } else if (scenario == 3u) {
            call.allele_index[1] = 0u; /* A retained single slot has no second allele. */
            ASSERT_EQ(expected, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
        } else if (scenario == 4u) {
            int32_t gt = 1;
            duckvep_haplotype_call_t decoded = {&gt, NULL, 0u, 1u, 1u, {0, 0u}, DUCKVEP_PHASE_STRICT};
            ASSERT_EQ(expected, duckvep_haplotype_stream_push_call(s, 0u, &decoded, NULL, 0u));
        } else {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
            expected = DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR;
            ASSERT_EQ(expected, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
            ASSERT_EQ(DUCKVEP_CARRIERS_DUPLICATE_CALL, s->carrier_error);
        }
        ASSERT_EQ(expected, duckvep_haplotype_stream_finish(s));
    }
    PASS();
}

TEST haplotype_stream_recycles_owned_alleles_and_projection_slots(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 64u);
    /* Two live transcripts suffice for all 64 inputs. Unequal allele sizes
     * force wrap padding while older events and projections are still live. */
    f.buffers.event_capacity = f.buffers.projection_capacity = 2u;
    f.buffers.allele_capacity = 8u;
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    uint32_t drained = 0u;
    for (uint32_t i = 0u; i <= 64u; i++) {
        uint8_t ref = 'A', alt[3] = {'A', 'A', 'A'};
        if (i % 2u) alt[0] = 'C';
        duckvep_haplotype_source_t source = {UINT64_MAX - i, &ref, alt,
            i < 64u ? f.starts[i] : 0u, 0u, 1u, (uint16_t)(i % 2u ? 1u : 3u), 0u, 0u, 0u, 0u};
        duckvep_haplotype_stream_status_t status;
        while ((status = i == 64u ? duckvep_haplotype_stream_finish(s)
            : haplotype_test_begin_candidates(s, &source, &i, 1u)) ==
                DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) {
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
            ASSERT_EQ(drained, leaf.carriers.transcript_index);
            ASSERT_EQ(1u, leaf.contributor_count);
            ASSERT_EQ(UINT64_MAX - drained, leaf.contributors[0].source.event_id);
            ASSERT_EQ('A', leaf.contributors[0].source.ref[0]);
            ASSERT_EQ(drained % 2u ? 'C' : 'A', leaf.contributors[0].source.alt[0]);
            const duckvep_haplotype_contributor_t *contributor = &leaf.contributors[0];
            const duckvep_haplotype_source_t *retained = &contributor->source;
            ASSERT(retained->ref >= f.alleles);
            ASSERT(retained->alt == retained->ref + retained->ref_len);
            ASSERT(retained->alt + retained->alt_len <= f.alleles + f.buffers.allele_capacity);
            duckvep_event_t prepared = {0};
            ASSERT(duckvep_event_prepare_small(retained->pos1, retained->ref, retained->ref_len,
                retained->alt, retained->alt_len, &prepared));
            prepared.chrom_id = retained->chrom_id;
            ASSERT(contributor->prepared != NULL);
            ASSERT_MEM_EQ(&prepared, contributor->prepared, sizeof prepared);
            ASSERT(contributor->projected != NULL);
            ASSERT_EQ(drained % 2u ? 1u : 2u, contributor->projected->cds_start);
            ASSERT_EQ(drained % 2u ? 1u : 0u, contributor->projected->ref_len);
            ASSERT_EQ(drained % 2u ? 1u : 2u, contributor->projected->alt_len);
            ASSERT(contributor->projected->ref == (drained % 2u ? retained->ref : NULL));
            ASSERT(contributor->projected->alt == retained->alt + prepared.alt_diff_offset);
            ASSERT_EQ(1, contributor->projected->variant_strand);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
            ASSERT_EQ(drained % 2u ? 12u : 14u, leaf.cds_length);
            uint8_t expected[14]; memset(expected, 'A', sizeof(expected));
            if (drained % 2u) expected[0] = 'C';
            ASSERT_EQ(0, memcmp(expected, leaf.cds, leaf.cds_length));
            ASSERT_EQ(4u, leaf.protein_length);
            ASSERT_EQ(0, memcmp(drained % 2u ? "QKKK" : "KKKK", leaf.protein, 4u));
            ASSERT_EQ(3u, leaf.carriers.call_count);
            uint32_t seen = 0u;
            for (uint32_t id = leaf.carriers.first_call; id;) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, id);
                ASSERT(call && call->key.sample_index < 3u);
                uint32_t sample = call->key.sample_index;
                ASSERT_EQ(0u, seen & (1u << sample)); seen |= 1u << sample;
                ASSERT_EQ(sample + 1u, call->key.lane);
                ASSERT_EQ(2u * sample + 1u, call->key.ploidy);
                ASSERT_EQ(sample % 2u, call->key.phase_set_present);
                ASSERT_EQ(sample % 2u ? -17 : 0, call->key.phase_set);
                id = call->next_leaf;
            }
            ASSERT_EQ(7u, seen);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
            drained++;
        }
        ASSERT_EQ(i == 64u ? DUCKVEP_HAPLOTYPE_STREAM_DONE : DUCKVEP_HAPLOTYPE_STREAM_OK, status);
        if (i == 64u) break;
        /* Destroy the input immediately, before any carrier is appended. */
        ref = 'N'; memset(alt, 'N', sizeof(alt));
        for (uint32_t sample = 0u; sample < 3u; sample++) {
            duckvep_carrier_key_t key = {sample, -17, (uint16_t)(sample + 1u),
                (uint16_t)(2u * sample + 1u), (uint8_t)(sample % 2u), 0u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
    }
    ASSERT_EQ(64u, drained); ASSERT_EQ(64u, s->input_events);
    ASSERT_EQ(64u, s->projected_events); ASSERT_EQ(64u, s->completed_leaves);
    ASSERT_EQ(64u * 13u, s->translated_bases);
    ASSERT_EQ(2u, s->peak_events); ASSERT_EQ(2u, s->peak_projections);
    ASSERT_EQ(8u, s->peak_alleles);
    ASSERT_EQ(0u, s->event_count); ASSERT_EQ(0u, s->projection_count);
    ASSERT_EQ(0u, s->allele_count); ASSERT_EQ(0u, s->carriers.call_count);
    PASS();
}

TEST haplotype_stream_projects_once_and_preserves_occupied_ancestors(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 2u);
    f.starts[1] = f.starts[0]; f.ends[1] = f.ends[0];
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    uint32_t candidates[] = {0u, 1u};
    /* An earlier upload can project downstream of the next uploaded allele
     * after unchanged REF/ALT prefixes are trimmed. Preserve the raw alleles. */
    duckvep_haplotype_source_t source = {200u, (const uint8_t *)"AAAAAA", (const uint8_t *)"AAAAAC",
        100u, 0u, 6u, 6u, 0u, 0u, 0u, 0u};
    duckvep_carrier_key_t keys[] = {{0u, 0, 1u, 2u, 0u, 0u}, {42u, 0, 2u, 4u, 1u, 0u}};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, candidates, 2u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &keys[0]));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &keys[1]));
    source.event_id++; source.pos1 = 102u; source.alt = (const uint8_t *)"G";
    source.ref = (const uint8_t *)"A"; source.ref_len = source.alt_len = 1u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, candidates, 2u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &keys[0]));
    uint32_t paths = 0u;
    duckvep_haplotype_stream_status_t status;
    while ((status = duckvep_haplotype_stream_finish(s)) == DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) {
        duckvep_haplotype_leaf_t leaf;
        while ((status = duckvep_haplotype_stream_next(s, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
            ASSERT_EQ(1u, leaf.carriers.call_count);
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
            ASSERT(call); ASSERT(leaf.carriers.transcript_index < 2u);
            uint32_t bit = leaf.carriers.transcript_index * 2u + (call->key.sample_index == 42u);
            ASSERT_EQ(0u, paths & (1u << bit)); paths |= 1u << bit;
            size_t count = call->key.sample_index == 42u ? 1u : 2u;
            ASSERT_EQ(count, leaf.contributor_count);
            ASSERT_EQ(200u, leaf.contributors[0].source.event_id);
            ASSERT_EQ(6u, leaf.contributors[0].source.alt_len);
            ASSERT_EQ(0, memcmp("AAAAAC", leaf.contributors[0].source.alt, 6u));
            if (count == 2u) ASSERT_EQ(201u, leaf.contributors[1].source.event_id);
            uint8_t expected[12]; memset(expected, 'A', sizeof(expected));
            expected[5] = 'C'; if (count == 2u) expected[2] = 'G';
            ASSERT_EQ(12u, leaf.cds_length);
            ASSERT_EQ(0, memcmp(expected, leaf.cds, sizeof(expected)));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
    ASSERT_EQ(15u, paths); ASSERT_EQ(4u, s->projected_events);
    ASSERT_EQ(4u, s->completed_leaves); ASSERT_EQ(48u, s->translated_bases);
    PASS();
}

TEST haplotype_stream_keeps_noncoding_contributors_without_poisoning_cds(void) {
    /* Two exons with three UTR bases at each transcript end. Every small
     * genomic span is classified independently by its literal CDS-base set;
     * strict replay withholds crossing spans, while raw compatibility retains
     * their omission as conditional evidence and replays the mapped source. */
    for (unsigned raw = 0u; raw < 2u; raw++) {
        for (int strand = -1; strand <= 1; strand += 2) {
            for (uint32_t pos = 100u; pos <= 151u; pos++) {
                for (uint16_t n = 1u; n <= 52u && pos + n - 1u <= 151u; n++) {
                    if (pos <= 105u && pos + n - 1u >= 105u) continue;
                    struct haplotype_stream_scene f;
                    haplotype_stream_scene_prepare(&f, 1u);
                    uint32_t cds_start = 103u, cds_end = 148u;
                    uint32_t exon_start[] = {100u, 143u}, exon_end[] = {108u, 151u};
                    uint32_t cdna_start[] = {1u, 10u}, cdna_end[] = {9u, 18u};
                    if (strand < 0) {
                        exon_start[0] = 143u; exon_start[1] = 100u;
                        exon_end[0] = 151u; exon_end[1] = 108u;
                    }
                    f.ends[0] = 151u; f.strands[0] = (int8_t)strand;
                    f.model.cds_start1 = &cds_start; f.model.cds_end1 = &cds_end;
                    f.exon_count[0] = 2u;
                    f.exons = (duckvep_exon_model_t){.exon_count = 2u,
                        .start1 = exon_start, .end1 = exon_end,
                        .cdna_start1 = cdna_start, .cdna_end1 = cdna_end};
                    memset(f.reference, strand > 0 ? 'A' : 'T', sizeof(f.reference));
                    duckvep_haplotype_stream_t *s = &f.stream;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                        s, &f.model, &f.exons, &f.sequences, &f.buffers));
                    uint8_t reference[52], alternate[52];
                    memset(reference, 'A', sizeof(reference));
                    memset(alternate, 'C', sizeof(alternate));
                    duckvep_haplotype_source_t events[] = {
                        {1u, reference, alternate, pos, 0u, n, n, raw, (uint8_t)raw, 0u, 0u},
                        {2u, (const uint8_t *)"A", (const uint8_t *)"G", 105u, 0u, 1u, 1u,
                            raw, (uint8_t)raw, 0u, 0u}
                    };
                    if (pos > 105u) {
                        duckvep_haplotype_source_t swap = events[0];
                        events[0] = events[1]; events[1] = swap;
                    }
                    duckvep_carrier_key_t key = {0u, 0, 1u, 1u, 0u, 0u};
                    uint32_t tx = 0u;
                    for (size_t i = 0u; i < 2u; i++) {
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                            haplotype_test_begin_candidates(s, events + i, &tx, 1u));
                        if (raw) {
                            duckvep_raw_gt_t call;
                            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                                (const uint8_t *)"1|1", 3u, 1u, &call));
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
                        } else {
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                haplotype_test_push_called(s, &key));
                        }
                    }
                    uint32_t coding_bases = 0u;
                    uint8_t expected[12];
                    memcpy(expected, f.reference, sizeof(expected));
                    expected[strand > 0 ? 2u : 9u] = strand > 0 ? 'G' : 'C';
                    for (uint32_t p = pos; p < pos + n; p++) {
                        if ((p >= 103u && p <= 108u) || (p >= 143u && p <= 148u)) {
                            uint32_t offset = p <= 108u ? p - 103u : p - 143u + 6u;
                            expected[strand > 0 ? offset : 11u - offset] = strand > 0 ? 'C' : 'G';
                            coding_bases++;
                        }
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
                        duckvep_haplotype_stream_finish(s));
                    duckvep_haplotype_leaf_t leaf;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
                    ASSERT_EQ(2u, leaf.contributor_count);
                    ASSERT_EQ(events[0].event_id, leaf.contributors[0].source.event_id);
                    ASSERT_EQ(events[1].event_id, leaf.contributors[1].source.event_id);
                    int omitted = raw && coding_bases && coding_bases != n;
                    ASSERT_EQ(coding_bases == n ? DUCKVEP_CDS_EDIT_OK : omitted
                        ? DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED : DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                        leaf.contributors[pos < 105u ? 0u : 1u].projection_status);
                    if (coding_bases == 0u || coding_bases == n || omitted) {
                        if (omitted) {
                            memcpy(expected, f.reference, sizeof(expected));
                            expected[strand > 0 ? 2u : 9u] = strand > 0 ? 'G' : 'C';
                            ASSERT(leaf.contributors[pos < 105u ? 0u : 1u].evidence_flags &
                                DUCKVEP_CARRIER_CONDITIONAL);
                        }
                        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                        ASSERT_EQ(omitted ? DUCKVEP_HAPLOTYPE_CONDITIONAL : DUCKVEP_HAPLOTYPE_OK,
                            leaf.sequence_status);
                        ASSERT_EQ(12u, leaf.cds_length);
                        ASSERT_EQ(0, memcmp(expected, leaf.cds, sizeof(expected)));
                        ASSERT_EQ(coding_bases && !omitted ? 2u : 1u, leaf.edit_count);
                    } else {
                        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS, leaf.projection_status);
                        ASSERT(leaf.cds == NULL && leaf.protein == NULL && leaf.blocks == NULL);
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
                }
            }
        }
    }
    PASS();
}

TEST haplotype_source_omission_checks_layout_sequence_and_ref(void) {
    for (int strand = -1; strand <= 1; strand += 2) {
        for (unsigned scenario = 0u; scenario < 13u; scenario++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            uint32_t starts[] = {100u, 146u}, ends[] = {105u, 151u};
            uint32_t cs[] = {1u, 7u}, ce[] = {6u, 12u};
            if (strand < 0) { starts[0] = 146u; starts[1] = 100u; ends[0] = 151u; ends[1] = 105u; }
            f.ends[0] = 151u; f.strands[0] = (int8_t)strand; f.exon_count[0] = 2u;
            f.exons = (duckvep_exon_model_t){.exon_count = 2u,
                .start1 = starts, .end1 = ends, .cdna_start1 = cs, .cdna_end1 = ce};
            memset(f.reference, strand > 0 ? 'A' : 'T', sizeof(f.reference));
            uint8_t ref[44]; memset(ref, 'A', sizeof(ref));
            uint8_t alt = 'C';
            duckvep_cds_edit_status_t expected = DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED;
            if (scenario == 1u || scenario == 2u) {
                ref[scenario == 1u ? 0u : 43u] = 'C'; expected = DUCKVEP_CDS_EDIT_REF_MISMATCH;
            }
            if (scenario == 3u) { alt = 'Z'; expected = DUCKVEP_CDS_EDIT_INVALID_ALLELE; }
            if (scenario == 4u) { f.sequence_offsets[0] = 13u; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 5u) { f.lengths[0] = 11u; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 6u) { ce[0]--; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 7u) { ends[0]++; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            if (scenario == 8u) { f.exons.cdna_start1 = NULL; expected = DUCKVEP_CDS_EDIT_INVALID_ARG; }
            uint32_t cached_start = 2u, cached_end = 13u, cached_exon = 0u;
            uint8_t cached_phase = 0u;
            if (scenario == 9u) {
                f.model.cds_cdna_start1 = &cached_start; f.model.cds_cdna_end1 = &cached_end;
                f.model.cds_start_exon_index = &cached_exon; f.model.cds_phase_offset = &cached_phase;
                expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
            }
            if (scenario == 10u) { f.lengths[0] = 0u; expected = DUCKVEP_CDS_EDIT_OUT_OF_CDS; }
            if (scenario == 11u) ref[20] = 'C'; /* Intronic REF is not represented in this CDS pool. */
            if (scenario == 12u) { ref[20] = 'N'; expected = DUCKVEP_CDS_EDIT_INVALID_ALLELE; }
            duckvep_event_t event;
            ASSERT(duckvep_event_prepare_replacement(104u, ref, sizeof(ref), &alt, 1u, &event));
            duckvep_prepared_cds_allele_t allele = {&event, ref, &alt, ref, sizeof(ref), 1u, 1};
            duckvep_haplotype_edit_t edit;
            duckvep_cds_edit_status_t status = duckvep_compat_vep_source_cds_edit_build(
                &f.model, &f.exons, &f.sequences, 0u, (int8_t)strand, &allele, &edit);
            if (status != expected) fprintf(stderr, "source mapping strand=%d scenario=%u\n", strand, scenario);
            ASSERT_EQ_FMT(expected, status, "%d");
            if (expected == DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED) {
                duckvep_haplotype_edit_t zero = {0};
                ASSERT_MEM_EQ(&zero, &edit, sizeof(edit));
            }
            if (scenario == 5u) {
                /* Intronic context cannot conceal an inconsistent CDS extent. */
                ASSERT(duckvep_event_prepare_replacement(120u, ref, 1u, &alt, 1u, &event));
                allele.ref_length = allele.alt_length = 1u;
                ASSERT_EQ(DUCKVEP_CDS_EDIT_INVALID_ARG, duckvep_compat_vep_source_cds_edit_build(
                    &f.model, &f.exons, &f.sequences, 0u, (int8_t)strand, &allele, &edit));
            }
        }
    }
    PASS();
}

TEST haplotype_source_skips_require_valid_model_ref_and_literal_geometry(void) {
    /* Haplosaurus _mutate_sequences tests raw /^[ACGT]*$/ before mutation.
     * Actual c,N,n,U,u observations: receipt
     * 286298b57d32972e0c9094ddfcdf98c7f0aa4e8d726830e51fcbb6ad9521a175.
     * Other lowercase bases and mixed strings exercise the same literal rule.
     * Native unsupported R/dash contrasts below do not claim upstream parity. */
    const char *skipped_alts[] = {"N", "n", "U", "u", "a", "c", "g", "t", "AN", "Ac"};
    enum {
        VALID_SOURCE, WRONG_REF, WRONG_POOL_OFFSET, WRONG_CDS_LENGTH,
        WRONG_EXON_CDNA, WRONG_CACHED_CDS, WRONG_RAW_START,
        WRONG_FEATURE_OFFSET, WRONG_REF_LENGTH, WRONG_VARIANT_STRAND,
        SOURCE_SCENARIO_COUNT
    };
    unsigned checked = 0u;
    for (int strand = -1; strand <= 1; strand += 2) {
        for (size_t i = 0u; i < sizeof(skipped_alts) / sizeof(skipped_alts[0]); i++) {
            for (unsigned scenario = 0u; scenario < SOURCE_SCENARIO_COUNT; scenario++) {
                struct haplotype_stream_scene f;
                haplotype_stream_scene_prepare(&f, 1u);
                f.starts[0] = 11u;
                f.ends[0] = 22u;
                f.strands[0] = (int8_t)strand;
                memcpy(f.reference, "ATGAAAGCCTAA", 12u);
                uint8_t ref = strand > 0 ? 'A' : 'T';
                if (scenario == WRONG_REF) ref = 'N'; /* N is not a wildcard. */
                const uint8_t *alt = (const uint8_t *)skipped_alts[i];
                uint16_t alt_length = (uint16_t)strlen(skipped_alts[i]);
                duckvep_event_t event;
                ASSERT(duckvep_event_prepare_replacement(strand > 0 ? 14u : 19u,
                    &ref, 1u, alt, alt_length, &event));
                duckvep_prepared_cds_allele_t allele = {
                    &event, &ref, alt, &ref, 1u, alt_length, 1};
                duckvep_cds_edit_status_t expected = DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED;
                uint32_t cached_start = 2u, cached_end = 13u, cached_exon = 0u;
                uint8_t cached_phase = 0u;
                switch (scenario) {
                case VALID_SOURCE:
                    break;
                case WRONG_REF:
                    expected = DUCKVEP_CDS_EDIT_REF_MISMATCH;
                    break;
                case WRONG_POOL_OFFSET:
                    f.sequence_offsets[0] = 13u;
                    expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
                    break;
                case WRONG_CDS_LENGTH:
                    f.lengths[0] = 11u;
                    expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
                    break;
                case WRONG_EXON_CDNA:
                    f.cdna_ends[0] = 11u;
                    expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
                    break;
                case WRONG_CACHED_CDS:
                    f.model.cds_cdna_start1 = &cached_start;
                    f.model.cds_cdna_end1 = &cached_end;
                    f.model.cds_start_exon_index = &cached_exon;
                    f.model.cds_phase_offset = &cached_phase;
                    expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
                    break;
                case WRONG_RAW_START:
                    event.raw_start1++;
                    expected = DUCKVEP_CDS_EDIT_INVALID_EVENT;
                    break;
                case WRONG_FEATURE_OFFSET:
                    event.feature_allele_offset = 1u;
                    expected = DUCKVEP_CDS_EDIT_INVALID_EVENT;
                    break;
                case WRONG_REF_LENGTH:
                    event.ref_diff_length = 2u;
                    expected = DUCKVEP_CDS_EDIT_INVALID_EVENT;
                    break;
                case WRONG_VARIANT_STRAND:
                    allele.variant_strand = -1;
                    expected = DUCKVEP_CDS_EDIT_INVALID_ARG;
                    break;
                }
                duckvep_haplotype_edit_t edit, zero = {0};
                memset(&edit, 0xff, sizeof(edit));
                duckvep_cds_edit_status_t status = duckvep_compat_vep_source_cds_edit_build(
                    &f.model, &f.exons, &f.sequences, 0u, (int8_t)strand, &allele, &edit);
                if (status != expected)
                    fprintf(stderr, "skipped ALT=%s strand=%d scenario=%u\n",
                        skipped_alts[i], strand, scenario);
                ASSERT_EQ_FMT(expected, status, "%d");
                ASSERT_MEM_EQ(&zero, &edit, sizeof(edit));
                ASSERT_MEM_EQ("ATGAAAGCCTAA", f.reference, 12u);
                checked++;
            }
        }
    }
    ASSERT_EQ(200u, checked);
    const char *unsupported_alts[] = {"R", "-", "N-", "NR"};
    for (size_t i = 0u; i < sizeof(unsupported_alts) / sizeof(unsupported_alts[0]); i++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        const uint8_t *alt = (const uint8_t *)unsupported_alts[i];
        duckvep_event_t event;
        ASSERT(duckvep_event_prepare_replacement(100u, (const uint8_t *)"A", 1u,
            alt, (uint16_t)strlen(unsupported_alts[i]), &event));
        duckvep_prepared_cds_allele_t allele = {&event, (const uint8_t *)"A", alt,
            (const uint8_t *)"A", 1u, (uint16_t)strlen(unsupported_alts[i]), 1};
        duckvep_haplotype_edit_t edit, zero = {0};
        memset(&edit, 0xff, sizeof(edit));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_INVALID_ALLELE, duckvep_compat_vep_source_cds_edit_build(
            &f.model, &f.exons, &f.sequences, 0u, 1, &allele, &edit));
        ASSERT_MEM_EQ(&zero, &edit, sizeof(edit));
    }
    PASS();
}

TEST haplotype_stream_decoded_calls_do_not_adopt_raw_allele_skipping(void) {
    const duckvep_phase_policy_t policies[] = {
        DUCKVEP_PHASE_STRICT, DUCKVEP_PHASE_VEP_COMPAT};
    for (size_t i = 0u; i < sizeof(policies) / sizeof(policies[0]); i++) {
        for (unsigned invalid_payload = 0u; invalid_payload < 2u; invalid_payload++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.starts[0] = 11u;
            f.ends[0] = 22u;
            memcpy(f.reference, "ATGGCNGCCTAA", 12u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            /* The independent feature removes the matching N anchor. A changed
             * ALT N still fails; it must not become a raw conditional no-op.
             * Decoded source metadata is zero; ALT index belongs to the call. */
            duckvep_haplotype_source_t source = {.event_id = 143856u,
                .pos1 = invalid_payload ? 15u : 16u,
                .ref = (const uint8_t *)(invalid_payload ? "C" : "N"),
                .alt = (const uint8_t *)(invalid_payload ? "N" : "NGCC"),
                .ref_len = 1u, .alt_len = invalid_payload ? 1u : 4u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_begin(&f.stream, &source));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_project(&f.stream, 0u));
            int32_t allele = 1;
            uint8_t phased = 1u;
            duckvep_haplotype_phase_set_t set = {0};
            duckvep_haplotype_call_t call = {&allele, &phased, 7u, 1u, 1u, set, policies[i]};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_call(&f.stream, 0u, &call, &set, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
                duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_next(&f.stream, &leaf));
            ASSERT_EQ(1u, leaf.carriers.call_count);
            ASSERT_EQ(1u, leaf.contributor_count);
            ASSERT_EQ(DUCKVEP_CARRIER_CALLED, leaf.evidence_flags);
            ASSERT_EQ(invalid_payload ? DUCKVEP_CDS_EDIT_INVALID_ALLELE : DUCKVEP_CDS_EDIT_OK,
                leaf.projection_status);
            ASSERT_EQ(leaf.projection_status, leaf.contributors[0].projection_status);
            if (invalid_payload) {
                ASSERT(leaf.cds == NULL);
                ASSERT(leaf.protein == NULL);
                ASSERT_EQ(0u, leaf.edit_count);
                ASSERT_EQ(0u, leaf.block_count);
            } else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
                ASSERT_EQ(15u, leaf.cds_length);
                ASSERT_MEM_EQ("ATGGCNGCCGCCTAA", leaf.cds, 15u);
                ASSERT_EQ(5u, leaf.protein_length);
                ASSERT_MEM_EQ("MAAA*", leaf.protein, 5u);
                ASSERT_EQ(1u, leaf.edit_count);
                ASSERT_EQ(1u, leaf.block_count);
                ASSERT_EQ(143856u, leaf.edit_event_ids[0]);
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE,
                duckvep_haplotype_stream_next(&f.stream, &leaf));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(&f.stream));
        }
    }
    PASS();
}

TEST haplotype_stream_unmapped_sources_retain_raw_slot_evidence(void) {
    const char *spellings[] = {"1|1", "0|1", ".|1", ".", "0|0"};
    for (int strand = -1; strand <= 1; strand += 2) {
        for (unsigned mode = 0u; mode < 5u; mode++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            uint32_t starts[] = {100u, 146u}, ends[] = {105u, 151u};
            uint32_t cs[] = {1u, 7u}, ce[] = {6u, 12u};
            if (strand < 0) { starts[0] = 146u; starts[1] = 100u; ends[0] = 151u; ends[1] = 105u; }
            f.ends[0] = 151u; f.strands[0] = (int8_t)strand; f.exon_count[0] = 2u;
            f.exons = (duckvep_exon_model_t){.exon_count = 2u,
                .start1 = starts, .end1 = ends, .cdna_start1 = cs, .cdna_end1 = ce};
            memset(f.reference, strand > 0 ? 'A' : 'T', sizeof(f.reference));
            /* Three interpretations of one record occupy 44+44+44+1+44 bytes. */
            uint8_t allele_storage[256];
            f.buffers.alleles = allele_storage; f.buffers.allele_capacity = sizeof(allele_storage);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            uint8_t ref[44]; memset(ref, 'A', sizeof(ref));
            duckvep_raw_gt_t call, anchor;
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                (const uint8_t *)spellings[mode], strlen(spellings[mode]), 1u, &call));
            ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                (const uint8_t *)"1|1", 3u, 1u, &anchor));
            for (unsigned interpretation = 0u; interpretation < 4u; interpretation++) {
                duckvep_haplotype_source_t source = {1u, ref,
                    interpretation == 0u ? ref : (const uint8_t *)(interpretation == 1u ? "C" : ""),
                    104u, 0u, sizeof(ref), interpretation == 0u ? sizeof(ref) : interpretation == 1u ? 1u : 0u,
                    interpretation < 2u ? interpretation : UINT32_MAX, 1u, 0u, 0u};
                if (interpretation == 3u) source = (duckvep_haplotype_source_t){
                    2u, (const uint8_t *)"A", (const uint8_t *)"G", 150u, 0u, 1u, 1u, 1u, 1u, 0u, 0u};
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                    &f.stream, 0u, 0u, interpretation == 3u ? &anchor : &call, 1u));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            unsigned carriers = 0u;
            duckvep_haplotype_stream_status_t status;
            while ((status = duckvep_haplotype_stream_next(&f.stream, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                ASSERT_EQ(mode == 4u ? DUCKVEP_HAPLOTYPE_OK : DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
                ASSERT_EQ(12u, leaf.cds_length);
                const char *expected_cds = strand > 0 ? "AAAAAAAAAAGA" : "TCTTTTTTTTTT";
                ASSERT_MEM_EQ(expected_cds, leaf.cds, 12u);
                ASSERT_EQ(1u, leaf.edit_count); ASSERT_EQ(2u, leaf.edit_event_ids[0]);
                ASSERT_EQ(mode == 4u ? 1u : 2u, leaf.contributor_count);
                if (mode != 4u) {
                    ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
                    ASSERT_EQ(DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED, leaf.contributors[0].projection_status);
                    ASSERT(leaf.contributors[0].evidence_flags & DUCKVEP_CARRIER_CONDITIONAL);
                }
                carriers += leaf.carriers.call_count;
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status); ASSERT_EQ(2u, carriers);
        }
    }
    PASS();
}

TEST haplotype_stream_reference_only_protein_uses_call_retention(void) {
    const char *references[] = {"ATGAAACCCTAA", "ATGAAATAACCC", "CTGAAACCCTAA", "ATGGCNCCCTAA"};
    const char *curated[] = {"MKP*", "MK*P", "MKP*", "MAP*"};
    const char *raw[] = {"MKP*", "MK*", "LKP*", "MAP*"};
    const char *alternate[] = {"MQP*", "MQ*", "LQP*", "MPP*"};
    const char *coding[] = {"MKP*", "MK*P", "LKP*", "MAP*"};
    for (unsigned reference = 0u; reference < 4u; reference++) {
        for (int strand = -1; strand <= 1; strand += 2) {
            for (unsigned peptide_edit = 0u; peptide_edit < 2u; peptide_edit++) {
                /* Missing-only, coding, shadowed, unmapped, UTR and intronic calls. */
                for (unsigned mode = 0u; mode < 7u; mode++) {
                    struct haplotype_stream_scene f;
                    haplotype_stream_scene_prepare(&f, 1u);
                    memcpy(f.reference, references[reference], 12u);
                    f.ends[0] = 117u; f.cdna_ends[0] = 18u; f.strands[0] = (int8_t)strand;
                    uint32_t cds_start = 103u, cds_end = 114u;
                    f.model.cds_start1 = &cds_start; f.model.cds_end1 = &cds_end;
                    uint32_t exon_starts[] = {100u, mode == 5u ? 121u : 133u};
                    uint32_t exon_ends[] = {108u, exon_starts[1] + 8u};
                    uint32_t cdna_starts[] = {1u, 10u}, cdna_ends[] = {9u, 18u};
                    if (mode >= 5u) {
                        f.ends[0] = exon_ends[1]; cds_end = exon_starts[1] + 5u;
                        if (strand < 0) {
                            exon_starts[0] = exon_starts[1]; exon_starts[1] = 100u;
                            exon_ends[0] = exon_ends[1]; exon_ends[1] = 108u;
                        }
                        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
                        f.exons.start1 = exon_starts; f.exons.end1 = exon_ends;
                        f.exons.cdna_start1 = cdna_starts; f.exons.cdna_end1 = cdna_ends;
                    }
                    uint32_t peptide_offsets[] = {0u, peptide_edit}, peptide_position = 2u;
                    uint8_t peptide_alt = 'W';
                    f.sequences.peptide_edit_offset = peptide_offsets;
                    f.sequences.peptide_edit_position1 = &peptide_position;
                    f.sequences.peptide_edit_alt = &peptide_alt;
                    f.sequences.peptide_edit_count = peptide_edit;
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
                    duckvep_raw_gt_t missing, call;
                    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                        (const uint8_t *)".", 1u, 1u, &missing));
                    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
                        (const uint8_t *)(mode ? "0|1" : "."), mode ? 3u : 1u, 1u, &call));
                    duckvep_haplotype_source_t source = {1u, (const uint8_t *)"C",
                        (const uint8_t *)"C", 100u, 0u, 1u, 1u, 0u, 1u, 0u, 0u};
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                        duckvep_haplotype_stream_push_raw_call(&f.stream, 0u, 0u, &missing, 1u));
                    uint8_t ref[3], alt[3];
                    uint32_t pos = mode >= 5u ? 115u : mode == 3u ? 102u
                        : mode == 4u ? 101u : strand > 0 ? 106u : 111u;
                    uint16_t length = mode == 3u ? 3u : 1u;
                    for (uint16_t i = 0u; i < length; i++) {
                        uint32_t p = pos + i;
                        ref[i] = mode >= 5u || p < cds_start ? 'C' : strand > 0 ? f.reference[p - cds_start]
                            : kprop_hgvs_complement(f.reference[cds_end - p]);
                        alt[i] = mode == 4u ? 'G' : strand > 0 ? 'C' : 'G';
                    }
                    for (uint32_t allele = 0u; allele < 2u; allele++) {
                        source = (duckvep_haplotype_source_t){2u, ref, allele ? alt : ref,
                            pos, 0u, length, length, allele, 1u, 0u, 0u};
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &source));
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
                            &f.stream, 0u, 0u, &call, mode != 2u));
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
                        duckvep_haplotype_stream_finish(&f.stream));
                    duckvep_haplotype_leaf_t leaf;
                    unsigned seen = 0u;
                    duckvep_haplotype_stream_status_t status;
                    char expected_reference[5]; memcpy(expected_reference, curated[reference], 5u);
                    if (peptide_edit) expected_reference[1] = 'W';
                    while ((status = duckvep_haplotype_stream_next(&f.stream, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_CONDITIONAL, leaf.sequence_status);
                        ASSERT_EQ(4u, leaf.reference_protein_length);
                        ASSERT_MEM_EQ(expected_reference, leaf.reference_protein, 4u);
                        ASSERT_MEM_EQ(coding[reference], leaf.reference_coding_protein, 5u);
                        ASSERT_EQ(4u, leaf.reference_coding_translation.length);
                        ASSERT_EQ(reference == 1u ? 3u : 4u,
                            leaf.reference_coding_translation.first_stop_position1);
                        ASSERT_EQ(reference != 3u, leaf.reference_coding_translation.unambiguous);
                        for (uint32_t id = leaf.carriers.first_call; id;) {
                            const duckvep_carrier_call_t *carrier = duckvep_carriers_call(&f.stream.carriers, id);
                            ASSERT(carrier); ASSERT(carrier->key.lane == 1u || carrier->key.lane == 2u);
                            unsigned bit = 1u << (carrier->key.lane - 1u);
                            ASSERT_EQ(0u, seen & bit); seen |= bit;
                            int edited = mode == 1u && carrier->key.lane == 2u;
                            int reference_route = !mode || mode >= 5u;
                            const char *expected = reference_route ? expected_reference
                                : edited ? alternate[reference] : raw[reference];
                            ASSERT_EQ(strlen(expected), leaf.protein_length);
                            ASSERT_MEM_EQ(expected, leaf.protein, leaf.protein_length);
                            ASSERT_EQ((size_t)edited, leaf.edit_count);
                            ASSERT_EQ(reference_route, leaf.protein == leaf.reference_protein);
                            if (reference_route) ASSERT_EQ(0u, leaf.flags & DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED);
                            id = carrier->next_leaf;
                        }
                    }
                    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status); ASSERT_EQ(3u, seen);
                }
            }
        }
    }
    PASS();
}

TEST haplotype_stream_retains_conflicts_and_latches_resource_errors(void) {
    struct haplotype_stream_scene f;
    duckvep_carrier_key_t key = {0u, 1, 1u, 2u, 1u, 0u};
    uint32_t tx = 0u;
    /* One bad REF, two conflicting edits, and each named capacity limit. */
    for (unsigned scenario = 0u; scenario < 13u; scenario++) {
        haplotype_stream_scene_prepare(&f, 1u);
        uint32_t peptide_offsets[] = {1u, 0u};
        if (scenario == 2u) f.buffers.event_capacity = 1u;
        if (scenario == 3u) f.buffers.projection_capacity = 1u;
        if (scenario == 4u) f.buffers.allele_capacity = 3u;
        if (scenario == 5u) f.buffers.leaf_capacity = 1u;
        if (scenario == 6u) f.buffers.cds_capacity = 11u;
        if (scenario == 7u) f.buffers.protein_capacity = 4u;
        if (scenario == 9u) f.buffers.edit_capacity = 1u;
        if (scenario == 10u) f.buffers.reference_protein_capacity = 5u;
        if (scenario == 11u) f.sequences.peptide_edit_offset = peptide_offsets;
        if (scenario == 12u) f.sequences.peptide_edit_count = 1u;
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)(scenario ? "A" : "C"),
            (const uint8_t *)"G", 100u, 0u, 1u, 1u, 0u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        duckvep_haplotype_stream_status_t status = DUCKVEP_HAPLOTYPE_STREAM_OK;
        if ((scenario >= 1u && scenario <= 5u) || scenario == 9u) {
            source.event_id = 2u; source.alt = (const uint8_t *)"C";
            if (scenario != 1u) source.pos1++;
            status = haplotype_test_begin_candidates(s, &source, &tx, 1u);
            if (scenario == 2u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL, status);
            else if (scenario == 3u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL, status);
            else if (scenario == 4u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL, status);
            else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, status);
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
            }
        } else if (scenario == 8u) {
            status = haplotype_test_begin_candidates(s, &source, &tx, 1u);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER, status);
        }
        duckvep_haplotype_leaf_t leaf;
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
            status = duckvep_haplotype_stream_next(s, &leaf);
            if (scenario == 5u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL, status);
            else if (scenario == 9u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL, status);
            else if (scenario == 6u || scenario == 7u || scenario == 10u)
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, status);
            else if (scenario >= 11u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG, status);
            else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, status);
                ASSERT_EQ(scenario == 0u ? DUCKVEP_CDS_EDIT_REF_MISMATCH : DUCKVEP_CDS_EDIT_OK,
                          leaf.projection_status);
                if (scenario == 1u) ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER, leaf.sequence_status);
                ASSERT_EQ(scenario == 0u ? 1u : 2u, leaf.contributor_count);
                ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
                if (scenario == 1u) ASSERT_EQ(2u, leaf.contributors[1].source.event_id);
                ASSERT(!leaf.cds && !leaf.protein);
                ASSERT_EQ(0u, leaf.cds_length); ASSERT_EQ(0u, leaf.protein_length);
                ASSERT_EQ(0, leaf.nominal_length_diff);
                ASSERT_EQ(1u, leaf.carriers.call_count);
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
            }
        }
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) {
            ASSERT_EQ(status, haplotype_test_push_called(s, &key));
            ASSERT_EQ(status, duckvep_haplotype_stream_finish(s));
            memset(&leaf, 0xa5, sizeof(leaf));
            ASSERT_EQ(status, duckvep_haplotype_stream_next(s, &leaf));
            duckvep_haplotype_leaf_t zero_leaf = {0};
            ASSERT_EQ(0, memcmp(&leaf, &zero_leaf, sizeof(leaf)));
        }
    }
    PASS();
}

TEST haplotype_stream_strict_overlap_union_is_reference_indexed_and_transactional(void) {
    for (unsigned strand = 0u; strand < 2u; strand++) for (unsigned scenario = 0u; scenario < 10u; scenario++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, scenario == 6u ? "atggctgaataa" : "ATGGCTGAATAA", 12u);
        f.strands[0] = strand ? -1 : 1;
        duckvep_haplotype_edit_t joint[8];
        uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = scenario == 3u ? 11u : 12u;
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_carrier_key_t key = {0u, 10, 1u, 2u, 1u, 0u};
        if (scenario == 7u) { key.phase_set = 0; key.phase_set_present = 0u; }
        uint32_t tx = 0u;
        duckvep_haplotype_source_t sources[2] = {
            {1u, (const uint8_t *)"GC", (const uint8_t *)"AT",
                strand ? 107u : 103u, 0u, 2u, 2u, 0u, 0u, 0u, 0u},
            {2u, (const uint8_t *)(strand ? "AG" : "CT"),
                (const uint8_t *)(scenario == 1u ? (strand ? "AC" : "GT") : (strand ? "AA" : "TT")),
                strand ? 106u : 104u, 0u, 2u, 2u, 0u, 0u, 0u, 0u}};
        if (scenario == 2u) { sources[1] = sources[0]; sources[1].event_id = 2u; }
        if (scenario == 4u) {
            sources[0].ref = (const uint8_t *)(strand ? "ATTCAGC" : "GCTGAAT");
            sources[0].alt = (const uint8_t *)(strand ? "GTTCAAT" : "ATTGAAC");
            sources[0].ref_len = sources[0].alt_len = 7u;
            sources[0].pos1 = strand ? 102u : 103u;
        }
        if (scenario == 5u) sources[1].alt = (const uint8_t *)(strand ? "TG" : "CA");
        if (scenario == 8u) {
            sources[0].ref = (const uint8_t *)(strand ? "AGC" : "GCT");
            sources[0].alt = (const uint8_t *)(strand ? "TGT" : "ACA");
            sources[0].ref_len = sources[0].alt_len = 3u;
            sources[0].pos1 = strand ? 106u : 103u;
            sources[1].alt = (const uint8_t *)"TA";
        }
        if (scenario == 9u) sources[0].alt = (const uint8_t *)(strand ? "GT" : "AC");
        for (unsigned i = 0u; i < 2u; i++) {
            unsigned at = strand && scenario != 4u ? 1u - i : i;
            duckvep_haplotype_source_t source = sources[at];
            source.event_id = i + 1u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        if (scenario == 3u) {
            memset(&leaf, 0xa5, sizeof(leaf));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, duckvep_haplotype_stream_next(s, &leaf));
            duckvep_haplotype_leaf_t zero = {0};
            ASSERT_EQ(0, memcmp(&zero, &leaf, sizeof(leaf)));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, duckvep_haplotype_stream_next(s, &leaf));
            f.buffers.joint_capacity = 12u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f.model, &f.exons, &f.sequences, &f.buffers));
            for (unsigned i = 0u; i < 2u; i++) {
                duckvep_haplotype_source_t source = sources[strand ? 1u - i : i];
                source.event_id = i + 1u;
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
            }
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(2u, leaf.contributor_count);
        if (scenario == 1u) {
            ASSERT(!leaf.cds && !leaf.protein && !leaf.joint_edits);
            ASSERT_EQ(DUCKVEP_REASON_CONTRADICTORY_EDITS, leaf.path_reason);
            ASSERT_EQ(0u, leaf.haplotype_so_mask); ASSERT(!leaf.nmd_stop_valid && !leaf.nmd_junction_valid);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_NMD_UNKNOWN, leaf.nmd);
        } else if (scenario == 4u) {
            ASSERT_EQ(0, memcmp("ATGATTGAACAA", leaf.cds, 12u));
            ASSERT_EQ(2u, leaf.joint_edit_count); ASSERT_EQ(2u, leaf.block_count);
            ASSERT_EQ(3u, leaf.edit_count);
        } else if (scenario == 5u || scenario == 8u || scenario == 9u) {
            ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
            ASSERT_EQ(0, memcmp(scenario == 9u ? "ATGATTGAATAA" : "ATGATAGAATAA", leaf.cds, 12u));
            ASSERT_EQ(0, memcmp("MIE*", leaf.protein, 4u));
            ASSERT_EQ(scenario == 8u ? 3u : 2u, leaf.edit_count);
            if (scenario == 8u) {
                /* c.4G>A only abuts the agreeing c.5_6 islands, so it stays its own edit as in
                 * ordinary replay; both still change codon 2 and form one block. */
                ASSERT_EQ(2u, leaf.joint_edit_count); ASSERT_EQ(1u, leaf.block_count);
                ASSERT_EQ(4u, leaf.joint_edits[0].cds_start); ASSERT_EQ(1u, leaf.joint_edits[0].ref_len);
                ASSERT_EQ(0, memcmp("A", leaf.joint_edits[0].alt, 1u));
                ASSERT_EQ(5u, leaf.joint_edits[1].cds_start); ASSERT_EQ(2u, leaf.joint_edits[1].ref_len);
                ASSERT_EQ(0, memcmp("TA", leaf.joint_edits[1].alt, 2u));
            } else ASSERT(!leaf.joint_edits);
        } else {
            ASSERT_EQ(0, memcmp("ATGATTGAATAA", leaf.cds, 12u));
            ASSERT_EQ(0, memcmp("MIE*", leaf.protein, 4u));
            ASSERT_EQ(1u, leaf.joint_edit_count); ASSERT_EQ(2u, leaf.edit_count);
            ASSERT_EQ(4u, leaf.joint_edits[0].cds_start);
            ASSERT_EQ(0, memcmp(scenario == 6u ? "gc" : "GC", leaf.joint_edits[0].ref, 2u));
            ASSERT_EQ(0, memcmp("AT", leaf.joint_edits[0].alt, 2u));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
        ASSERT_EQ(0u, s->event_count); ASSERT_EQ(0u, s->projection_count);
    }
    PASS();
}

TEST haplotype_stream_overlap_post_stop_roles_follow_source_islands(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) for (unsigned mixed = 0u; mixed < 2u; mixed++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u); f.strands[0] = reverse ? -1 : 1;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        duckvep_haplotype_source_t sources[3] = {
            {1u, (const uint8_t *)(mixed ? (reverse ? "TCAGC" : "GCTGA") : (reverse ? "AGC" : "GCT")),
                (const uint8_t *)(mixed ? (reverse ? "GCTTA" : "TAAGC") : (reverse ? "TTA" : "TAA")),
                reverse ? (mixed ? 104u : 106u) : 103u, 0u, mixed ? 5u : 3u, mixed ? 5u : 3u, 0u, 0u, 0u, 1u},
            {2u, (const uint8_t *)(reverse ? "AGC" : "GCT"), (const uint8_t *)(reverse ? "TTA" : "TAA"),
                reverse ? 106u : 103u, 0u, 3u, 3u, 0u, 0u, 0u, 1u},
            {3u, (const uint8_t *)(reverse ? "T" : "A"), (const uint8_t *)(reverse ? "G" : "C"),
                reverse ? 104u : 107u, 0u, 1u, 1u, 0u, 0u, 0u, 1u}};
        unsigned order[3] = {0u, 1u, 2u};
        if (reverse) { order[0] = mixed ? 0u : 2u; order[1] = mixed ? 2u : 0u; order[2] = 1u; }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        for (unsigned i = 0u; i < 3u; i++) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, sources+order[i], &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(2u, leaf.joint_edit_count); ASSERT_EQ(3u+mixed, leaf.edit_count);
        ASSERT_MEM_EQ("ATGTAAGCATAA", leaf.cds, 12u); ASSERT_MEM_EQ("M*", leaf.protein, 2u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_NMD_ESCAPE, leaf.nmd);
        for (size_t i = 0u; i < leaf.contributor_count; i++) {
            const duckvep_haplotype_contributor_t *c = &leaf.contributors[i];
            ASSERT_EQ(c->source.event_id == 3u ? DUCKVEP_ROLE_POST_STOP : DUCKVEP_ROLE_APPLIED, c->role);
            ASSERT_EQ(c->source.event_id == 3u || (mixed && c->source.event_id == 1u) ? 1u : 0u, c->post_stop_edits);
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }
    PASS();
}

/* Intersecting literal indel islands on ATG GCT GAA TAA. Each source is written as an anchored
 * CDS-orientation REF/ALT at a 1-based CDS position; the reverse-strand scene carries its reverse
 * complement at the mirrored genomic position. Expected alleles are derived by hand from the
 * immutable reference: duplicate assertions are one edit, uniquely placed deletions remove the
 * union of their bases, and an unslidable insertion precedes the substitution at its site. */
struct literal_overlap_case {
    const char *ref[3], *alt[3];
    uint32_t cds_pos[3];
    unsigned sources;
    const char *cds, *protein;
    uint32_t joint_start, joint_ref_len;
    const char *joint_alt;
    unsigned post_stop_mask;
};

static char literal_overlap_complement(char base) {
    return base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
}

static void literal_overlap_source(const struct literal_overlap_case *c, unsigned i, int reverse,
    char ref[8], char alt[8], duckvep_haplotype_source_t *source) {
    size_t rl = strlen(c->ref[i]), al = strlen(c->alt[i]);
    for (size_t j = 0u; j < rl; j++)
        ref[j] = reverse ? literal_overlap_complement(c->ref[i][rl - 1u - j]) : c->ref[i][j];
    for (size_t j = 0u; j < al; j++)
        alt[j] = reverse ? literal_overlap_complement(c->alt[i][al - 1u - j]) : c->alt[i][j];
    *source = (duckvep_haplotype_source_t){i + 1u, (const uint8_t *)ref, (const uint8_t *)alt,
        reverse ? 112u - (c->cds_pos[i] + (uint32_t)rl - 1u) : 99u + c->cds_pos[i],
        0u, (uint32_t)rl, (uint32_t)al, 0u, 0u, 0u, 0u};
}

TEST haplotype_stream_literal_indel_overlaps_compose_or_refuse(void) {
    static const struct literal_overlap_case composed[] = {
        /* duplicate deletion of GCT (c.4_6), itself slidable to CTG */
        {{"GGCT", "GGCT"}, {"G", "G"}, {3u, 3u}, 2u, "ATGGAATAA", "ME*", 4u, 3u, "", 0u},
        /* the same deletion asserted three times */
        {{"GGCT", "GGCT", "GGCT"}, {"G", "G", "G"}, {3u, 3u, 3u}, 3u, "ATGGAATAA", "ME*", 4u, 3u, "", 0u},
        /* c.5_6 nested inside c.4_9, both unslidable */
        {{"GGCTGAA", "GCT"}, {"G", "G"}, {3u, 4u}, 2u, "ATGTAA", "M*", 4u, 6u, "", 0u},
        /* c.5_6 and c.6_7 overlap at c.6, both unslidable */
        {{"GCT", "CTG"}, {"G", "C"}, {4u, 5u}, 2u, "ATGGAATAA", "ME*", 5u, 3u, "", 0u},
        /* identical CCC insertions between c.6 and c.7 */
        {{"T", "T"}, {"TCCC", "TCCC"}, {6u, 6u}, 2u, "ATGGCTCCCGAATAA", "MAPE*", 7u, 0u, "CCC", 0u},
        /* identical insertions at the terminal interbase site */
        {{"A", "A"}, {"AC", "AC"}, {12u, 12u}, 2u, "ATGGCTGAATAAC", "MAE*", 13u, 0u, "C", 3u},
        /* unslidable CCC insertion before c.7G>A */
        {{"T", "G"}, {"TCCC", "A"}, {6u, 7u}, 2u, "ATGGCTCCCAAATAA", "MAPK*", 7u, 1u, "CCCA", 0u},
        /* c.7G>T stop, then a duplicated c.11_12 deletion after it */
        {{"G", "TAA", "TAA"}, {"T", "T", "T"}, {7u, 10u, 10u}, 3u, "ATGGCTTAAT", "MA*", 7u, 1u, "T", 6u},
        /* Opaque replacement payloads have no per-base ancestry; exact aliases apply once. */
        {{"GCT", "GCT"}, {"TTTCCC", "TTTCCC"}, {4u, 4u}, 2u, "ATGTTTCCCGAATAA", "MFPE*", 4u, 3u, "TTTCCC", 0u},
        {{"GCTGAA", "GCTGAA"}, {"CCC", "CCC"}, {4u, 4u}, 2u, "ATGCCCTAA", "MP*", 4u, 6u, "CCC", 0u},
        {{"GCT", "GCT"}, {"TAAGAA", "TAAGAA"}, {4u, 4u}, 2u, "ATGTAAGAAGAATAA", "M*", 4u, 3u, "TAAGAA", 0u},
        {{"GCT", "GCT", "A"}, {"TAAGAA", "TAAGAA", "AC"}, {4u, 4u, 12u}, 3u, "ATGTAAGAAGAATAAC", "M*", 4u, 3u, "TAAGAA", 4u},
        /* Records apply at their own coordinates, also where an equivalent placement exists in a
         * repeat: c.4_6del and c.5_6del remove c.4_6; CCC before c.7G>A; AAA inserted before the
         * deleted c.5_6. */
        {{"GGCT", "GCT"}, {"G", "G"}, {3u, 4u}, 2u, "ATGGAATAA", "ME*", 4u, 3u, "", 0u},
        {{"T", "G"}, {"TGGG", "A"}, {6u, 7u}, 2u, "ATGGCTGGGAAATAA", "MAGK*", 7u, 1u, "GGGA", 0u},
        {{"G", "GCT"}, {"GAAA", "G"}, {4u, 4u}, 2u, "ATGGAAAGAATAA", "MERI", 5u, 2u, "AAA", 0u},
    };
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned w = 0u; w < sizeof(composed) / sizeof(composed[0]); w++)
    for (unsigned swap = 0u; swap < 2u; swap++) {
        const struct literal_overlap_case *c = &composed[w];
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u);
        f.strands[0] = reverse ? -1 : 1;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        char refs[3][8], alts[3][8];
        duckvep_haplotype_source_t sources[3];
        for (unsigned i = 0u; i < c->sources; i++) literal_overlap_source(c, i, (int)reverse, refs[i], alts[i], &sources[i]);
        /* Genomic order is required; sources at one position take event ids in either order. */
        unsigned order[3] = {0u, 1u, 2u};
        for (unsigned i = 1u; i < c->sources; i++)
            for (unsigned j = i; j > 0u; j--) {
                const duckvep_haplotype_source_t *a = &sources[order[j - 1u]], *b = &sources[order[j]];
                if (a->pos1 > b->pos1 || (swap && a->pos1 == b->pos1)) {
                    unsigned t = order[j]; order[j] = order[j - 1u]; order[j - 1u] = t;
                }
            }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned i = 0u; i < c->sources; i++) {
            duckvep_haplotype_source_t source = sources[order[i]];
            source.event_id = i + 1u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        size_t length = strlen(c->cds), residues = strlen(c->protein);
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(length, leaf.cds_length); ASSERT_MEM_EQ(c->cds, leaf.cds, length);
        ASSERT_EQ(residues, leaf.protein_length); ASSERT_MEM_EQ(c->protein, leaf.protein, residues);
        ASSERT_EQ((int64_t)length - 12, leaf.nominal_length_diff);
        ASSERT_EQ(length != 12u, (leaf.flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u);
        ASSERT_EQ(c->sources, leaf.contributor_count); ASSERT_EQ(c->sources, leaf.edit_count);
        const duckvep_haplotype_edit_t *last = &leaf.joint_edits[leaf.joint_edit_count - 1u];
        unsigned joint_count = c->post_stop_mask && c->sources == 3u ? 2u : 1u;
        ASSERT_EQ(joint_count, leaf.joint_edit_count);
        ASSERT_EQ(joint_count, leaf.block_count);
        const duckvep_haplotype_edit_t *first = &leaf.joint_edits[0];
        ASSERT_EQ(c->joint_start, first->cds_start); ASSERT_EQ(c->joint_ref_len, first->ref_len);
        ASSERT_EQ(strlen(c->joint_alt), first->alt_len);
        ASSERT_MEM_EQ(c->joint_alt, first->alt, first->alt_len);
        if (c->post_stop_mask == 6u) {
            ASSERT_EQ(11u, last->cds_start); ASSERT_EQ(2u, last->ref_len); ASSERT_EQ(0u, last->alt_len);
        } else if (c->post_stop_mask) {
            ASSERT_EQ(13u, last->cds_start); ASSERT_EQ(0u, last->ref_len); ASSERT_EQ(1u, last->alt_len);
            ASSERT_MEM_EQ("C", last->alt, 1u);
        }
        for (size_t i = 0u; i < leaf.contributor_count; i++) {
            const duckvep_haplotype_contributor_t *k = &leaf.contributors[i];
            unsigned post = (c->post_stop_mask >> order[k->source.event_id - 1u]) & 1u;
            ASSERT_EQ(1u, k->edit_count); ASSERT_EQ(post, k->post_stop_edits);
            ASSERT_EQ(post ? DUCKVEP_ROLE_POST_STOP : DUCKVEP_ROLE_APPLIED, k->role);
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(0u, f.stream.event_count);
    }

    /* Contradictory records and reference mismatches publish no sequence. */
    static const struct {
        struct literal_overlap_case c;
        duckvep_prediction_status_t status;
    } refused[] = {
        /* terminal insertion aliases must match their original reference anchor */
        {{{"T", "T"}, {"TC", "TC"}, {12u, 12u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT},
        /* different insertions in one interbase site have no order */
        {{{"T", "T"}, {"TCCC", "TAAA"}, {6u, 6u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_EDIT_CONFLICT},
        /* the c.5_6 deletion removes the base c.5C>A changes */
        {{{"GCT", "C"}, {"G", "A"}, {4u, 5u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_EDIT_CONFLICT},
        /* so does c.4_6del, though it could be written c.5_7del */
        {{{"GGCT", "C"}, {"G", "A"}, {3u, 5u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_EDIT_CONFLICT},
        /* c.8del removes the base c.8A>T changes; c.9del would be an equivalent deletion */
        {{{"GA", "A"}, {"G", "T"}, {7u, 8u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_EDIT_CONFLICT},
        /* a distinct mutation inside an opaque replacement claims one of its bases */
        {{{"CT", "C"}, {"AAA", "G"}, {5u, 5u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_EDIT_CONFLICT},
        /* normalized aliases cannot conceal an incorrect unchanged padding base */
        {{{"GCT", "AGCT"}, {"TTTCCC", "ATTTCCC"}, {4u, 3u}, 2u, NULL, NULL, 0u, 0u, NULL, 0u}, DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT},
    };
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned w = 0u; w < sizeof(refused) / sizeof(refused[0]); w++) {
        const struct literal_overlap_case *c = &refused[w].c;
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u);
        f.strands[0] = reverse ? -1 : 1;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        char refs[2][8], alts[2][8];
        duckvep_haplotype_source_t sources[2];
        for (unsigned i = 0u; i < 2u; i++) literal_overlap_source(c, i, (int)reverse, refs[i], alts[i], &sources[i]);
        unsigned lo = sources[0].pos1 > sources[1].pos1;
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned i = 0u; i < 2u; i++) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &sources[i ? 1u - lo : lo], &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(refused[w].status, leaf.path_status);
        if (refused[w].status == DUCKVEP_PREDICTION_EDIT_CONFLICT)
            ASSERT_EQ(DUCKVEP_REASON_CONTRADICTORY_EDITS, leaf.path_reason);
        if (refused[w].status == DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT)
            ASSERT(leaf.path_reason == DUCKVEP_REASON_REFERENCE_MISMATCH ||
                (leaf.path_reason == DUCKVEP_REASON_PROJECTION &&
                 leaf.prediction_projection == DUCKVEP_CDS_EDIT_REF_MISMATCH));
        ASSERT(!leaf.cds && !leaf.protein && !leaf.joint_edits);
        ASSERT_EQ(0u, leaf.cds_length); ASSERT_EQ(0u, leaf.block_count);
        for (size_t i = 0u; i < leaf.contributor_count; i++) {
            const duckvep_haplotype_contributor_t *k = &leaf.contributors[i];
            ASSERT_EQ(k->edit_count || k->source_replaced ? DUCKVEP_ROLE_UNAPPLIED : DUCKVEP_ROLE_OMITTED, k->role);
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }

    /* The grown CDS is bounded by the charged CDS buffer; overflow fails the stream atomically. */
    {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u);
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        f.buffers.cds_capacity = 14u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned i = 0u; i < 2u; i++) {
            duckvep_haplotype_source_t source = {i + 1u, (const uint8_t *)"T", (const uint8_t *)"TCCC",
                105u, 0u, 1u, 4u, 0u, 0u, 0u, 0u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        memset(&leaf, 0xa5, sizeof(leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, duckvep_haplotype_stream_next(&f.stream, &leaf));
        duckvep_haplotype_leaf_t zero = {0};
        ASSERT_MEM_EQ(&zero, &leaf, sizeof(leaf));
    }
    PASS();
}

TEST haplotype_stream_opaque_endpoint_order(void) {
    const struct literal_overlap_case c = {{"G", "GCT"}, {"GCCC", "TTTCCC"},
        {3u, 4u}, 2u, "ATGCCCTTTCCCGAATAA", NULL, 0u, 0u, NULL, 0u};
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned swap = 0u; swap < 2u; swap++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u);
        f.strands[0] = reverse ? -1 : 1;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream,
            &f.model, &f.exons, &f.sequences, &f.buffers));
        char refs[2][8], alts[2][8]; duckvep_haplotype_source_t sources[2];
        for (unsigned i = 0u; i < 2u; i++)
            literal_overlap_source(&c, i, (int)reverse, refs[i], alts[i], &sources[i]);
        if (swap) { sources[0].event_id = 2u; sources[1].event_id = 1u; }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned i = 0u; i < 2u; i++) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream,
                &sources[reverse ? 1u - i : i], &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(18u, leaf.cds_length); ASSERT_MEM_EQ(c.cds, leaf.cds, 18u);
        ASSERT_EQ(6, leaf.nominal_length_diff);
        ASSERT_EQ(1u, leaf.joint_edit_count);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }
    PASS();
}

TEST haplotype_stream_boundary_insertions_preserve_source_exons(void) {
    const struct literal_overlap_case c = {{"T", "C"}, {"TCCC", "CCCC"},
        {6u, 7u}, 2u, "ATGGCTCCCCCCCAATAA", NULL, 0u, 0u, NULL, 0u};
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned swap = 0u; swap < 2u; swap++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTCAATAA", 12u);
        f.strands[0] = reverse ? -1 : 1;
        uint32_t exon_start[2] = {reverse ? 106u : 100u, reverse ? 100u : 106u};
        uint32_t exon_end[2] = {reverse ? 111u : 105u, reverse ? 105u : 111u};
        uint32_t cdna_start[2] = {1u, 7u}, cdna_end[2] = {6u, 12u};
        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
        f.exons.start1 = exon_start; f.exons.end1 = exon_end;
        f.exons.cdna_start1 = cdna_start; f.exons.cdna_end1 = cdna_end;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream,
            &f.model, &f.exons, &f.sequences, &f.buffers));
        char refs[2][8], alts[2][8]; duckvep_haplotype_source_t sources[2];
        for (unsigned i = 0u; i < 2u; i++)
            literal_overlap_source(&c, i, (int)reverse, refs[i], alts[i], &sources[i]);
        unsigned order[2] = {reverse ? 1u : 0u, reverse ? 0u : 1u};
        if (swap) { sources[0].event_id = 2u; sources[1].event_id = 1u; }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned i = 0u; i < 2u; i++) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream,
                &sources[order[i]], &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        /* The insertions sit in different genomic gaps on either side of the intron: both apply, in
         * transcript order, on either strand. */
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(18u, leaf.cds_length); ASSERT_MEM_EQ(c.cds, leaf.cds, 18u);
        ASSERT_EQ(6, leaf.nominal_length_diff);
        ASSERT_EQ(2u, leaf.contributor_count);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }
    PASS();
}

TEST haplotype_stream_cross_exon_opaque_alias_keeps_atomic_geometry(void) {
    const struct literal_overlap_case c = {{"GCT", "GCT", "G"}, {"TTTCCC", "TTTCCC", "A"},
        {4u, 4u, 7u}, 3u, "ATGTTTCCCAAATAA", "MFPK*", 4u, 4u, "TTTCCCA", 0u};
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned swap = 0u; swap < 2u; swap++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u);
        f.strands[0] = reverse ? -1 : 1;
        uint32_t exon_start[2] = {reverse ? 106u : 100u, reverse ? 100u : 106u};
        uint32_t exon_end[2] = {reverse ? 111u : 105u, reverse ? 105u : 111u};
        uint32_t cdna_start[2] = {1u, 7u}, cdna_end[2] = {6u, 12u};
        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
        f.exons.start1 = exon_start; f.exons.end1 = exon_end;
        f.exons.cdna_start1 = cdna_start; f.exons.cdna_end1 = cdna_end;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream,
            &f.model, &f.exons, &f.sequences, &f.buffers));
        char refs[3][8], alts[3][8]; duckvep_haplotype_source_t sources[3];
        for (unsigned i = 0u; i < 3u; i++)
            literal_overlap_source(&c, i, (int)reverse, refs[i], alts[i], &sources[i]);
        unsigned order[3] = {0u, 1u, 2u};
        if (swap) { sources[0].event_id = 2u; sources[1].event_id = 1u;
            order[0] = 1u; order[1] = 0u; }
        if (reverse) { order[0] = 2u; order[1] = swap ? 1u : 0u; order[2] = swap ? 0u : 1u; }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned i = 0u; i < 3u; i++) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream,
                &sources[order[i]], &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(15u, leaf.cds_length); ASSERT_MEM_EQ(c.cds, leaf.cds, 15u);
        ASSERT_EQ(3u, leaf.edit_count);
        ASSERT(leaf.joint_edit_count >= 1u && leaf.joint_edit_count <= 2u);
        ASSERT_EQ(4u, leaf.joint_edits[0].cds_start);
        ASSERT_EQ(3, leaf.nominal_length_diff);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }
    PASS();
}

TEST haplotype_stream_opaque_alias_capacity_and_reuse(void) {
    const struct literal_overlap_case c = {{"GCT", "GCT", "A"}, {"TTTCCC", "TTTCCC", "AC"},
        {4u, 4u, 12u}, 3u, "ATGTTTCCCGAATAAC", "MFPE*", 4u, 3u, "TTTCCC", 0u};
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned terminal = 0u; terminal < 2u; terminal++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGGCTGAATAA", 12u); f.strands[0] = reverse ? -1 : 1;
        duckvep_haplotype_edit_t joint[8]; uint8_t alt[12], assigned[12];
        f.buffers.joint_edits = joint; f.buffers.joint_alt = alt;
        f.buffers.joint_assigned = assigned; f.buffers.joint_capacity = 12u;
        char refs[3][8] = {{0}}, alts[3][8] = {{0}};
        char refs_before[3][8], alts_before[3][8];
        duckvep_haplotype_source_t sources[3];
        for (unsigned i = 0u; i < 3u; i++) literal_overlap_source(&c, i, (int)reverse, refs[i], alts[i], &sources[i]);
        memcpy(refs_before, refs, sizeof(refs)); memcpy(alts_before, alts, sizeof(alts));
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        size_t length = 15u + terminal;
        for (unsigned attempt = 0u; attempt < 3u; attempt++) {
            f.buffers.cds_capacity = attempt ? length : length - 1u;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            for (unsigned i = 0u; i < 2u + terminal; i++) {
                unsigned source = reverse && terminal ? (i + 2u) % 3u : i;
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &sources[source], &tx, 1u));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
            }
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            memset(&leaf, 0xa5, sizeof(leaf));
            if (!attempt) {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, duckvep_haplotype_stream_next(&f.stream, &leaf));
                duckvep_haplotype_leaf_t zero = {0};
                ASSERT_MEM_EQ(&zero, &leaf, sizeof(leaf));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL, duckvep_haplotype_stream_next(&f.stream, &leaf));
            } else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
                ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
                ASSERT_EQ(length, leaf.cds_length); ASSERT_MEM_EQ(c.cds, leaf.cds, length);
                ASSERT_EQ(1u + terminal, leaf.joint_edit_count);
                ASSERT_EQ(f.cds + 3u, leaf.joint_edits[0].alt);
                ASSERT_MEM_EQ("TTTCCC", leaf.joint_edits[0].alt, 6u);
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(&f.stream));
                ASSERT_EQ(0u, f.stream.event_count); ASSERT_EQ(0u, f.stream.projection_count);
            }
            ASSERT_MEM_EQ("ATGGCTGAATAA", f.reference, 12u);
            ASSERT_MEM_EQ(refs_before, refs, sizeof(refs));
            ASSERT_MEM_EQ(alts_before, alts, sizeof(alts));
        }
    }
    PASS();
}

TEST haplotype_stream_mnv_context_does_not_conflict_with_another_edit(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        if (reverse) {
            memset(f.reference, 'T', sizeof(f.reference));
            f.strands[0] = -1;
        }
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_carrier_key_t key = {0u, 10, 1u, 2u, 1u, 0u};
        uint32_t tx = 0u;
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)"AAA", (const uint8_t *)"CAC",
            100u, 0u, 3u, 3u, 0u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        source.event_id = 2u; source.pos1 = 101u; source.ref_len = source.alt_len = 1u;
        source.ref = (const uint8_t *)"A"; source.alt = (const uint8_t *)"G";
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
        ASSERT_EQ(2u, leaf.contributor_count); ASSERT_EQ(3u, leaf.edit_count);
        ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
        ASSERT_EQ(0, memcmp("CAC", leaf.contributors[0].source.alt, 3u));
        ASSERT_EQ(2u, leaf.contributors[1].source.event_id);
        ASSERT_EQ(12u, leaf.cds_length); ASSERT_EQ(4u, leaf.protein_length);
        ASSERT_EQ(0, memcmp(reverse ? "TTTTTTTTTGCG" : "CGCAAAAAAAAA", leaf.cds, 12u));
        ASSERT_EQ(0, memcmp(reverse ? "FFFA" : "RKKK", leaf.protein, 4u));
        ASSERT_EQ(2u, s->projected_events); ASSERT_EQ(1u, s->completed_leaves);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    }
    PASS();
}

TEST haplotype_stream_predicts_complete_utr_substitution_from_joint_cdna(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    f.ends[0] = 119u; f.cdna_ends[0] = 20u;
    f.model.cds_start1 = (uint32_t[]){106u}; f.model.cds_end1 = (uint32_t[]){114u};
    f.lengths[0] = 9u;
    memcpy(f.reference, "ATGAAATAA", 9u);
    memcpy(f.flank, "AAAAAACCCCC", 11u);
    f.sequences.flank_bytes = f.flank; f.sequences.flank_bytes_len = 11u;
    f.sequences.pre_cds_offset = f.pre_flank_offsets; f.sequences.pre_cds_length = f.pre_flank_lengths;
    f.sequences.post_cds_offset = f.post_flank_offsets; f.sequences.post_cds_length = f.post_flank_lengths;
    f.pre_flank_lengths[0] = 6u; f.post_flank_offsets[0] = 6u; f.post_flank_lengths[0] = 5u;
    f.sequences.flanks_complete = 1u;
    f.buffers.noncoding = f.noncoding; f.buffers.noncoding_capacity = sizeof(f.noncoding);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
    uint32_t tx = 0u;
    duckvep_haplotype_source_t source = {1u, (const uint8_t *)"A", (const uint8_t *)"C",
        102u, 0u, 1u, 1u, 0u, 0u, 0u, 1u};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
    duckvep_carrier_key_t key = carrier_test_key(1u);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
    duckvep_haplotype_leaf_t leaf;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
    ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR), leaf.haplotype_so_mask);
    ASSERT_EQ((int)DUCKVEP_CDS_EDIT_OUT_OF_CDS, (int)leaf.contributors[0].projection_status);
    ASSERT_EQ(1u, leaf.contributors[0].source.event_id);
    ASSERT_EQ(0, memcmp(f.cds, "ATGAAATAA", 9u));
    uint32_t spliced_end[] = {125u}, spliced_cds_start[] = {112u}, spliced_cds_end[] = {120u};
    for (unsigned control = 0u; control < 3u; control++) {
        f.pre_flank_lengths[0] = control == 0u ? 7u : 6u;
        f.post_flank_offsets[0] = control == 0u ? 7u : 6u;
        f.post_flank_lengths[0] = control < 2u ? 4u : 5u;
        if (control == 2u) {
            f.model.end1 = spliced_end;
            f.model.cds_start1 = spliced_cds_start; f.model.cds_end1 = spliced_cds_end;
            f.exon_count[0] = 2u; f.exons.exon_count = 2u;
            f.starts[1] = 110u; f.ends[0] = 103u; f.ends[1] = 125u;
            f.cdna_starts[1] = 5u; f.cdna_ends[0] = 4u; f.cdna_ends[1] = 20u;
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ((int)DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, (int)leaf.path_status);
        ASSERT_EQ(0u, leaf.haplotype_so_mask);
        ASSERT_EQ(0, memcmp(f.cds, "ATGAAATAA", 9u));
    }
    PASS();
}

TEST haplotype_stream_conditional_rna_on_both_strands(void) {
    static const char raw[] = "AAAAAAATGAAAGGGTAA CCC";
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned scenario = 0u; scenario < 21u; scenario++) {
        struct haplotype_stream_scene f;
        uint8_t prediction[128] = {0}, alt[] = {'G', 'A'};
        uint8_t presence[] = {2u}, boundary[] = {0u}, kind[] = {1u, 1u};
        uint32_t offsets[] = {0u, 1u}, starts[] = {10u, 7u}, ends[] = {10u, 7u};
        uint32_t alt_lengths[] = {1u, 1u}, cds_starts[] = {7u}, cds_ends[] = {18u};
        uint64_t alt_offsets[] = {0u, 1u};
        uint8_t table[] = {1u}, compatibility_table[] = {0u};
        uint64_t flags[] = {DUCKVEP_TX_RNA_EDIT};
        uint8_t ref[] = {reverse ? 'C' : 'G'}, genomic_alt[] = {reverse ? 'T' : 'A', 'A'};
        unsigned position = 13u;
        if (scenario == 1u || scenario == 3u) position = 10u;
        if (scenario == 2u) position = 11u;
        if (scenario == 13u) position = 7u;
        if (scenario >= 19u) {
            position = 16u;
            ref[0] = reverse ? 'A' : 'T';
            genomic_alt[0] = reverse ? 'G' : 'C';
        }
        if ((scenario >= 1u && scenario <= 3u) || scenario == 13u) {
            ref[0] = reverse ? 'T' : 'A';
            genomic_alt[0] = scenario == 3u || scenario == 13u
                ? (reverse ? 'C' : 'G') : (reverse ? 'G' : 'C');
        }
        if (scenario == 4u) boundary[0] = 1u;
        if (scenario == 5u) flags[0] |= DUCKVEP_TX_CDS_START_NF;
        if (scenario == 6u) presence[0] = 1u;
        if (scenario == 7u) cds_ends[0] = 17u;
        if (scenario == 8u) alt[0] = 'A';
        if (scenario == 9u) alt_lengths[0] = 0u;
        if (scenario == 10u) ref[0] = 'A';
        if (scenario == 12u) table[0] = 2u;
        if (scenario == 14u) { starts[0] = ends[0] = 7u; alt[0] = 'G'; }
        if (scenario == 16u) flags[0] |= DUCKVEP_TX_INITIAL_METHIONINE;
        if (scenario == 17u) { starts[0] = ends[0] = 7u; alt[0] = 'A'; }
        if (scenario == 18u) offsets[1] = 2u;
        haplotype_stream_scene_prepare(&f, 1u);
        f.strands[0] = reverse ? -1 : 1;
        f.model.flags = flags;
        f.ends[0] = 120u; f.cdna_ends[0] = 21u;
        f.model.cds_start1 = (uint32_t[]){reverse ? 103u : 106u};
        f.model.cds_end1 = (uint32_t[]){reverse ? 114u : 117u};
        f.lengths[0] = 0u; f.sequences.cds_bytes = NULL;
        f.sequences.cds_bytes_len = 0u;
        f.sequences.codon_table = compatibility_table;
        memcpy(f.cdna, raw, 18u); memcpy(f.cdna + 18u, "CCC", 3u);
        if (scenario == 15u || scenario == 17u) f.cdna[6] = 'G';
        f.cdna_lengths[0] = scenario == 20u ? 20u : 21u;
        if (scenario == 20u) f.cdna_ends[0] = 20u;
        f.sequences.cdna_bytes = f.cdna; f.sequences.cdna_bytes_len = f.cdna_lengths[0];
        f.sequences.cdna_offset = f.cdna_offsets; f.sequences.cdna_length = f.cdna_lengths;
        f.sequences.cdna_provided = 1u;
        f.buffers.noncoding = f.noncoding; f.buffers.noncoding_capacity = sizeof(f.noncoding);
        f.buffers.prediction_protein = prediction; f.buffers.prediction_protein_capacity = sizeof(prediction);
        duckvep_haplotype_source_t source = {101u, ref, genomic_alt,
            reverse ? 121u - position : 99u + position, 0u, 1u,
            scenario == 11u ? 2u : 1u, 0u, 0u, 0u, 1u};
        uint32_t tx = 0u;
        duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned replay = 0u; replay < (scenario == 0u ? 2u : 1u); replay++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
            duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        f.stream.rna = (duckvep_rna_programme_view_t){.presence = presence, .boundary_overrides = boundary,
            .offsets = offsets, .starts = starts, .ends = ends, .kinds = kind,
            .alt_offsets = alt_offsets, .alt_lengths = alt_lengths, .alt_bytes = alt,
            .edit_count = scenario == 18u ? 2u : 1u,
            .alt_bytes_length = scenario == 18u ? 2u : 1u,
            .cds_starts = cds_starts, .cds_ends = cds_ends, .codon_tables = table};
        duckvep_haplotype_stream_status_t begun = haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u);
        if (scenario == 10u && begun == DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG) continue;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, begun);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        if (scenario && scenario != 12u && scenario != 17u && scenario != 18u && scenario < 19u) {
            ASSERT_EQ(1, leaf.path_status != DUCKVEP_PREDICTION_CONDITIONAL_RNA);
            ASSERT_EQ(0u, leaf.haplotype_so_mask);
            continue;
        }
        ASSERT_EQ(DUCKVEP_PREDICTION_CONDITIONAL_RNA, leaf.path_status);
        ASSERT_EQ(DUCKVEP_REASON_ASSUMED_RNA_PROGRAMME, leaf.path_reason);
        ASSERT_EQ(DUCKVEP_PREDICTION_CONDITIONAL_RNA, leaf.prediction_status);
        ASSERT_EQ(DUCKVEP_REASON_ASSUMED_RNA_PROGRAMME, leaf.prediction_reason);
        ASSERT_EQ(4u, leaf.prediction_reference_protein_length);
        ASSERT_EQ(scenario == 19u ? 5u : scenario == 12u ? 3u : 4u,
            leaf.prediction_protein_length);
        ASSERT_MEM_EQ((scenario == 17u ? "MKG*" : "MEG*"), leaf.prediction_reference_protein, 4u);
        ASSERT_MEM_EQ((scenario == 19u ? "MEGQP" : scenario == 20u ? "MEGQ" :
            scenario == 12u ? "ME*" : scenario == 17u ? "MKR*" : "MER*"),
            leaf.prediction_protein, scenario == 19u ? 5u : scenario == 12u ? 3u : 4u);
        ASSERT_EQ(0u, leaf.haplotype_so_mask);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_NMD_UNKNOWN, leaf.nmd);
        ASSERT_MEM_EQ((scenario == 17u ? "AAAAAAGTGAAAGGGTAA" : "AAAAAAATGAAAGGGTAA"), f.cdna, 18u);
        }
    }
    return 0;
}

TEST haplotype_stream_replays_noncoding_indels_on_both_strands(void) {
    static const char reference[] = "AAAAAAATGAAATAACCCCC";
    for (unsigned lower = 0u; lower < 2u; lower++)
    for (unsigned feature = 0u; feature < 3u; feature++)
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned kind = 0u; kind < 5u; kind++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.strands[0] = reverse ? -1 : 1;
        f.model.end1 = (uint32_t[]){119u};
        f.ends[0] = 119u; f.cdna_ends[0] = 20u;
        f.model.cds_start1 = (uint32_t[]){feature ? (reverse ? 105u : 106u) : 0u};
        f.model.cds_end1 = (uint32_t[]){feature ? (reverse ? 113u : 114u) : 0u};
        f.lengths[0] = 9u; memcpy(f.reference, "ATGAAATAA", 9u);
        memcpy(f.flank, "AAAAAACCCCC", 11u);
        f.sequences.flank_bytes = f.flank; f.sequences.flank_bytes_len = 11u;
        f.sequences.pre_cds_offset = f.pre_flank_offsets; f.sequences.pre_cds_length = f.pre_flank_lengths;
        f.sequences.post_cds_offset = f.post_flank_offsets; f.sequences.post_cds_length = f.post_flank_lengths;
        f.pre_flank_lengths[0] = 6u; f.post_flank_offsets[0] = 6u; f.post_flank_lengths[0] = 5u;
        f.sequences.flanks_complete = 1u;
        memcpy(f.cdna, reference, 20u); f.cdna_lengths[0] = 20u;
        if (lower) {
            memcpy(f.reference, "atgaaataa", 9u);
            memcpy(f.flank, "aaaaaaccccc", 11u);
            memcpy(f.cdna, "aaaaaaatgaaataaccccc", 20u);
        }
        f.sequences.cdna_bytes = f.cdna; f.sequences.cdna_bytes_len = 20u;
        f.sequences.cdna_offset = f.cdna_offsets; f.sequences.cdna_length = f.cdna_lengths;
        f.sequences.cdna_provided = 1u;
        f.buffers.noncoding = f.noncoding; f.buffers.noncoding_capacity = sizeof(f.noncoding);
        unsigned position = feature == 2u ? 16u : 1u;
        unsigned separation = feature == 2u ? 2u : 3u; /* zero-based source REF position */
        char base = reference[position];
        uint8_t refs[2][3] = {{(uint8_t)base, (uint8_t)base, 0}, {(uint8_t)base, (uint8_t)base, 0}};
        uint8_t alts[2][4] = {{(uint8_t)base, 'G', 'T', 0}, {(uint8_t)base, 0, 0, 0}};
        size_t nr = kind == 1u || kind == 2u ? 2u : 1u;
        size_t na = kind == 1u ? 1u : kind == 2u ? 3u : 2u;
        if (kind >= 3u) alts[0][1] = kind == 4u ? (uint8_t)base : 'G';
        size_t count = kind >= 3u ? 2u : 1u;
        duckvep_haplotype_source_t sources[2] = {{0}};
        char expected[32] = {0};
        memcpy(expected, reference, position);
        memcpy(expected + position, alts[0], na);
        memcpy(expected + position + na, reference + position + nr, 20u - position - nr);
        size_t expected_length = 20u + na - nr;
        if (count == 2u) {
            size_t deletion = position + separation + 1u; /* VCF anchor then one deleted base */
            memmove(expected + deletion, expected + deletion + 1u, expected_length - deletion - 1u);
            expected_length--;
        }
        for (size_t i = 0u; i < count; i++) {
            size_t ref_len = i ? 2u : nr, alt_len = i ? 1u : na;
            unsigned p = position + (i ? separation : 0u);
            if (reverse) {
                for (size_t j = 0u; j < ref_len / 2u; j++) { uint8_t c = refs[i][j]; refs[i][j] = refs[i][ref_len-j-1u]; refs[i][ref_len-j-1u] = c; }
                for (size_t j = 0u; j < alt_len / 2u; j++) { uint8_t c = alts[i][j]; alts[i][j] = alts[i][alt_len-j-1u]; alts[i][alt_len-j-1u] = c; }
                for (size_t j = 0u; j < ref_len; j++) refs[i][j] = refs[i][j] == 'A' ? 'T' : refs[i][j] == 'C' ? 'G' : refs[i][j] == 'G' ? 'C' : 'A';
                for (size_t j = 0u; j < alt_len; j++) alts[i][j] = alts[i][j] == 'A' ? 'T' : alts[i][j] == 'C' ? 'G' : alts[i][j] == 'G' ? 'C' : 'A';
            }
            sources[i] = (duckvep_haplotype_source_t){101u+i, refs[i], alts[i],
                reverse ? 120u - p - (uint32_t)ref_len : 100u + p, 0u,
                (uint16_t)ref_len, (uint16_t)alt_len, 0u, 0u, 0u, 1u};
        }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        for (size_t k = 0u; k < count; k++) {
            size_t i = reverse ? count-k-1u : k;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, sources+i, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        if (leaf.path_status != DUCKVEP_PREDICTION_PREDICTED)
            fprintf(stderr, "RNA witness feature=%u reverse=%u kind=%u status=%d reason=%d edits=%zu p0=%d p1=%d\n", feature, reverse, kind, leaf.path_status, leaf.path_reason, leaf.edit_count, leaf.contributors[0].projection_status, count>1u ? leaf.contributors[1].projection_status : -1);
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_MEM_EQ(expected, f.noncoding + 20u, expected_length);
        ASSERT_MEM_EQ((lower ? "atgaaataa" : "ATGAAATAA"), f.reference, 9u);
        ASSERT_MEM_EQ((lower ? "aaaaaaccccc" : "AAAAAACCCCC"), f.flank, 11u);
        ASSERT_MEM_EQ((lower ? "aaaaaaatgaaataaccccc" : reference), f.cdna, 20u);
        if (kind == 4u) ASSERT_MEM_EQ(reference, f.noncoding + 20u, 20u);
        ASSERT_EQ(kind == 4u ? 0u : feature == 0u ? DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON) :
            feature == 1u ? DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) : DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR), leaf.haplotype_so_mask);
        ASSERT_EQ(0u, leaf.edit_count); ASSERT_EQ(count, leaf.contributor_count);
        if (feature) { ASSERT_MEM_EQ("ATGAAATAA", leaf.cds, 9u); ASSERT_MEM_EQ("MK*", leaf.protein, 3u); }
        else { ASSERT_EQ(NULL, leaf.cds); ASSERT_EQ(NULL, leaf.protein); }
    }
    PASS();
}

TEST haplotype_stream_projects_utr_indel_ejc50_and_capacity(void) {
    for (unsigned deletion = 0u; deletion < 2u; deletion++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.model.end1 = (uint32_t[]){239u};
        f.model.cds_start1 = (uint32_t[]){106u}; f.model.cds_end1 = (uint32_t[]){114u};
        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
        f.ends[0] = 179u; f.starts[1] = 200u; f.ends[1] = 239u;
        f.cdna_ends[0] = 80u; f.cdna_starts[1] = 81u; f.cdna_ends[1] = 120u;
        f.lengths[0] = 9u; memcpy(f.reference, "ATGAAATAA", 9u);
        memcpy(f.flank, "AAAAAA", 6u); memset(f.flank + 6u, 'C', 105u);
        f.sequences.flank_bytes = f.flank; f.sequences.flank_bytes_len = 111u;
        f.sequences.pre_cds_offset = f.pre_flank_offsets; f.sequences.pre_cds_length = f.pre_flank_lengths;
        f.sequences.post_cds_offset = f.post_flank_offsets; f.sequences.post_cds_length = f.post_flank_lengths;
        f.pre_flank_lengths[0] = 6u; f.post_flank_offsets[0] = 6u; f.post_flank_lengths[0] = 105u;
        f.sequences.flanks_complete = 1u;
        f.buffers.noncoding = f.noncoding;
        uint8_t ref[21]; memset(ref, 'C', sizeof ref);
        const uint8_t *alt = (const uint8_t *)(deletion ? "C" : "CAAAAAAAA");
        duckvep_haplotype_source_t source = {1u, ref, alt, 120u, 0u,
            deletion ? 21u : 1u, deletion ? 1u : 9u, 0u, 0u, 0u, 1u};
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        for (unsigned capacity = 0u; capacity < 3u; capacity++) {
            f.buffers.noncoding_capacity = capacity == 1u ?
                2u * (deletion ? 119u : 127u) : sizeof(f.noncoding);
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
            if (capacity == 1u) {
                ASSERT_EQ(DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, leaf.path_status);
                ASSERT_EQ(0u, leaf.haplotype_so_mask);
            } else {
                ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
                ASSERT_EQ(deletion ? DUCKVEP_HAPLOTYPE_NMD_ESCAPE : DUCKVEP_HAPLOTYPE_NMD_TRIGGER, leaf.nmd);
                ASSERT_EQ(15u, leaf.nmd_stop_position1);
                ASSERT_EQ(deletion ? 60u : 88u, leaf.nmd_junction_position1);
                ASSERT_MEM_EQ("ATGAAATAA", leaf.cds, 9u);
                ASSERT_MEM_EQ("MK*", leaf.protein, 3u);
            }
        }
    }
    PASS();
}

TEST haplotype_stream_preserves_restored_rna_and_tests_exon_geometry(void) {
    for (unsigned lower = 0u; lower < 2u; lower++)
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
    for (unsigned kind = 0u; kind < 3u; kind++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.strands[0] = reverse ? -1 : 1;
        f.model.start1 = (uint32_t[]){100u};
        f.model.end1 = (uint32_t[]){reverse ? 264u : 229u};
        f.model.cds_start1 = (uint32_t[]){reverse ? 250u : 106u};
        f.model.cds_end1 = (uint32_t[]){reverse ? 258u : 114u};
        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
        f.starts[0] = reverse ? 200u : 100u; f.starts[1] = reverse ? 100u : 200u;
        f.ends[0] = f.starts[0] + 64u; f.ends[1] = f.starts[1] + 29u;
        f.cdna_starts[1] = 66u; f.cdna_ends[0] = 65u; f.cdna_ends[1] = 95u;
        f.lengths[0] = 9u;
        memcpy(f.reference, lower ? "atgaaataa" : "ATGAAATAA", 9u);
        memcpy(f.flank, lower ? "aaaaaa" : "AAAAAA", 6u);
        memset(f.flank + 6u, lower ? 'c' : 'C', 80u);
        f.sequences.flank_bytes = f.flank; f.sequences.flank_bytes_len = 86u;
        f.sequences.pre_cds_offset = f.pre_flank_offsets; f.sequences.pre_cds_length = f.pre_flank_lengths;
        f.sequences.post_cds_offset = f.post_flank_offsets; f.sequences.post_cds_length = f.post_flank_lengths;
        f.pre_flank_lengths[0] = 6u; f.post_flank_offsets[0] = 6u; f.post_flank_lengths[0] = 80u;
        f.sequences.flanks_complete = 1u;
        f.buffers.noncoding = f.noncoding; f.buffers.noncoding_capacity = sizeof(f.noncoding);
        duckvep_haplotype_source_t sources[2] = {{0}};
        for (unsigned i = 0u; i < 2u; i++) {
            unsigned p = i ? (kind ? 75u : 35u) : 25u;
            unsigned exon = p <= 65u ? 0u : 1u;
            unsigned offset = p - (exon ? 66u : 1u);
            unsigned ref_len = (i == (kind == 2u ? 0u : 1u)) ? 2u : 1u;
            sources[i] = (duckvep_haplotype_source_t){101u+i,
                (const uint8_t *)(reverse ? "GG" : "CC"), (const uint8_t *)(reverse ? "GG" : "CC"),
                reverse ? f.ends[exon]-offset-ref_len+1u : f.starts[exon]+offset, 0u,
                (uint16_t)ref_len, (uint16_t)(3u-ref_len), 0u, 0u, 0u, 1u};
        }
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        for (unsigned k = 0u; k < 2u; k++) {
            unsigned i = reverse ? 1u-k : k;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, sources+i, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        }
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(0u, leaf.haplotype_so_mask);
        ASSERT_MEM_EQ("AAAAAAATGAAATAA", f.noncoding + 95u, 15u);
        for (size_t i = 15u; i < 95u; i++) ASSERT_EQ('C', f.noncoding[95u+i]);
        ASSERT_MEM_EQ((lower ? "atgaaataa" : "ATGAAATAA"), f.reference, 9u);
        ASSERT_MEM_EQ((lower ? "aaaaaa" : "AAAAAA"), f.flank, 6u);
        for (size_t i = 6u; i < 86u; i++) ASSERT_EQ(lower ? 'c' : 'C', f.flank[i]);
        ASSERT_MEM_EQ("ATGAAATAA", leaf.cds, 9u); ASSERT_MEM_EQ("MK*", leaf.protein, 3u);
        ASSERT_EQ(0u, leaf.edit_count); ASSERT_EQ(2u, leaf.contributor_count);
        for (size_t i = 0u; i < 2u; i++) {
            ASSERT_EQ(DUCKVEP_ROLE_APPLIED, leaf.contributors[i].role);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS, leaf.contributors[i].projection_status);
            ASSERT_EQ(101u <= leaf.contributors[i].source.event_id && leaf.contributors[i].source.event_id <= 102u, 1u);
            ASSERT_MEM_EQ((reverse ? "GG" : "CC"), leaf.contributors[i].source.ref, leaf.contributors[i].source.ref_len);
        }
        if (!kind) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE, leaf.nmd);
            ASSERT_EQ(0u, leaf.nmd_stop_valid); ASSERT_EQ(0u, leaf.nmd_junction_valid);
        } else {
            /* Stop = 6 pre-CDS bases + 3 codons. J = exon 1 length + its net delta. */
            unsigned junction = kind == 1u ? 66u : 64u;
            ASSERT_EQ(15u, leaf.nmd_stop_position1); ASSERT_EQ(junction, leaf.nmd_junction_position1);
            ASSERT_EQ(junction - 15u > 50u ? DUCKVEP_HAPLOTYPE_NMD_TRIGGER : DUCKVEP_HAPLOTYPE_NMD_ESCAPE, leaf.nmd);
        }
    }
    PASS();
}

TEST haplotype_stream_terminal_rna_insertion_updates_long_exon_exception(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) for (unsigned last = 406u; last <= 408u; last++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.strands[0] = reverse ? -1 : 1;
        f.starts[0] = reverse ? 700u : 100u; f.ends[0] = f.starts[0] + 79u;
        f.starts[1] = reverse ? 100u : 200u; f.ends[1] = f.starts[1] + last - 1u;
        f.model.start1 = (uint32_t[]){100u};
        f.model.end1 = (uint32_t[]){reverse ? 779u : f.ends[1]};
        f.model.cds_start1 = (uint32_t[]){reverse ? f.ends[1]-14u : 206u};
        f.model.cds_end1 = (uint32_t[]){reverse ? f.ends[1]-6u : 214u};
        f.exon_count[0] = 2u; f.exons.exon_count = 2u;
        f.cdna_starts[1] = 81u; f.cdna_ends[0] = 80u; f.cdna_ends[1] = 80u+last;
        memcpy(f.reference, "ATGAAATAA", 9u); f.lengths[0] = 9u;
        uint8_t flank[512], rna[1024];
        memset(flank, 'A', 86u); memset(flank+86u, 'C', last-15u);
        f.sequences.flank_bytes = flank; f.sequences.flank_bytes_len = 71u+last;
        f.sequences.pre_cds_offset = f.pre_flank_offsets; f.sequences.pre_cds_length = f.pre_flank_lengths;
        f.sequences.post_cds_offset = f.post_flank_offsets; f.sequences.post_cds_length = f.post_flank_lengths;
        f.pre_flank_lengths[0] = 86u; f.post_flank_offsets[0] = 86u; f.post_flank_lengths[0] = last-15u;
        f.sequences.flanks_complete = 1u; f.buffers.noncoding = rna; f.buffers.noncoding_capacity = sizeof(rna);
        duckvep_haplotype_source_t source = {101u, (const uint8_t *)(reverse ? "G" : "C"),
            (const uint8_t *)"CG", reverse ? f.starts[1] : f.ends[1], 0u, 1u, 2u, 0u, 0u, 0u, 1u};
        uint32_t tx = 0u; duckvep_carrier_key_t key = carrier_test_key(1u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_NMD_ESCAPE, leaf.nmd);
        ASSERT_EQ(95u, leaf.nmd_stop_position1); ASSERT_EQ(80u, leaf.nmd_junction_position1);
        ASSERT_EQ(last+1u > DUCKVEP_HAPLOTYPE_NMD_LONG_EXON_BASES,
            !!(leaf.nmd_exceptions & DUCKVEP_HAPLOTYPE_NMD_EXCEPTION_LONG_EXON));
        ASSERT_MEM_EQ("ATGAAATAA", leaf.cds, 9u); ASSERT_MEM_EQ("MK*", leaf.protein, 3u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
    }
    PASS();
}

TEST haplotype_stream_keeps_noncoding_raw_empty_slot_outside_strict_prediction(void) {
    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.strands[0] = reverse ? -1 : 1;
        f.model.end1 = (uint32_t[]){119u}; f.ends[0] = 119u; f.cdna_ends[0] = 20u;
        f.model.cds_start1 = (uint32_t[]){0u}; f.model.cds_end1 = (uint32_t[]){0u}; f.lengths[0] = 0u;
        memcpy(f.cdna, "AAAAAAATGAAATAACCCCC", 20u); f.cdna_lengths[0] = 20u;
        f.sequences.cdna_bytes = f.cdna; f.sequences.cdna_bytes_len = 20u;
        f.sequences.cdna_offset = f.cdna_offsets; f.sequences.cdna_length = f.cdna_lengths;
        f.sequences.cdna_provided = 1u;
        f.buffers.noncoding = f.noncoding; f.buffers.noncoding_capacity = sizeof(f.noncoding);
        duckvep_haplotype_source_t source = {101u,
            (const uint8_t *)(reverse ? "GGGGGTTATTTCATTTTTTT" : "AAAAAAATGAAATAACCCCC"),
            (const uint8_t *)"", 100u, 0u, 20u, 0u, UINT32_MAX, 1u, 1u, 0u};
        uint32_t tx = 0u;
        duckvep_raw_gt_t call;
        ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw((const uint8_t *)"1", 1u, 1u, &call));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(&f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(&f.stream, tx, 0u, &call, 1u));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        ASSERT_EQ_FMT((int)DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, (int)leaf.path_status, "%d");
        ASSERT_EQ(DUCKVEP_REASON_NON_STRICT_PHASE_POLICY, leaf.path_reason);
        ASSERT_EQ(0u, leaf.haplotype_so_mask);
        ASSERT_EQ(UINT32_MAX, leaf.contributors[0].source.allele_index);
        ASSERT_EQ(NULL, leaf.cds); ASSERT_EQ(NULL, leaf.protein);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_NMD_UNKNOWN, leaf.nmd);
        ASSERT_EQ(1u, leaf.contributor_count); ASSERT_EQ(101u, leaf.contributors[0].source.event_id);
    }
    PASS();
}

/* One raw record's interpretations (REF, ALTs, empty slot) for each sample, as the host plan routes them. */
static duckvep_haplotype_stream_status_t raw_proof_record(duckvep_haplotype_stream_t *s, uint64_t id,
    uint32_t pos1, const char *ref, const char *const *alts, uint32_t alt_count,
    const char *const *gts, size_t samples) {
    for (uint32_t a = 0u; a <= alt_count + 1u; a++) {
        uint32_t allele = a > alt_count ? UINT32_MAX : a;
        const char *alt = !a ? ref : a > alt_count ? "" : alts[a - 1u];
        duckvep_haplotype_source_t source = {id, (const uint8_t *)ref, (const uint8_t *)alt, pos1, 0u,
            (uint16_t)strlen(ref), (uint16_t)strlen(alt), allele, 1u, 0u, allele};
        duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_begin(s, &source);
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) status = duckvep_haplotype_stream_project(s, 0u);
        for (size_t i = 0u; status == DUCKVEP_HAPLOTYPE_STREAM_OK && i < samples; i++) {
            duckvep_raw_gt_t call;
            if (!gts[i]) continue;
            if (duckvep_phase_parse_vep_raw((const uint8_t *)gts[i], strlen(gts[i]), alt_count, &call) !=
                DUCKVEP_RAW_GT_OK) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
            status = duckvep_haplotype_stream_push_raw_call(s, 0u, (uint32_t)i, &call, 1u);
        }
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    return DUCKVEP_HAPLOTYPE_STREAM_OK;
}

/* Complete original diploid calls identify one unordered copy pair. */
TEST haplotype_stream_raw_single_het_has_two_real_copies(void) {
    static const char *const alts_a[] = {"G"}, *const alts_b[] = {"T"};
    static const char *const het[] = {"0/1", "1/0", "0|1", "1|0"};
    for (unsigned reverse = 0u; reverse < 2u; reverse++)
      for (unsigned order = 0u; order < 2u; order++)
        for (unsigned spelling = 0u; spelling < 4u; spelling++) {
            struct haplotype_stream_scene f;
            haplotype_stream_scene_prepare(&f, 1u);
            f.strands[0] = reverse ? -1 : 1;
            memcpy(f.reference, "ATGAAACCCTAA", 12u);
            const char *a[] = {order ? "1/1" : het[spelling]};
            const char *b[] = {order ? het[spelling] : "1|1"};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
            const char *const minus_a[] = {"A"}, *const minus_b[] = {"C"};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(&f.stream, 11u,
                104u, reverse ? "G" : "A", reverse ? minus_a : alts_a, 1u, a, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(&f.stream, 12u,
                107u, reverse ? "T" : "C", reverse ? minus_b : alts_b, 1u, b, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
            duckvep_haplotype_leaf_t leaf;
            unsigned lanes = 0u;
            while (duckvep_haplotype_stream_next(&f.stream, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.prediction_status);
                uint32_t id = leaf.carriers.first_call;
                for (uint32_t i = 0u; i < leaf.carriers.call_count; i++) {
                    const duckvep_carrier_call_t *c = duckvep_carriers_call(&f.stream.carriers, id);
                    lanes |= 1u << (c->key.lane - 1u);
                    unsigned alt_lane = spelling % 2u ? 1u : 2u;
                    ASSERT_MEM_EQ((c->key.lane == alt_lane ? "ATGAGACTCTAA" :
                        (order != reverse) ? "ATGAGACCCTAA" : "ATGAAACTCTAA"), leaf.cds, 12u);
                    id = c->next_leaf;
                }
                ASSERT_EQ(0u, leaf.ordered_replacements);
                ASSERT_MEM_EQ("ATGAAACCCTAA", f.reference, 12u);
            }
            ASSERT_EQ(3u, lanes);
        }
    PASS();
}

TEST haplotype_stream_raw_single_het_multialt_is_one_original_record(void) {
    static const char *const alts_a[] = {"T", "G"}, *const alts_b[] = {"T"};
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    memcpy(f.reference, "ATGAAACCCTAA", 12u);
    const char *a[] = {"0/2", "2/0", "1/2", "2/1"};
    const char *b[] = {"1/1", "1/1", "1/1", "1/1"};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(&f.stream, 11u, 104u, "A", alts_a, 2u, a, 4u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(&f.stream, 12u, 107u, "C", alts_b, 1u, b, 4u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
    duckvep_haplotype_leaf_t leaf;
    unsigned seen[4] = {0};
    while (duckvep_haplotype_stream_next(&f.stream, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        uint32_t id = leaf.carriers.first_call;
        for (uint32_t i = 0u; i < leaf.carriers.call_count; i++) {
            const duckvep_carrier_call_t *c = duckvep_carriers_call(&f.stream.carriers, id);
            duckvep_prediction_status_t status; duckvep_prediction_reason_t reason;
            duckvep_haplotype_carrier_prediction(&f.stream.carriers, &leaf, c, &status, &reason);
            ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, status);
            {
                static const char *const expected[4][2] = {
                    {"ATGAAACTCTAA", "ATGAGACTCTAA"}, {"ATGAGACTCTAA", "ATGAAACTCTAA"},
                    {"ATGATACTCTAA", "ATGAGACTCTAA"}, {"ATGAGACTCTAA", "ATGATACTCTAA"}
                };
                ASSERT_MEM_EQ(expected[c->key.sample_index][c->key.lane - 1u], leaf.cds, 12u);
                ASSERT_EQ(1u, c->copy_has_het); ASSERT_EQ(11u, c->copy_het_record);
            }
            seen[c->key.sample_index] |= 1u << (c->key.lane - 1u);
            id = c->next_leaf;
        }
    }
    for (unsigned i = 0u; i < 4u; i++) ASSERT_EQ(3u, seen[i]);
    PASS();
}

TEST haplotype_stream_raw_alt_alt_requires_both_siblings(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    memcpy(f.reference, "ATGAAACCCTAA", 12u);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
    duckvep_raw_gt_t call;
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw(
        (const uint8_t *)"1/2", 3u, 2u, &call));
    /* A direct caller projects ALT1 but never supplies the original ALT2. */
    duckvep_haplotype_source_t only = {11u, (const uint8_t *)"A", (const uint8_t *)"T",
        104u, 0u, 1u, 1u, 1u, 1u, 0u, 1u};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(&f.stream, &only));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(&f.stream, 0u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(
        &f.stream, 0u, 0u, &call, 1u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
    duckvep_haplotype_leaf_t leaf;
    unsigned lanes = 0u;
    while (duckvep_haplotype_stream_next(&f.stream, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        uint32_t id = leaf.carriers.first_call;
        for (uint32_t i = 0u; i < leaf.carriers.call_count; i++) {
            const duckvep_carrier_call_t *c = duckvep_carriers_call(&f.stream.carriers, id);
            duckvep_prediction_status_t status; duckvep_prediction_reason_t reason;
            duckvep_haplotype_carrier_prediction(&f.stream.carriers, &leaf, c, &status, &reason);
            ASSERT_EQ(DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, status);
            ASSERT_EQ(DUCKVEP_REASON_NON_STRICT_PHASE_POLICY, reason);
            lanes |= 1u << (c->key.lane - 1u);
            id = c->next_leaf;
        }
    }
    ASSERT_EQ(1u, lanes);
    PASS();
}

TEST haplotype_stream_raw_single_het_original_evidence_is_sample_local(void) {
    static const char *const bad[] = {"0/1", "1", "1/1/1", "./1", "|0/1", "/1|0", "0", "0/0/0", "./."};
    static const char *const alts_a[] = {"G"}, *const alts_b[] = {"T"};
    for (unsigned order = 0u; order < 2u; order++)
      for (size_t scenario = 0u; scenario < sizeof(bad) / sizeof(bad[0]); scenario++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        memcpy(f.reference, "ATGAAACCCTAA", 12u);
        const char *a[] = {"0/1", order ? bad[scenario] : "0/1", "0/0", "1/1"};
        const char *b[] = {"1/1", order ? "0/1" : bad[scenario], "1/1", "1/1"};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(&f.stream, 11u, 104u, "A", alts_a, 1u, a, 4u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(&f.stream, 12u, 107u, "C", alts_b, 1u, b, 4u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        duckvep_haplotype_leaf_t leaf;
        unsigned seen[4] = {0};
        while (duckvep_haplotype_stream_next(&f.stream, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            uint32_t id = leaf.carriers.first_call;
            for (uint32_t i = 0u; i < leaf.carriers.call_count; i++) {
                const duckvep_carrier_call_t *c = duckvep_carriers_call(&f.stream.carriers, id);
                duckvep_prediction_status_t status; duckvep_prediction_reason_t reason;
                duckvep_haplotype_carrier_prediction(&f.stream.carriers, &leaf, c, &status, &reason);
                ASSERT_EQ(c->key.sample_index == 1u ? DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT :
                    DUCKVEP_PREDICTION_PREDICTED, status);
                seen[c->key.sample_index] |= 1u << (c->key.lane - 1u);
                id = c->next_leaf;
            }
            ASSERT_MEM_EQ("ATGAAACCCTAA", f.reference, 12u);
        }
        for (unsigned i = 0u; i < 4u; i++) ASSERT_EQ(3u, seen[i]);
        ASSERT_EQ(0u, f.stream.carriers.call_count);
        ASSERT_EQ(0u, f.stream.carriers.prefix_count);
        ASSERT_EQ(0u, f.stream.carriers.transcript_count);
      }
    PASS();
}

/* Homozygous original calls prove two copies of the same raw path. */
TEST haplotype_stream_raw_homozygous_alt_paths_predict_only_proven_copies(void) {
    static const char *const a_alts[] = {"G", "T"}, *const b_alts[] = {"T"};
    /* Samples: 0 proof; 1 unphased het pair joins lane 2; 2 triploid 1/1/0; 3 missing slot; 4 haploid
     * one-slot GT; 5 prefixed spelling hides an ALT copy; 6 heterogeneous ALTs with a missing background; 7 other sample, site A only. */
    static const char *const at_a[] = {"1/1", "0/1", "1/1/0", "1/1", "1", "|0/1", "1/2", "1|1"};
    static const char *const at_b[] = {"1|1", "0/1", "1/1", "1/.", "1", "1/1", "1/.", "0/0"};
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    memcpy(f.reference, "ATGAAACCCTAA", 12u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f.model, &f.exons, &f.sequences, &f.buffers));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(s, 11u, 104u, "A", a_alts, 2u, at_a, 8u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(s, 12u, 107u, "C", b_alts, 1u, at_b, 8u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    unsigned predicted = 0u, refused = 0u, both = 0u, single = 0u, ab_views = 0u;
    uint32_t ab_leaf = 0u;
    while (duckvep_haplotype_stream_next(s, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        /* The conditional 1/. lane replays the same bytes; only the called path is the proof leaf. */
        int edited_ab = leaf.cds && !memcmp(leaf.cds, "ATGAGACTCTAA", 12u) &&
            leaf.evidence_flags == DUCKVEP_CARRIER_CALLED;
        int edited_a = leaf.cds && !memcmp(leaf.cds, "ATGAGACCCTAA", 12u);
        uint32_t id = leaf.carriers.first_call;
        for (uint32_t i = 0u; i < leaf.carriers.call_count; i++) {
            const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, id);
            ASSERT(c != NULL);
            duckvep_prediction_status_t status;
            duckvep_prediction_reason_t reason;
            duckvep_haplotype_carrier_prediction(&s->carriers, &leaf, c, &status, &reason);
            uint32_t sample = c->key.sample_index;
            if (status == DUCKVEP_PREDICTION_PREDICTED) {
                predicted |= 1u << (sample * 2u + c->key.lane - 1u);
                ASSERT(sample == 0u || sample == 7u);
                ASSERT_EQ((int)DUCKVEP_REASON_SUPPORTED_DOMAIN, (int)reason);
                ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), leaf.haplotype_so_mask);
                ASSERT(sample ? edited_a : edited_ab);
                ASSERT_EQ(0u, leaf.ordered_replacements);
                if (!sample) { ASSERT_MEM_EQ("MRL", leaf.protein, 3u); both++; }
                else { ASSERT_MEM_EQ("MRP", leaf.protein, 3u); single++; }
            } else {
                refused |= 1u << (sample * 2u + c->key.lane - 1u);
                ASSERT(sample != 0u && sample != 7u);
                if (status == DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT)
                    ASSERT_EQ((int)DUCKVEP_REASON_NON_STRICT_PHASE_POLICY, (int)reason);
            }
            /* The unproven carriers that end on the proven path keep the compatibility replay. */
            if (edited_ab && sample != 0u) ASSERT(status != DUCKVEP_PREDICTION_PREDICTED);
            id = c->next_leaf;
        }
        if (edited_ab) {
            /* Unproven carriers sharing the same bytes remain a separate view;
             * they neither qualify nor censor the observed pair. */
            ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
            ab_views++;
            if (ab_views == 1u) {
                ab_leaf = leaf.carriers.id;
                ASSERT_EQ(2u, leaf.carriers.call_count);
                ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.prediction_status);
            } else {
                ASSERT_EQ(ab_leaf, leaf.carriers.id);
                ASSERT_EQ(4u, leaf.carriers.call_count);
                ASSERT_EQ((int)DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, (int)leaf.prediction_status);
            }
        }
        if (edited_a) ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.prediction_status);
    }
    ASSERT_EQ(3u << 0u | 3u << 14u, predicted);
    ASSERT_EQ(0xffffu & ~(3u << 0u | 3u << 14u), refused);
    ASSERT_EQ(2u, both); ASSERT_EQ(2u, single); ASSERT_EQ(2u, ab_views);

    /* Strict decoded replay of the same events reaches the same CDS, protein and SO mask. */
    haplotype_stream_scene_prepare(&f, 1u);
    memcpy(f.reference, "ATGAAACCCTAA", 12u);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f.model, &f.exons, &f.sequences, &f.buffers));
    duckvep_haplotype_source_t strict[] = {
        {11u, (const uint8_t *)"A", (const uint8_t *)"G", 104u, 0u, 1u, 1u, 0u, 0u, 0u, 1u},
        {12u, (const uint8_t *)"C", (const uint8_t *)"T", 107u, 0u, 1u, 1u, 0u, 0u, 0u, 1u}};
    duckvep_carrier_key_t key = {0u, 0, 1u, 2u, 0u, 0u};
    uint32_t tx = 0u;
    for (unsigned i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, strict + i, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
    ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.prediction_status);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), leaf.haplotype_so_mask);
    ASSERT_MEM_EQ("ATGAGACTCTAA", leaf.cds, 12u); ASSERT_MEM_EQ("MRL", leaf.protein, 3u);
    PASS();
}

/* The proof sample's leaf view after two raw records at 104 (A>G) and 107 (C>T); *views counts the views of its path.
 * The view's borrowed CDS is copied to cds. */
static int raw_proof_view(struct haplotype_stream_scene *f, const char *const *at_a,
    const char *const *at_b, size_t samples, uint32_t sample, duckvep_haplotype_leaf_t *out, unsigned *views, uint8_t cds[12]) {
    static const char *const a_alts[] = {"G"}, *const b_alts[] = {"T"};
    duckvep_haplotype_stream_t *s = &f->stream;
    haplotype_stream_scene_prepare(f, 1u);
    memcpy(f->reference, "ATGAAACCCTAA", 12u);
    if (duckvep_haplotype_stream_init(s, &f->model, &f->exons, &f->sequences, &f->buffers) ||
        raw_proof_record(s, 11u, 104u, "A", a_alts, 1u, at_a, samples) ||
        raw_proof_record(s, 12u, 107u, "C", b_alts, 1u, at_b, samples) ||
        duckvep_haplotype_stream_finish(s) != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) return 0;
    int found = 0;
    uint32_t ids[16], n = 0u;
    duckvep_haplotype_leaf_t leaf;
    while (duckvep_haplotype_stream_next(s, &leaf) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        uint32_t id = leaf.carriers.first_call, zero = 0u;
        for (uint32_t i = 0u; i < leaf.carriers.call_count; i++) {
            const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, id);
            if (!c) return 0;
            zero += c->key.sample_index == sample;
            id = c->next_leaf;
        }
        /* A view is never mixed: the sample's view holds only its two lanes. */
        if (zero && (zero != 2u || leaf.carriers.call_count != 2u || found)) return 0;
        if (zero) {
            if (!leaf.cds) return 0;
            *out = leaf;
            memcpy(cds, leaf.cds, 12u);
            found = 1;
        }
        if (n == 16u) return 0;
        ids[n++] = leaf.carriers.id;
    }
    *views = 0u;
    for (uint32_t i = 0u; found && i < n; i++) *views += ids[i] == out->carriers.id;
    return found;
}

/* Context independence: adding another sample's unproven calls to the proof path, before or after, never
 * changes the complete sample's row; any omitted REF call that is not diploid leaves a copy unobserved
 * whichever side of the ALT calls it arrives, while diploid hom-REF does not. Independent witness: the
 * sample alone, whose prediction is the strict decoded replay checked above (MKP* -> MRL*). */
TEST haplotype_stream_raw_copy_proof_is_context_independent(void) {
    static const char *const hom_a[] = {"1/1", "0/1"}, *const hom_b[] = {"1/1", "0/1"};
    static const char *const het_a[] = {"0/1", "1/1"}, *const het_b[] = {"0/1", "1/1"};
    struct haplotype_stream_scene f;
    duckvep_haplotype_leaf_t alone, mixed;
    uint8_t cds[12];
    unsigned views;
    ASSERT(raw_proof_view(&f, hom_a, hom_b, 1u, 0u, &alone, &views, cds));
    ASSERT_EQ(1u, views);
    ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)alone.prediction_status);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), alone.haplotype_so_mask);
    ASSERT_MEM_EQ("ATGAGACTCTAA", cds, 12u);
    duckvep_haplotype_nmd_t nmd = alone.nmd;
    ASSERT(raw_proof_view(&f, hom_a, hom_b, 2u, 0u, &mixed, &views, cds));
    ASSERT_EQ(2u, views); /* Sample 1's lane 2 is the unproven view of the same path. */
    ASSERT_EQ((int)alone.prediction_status, (int)mixed.prediction_status);
    ASSERT_EQ((int)alone.prediction_reason, (int)mixed.prediction_reason);
    ASSERT_EQ(alone.haplotype_so_mask, mixed.haplotype_so_mask);
    ASSERT_EQ((int)nmd, (int)mixed.nmd);
    ASSERT_MEM_EQ("ATGAGACTCTAA", cds, 12u);
    /* Sample order does not matter: here the unproven sample has index 0 and the proof is sample 1. */
    ASSERT(raw_proof_view(&f, het_a, het_b, 2u, 1u, &mixed, &views, cds));
    ASSERT_EQ(2u, views);
    ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)mixed.prediction_status);
    ASSERT_EQ(alone.haplotype_so_mask, mixed.haplotype_so_mask);
    ASSERT_MEM_EQ("ATGAGACTCTAA", cds, 12u);

    static const char *const omitted[] = {"0", "0/0/0", "0/0"};
    static const char *const alt[] = {"1/1"};
    for (unsigned g = 0u; g < 3u; g++)
        for (unsigned order = 0u; order < 2u; order++) {
            const char *const *a = order ? omitted + g : alt, *const *b = order ? alt : omitted + g;
            ASSERT(raw_proof_view(&f, a, b, 1u, 0u, &mixed, &views, cds));
            ASSERT_EQ(1u, views);
            ASSERT_MEM_EQ((order ? "ATGAAACTCTAA" : "ATGAGACCCTAA"), cds, 12u);
            ASSERT_EQ((int)(g == 2u ? DUCKVEP_PREDICTION_PREDICTED : DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT),
                      (int)mixed.prediction_status);
        }

    /* A qualification-only push charges a carrier slot before any mutation and adds no event. */
    uint16_t chrom[] = {0u};
    uint32_t start[] = {1u}, end[] = {100u};
    duckvep_transcript_model_t model = {.chrom_id = chrom, .start1 = start, .end1 = end, .transcript_count = 1u};
    struct carrier_test_pool pool;
    duckvep_carrier_buffers_t buffers = carrier_test_buffers(&pool);
    buffers.call_capacity = 1u;
    duckvep_carriers_t carriers;
    uint32_t done;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&carriers, &model, &buffers));
    duckvep_carrier_key_t lane1 = {0u, 0, 1u, 2u, 0u, 0u, 0u}, lane2 = lane1;
    lane2.lane = 2u;
    lane2.copy_unproven = 1u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&carriers, 0u, 10u, 1u, &done));
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_push(&carriers, 0u, &lane1, 0u));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&carriers, 0u, &lane1, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, duckvep_carriers_push(&carriers, 0u, &lane2, 0u));
    ASSERT_EQ(1u, carriers.call_count);
    ASSERT_EQ(1u, carriers.prefix_count);
    lane1.copy_unproven = 1u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&carriers, 0u, &lane1, 0u));
    ASSERT_EQ(1u, carriers.prefix_count);
    ASSERT_EQ(1u, pool.calls[0].key.copy_unproven);
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, duckvep_carriers_finish(&carriers, &done));
    duckvep_carrier_leaf_t view;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&carriers, &view));
    ASSERT_EQ(1u, view.call_count);
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&carriers, &view));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&carriers));
    ASSERT_EQ(0u, carriers.call_count);
    ASSERT_EQ(0u, carriers.prefix_count);
    PASS();
}

/* Overlapping raw records keep ordered compatibility replay and its diagnostic even when homozygous;
 * a repeated call is an error; a refused capacity push does not change a carrier's copy qualification. */
TEST haplotype_stream_raw_homozygous_alt_near_neighbours(void) {
    static const char *const hom[] = {"1/1"};
    static const char *const ga[] = {"GA"}, *const t[] = {"T"};
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    memcpy(f.reference, "ATGAAACCCTAA", 12u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f.model, &f.exons, &f.sequences, &f.buffers));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(s, 11u, 104u, "AA", ga, 1u, hom, 1u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, raw_proof_record(s, 12u, 105u, "A", t, 1u, hom, 1u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
    ASSERT_EQ(1u, leaf.ordered_replacements);
    ASSERT_EQ(2u, leaf.carriers.call_count);
    ASSERT_EQ((int)DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, (int)leaf.prediction_status);
    ASSERT_EQ((int)DUCKVEP_REASON_NON_STRICT_PHASE_POLICY, (int)leaf.prediction_reason);
    ASSERT_EQ(0u, leaf.haplotype_so_mask);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));

    haplotype_stream_scene_prepare(&f, 1u);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f.model, &f.exons, &f.sequences, &f.buffers));
    duckvep_raw_gt_t call;
    ASSERT_EQ(DUCKVEP_RAW_GT_OK, duckvep_phase_parse_vep_raw((const uint8_t *)"1/1", 3u, 1u, &call));
    duckvep_haplotype_source_t source = {11u, (const uint8_t *)"A", (const uint8_t *)"G", 104u, 0u,
        1u, 1u, 1u, 1u, 0u, 1u};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_begin(s, &source));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_project(s, 0u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR, duckvep_haplotype_stream_push_raw_call(s, 0u, 0u, &call, 1u));
    ASSERT_EQ(DUCKVEP_CARRIERS_DUPLICATE_CALL, s->carrier_error);

    uint16_t chrom[] = {0u};
    uint32_t start[] = {1u}, end[] = {100u};
    duckvep_transcript_model_t model = {.chrom_id = chrom, .start1 = start, .end1 = end, .transcript_count = 1u};
    struct carrier_test_pool pool;
    duckvep_carrier_buffers_t buffers = carrier_test_buffers(&pool);
    buffers.prefix_capacity = 1u;
    duckvep_carriers_t carriers;
    uint32_t done;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&carriers, &model, &buffers));
    duckvep_carrier_key_t proven = {0u, 0, 1u, 2u, 0u, 0u, 0u}, unproven = proven;
    unproven.copy_unproven = 1u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&carriers, 0u, 10u, 1u, &done));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&carriers, 0u, &proven, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&carriers, 0u, 20u, 2u, &done));
    ASSERT_EQ(DUCKVEP_CARRIERS_PREFIX_FULL, duckvep_carriers_push(&carriers, 0u, &unproven, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0u, pool.calls[0].key.copy_unproven);
    unproven.copy_unproven = 2u;
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_push(&carriers, 0u, &unproven, DUCKVEP_CARRIER_CALLED));
    PASS();
}

TEST haplotype_stream_replays_joint_noncoding_cdna_without_cds(void) {
    struct haplotype_stream_scene f;
    uint32_t cds_none[] = {0u};
    uint32_t tx = 0u;
    duckvep_carrier_key_t key = carrier_test_key(1u);
    duckvep_haplotype_leaf_t leaf;
    duckvep_haplotype_source_t source[] = {
        {101u, (const uint8_t *)"C", (const uint8_t *)"A", 105u, 0u, 1u, 1u, 0u, 0u, 0u, 1u},
        {102u, (const uint8_t *)"T", (const uint8_t *)"C", 125u, 0u, 1u, 1u, 0u, 0u, 0u, 1u}
    };

    haplotype_stream_scene_prepare(&f, 1u);
    f.model.cds_start1 = f.model.cds_end1 = cds_none;
    f.model.end1 = (uint32_t[]){129u};
    f.exon_count[0] = 2u;
    f.exons.exon_count = 2u;
    f.starts[0] = 100u; f.ends[0] = 109u;
    f.starts[1] = 120u; f.ends[1] = 129u;
    f.cdna_starts[0] = 1u; f.cdna_ends[0] = 10u;
    f.cdna_starts[1] = 11u; f.cdna_ends[1] = 20u;
    memcpy(f.cdna, "AAAAACCCCCGGGGGTTTTT", 20u);
    f.sequences.cdna_bytes = f.cdna;
    f.sequences.cdna_bytes_len = 20u;
    f.cdna_lengths[0] = 20u;
    f.sequences.cdna_offset = f.cdna_offsets;
    f.sequences.cdna_length = f.cdna_lengths;
    f.sequences.cdna_provided = 1u;
    f.buffers.noncoding = f.noncoding;
    f.buffers.noncoding_capacity = sizeof(f.noncoding);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
    for (size_t i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
            haplotype_test_begin_candidates(&f.stream, source + i, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
            haplotype_test_push_called(&f.stream, &key));
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
        duckvep_haplotype_stream_finish(&f.stream));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
        duckvep_haplotype_stream_next(&f.stream, &leaf));
    ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON),
        leaf.haplotype_so_mask);
    ASSERT_EQ(NULL, leaf.cds);
    ASSERT_EQ(NULL, leaf.reference_cds);
    ASSERT_EQ(2u, leaf.contributor_count);
    ASSERT_EQ(101u, leaf.contributors[0].source.event_id);
    ASSERT_EQ(102u, leaf.contributors[1].source.event_id);
    ASSERT_EQ((int)DUCKVEP_ROLE_APPLIED, (int)leaf.contributors[0].role);
    ASSERT_EQ((int)DUCKVEP_ROLE_APPLIED, (int)leaf.contributors[1].role);
    ASSERT_MEM_EQ("AAAAACCCCCGGGGGTTTTT", f.noncoding, 20u);
    ASSERT_MEM_EQ("AAAAAACCCCGGGGGCTTTT", f.noncoding + 20u, 20u);
    PASS();
}

TEST haplotype_stream_prediction_is_keyed_and_retains_every_contributor(void) {
    /* Complete and split-domain carriers retain identical source operands in separate leaves. */
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    uint32_t tx = 0u;
    duckvep_haplotype_source_t source = {1u, (const uint8_t *)"A", (const uint8_t *)"C",
        100u, 0u, 1u, 1u, 0u, 0u, 0u, 2u};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
    int32_t gt[2] = {2, 0};
    uint8_t phase[2] = {0u, 1u};
    const duckvep_haplotype_phase_set_t one[1] = {{10, 1u}}, two[2] = {{5, 1u}, {10, 1u}};
    duckvep_haplotype_call_t call = {gt, phase, 0u, 2u, 2u, {10, 1u}, DUCKVEP_PHASE_STRICT};
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_call(s, tx, &call, one, 1u));
    call.sample_index = 1u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_push_call(s, tx, &call, two, 2u));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    unsigned split = 0u;
    for (unsigned i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(1u, leaf.carriers.call_count);
        ASSERT_EQ(1u, leaf.contributor_count);
        ASSERT_EQ(2u, leaf.contributors[0].source.alt_ordinal);
        ASSERT_EQ((int)DUCKVEP_ROLE_APPLIED, (int)leaf.contributors[0].role);
        ASSERT_EQ(1u, leaf.contributors[0].edit_count);
        ASSERT_EQ(1u, leaf.listed_edit_count);
        /* The scene CDS has no start codon, so it is read as if its start were not annotated. */
        ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
        ASSERT_EQ((int)DUCKVEP_REASON_SUPPORTED_DOMAIN, (int)leaf.path_reason);
        const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, leaf.carriers.first_call);
        duckvep_prediction_status_t status;
        duckvep_prediction_reason_t reason;
        ASSERT(c != NULL);
        duckvep_haplotype_carrier_prediction(&s->carriers, &leaf, c, &status, &reason);
        if (c->key.sample_index == 1u) {
            split++;
            ASSERT_EQ(1u, c->key.domain_split);
            ASSERT_EQ((int)DUCKVEP_PREDICTION_INCOMPLETE_INPUT, (int)status);
            ASSERT_EQ((int)DUCKVEP_REASON_CROSS_PS_UNRESOLVED, (int)reason);
        } else {
            ASSERT_EQ(0u, c->key.domain_split);
            ASSERT_EQ((int)leaf.path_status, (int)status);
            ASSERT_EQ((int)leaf.path_reason, (int)reason);
        }
        ASSERT_EQ(0u, c->next_leaf);
    }
    ASSERT_EQ(1u, split);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
    PASS();
}

/* ---- same-codon classifier ------------------------------------------------ */

static const char haplo_class_cds[] = "ATGGCTGAACTGAAACGCAGCTAA"; /* M A E L K R S *; codon 7 (AGC) has AGC->TCC pairs */

static char haplo_class_complement(char c) { return c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A'; }

/* Standard code, TCAG order, written independently of the kernel's tables. */
static char haplo_class_amino(const char *codon) {
    static const char table[] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    static const char order[] = "TCAG";
    unsigned index = 0u;
    for (unsigned i = 0u; i < 3u; i++) index = index * 4u + (unsigned)(strchr(order, codon[i]) - order);
    return table[index];
}

struct haplo_class_edit { unsigned pos; const char *ref, *alt; }; /* transcript orientation, 1-based */

/* Independent expectation for frame-preserving substitutions, decided on the edited cDNA alone:
 * 0 synonymous, 1 missense, 3 stop_gained (a new first stop before the terminator), 4 start_lost (the edited
 * CDS does not begin with ATG), 5 stop_lost (no stop left), 6 stop_retained (the terminator changed to another
 * stop, peptide identical). Edits are equal-length substitutions. */
static int haplo_class_expect_substitutions(const struct haplo_class_edit *edits, size_t count) {
    char alt[sizeof(haplo_class_cds)];
    memcpy(alt, haplo_class_cds, sizeof(alt));
    for (size_t i = 0u; i < count; i++) memcpy(alt + edits[i].pos - 1u, edits[i].alt, strlen(edits[i].alt));
    if (memcmp(alt, "ATG", 3u)) return 4;
    int same = 1;
    unsigned first_stop = 0u;
    for (unsigned codon = 0u; codon < 8u; codon++) {
        char a = haplo_class_amino(alt + 3u * codon), r = haplo_class_amino(haplo_class_cds + 3u * codon);
        if (a == '*' && !first_stop) first_stop = codon + 1u;
        if (codon != 7u) same &= a == r;
    }
    if (!first_stop) return 5;
    if (first_stop < 8u) return 3;
    if (same) return memcmp(alt + 21u, haplo_class_cds + 21u, 3u) ? 6 : 0;
    return 1;
}

/* Runs one lane of records (transcript-orientation edits) through the stream and returns the leaf. */
static duckvep_haplotype_stream_status_t haplo_class_run(struct haplotype_stream_scene *f, int reverse,
    const struct haplo_class_edit *edits, size_t count, uint16_t ploidy, duckvep_haplotype_leaf_t *leaf) {
    haplotype_stream_scene_prepare(f, 1u);
    f->sequences.cds_bytes = (const uint8_t *)haplo_class_cds;
    f->sequences.cds_bytes_len = 24u;
    f->lengths[0] = f->cdna_ends[0] = 24u;
    f->ends[0] = f->starts[0] + 23u;
    f->strands[0] = reverse ? -1 : 1;
    duckvep_haplotype_stream_t *s = &f->stream;
    duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_init(
        s, &f->model, &f->exons, &f->sequences, &f->buffers);
    if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    /* Genomic records ascend; a minus-strand transcript reverses the transcript-ordered list. */
    char ref[4][64], alt[4][64];
    unsigned pos[4], ref_len[4], alt_len[4];
    if (count > 4u) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    for (size_t k = 0u; k < count; k++) {
        const struct haplo_class_edit *e = &edits[reverse ? count - 1u - k : k];
        size_t rl = strlen(e->ref), al = strlen(e->alt);
        ref_len[k] = (unsigned)rl; alt_len[k] = (unsigned)al;
        if (!reverse) {
            memcpy(ref[k], e->ref, rl); memcpy(alt[k], e->alt, al); pos[k] = 99u + e->pos;
        } else {
            for (size_t i = 0u; i < rl; i++) ref[k][i] = haplo_class_complement(e->ref[rl - 1u - i]);
            for (size_t i = 0u; i < al; i++) alt[k][i] = haplo_class_complement(e->alt[al - 1u - i]);
            pos[k] = 124u - (e->pos + (unsigned)rl - 1u);
        }
    }
    duckvep_carrier_key_t key = {0u, 10, 1u, ploidy, 1u, 0u};
    uint32_t tx = 0u;
    for (size_t k = 0u; k < count; k++) {
        duckvep_haplotype_source_t source = {k + 1u, (const uint8_t *)ref[k], (const uint8_t *)alt[k],
            pos[k], 0u, (uint16_t)ref_len[k], (uint16_t)alt_len[k], 0u, 0u, 0u, 1u};
        status = haplotype_test_begin_candidates(s, &source, &tx, 1u);
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) status = haplotype_test_push_called(s, &key);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    status = duckvep_haplotype_stream_finish(s);
    if (status != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) return status;
    return duckvep_haplotype_stream_next(s, leaf);
}

TEST haplotype_unchanged_ambiguous_first_codon_is_not_a_changed_effect(void) {
    static const char cds[] = "NAGGCTGAACTGAAACGCAGCTAA";
    static const struct {
        unsigned pos[2];
        const char *ref[2], *alt[2];
        size_t count;
        uint64_t mask;
    } cases[] = {
        {{9u, 0u}, {"A", NULL}, {"G", NULL}, 1u, DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)},
        {{2u, 0u}, {"A", NULL}, {"G", NULL}, 1u, DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE)},
        {{2u, 9u}, {"A", "A"}, {"G", "G"}, 2u,
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) | DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)},
        /* The insertion yields NGA GGC TGA: an unknown codon and a frame-displaced early stop. */
        {{2u, 0u}, {"A", NULL}, {"GA", NULL}, 1u,
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) | DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
            DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)}
    };
    for (int reverse = 0; reverse < 2; reverse++) for (size_t n = 0u; n < sizeof(cases)/sizeof(cases[0]); n++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        f.sequences.cds_bytes = (const uint8_t *)cds;
        f.sequences.cds_bytes_len = sizeof(cds) - 1u;
        f.lengths[0] = f.cdna_ends[0] = 24u;
        f.ends[0] = 123u;
        f.strands[0] = reverse ? -1 : 1;
        uint64_t flags = DUCKVEP_TX_CDS_START_NF;
        uint8_t phase = 1u;
        f.model.flags = &flags;
        f.model.cds_phase_offset = &phase;
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        uint32_t tx = 0u;
        duckvep_carrier_key_t key = carrier_test_key(1u);
        for (size_t k = 0u; k < cases[n].count; k++) {
            size_t i = reverse ? cases[n].count - 1u - k : k;
            size_t rl = strlen(cases[n].ref[i]), al = strlen(cases[n].alt[i]);
            char ref[4], alt[4];
            for (size_t j = 0u; j < rl; j++) ref[j] = reverse ?
                haplo_class_complement(cases[n].ref[i][rl - 1u - j]) : cases[n].ref[i][j];
            for (size_t j = 0u; j < al; j++) alt[j] = reverse ?
                haplo_class_complement(cases[n].alt[i][al - 1u - j]) : cases[n].alt[i][j];
            duckvep_haplotype_source_t source = {k + 1u, (const uint8_t *)ref, (const uint8_t *)alt,
                reverse ? 124u - cases[n].pos[i] : 99u + cases[n].pos[i],
                0u, (uint16_t)rl, (uint16_t)al, 0u, 0u, 0u, 1u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        if (leaf.haplotype_so_mask != cases[n].mask ||
            leaf.path_status != DUCKVEP_PREDICTION_PREDICTED)
            fprintf(stderr, "[ambiguous codon] reverse=%d case=%zu status=%d reason=%d mask=%llu expected=%llu cds=%.*s protein=%.*s\n",
                reverse, n, leaf.path_status, leaf.path_reason, (unsigned long long)leaf.haplotype_so_mask,
                (unsigned long long)cases[n].mask, (int)leaf.cds_length, leaf.cds,
                (int)leaf.protein_length, leaf.protein);
        ASSERT_EQ(DUCKVEP_PREDICTION_PREDICTED, leaf.path_status);
        ASSERT_EQ(cases[n].mask, leaf.haplotype_so_mask);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
    }
    PASS();
}

TEST haplotype_annotated_non_atg_start_retains_source_and_product_effects(void) {
    static const struct {
        const char *start, *replacement;
        unsigned position;
        uint8_t table, open_start;
        duckvep_prediction_status_t status;
        uint64_t mask;
    } cases[] = {
        {"CTG", "T", 1u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED,
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED) | DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)},
        {"TTG", "C", 1u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED,
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED) | DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)},
        {"CTG", "A", 4u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_MISSENSE)},
        {"CTG", "T", 7u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)},
        {"CTG", "C", 2u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_START_LOST)},
        {"CTG", "A", 4u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_MISSENSE)},
        /* Without a declared code CTG is not a start, so the start is read as not annotated. */
        {"CTG", "A", 4u, 0u, 0u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_MISSENSE)},
        {"CTG", "A", 4u, 1u, 1u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_MISSENSE)},
        {"NTG", "A", 4u, 1u, 0u, DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, 0u},
        /* CCG is no start codon: the start is read as not annotated. */
        {"CCG", "A", 4u, 1u, 0u, DUCKVEP_PREDICTION_PREDICTED, DUCKVEP_SO(DUCKVEP_SO_MISSENSE)},
        {"CTG", "T", 1u, 1u, 0u, DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT, 0u}
    };
    for (int reverse = 0; reverse < 2; reverse++) for (size_t n = 0u; n < sizeof(cases)/sizeof(cases[0]); n++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        uint8_t cds[24];
        memcpy(cds, cases[n].start, 3u);
        memcpy(cds + 3u, "GCTGAACTGAAACGCAGCTAA", 21u);
        f.sequences.cds_bytes = cds;
        f.sequences.cds_bytes_len = sizeof(cds);
        f.lengths[0] = f.cdna_ends[0] = 24u;
        f.ends[0] = 123u;
        f.strands[0] = reverse ? -1 : 1;
        uint8_t table = cases[n].table;
        uint64_t flags = cases[n].open_start ? DUCKVEP_TX_CDS_START_NF : 0u;
        f.model.flags = &flags;
        f.sequences.codon_table = cases[n].table ? &table : NULL;
        uint32_t tx = 0u;
        char ref = (char)cds[cases[n].position - 1u], alt = cases[n].replacement[0];
        if (n == 10u) ref = 'A';
        if (reverse) { ref = haplo_class_complement(ref); alt = haplo_class_complement(alt); }
        duckvep_haplotype_source_t source = {n + 1u, (const uint8_t *)&ref, (const uint8_t *)&alt,
            reverse ? 124u - cases[n].position : 99u + cases[n].position,
            0u, 1u, 1u, 0u, 0u, 0u, 1u};
        duckvep_carrier_key_t key = carrier_test_key(1u);
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(&f.stream, &source, &tx, 1u));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(&f.stream, &key));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(&f.stream));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(&f.stream, &leaf));
        if (cases[n].status != leaf.path_status)
            fprintf(stderr, "[annotated initiation] reverse=%d case=%zu expected=%d actual=%d reason=%d cds=%.*s protein=%.*s\n",
                reverse, n, cases[n].status, leaf.path_status, leaf.path_reason,
                (int)leaf.cds_length, leaf.cds, (int)leaf.protein_length, leaf.protein);
        ASSERT_EQ((int)cases[n].status, (int)leaf.path_status);
        if (cases[n].status == DUCKVEP_PREDICTION_PREDICTED)
            ASSERT_EQ(cases[n].mask, leaf.haplotype_so_mask);
        if (n == 10u) ASSERT_EQ(NULL, leaf.reference_protein);
        else {
            /* Compatibility and coding references are distinct translation operands. */
            char compatibility_first = cases[n].start[1] == 'C' ? 'P' : 'M';
            char coding_first = cases[n].start[0] == 'N' ? 'X' : cases[n].start[1] == 'C' ? 'P' : 'L';
            ASSERT_EQ(compatibility_first, leaf.reference_protein[0]);
            ASSERT_EQ(coding_first, leaf.reference_coding_protein[0]);
        }
        ASSERT_MEM_EQ(cds, cases[n].start, 3u);
        ASSERT_EQ(n + 1u, leaf.contributors[0].source.event_id);
        ASSERT_EQ(ref, (char)leaf.contributors[0].source.ref[0]);
        ASSERT_EQ(alt, (char)leaf.contributors[0].source.alt[0]);
    }
    PASS();
}

TEST haplotype_same_codon_classifier_matches_independent_translation(void) {
    /* Every single and pair of substitutions over all 24 CDS positions, on both strands: identical peptide
     * is synonymous, a changed peptide missense, a new first stop before the terminator stop_gained, a lost
     * ATG start_lost, a terminator that is no longer a stop stop_lost and one that is another stop with an
     * identical peptide stop_retained. Nothing stays pending. The combined sequence is classified, never a
     * member alone. */
    static const char bases[] = "ACGT";
    unsigned decided[7] = {0u, 0u, 0u, 0u, 0u, 0u, 0u}, seen_cis_syn_of_missense = 0u;
    for (int reverse = 0; reverse < 2; reverse++) {
        for (unsigned p1 = 1u; p1 <= 24u; p1++) for (unsigned b1 = 0u; b1 < 4u; b1++) {
            if (bases[b1] == haplo_class_cds[p1 - 1u]) continue;
            char one[2] = {bases[b1], 0}, r1[2] = {haplo_class_cds[p1 - 1u], 0};
            for (unsigned p2 = p1; p2 <= 24u; p2++) for (unsigned b2 = 0u; b2 < 4u; b2++) {
                if (p2 == p1 && b2) break;
                if (p2 > p1 && bases[b2] == haplo_class_cds[p2 - 1u]) continue;
                char two[2] = {bases[b2], 0}, r2[2] = {haplo_class_cds[p2 - 1u], 0};
                struct haplo_class_edit edits[2] = {{p1, r1, one}, {p2, r2, two}};
                size_t count = p2 == p1 ? 1u : 2u;
                struct haplotype_stream_scene f;
                duckvep_haplotype_leaf_t leaf;
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplo_class_run(&f, reverse, edits, count, 2u, &leaf));
                int expected = haplo_class_expect_substitutions(edits, count);
                ASSERT_EQ(count, leaf.edit_count);
                static const int so_of[7] = {DUCKVEP_SO_SYNONYMOUS, DUCKVEP_SO_MISSENSE, 0, DUCKVEP_SO_STOP_GAINED,
                    DUCKVEP_SO_START_LOST, DUCKVEP_SO_STOP_LOST, DUCKVEP_SO_STOP_RETAINED};
                ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
                ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.prediction_status);
                uint64_t expected_mask = DUCKVEP_SO(so_of[expected]);
                char final[sizeof(haplo_class_cds)];
                memcpy(final, haplo_class_cds, sizeof(final));
                for (size_t e = 0u; e < count; e++) final[edits[e].pos - 1u] = edits[e].alt[0];
                unsigned stop = 8u;
                for (unsigned c = 0u; c < 8u; c++)
                    if (haplo_class_amino(final + c * 3u) == '*') { stop = c; break; }
                if (expected != 4) {
                    if (expected == 0) expected_mask = 0u;
                    for (unsigned c = 1u; c < stop && c < 7u; c++) {
                        const char *r = haplo_class_cds + 3u * c, *a = final + 3u * c;
                        if (!memcmp(r, a, 3u)) continue;
                        char ra = haplo_class_amino(r), aa = haplo_class_amino(a);
                        if (aa == '*' || ra == '*') continue;
                        expected_mask |= DUCKVEP_SO(ra == aa ? DUCKVEP_SO_SYNONYMOUS : DUCKVEP_SO_MISSENSE);
                    }
                }
                ASSERT_EQ(expected_mask, leaf.haplotype_so_mask);
                ASSERT_EQ((int)duckvep_so_impact(expected_mask), (int)duckvep_so_impact(leaf.haplotype_so_mask));
                decided[expected]++;
                /* A pair whose members are each missense but whose combination is synonymous. */
                if (count == 2u && !expected && (p1 - 1u) / 3u == (p2 - 1u) / 3u) {
                    struct haplo_class_edit solo[1] = {edits[0]};
                    if (haplo_class_expect_substitutions(solo, 1u) == 1) seen_cis_syn_of_missense++;
                }
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(&f.stream, &leaf));
            }
        }
    }
    ASSERT(decided[0] > 0u && decided[1] > 0u && decided[3] > 0u && decided[4] > 0u && decided[5] > 0u && decided[6] > 0u);
    ASSERT(seen_cis_syn_of_missense > 0u);
    PASS();
}

TEST haplotype_same_codon_classifier_names_pure_and_mixed_indels_and_frames(void) {
    /* CDS ATG GCT GAA CTG AAA CGC AGC TAA. mask 0 with status PREDICTED never occurs; every case is decided
     * (slice 5 decides the start/terminal cases that were eligible_classifier_pending). */
    const uint64_t FS = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), SG = DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        SL = DUCKVEP_SO(DUCKVEP_SO_STOP_LOST), ST = DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        SR = DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
    const struct { const char *name; size_t count; struct haplo_class_edit edits[3]; uint64_t mask; } cases[] = {
        {"insertion at a codon boundary", 1u, {{12u, "G", "GTCT"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)},
        {"reached inserted codon before gained stop", 2u, {{6u, "T", "TAAA"}, {13u, "A", "T"}},
            SG | DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)},
        {"insertion inside a codon", 1u, {{4u, "G", "GAAA"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)},
        {"whole-codon deletion", 1u, {{12u, "GAAA", "G"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION)},
        {"two insertions", 2u, {{6u, "T", "TGGG"}, {18u, "C", "CGGG"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)},
        {"two deletions", 2u, {{6u, "TGAA", "T"}, {15u, "ACGC", "A"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION)},
        {"insertion and deletion", 2u, {{12u, "G", "GTCT"}, {15u, "ACGC", "A"}}, DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING)},
        {"insertion and substitution", 2u, {{5u, "C", "G"}, {12u, "G", "GTCT"}},
            DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING) | DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) |
            DUCKVEP_SO(DUCKVEP_SO_MISSENSE)},
        {"in-frame replacement", 1u, {{7u, "GAACTG", "GGG"}}, DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING)},
        /* Frame opening: the displaced frame reads T AAA (stop at alt codon 5) inside the interval. */
        {"one-base insertion, early stop", 1u, {{12u, "G", "GT"}}, SG | FS},
        {"frame restored after that stop", 2u, {{12u, "G", "GT"}, {17u, "GC", "G"}}, SG | FS},
        /* One base deleted: the -1 frame has no stop before the CDS ends. */
        {"one-base deletion runs off the CDS (frame and abolished termination)", 1u, {{15u, "AC", "A"}}, FS | SL},
        /* Deletion then insertion: restored before any stop, peptide changed. */
        {"restored before termination", 2u, {{5u, "CT", "C"}, {8u, "A", "AG"}}, DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING)},
        /* A stop outside every displaced interval is stop_gained alone; the later frame edit is post-stop. */
        {"stop gained, frame edit after it", 2u, {{7u, "G", "T"}, {15u, "AC", "A"}}, SG},
        {"insertion before the stop codon", 1u, {{21u, "C", "CCCC"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)},
        {"insertion inside the stop codon abolishes it", 1u, {{22u, "T", "TGGG"}}, SL},
        {"start codon deletion", 1u, {{1u, "ATGG", "A"}}, ST},
        {"frame edit in the terminal codon reads the stop in the displaced frame", 1u, {{22u, "T", "TA"}}, FS},
        {"start codon SNV", 1u, {{1u, "A", "C"}}, ST},
        {"start codon suppresses a stop gained and a frame edit", 3u, {{1u, "A", "C"}, {7u, "G", "T"}, {15u, "AC", "A"}}, ST},
        {"terminator SNV to a sense codon", 1u, {{24u, "A", "C"}}, SL},
        {"terminator SNV to another stop", 1u, {{23u, "A", "G"}}, SR},
        {"terminator terminated earlier by a stop, then abolished", 2u, {{7u, "G", "T"}, {24u, "A", "C"}}, SG},
        {"stop codon inserted before the terminator", 1u, {{21u, "C", "CTAA"}}, SR},
        {"terminator deleted", 1u, {{21u, "CTAA", "C"}}, SL},
        {"last codon and terminator deleted", 1u, {{18u, "CAGCTAA", "C"}}, SL},
        {"start rescued by an insertion inside the codon", 1u, {{1u, "A", "ATGA"}}, DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)}};
    for (size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); i++) {
        struct haplotype_stream_scene f;
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ_FMT(DUCKVEP_HAPLOTYPE_STREAM_OK, haplo_class_run(&f, 0, cases[i].edits, cases[i].count, 2u, &leaf), "%d");
        if (cases[i].mask) {
            ASSERT_EQ_FMT((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status, "%d");
            if (cases[i].mask != leaf.haplotype_so_mask)
                fprintf(stderr, "[coding case] %s: expected=%llu observed=%llu dna=%.*s protein=%.*s\n",
                    cases[i].name, (unsigned long long)cases[i].mask,
                    (unsigned long long)leaf.haplotype_so_mask, (int)leaf.cds_length, leaf.cds,
                    (int)leaf.protein_length, leaf.protein);
            ASSERT_EQ((long long)cases[i].mask, (long long)leaf.haplotype_so_mask);
            ASSERT_EQ((int)(cases[i].mask & (FS | SG | SL | ST) ? DUCKVEP_IMPACT_HIGH :
                            cases[i].mask == SR ? DUCKVEP_IMPACT_LOW : DUCKVEP_IMPACT_MODERATE),
                      (int)duckvep_so_impact(leaf.haplotype_so_mask));
        } else {
            FAILm("every case is decided");
        }
        /* Every contributor stays attached, whatever its role. */
        ASSERT_EQ(cases[i].count, leaf.contributor_count);
    }
    PASS();
}

/* Frame classifier oracle: edits are single-base deletions, single-base insertions and substitutions on the
 * 24-base CDS. The expectation edits the cDNA, translates it with the standard code and walks the edits
 * base by base to flag displaced bases; it uses no production structure. */
enum { HF_DEL, HF_INS, HF_SUB };
struct hf_edit { int kind; unsigned p; char base; }; /* DEL: base p; INS: base after p; SUB: base p */

static int hf_expect(const struct hf_edit *edits, size_t count, uint64_t *mask, bool *post_stop) {
    /* Always decides (slice 5). Besides the displaced flag every edited base carries whether it holds the
     * reference terminator: retained bases of cDNA 22-24, and every base of an edit that reaches them (an
     * insertion after base 22 or 23, a deletion or substitution of base 22-24). */
    char alt[64]; int displaced[64], window[64], ancestry[64]; size_t n = 0u; int offset = 0; unsigned at = 1u;
    bool net_frame = false;
    size_t edit_start[4] = {0u, 0u, 0u, 0u}; /* 0-based position of each edit's first edited base in the alt CDS */
    for (size_t i = 0u; i < count; i++) {
        const struct hf_edit *e = &edits[i];
        unsigned start = e->kind == HF_INS ? e->p + 1u : e->p;
        for (; at < start; at++) {
            alt[n] = haplo_class_cds[at - 1u]; displaced[n] = offset % 3 != 0;
            window[n] = at >= 22u; ancestry[n++] = (int)at - 1;
        }
        edit_start[i] = n;
        if (e->kind == HF_INS) {
            alt[n] = e->base; displaced[n] = offset % 3 != 0 || (offset + 1) % 3 != 0;
            window[n] = e->p == 22u || e->p == 23u; ancestry[n++] = -1; offset++; net_frame = true;
        } else if (e->kind == HF_DEL) {
            offset--; at = e->p + 1u; net_frame = true;
        } else {
            alt[n] = e->base; displaced[n] = offset % 3 != 0;
            window[n] = e->p >= 22u; ancestry[n++] = (int)e->p - 1; at = e->p + 1u;
        }
    }
    for (; at <= 24u; at++) {
        alt[n] = haplo_class_cds[at - 1u]; displaced[n] = offset % 3 != 0;
        window[n] = at >= 22u; ancestry[n++] = (int)at - 1;
    }
    /* A deletion that removes the whole terminator with nothing flagged leaves an empty window at the end. */
    size_t window_begin = n + 1u;
    for (size_t i = 0u; i < n && window_begin == n + 1u; i++) if (window[i]) window_begin = i + 1u;
    size_t codons = n / 3u, first = 0u;
    for (size_t c = 1u; c <= codons && !first; c++) if (haplo_class_amino(alt + 3u * (c - 1u)) == '*') first = c;
    /* Roles are per edit: post-stop when the edit starts after the first stop codon of the edited sequence,
     * whatever interaction block it shares with earlier edits. */
    for (size_t i = 0u; i < count; i++) post_stop[i] = first && edit_start[i] >= 3u * first;
    if (n < 3u || memcmp(alt, "ATG", 3u)) { *mask = DUCKVEP_SO(DUCKVEP_SO_START_LOST); return 1; }
    if (!first) {
        *mask = DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
        if (n % 3u) *mask |= DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
        goto reached_codons;
    }
    size_t stop_start = 3u * first - 2u;
    bool same = first == 8u;
    for (size_t c = 0u; same && c + 1u < 8u; c++) same = haplo_class_amino(alt + 3u * c) == haplo_class_amino(haplo_class_cds + 3u * c);
    bool in_displaced = displaced[stop_start - 1u] || displaced[stop_start] || displaced[stop_start + 1u];
    if (stop_start < window_begin && !same) {
        *mask = DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
        if (in_displaced) *mask |= DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
    } else if (in_displaced) {
        *mask = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
    } else if (same) {
        *mask = DUCKVEP_SO(stop_start != window_begin || memcmp(alt + stop_start - 1u, haplo_class_cds + 21u, 3u) ?
            DUCKVEP_SO_STOP_RETAINED : DUCKVEP_SO_SYNONYMOUS);
    } else {
        *mask = net_frame ? DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING) : DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
    }
reached_codons:
    *mask &= ~DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    for (size_t c = 1u; c < (first ? first - 1u : codons) && c < 7u; c++) {
        size_t a = 3u * c;
        int r = ancestry[a];
        if (r < 0 || r % 3 || r + 2 >= 24 || ancestry[a + 1u] != r + 1 ||
            ancestry[a + 2u] != r + 2 || !memcmp(alt + a, haplo_class_cds + r, 3u)) continue;
        char before = haplo_class_amino(haplo_class_cds + r), after = haplo_class_amino(alt + a);
        if (before != '*' && after != '*')
            *mask |= DUCKVEP_SO(before == after ? DUCKVEP_SO_SYNONYMOUS : DUCKVEP_SO_MISSENSE);
    }
    return 1;
}

/* An indel next to a repeated base has several equivalent placements (VCF left alignment differs by strand
 * orientation), and the frame SO is normalized-edit-path-sensitive by contract. The exhaustive oracle uses only
 * edits with a unique placement: a deleted base differs from both neighbours, an inserted base from both flanks. */
static int hf_unambiguous(const struct hf_edit *e) {
    if (e->kind == HF_SUB) return 1;
    if (e->kind == HF_DEL) return haplo_class_cds[e->p - 1u] != haplo_class_cds[e->p - 2u] &&
        haplo_class_cds[e->p - 1u] != haplo_class_cds[e->p];
    return e->base != haplo_class_cds[e->p - 1u] && e->base != haplo_class_cds[e->p];
}

static void hf_to_stream(const struct hf_edit *e, char ref[3], char alt[3], struct haplo_class_edit *out) {
    char cds_base = haplo_class_cds[e->p - 1u];
    if (e->kind == HF_DEL) {
        ref[0] = haplo_class_cds[e->p - 2u]; ref[1] = cds_base; ref[2] = 0; alt[0] = ref[0]; alt[1] = 0;
        out->pos = e->p - 1u;
    } else if (e->kind == HF_INS) {
        ref[0] = cds_base; ref[1] = 0; alt[0] = cds_base; alt[1] = e->base; alt[2] = 0; out->pos = e->p;
    } else {
        ref[0] = cds_base; ref[1] = 0; alt[0] = e->base; alt[1] = 0; out->pos = e->p;
    }
    out->ref = ref; out->alt = alt;
}

static int hf_check(const struct hf_edit *edits, size_t count, int reverse, unsigned long *decided, unsigned long *pending,
                    unsigned long *seen_masks, unsigned long *post_stops) {
    char ref[4][3], alt[4][3];
    struct haplo_class_edit stream[4];
    for (size_t i = 0u; i < count; i++) hf_to_stream(&edits[i], ref[i], alt[i], &stream[i]);
    struct haplotype_stream_scene f;
    duckvep_haplotype_leaf_t leaf;
    if (haplo_class_run(&f, reverse, stream, count, 2u, &leaf) != DUCKVEP_HAPLOTYPE_STREAM_OK) return 0;
    uint64_t mask = 0u;
    bool post_stop[4] = {false, false, false, false};
    int decides = hf_expect(edits, count, &mask, post_stop);
    if (leaf.contributor_count != count) return 0;
    for (size_t i = 0u; i < count; i++) {
        uint64_t id = (reverse ? count - 1u - i : i) + 1u;
        for (size_t j = 0u; j < count; j++)
            if (leaf.contributors[j].source.event_id == id &&
                leaf.contributors[j].role != (post_stop[i] ? DUCKVEP_ROLE_POST_STOP : DUCKVEP_ROLE_APPLIED)) return 0;
        (*post_stops) += post_stop[i];
    }
    if (!decides || leaf.path_status != DUCKVEP_PREDICTION_PREDICTED || leaf.haplotype_so_mask != mask) {
        fprintf(stderr, "[frame witness] reverse=%d count=%zu edits=(%d,%u,%c),(%d,%u,%c),(%d,%u,%c) "
                "dna=%.*s expected=%llu observed=%llu status=%d\n", reverse, count,
                edits[0].kind, edits[0].p, edits[0].base,
                count > 1u ? edits[1].kind : -1, count > 1u ? edits[1].p : 0u, count > 1u ? edits[1].base : '-',
                count > 2u ? edits[2].kind : -1, count > 2u ? edits[2].p : 0u, count > 2u ? edits[2].base : '-',
                (int)leaf.cds_length, leaf.cds, (unsigned long long)mask,
                (unsigned long long)leaf.haplotype_so_mask, (int)leaf.path_status);
        fprintf(stderr, "[frame peptide] ref=%.*s alt=%.*s refcoding=%.*s joint=%zu edits=%zu\n",
                (int)leaf.reference_protein_length, leaf.reference_protein,
                (int)leaf.protein_length, leaf.protein,
                (int)leaf.reference_protein_length, leaf.reference_coding_protein,
                leaf.joint_edit_count, leaf.edit_count);
        return 0;
    }
    (*decided)++;
    const uint64_t FS = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT);
    if (mask == (FS | DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED))) *seen_masks |= 1u;
    else if (mask == FS) *seen_masks |= 2u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING)) *seen_masks |= 4u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)) *seen_masks |= 8u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)) *seen_masks |= 16u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_START_LOST)) *seen_masks |= 32u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_STOP_LOST)) *seen_masks |= 64u;
    else if (mask == (FS | DUCKVEP_SO(DUCKVEP_SO_STOP_LOST))) *seen_masks |= 128u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED)) *seen_masks |= 256u;
    else if (mask == DUCKVEP_SO(DUCKVEP_SO_MISSENSE)) *seen_masks |= 512u;
    (void)pending;
    return 1;
}

TEST haplotype_frame_classifier_matches_independent_translation_and_frame_walk(void) {
    /* Every single edit and every pair of edits (at least two bases apart) drawn from single-base deletions,
     * insertions and substitutions over every CDS position (the start codon, the body and the terminator), on
     * both strands, plus every triple of one insertion/deletion family. The classification follows the
     * translated edited sequence, the displaced-frame walk and the position of the reference terminator, never
     * the nominal net offset: the same +1/-1 history is rescued (protein altering), terminates early
     * (frameshift with stop_gained), is a stop outside the displaced frame, or runs off the CDS (frameshift
     * with stop_lost). Every path is decided. */
    struct hf_edit universe[128];
    size_t n = 0u;
    for (unsigned p = 1u; p <= 24u; p++) {
        if (p >= 2u) universe[n++] = (struct hf_edit){HF_DEL, p, 0};
        if (p < 24u) universe[n++] = (struct hf_edit){HF_INS, p, "GT"[p & 1u]};
        char pick = "ACGT"[(p * 7u) & 3u];
        universe[n++] = (struct hf_edit){HF_SUB, p, pick == haplo_class_cds[p - 1u] ? "CGTA"[(p * 7u) & 3u] : pick};
        if (p > 1u && p < 24u) universe[n++] = (struct hf_edit){HF_INS, p, "AC"[p & 1u]};
        if (p == 23u || p == 24u) universe[n++] = (struct hf_edit){HF_SUB, p, 'G'}; /* TAA -> TGA / TAG */
    }
    size_t kept = 0u;
    for (size_t i = 0u; i < n; i++) if (hf_unambiguous(&universe[i])) universe[kept++] = universe[i];
    n = kept;
    unsigned long decided = 0u, pending = 0u, seen = 0u, post_stops = 0u;
    for (int reverse = 0; reverse < 2; reverse++) {
        for (size_t a = 0u; a < n; a++) {
            ASSERT(hf_check(&universe[a], 1u, reverse, &decided, &pending, &seen, &post_stops));
            for (size_t b = 0u; b < n; b++) {
                const struct hf_edit pair[2] = {universe[a], universe[b]};
                /* Ascending order, at least three positions apart (closer indels at one site are the slice-2 overlap cases). */
                unsigned lo = pair[0].p, hi = pair[1].p;
                if (hi < lo + 3u) continue;
                ASSERT(hf_check(pair, 2u, reverse, &decided, &pending, &seen, &post_stops));
            }
        }
        for (unsigned a = 4u; a <= 20u; a += 2u) for (unsigned b = a + 3u; b <= 20u; b += 2u)
            for (unsigned c = b + 3u; c <= 20u; c += 2u) for (unsigned k = 0u; k < 27u; k++) {
                struct hf_edit triple[3];
                unsigned kinds[3] = {k % 3u, (k / 3u) % 3u, k / 9u}, pos[3] = {a, b, c};
                for (unsigned i = 0u; i < 3u; i++)
                    triple[i] = (struct hf_edit){(int)kinds[i], pos[i], kinds[i] == HF_SUB ? (haplo_class_cds[pos[i] - 1u] == 'T' ? 'G' : 'T') : "GC"[(a + i) & 1u]};
                if (!hf_unambiguous(&triple[0]) || !hf_unambiguous(&triple[1]) || !hf_unambiguous(&triple[2])) continue;
                ASSERT(hf_check(triple, 3u, reverse, &decided, &pending, &seen, &post_stops));
            }
    }
    /* Every outcome must occur: both frame outcomes, protein_altering, synonymous, stop_gained alone,
     * start_lost, stop_lost with and without frameshift, stop_retained and missense. Nothing is pending. */
    if (seen != 1023u) fprintf(stderr, "[frame classifier coverage] seen=%lu\n", seen);
    ASSERT_EQ(1023u, (unsigned)seen);
    ASSERT(decided > 1000u && pending == 0u && post_stops > 100u);
    fprintf(stderr, "[frame classifier coverage] decided=%lu pending=%lu post_stop_edits=%lu\n", decided, pending, post_stops);
    PASS();
}

TEST haplotype_same_codon_classifier_is_keyed_per_carrier(void) {
    struct haplotype_stream_scene f;
    duckvep_haplotype_leaf_t leaf;
    struct haplo_class_edit edit = {5u, "C", "G"};
    /* A lane of a complete phased call is one haplotype at any ploidy, so the carrier of a haploid,
     * triploid or tetraploid call takes the path's prediction, exactly as a diploid carrier does. */
    for (uint16_t ploidy = 1u; ploidy <= 4u; ploidy++) {
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplo_class_run(&f, 0, &edit, 1u, ploidy, &leaf));
        ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
        ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.prediction_status);
        ASSERT_EQ((int)DUCKVEP_REASON_SUPPORTED_DOMAIN, (int)leaf.prediction_reason);
        const duckvep_carrier_call_t *c = duckvep_carriers_call(&f.stream.carriers, leaf.carriers.first_call);
        ASSERT(c != NULL);
        ASSERT_EQ(ploidy, c->key.ploidy);
        duckvep_prediction_status_t status;
        duckvep_prediction_reason_t reason;
        duckvep_haplotype_carrier_prediction(&f.stream.carriers, &leaf, c, &status, &reason);
        ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)status);
        ASSERT_EQ((int)DUCKVEP_REASON_SUPPORTED_DOMAIN, (int)reason);
    }
    PASS();
}

TEST haplotype_stream_edit_provenance_survives_islands_sorting_and_blocks(void) {
    /* Exhaust eight genomic positions independently: unchanged, part of one
     * uploaded MNV, or a separate SNV. The direct digit/run oracle never uses
     * the production splitter or sort. Repeat on both transcript strands. */
    for (unsigned reverse = 0u; reverse < 2u; reverse++) for (unsigned code = 1u; code < 6561u; code++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        if (reverse) { memset(f.reference, 'T', sizeof(f.reference)); f.strands[0] = -1; }
        uint8_t digits[8], alternate[8], expected_cds[12];
        unsigned value = code, has_mnv = 0u;
        size_t source_count = 0u;
        memset(expected_cds, reverse ? 'T' : 'A', sizeof(expected_cds));
        for (unsigned i = 0u; i < 8u; i++) {
            digits[i] = (uint8_t)(value % 3u); value /= 3u;
            alternate[i] = digits[i] == 1u ? 'C' : 'A';
            has_mnv |= digits[i] == 1u;
            source_count += digits[i] == 2u;
            if (digits[i]) expected_cds[reverse ? 11u - i : i] =
                digits[i] == 1u ? (reverse ? 'G' : 'C') : (reverse ? 'C' : 'G');
        }
        source_count += has_mnv;
        f.edit_event_ids[8] = UINT64_C(0xface0123456789ab);
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        duckvep_carrier_key_t key = {0u, 0, 1u, 1u, 0u, 0u};
        uint32_t tx = 0u;
        if (has_mnv) {
            duckvep_haplotype_source_t source = {17u, (const uint8_t *)"AAAAAAAA", alternate,
                100u, 0u, 8u, 8u, 0u, 0u, 0u, 0u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
        for (unsigned i = 0u; i < 8u; i++) if (digits[i] == 2u) {
            duckvep_haplotype_source_t source = {UINT64_MAX - i, (const uint8_t *)"A",
                (const uint8_t *)"G", 100u + i, 0u, 1u, 1u, 0u, 0u, 0u, 0u};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_push_called(s, &key));
        }
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
        ASSERT_EQ(source_count, leaf.contributor_count);
        ASSERT_EQ(12u, leaf.cds_length);
        ASSERT_EQ(0, memcmp(expected_cds, leaf.cds, sizeof(expected_cds)));
        size_t edit = 0u;
        for (unsigned i = 0u; i < 8u;) {
            if (!digits[i]) { i++; continue; }
            unsigned begin = i++;
            if (digits[begin] == 1u) while (i < 8u && digits[i] == 1u) i++;
            ASSERT(edit < leaf.edit_count);
            size_t at = reverse ? leaf.edit_count - 1u - edit : edit;
            ASSERT_EQ(digits[begin] == 1u ? 17u : UINT64_MAX - begin, leaf.edit_event_ids[at]);
            ASSERT_EQ(reverse ? 13u - i : begin + 1u, f.edits[at].cds_start);
            ASSERT_EQ(i - begin, f.edits[at].ref_len);
            ASSERT_EQ(i - begin, f.edits[at].alt_len);
            edit++;
        }
        ASSERT_EQ(edit, leaf.edit_count);
        size_t covered = 0u;
        for (size_t i = 0u; i < leaf.block_count; i++) {
            ASSERT_EQ(covered, leaf.blocks[i].edit_begin);
            ASSERT(leaf.blocks[i].edit_count > 0u && leaf.blocks[i].edit_count <= edit - covered);
            covered += leaf.blocks[i].edit_count;
        }
        ASSERT_EQ(edit, covered);
        ASSERT_EQ(UINT64_C(0xface0123456789ab), f.edit_event_ids[8]);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    }
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    f.buffers.edit_event_ids = NULL;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG, duckvep_haplotype_stream_init(
        &f.stream, &f.model, &f.exons, &f.sequences, &f.buffers));
    PASS();
}

TEST haplotype_stream_calls_match_permitted_genotypes(void) {
    uint32_t cases = 0u;
    /* This is an additional finite lane, independent of the existing generated
     * carrier/edit properties. Enumerate genotype permutations rather than
     * calling the production phase reducer to construct expected evidence. */
    for (uint16_t ploidy = 1u; ploidy <= 4u; ploidy++) {
        for (uint32_t genotype = 0u; genotype < (1u << (2u * ploidy)); genotype++) {
            for (uint32_t bits = 0u; bits < (1u << ploidy); bits++) {
                int32_t alleles[4];
                uint8_t phase[4], order[4], possible[4] = {0}, missing_at[4] = {0};
                unsigned missing = 0u, unphased = 0u, homogeneous = 1u;
                for (uint16_t i = 0u; i < ploidy; i++) {
                    alleles[i] = (int32_t)((genotype >> (2u * i)) & 3u) - 1;
                    phase[i] = (uint8_t)((bits >> i) & 1u);
                    order[i] = (uint8_t)i;
                    if (alleles[i] < 0) missing++;
                    if (alleles[i] < 0 || alleles[i] != alleles[0]) homogeneous = 0u;
                    if (!phase[i]) unphased++;
                }
                do {
                    int allowed = 1;
                    for (uint16_t i = 0u; i < ploidy; i++)
                        if (phase[i] && order[i] != i) allowed = 0;
                    if (!allowed) continue;
                    for (uint16_t i = 0u; i < ploidy; i++) {
                        int32_t allele = alleles[order[i]];
                        possible[i] |= allele < 0 ? 7u : (uint8_t)(1u << allele);
                        if (allele < 0) missing_at[i] = 1u;
                    }
                } while (phase_next_permutation(order, ploidy));
                for (unsigned compat = 0u; compat < 2u; compat++) {
                    for (unsigned reverse = 0u; reverse < 2u; reverse++) {
                        for (uint32_t alt = 1u; alt <= 2u; alt++) {
                            cases++;
                            struct haplotype_stream_scene f;
                            haplotype_stream_scene_prepare(&f, 1u);
                            if (reverse) {
                                memset(f.reference, 'T', sizeof(f.reference));
                                f.strands[0] = -1;
                            }
                            duckvep_haplotype_stream_t *s = &f.stream;
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
                                s, &f.model, &f.exons, &f.sequences, &f.buffers));
                            duckvep_haplotype_phase_set_t sets[] = {{0, 0u}, {-3, 1u}, {0, 1u}, {17, 1u}};
                            size_t set_count = compat ? 1u : 4u;
                            duckvep_haplotype_call_t call = {alleles, bits ? phase : NULL, 0u,
                                alt, ploidy, {0, 1u}, compat ? DUCKVEP_PHASE_VEP_COMPAT : DUCKVEP_PHASE_STRICT};
                            duckvep_haplotype_source_t source = {71u, (const uint8_t *)"A",
                                (const uint8_t *)(alt == 1u ? "C" : "G"), 103u, 0u, 1u, 1u, 0u, 0u, 0u, 0u};
                            uint32_t tx = 0u;
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                haplotype_test_begin_candidates(s, &source, &tx, 1u));
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                                duckvep_haplotype_stream_push_call(s, tx, &call, sets, set_count));
                            uint8_t expected[4][4] = {{0}}, seen[4][4] = {{0}};
                            unsigned expected_calls = 0u;
                            if (compat) {
                                unsigned lane = 0u;
                                if (missing) for (unsigned i = 0u; i < ploidy; i++)
                                    expected[0][i] = DUCKVEP_CARRIER_MISSING;
                                for (unsigned i = 0u; i < ploidy; i++) {
                                    if (alleles[i] < 0) continue;
                                    if (alleles[i] == (int32_t)alt) expected[0][lane] |= DUCKVEP_CARRIER_CALLED;
                                    lane++;
                                }
                            } else {
                                for (unsigned i = 0u; i < ploidy; i++) {
                                    int resolved = ploidy == 1u || homogeneous || phase[i] || unphased == 1u ||
                                        (possible[i] & (possible[i] - 1u)) == 0u;
                                    uint8_t evidence = 0u;
                                    if (resolved) {
                                        if (alleles[i] < 0) evidence = DUCKVEP_CARRIER_MISSING;
                                        else if (alleles[i] == (int32_t)alt) evidence = DUCKVEP_CARRIER_CALLED;
                                    } else {
                                        if (possible[i] & (1u << alt)) evidence = DUCKVEP_CARRIER_UNPHASED;
                                        if (missing_at[i]) evidence |= DUCKVEP_CARRIER_MISSING;
                                    }
                                    for (unsigned set = 0u; set < 4u; set++) {
                                        if (set == 2u || ploidy == 1u || homogeneous || !bits)
                                            expected[set][i] = evidence;
                                    }
                                }
                            }
                            for (size_t set = 0u; set < set_count; set++)
                                for (unsigned i = 0u; i < ploidy; i++)
                                    if (expected[set][i]) expected_calls++;
                            ASSERT_EQ(expected_calls, s->carriers.call_count);
                            duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_finish(s);
                            unsigned seen_flags = 0u;
                            if (expected_calls) {
                                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, status);
                                duckvep_haplotype_leaf_t leaf;
                                while ((status = duckvep_haplotype_stream_next(s, &leaf)) ==
                                        DUCKVEP_HAPLOTYPE_STREAM_OK) {
                                    ASSERT(leaf.evidence_flags > 0u && leaf.evidence_flags < 8u);
                                    ASSERT_EQ(0u, seen_flags & (1u << leaf.evidence_flags));
                                    seen_flags |= 1u << leaf.evidence_flags;
                                    ASSERT_EQ(1u, leaf.contributor_count);
                                    ASSERT_EQ(71u, leaf.contributors[0].source.event_id);
                                    ASSERT_EQ(leaf.evidence_flags, leaf.contributors[0].evidence_flags);
                                    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
                                    ASSERT_EQ(!!(leaf.evidence_flags & DUCKVEP_CARRIER_CALLED), leaf.edit_count);
                                    if (leaf.evidence_flags != DUCKVEP_CARRIER_CALLED) {
                                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE, leaf.sequence_status);
                                        ASSERT(!leaf.cds && !leaf.protein);
                                        ASSERT_EQ(0u, leaf.cds_length); ASSERT_EQ(0u, leaf.protein_length);
                                    } else {
                                        uint8_t cds[12]; memset(cds, reverse ? 'T' : 'A', sizeof(cds));
                                        cds[reverse ? 8u : 3u] = (uint8_t)(reverse ? (alt == 1u ? 'G' : 'C') :
                                                                                                 (alt == 1u ? 'C' : 'G'));
                                        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
                                        ASSERT_EQ(12u, leaf.cds_length); ASSERT_EQ(4u, leaf.protein_length);
                                        ASSERT_EQ(0, memcmp(cds, leaf.cds, sizeof(cds)));
                                        const char *protein = reverse ? (alt == 1u ? "FFLF" : "FFFF") :
                                                                                        (alt == 1u ? "KQKK" : "KEKK");
                                        ASSERT_EQ(0, memcmp(protein, leaf.protein, 4u));
                                    }
                                    unsigned carriers = 0u;
                                    for (uint32_t id = leaf.carriers.first_call; id;) {
                                        const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, id);
                                        ASSERT(c); ASSERT_EQ(0u, c->key.sample_index);
                                        ASSERT_EQ(ploidy, c->key.ploidy);
                                        ASSERT(c->key.lane > 0u && c->key.lane <= ploidy);
                                        size_t set = 0u;
                                        while (set < set_count && (sets[set].present != c->key.phase_set_present ||
                                            (sets[set].present && sets[set].value != c->key.phase_set))) set++;
                                        ASSERT(set < set_count);
                                        ASSERT_EQ(0u, seen[set][c->key.lane - 1u]);
                                        seen[set][c->key.lane - 1u] = leaf.evidence_flags;
                                        carriers++; id = c->next_leaf;
                                    }
                                    ASSERT_EQ(carriers, leaf.carriers.call_count);
                                }
                                ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
                                status = duckvep_haplotype_stream_finish(s);
                            }
                            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
                            ASSERT_EQ(0, memcmp(expected, seen, sizeof(expected)));
                            ASSERT_EQ((seen_flags & (1u << DUCKVEP_CARRIER_CALLED)) ? 12u : 0u,
                                      s->translated_bases);
                            ASSERT_EQ(1u, s->projected_events);
                        }
                    }
                }
            }
        }
    }
    ASSERT_EQ(37440u, cases);
    PASS();
}

TEST haplotype_stream_preserves_later_phase_sets_and_uncertain_prefixes(void) {
    struct haplotype_stream_scene f;
    haplotype_stream_scene_prepare(&f, 1u);
    duckvep_haplotype_stream_t *s = &f.stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
        s, &f.model, &f.exons, &f.sequences, &f.buffers));
    const duckvep_haplotype_phase_set_t sets[] = {{10, 1u}, {20, 1u}};
    uint8_t expected[3][2][2][3] = {{{{0}}}}, seen[3][2][2] = {{{0}}};
    uint32_t tx = 0u;
    for (uint32_t event = 0u; event < 3u; event++) {
        duckvep_haplotype_source_t source = {event + 1u, (const uint8_t *)"A",
            (const uint8_t *)(event == 1u ? "G" : "C"), 100u + event, 0u, 1u, 1u, 0u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        for (uint32_t sample = 0u; sample < 3u; sample++) {
            int32_t gt[2] = {event == 2u ? 0 : 1, event == 1u ? 0 : 1};
            uint8_t phase[2] = {1u, 1u};
            if (event == 1u && sample == 2u) gt[1] = -1;
            if (event == 2u && sample == 1u) phase[0] = phase[1] = 0u;
            duckvep_haplotype_call_t call = {gt, phase, sample, 1u, 2u,
                event ? sets[event - 1u] : (duckvep_haplotype_phase_set_t){0, 0u}, DUCKVEP_PHASE_STRICT};
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_call(s, tx, &call, sets, 2u));
            for (unsigned set = 0u; set < 2u; set++) {
                for (unsigned lane = 0u; lane < 2u; lane++) {
                    if (!event) expected[sample][set][lane][event] = DUCKVEP_CARRIER_CALLED;
                    else if (event == 1u && set == 0u) {
                        if (lane == 0u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_CALLED;
                        else if (sample == 2u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_MISSING;
                    } else if (event == 2u) {
                        if (sample == 1u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_UNPHASED;
                        else if (set == 1u && lane == 1u) expected[sample][set][lane][event] = DUCKVEP_CARRIER_CALLED;
                    }
                }
            }
        }
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY, duckvep_haplotype_stream_finish(s));
    duckvep_haplotype_leaf_t leaf;
    duckvep_haplotype_stream_status_t status;
    unsigned leaves = 0u, calls = 0u;
    while ((status = duckvep_haplotype_stream_next(s, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        leaves++;
        for (uint32_t id = leaf.carriers.first_call; id;) {
            const duckvep_carrier_call_t *c = duckvep_carriers_call(&s->carriers, id);
            ASSERT(c && c->key.sample_index < 3u && c->key.phase_set_present);
            ASSERT(c->key.phase_set == 10 || c->key.phase_set == 20);
            ASSERT(c->key.lane == 1u || c->key.lane == 2u);
            ASSERT_EQ(2u, c->key.ploidy);
            unsigned sample = c->key.sample_index, set = c->key.phase_set == 20, lane = c->key.lane - 1u;
            ASSERT_EQ(0u, seen[sample][set][lane]); seen[sample][set][lane] = 1u;
            unsigned count = 0u;
            uint8_t evidence = 0u, cds[12]; memset(cds, 'A', sizeof(cds));
            for (unsigned event = 0u; event < 3u; event++) {
                uint8_t flags = expected[sample][set][lane][event];
                if (!flags) continue;
                ASSERT(count < leaf.contributor_count);
                ASSERT_EQ(event + 1u, leaf.contributors[count].source.event_id);
                ASSERT_EQ(flags, leaf.contributors[count].evidence_flags);
                evidence |= flags; count++;
                if (flags & DUCKVEP_CARRIER_CALLED) cds[event] = event == 1u ? 'G' : 'C';
            }
            ASSERT_EQ(count, leaf.contributor_count); ASSERT_EQ(evidence, leaf.evidence_flags);
            ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, leaf.projection_status);
            if (evidence == DUCKVEP_CARRIER_CALLED) {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK, leaf.sequence_status);
                ASSERT_EQ(12u, leaf.cds_length);
                ASSERT_EQ(0, memcmp(cds, leaf.cds, sizeof(cds)));
            } else {
                ASSERT_EQ(DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE, leaf.sequence_status);
                ASSERT(!leaf.cds && !leaf.protein);
            }
            calls++; id = c->next_leaf;
        }
    }
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, status);
    ASSERT_EQ(12u, calls); ASSERT_EQ(6u, leaves);
    ASSERT_EQ(3u, s->projected_events); ASSERT_EQ(36u, s->translated_bases);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_DONE, duckvep_haplotype_stream_finish(s));
    PASS();
}

TEST haplotype_stream_rejects_invalid_calls_and_latches_partial_broadcast(void) {
    for (unsigned scenario = 0u; scenario < 15u; scenario++) {
        struct haplotype_stream_scene f;
        haplotype_stream_scene_prepare(&f, 1u);
        if (scenario == 14u) f.buffers.carriers.call_capacity = 1u;
        duckvep_haplotype_stream_t *s = &f.stream;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(
            s, &f.model, &f.exons, &f.sequences, &f.buffers));
        uint32_t tx = 0u;
        duckvep_haplotype_source_t source = {1u, (const uint8_t *)"A", (const uint8_t *)"C",
            100u, 0u, 1u, 1u, 0u, 0u, 0u, 0u};
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, haplotype_test_begin_candidates(s, &source, &tx, 1u));
        int32_t gt[2] = {1, 0};
        uint8_t phase[2] = {1u, 1u};
        duckvep_haplotype_phase_set_t sets[] = {{0, 0u}, {10, 1u}};
        duckvep_haplotype_call_t call = {gt, phase, 0u, 1u, 2u, {10, 1u}, DUCKVEP_PHASE_STRICT};
        if (scenario == 0u) call.alleles = NULL;
        if (scenario == 1u) call.ploidy = 0u;
        if (scenario == 2u) call.alt_index = 0u;
        if (scenario == 3u) call.alt_index = UINT32_MAX;
        if (scenario == 4u) call.phase_set.present = 2u;
        if (scenario == 5u) call.policy = (duckvep_phase_policy_t)99;
        if (scenario == 6u) sets[1].present = 0u; /* Repeated absent set. */
        if (scenario == 7u) { sets[0] = sets[1]; sets[1].value = 5; }
        if (scenario == 8u) call.phase_set.value = 11;
        if (scenario == 9u) gt[1] = -2;
        if (scenario == 10u) phase[1] = 2u;
        if (scenario == 11u) call.policy = DUCKVEP_PHASE_VEP_COMPAT;
        if (scenario == 12u) tx = 1u; /* Not projected by begin. */
        if (scenario == 13u) {
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK,
                duckvep_haplotype_stream_push_call(s, tx, &call, sets, 2u));
            call.policy = DUCKVEP_PHASE_VEP_COMPAT;
        }
        if (scenario == 14u) gt[1] = 1; /* First lane lands, second exceeds the call pool. */
        duckvep_haplotype_stream_status_t expected = scenario == 14u
            ? DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR : DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        ASSERT_EQ(expected, duckvep_haplotype_stream_push_call(s, tx, &call, sets,
            scenario == 13u ? 1u : 2u));
        ASSERT_EQ(scenario >= 13u ? 1u : 0u, s->carriers.call_count);
        if (scenario == 14u) ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, s->carrier_error);
        ASSERT_EQ(expected, duckvep_haplotype_stream_push_call(s, tx, &call, NULL, 0u));
        ASSERT_EQ(expected, duckvep_haplotype_stream_finish(s));
        duckvep_haplotype_leaf_t leaf, zero = {0}; memset(&leaf, 0xa5, sizeof(leaf));
        ASSERT_EQ(expected, duckvep_haplotype_stream_next(s, &leaf));
        ASSERT_EQ(0, memcmp(&leaf, &zero, sizeof(leaf)));
    }
    PASS();
}

/* An independent dense bit matrix is the oracle, not the prefix index. Each
 * bit names a carried event; the empty mask is an implicit reference path. */
static int carrier_test_drain(duckvep_carriers_t *s, uint32_t tx,
                              const uint8_t masks[8], uint64_t event_base) {
    uint64_t seen_paths = 0u;
    unsigned seen_samples = 0u;
    duckvep_carrier_leaf_t leaf;
    duckvep_carriers_status_t status;
    while ((status = duckvep_carriers_next_leaf(s, &leaf)) == DUCKVEP_CARRIERS_OK) {
        duckvep_carrier_event_t events[6], before[6];
        size_t required = 0u;
        memset(events, 0xa5, sizeof events);
        memcpy(before, events, sizeof before);
        if (leaf.transcript_index != tx || !leaf.event_count || !leaf.call_count ||
            duckvep_carriers_leaf_events(s, leaf.id, events, 0u, &required) !=
                DUCKVEP_CARRIERS_OUTPUT_FULL || required != leaf.event_count ||
            memcmp(events, before, sizeof before) != 0 ||
            duckvep_carriers_leaf_events(s, leaf.id, events, 6u, &required) !=
                DUCKVEP_CARRIERS_OK) return 0;
        unsigned mask = 0u;
        for (size_t i = 0u; i < required; i++) {
            if (events[i].evidence_flags != DUCKVEP_CARRIER_CALLED ||
                events[i].event_id < event_base || events[i].event_id >= event_base + 6u ||
                (i && events[i].event_id <= events[i - 1u].event_id)) return 0;
            mask |= 1u << (unsigned)(events[i].event_id - event_base);
        }
        if (!mask || (seen_paths & (UINT64_C(1) << mask))) return 0;
        seen_paths |= UINT64_C(1) << mask;
        uint32_t count = 0u;
        for (uint32_t id = leaf.first_call; id;) {
            const duckvep_carrier_call_t *call = duckvep_carriers_call(s, id);
            if (!call || ++count > 8u || call->key.sample_index >= 8u) return 0;
            uint32_t sample = call->key.sample_index;
            duckvep_carrier_key_t key = carrier_test_key(sample);
            if ((seen_samples & (1u << sample)) || masks[sample] != mask ||
                call->key.lane != key.lane || call->key.ploidy != key.ploidy ||
                call->key.phase_set != key.phase_set ||
                call->key.phase_set_present != key.phase_set_present) return 0;
            seen_samples |= 1u << sample;
            id = call->next_leaf;
        }
        if (count != leaf.call_count) return 0;
    }
    if (status != DUCKVEP_CARRIERS_DONE) return 0;
    for (unsigned sample = 0u; sample < 8u; sample++) {
        if (!!(seen_samples & (1u << sample)) != !!masks[sample]) return 0;
    }
    return duckvep_carriers_release(s) == DUCKVEP_CARRIERS_OK;
}

TEST carrier_stream_lifetime_keys_and_capacity(void) {
    uint16_t chrom[] = {1u, 1u, 2u};
    uint32_t end[] = {3u, 8u, 4u};
    duckvep_transcript_model_t model = {0};
    model.transcript_count = 3u; model.chrom_id = chrom; model.end1 = end;
    struct carrier_test_pool pool, before_pool;
    duckvep_carrier_buffers_t b = carrier_test_buffers(&pool);
    duckvep_carriers_t s, before;
    uint32_t completed;
    duckvep_carrier_key_t key = carrier_test_key(0u);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 1u, 10u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 1u, &key, DUCKVEP_CARRIER_CALLED));
    /* NULL PS and an explicitly present zero are different carrier keys. */
    key.phase_set_present = 1u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(3u, s.call_count);
    ASSERT_EQ(2u, s.prefix_count);
    before = s; memcpy(&before_pool, &pool, sizeof pool);
    ASSERT_EQ(DUCKVEP_CARRIERS_DUPLICATE_CALL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    key.ploidy = 2u;
    ASSERT_EQ(DUCKVEP_CARRIERS_PLOIDY_CONFLICT, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0, memcmp(&before, &s, sizeof s));
    ASSERT_EQ(0, memcmp(&before_pool, &pool, sizeof pool));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 3u, 11u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY,
              duckvep_carriers_advance(&s, 1u, 4u, 12u, &completed));
    ASSERT_EQ(0u, completed);
    ASSERT_EQ(DUCKVEP_CARRIERS_INPUT_ORDER,
              duckvep_carriers_advance(&s, 1u, 5u, 12u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_release(&s));
    duckvep_carrier_leaf_t leaf;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(2u, leaf.call_count);
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    ASSERT_EQ(1u, s.transcript_count); ASSERT_EQ(1u, s.call_count); ASSERT_EQ(1u, s.prefix_count);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 4u, 12u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_INPUT_ORDER,
              duckvep_carriers_advance(&s, 1u, 4u, 11u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY,
              duckvep_carriers_advance(&s, 2u, 1u, 13u, &completed));
    ASSERT_EQ(1u, completed);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 2u, 1u, 13u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 2u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, duckvep_carriers_finish(&s, &completed));
    ASSERT_EQ(2u, completed);
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_finish(&s, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_finish(&s, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_INPUT_ORDER,
              duckvep_carriers_advance(&s, 2u, 2u, 14u, &completed));
    ASSERT_EQ(0u, s.transcript_count); ASSERT_EQ(0u, s.call_count); ASSERT_EQ(0u, s.prefix_count);

    b.transcript_capacity = 1u; b.transcript_buckets = 2u;
    b.call_capacity = 2u; b.call_buckets = 4u;
    b.prefix_capacity = 2u; b.prefix_buckets = 4u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 1u, 1u, &completed));
    for (uint32_t sample = 0u; sample < 2u; sample++) {
        key = carrier_test_key(sample);
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    }
    before = s; memcpy(&before_pool, &pool, sizeof pool);
    key = carrier_test_key(2u);
    ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_FULL, duckvep_carriers_push(&s, 1u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0, memcmp(&before, &s, sizeof s));
    ASSERT_EQ(0, memcmp(&before_pool, &pool, sizeof pool));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 2u, 2u, &completed));
    for (uint32_t sample = 0u; sample < 2u; sample++) {
        key = carrier_test_key(sample);
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    }
    ASSERT_EQ(2u, s.prefix_count); /* Reuse succeeds even when no slot is free. */
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 3u, 3u, &completed));
    before = s; memcpy(&before_pool, &pool, sizeof pool);
    ASSERT_EQ(DUCKVEP_CARRIERS_PREFIX_FULL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0, memcmp(&before, &s, sizeof s));
    ASSERT_EQ(0, memcmp(&before_pool, &pool, sizeof pool));
    b.call_capacity = 0u; b.call_buckets = 0u;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, 1u, 1u, &completed));
    ASSERT_EQ(DUCKVEP_CARRIERS_CALL_FULL, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
    ASSERT_EQ(0u, s.transcript_count); ASSERT_EQ(0u, s.prefix_count);
    b.transcript_buckets = 3u;
    ASSERT_EQ(DUCKVEP_CARRIERS_INVALID_ARG, duckvep_carriers_init(&s, &model, &b));
    PASS();
}

TEST carrier_stream_retains_split_phase_evidence(void) {
    uint16_t chrom[] = {1u};
    uint32_t end[] = {20u};
    duckvep_transcript_model_t model = {0};
    model.transcript_count = 1u; model.chrom_id = chrom; model.end1 = end;
    struct carrier_test_pool pool;
    duckvep_carrier_buffers_t b = carrier_test_buffers(&pool);
    for (unsigned initial = 0u; initial < 2u; initial++) {
        duckvep_carriers_t s;
        uint32_t completed;
        duckvep_carrier_key_t key = carrier_test_key(0u);
        key.ploidy = 2u; key.phase_set_present = 1u; key.phase_set = 7;
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
        for (unsigned i = 0u; i < 3u; i++) {
            key.domain_split = initial ? i == 0u : i == 1u;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, i + 1u, i + 1u, &completed));
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
        }
        ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, duckvep_carriers_finish(&s, &completed));
        duckvep_carrier_leaf_t leaf;
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
        ASSERT_EQ(1u, leaf.call_count);
        const duckvep_carrier_key_t *retained = &b.calls[leaf.first_call - 1u].key;
        ASSERT_EQ(1u, retained->domain_split);
        ASSERT_EQ(7, retained->phase_set);
        ASSERT_EQ(2u, retained->ploidy);
        ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
    }
    PASS();
}

TEST carrier_stream_partitions_phase_completeness(void) {
    uint16_t chrom[] = {1u};
    uint32_t end[] = {20u};
    duckvep_transcript_model_t model = {0};
    model.transcript_count = 1u; model.chrom_id = chrom; model.end1 = end;
    struct carrier_test_pool pool;
    duckvep_carrier_buffers_t b = carrier_test_buffers(&pool);
    for (unsigned initial = 0u; initial < 2u; initial++) {
        for (unsigned reverse = 0u; reverse < 2u; reverse++) {
            duckvep_carriers_t s;
            uint32_t completed;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
            for (unsigned i = 0u; i < 3u; i++) {
                ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 1u, i + 1u, i + 1u, &completed));
                for (unsigned j = 0u; j < 3u; j++) {
                    unsigned sample = reverse ? 2u - j : j;
                    duckvep_carrier_key_t key = carrier_test_key(sample);
                    key.domain_split = sample && (initial ? i == 0u : i == 1u);
                    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, 0u, &key, DUCKVEP_CARRIER_CALLED));
                }
            }
            ASSERT_EQ(initial ? 6u : 5u, s.prefix_count);
            ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, duckvep_carriers_finish(&s, &completed));
            for (unsigned i = 0u; i < 2u; i++) {
                duckvep_carrier_leaf_t leaf;
                duckvep_carrier_event_t events[3];
                size_t count;
                ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
                ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_leaf_events(&s, leaf.id, events, 3u, &count));
                ASSERT_EQ(3u, count);
                for (unsigned j = 0u; j < 3u; j++) {
                    ASSERT_EQ(j + 1u, events[j].event_id);
                    ASSERT_EQ(DUCKVEP_CARRIER_CALLED, events[j].evidence_flags);
                }
                uint32_t id = leaf.first_call;
                uint8_t split = b.calls[id - 1u].key.domain_split;
                ASSERT_EQ(split ? 2u : 1u, leaf.call_count);
                for (unsigned j = 0u; j < leaf.call_count; j++) {
                    const duckvep_carrier_call_t *call = duckvep_carriers_call(&s, id);
                    ASSERT(call != NULL);
                    ASSERT_EQ(split, call->key.domain_split);
                    ASSERT_EQ(split, call->key.sample_index != 0u);
                    id = call->next_leaf;
                }
                ASSERT_EQ(0u, id);
            }
            duckvep_carrier_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
            ASSERT_EQ(0u, s.prefix_count);
            ASSERT_EQ(0u, s.call_count);
            for (uint32_t i = 0u; i < b.prefix_buckets; i++) ASSERT_EQ(0u, b.prefix_index[i].id);
            for (uint32_t i = 0u; i < b.call_buckets; i++) ASSERT_EQ(0u, b.call_index[i].id);
        }
    }
    PASS();
}

struct carrier_matrix_case { uint8_t masks[2][8]; };

TEST carrier_stream_reuses_slots_after_many_transcripts(void) {
    enum { TRANSCRIPTS = 1000 };
    uint16_t chrom[TRANSCRIPTS] = {0};
    uint32_t end[TRANSCRIPTS];
    for (unsigned i = 0u; i < TRANSCRIPTS; i++) end[i] = i + 1u;
    duckvep_transcript_model_t model = {0};
    model.transcript_count = TRANSCRIPTS; model.chrom_id = chrom; model.end1 = end;
    struct carrier_test_pool pool;
    duckvep_carrier_buffers_t b = carrier_test_buffers(&pool);
    b.transcript_capacity = 1u; b.transcript_buckets = 2u;
    b.call_capacity = 8u; b.call_buckets = 16u;
    b.prefix_capacity = 1u; b.prefix_buckets = 2u;
    duckvep_carriers_t s;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    for (uint32_t tx = 0u; tx <= TRANSCRIPTS; tx++) {
        uint32_t completed;
        duckvep_carriers_status_t status = tx == TRANSCRIPTS
            ? duckvep_carriers_finish(&s, &completed)
            : duckvep_carriers_advance(&s, 0u, tx + 1u, UINT64_MAX - tx, &completed);
        if (tx) {
            ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, status);
            ASSERT_EQ(tx - 1u, completed);
            duckvep_carrier_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(8u, leaf.call_count);
            duckvep_carrier_event_t event;
            size_t required;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK,
                      duckvep_carriers_leaf_events(&s, leaf.id, &event, 1u, &required));
            ASSERT_EQ(1u, required);
            ASSERT_EQ(UINT64_MAX - (tx - 1u), event.event_id);
            ASSERT_EQ(DUCKVEP_CARRIER_CALLED, event.evidence_flags);
            unsigned seen = 0u;
            for (uint32_t id = leaf.first_call; id;) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s, id);
                ASSERT(call != NULL);
                uint32_t sample = call->key.sample_index - (tx - 1u) * 8u;
                ASSERT(sample < 8u);
                ASSERT(!(seen & (1u << sample)));
                seen |= 1u << sample;
                ASSERT_EQ(sample == 7u ? UINT16_MAX : sample + 1u, call->key.ploidy);
                ASSERT_EQ(call->key.ploidy, call->key.lane);
                id = call->next_leaf;
            }
            ASSERT_EQ(255u, seen);
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
            ASSERT_EQ(0u, s.transcript_count); ASSERT_EQ(0u, s.call_count); ASSERT_EQ(0u, s.prefix_count);
            status = tx == TRANSCRIPTS ? duckvep_carriers_finish(&s, &completed)
                : duckvep_carriers_advance(&s, 0u, tx + 1u, UINT64_MAX - tx, &completed);
        }
        if (tx == TRANSCRIPTS) {
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, status);
            break;
        }
        ASSERT_EQ(DUCKVEP_CARRIERS_OK, status);
        for (uint32_t sample = 0u; sample < 8u; sample++) {
            duckvep_carrier_key_t key = carrier_test_key(tx * 8u + sample);
            key.ploidy = sample == 7u ? UINT16_MAX : (uint16_t)(sample + 1u);
            key.lane = key.ploidy;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, tx, &key, DUCKVEP_CARRIER_CALLED));
        }
    }
    ASSERT_EQ(1u, s.peak_transcripts); ASSERT_EQ(8u, s.peak_calls); ASSERT_EQ(1u, s.peak_prefixes);
    PASS();
}

TEST carrier_stream_expiry_heap_matches_dense_active_set(void) {
    enum { TRANSCRIPTS = 256, SLOTS = 64 };
    uint16_t chrom[TRANSCRIPTS] = {0};
    uint32_t end[TRANSCRIPTS], start[TRANSCRIPTS], active[SLOTS];
    uint8_t live[TRANSCRIPTS] = {0};
    duckvep_carrier_transcript_t transcripts[SLOTS];
    duckvep_carrier_call_t calls[SLOTS];
    duckvep_carrier_prefix_t prefixes[SLOTS];
    duckvep_carrier_bucket_t tx_index[2 * SLOTS], call_index[2 * SLOTS], prefix_index[2 * SLOTS];
    duckvep_carrier_buffers_t b = {transcripts, calls, prefixes, active,
        tx_index, call_index, prefix_index, SLOTS, SLOTS, SLOTS,
        2 * SLOTS, 2 * SLOTS, 2 * SLOTS};
    /* Model ordinals are not opening order. First fill every slot with tied
     * ends; then interleave openings and expiry with up to 64 live positions. */
    for (uint32_t pos = 1u; pos <= TRANSCRIPTS; pos++) {
        uint32_t tx = (pos * 73u) % TRANSCRIPTS;
        start[tx] = pos;
        end[tx] = pos <= SLOTS ? SLOTS : pos + (tx * 29u) % SLOTS;
    }
    duckvep_transcript_model_t model = {0};
    model.transcript_count = TRANSCRIPTS; model.chrom_id = chrom; model.end1 = end;
    duckvep_carriers_t s;
    ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_init(&s, &model, &b));
    uint32_t drained = 0u;
    for (uint32_t pos = 1u; pos <= TRANSCRIPTS + 1u; pos++) {
        int eof = pos > TRANSCRIPTS;
        for (;;) {
            uint32_t expected = UINT32_MAX, completed;
            for (uint32_t tx = 0u; tx < TRANSCRIPTS; tx++) {
                if (live[tx] && (expected == UINT32_MAX || end[tx] < end[expected])) expected = tx;
            }
            duckvep_carriers_status_t status = eof ? duckvep_carriers_finish(&s, &completed)
                : duckvep_carriers_advance(&s, 0u, pos, pos, &completed);
            if (expected == UINT32_MAX || (!eof && end[expected] >= pos)) {
                ASSERT_EQ(eof ? DUCKVEP_CARRIERS_DONE : DUCKVEP_CARRIERS_OK, status);
                break;
            }
            ASSERT_EQ(DUCKVEP_CARRIERS_TRANSCRIPT_READY, status);
            ASSERT_EQ(expected, completed);
            duckvep_carrier_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(expected, leaf.transcript_index);
            ASSERT_EQ(1u, leaf.call_count);
            duckvep_carrier_event_t event;
            size_t required;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK,
                      duckvep_carriers_leaf_events(&s, leaf.id, &event, 1u, &required));
            ASSERT_EQ(1u, required); ASSERT_EQ(start[expected], event.event_id);
            ASSERT_EQ(DUCKVEP_CARRIER_CALLED, event.evidence_flags);
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&s, leaf.first_call);
            ASSERT(call != NULL); ASSERT_EQ(expected, call->key.sample_index);
            ASSERT_EQ(0u, call->next_leaf);
            ASSERT_EQ(DUCKVEP_CARRIERS_DONE, duckvep_carriers_next_leaf(&s, &leaf));
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_release(&s));
            live[expected] = 0u;
            drained++;
        }
        if (!eof) {
            uint32_t tx = (pos * 73u) % TRANSCRIPTS, completed;
            duckvep_carrier_key_t key = carrier_test_key(tx);
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_push(&s, tx, &key, DUCKVEP_CARRIER_CALLED));
            live[tx] = 1u;
            ASSERT_EQ(DUCKVEP_CARRIERS_OK, duckvep_carriers_advance(&s, 0u, pos, pos, &completed));
        }
        uint32_t live_count = 0u, earliest = UINT32_MAX;
        for (uint32_t tx = 0u; tx < TRANSCRIPTS; tx++) if (live[tx]) {
            live_count++;
            if (earliest == UINT32_MAX || end[tx] < end[earliest]) earliest = tx;
        }
        ASSERT_EQ(live_count, s.transcript_count);
        ASSERT_EQ(live_count, s.call_count); ASSERT_EQ(live_count, s.prefix_count);
        if (live_count) {
            ASSERT_EQ(earliest, transcripts[active[0] - 1u].transcript_index);
            for (uint32_t child = 1u; child < live_count; child++) {
                uint32_t a = transcripts[active[(child - 1u) / 2u] - 1u].transcript_index;
                uint32_t z = transcripts[active[child] - 1u].transcript_index;
                ASSERT(end[a] < end[z] || (end[a] == end[z] && a < z));
            }
        }
    }
    ASSERT_EQ(TRANSCRIPTS, drained);
    ASSERT_EQ(SLOTS, s.peak_transcripts); ASSERT_EQ(SLOTS, s.peak_calls); ASSERT_EQ(SLOTS, s.peak_prefixes);
    PASS();
}

static enum theft_alloc_res carrier_matrix_alloc(struct theft *t, void *env, void **instance) {
    (void)env;
    struct carrier_matrix_case *c = malloc(sizeof(*c));
    if (!c) return THEFT_ALLOC_ERROR;
    for (unsigned tx = 0u; tx < 2u; tx++) for (unsigned sample = 0u; sample < 8u; sample++) {
        c->masks[tx][sample] = (uint8_t)kprop_bounded(t, tx ? 64u : 16u);
    }
    /* Include cohort prefix sharing on every trial, not just by chance. */
    for (unsigned tx = 0u; tx < 2u; tx++) c->masks[tx][7] = c->masks[tx][0];
    *instance = c;
    return THEFT_ALLOC_OK;
}

static void carrier_matrix_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static enum theft_trial_res prop_carrier_stream_matches_dense_matrix(struct theft *t, void *arg) {
    (void)t;
    const struct carrier_matrix_case *c = arg;
    uint16_t chrom[2][2] = {{1u, 1u}, {2u, 2u}};
    uint32_t end[2][2] = {{3u, 6u}, {4u, 6u}};
    duckvep_transcript_model_t models[2] = {{0}, {0}};
    uint8_t expected[2][2][8];
    struct carrier_test_pool pools[2];
    duckvep_carriers_t streams[2];
    const uint64_t base = UINT64_MAX - 6u;
    for (unsigned run = 0u; run < 2u; run++) {
        models[run].transcript_count = 2u;
        models[run].chrom_id = chrom[run]; models[run].end1 = end[run];
        memcpy(expected[run], c->masks, sizeof c->masks);
        for (unsigned sample = 0u; sample < 8u; sample++) expected[run][0][sample] &= run ? 15u : 7u;
        duckvep_carrier_buffers_t b = carrier_test_buffers(&pools[run]);
        if (duckvep_carriers_init(&streams[run], &models[run], &b) != DUCKVEP_CARRIERS_OK)
            return THEFT_TRIAL_ERROR;
    }
    /* Interleave different models with identical local ordinals but different
     * contigs/transcript endpoints. One also resumes between every carrier row. */
    for (uint32_t pos = 1u; pos <= 6u; pos++) for (unsigned run = 0u; run < 2u; run++) {
        duckvep_carriers_t *s = &streams[run];
        uint32_t tx;
        duckvep_carriers_status_t status;
        while ((status = duckvep_carriers_advance(s, (uint16_t)(run + 1u), pos, base + pos - 1u, &tx)) ==
                   DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
            if (!carrier_test_drain(s, tx, expected[run][tx], base)) return THEFT_TRIAL_FAIL;
        }
        if (status != DUCKVEP_CARRIERS_OK) return THEFT_TRIAL_FAIL;
        for (tx = 0u; tx < 2u; tx++) for (uint32_t sample = 0u; sample < 8u; sample++) {
            if (!(expected[run][tx][sample] & (1u << (pos - 1u)))) continue;
            uint32_t completed;
            if (run && duckvep_carriers_advance(s, (uint16_t)(run + 1u), pos, base + pos - 1u, &completed) !=
                           DUCKVEP_CARRIERS_OK) return THEFT_TRIAL_FAIL;
            duckvep_carrier_key_t key = carrier_test_key(sample);
            if (duckvep_carriers_push(s, tx, &key, DUCKVEP_CARRIER_CALLED) != DUCKVEP_CARRIERS_OK) return THEFT_TRIAL_FAIL;
        }
    }
    for (unsigned run = 0u; run < 2u; run++) {
        duckvep_carriers_t *s = &streams[run];
        uint32_t tx;
        duckvep_carriers_status_t status;
        while ((status = duckvep_carriers_finish(s, &tx)) == DUCKVEP_CARRIERS_TRANSCRIPT_READY) {
            if (!carrier_test_drain(s, tx, expected[run][tx], base)) return THEFT_TRIAL_FAIL;
        }
        if (status != DUCKVEP_CARRIERS_DONE || s->transcript_count || s->call_count || s->prefix_count ||
            s->peak_transcripts > 2u || s->peak_calls > 16u || s->peak_prefixes > 64u)
            return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST carrier_stream_matches_dense_matrix_across_batches(void) {
    struct theft_type_info type = {.alloc = carrier_matrix_alloc, .free = carrier_matrix_free};
    struct theft_run_config cfg = {0};
    cfg.name = "sparse carrier paths == dense event matrix across input batches";
    cfg.prop1 = prop_carrier_stream_matches_dense_matrix;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* The dense genotype masks and direct genomic-to-CDS index formula are
 * independent of the replay store, projector and prefix traversal. The same
 * model-local ordinal and event IDs coexist on opposite strands/contigs. */
static enum theft_trial_res prop_haplotype_stream_matches_dense_models(struct theft *t, void *arg) {
    (void)t;
    const struct carrier_matrix_case *c = arg;
    struct haplotype_stream_scene f[2];
    for (unsigned run = 0u; run < 2u; run++) {
        haplotype_stream_scene_prepare(&f[run], 1u);
        f[run].chrom[0] = (uint16_t)(run + 1u);
        if (run) {
            memset(f[run].reference, 'T', sizeof(f[run].reference));
            f[run].strands[0] = -1;
        }
        if (duckvep_haplotype_stream_init(&f[run].stream, &f[run].model,
            &f[run].exons, &f[run].sequences, &f[run].buffers) != DUCKVEP_HAPLOTYPE_STREAM_OK)
            return THEFT_TRIAL_ERROR;
    }
    uint32_t tx = 0u;
    for (uint32_t event = 0u; event < 6u; event++) for (unsigned run = 0u; run < 2u; run++) {
        uint8_t ref = 'A', alt = 'C';
        duckvep_haplotype_source_t source = {UINT64_MAX - 5u + event, &ref, &alt,
            100u + event, (uint16_t)(run + 1u), 1u, 1u, 0u, 0u, 0u, 0u};
        duckvep_haplotype_stream_t *s = &f[run].stream;
        if (haplotype_test_begin_candidates(s, &source, &tx, 1u) != DUCKVEP_HAPLOTYPE_STREAM_OK)
            return THEFT_TRIAL_FAIL;
        ref = alt = 'N';
        for (uint32_t sample = 0u; sample < 8u; sample++) {
            if (!(c->masks[run][sample] & (1u << event))) continue;
            duckvep_carrier_key_t key = carrier_test_key(sample);
            if (haplotype_test_push_called(s, &key) != DUCKVEP_HAPLOTYPE_STREAM_OK)
                return THEFT_TRIAL_FAIL;
        }
    }
    for (unsigned run = 0u; run < 2u; run++) {
        duckvep_haplotype_stream_t *s = &f[run].stream;
        unsigned seen_samples = 0u, path_count = 0u;
        uint64_t seen_paths = 0u;
        uint8_t reference_peptide[5];
        duckvep_translation_t reference_translation;
        if (duckvep_translate_cds(f[run].reference, sizeof(f[run].reference),
                DUCKVEP_CODON_TABLE_STANDARD,
                reference_peptide, sizeof(reference_peptide), &reference_translation) !=
            DUCKVEP_TRANSLATION_OK) return THEFT_TRIAL_ERROR;
        duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_finish(s);
        if (s->carriers.call_count) {
            if (status != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) return THEFT_TRIAL_FAIL;
            duckvep_haplotype_leaf_t leaf;
            while ((status = duckvep_haplotype_stream_next(s, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                if (leaf.carriers.transcript_index || leaf.projection_status != DUCKVEP_CDS_EDIT_OK ||
                    leaf.sequence_status != DUCKVEP_HAPLOTYPE_OK || leaf.cds_length != 12u ||
                    leaf.protein_length != 4u || leaf.flags || !leaf.carriers.call_count)
                    return THEFT_TRIAL_FAIL;
                unsigned mask = 0u;
                uint8_t expected[12]; memset(expected, run ? 'T' : 'A', sizeof(expected));
                for (size_t i = 0u; i < leaf.contributor_count; i++) {
                    const duckvep_haplotype_source_t *source = &leaf.contributors[i].source;
                    if (source->event_id < UINT64_MAX - 5u || source->chrom_id != run + 1u ||
                        (i && source->event_id <= leaf.contributors[i - 1u].source.event_id) ||
                        source->ref[0] != 'A' || source->alt[0] != 'C') return THEFT_TRIAL_FAIL;
                    unsigned event = (unsigned)(source->event_id - (UINT64_MAX - 5u));
                    if (source->pos1 != 100u + event) return THEFT_TRIAL_FAIL;
                    mask |= 1u << event;
                    expected[run ? 11u - event : event] = run ? 'G' : 'C';
                }
                if (!mask || memcmp(expected, leaf.cds, sizeof(expected)) ||
                    (seen_paths & (UINT64_C(1) << mask))) return THEFT_TRIAL_FAIL;
                /* Consume the stream's actual completed leaf, including its
                 * ascending edits, without a second apply or translation. */
                duckvep_edit_set_t set = {s->buffers.edits, leaf.edit_count};
                duckvep_haplotype_result_t applied = {leaf.cds_length, 0, leaf.flags, leaf.edit_count};
                duckvep_coding_context_t context;
                if (duckvep_coding_context_open_replay(leaf.reference_cds, 12u, &set,
                        f[run].strands[0], DUCKVEP_CODON_TABLE_STANDARD, leaf.cds, &applied,
                        reference_peptide, &reference_translation, leaf.protein, &leaf.translation,
                        &context) != DUCKVEP_CODING_CONTEXT_OK ||
                    context.alt_cds != leaf.cds || context.alt_peptide != leaf.protein ||
                    context.ref_peptide != reference_peptide || context.applied_edits != leaf.edit_count)
                    return THEFT_TRIAL_FAIL;
                for (size_t block = 0u; block < leaf.block_count; block++) {
                    duckvep_sequence_delta_t facts;
                    if (duckvep_coding_context_block_delta_fill(&context, set.edits, set.count,
                            leaf.blocks + block, 0u, &facts) != DUCKVEP_CONTEXT_DELTA_OK || !facts.valid)
                        return THEFT_TRIAL_FAIL;
                }
                seen_paths |= UINT64_C(1) << mask;
                path_count++;
                unsigned count = 0u;
                for (uint32_t id = leaf.carriers.first_call; id;) {
                    const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->carriers, id);
                    if (!call || call->key.sample_index >= 8u || ++count > 8u) return THEFT_TRIAL_FAIL;
                    unsigned sample = call->key.sample_index;
                    duckvep_carrier_key_t key = carrier_test_key(sample);
                    if ((seen_samples & (1u << sample)) || c->masks[run][sample] != mask ||
                        call->key.lane != key.lane || call->key.ploidy != key.ploidy ||
                        call->key.phase_set_present != key.phase_set_present ||
                        call->key.phase_set != (key.phase_set_present ? key.phase_set : 0))
                        return THEFT_TRIAL_FAIL;
                    seen_samples |= 1u << sample;
                    id = call->next_leaf;
                }
                if (count != leaf.carriers.call_count) return THEFT_TRIAL_FAIL;
            }
            if (status != DUCKVEP_HAPLOTYPE_STREAM_DONE) return THEFT_TRIAL_FAIL;
            status = duckvep_haplotype_stream_finish(s);
        }
        if (status != DUCKVEP_HAPLOTYPE_STREAM_DONE || s->event_count || s->projection_count ||
            s->allele_count || s->completed_leaves != path_count || s->projected_events != 6u ||
            s->translated_bases != 12u * path_count) return THEFT_TRIAL_FAIL;
        for (unsigned sample = 0u; sample < 8u; sample++)
            if (!!(seen_samples & (1u << sample)) != !!c->masks[run][sample]) return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST haplotype_stream_matches_dense_models_across_batches(void) {
    struct theft_type_info type = {.alloc = carrier_matrix_alloc, .free = carrier_matrix_free};
    struct theft_run_config cfg = {0};
    cfg.name = "owned haplotype replay == dense genomic edits in coexisting models";
    cfg.prop1 = prop_haplotype_stream_matches_dense_models;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ---- NMD attribution, rule ejc50 ------------------------------------------------------
 * One transcript is laid out on either strand with a random UTR and 1 to 4 exons (a UTR-only exon is allowed) around a
 * random CDS. The oracle edits the spliced mRNA in transcript orientation, following every base's exon through the
 * edits, translates the edited CDS with an independent standard-code table and applies the rule; it shares nothing
 * with the kernel's exon walk. Records are written in genome orientation and VCF form, so the stream projects,
 * differences and orders them itself. Edits stay off the start and terminator codons (start_lost and terminator
 * windows are covered elsewhere), a margin from exon boundaries (an indel next to a boundary belongs to neither exon)
 * and away from one-step placement shifts (a normalized island must not move). */
#define NM_CDS_MAX 240u
#define NM_FLANK 8u
#define NM_INTRON 12u
#define NM_GENOME_MAX (NM_CDS_MAX + 2u * 40u + 4u * NM_INTRON + 2u * NM_FLANK + 64u)

enum { NM_SUB, NM_INS, NM_DEL };
struct nm_edit { int kind; unsigned p, n; char bases[4]; }; /* tx orientation; INS after base p; DEL bases p+1..p+n; p is a CDS position */
struct nm_case {
    unsigned reverse, cache, u5, u3, cds_len, exons, edit_count;
    unsigned ends[4]; /* cumulative exon ends in mRNA coordinates; the last is the mRNA length */
    char cds[NM_CDS_MAX + 1u];
    struct nm_edit edits[3];
};
struct nm_expect {
    unsigned premature, stop_valid, junction_valid, applied[3];
    int nmd; uint64_t stop, junction;
};

struct nmd_scene {
    duckvep_transcript_model_t model;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t sequences;
    uint16_t chrom[1], exon_count[1];
    uint32_t tstart[1], tend[1], cstart[1], cend[1], exon_offset[1], cds_cdna_start[1], cds_cdna_end[1], length[1];
    uint32_t ex_start[4], ex_end[4], ex_cdna_start[4], ex_cdna_end[4];
    uint64_t sequence_offset[1];
    int8_t strand[1];
    uint8_t genome[NM_GENOME_MAX], reference[NM_CDS_MAX];
    unsigned genome_length;
    struct carrier_test_pool carrier_pool;
    duckvep_haplotype_stored_event_t events[8];
    duckvep_haplotype_projection_t projections[16];
    duckvep_haplotype_contributor_t contributors[8];
    duckvep_haplotype_edit_t edits[8];
    uint64_t edit_event_ids[9];
    duckvep_haplotype_block_t blocks[8];
    duckvep_carrier_event_t leaf_events[8];
    uint8_t alleles[128], cds[NM_CDS_MAX + 64u], protein[NM_CDS_MAX + 64u];
    uint8_t reference_protein[NM_CDS_MAX + 64u], reference_coding_protein[NM_CDS_MAX + 64u];
    duckvep_haplotype_stream_buffers_t buffers;
    duckvep_haplotype_stream_t stream;
};

static unsigned nm_exon_of(const struct nm_case *c, unsigned t) { /* mRNA position t (1-based) -> exon 1..n */
    unsigned k = 0u;
    while (k + 1u < c->exons && t > c->ends[k]) k++;
    return k + 1u;
}
static unsigned nm_index(const struct nm_case *c, unsigned t) { /* transcript-orientation genome index of mRNA base t */
    return NM_FLANK + t + (nm_exon_of(c, t) - 1u) * NM_INTRON;
}
static unsigned nm_genomic(const struct nm_case *c, struct nmd_scene *f, unsigned t) {
    return c->reverse ? f->genome_length + 1u - nm_index(c, t) : nm_index(c, t);
}
static char nm_comp(char b) { return haplo_class_complement(b); }

/* Deterministic filler so an intron or flank never matches the transcript by accident. */
static char nm_filler(unsigned i) { return "GTCAAGCT"[(i * 5u + 3u) & 7u]; }

static void nm_scene_prepare(struct nmd_scene *f, const struct nm_case *c) {
    memset(f, 0, sizeof(*f));
    unsigned mrna = c->ends[c->exons - 1u];
    char tx[NM_GENOME_MAX];
    unsigned n = 0u;
    for (unsigned i = 0u; i < NM_FLANK; i++) tx[n++] = nm_filler(i);
    for (unsigned t = 1u; t <= mrna; t++) {
        if (t > 1u && nm_exon_of(c, t) != nm_exon_of(c, t - 1u))
            for (unsigned i = 0u; i < NM_INTRON; i++) tx[n++] = nm_filler(100u + t + i);
        tx[n++] = t <= c->u5 || t > c->u5 + c->cds_len ? nm_filler(200u + t) : c->cds[t - c->u5 - 1u];
    }
    for (unsigned i = 0u; i < NM_FLANK; i++) tx[n++] = nm_filler(300u + i);
    f->genome_length = n;
    for (unsigned i = 0u; i < n; i++) f->genome[i] = (uint8_t)(c->reverse ? nm_comp(tx[n - 1u - i]) : tx[i]);
    memcpy(f->reference, c->cds, c->cds_len);
    f->strand[0] = c->reverse ? -1 : 1;
    f->exon_count[0] = (uint16_t)c->exons;
    unsigned lo = nm_genomic(c, f, 1u), hi = nm_genomic(c, f, mrna);
    f->tstart[0] = lo < hi ? lo : hi; f->tend[0] = lo < hi ? hi : lo;
    lo = nm_genomic(c, f, c->u5 + 1u); hi = nm_genomic(c, f, c->u5 + c->cds_len);
    f->cstart[0] = lo < hi ? lo : hi; f->cend[0] = lo < hi ? hi : lo;
    for (unsigned k = 0u; k < c->exons; k++) {
        unsigned first = k ? c->ends[k - 1u] + 1u : 1u, last = c->ends[k];
        unsigned a = nm_genomic(c, f, first), b = nm_genomic(c, f, last);
        f->ex_start[k] = a < b ? a : b; f->ex_end[k] = a < b ? b : a;
        f->ex_cdna_start[k] = first; f->ex_cdna_end[k] = last;
    }
    f->length[0] = c->cds_len;
    f->cds_cdna_start[0] = c->u5 + 1u; f->cds_cdna_end[0] = c->u5 + c->cds_len;
    f->model.transcript_count = 1u; f->model.chrom_id = f->chrom;
    f->model.start1 = f->tstart; f->model.end1 = f->tend; f->model.cds_start1 = f->cstart; f->model.cds_end1 = f->cend;
    f->model.strand = f->strand; f->model.exon_offset = f->exon_offset; f->model.exon_count = f->exon_count;
    if (c->cache) { f->model.cds_cdna_start1 = f->cds_cdna_start; f->model.cds_cdna_end1 = f->cds_cdna_end; }
    f->exons.exon_count = c->exons; f->exons.start1 = f->ex_start; f->exons.end1 = f->ex_end;
    f->exons.cdna_start1 = f->ex_cdna_start; f->exons.cdna_end1 = f->ex_cdna_end;
    f->sequences.transcript_count = 1u; f->sequences.cds_bytes = f->reference; f->sequences.cds_bytes_len = c->cds_len;
    f->sequences.cds_offset = f->sequence_offset; f->sequences.cds_length = f->length;
    f->buffers = (duckvep_haplotype_stream_buffers_t){
        .carriers = carrier_test_buffers(&f->carrier_pool),
        .events = f->events, .projections = f->projections, .alleles = f->alleles,
        .event_capacity = 8u, .projection_capacity = 16u, .allele_capacity = sizeof(f->alleles),
        .leaf_events = f->leaf_events, .contributors = f->contributors, .edits = f->edits,
        .blocks = f->blocks, .edit_event_ids = f->edit_event_ids,
        .leaf_capacity = 8u, .edit_capacity = 8u, .cds = f->cds, .protein = f->protein,
        .cds_capacity = sizeof(f->cds), .protein_capacity = sizeof(f->protein),
        .reference_protein = f->reference_protein, .reference_coding_protein = f->reference_coding_protein,
        .reference_protein_capacity = sizeof(f->reference_protein)};
}

/* Genome-orientation VCF record of an edit (the same forms the R generator writes). */
struct nm_record { unsigned pos, ref_len, alt_len; char ref[8], alt[8]; unsigned edit; };
static void nm_record_of(const struct nm_case *c, struct nmd_scene *f, const struct nm_edit *e, struct nm_record *r) {
    unsigned p = c->u5 + e->p; /* mRNA coordinate of CDS base p */
    memset(r, 0, sizeof(*r));
    if (e->kind == NM_SUB) {
        unsigned g = nm_genomic(c, f, p);
        r->pos = g; r->ref_len = r->alt_len = 1u;
        r->ref[0] = (char)f->genome[g - 1u]; r->alt[0] = c->reverse ? nm_comp(e->bases[0]) : e->bases[0];
    } else if (e->kind == NM_DEL) {
        unsigned first = nm_genomic(c, f, p + 1u), last = nm_genomic(c, f, p + e->n);
        unsigned low = c->reverse ? last : first, high = c->reverse ? first : last;
        r->pos = c->reverse ? low - 1u : nm_genomic(c, f, p);
        unsigned from = r->pos, to = high; /* ref = anchor + deleted bases */
        r->ref_len = to - from + 1u; r->alt_len = 1u;
        for (unsigned i = 0u; i < r->ref_len; i++) r->ref[i] = (char)f->genome[from - 1u + i];
        r->alt[0] = r->ref[0];
    } else { /* insertion after CDS base p */
        unsigned g = c->reverse ? nm_genomic(c, f, p + 1u) : nm_genomic(c, f, p);
        r->pos = g; r->ref_len = 1u; r->ref[0] = (char)f->genome[g - 1u];
        if (!c->reverse) { r->alt[0] = r->ref[0]; for (unsigned i = 0u; i < e->n; i++) r->alt[1u + i] = e->bases[i]; }
        else { r->alt[0] = r->ref[0]; for (unsigned i = 0u; i < e->n; i++) r->alt[1u + i] = nm_comp(e->bases[e->n - 1u - i]); }
        r->alt_len = 1u + e->n;
    }
}

static unsigned nm_edit_start(const struct nm_case *c, const struct nm_edit *e) { /* reference mRNA coordinate of the first base or insertion point */
    return c->u5 + e->p + (e->kind == NM_SUB ? 0u : 1u);
}
static int nm_edit_change(const struct nm_edit *e) { return e->kind == NM_SUB ? 0 : e->kind == NM_INS ? (int)e->n : -(int)e->n; }

/* Independent oracle: edit the mRNA in transcript orientation, keep each base's exon and reference-CDS origin. */
static void nm_oracle(const struct nm_case *c, struct nm_expect *x, unsigned *first_stop_codon, int post[3]) {
    memset(x, 0, sizeof(*x));
    unsigned mrna = c->ends[c->exons - 1u], out_n = 0u, at = 1u;
    char base[NM_GENOME_MAX], out[NM_GENOME_MAX];
    unsigned lab[NM_GENOME_MAX], out_lab[NM_GENOME_MAX], origin[NM_GENOME_MAX], out_origin[NM_GENOME_MAX], start0[3];
    for (unsigned t = 1u; t <= mrna; t++) {
        base[t - 1u] = t <= c->u5 || t > c->u5 + c->cds_len ? 'N' : c->cds[t - c->u5 - 1u];
        lab[t - 1u] = nm_exon_of(c, t);
        origin[t - 1u] = t > c->u5 && t <= c->u5 + c->cds_len ? t - c->u5 - 1u : UINT32_MAX;
    }
    for (unsigned k = 0u; k < c->edit_count; k++) { /* edits are ascending by construction */
        const struct nm_edit *e = &c->edits[k];
        unsigned s = nm_edit_start(c, e), rl = e->kind == NM_INS ? 0u : e->kind == NM_SUB ? 1u : e->n;
        for (; at < s; at++) { out[out_n] = base[at - 1u]; out_lab[out_n] = lab[at - 1u]; out_origin[out_n++] = origin[at - 1u]; }
        start0[k] = out_n;
        if (e->kind != NM_DEL) {
            unsigned m = e->kind == NM_SUB ? 1u : e->n;
            for (unsigned i = 0u; i < m; i++) { out[out_n] = e->bases[i]; out_lab[out_n] = lab[s - 1u]; out_origin[out_n++] = UINT32_MAX; }
        }
        at = s + rl;
    }
    for (; at <= mrna; at++) { out[out_n] = base[at - 1u]; out_lab[out_n] = lab[at - 1u]; out_origin[out_n++] = origin[at - 1u]; }
    unsigned cds_n = out_n - c->u5 - c->u3, first = 0u;
    const char *edited = out + c->u5;
    for (unsigned codon = 0u; 3u * codon + 3u <= cds_n && !first; codon++)
        if (haplo_class_amino(edited + 3u * codon) == '*') first = codon + 1u;
    *first_stop_codon = first;
    for (unsigned k = 0u; k < c->edit_count; k++) post[k] = first && start0[k] - c->u5 >= 3u * first;
    if (cds_n < 3u || memcmp(edited, "ATG", 3u) || !first) { x->nmd = DUCKVEP_HAPLOTYPE_NMD_UNKNOWN; return; }
    unsigned term_start = cds_n;
    for (unsigned i = 0u; i < cds_n; i++)
        if (out_origin[c->u5 + i] != UINT32_MAX && out_origin[c->u5 + i] >= c->cds_len - 3u) { term_start = i; break; }
    bool same = first == c->cds_len / 3u;
    for (unsigned codon = 0u; same && codon + 1u < first; codon++)
        same = haplo_class_amino(edited + 3u * codon) == haplo_class_amino(c->cds + 3u * codon);
    x->premature = 3u * first - 3u < term_start && !same;
    if (!x->premature) { x->nmd = DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE; return; }
    x->stop_valid = 1u; x->stop = c->u5 + 3u * first;
    if (c->exons == 1u) { x->nmd = DUCKVEP_HAPLOTYPE_NMD_ESCAPE; return; }
    unsigned junction = 0u;
    for (unsigned i = 0u; i < out_n; i++) if (out_lab[i] == c->exons - 1u) junction = i + 1u;
    x->junction_valid = 1u; x->junction = junction;
    x->nmd = (int64_t)junction - (int64_t)x->stop > 50 ? DUCKVEP_HAPLOTYPE_NMD_TRIGGER : DUCKVEP_HAPLOTYPE_NMD_ESCAPE;
}

static int nm_compare_records(const void *a, const void *b) {
    const struct nm_record *x = a, *y = b;
    return x->pos < y->pos ? -1 : x->pos > y->pos;
}

static duckvep_haplotype_stream_status_t nm_run(const struct nm_case *c, struct nmd_scene *f, duckvep_haplotype_leaf_t *leaf) {
    nm_scene_prepare(f, c);
    duckvep_haplotype_stream_t *s = &f->stream;
    duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_init(s, &f->model, &f->exons, &f->sequences, &f->buffers);
    if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    struct nm_record records[3];
    for (unsigned k = 0u; k < c->edit_count; k++) { nm_record_of(c, f, &c->edits[k], &records[k]); records[k].edit = k; }
    qsort(records, c->edit_count, sizeof(records[0]), nm_compare_records);
    duckvep_carrier_key_t key = {0u, 10, 1u, 2u, 1u, 0u};
    uint32_t tx = 0u;
    for (unsigned k = 0u; k < c->edit_count; k++) {
        duckvep_haplotype_source_t source = {records[k].edit + 1u, (const uint8_t *)records[k].ref, (const uint8_t *)records[k].alt,
            records[k].pos, 0u, (uint16_t)records[k].ref_len, (uint16_t)records[k].alt_len, 0u, 0u, 0u, 1u};
        status = haplotype_test_begin_candidates(s, &source, &tx, 1u);
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) status = haplotype_test_push_called(s, &key);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    status = duckvep_haplotype_stream_finish(s);
    if (status != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) return status;
    return duckvep_haplotype_stream_next(s, leaf);
}

static const char *const nm_sense[] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG",
    "CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA",
    "AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};

/* ATG, sense codons with no run of three equal bases, and a stop: the domain of the v1 classifier. */
static void nm_random_cds(struct theft *t, unsigned codons, char *cds) {
    memcpy(cds, "ATG", 3u);
    for (unsigned i = 1u; i + 1u < codons; i++) {
        for (unsigned attempt = 0u;; attempt++) {
            const char *codon = nm_sense[kprop_bounded(t, sizeof(nm_sense) / sizeof(nm_sense[0]))];
            char *at = cds + 3u * i;
            memcpy(at, codon, 3u);
            bool run = (at[-1] == at[0] && at[-2] == at[-1]) || (at[-1] == at[0] && at[0] == at[1]) || (at[0] == at[1] && at[1] == at[2]);
            if (!run || attempt > 20u) break;
        }
    }
    memcpy(cds + 3u * (codons - 1u), (const char *[]){"TAA", "TAG", "TGA"}[kprop_bounded(t, 3u)], 3u);
    cds[3u * codons] = 0;
}

/* One indel is unambiguous when neither one-step shift of its placement gives the same edited sequence. */
static bool nm_unambiguous(const struct nm_case *c, const struct nm_edit *e) {
    if (e->kind == NM_SUB) return true;
    const char *s = c->cds; unsigned p = e->p; /* CDS base p is s[p - 1] */
    if (e->kind == NM_INS) return e->bases[0] != s[p] && e->bases[e->n - 1u] != s[p - 1u];
    return s[p - 1u] != s[p + e->n - 1u] && s[p] != s[p + e->n];
}
/* An indel keeps three mRNA bases of margin from every exon boundary, so it cannot belong to either exon. */
static bool nm_inside_one_exon(const struct nm_case *c, const struct nm_edit *e) {
    unsigned lo = c->u5 + e->p - 2u, hi = c->u5 + e->p + (e->kind == NM_SUB ? 0u : e->n) + 3u;
    return nm_exon_of(c, lo) == nm_exon_of(c, hi);
}

/* Draws one edit into e after CDS base `after` (the previous edit's last base), or fails when none fits. */
static bool nm_random_edit(struct theft *t, struct nm_case *c, unsigned after, unsigned limit, unsigned stop_bias, struct nm_edit *e) {
    unsigned lo = after + 6u;
    if (lo + 2u > limit) return false;
    unsigned kind = (unsigned)kprop_bounded(t, 4u);
    memset(e, 0, sizeof(*e));
    e->p = lo + (unsigned)kprop_bounded(t, limit - lo);
    if (kind <= 1u) {
        e->kind = NM_SUB; e->n = 1u;
        char ref = c->cds[e->p - 1u];
        if (kind == 0u || stop_bias) { /* try for a stop codon: move to a codon where one substitution makes TAA/TAG/TGA */
            for (unsigned tries = 0u; tries < 30u; tries++) {
                unsigned codon = (e->p - 1u) / 3u, at = 3u * codon;
                if (codon < 1u || at + 6u > c->cds_len) break;
                static const char *const stops[3] = {"TAA", "TAG", "TGA"};
                bool found = false;
                for (unsigned s = 0u; s < 3u && !found; s++) {
                    unsigned diff = 0u, where = 0u;
                    for (unsigned k = 0u; k < 3u; k++) if (c->cds[at + k] != stops[s][k]) { diff++; where = k; }
                    if (diff == 1u) { e->p = at + where + 1u; e->bases[0] = stops[s][where]; found = true; }
                }
                if (found) { ref = c->cds[e->p - 1u]; break; }
                e->p = lo + (unsigned)kprop_bounded(t, limit - lo);
            }
        }
        if (!e->bases[0] || e->bases[0] == ref) { e->bases[0] = "ACGT"[kprop_bounded(t, 4u)]; if (e->bases[0] == c->cds[e->p - 1u]) e->bases[0] = c->cds[e->p - 1u] == 'A' ? 'C' : 'A'; }
        return e->p >= lo && e->p <= limit && e->p > 3u;
    }
    if (kind == 2u) {
        e->kind = NM_INS; e->n = 1u + (unsigned)kprop_bounded(t, 3u);
        for (unsigned i = 0u; i < e->n; i++) e->bases[i] = "ACGT"[kprop_bounded(t, 4u)];
    } else {
        e->kind = NM_DEL; e->n = kprop_bounded(t, 2u) ? 3u : 1u + (unsigned)kprop_bounded(t, 3u);
        if (e->p + e->n + 3u > limit) return false;
    }
    return e->p > 4u && e->p <= limit && nm_unambiguous(c, e) && nm_inside_one_exon(c, e);
}

static unsigned nm_seen[16];
enum { NM_SEEN_TRIGGER, NM_SEEN_ESCAPE, NM_SEEN_UNKNOWN, NM_SEEN_NA, NM_SEEN_INTRONLESS, NM_SEEN_REVERSE, NM_SEEN_UTR, NM_SEEN_INDEL_SHIFT,
       NM_SEEN_POST, NM_SEEN_FLIP, NM_SEEN_EXONS4, NM_SEEN_RUN, NM_SEEN_ATTRIBUTED_POST };

/* Fills the mRNA/exon layout of a case: UTR lengths, cuts with at least eight bases per exon. */
static bool nm_random_layout(struct theft *t, struct nm_case *c, unsigned codons) {
    c->reverse = (unsigned)kprop_bounded(t, 2u); c->cache = (unsigned)kprop_bounded(t, 2u);
    c->cds_len = 3u * codons;
    c->u5 = kprop_bounded(t, 3u) ? (unsigned)kprop_bounded(t, 22u) : 0u;
    c->u3 = kprop_bounded(t, 3u) ? (unsigned)kprop_bounded(t, 22u) : 0u;
    c->exons = 1u + (unsigned)kprop_bounded(t, 4u);
    unsigned mrna = c->u5 + c->cds_len + c->u3, previous = 0u;
    for (unsigned k = 0u; k + 1u < c->exons; k++) {
        unsigned lo = previous + 8u, hi = mrna - 8u * (c->exons - 1u - k) - 1u;
        if (lo > hi) return false;
        c->ends[k] = lo + (unsigned)kprop_bounded(t, hi - lo + 1u);
        previous = c->ends[k];
    }
    c->ends[c->exons - 1u] = mrna;
    return true;
}

static enum theft_alloc_res nm_case_alloc(struct theft *t, void *env, void **instance) {
    (void)env;
    struct nm_case *c = calloc(1u, sizeof(*c));
    if (!c) return THEFT_ALLOC_ERROR;
    unsigned codons = 14u + (unsigned)kprop_bounded(t, 54u);
    if (!nm_random_layout(t, c, codons)) { c->cds_len = 0u; *instance = c; return THEFT_ALLOC_OK; }
    nm_random_cds(t, codons, c->cds);
    unsigned count = 1u + (unsigned)kprop_bounded(t, 3u), after = 3u, limit = c->cds_len - 12u;
    for (unsigned k = 0u; k < count; k++) {
        struct nm_edit e;
        unsigned bias = !k && kprop_bounded(t, 2u);
        bool ok = false;
        for (unsigned tries = 0u; tries < 20u && !ok; tries++) ok = nm_random_edit(t, c, after, limit, bias, &e);
        if (!ok) break;
        c->edits[c->edit_count++] = e;
        after = e.p + (e.kind == NM_DEL ? e.n : 1u);
    }
    *instance = c;
    return THEFT_ALLOC_OK;
}
static void nm_case_free(void *instance, void *env) { (void)env; free(instance); }

static enum theft_trial_res prop_nmd_matches_edited_geometry_oracle(struct theft *t, void *arg) {
    (void)t;
    const struct nm_case *c = arg;
    if (!c->cds_len || !c->edit_count) return THEFT_TRIAL_SKIP;
    struct nmd_scene *f = malloc(sizeof(*f));
    if (!f) return THEFT_TRIAL_ERROR;
    duckvep_haplotype_leaf_t leaf;
    struct nm_expect x;
    unsigned first;
    int post[3];
    nm_oracle(c, &x, &first, post);
    enum theft_trial_res result = THEFT_TRIAL_FAIL;
    if (nm_run(c, f, &leaf) != DUCKVEP_HAPLOTYPE_STREAM_OK || leaf.path_status != DUCKVEP_PREDICTION_PREDICTED ||
        leaf.contributor_count != c->edit_count || (int)leaf.nmd != x.nmd ||
        leaf.nmd_stop_valid != x.stop_valid || leaf.nmd_junction_valid != x.junction_valid ||
        (x.stop_valid && leaf.nmd_stop_position1 != x.stop) || (x.junction_valid && leaf.nmd_junction_position1 != x.junction) ||
        ((leaf.haplotype_so_mask & DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)) != 0u) != (x.premature != 0u)) goto done;
    for (unsigned k = 0u; k < c->edit_count; k++) {
        bool found = false;
        for (size_t j = 0u; j < leaf.contributor_count; j++) if (leaf.contributors[j].source.event_id == k + 1u) {
            found = true;
            if (leaf.contributors[j].role != (post[k] ? DUCKVEP_ROLE_POST_STOP : DUCKVEP_ROLE_APPLIED)) goto done;
            /* Attribution (applied, or post_stop and moved J) equals the oracle's: applied edits plus the post-stop
             * indels that change length at or before the penultimate exon's last base. Roles are unchanged. */
            if (x.junction_valid) {
                int moved = nm_edit_change(&c->edits[k]) && nm_edit_start(c, &c->edits[k]) <= c->ends[c->exons - 2u];
                int attributed = leaf.contributors[j].role == DUCKVEP_ROLE_APPLIED ||
                    (leaf.contributors[j].role == DUCKVEP_ROLE_POST_STOP && leaf.contributors[j].nmd_moved_junction);
                if (attributed != (!post[k] || moved)) goto done;
                nm_seen[NM_SEEN_ATTRIBUTED_POST] += post[k] && moved;
            }
        }
        if (!found) goto done;
    }
    nm_seen[x.nmd == DUCKVEP_HAPLOTYPE_NMD_TRIGGER ? NM_SEEN_TRIGGER : x.nmd == DUCKVEP_HAPLOTYPE_NMD_ESCAPE ? NM_SEEN_ESCAPE :
            x.nmd == DUCKVEP_HAPLOTYPE_NMD_UNKNOWN ? NM_SEEN_UNKNOWN : NM_SEEN_NA]++;
    nm_seen[NM_SEEN_INTRONLESS] += c->exons == 1u && x.stop_valid;
    nm_seen[NM_SEEN_REVERSE] += c->reverse; nm_seen[NM_SEEN_UTR] += c->u5 && c->u3; nm_seen[NM_SEEN_EXONS4] += c->exons == 4u && x.junction_valid;
    for (unsigned k = 0u; k < c->edit_count; k++) {
        nm_seen[NM_SEEN_POST] += post[k];
        if (x.junction_valid && c->exons > 1u && nm_edit_start(c, &c->edits[k]) <= c->ends[c->exons - 2u] && nm_edit_change(&c->edits[k]))
            nm_seen[NM_SEEN_INDEL_SHIFT]++;
    }
    result = THEFT_TRIAL_PASS;
done:
    free(f);
    return result;
}

TEST haplotype_nmd_matches_independent_edited_transcript_geometry(void) {
    /* Random layouts (either strand, UTR or none, 1 to 4 exons, UTR-only exons) and 1 to 3 edits (substitutions biased to
     * stops, and indels of 1 to 6 bases) recompute S and J from the edited sequence and exon geometry by an independent
     * oracle: the same prediction, S, J, the stop_gained agreement and the applied/post_stop role of every edit. */
    memset(nm_seen, 0, sizeof(nm_seen));
    struct theft_type_info type = {.alloc = nm_case_alloc, .free = nm_case_free};
    struct theft_run_config cfg = {0};
    cfg.name = "ejc50 == edited-transcript oracle";
    cfg.prop1 = prop_nmd_matches_edited_geometry_oracle;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    fprintf(stderr, "[nmd oracle coverage] attributed_post_stop=%u trigger=%u escape=%u unknown=%u not_applicable=%u intronless=%u reverse=%u utr=%u four_exon=%u "
        "indel_before_junction=%u post_stop_edits=%u\n", nm_seen[NM_SEEN_ATTRIBUTED_POST], nm_seen[NM_SEEN_TRIGGER], nm_seen[NM_SEEN_ESCAPE], nm_seen[NM_SEEN_UNKNOWN],
        nm_seen[NM_SEEN_NA], nm_seen[NM_SEEN_INTRONLESS], nm_seen[NM_SEEN_REVERSE], nm_seen[NM_SEEN_UTR], nm_seen[NM_SEEN_EXONS4],
        nm_seen[NM_SEEN_INDEL_SHIFT], nm_seen[NM_SEEN_POST]);
    if (cfg.trials >= KPROP_DEFAULT_TRIALS)
        ASSERT(nm_seen[NM_SEEN_TRIGGER] > 20u && nm_seen[NM_SEEN_ESCAPE] > 20u && nm_seen[NM_SEEN_UNKNOWN] > 5u &&
               nm_seen[NM_SEEN_NA] > 5u && nm_seen[NM_SEEN_INTRONLESS] > 5u && nm_seen[NM_SEEN_INDEL_SHIFT] > 20u &&
               nm_seen[NM_SEEN_UTR] > 20u && nm_seen[NM_SEEN_EXONS4] > 5u && nm_seen[NM_SEEN_POST] > 10u &&
               nm_seen[NM_SEEN_REVERSE] > 100u && nm_seen[NM_SEEN_ATTRIBUTED_POST] > 10u);
    PASS();
}

/* 240-base CDS: ATG AGC (CAA)x77 TAA, so a stop codon at CDS 49-51 (S = u5 + 51) leaves room for any junction distance. */
static void nm_caa_cds(char *cds) {
    memcpy(cds, "ATGAGC", 6u);
    for (unsigned i = 0u; i < 77u; i++) memcpy(cds + 6u + 3u * i, "CAA", 3u);
    memcpy(cds + 237u, "TAA", 3u);
    cds[240] = 0;
}

TEST haplotype_nmd_threshold_boundary_is_exactly_fifty(void) {
    /* J - S from -5 to 80 for a stop codon at CDS 17, on both strands, with and without a UTR and with and without the
     * model's cDNA origin cache: 50 escapes, 51 triggers, and S and J are the reported edited coordinates. A single
     * exon escapes at every distance and reports no junction. */
    unsigned checked = 0u;
    for (unsigned reverse = 0u; reverse < 2u; reverse++) for (unsigned cache = 0u; cache < 2u; cache++)
        for (unsigned utr = 0u; utr < 2u; utr++) for (int d = -5; d <= 80; d++) for (unsigned exons = 1u; exons <= 2u; exons++) {
            struct nm_case *c = calloc(1u, sizeof(*c));
            struct nmd_scene *f = malloc(sizeof(*f));
            ASSERT(c && f);
            c->reverse = reverse; c->cache = cache; c->u5 = utr ? 9u : 0u; c->u3 = utr ? 5u : 0u; c->cds_len = 240u; c->exons = exons;
            nm_caa_cds(c->cds);
            unsigned stop = c->u5 + 51u;
            c->ends[exons - 1u] = c->u5 + 240u + c->u3;
            if (exons == 2u) c->ends[0] = (unsigned)((int)stop + d);
            c->edits[0] = (struct nm_edit){NM_SUB, 49u, 1u, "T"};
            c->edit_count = 1u;
            duckvep_haplotype_leaf_t leaf;
            ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, nm_run(c, f, &leaf));
            ASSERT_EQ((int)DUCKVEP_PREDICTION_PREDICTED, (int)leaf.path_status);
            ASSERT(leaf.haplotype_so_mask & DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED));
            ASSERT_EQ(1u, leaf.nmd_stop_valid);
            ASSERT_EQ(stop, leaf.nmd_stop_position1);
            if (exons == 1u) {
                ASSERT_EQ((int)DUCKVEP_HAPLOTYPE_NMD_ESCAPE, (int)leaf.nmd);
                ASSERT_EQ(0u, leaf.nmd_junction_valid);
            } else {
                ASSERT_EQ(1u, leaf.nmd_junction_valid);
                ASSERT_EQ((unsigned)((int)stop + d), leaf.nmd_junction_position1);
                ASSERT_EQ((int)(d > 50 ? DUCKVEP_HAPLOTYPE_NMD_TRIGGER : DUCKVEP_HAPLOTYPE_NMD_ESCAPE), (int)leaf.nmd);
            }
            free(c); free(f);
            checked++;
        }
    ASSERT_EQ(2u * 2u * 2u * 86u * 2u, checked);
    PASS();
}

struct nm_post_case { struct nm_case base; struct nm_edit extra; unsigned valid; };

static enum theft_alloc_res nm_post_alloc(struct theft *t, void *env, void **instance) {
    (void)env;
    struct nm_post_case *p = calloc(1u, sizeof(*p));
    if (!p) return THEFT_ALLOC_ERROR;
    struct nm_case *c = &p->base;
    unsigned codons = 40u + (unsigned)kprop_bounded(t, 30u);
    c->reverse = (unsigned)kprop_bounded(t, 2u); c->cache = (unsigned)kprop_bounded(t, 2u);
    c->cds_len = 3u * codons;
    c->u5 = kprop_bounded(t, 2u) ? (unsigned)kprop_bounded(t, 16u) : 0u; c->u3 = kprop_bounded(t, 2u) ? (unsigned)kprop_bounded(t, 16u) : 0u;
    c->exons = 2u + (unsigned)kprop_bounded(t, 2u);
    nm_random_cds(t, codons, c->cds);
    /* One stop-creating substitution early enough that the junction can sit 47 to 53 bases after it. */
    struct nm_edit stop = {0};
    for (unsigned tries = 0u; tries < 40u && !stop.n; tries++) {
        unsigned codon = 2u + (unsigned)kprop_bounded(t, 12u), at = 3u * codon;
        static const char *const stops[3] = {"TAA", "TAG", "TGA"};
        for (unsigned s = 0u; s < 3u && !stop.n; s++) {
            unsigned diff = 0u, where = 0u;
            for (unsigned k = 0u; k < 3u; k++) if (c->cds[at + k] != stops[s][k]) { diff++; where = k; }
            if (diff == 1u) stop = (struct nm_edit){NM_SUB, at + where + 1u, 1u, {stops[s][where]}};
        }
    }
    if (!stop.n) { *instance = p; return THEFT_ALLOC_OK; }
    c->edits[0] = stop; c->edit_count = 1u;
    unsigned s_cds = 3u * ((stop.p - 1u) / 3u + 1u), s_mrna = c->u5 + s_cds; /* S as an mRNA coordinate */
    unsigned mrna = c->u5 + c->cds_len + c->u3, junction = s_mrna + 47u + (unsigned)kprop_bounded(t, 7u);
    if (junction + 30u > mrna) { *instance = p; return THEFT_ALLOC_OK; }
    if (c->exons == 2u) c->ends[0] = junction;
    else {
        unsigned cut = 10u + (unsigned)kprop_bounded(t, junction > 30u ? junction - 30u : 1u);
        if (cut + 8u > junction) { *instance = p; return THEFT_ALLOC_OK; }
        c->ends[0] = cut; c->ends[1] = junction;
    }
    c->ends[c->exons - 1u] = mrna;
    /* The post-stop edit: after the stop with margin, before the last 12 CDS bases, anywhere in the penultimate or last exon. */
    unsigned from = s_cds + 6u, limit = c->cds_len - 12u;
    if (from + 8u > limit) { *instance = p; return THEFT_ALLOC_OK; }
    struct nm_edit e = {0};
    unsigned kind = (unsigned)kprop_bounded(t, 3u);
    e.p = from + (unsigned)kprop_bounded(t, limit - from);
    if (kind == 0u) { e.kind = NM_SUB; e.n = 1u; e.bases[0] = c->cds[e.p - 1u] == 'A' ? 'G' : 'A'; }
    else if (kind == 1u) { e.kind = NM_INS; e.n = 1u + (unsigned)kprop_bounded(t, 3u); for (unsigned i = 0u; i < e.n; i++) e.bases[i] = "ACGT"[kprop_bounded(t, 4u)]; }
    else { e.kind = NM_DEL; e.n = 1u + (unsigned)kprop_bounded(t, 3u); }
    if (e.p + e.n + 3u > limit || !nm_unambiguous(c, &e) || !nm_inside_one_exon(c, &e)) { *instance = p; return THEFT_ALLOC_OK; }
    p->extra = e; p->valid = 1u;
    *instance = p;
    return THEFT_ALLOC_OK;
}

static enum theft_trial_res prop_nmd_ignores_post_stop_edits_except_through_junction(struct theft *t, void *arg) {
    (void)t;
    const struct nm_post_case *p = arg;
    if (!p->valid) return THEFT_TRIAL_SKIP;
    struct nm_case with = p->base;
    with.edits[1] = p->extra; with.edit_count = 2u;
    struct nmd_scene *f = malloc(sizeof(*f)), *g = malloc(sizeof(*g));
    if (!f || !g) { free(f); free(g); return THEFT_TRIAL_ERROR; }
    duckvep_haplotype_leaf_t before, after;
    enum theft_trial_res result = THEFT_TRIAL_FAIL;
    if (nm_run(&p->base, f, &before) != DUCKVEP_HAPLOTYPE_STREAM_OK || nm_run(&with, g, &after) != DUCKVEP_HAPLOTYPE_STREAM_OK ||
        before.path_status != DUCKVEP_PREDICTION_PREDICTED || after.path_status != DUCKVEP_PREDICTION_PREDICTED ||
        !before.nmd_stop_valid || !after.nmd_stop_valid || !before.nmd_junction_valid || !after.nmd_junction_valid) goto done;
    /* The stop, its attribution and the extra edit's role never depend on an edit after the stop. */
    if (before.nmd_stop_position1 != after.nmd_stop_position1 || after.contributor_count != 2u) goto done;
    for (size_t j = 0u; j < after.contributor_count; j++) {
        uint64_t id = after.contributors[j].source.event_id;
        if (id == 1u && after.contributors[j].role != DUCKVEP_ROLE_APPLIED) goto done;
        if (id == 2u && after.contributors[j].role != DUCKVEP_ROLE_POST_STOP) goto done;
        /* Its attribution follows J: attributed exactly when it is an indel at or before the penultimate exon's last base. */
        if (id == 2u && (after.contributors[j].nmd_moved_junction != 0u) !=
            (nm_edit_change(&p->extra) != 0 && nm_edit_start(&with, &p->extra) <= p->base.ends[p->base.exons - 2u])) goto done;
        if (id == 1u && after.contributors[j].nmd_moved_junction) goto done;
    }
    /* J moves by exactly the extra edit's length change when it starts at or before the penultimate exon's last base. */
    bool moves = nm_edit_start(&with, &p->extra) <= p->base.ends[p->base.exons - 2u];
    int64_t expect = (int64_t)before.nmd_junction_position1 + (moves ? nm_edit_change(&p->extra) : 0);
    if ((int64_t)after.nmd_junction_position1 != expect) goto done;
    int64_t distance = (int64_t)after.nmd_junction_position1 - (int64_t)after.nmd_stop_position1;
    if ((int)after.nmd != (int)(distance > 50 ? DUCKVEP_HAPLOTYPE_NMD_TRIGGER : DUCKVEP_HAPLOTYPE_NMD_ESCAPE)) goto done;
    /* An edit that leaves the distance alone leaves the whole prediction alone. */
    if (!moves || !nm_edit_change(&p->extra)) {
        if (after.nmd != before.nmd || after.nmd_junction_position1 != before.nmd_junction_position1) goto done;
    }
    nm_seen[NM_SEEN_POST]++;
    if (moves && nm_edit_change(&p->extra)) { nm_seen[NM_SEEN_INDEL_SHIFT]++; nm_seen[NM_SEEN_FLIP] += after.nmd != before.nmd; }
    else nm_seen[NM_SEEN_RUN]++;
    result = THEFT_TRIAL_PASS;
done:
    free(f); free(g);
    return result;
}

TEST haplotype_nmd_is_invariant_to_post_stop_edits_that_leave_the_junction_distance(void) {
    /* A stop-creating substitution with the junction 47 to 53 bases after it, then one more edit after the stop codon.
     * S, the applied set and the prediction never depend on it, except through the edited junction: a substitution or
     * an indel in the last exon leaves the prediction and J untouched, while an indel up to the penultimate exon's last
     * base moves J by exactly its length change (and can flip the prediction across the threshold). */
    memset(nm_seen, 0, sizeof(nm_seen));
    struct theft_type_info type = {.alloc = nm_post_alloc, .free = nm_case_free};
    struct theft_run_config cfg = {0};
    cfg.name = "post-stop edits move J only";
    cfg.prop1 = prop_nmd_ignores_post_stop_edits_except_through_junction;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    fprintf(stderr, "[nmd post-stop coverage] checked=%u junction_shifts=%u flips=%u unchanged=%u\n", nm_seen[NM_SEEN_POST],
        nm_seen[NM_SEEN_INDEL_SHIFT], nm_seen[NM_SEEN_FLIP], nm_seen[NM_SEEN_RUN]);
    if (cfg.trials >= KPROP_DEFAULT_TRIALS)
        ASSERT(nm_seen[NM_SEEN_INDEL_SHIFT] > 30u && nm_seen[NM_SEEN_FLIP] > 5u && nm_seen[NM_SEEN_RUN] > 30u);
    PASS();
}

TEST haplotype_nmd_no_op_alt_is_rejected_so_no_reference_lane_leaf_exists(void) {
    /* Reference-only lanes are not_applicable (a path with no edit has no termination change; the kernel labels such a
     * leaf so), but no such leaf can be built from records: an ALT equal to REF is rejected before any leaf exists, and
     * an ALT-free lane of a strict decoded call has no carrier. The one-substitution control is decided by the rule. */
    struct nm_case *c = calloc(1u, sizeof(*c));
    struct nmd_scene *f = malloc(sizeof(*f));
    ASSERT(c && f);
    c->cds_len = 240u; c->exons = 2u; c->u5 = 3u; c->u3 = 4u; c->ends[0] = 90u; c->ends[1] = 247u;
    nm_caa_cds(c->cds);
    c->edits[0] = (struct nm_edit){NM_SUB, 49u, 1u, "T"}; c->edit_count = 1u;
    duckvep_haplotype_leaf_t leaf;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, nm_run(c, f, &leaf));
    ASSERT_EQ((int)DUCKVEP_HAPLOTYPE_NMD_ESCAPE, (int)leaf.nmd); /* J - S = 90 - 54 = 36 */
    ASSERT_EQ(54u, (unsigned)leaf.nmd_stop_position1);
    ASSERT_EQ(90u, (unsigned)leaf.nmd_junction_position1);
    nm_scene_prepare(f, c);
    duckvep_haplotype_stream_t *s = &f->stream;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_OK, duckvep_haplotype_stream_init(s, &f->model, &f->exons, &f->sequences, &f->buffers));
    struct nm_record r;
    nm_record_of(c, f, &c->edits[0], &r);
    char same_alt[2] = {r.ref[0], 0};
    duckvep_haplotype_source_t source = {1u, (const uint8_t *)r.ref, (const uint8_t *)same_alt, r.pos, 0u, 1u, 1u, 0u, 0u, 0u, 1u};
    uint32_t tx = 0u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG, haplotype_test_begin_candidates(s, &source, &tx, 1u));
    free(c); free(f);
    PASS();
}
