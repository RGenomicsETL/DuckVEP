/* Bounded hypothetical diploid-arrangement replay shared by both hosts. */
#include "core/duckvep_core_arrangements.h"

#include "kernel/src/duckvep_budget.h"
#include "kernel/src/duckvep_phase.h"

#include <stdio.h>
#include <string.h>

/* Every staged REF/ALT pair is bounded independently, and prediction copies are
 * bounded before input is read. These limits keep enumeration storage auditable. */
enum { ARRANGEMENT_CALL_BYTES = 65536u, ARRANGEMENT_SEQUENCE_BYTES = 65536u };

typedef struct {
    uint64_t event_id, position;
    uint32_t chrom, allele_index, transcript, sample, site;
    uint32_t ref_offset, alt_offset, ref_length, alt_length;
    int32_t allele[2];
    int64_t phase_set;
    uint8_t phase[2], phase_present, phase_set_present;
} arrangement_call_t;

typedef struct {
    uint32_t transcript;
    size_t cds_length, protein_length;
    int64_t nominal_length_diff;
    uint64_t consequence_mask;
    uint8_t prediction_status, prediction_reason, hypothetical, valid;
} arrangement_prediction_t;

const char *const duckvep_arrangement_limit_names[DUCKVEP_ARRANGEMENT_LIMIT_COUNT] = {
    "max_sites", "max_calls", "max_arrangements", "max_replays"
};

struct duckvep_arrangement_state {
    duckvep_arrangement_config_t config;
    duckvep_hap_state_t *replayer;
    arrangement_call_t *calls;
    duckvep_phase_arrangement_site_t *sites;
    uint16_t *lanes;
    uint32_t *transcripts;
    arrangement_prediction_t *predictions;
    uint8_t *call_bytes, *cds_bytes, *protein_bytes;
    size_t call_byte_used, call_byte_capacity, call_count, site_count, transcript_count;
    size_t arrangement_count, prediction_count, row_count, output_index, workspace_bytes;
    uint32_t sample;
    uint8_t loaded;
};

static int mul_size(size_t left, size_t right, size_t *out) {
    if (left && right > SIZE_MAX / left) return 0;
    *out = left * right;
    return 1;
}

static int add_size(size_t left, size_t right, size_t *out) {
    if (right > SIZE_MAX - left) return 0;
    *out = left + right;
    return 1;
}

static int bit_valid(const uint64_t *bits, size_t bit) {
    return !bits || ((bits[bit >> 6u] >> (bit & 63u)) & UINT64_C(1));
}

static void arrangement_error(char *error, size_t error_size, const char *message) {
    if (error_size) snprintf(error, error_size, "%s", message);
}

void duckvep_arrangement_config_defaults(duckvep_arrangement_config_t *config,
    const duckvep_owned_model_t *model) {
    memset(config, 0, sizeof(*config));
    config->model = model;
    config->max_sites = DUCKVEP_ARRANGEMENT_DEFAULT_SITES;
    config->max_calls = DUCKVEP_ARRANGEMENT_DEFAULT_CALLS;
    config->max_arrangements = DUCKVEP_ARRANGEMENT_DEFAULT_ALTERNATIVES;
    config->max_replays = DUCKVEP_ARRANGEMENT_DEFAULT_REPLAYS;
}

void duckvep_arrangement_close(duckvep_arrangement_state_t *state) {
    if (!state) return;
    duckvep_hap_close(state->replayer);
    duckvep_budget_free(state->protein_bytes);
    duckvep_budget_free(state->cds_bytes);
    duckvep_budget_free(state->predictions);
    duckvep_budget_free(state->transcripts);
    duckvep_budget_free(state->lanes);
    duckvep_budget_free(state->sites);
    duckvep_budget_free(state->call_bytes);
    duckvep_budget_free(state->calls);
    duckvep_budget_free(state);
}

static int allocate_array(void **out, size_t count, size_t size) {
    if (!count || count > SIZE_MAX / size) return 0;
    *out = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, count, size);
    return *out != NULL;
}

int duckvep_arrangement_config_check(const duckvep_arrangement_config_t *config,
    char *error, size_t error_size) {
    enum { WORKSPACE_LIMIT = 134217728u };
    size_t lane_count, call_bytes, carrier_count, sequence_bytes, calls_bytes, sites_bytes,
        lanes_bytes, transcripts_bytes, predictions_bytes, outer_bytes;
    if (error_size) error[0] = '\0';
    if (!config || !config->model || !config->max_sites || !config->max_calls ||
        !config->max_arrangements || !config->max_replays ||
        config->max_sites > DUCKVEP_PHASE_ARRANGEMENT_MAX_SITES ||
        config->max_calls > UINT32_MAX / 2u || config->max_arrangements > (UINT32_MAX - 1u) / 2u) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: invalid bounded limits");
        return 0;
    }
    if (!mul_size(config->max_arrangements, config->max_sites, &lane_count) ||
        !mul_size(config->max_calls, 2u, &carrier_count) ||
        !mul_size(config->max_calls, ARRANGEMENT_CALL_BYTES, &call_bytes) ||
        !mul_size(config->max_replays, ARRANGEMENT_SEQUENCE_BYTES, &sequence_bytes) ||
        !mul_size(config->max_calls, sizeof(arrangement_call_t), &calls_bytes) ||
        !mul_size(config->max_sites, sizeof(duckvep_phase_arrangement_site_t), &sites_bytes) ||
        !mul_size(lane_count, sizeof(uint16_t), &lanes_bytes) ||
        !mul_size(config->max_calls, sizeof(uint32_t), &transcripts_bytes) ||
        !mul_size(config->max_replays, sizeof(arrangement_prediction_t), &predictions_bytes) ||
        !add_size(sizeof(duckvep_arrangement_state_t), calls_bytes, &outer_bytes) ||
        !add_size(outer_bytes, call_bytes, &outer_bytes) || !add_size(outer_bytes, sites_bytes, &outer_bytes) ||
        !add_size(outer_bytes, lanes_bytes, &outer_bytes) || !add_size(outer_bytes, transcripts_bytes, &outer_bytes) ||
        !add_size(outer_bytes, predictions_bytes, &outer_bytes) || !add_size(outer_bytes, sequence_bytes, &outer_bytes) ||
        !add_size(outer_bytes, sequence_bytes, &outer_bytes) || outer_bytes >= WORKSPACE_LIMIT) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: bounded workspace overflows size_t or capacity");
        return 0;
    }
    return 1;
}

duckvep_arrangement_state_t *duckvep_arrangement_open(const duckvep_arrangement_config_t *config,
    char *error, size_t error_size) {
    enum { WORKSPACE_LIMIT = 134217728u };
    duckvep_arrangement_state_t *state;
    duckvep_hap_config_t replay;
    size_t lane_count, call_bytes, carrier_count, sequence_bytes, calls_bytes, sites_bytes,
        lanes_bytes, transcripts_bytes, predictions_bytes, outer_bytes, replay_bytes;
    if (!duckvep_arrangement_config_check(config, error, error_size)) return NULL;
    /* Compute every charged allocation, including both lane sequence stores and
     * the reusable native replay workspace, before allocating any of them. */
    if (!mul_size(config->max_arrangements, config->max_sites, &lane_count) ||
        !mul_size(config->max_calls, 2u, &carrier_count) ||
        !mul_size(config->max_calls, ARRANGEMENT_CALL_BYTES, &call_bytes) ||
        !mul_size(config->max_replays, ARRANGEMENT_SEQUENCE_BYTES, &sequence_bytes) ||
        !mul_size(config->max_calls, sizeof(*state->calls), &calls_bytes) ||
        !mul_size(config->max_sites, sizeof(*state->sites), &sites_bytes) ||
        !mul_size(lane_count, sizeof(*state->lanes), &lanes_bytes) ||
        !mul_size(config->max_calls, sizeof(*state->transcripts), &transcripts_bytes) ||
        !mul_size(config->max_replays, sizeof(*state->predictions), &predictions_bytes) ||
        !add_size(sizeof(*state), calls_bytes, &outer_bytes) ||
        !add_size(outer_bytes, call_bytes, &outer_bytes) ||
        !add_size(outer_bytes, sites_bytes, &outer_bytes) ||
        !add_size(outer_bytes, lanes_bytes, &outer_bytes) ||
        !add_size(outer_bytes, transcripts_bytes, &outer_bytes) ||
        !add_size(outer_bytes, predictions_bytes, &outer_bytes) ||
        !add_size(outer_bytes, sequence_bytes, &outer_bytes) || !add_size(outer_bytes, sequence_bytes, &outer_bytes) ||
        outer_bytes >= WORKSPACE_LIMIT) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: bounded workspace overflows size_t or capacity");
        return NULL;
    }
    duckvep_hap_config_defaults(&replay, config->model);
    replay.policy = DUCKVEP_PHASE_STRICT;
    replay.hgvs = 0;
    replay.limits[DUCKVEP_HAP_LIMIT_EVENTS] = config->max_calls;
    replay.limits[DUCKVEP_HAP_LIMIT_TRANSCRIPTS] = config->max_calls;
    replay.limits[DUCKVEP_HAP_LIMIT_CARRIERS] = carrier_count;
    replay.limits[DUCKVEP_HAP_LIMIT_PREFIXES] = carrier_count;
    replay.limits[DUCKVEP_HAP_LIMIT_PROJECTIONS] = carrier_count;
    replay.limits[DUCKVEP_HAP_LIMIT_ALLELES] = call_bytes;
    replay.limits[DUCKVEP_HAP_LIMIT_LEAF_EVENTS] = config->max_calls;
    replay.limits[DUCKVEP_HAP_LIMIT_LEAF_EDITS] = config->max_calls;
    replay.limits[DUCKVEP_HAP_LIMIT_SEQUENCE] = ARRANGEMENT_SEQUENCE_BYTES;
    replay.limits[DUCKVEP_HAP_LIMIT_PLOIDY] = 2u;
    replay.limits[DUCKVEP_HAP_LIMIT_PHASE_SETS] = 1u;
    replay.limits[DUCKVEP_HAP_LIMIT_DIFFERENCES] = config->max_calls;
    replay.limits[DUCKVEP_HAP_LIMIT_HGVS_OPERATIONS] = config->max_calls;
    replay.limits[DUCKVEP_HAP_LIMIT_WORKSPACE] = WORKSPACE_LIMIT - outer_bytes;
    if (!duckvep_hap_workspace_estimate(&replay, &replay_bytes) ||
        !add_size(outer_bytes, replay_bytes, &outer_bytes) || outer_bytes > WORKSPACE_LIMIT) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: replay workspace exceeds capacity");
        return NULL;
    }
    state = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1u, sizeof(*state));
    if (!state) goto exhausted;
    state->config = *config;
    state->call_byte_capacity = call_bytes;
    if (!allocate_array((void **)&state->calls, config->max_calls, sizeof(*state->calls)) ||
        !allocate_array((void **)&state->call_bytes, call_bytes, sizeof(*state->call_bytes)) ||
        !allocate_array((void **)&state->sites, config->max_sites, sizeof(*state->sites)) ||
        !allocate_array((void **)&state->lanes, lane_count, sizeof(*state->lanes)) ||
        !allocate_array((void **)&state->transcripts, config->max_calls, sizeof(*state->transcripts)) ||
        !allocate_array((void **)&state->predictions, config->max_replays, sizeof(*state->predictions)) ||
        !allocate_array((void **)&state->cds_bytes, config->max_replays, ARRANGEMENT_SEQUENCE_BYTES) ||
        !allocate_array((void **)&state->protein_bytes, config->max_replays, ARRANGEMENT_SEQUENCE_BYTES)) goto exhausted;
    state->replayer = duckvep_hap_open(&replay, error, error_size);
    if (!state->replayer) goto failed;
    state->workspace_bytes = outer_bytes;
    return state;
exhausted:
    arrangement_error(error, error_size, "duckvep_haplotype_arrangements: bounded workspace allocation failed");
failed:
    duckvep_arrangement_close(state);
    return NULL;
}

static int find_transcript(duckvep_arrangement_state_t *state, uint32_t transcript, size_t *index) {
    for (size_t i = 0u; i < state->transcript_count; i++) {
        if (state->transcripts[i] == transcript) { *index = i; return 1; }
    }
    if (state->transcript_count == state->config.max_calls) return 0;
    *index = state->transcript_count;
    state->transcripts[state->transcript_count++] = transcript;
    return 1;
}

static int find_site(const duckvep_arrangement_state_t *state, uint64_t event_id, size_t *index) {
    for (size_t i = 0u; i < state->site_count; i++) {
        if (state->sites[i].source_id == event_id) { *index = i; return 1; }
    }
    return 0;
}

static int stage_row(duckvep_arrangement_state_t *state, const duckvep_hap_row_t *row,
    char *error, size_t error_size) {
    arrangement_call_t *call;
    duckvep_phase_arrangement_site_t *site;
    size_t site_index, transcript_index;
    if (state->call_count == state->config.max_calls) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: max_calls exhausted before enumeration");
        return 0;
    }
    if ((row->null_mask & UINT32_C(0x000071ff)) || !row->event_id || row->allele_index != 1u || row->copies != 1 ||
        row->versions != 1 || row->ploidies != 1 || row->gt_length != 2u ||
        !bit_valid(row->gt_validity, row->gt_offset) || !bit_valid(row->gt_validity, row->gt_offset + 1u) ||
        !row->ref || !row->alt || !row->ref_len || !row->alt_len ||
        row->ref_len > ARRANGEMENT_CALL_BYTES || row->alt_len > ARRANGEMENT_CALL_BYTES ||
        row->gt[row->gt_offset] < 0 || row->gt[row->gt_offset] > 1 ||
        row->gt[row->gt_offset + 1u] < 0 || row->gt[row->gt_offset + 1u] > 1 ||
        row->gt[row->gt_offset] == row->gt[row->gt_offset + 1u]) {
        snprintf(error, error_size,
            "duckvep_haplotype_arrangements: requires non-NULL diploid biallelic heterozygous alt_events rows (null=%u,event=%llu,alt=%u,copies=%lld,versions=%lld,ploidies=%lld,gt=%zu)",
            row->null_mask, (unsigned long long)row->event_id, row->allele_index, (long long)row->copies,
            (long long)row->versions, (long long)row->ploidies, row->gt_length);
        return 0;
    }
    if (row->have_phase && (row->phase_length != 2u || !bit_valid(row->phase_validity, row->phase_offset) ||
        !bit_valid(row->phase_validity, row->phase_offset + 1u))) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: phase_before must be a complete diploid call");
        return 0;
    }
    if (row->ref_len > state->call_byte_capacity - state->call_byte_used ||
        row->alt_len > state->call_byte_capacity - state->call_byte_used - row->ref_len) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: bounded call storage exhausted before enumeration");
        return 0;
    }
    if (state->call_count && row->sample != state->sample) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: input contains more than one sample");
        return 0;
    }
    if (!find_transcript(state, row->transcript, &transcript_index)) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: max_calls exhausted by transcript identities");
        return 0;
    }
    call = &state->calls[state->call_count];
    memset(call, 0, sizeof(*call));
    call->event_id = row->event_id; call->position = row->pos; call->chrom = row->chrom;
    call->allele_index = row->allele_index; call->transcript = row->transcript; call->sample = row->sample;
    call->ref_offset = (uint32_t)state->call_byte_used; call->ref_length = row->ref_len;
    memcpy(state->call_bytes + state->call_byte_used, row->ref, row->ref_len);
    state->call_byte_used += row->ref_len;
    call->alt_offset = (uint32_t)state->call_byte_used; call->alt_length = row->alt_len;
    memcpy(state->call_bytes + state->call_byte_used, row->alt, row->alt_len);
    state->call_byte_used += row->alt_len;
    call->allele[0] = row->gt[row->gt_offset]; call->allele[1] = row->gt[row->gt_offset + 1u];
    call->phase_present = (uint8_t)(row->have_phase != 0);
    if (call->phase_present) {
        call->phase[0] = row->phase[row->phase_offset] != 0;
        call->phase[1] = row->phase[row->phase_offset + 1u] != 0;
    }
    call->phase_set_present = (uint8_t)(row->phase_set_present != 0);
    call->phase_set = row->phase_set;
    if (!find_site(state, row->event_id, &site_index)) {
        if (state->site_count == state->config.max_sites) {
            arrangement_error(error, error_size, "duckvep_haplotype_arrangements: max_sites exhausted before enumeration");
            return 0;
        }
        site_index = state->site_count++;
        site = &state->sites[site_index];
        memset(site, 0, sizeof(*site));
        site->source_id = row->event_id; site->allele[0] = call->allele[0]; site->allele[1] = call->allele[1];
        site->sample_index = row->sample; site->source_alt_count = 1u;
        site->phase_before[0] = call->phase[0]; site->phase_before[1] = call->phase[1];
        site->phase_set_present = call->phase_set_present; site->phase_set = call->phase_set;
    } else {
        site = &state->sites[site_index];
        if (site->sample_index != row->sample || site->allele[0] != call->allele[0] ||
            site->allele[1] != call->allele[1] || site->phase_before[0] != call->phase[0] ||
            site->phase_before[1] != call->phase[1] || site->phase_set_present != call->phase_set_present ||
            (site->phase_set_present && site->phase_set != call->phase_set)) {
            arrangement_error(error, error_size,
                "duckvep_haplotype_arrangements: source call differs between transcript projections");
            return 0;
        }
    }
    call->site = (uint32_t)site_index;
    state->sample = row->sample;
    state->call_count++;
    (void)transcript_index;
    return 1;
}

typedef struct {
    duckvep_arrangement_state_t *state;
    size_t hypothesis, transcript_index, lane, next_call;
    duckvep_hap_replay_context_t context;
    duckvep_hap_row_t row;
} replay_input_t;

static int replay_next(void *context, duckvep_hap_row_t *row, char *error, size_t error_size) {
    replay_input_t *input = context;
    duckvep_arrangement_state_t *state = input->state;
    const arrangement_call_t *call;
    while (input->next_call < state->call_count) {
        call = &state->calls[input->next_call++];
        if (call->transcript != state->transcripts[input->transcript_index] ||
            state->lanes[input->hypothesis * state->site_count + call->site] != input->lane + 1u) continue;
        memset(&input->row, 0, sizeof(input->row));
        input->row.event_id = call->event_id; input->row.pos = call->position; input->row.chrom = call->chrom;
        input->row.allele_index = call->allele_index; input->row.transcript = call->transcript; input->row.sample = call->sample;
        input->row.ref = state->call_bytes + call->ref_offset; input->row.ref_len = call->ref_length;
        input->row.alt = state->call_bytes + call->alt_offset; input->row.alt_len = call->alt_length;
        input->row.copies = input->row.versions = input->row.ploidies = 1u;
        input->row.replay = &input->context;
        *row = input->row;
        return 1;
    }
    (void)error; (void)error_size;
    return 0;
}

typedef struct { duckvep_arrangement_state_t *state; size_t index, count; } prediction_sink_t;

static int prediction_sink(void *context, const duckvep_hap_prediction_t *prediction,
    char *error, size_t error_size) {
    prediction_sink_t *sink = context;
    duckvep_arrangement_state_t *state = sink->state;
    arrangement_prediction_t *stored;
    if (sink->count++ || sink->index >= state->config.max_replays ||
        prediction->cds_length > ARRANGEMENT_SEQUENCE_BYTES || prediction->protein_length > ARRANGEMENT_SEQUENCE_BYTES) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: replay prediction exceeds bounded storage");
        return 0;
    }
    stored = &state->predictions[sink->index];
    stored->transcript = prediction->transcript; stored->cds_length = prediction->cds_length;
    stored->protein_length = prediction->protein_length; stored->nominal_length_diff = prediction->nominal_length_diff;
    stored->consequence_mask = prediction->consequence_mask;
    stored->prediction_status = (uint8_t)prediction->prediction_status;
    stored->prediction_reason = (uint8_t)prediction->prediction_reason;
    stored->hypothetical = prediction->hypothetical;
    stored->valid = 1u;
    if (prediction->cds_length) memcpy(state->cds_bytes + sink->index * ARRANGEMENT_SEQUENCE_BYTES,
        prediction->cds, prediction->cds_length);
    if (prediction->protein_length) memcpy(state->protein_bytes + sink->index * ARRANGEMENT_SEQUENCE_BYTES,
        prediction->protein, prediction->protein_length);
    return 1;
}

static int replay_hypotheses(duckvep_arrangement_state_t *state, char *error, size_t error_size) {
    for (size_t hypothesis = 0u; hypothesis < state->arrangement_count; hypothesis++) {
        for (size_t transcript = 0u; transcript < state->transcript_count; transcript++) {
            for (size_t lane = 0u; lane < 2u; lane++) {
                replay_input_t input;
                prediction_sink_t sink;
                size_t prediction_index = (hypothesis * state->transcript_count + transcript) * 2u + lane;
                int has_contributor = 0;
                for (size_t call = 0u; call < state->call_count; call++) {
                    if (state->calls[call].transcript == state->transcripts[transcript] &&
                        state->lanes[hypothesis * state->site_count + state->calls[call].site] == lane + 1u) {
                        has_contributor = 1;
                        break;
                    }
                }
                sink.state = state; sink.index = prediction_index; sink.count = 0u;
                if (!has_contributor) {
                    duckvep_hap_prediction_t reference;
                    if (!duckvep_hap_reference(state->replayer, state->transcripts[transcript], &reference,
                            error, error_size) || !prediction_sink(&sink, &reference, error, error_size)) return 0;
                } else {
                    duckvep_hap_input_t replay = {&input, replay_next};
                    memset(&input, 0, sizeof(input));
                    input.state = state; input.hypothesis = hypothesis; input.transcript_index = transcript;
                    input.lane = lane;
                    input.context.kind = DUCKVEP_HAP_REPLAY_HYPOTHETICAL_ALT;
                    input.context.partition = (uint32_t)(hypothesis * 2u + lane + 1u);
                    input.context.lane = (uint16_t)(lane + 1u);
                    if (!duckvep_hap_replay(state->replayer, &replay, prediction_sink, &sink, error, error_size)) return 0;
                }
                if (sink.count != 1u) {
                    arrangement_error(error, error_size,
                        "duckvep_haplotype_arrangements: replay did not yield one lane prediction");
                    return 0;
                }
            }
        }
    }
    return 1;
}

int duckvep_arrangement_load(duckvep_arrangement_state_t *state, const duckvep_hap_input_t *input,
    char *error, size_t error_size) {
    duckvep_phase_arrangement_status_t phase_status;
    size_t alternatives = 0u, predictions, rows;
    if (error_size) error[0] = '\0';
    if (!state || !input || !input->next || state->loaded) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: invalid or reused arrangement state");
        return 0;
    }
    for (;;) {
        duckvep_hap_row_t row;
        int status = input->next(input->context, &row, error, error_size);
        if (status < 0) return 0;
        if (!status) break;
        if (!stage_row(state, &row, error, error_size)) return 0;
    }
    if (!state->site_count || !state->transcript_count) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: input has no eligible calls");
        return 0;
    }
    phase_status = duckvep_phase_arrange_diploid(state->sites, state->site_count, state->config.max_arrangements,
        NULL, 0u, &alternatives);
    if (phase_status != DUCKVEP_PHASE_ARRANGEMENT_OK) {
        snprintf(error, error_size, "duckvep_haplotype_arrangements: phase arrangement status %u before publication",
            (unsigned)phase_status);
        return 0;
    }
    if (!alternatives || alternatives > state->config.max_arrangements ||
        !mul_size(alternatives, state->transcript_count, &predictions) ||
        !mul_size(predictions, 2u, &predictions) || predictions > state->config.max_replays ||
        !mul_size(alternatives, state->call_count, &rows) || !mul_size(rows, 2u, &rows)) {
        arrangement_error(error, error_size, "duckvep_haplotype_arrangements: alternative/replay capacity exceeded before publication");
        return 0;
    }
    phase_status = duckvep_phase_arrange_diploid(state->sites, state->site_count, state->config.max_arrangements,
        state->lanes, state->config.max_arrangements * state->site_count, &alternatives);
    if (phase_status != DUCKVEP_PHASE_ARRANGEMENT_OK) {
        snprintf(error, error_size, "duckvep_haplotype_arrangements: phase arrangement status %u before publication",
            (unsigned)phase_status);
        return 0;
    }
    state->arrangement_count = alternatives;
    state->prediction_count = predictions;
    state->row_count = rows;
    if (!replay_hypotheses(state, error, error_size)) return 0;
    state->loaded = 1u;
    return 1;
}

int duckvep_arrangement_next(duckvep_arrangement_state_t *state, duckvep_arrangement_row_t *row) {
    while (state && state->output_index < state->row_count) {
        size_t index = state->output_index++;
        size_t hypothesis = index / (state->call_count * 2u);
        size_t remainder = index % (state->call_count * 2u);
        size_t lane = remainder / state->call_count;
        arrangement_call_t *call = &state->calls[remainder % state->call_count];
        size_t transcript_index;
        arrangement_prediction_t *prediction;
        uint16_t assigned;
        int lane_has_contributor = 0;
        for (transcript_index = 0u; transcript_index < state->transcript_count; transcript_index++) {
            if (state->transcripts[transcript_index] == call->transcript) break;
        }
        if (transcript_index == state->transcript_count) return 0;
        prediction = &state->predictions[(hypothesis * state->transcript_count + transcript_index) * 2u + lane];
        for (size_t source = 0u; source < state->call_count; source++) {
            if (state->calls[source].transcript == call->transcript &&
                state->lanes[hypothesis * state->site_count + state->calls[source].site] == lane + 1u) {
                lane_has_contributor = 1;
                break;
            }
        }
        assigned = state->lanes[hypothesis * state->site_count + call->site];
        memset(row, 0, sizeof(*row));
        row->hypothesis_id = hypothesis + 1u; row->lane = (uint16_t)(lane + 1u);
        row->event_id = call->event_id; row->transcript = call->transcript; row->sample = call->sample;
        row->chrom = call->chrom; row->position = call->position; row->allele_index = call->allele_index;
        row->reference = state->call_bytes + call->ref_offset; row->reference_length = call->ref_length;
        row->alternate = state->call_bytes + call->alt_offset; row->alternate_length = call->alt_length;
        row->allele0 = call->allele[0]; row->allele1 = call->allele[1];
        row->phase0 = call->phase[0]; row->phase1 = call->phase[1]; row->phase_present = call->phase_present;
        row->phase_set_present = call->phase_set_present; row->phase_set = call->phase_set;
        row->assigned_allele = assigned == lane + 1u ? 1 : 0;
        row->contributes = (uint8_t)(row->assigned_allele != 0);
        row->reference_lane = (uint8_t)(!lane_has_contributor);
        row->cds = state->cds_bytes + ((hypothesis * state->transcript_count + transcript_index) * 2u + lane) * ARRANGEMENT_SEQUENCE_BYTES;
        row->protein = state->protein_bytes + ((hypothesis * state->transcript_count + transcript_index) * 2u + lane) * ARRANGEMENT_SEQUENCE_BYTES;
        row->cds_length = prediction->cds_length; row->protein_length = prediction->protein_length;
        row->nominal_length_diff = prediction->nominal_length_diff;
        row->consequence_mask = prediction->consequence_mask;
        row->prediction_status = prediction->prediction_status;
        row->prediction_reason = prediction->prediction_reason;
        row->hypothetical = prediction->hypothetical;
        return 1;
    }
    return 0;
}

size_t duckvep_arrangement_count(const duckvep_arrangement_state_t *state) {
    return state ? state->arrangement_count : 0u;
}

size_t duckvep_arrangement_workspace_bytes(const duckvep_arrangement_state_t *state) {
    return state ? state->workspace_bytes : 0u;
}

const char *duckvep_arrangement_prediction_status(uint8_t status) {
    switch ((duckvep_prediction_status_t)status) {
    case DUCKVEP_PREDICTION_ELIGIBLE: return "eligible_classifier_pending";
    case DUCKVEP_PREDICTION_PREDICTED: return "predicted";
    case DUCKVEP_PREDICTION_INCOMPLETE_INPUT: return "incomplete_input";
    case DUCKVEP_PREDICTION_EDIT_CONFLICT: return "edit_conflict";
    case DUCKVEP_PREDICTION_UNSUPPORTED_OVERLAP: return "unsupported_overlap";
    case DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT: return "unsupported_context";
    default: return "unsupported";
    }
}
