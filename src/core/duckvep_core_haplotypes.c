/* The phased native replay stream behind duckvep_haplotypes, host-neutral (see duckvep_core_haplotypes.h). The workspace is
 * allocated once and retains only the active window. DuckDB is reached only through the duckvep_h_* host layer of
 * "duckvep_host.h"; both the v1 and the v2 host compile this file unchanged. */
#include "duckvep_host.h"

#include "core/duckvep_core_haplotypes.h"
#include "kernel/src/duckvep_budget.h"
#include "kernel/src/duckvep_haplotype_stream.h"
#include "kernel/src/duckvep_sequence_diff.h"
#include "kernel/src/duckvep_dna.h"
#include "kernel/src/duckvep_effect.h"
#include "kernel/src/duckvep_hgvs.h"
#include "kernel/src/duckvep_annotation_internal.h"
#include "kernel/src/duckvep_event.h"
#include "duckvep_reference.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The last existing column stays last: nominal_length_diff is documented and tested as the
 * final field, so slice-2 columns are inserted before it. */
enum { HAPLOTYPE_LIST_COLUMN = 9, HAPLOTYPE_STOP_COLUMN = 14,
    HAPLOTYPE_HGVSP_COLUMN = 15, HAPLOTYPE_HGVSP_STATUS_COLUMN = 16,
    HAPLOTYPE_POLICY_COLUMN = 17, HAPLOTYPE_STATUS_COLUMN = 18, HAPLOTYPE_REASON_COLUMN = 19,
    HAPLOTYPE_PROVENANCE_COLUMN = 20, HAPLOTYPE_EDITS_COLUMN = 21,
    HAPLOTYPE_CARRIER_PREDICTION_COLUMN = 22, HAPLOTYPE_CONSEQUENCES_COLUMN = 23,
    HAPLOTYPE_IMPACT_COLUMN = 24,
    HAPLOTYPE_NMD_RULE_COLUMN = 25, HAPLOTYPE_NMD_COLUMN = 26, HAPLOTYPE_NMD_STOP_COLUMN = 27,
    HAPLOTYPE_NMD_JUNCTION_COLUMN = 28, HAPLOTYPE_NMD_CONTRIBUTORS_COLUMN = 29,
    HAPLOTYPE_NOMINAL_LENGTH_COLUMN = 30, HAPLOTYPE_OUTPUT_COLUMNS = 31 };
enum { HAPLOTYPE_PROVENANCE_FIELDS = 10, HAPLOTYPE_EDIT_FIELDS = 7, HAPLOTYPE_CARRIER_PREDICTION_FIELDS = 10 };
#define HAPLOTYPE_POLICY_VERSION "duckvep-coding-v1"
enum { HAPLOTYPE_BLOCK_EVENT_FIELD = 9, HAPLOTYPE_BLOCK_FIELDS = 10 };

const char *const duckvep_hap_limit_names[DUCKVEP_HAP_LIMIT_COUNT] = {"max_active_events", "max_active_transcripts",
    "max_active_carriers", "max_active_prefixes", "max_active_projections", "max_allele_bytes",
    "max_leaf_events", "max_leaf_edits", "max_sequence_bases", "max_ploidy", "max_phase_sets",
    "max_alignment_cells", "max_leaf_differences", "max_hgvs_operations", "max_hgvs_bytes",
    "max_hgvs_reference_bytes", "workspace_limit"};
const uint64_t duckvep_hap_limit_defaults[DUCKVEP_HAP_LIMIT_COUNT] = {16384, 4096, 65536, 262144, 262144, 8388608,
    4096, 65536, 1048576, 64, 1024, 16777216, 65536, 65536, 1048576,
    DUCKVEP_REFERENCE_DEFAULT_BYTES, 268435456};

int duckvep_hap_limit_valid(unsigned index, uint64_t n) {
    if (index >= DUCKVEP_HAP_LIMIT_COUNT || !n || n > SIZE_MAX) return 0;
    if (index <= DUCKVEP_HAP_LIMIT_PROJECTIONS && n > (UINT32_C(1) << 29)) return 0;
    if (index == DUCKVEP_HAP_LIMIT_PLOIDY && n > UINT16_MAX) return 0;
    return 1;
}

void duckvep_hap_config_defaults(duckvep_hap_config_t *config, const duckvep_owned_model_t *model) {
    memset(config, 0, sizeof *config);
    config->model = model;
    config->policy = DUCKVEP_PHASE_STRICT;
    for (unsigned i = 0u; i < DUCKVEP_HAP_LIMIT_COUNT; i++) config->limits[i] = (size_t)duckvep_hap_limit_defaults[i];
}

int duckvep_hap_config_check(const duckvep_hap_config_t *config, char *error, size_t error_size) {
    if (!config->model) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: require a nonempty calls query and loaded model name");
        return 0;
    }
    if (config->model->lifted) {
        duckvep_sql_set_error(error, error_size,
            "duckvep_haplotypes: phased edit sets are not supported for models with wrapped circular objects");
        return 0;
    }
    if (config->source_records && config->policy != DUCKVEP_PHASE_VEP116_COMPAT) {
        duckvep_sql_set_error(error, error_size,
            "duckvep_haplotypes: input_mode must be 'alt_events' or 'source_records'; source_records requires phase_policy='vep_compat'");
        return 0;
    }
    for (unsigned i = 0u; i < DUCKVEP_HAP_LIMIT_COUNT; i++) {
        if (!duckvep_hap_limit_valid(i, config->limits[i])) {
            char text[128];
            snprintf(text, sizeof text, "duckvep_haplotypes: invalid %s", duckvep_hap_limit_names[i]);
            duckvep_sql_set_error(error, error_size, text);
            return 0;
        }
    }
    return 1;
}

/* The output schema: names and SQL types. Record types are shared by the list columns. */
#define HAP_CARRIER "STRUCT(sample_index UINTEGER, phase_set BIGINT, haplotype_lane USMALLINT, ploidy USMALLINT)[]"
#define HAP_CONTRIBUTOR_7 "event_index UBIGINT, seq_region UINTEGER, position UBIGINT, reference VARCHAR, " \
    "alternate VARCHAR, evidence_flags UTINYINT, projection_status VARCHAR"
#define HAP_DIFFERENCE "STRUCT(ref_start0 UBIGINT, alt_start0 UBIGINT, reference VARCHAR, alternate VARCHAR, " \
    "alignment_start0 UBIGINT)[]"
static const char *const column_names[HAPLOTYPE_OUTPUT_COLUMNS] = {"transcript_index", "cds", "protein",
    "sequence_flags", "evidence_flags", "projection_status", "sequence_status", "edit_count", "carrier_count",
    "carriers", "contributors", "coding_blocks", "cds_differences", "protein_differences",
    "stop_in_displaced_frame", "hgvsp", "hgvsp_status", "prediction_policy", "prediction_status",
    "prediction_reason", "contributor_provenance", "normalized_edits", "carrier_predictions",
    "haplotype_consequences", "haplotype_impact", "nmd_rule", "nmd_prediction", "nmd_stop_position",
    "nmd_junction_position", "nmd_contributors", "nominal_length_diff"};
static const char *const column_types[HAPLOTYPE_OUTPUT_COLUMNS] = {"UINTEGER", "VARCHAR", "VARCHAR",
    "UINTEGER", "UTINYINT", "VARCHAR", "VARCHAR", "UBIGINT", "UINTEGER",
    HAP_CARRIER, NULL,
    "STRUCT(cds_start UINTEGER, reference VARCHAR, alternate VARCHAR, alt_start0 UBIGINT, length_change BIGINT, "
    "sequence_flags UINTEGER, coding_status VARCHAR, local_consequence_mask UBIGINT, after_first_stop BOOLEAN, "
    "event_indices UBIGINT[])[]",
    HAP_DIFFERENCE, HAP_DIFFERENCE,
    "BOOLEAN", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR",
    "STRUCT(event_index UBIGINT, alt_index UINTEGER, seq_region UINTEGER, position UBIGINT, reference VARCHAR, "
    "alternate VARCHAR, evidence_flags UTINYINT, projection_status VARCHAR, role VARCHAR, edit_count UINTEGER)[]",
    "STRUCT(edit_index UBIGINT, event_index UBIGINT, block_index UBIGINT, cds_start UINTEGER, reference VARCHAR, "
    "alternate VARCHAR, variant_strand TINYINT)[]",
    "STRUCT(sample_index UINTEGER, phase_set BIGINT, haplotype_lane USMALLINT, prediction_status VARCHAR, "
    "prediction_reason VARCHAR, haplotype_impact VARCHAR, haplotype_consequences VARCHAR[], nmd_prediction VARCHAR, "
    "nmd_stop_position UBIGINT, nmd_junction_position UBIGINT)[]",
    "VARCHAR[]", "VARCHAR", "VARCHAR", "VARCHAR", "UBIGINT", "UBIGINT", "UBIGINT[]", "BIGINT"};

const char *duckvep_hap_column_name(unsigned index) {
    return index < HAPLOTYPE_OUTPUT_COLUMNS ? column_names[index] : NULL;
}

const char *duckvep_hap_column_type(unsigned index, int source_records) {
    if (index >= HAPLOTYPE_OUTPUT_COLUMNS) return NULL;
    if (index == 10u)
        return source_records ? "STRUCT(" HAP_CONTRIBUTOR_7 ", alt_index UINTEGER)[]" : "STRUCT(" HAP_CONTRIBUTOR_7 ")[]";
    return column_types[index];
}

typedef duckvep_hap_config_t haplotype_bind_t;

struct duckvep_hap_state {
    duckvep_hap_config_t bind;
    duckvep_hap_row_t row;
    int have_row, eof, have_call;
    uint32_t last_tx;
    duckvep_haplotype_stream_t stream;
    duckvep_haplotype_stream_buffers_t buffers;
    int32_t *gt;
    uint8_t *phase;
    duckvep_haplotype_phase_set_t *sets;
    duckvep_sequence_diff_scratch_t difference_scratch;
    duckvep_sequence_difference_t *differences;
    uint8_t *difference_reference;
    uint32_t difference_transcript;
    int have_difference_reference;
    duckvep_hgvs_protein_operation_t *protein_operations;
    duckvep_reference_reader_t reference;
    duckvep_delta_scratch_t hgvs_scratch;
    uint8_t *hgvs_shifted_allele;
    size_t hgvs_shifted_capacity;
    char *hgvsp;
    size_t workspace_bytes;
};

static void assign_text(duckvep_h_vector vector, size_t row, const char *text) {
    duckvep_h_assign_string(vector, row, text, strlen(text));
}

void duckvep_hap_close(duckvep_hap_state_t *s) {
    if (!s) return;
    duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    duckvep_budget_free(b->carriers.transcripts); duckvep_budget_free(b->carriers.calls); duckvep_budget_free(b->carriers.prefixes);
    duckvep_budget_free(b->carriers.active_transcripts); duckvep_budget_free(b->carriers.transcript_index);
    duckvep_budget_free(b->carriers.call_index); duckvep_budget_free(b->carriers.prefix_index);
    duckvep_budget_free(b->events); duckvep_budget_free(b->projections); duckvep_budget_free(b->alleles); duckvep_budget_free(b->leaf_events);
    duckvep_budget_free(b->contributors); duckvep_budget_free(b->edits); duckvep_budget_free(b->blocks); duckvep_budget_free(b->cds); duckvep_budget_free(b->protein);
    duckvep_budget_free(b->edit_event_ids);
    duckvep_budget_free(s->difference_scratch.scores); duckvep_budget_free(s->difference_scratch.trace); duckvep_budget_free(s->differences);
    duckvep_budget_free(s->difference_reference); duckvep_budget_free(b->reference_protein);
    duckvep_budget_free(b->reference_coding_protein);
    duckvep_budget_free(s->protein_operations); duckvep_budget_free(s->hgvsp);
    duckvep_reference_reader_close(&s->reference);
    duckvep_budget_free(s->reference.bases);
    duckvep_budget_free(s->hgvs_scratch.edits); duckvep_budget_free(s->hgvs_scratch.alt_cds);
    duckvep_budget_free(s->hgvs_scratch.ref_peptide); duckvep_budget_free(s->hgvs_scratch.alt_peptide);
    duckvep_budget_free(s->hgvs_shifted_allele);
    duckvep_budget_free(s->gt); duckvep_budget_free(s->phase); duckvep_budget_free(s->sets); duckvep_budget_free(s);
}

static uint32_t bucket_count(uint32_t capacity) {
    uint32_t n = 2u;
    while (n < capacity * 2u) n *= 2u;
    return n;
}

static int workspace_allocate(duckvep_hap_state_t *s, const duckvep_hap_config_t *bind) {
    const size_t *n = bind->limits;
    duckvep_haplotype_stream_buffers_t *b = &s->buffers;
    b->carriers.transcript_capacity = (uint32_t)n[DUCKVEP_HAP_LIMIT_TRANSCRIPTS];
    b->carriers.call_capacity = (uint32_t)n[DUCKVEP_HAP_LIMIT_CARRIERS];
    b->carriers.prefix_capacity = (uint32_t)n[DUCKVEP_HAP_LIMIT_PREFIXES];
    b->carriers.transcript_buckets = bucket_count(b->carriers.transcript_capacity);
    b->carriers.call_buckets = bucket_count(b->carriers.call_capacity);
    b->carriers.prefix_buckets = bucket_count(b->carriers.prefix_capacity);
    b->event_capacity = (uint32_t)n[DUCKVEP_HAP_LIMIT_EVENTS]; b->projection_capacity = (uint32_t)n[DUCKVEP_HAP_LIMIT_PROJECTIONS];
    b->allele_capacity = n[DUCKVEP_HAP_LIMIT_ALLELES]; b->leaf_capacity = n[DUCKVEP_HAP_LIMIT_LEAF_EVENTS];
    b->edit_capacity = n[DUCKVEP_HAP_LIMIT_LEAF_EDITS];
    if (n[DUCKVEP_HAP_LIMIT_SEQUENCE] == SIZE_MAX || (bind->hgvs && n[DUCKVEP_HAP_LIMIT_HGVS_BYTES] == SIZE_MAX)) return 0;
    b->cds_capacity = n[DUCKVEP_HAP_LIMIT_SEQUENCE];
    b->protein_capacity = n[DUCKVEP_HAP_LIMIT_SEQUENCE] + 1u;
    if (b->protein_capacity > SIZE_MAX / 2u) return 0;
    s->difference_scratch.score_capacity = b->protein_capacity * 2u;
    s->difference_scratch.trace_capacity = n[DUCKVEP_HAP_LIMIT_ALIGNMENT];
    b->reference_protein_capacity = n[DUCKVEP_HAP_LIMIT_SEQUENCE] / 3u + 2u;
    if (bind->hgvs) {
        s->reference.capacity = bind->model->reference_fasta_path ? n[DUCKVEP_HAP_LIMIT_HGVS_REFERENCE] : 0u;
        /* Alternating differing/retained bases maximize one uint16-length MNV's islands. */
        size_t islands = ((size_t)UINT16_MAX + 1u) / 2u;
        s->hgvs_scratch.edits_cap = n[DUCKVEP_HAP_LIMIT_LEAF_EDITS] < islands ? n[DUCKVEP_HAP_LIMIT_LEAF_EDITS] : islands;
        s->hgvs_scratch.alt_cds_cap = n[DUCKVEP_HAP_LIMIT_SEQUENCE];
        s->hgvs_scratch.ref_peptide_cap = s->hgvs_scratch.alt_peptide_cap = n[DUCKVEP_HAP_LIMIT_SEQUENCE] / 3u + 1u;
        s->hgvs_shifted_capacity = n[DUCKVEP_HAP_LIMIT_ALLELES] < UINT16_MAX + 1u ? n[DUCKVEP_HAP_LIMIT_ALLELES] : UINT16_MAX + 1u;
    }
    /* Count every byte before allocating any of the arrays. Each pointer has
     * one owner and one cleanup site; no allocator is used by the scan loop. */
#define ARRAYS(X) \
    X(b->carriers.transcripts, b->carriers.transcript_capacity) \
    X(b->carriers.calls, b->carriers.call_capacity) \
    X(b->carriers.prefixes, b->carriers.prefix_capacity) \
    X(b->carriers.active_transcripts, b->carriers.transcript_capacity) \
    X(b->carriers.transcript_index, b->carriers.transcript_buckets) \
    X(b->carriers.call_index, b->carriers.call_buckets) \
    X(b->carriers.prefix_index, b->carriers.prefix_buckets) \
    X(b->events, b->event_capacity) X(b->projections, b->projection_capacity) \
    X(b->alleles, b->allele_capacity) X(b->leaf_events, b->leaf_capacity) \
    X(b->contributors, b->leaf_capacity) X(b->edits, b->edit_capacity) \
    X(b->blocks, b->edit_capacity) X(b->edit_event_ids, b->edit_capacity) \
    X(b->cds, b->cds_capacity) X(b->protein, b->protein_capacity) \
    X(s->difference_scratch.scores, s->difference_scratch.score_capacity) \
    X(s->difference_scratch.trace, s->difference_scratch.trace_capacity) \
    X(s->differences, n[DUCKVEP_HAP_LIMIT_DIFFERENCES]) \
    X(s->difference_reference, n[DUCKVEP_HAP_LIMIT_SEQUENCE]) \
    X(b->reference_protein, b->reference_protein_capacity) \
    X(b->reference_coding_protein, b->reference_protein_capacity) \
    X(s->protein_operations, bind->hgvs ? n[DUCKVEP_HAP_LIMIT_HGVS_OPERATIONS] : 0u) \
    X(s->hgvsp, bind->hgvs ? n[DUCKVEP_HAP_LIMIT_HGVS_BYTES] + 1u : 0u) \
    X(s->reference.bases, s->reference.capacity) \
    X(s->hgvs_scratch.edits, s->hgvs_scratch.edits_cap) \
    X(s->hgvs_scratch.alt_cds, s->hgvs_scratch.alt_cds_cap) \
    X(s->hgvs_scratch.ref_peptide, s->hgvs_scratch.ref_peptide_cap) \
    X(s->hgvs_scratch.alt_peptide, s->hgvs_scratch.alt_peptide_cap) \
    X(s->hgvs_shifted_allele, s->hgvs_shifted_capacity) \
    X(s->gt, bind->source_records ? 0u : n[DUCKVEP_HAP_LIMIT_PLOIDY]) \
    X(s->phase, bind->source_records ? 0u : n[DUCKVEP_HAP_LIMIT_PLOIDY]) \
    X(s->sets, bind->source_records ? 0u : n[DUCKVEP_HAP_LIMIT_PHASE_SETS])
#define COUNT(p, count) \
    if ((count) > (n[DUCKVEP_HAP_LIMIT_WORKSPACE] - s->workspace_bytes) / sizeof(*(p))) return 0; \
    s->workspace_bytes += (count) * sizeof(*(p));
    s->workspace_bytes = sizeof(*s);
    if (s->workspace_bytes > n[DUCKVEP_HAP_LIMIT_WORKSPACE]) return 0;
    ARRAYS(COUNT)
#undef COUNT
#define ALLOCATE(p, count) if ((count) && !((p) = duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE, (count) * sizeof(*(p))))) return 0;
    ARRAYS(ALLOCATE)
#undef ALLOCATE
#undef ARRAYS
    return 1;
}
static const char *projection_name(duckvep_cds_edit_status_t status) {
    switch (status) {
    case DUCKVEP_CDS_EDIT_OK: return "ok";
    case DUCKVEP_CDS_EDIT_REF_MISMATCH: return "reference_mismatch";
    case DUCKVEP_CDS_EDIT_SOURCE_SHADOWED: return "shadowed_duplicate";
    case DUCKVEP_CDS_EDIT_SOURCE_UNMAPPED: return "source_unmapped";
    case DUCKVEP_CDS_EDIT_SOURCE_ALLELE_SKIPPED: return "source_allele_skipped";
    case DUCKVEP_CDS_EDIT_INVALID_ARG: return "invalid_argument";
    case DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND: return "unsupported_kind";
    case DUCKVEP_CDS_EDIT_INVALID_EVENT: return "invalid_event";
    case DUCKVEP_CDS_EDIT_OUT_OF_CDS: return "outside_cds";
    case DUCKVEP_CDS_EDIT_NON_CONTIGUOUS: return "non_contiguous";
    case DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL: return "projection_capacity";
    case DUCKVEP_CDS_EDIT_INVALID_ALLELE: return "invalid_allele";
    default: return "unavailable_projection";
    }
}

static const char *sequence_name(duckvep_haplotype_status_t status) {
    switch (status) {
    case DUCKVEP_HAPLOTYPE_OK: return "ok";
    case DUCKVEP_HAPLOTYPE_CONDITIONAL: return "conditional";
    case DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE: return "incomplete_input";
    case DUCKVEP_HAPLOTYPE_EDIT_ORDER: return "edit_conflict";
    case DUCKVEP_HAPLOTYPE_REF_MISMATCH: return "reference_mismatch";
    case DUCKVEP_HAPLOTYPE_INVALID_BASE: return "invalid_base";
    default: return "invalid_sequence";
    }
}

static int exhausted_limit(duckvep_haplotype_stream_status_t status, duckvep_carriers_status_t carrier) {
    switch (status) {
    case DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL: return DUCKVEP_HAP_LIMIT_EVENTS;
    case DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL: return DUCKVEP_HAP_LIMIT_PROJECTIONS;
    case DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL: return DUCKVEP_HAP_LIMIT_ALLELES;
    case DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL: return DUCKVEP_HAP_LIMIT_LEAF_EVENTS;
    case DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL: return DUCKVEP_HAP_LIMIT_LEAF_EDITS;
    case DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL: return DUCKVEP_HAP_LIMIT_SEQUENCE;
    case DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR:
        if (carrier == DUCKVEP_CARRIERS_TRANSCRIPT_FULL) return DUCKVEP_HAP_LIMIT_TRANSCRIPTS;
        if (carrier == DUCKVEP_CARRIERS_CALL_FULL) return DUCKVEP_HAP_LIMIT_CARRIERS;
        if (carrier == DUCKVEP_CARRIERS_PREFIX_FULL) return DUCKVEP_HAP_LIMIT_PREFIXES;
        return -1;
    default: return -1;
    }
}

static void null_cell(duckvep_h_vector vector, size_t row) {
    duckvep_h_set_null(vector, row);
}

static int prepare_difference_reference(duckvep_hap_state_t *s, const duckvep_hap_config_t *bind,
    const duckvep_haplotype_leaf_t *leaf, char *error, size_t error_size) {
    uint32_t tx = leaf->carriers.transcript_index;
    if (!leaf->cds || (s->have_difference_reference && tx == s->difference_transcript)) return 1;
    const duckvep_sequence_pool_t *seq = s->stream.sequences;
    size_t length = seq->cds_length[tx];
    if (length > bind->limits[DUCKVEP_HAP_LIMIT_SEQUENCE]) {
        snprintf(error, error_size,
            "duckvep_haplotypes: max_sequence_bases=%zu, reference requires=%zu at transcript %u",
            bind->limits[DUCKVEP_HAP_LIMIT_SEQUENCE], length, tx);
        return 0;
    }
    /* CDS alignment uses replay's canonical spelling; reference protein
     * preparation retains the model bytes for Ensembl's exact stop convention. */
    for (size_t i = 0u; i < length; i++) {
        char base = duckvep_dna_normalize_n((char)leaf->reference_cds[i]);
        if (!base) {
            duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid reference CDS base");
            return 0;
        }
        s->difference_reference[i] = (uint8_t)base;
    }
    s->difference_transcript = tx; s->have_difference_reference = 1;
    return 1;
}

/* Both sequence axes reuse the same bounded traceback and descriptor storage.
 * DuckDB copies one list's spans before the next axis resets those descriptors. */
static int append_sequence_differences(duckvep_h_vector vector, size_t row, duckvep_hap_state_t *s,
    const duckvep_hap_config_t *bind, const duckvep_haplotype_leaf_t *leaf, int protein,
    char *error, size_t error_size) {
    int known = leaf->cds && (!protein || leaf->reference_protein);
    const uint8_t *reference = protein ? leaf->reference_protein : s->difference_reference;
    const uint8_t *alternate = protein ? leaf->protein : leaf->cds;
    duckvep_sequence_diff_result_t result = {0};
    if (known) {
        size_t ref_length = protein ? leaf->reference_protein_length
            : s->stream.sequences->cds_length[leaf->carriers.transcript_index];
        size_t alt_length = protein ? leaf->protein_length : leaf->cds_length;
        duckvep_sequence_diff_status_t status = duckvep_sequence_differences(reference, ref_length,
            alternate, alt_length, (leaf->flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u,
            &s->difference_scratch, s->differences, bind->limits[DUCKVEP_HAP_LIMIT_DIFFERENCES], &result);
        if (status != DUCKVEP_SEQUENCE_DIFF_OK) {
            int limit = status == DUCKVEP_SEQUENCE_DIFF_TRACE_FULL ? DUCKVEP_HAP_LIMIT_ALIGNMENT :
                status == DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL ? DUCKVEP_HAP_LIMIT_DIFFERENCES : DUCKVEP_HAP_LIMIT_SEQUENCE;
            size_t required = limit == DUCKVEP_HAP_LIMIT_ALIGNMENT ? result.trace_cells :
                limit == DUCKVEP_HAP_LIMIT_DIFFERENCES ? result.count : alt_length;
            snprintf(error, error_size,
                "duckvep_haplotypes: %s difference status %u, %s=%zu, required=%zu at transcript %u",
                protein ? "protein" : "CDS", (unsigned)status, duckvep_hap_limit_names[limit],
                bind->limits[limit], required, leaf->carriers.transcript_index);
            return 0;
        }
    }
    duckvep_h_list_entry entry;
    if (!duckvep_h_list_extend(vector, result.count, &entry)) return 0;
    size_t base = entry.offset;
    ((duckvep_h_list_entry *)duckvep_h_data(vector))[row] = entry;
    if (!known) null_cell(vector, row);
    duckvep_h_vector records = duckvep_h_list_child(vector), fields[5];
    duckvep_h_ensure_validity(records);
    for (unsigned j = 0u; j < 5u; j++) {
        fields[j] = duckvep_h_struct_child(records, j);
        duckvep_h_ensure_validity(fields[j]);
    }
    for (size_t i = 0u; i < result.count; i++) {
        size_t at = base + i;
        duckvep_h_set_valid(records, at);
        for (unsigned j = 0u; j < 5u; j++)
            duckvep_h_set_valid(fields[j], at);
        const duckvep_sequence_difference_t *d = &s->differences[i];
        ((uint64_t *)duckvep_h_data(fields[0]))[at] = d->ref_start0;
        ((uint64_t *)duckvep_h_data(fields[1]))[at] = d->alt_start0;
        duckvep_h_assign_string(fields[2], at,
            (const char *)reference + d->ref_start0, d->ref_length);
        duckvep_h_assign_string(fields[3], at,
            (const char *)alternate + d->alt_start0, d->alt_length);
        ((uint64_t *)duckvep_h_data(fields[4]))[at] = d->alignment_start0;
    }
    return 1;
}

typedef struct {
    duckvep_hap_state_t *state;
    const duckvep_hap_config_t *bind;
    duckvep_hgvs_status_t status;
    size_t length;
    char *error;
    size_t error_size;
} haplotype_hgvs_observer_t;

static int haplotype_hgvs_observe(void *pointer, const duckvep_variant_batch_t *variants,
    const duckvep_consequence_t *row, const duckvep_pair_facts_t *facts) {
    haplotype_hgvs_observer_t *o = pointer;
    duckvep_hap_state_t *s = o->state;
    const duckvep_owned_model_t *m = o->bind->model;
    if (!facts || !facts->event) return 0;
    if (facts->transcript_edit_status != DUCKVEP_TRANSCRIPT_EDIT_OK) {
        o->status = facts->transcript_edit_status == DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT
            ? DUCKVEP_HGVS_NOT_APPLICABLE : DUCKVEP_HGVS_INVALID_PROJECTION;
        return 1;
    }
    int available;
    duckvep_hgvs_reference_window_t shift, lookup;
    if (!duckvep_reference_reader_windows(&s->reference, facts->event, &available,
            &shift, &lookup, o->error, o->error_size)) return 0;
    o->status = available ? duckvep_hgvs_uploaded_reference_validate(&lookup, facts->event,
        variants->allele_bytes + variants->ref_offset[0], variants->ref_length[0]) : DUCKVEP_HGVS_OK;
    if (o->status != DUCKVEP_HGVS_OK) return 1;
    duckvep_transcript_edit_t edit;
    duckvep_hgvs_dna_fact_t dna;
    o->status = duckvep_hgvs_dna_pair_build(&m->transcripts, &m->exons, &m->sequences,
        variants, facts, available ? &shift : NULL, available ? &lookup : NULL,
        &s->hgvs_scratch, &edit, &dna);
    if (o->status != DUCKVEP_HGVS_OK) return 1;
    duckvep_pair_facts_t protein_facts = *facts;
    protein_facts.transcript_edit = &edit;
    duckvep_hgvs_protein_pair_t protein;
    size_t required;
    o->status = duckvep_hgvs_protein_pair_build(&m->transcripts, &m->exons, &m->sequences,
        variants, row, &protein_facts, &dna, available ? &lookup : NULL,
        &s->hgvs_scratch, s->hgvs_shifted_allele, s->hgvs_shifted_capacity, &required, &protein);
    if (o->status != DUCKVEP_HGVS_OK) return 1;
    o->status = duckvep_hgvs_protein_render(&protein.fact, 1, s->hgvsp,
        o->bind->limits[DUCKVEP_HAP_LIMIT_HGVS_BYTES] + 1u, &o->length);
    if (o->status == DUCKVEP_HGVS_BUFFER_TOO_SMALL) {
        snprintf(o->error, o->error_size,
            "duckvep_haplotypes: max_hgvs_bytes=%zu, required=%zu at transcript %u",
            o->bind->limits[DUCKVEP_HAP_LIMIT_HGVS_BYTES], o->length, row->tx_idx);
        return 0;
    }
    return 1;
}

static int append_single_event_hgvsp(duckvep_h_vector text, duckvep_h_vector status_vector, size_t row,
    duckvep_hap_state_t *s, const duckvep_hap_config_t *bind, const duckvep_haplotype_leaf_t *leaf,
    char *error, size_t error_size) {
    const duckvep_haplotype_contributor_t *contributor = &leaf->contributors[0];
    const duckvep_haplotype_source_t *source = &contributor->source;
    size_t bytes = (size_t)source->ref_len + source->alt_len;
    const duckvep_event_t *event = contributor->prepared;
    duckvep_event_t normalized = {0};
    if (source->source_record) {
        if (!duckvep_event_prepare_small(source->pos1, source->ref, source->ref_len,
                source->alt, source->alt_len, &normalized)) {
            duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid source allele for HGVS");
            return 0;
        }
        normalized.chrom_id = source->chrom_id;
        event = &normalized;
    }
    uint32_t ref_offset = 0u, alt_offset = source->ref_len;
    duckvep_variant_batch_t variant = {.chrom_id = &source->chrom_id, .pos1 = &source->pos1,
        .end1 = &event->raw_end1, .ref_offset = &ref_offset, .alt_offset = &alt_offset,
        .ref_length = &source->ref_len, .alt_length = &source->alt_len,
        .allele_bytes = source->ref, .allele_bytes_len = bytes, .variant_kind = &event->kind,
        .count = 1u};
    haplotype_hgvs_observer_t observer = {.state = s, .bind = bind,
        .status = DUCKVEP_HGVS_NOT_APPLICABLE, .error = error, .error_size = error_size};
    duckvep_error_t native_error = {0};
    if (duckvep_annotate_pair_observed(bind->model->kernel, &variant, event,
            leaf->carriers.transcript_index, &s->hgvs_scratch,
            source->source_record ? NULL : contributor->projected, haplotype_hgvs_observe,
            &observer, &native_error) != DUCKVEP_OK) {
        if (!error[0]) duckvep_sql_set_error(error, error_size, native_error.message);
        return 0;
    }
    const char *name;
    switch (observer.status) {
        case DUCKVEP_HGVS_OK:
            duckvep_h_set_valid(text, row);
            duckvep_h_assign_string(text, row, s->hgvsp, observer.length);
            name = "ok"; break;
        case DUCKVEP_HGVS_NOT_APPLICABLE: name = "not_applicable"; break;
        case DUCKVEP_HGVS_MISSING_REFERENCE: name = "missing_reference"; break;
        case DUCKVEP_HGVS_REFERENCE_MISMATCH: name = "reference_mismatch"; break;
        case DUCKVEP_HGVS_INVALID_ALLELE: name = "invalid_allele"; break;
        case DUCKVEP_HGVS_MISSING_PEPTIDE: name = "missing_peptide"; break;
        case DUCKVEP_HGVS_MISSING_TRANSCRIPT_TAIL: name = "missing_transcript_tail"; break;
        case DUCKVEP_HGVS_MISSING_TRANSCRIPT_FLANK: name = "missing_transcript_flank"; break;
        case DUCKVEP_HGVS_UNSUPPORTED_PROTEIN: name = "unsupported_protein"; break;
        case DUCKVEP_HGVS_UNSUPPORTED_EDIT: name = "unsupported_coding_context"; break;
        case DUCKVEP_HGVS_INVALID_PROJECTION: name = "invalid_projection"; break;
        default:
            snprintf(error, error_size,
                "duckvep_haplotypes: HGVS status %u at event %llu; max_sequence_bases=%zu, max_leaf_edits=%zu, max_allele_bytes=%zu",
                (unsigned)observer.status, (unsigned long long)source->event_id,
                bind->limits[DUCKVEP_HAP_LIMIT_SEQUENCE], bind->limits[DUCKVEP_HAP_LIMIT_LEAF_EDITS], bind->limits[DUCKVEP_HAP_LIMIT_ALLELES]);
            return 0;
    }
    assign_text(status_vector, row, name);
    return 1;
}

static int append_hgvsp(duckvep_h_vector text, duckvep_h_vector status_vector, size_t row,
    duckvep_hap_state_t *s, const duckvep_hap_config_t *bind, const duckvep_haplotype_leaf_t *leaf,
    const duckvep_coding_context_t *coding, char *error, size_t error_size) {
    null_cell(text, row);
    const char *name = !bind->hgvs ? "not_requested" :
        !leaf->cds || !leaf->protein ? "unavailable_sequence" :
        leaf->sequence_status != DUCKVEP_HAPLOTYPE_OK ? "incomplete_input" :
        !leaf->reference_protein || !leaf->reference_protein_length ? "missing_reference_protein" :
        leaf->ordered_replacements ? "unsupported_ordered_replacements" : NULL;
    /* A singleton source allele uses the independent-event VEP presentation,
     * including its original anchor and cached predicates. Physical MNV islands
     * do not turn that one source into a compound event. */
    if (!name && leaf->contributor_count == 1u) {
        const duckvep_haplotype_contributor_t *c = &leaf->contributors[0];
        if ((!c->source.source_record || (c->source.allele_index && c->source.allele_index != UINT32_MAX)) &&
            (c->projection_status == DUCKVEP_CDS_EDIT_OK || c->projection_status == DUCKVEP_CDS_EDIT_OUT_OF_CDS))
            return append_single_event_hgvsp(text, status_vector, row, s, bind, leaf, error, error_size);
    }
    size_t count = 0u;
    /* The raw reference-only route already borrows the prepared reference.
     * No codon replay is needed to establish equality of these exact operands. */
    int reference_only = !name && !leaf->edit_count && leaf->protein == leaf->reference_protein;
    if (!name && !reference_only) {
        duckvep_hgvs_protein_reference_t reference = {
            leaf->reference_protein, leaf->reference_protein_length};
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_build(coding,
            &reference, s->buffers.edits, leaf->edit_count, leaf->blocks, leaf->block_count,
            bind->model->transcripts.flags[leaf->carriers.transcript_index],
            s->protein_operations, bind->limits[DUCKVEP_HAP_LIMIT_HGVS_OPERATIONS], &count);
        if (status == DUCKVEP_HGVS_BUFFER_TOO_SMALL) {
            snprintf(error, error_size,
                "duckvep_haplotypes: max_hgvs_operations=%zu exhausted at transcript %u",
                bind->limits[DUCKVEP_HAP_LIMIT_HGVS_OPERATIONS], leaf->carriers.transcript_index);
            return 0;
        }
        switch (status) {
            case DUCKVEP_HGVS_OK: break;
            case DUCKVEP_HGVS_NOT_APPLICABLE: name = "not_applicable"; break;
            case DUCKVEP_HGVS_MISSING_PEPTIDE: name = "missing_peptide"; break;
            case DUCKVEP_HGVS_UNSUPPORTED_PROTEIN: name = "unsupported_protein"; break;
            case DUCKVEP_HGVS_UNSUPPORTED_EDIT: name = "unsupported_coding_context"; break;
            default:
                snprintf(error, error_size, "duckvep_haplotypes: HGVS status %u at transcript %u",
                    (unsigned)status, leaf->carriers.transcript_index);
                return 0;
        }
    }
    if (!name) {
        size_t required;
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_haplotype_render(s->protein_operations,
            count, 1, s->hgvsp, bind->limits[DUCKVEP_HAP_LIMIT_HGVS_BYTES] + 1u, &required);
        if (status != DUCKVEP_HGVS_OK) {
            snprintf(error, error_size,
                "duckvep_haplotypes: HGVS render status %u, max_hgvs_bytes=%zu, required=%zu at transcript %u",
                (unsigned)status, bind->limits[DUCKVEP_HAP_LIMIT_HGVS_BYTES], required, leaf->carriers.transcript_index);
            return 0;
        }
        duckvep_h_set_valid(text, row);
        duckvep_h_assign_string(text, row, s->hgvsp, required);
        name = "ok";
    }
    assign_text(status_vector, row, name);
    return 1;
}

static const char *status_name(duckvep_prediction_status_t status) {
    switch (status) {
    case DUCKVEP_PREDICTION_ELIGIBLE: return "eligible_classifier_pending";
    case DUCKVEP_PREDICTION_PREDICTED: return "predicted";
    case DUCKVEP_PREDICTION_INCOMPLETE_INPUT: return "incomplete_input";
    case DUCKVEP_PREDICTION_EDIT_CONFLICT: return "edit_conflict";
    case DUCKVEP_PREDICTION_UNSUPPORTED_OVERLAP: return "unsupported_overlap";
    default: return "unsupported_context";
    }
}

static const char *reason_name(duckvep_prediction_reason_t reason, duckvep_cds_edit_status_t projection) {
    switch (reason) {
    case DUCKVEP_REASON_SUPPORTED_DOMAIN: return "supported_domain";
    case DUCKVEP_REASON_MISSING_CALL: return "missing_call";
    case DUCKVEP_REASON_UNPHASED_HETEROZYGOUS: return "unphased_heterozygous";
    case DUCKVEP_REASON_CROSS_PS_UNRESOLVED: return "unresolved_cross_ps_phase";
    case DUCKVEP_REASON_CONTRADICTORY_EDITS: return "contradictory_edits";
    case DUCKVEP_REASON_OVERLAPPING_EDITS: return "overlapping_edits";
    case DUCKVEP_REASON_DUPLICATE_EDITS: return "duplicate_edits";
    case DUCKVEP_REASON_SAME_GAP_INSERTIONS: return "ambiguous_same_gap_insertions";
    case DUCKVEP_REASON_NON_STRICT_PHASE_POLICY: return "non_strict_phase_policy";
    case DUCKVEP_REASON_NON_DIPLOID_CALL: return "non_diploid_call";
    case DUCKVEP_REASON_PROJECTION: return projection_name(projection);
    case DUCKVEP_REASON_REFERENCE_MISMATCH: return "reference_mismatch";
    case DUCKVEP_REASON_INVALID_BASE: return "invalid_base";
    case DUCKVEP_REASON_TRANSCRIPT_NOT_CODING: return "transcript_not_coding";
    case DUCKVEP_REASON_NON_STANDARD_CODON_TABLE: return "non_standard_codon_table";
    case DUCKVEP_REASON_CURATED_TRANSCRIPT: return "curated_transcript";
    case DUCKVEP_REASON_INCOMPLETE_CDS: return "incomplete_cds";
    case DUCKVEP_REASON_NONCANONICAL_START: return "noncanonical_start";
    case DUCKVEP_REASON_NONCANONICAL_STOP: return "noncanonical_stop";
    case DUCKVEP_REASON_INTERNAL_STOP: return "internal_stop";
    case DUCKVEP_REASON_NON_LITERAL_ALLELE: return "non_literal_allele";
    case DUCKVEP_REASON_ALLELE_OVER_50: return "allele_over_50_bases";
    case DUCKVEP_REASON_START_STOP_CLASSIFIER_PENDING: return "start_stop_classifier_pending";
    default: return "invalid_sequence";
    }
}

static const char *role_name(uint8_t role) {
    switch (role) {
    case DUCKVEP_ROLE_SHADOWED: return "shadowed";
    case DUCKVEP_ROLE_UNAPPLIED: return "unapplied";
    case DUCKVEP_ROLE_APPLIED: return "applied";
    case DUCKVEP_ROLE_POST_STOP: return "post_stop";
    default: return "omitted";
    }
}

/* NMD attribution: the applied edits (up to and including the stop) plus the post_stop edits that changed length at or
 * before the penultimate exon's last base, i.e. exactly the edits that moved J. Roles are unchanged. */
static int nmd_attributed(const duckvep_haplotype_contributor_t *c) {
    return c->role == DUCKVEP_ROLE_APPLIED || (c->role == DUCKVEP_ROLE_POST_STOP && c->nmd_moved_junction);
}

static const char *nmd_name(duckvep_haplotype_nmd_t nmd) {
    switch (nmd) {
    case DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE: return "not_applicable";
    case DUCKVEP_HAPLOTYPE_NMD_ESCAPE: return "escape";
    case DUCKVEP_HAPLOTYPE_NMD_TRIGGER: return "trigger";
    default: return "unknown";
    }
}

/* Versioned coding-v1 status/reason plus complete contributor and normalized-edit
 * provenance. Every contributor of the leaf is listed, whatever its role. */
static int append_prediction(duckvep_h_chunk output, size_t row, duckvep_hap_state_t *s,
    const duckvep_hap_config_t *bind, const duckvep_haplotype_leaf_t *leaf) {
    duckvep_h_vector v[HAPLOTYPE_OUTPUT_COLUMNS];
    for (unsigned i = HAPLOTYPE_POLICY_COLUMN; i <= HAPLOTYPE_CARRIER_PREDICTION_COLUMN; i++) {
        v[i] = duckvep_h_chunk_vector(output, i);
        duckvep_h_set_valid(v[i], row);
    }
    assign_text(v[HAPLOTYPE_POLICY_COLUMN], row, HAPLOTYPE_POLICY_VERSION);
    assign_text(v[HAPLOTYPE_STATUS_COLUMN], row, status_name(leaf->prediction_status));
    assign_text(v[HAPLOTYPE_REASON_COLUMN], row, reason_name(leaf->prediction_reason, leaf->prediction_projection));
    /* The shared edited sequence has one reduced SO set, ordered by severity rank; it is emitted for
     * every carrier whose own keyed status is predicted, and for the row only when all are. */
    unsigned bits[DUCKVEP_SO_BIT_COUNT], n = 0u;
    if (leaf->path_status == DUCKVEP_PREDICTION_PREDICTED) {
        for (unsigned bit = 0u; bit < DUCKVEP_SO_BIT_COUNT; bit++)
            if (leaf->haplotype_so_mask & DUCKVEP_SO(bit)) {
                unsigned at = n++;
                while (at && duckvep_so_rank(bits[at - 1u]) > duckvep_so_rank(bit)) { bits[at] = bits[at - 1u]; at--; }
                bits[at] = bit;
            }
    }
    const size_t counts[] = {leaf->contributor_count, leaf->listed_edit_count, leaf->carriers.call_count};
    const unsigned field_counts[] = {HAPLOTYPE_PROVENANCE_FIELDS, HAPLOTYPE_EDIT_FIELDS,
        HAPLOTYPE_CARRIER_PREDICTION_FIELDS};
    uint32_t call_id = leaf->carriers.first_call;
    for (unsigned list = 0u; list < 3u; list++) {
        duckvep_h_vector vector = v[HAPLOTYPE_PROVENANCE_COLUMN + list];
        duckvep_h_list_entry entry;
        if (!duckvep_h_list_extend(vector, counts[list], &entry)) return 0;
        ((duckvep_h_list_entry *)duckvep_h_data(vector))[row] = entry;
        duckvep_h_vector records = duckvep_h_list_child(vector), fields[HAPLOTYPE_PROVENANCE_FIELDS];
        duckvep_h_ensure_validity(records);
        for (unsigned j = 0u; j < field_counts[list]; j++) {
            fields[j] = duckvep_h_struct_child(records, j);
            duckvep_h_ensure_validity(fields[j]);
        }
        size_t block = 0u;
        duckvep_h_list_entry carrier_terms = {0u, 0u};
        duckvep_h_vector carrier_term_values = NULL;
        if (list == 2u && n) {
            /* One copy of the shared term list serves every predicted carrier of this leaf. */
            if (!duckvep_h_list_extend(fields[6], n, &carrier_terms)) return 0;
            carrier_term_values = duckvep_h_list_values(fields[6]);
            duckvep_h_ensure_validity(carrier_term_values);
            for (unsigned k = 0u; k < n; k++) {
                duckvep_h_set_valid(carrier_term_values, carrier_terms.offset + k);
                assign_text(carrier_term_values, carrier_terms.offset + k,
                    duckvep_so_name((duckvep_so_bit_t)bits[k]));
            }
        }
        for (size_t i = 0u; i < counts[list]; i++) {
            size_t at = entry.offset + i;
            duckvep_h_set_valid(records, at);
            for (unsigned j = 0u; j < field_counts[list]; j++)
                duckvep_h_set_valid(fields[j], at);
            if (!list) {
                const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
                ((uint64_t *)duckvep_h_data(fields[0]))[at] = c->source.event_id;
                if (c->source.alt_ordinal == UINT32_MAX) null_cell(fields[1], at);
                else ((uint32_t *)duckvep_h_data(fields[1]))[at] = c->source.alt_ordinal;
                ((uint32_t *)duckvep_h_data(fields[2]))[at] = c->source.chrom_id;
                ((uint64_t *)duckvep_h_data(fields[3]))[at] = c->source.pos1;
                duckvep_h_assign_string(fields[4], at, (const char *)c->source.ref, c->source.ref_len);
                duckvep_h_assign_string(fields[5], at, (const char *)c->source.alt, c->source.alt_len);
                ((uint8_t *)duckvep_h_data(fields[6]))[at] = c->evidence_flags;
                assign_text(fields[7], at, projection_name(c->projection_status));
                assign_text(fields[8], at, role_name(c->role));
                ((uint32_t *)duckvep_h_data(fields[9]))[at] = c->edit_count;
            } else if (list == 2u) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->stream.carriers, call_id);
                if (!call) return 0;
                duckvep_prediction_status_t cs;
                duckvep_prediction_reason_t cr;
                duckvep_haplotype_carrier_prediction(leaf, call, &cs, &cr);
                ((uint32_t *)duckvep_h_data(fields[0]))[at] = call->key.sample_index;
                ((int64_t *)duckvep_h_data(fields[1]))[at] = call->key.phase_set;
                if (!call->key.phase_set_present) null_cell(fields[1], at);
                ((uint16_t *)duckvep_h_data(fields[2]))[at] = call->key.lane;
                assign_text(fields[3], at, status_name(cs));
                assign_text(fields[4], at, reason_name(cr, leaf->prediction_projection));
                /* An empty list with NULL IMPACT is a lane equal to the reference; NULL both when the
                 * carrier itself is not predicted. */
                int carrier_decided = cs == DUCKVEP_PREDICTION_PREDICTED;
                assign_text(fields[7], at,
                    nmd_name(carrier_decided ? leaf->nmd : DUCKVEP_HAPLOTYPE_NMD_UNKNOWN));
                if (carrier_decided && leaf->nmd_stop_valid) ((uint64_t *)duckvep_h_data(fields[8]))[at] = leaf->nmd_stop_position1;
                else null_cell(fields[8], at);
                if (carrier_decided && leaf->nmd_junction_valid) ((uint64_t *)duckvep_h_data(fields[9]))[at] = leaf->nmd_junction_position1;
                else null_cell(fields[9], at);
                if (cs == DUCKVEP_PREDICTION_PREDICTED) {
                    if (n) assign_text(fields[5], at,
                        duckvep_impact_name(duckvep_so_impact(leaf->haplotype_so_mask)));
                    else null_cell(fields[5], at);
                    ((duckvep_h_list_entry *)duckvep_h_data(fields[6]))[at] = carrier_terms;
                } else {
                    null_cell(fields[5], at);
                    ((duckvep_h_list_entry *)duckvep_h_data(fields[6]))[at] = (duckvep_h_list_entry){0u, 0u};
                    null_cell(fields[6], at);
                }
                call_id = call->next_leaf;
            } else {
                const duckvep_haplotype_edit_t *e = &s->buffers.edits[i];
                ((uint64_t *)duckvep_h_data(fields[0]))[at] = i;
                ((uint64_t *)duckvep_h_data(fields[1]))[at] = s->buffers.edit_event_ids[i];
                while (leaf->cds && block < leaf->block_count &&
                       i >= leaf->blocks[block].edit_begin + leaf->blocks[block].edit_count) block++;
                if (leaf->cds && block < leaf->block_count && i >= leaf->blocks[block].edit_begin)
                    ((uint64_t *)duckvep_h_data(fields[2]))[at] = block;
                else null_cell(fields[2], at);
                ((uint32_t *)duckvep_h_data(fields[3]))[at] = e->cds_start;
                duckvep_h_assign_string(fields[4], at, e->ref_len ? (const char *)e->ref : "", e->ref_len);
                duckvep_h_assign_string(fields[5], at, e->alt_len ? (const char *)e->alt : "", e->alt_len);
                ((int8_t *)duckvep_h_data(fields[6]))[at] = e->variant_strand;
            }
        }
    }
    (void)bind;
    /* Whole-haplotype reduced SO set: NULL unless the same-codon classifier decided every carrier
     * of this row; an empty list (with NULL IMPACT) is a lane equal to the reference. */
    duckvep_h_vector consequences =
        duckvep_h_chunk_vector(output, HAPLOTYPE_CONSEQUENCES_COLUMN);
    duckvep_h_vector impact = duckvep_h_chunk_vector(output, HAPLOTYPE_IMPACT_COLUMN);
    duckvep_h_list_entry entry = {0u, 0u};
    int decided = leaf->prediction_status == DUCKVEP_PREDICTION_PREDICTED;
    if (!duckvep_h_list_extend(consequences, decided ? n : 0u, &entry)) return 0;
    ((duckvep_h_list_entry *)duckvep_h_data(consequences))[row] = entry;
    if (!decided) null_cell(consequences, row);
    else {
        duckvep_h_set_valid(consequences, row);
        duckvep_h_vector terms = duckvep_h_list_values(consequences);
        duckvep_h_ensure_validity(terms);
        for (unsigned i = 0u; i < n; i++) {
            duckvep_h_set_valid(terms, entry.offset + i);
            assign_text(terms, entry.offset + i, duckvep_so_name((duckvep_so_bit_t)bits[i]));
        }
    }
    if (decided && n) {
        duckvep_h_set_valid(impact, row);
        assign_text(impact, row, duckvep_impact_name(duckvep_so_impact(leaf->haplotype_so_mask)));
    } else null_cell(impact, row);
    /* ejc50-v1 row summary, gated like the consequences: unknown unless every carrier is decided.
     * The contributors that put the stop there are the APPLIED ones (the edits up to and including the stop),
     * by event_index; post_stop, shadowed and omitted sources are not listed. */
    assign_text(duckvep_h_chunk_vector(output, HAPLOTYPE_NMD_RULE_COLUMN), row,
        DUCKVEP_HAPLOTYPE_NMD_RULE);
    assign_text(duckvep_h_chunk_vector(output, HAPLOTYPE_NMD_COLUMN), row,
        nmd_name(decided ? leaf->nmd : DUCKVEP_HAPLOTYPE_NMD_UNKNOWN));
    duckvep_h_vector nmd_stop = duckvep_h_chunk_vector(output, HAPLOTYPE_NMD_STOP_COLUMN);
    duckvep_h_vector nmd_junction = duckvep_h_chunk_vector(output, HAPLOTYPE_NMD_JUNCTION_COLUMN);
    if (decided && leaf->nmd_stop_valid) ((uint64_t *)duckvep_h_data(nmd_stop))[row] = leaf->nmd_stop_position1;
    else null_cell(nmd_stop, row);
    if (decided && leaf->nmd_junction_valid) ((uint64_t *)duckvep_h_data(nmd_junction))[row] = leaf->nmd_junction_position1;
    else null_cell(nmd_junction, row);
    duckvep_h_vector attribution = duckvep_h_chunk_vector(output, HAPLOTYPE_NMD_CONTRIBUTORS_COLUMN);
    size_t attributed = 0u;
    if (decided && leaf->nmd_stop_valid)
        for (size_t i = 0u; i < leaf->contributor_count; i++) attributed += nmd_attributed(&leaf->contributors[i]);
    duckvep_h_list_entry attribution_entry;
    if (!duckvep_h_list_extend(attribution, attributed, &attribution_entry)) return 0;
    ((duckvep_h_list_entry *)duckvep_h_data(attribution))[row] = attribution_entry;
    if (!(decided && leaf->nmd_stop_valid)) null_cell(attribution, row);
    else {
        duckvep_h_vector values = duckvep_h_list_values(attribution);
        duckvep_h_ensure_validity(values);
        size_t at = 0u;
        for (size_t i = 0u; i < leaf->contributor_count; i++) {
            if (!nmd_attributed(&leaf->contributors[i])) continue;
            duckvep_h_set_valid(values, attribution_entry.offset + at);
            ((uint64_t *)duckvep_h_data(values))[attribution_entry.offset + at++] = leaf->contributors[i].source.event_id;
        }
    }
    return 1;
}

static int append_leaf(duckvep_h_chunk output, size_t row, duckvep_hap_state_t *s,
    const duckvep_hap_config_t *bind, const duckvep_haplotype_leaf_t *leaf,
    char *error, size_t error_size) {
    if (!prepare_difference_reference(s, bind, leaf, error, error_size)) return 0;
    duckvep_coding_context_t coding;
    uint32_t tx = leaf->carriers.transcript_index;
    if (!leaf->ordered_replacements && (leaf->block_count ||
            (bind->hgvs && leaf->cds && leaf->reference_protein &&
             leaf->reference_coding_translation.length))) {
        const duckvep_owned_model_t *model = bind->model;
        duckvep_edit_set_t edits = {s->buffers.edits, leaf->edit_count};
        duckvep_haplotype_result_t applied = {leaf->cds_length,
            leaf->nominal_length_diff,
            leaf->flags & ~(uint32_t)DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED, leaf->edit_count};
        duckvep_codon_table_t table = model->sequences.codon_table
            ? (duckvep_codon_table_t)model->sequences.codon_table[tx] : DUCKVEP_CODON_TABLE_STANDARD;
        const duckvep_event_t *event = NULL;
        if (leaf->edit_count == 1u) {
            for (size_t i = 0u; i < leaf->contributor_count; i++)
                if (leaf->contributors[i].source.event_id == leaf->edit_event_ids[0])
                    event = leaf->contributors[i].prepared;
        }
        if (duckvep_coding_context_open_replay(leaf->reference_cds, model->sequences.cds_length[tx],
                &edits, model->transcripts.strand[tx], table, leaf->cds, &applied,
                leaf->reference_coding_protein, &leaf->reference_coding_translation,
                s->buffers.protein, &leaf->translation, &coding) != DUCKVEP_CODING_CONTEXT_OK ||
            duckvep_coding_context_attach_model(&model->transcripts, &model->exons, &model->sequences,
                tx, event, leaf->edit_count == 1u ? edits.edits[0].cds_start : 0u, &coding) !=
                DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
            duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid completed coding context");
            return 0;
        }
    }
    duckvep_h_vector v[HAPLOTYPE_OUTPUT_COLUMNS];
    for (unsigned i = 0u; i < HAPLOTYPE_OUTPUT_COLUMNS; i++) {
        v[i] = duckvep_h_chunk_vector(output, i);
        duckvep_h_set_valid(v[i], row);
    }
    ((uint32_t *)duckvep_h_data(v[0]))[row] = leaf->carriers.transcript_index;
    if (leaf->cds) duckvep_h_assign_string(v[1], row, (const char *)leaf->cds, leaf->cds_length);
    else null_cell(v[1], row);
    if (leaf->protein) duckvep_h_assign_string(v[2], row, (const char *)leaf->protein, leaf->protein_length);
    else null_cell(v[2], row);
    ((uint32_t *)duckvep_h_data(v[3]))[row] = leaf->flags;
    if (leaf->cds)
        ((int64_t *)duckvep_h_data(v[HAPLOTYPE_NOMINAL_LENGTH_COLUMN]))[row] =
            leaf->nominal_length_diff;
    else null_cell(v[HAPLOTYPE_NOMINAL_LENGTH_COLUMN], row);
    ((uint8_t *)duckvep_h_data(v[4]))[row] = leaf->evidence_flags;
    assign_text(v[5], row, projection_name(leaf->projection_status));
    assign_text(v[6], row, leaf->projection_status == DUCKVEP_CDS_EDIT_OK
        ? sequence_name(leaf->sequence_status) : "unavailable_projection");
    ((uint64_t *)duckvep_h_data(v[7]))[row] = leaf->edit_count;
    ((uint32_t *)duckvep_h_data(v[8]))[row] = leaf->carriers.call_count;
    if (leaf->cds && !leaf->ordered_replacements)
        ((bool *)duckvep_h_data(v[HAPLOTYPE_STOP_COLUMN]))[row] = leaf->stop_in_displaced_frame != 0u;
    else null_cell(v[HAPLOTYPE_STOP_COLUMN], row);
    const size_t counts[] = {leaf->carriers.call_count, leaf->contributor_count, leaf->block_count};
    const unsigned field_counts[] = {4u, bind->source_records ? 8u : 7u, HAPLOTYPE_BLOCK_FIELDS};
    for (unsigned list = 0u; list < 3u; list++) {
        duckvep_h_vector vector = v[HAPLOTYPE_LIST_COLUMN + list];
        size_t count = counts[list];
        duckvep_h_list_entry entry;
        if (!duckvep_h_list_extend(vector, count, &entry)) return 0;
        size_t base = entry.offset;
        ((duckvep_h_list_entry *)duckvep_h_data(vector))[row] = entry;
        if (list >= 2u && !leaf->cds) null_cell(vector, row);
        duckvep_h_vector records = duckvep_h_list_child(vector), fields[HAPLOTYPE_BLOCK_FIELDS];
        duckvep_h_ensure_validity(records);
        for (unsigned j = 0u; j < field_counts[list]; j++) {
            fields[j] = duckvep_h_struct_child(records, j);
            duckvep_h_ensure_validity(fields[j]);
        }
        size_t event_base = 0u;
        if (list == 2u && leaf->cds) {
            duckvep_h_vector event_vector = fields[HAPLOTYPE_BLOCK_EVENT_FIELD];
            duckvep_h_list_entry events;
            if (!duckvep_h_list_extend(event_vector, leaf->edit_count, &events)) return 0;
            event_base = events.offset;
            duckvep_h_vector ids = duckvep_h_list_values(event_vector);
            duckvep_h_ensure_validity(ids);
            uint64_t *data = duckvep_h_data(ids);
            for (size_t i = 0u; i < leaf->edit_count; i++) {
                data[event_base + i] = leaf->edit_event_ids[i];
                duckvep_h_set_valid(ids, event_base + i);
            }
        }
        uint32_t call_id = leaf->carriers.first_call;
        for (size_t i = 0u; i < count; i++) {
            size_t at = base + i;
            duckvep_h_set_valid(records, at);
            for (unsigned j = 0u; j < field_counts[list]; j++)
                duckvep_h_set_valid(fields[j], at);
            if (!list) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&s->stream.carriers, call_id);
                if (!call) return 0;
                ((uint32_t *)duckvep_h_data(fields[0]))[at] = call->key.sample_index;
                ((int64_t *)duckvep_h_data(fields[1]))[at] = call->key.phase_set;
                if (!call->key.phase_set_present) null_cell(fields[1], at);
                ((uint16_t *)duckvep_h_data(fields[2]))[at] = call->key.lane;
                ((uint16_t *)duckvep_h_data(fields[3]))[at] = call->key.ploidy;
                call_id = call->next_leaf;
            } else if (list == 1u) {
                const duckvep_haplotype_contributor_t *c = &leaf->contributors[i];
                ((uint64_t *)duckvep_h_data(fields[0]))[at] = c->source.event_id;
                ((uint32_t *)duckvep_h_data(fields[1]))[at] = c->source.chrom_id;
                ((uint64_t *)duckvep_h_data(fields[2]))[at] = c->source.pos1;
                duckvep_h_assign_string(fields[3], at, (const char *)c->source.ref, c->source.ref_len);
                duckvep_h_assign_string(fields[4], at, (const char *)c->source.alt, c->source.alt_len);
                ((uint8_t *)duckvep_h_data(fields[5]))[at] = c->evidence_flags;
                assign_text(fields[6], at, projection_name(c->projection_status));
                if (bind->source_records) {
                    if (c->source.allele_index == UINT32_MAX) null_cell(fields[7], at);
                    else ((uint32_t *)duckvep_h_data(fields[7]))[at] = c->source.allele_index;
                }
            } else if (list == 2u) {
                const duckvep_haplotype_block_t *block = &leaf->blocks[i];
                ((uint32_t *)duckvep_h_data(fields[0]))[at] = block->cds_start;
                duckvep_h_assign_string(fields[1], at,
                    (const char *)leaf->reference_cds + block->cds_start - 1u, block->ref_len);
                duckvep_h_assign_string(fields[2], at,
                    (const char *)leaf->cds + block->alt_start0, block->alt_len);
                ((uint64_t *)duckvep_h_data(fields[3]))[at] = block->alt_start0;
                ((int64_t *)duckvep_h_data(fields[4]))[at] = block->length_diff;
                ((uint32_t *)duckvep_h_data(fields[5]))[at] = block->flags;
                duckvep_sequence_delta_t delta;
                duckvep_context_delta_status_t status = leaf->ordered_replacements
                    ? DUCKVEP_CONTEXT_DELTA_UNSUPPORTED : duckvep_coding_context_block_delta_fill(
                    &coding, s->buffers.edits, leaf->edit_count, block,
                    bind->model->transcripts.flags[tx], &delta);
                const char *name = leaf->ordered_replacements ? "unsupported_ordered_replacements" :
                    status == DUCKVEP_CONTEXT_DELTA_OK ? "ok" :
                    status == DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL ? "missing_transcript_tail" :
                    status == DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK ? "missing_transcript_flank" :
                    status == DUCKVEP_CONTEXT_DELTA_UNSUPPORTED ? "unsupported" : "invalid_argument";
                assign_text(fields[6], at, name);
                if (status == DUCKVEP_CONTEXT_DELTA_OK)
                    ((uint64_t *)duckvep_h_data(fields[7]))[at] = duckvep_effect_eval_coding_delta(&delta);
                else null_cell(fields[7], at);
                ((bool *)duckvep_h_data(fields[8]))[at] =
                    leaf->translation.first_stop_position1 &&
                    block->alt_start0 / 3u >= leaf->translation.first_stop_position1;
                ((duckvep_h_list_entry *)duckvep_h_data(fields[HAPLOTYPE_BLOCK_EVENT_FIELD]))[at] =
                    (duckvep_h_list_entry){event_base + block->edit_begin, block->edit_count};
            }
        }
    }
    if (!append_prediction(output, row, s, bind, leaf)) return 0;
    return append_sequence_differences(v[12], row, s, bind, leaf, 0, error, error_size) &&
        append_sequence_differences(v[13], row, s, bind, leaf, 1, error, error_size) &&
        append_hgvsp(v[HAPLOTYPE_HGVSP_COLUMN], v[HAPLOTYPE_HGVSP_STATUS_COLUMN], row,
            s, bind, leaf, &coding, error, error_size);
}

static int bit_valid(const uint64_t *validity, size_t index) {
    return !validity || (validity[index >> 6] >> (index & 63u) & 1u);
}

static duckvep_haplotype_stream_status_t consume_call(duckvep_hap_state_t *s, char *error, size_t error_size) {
    const duckvep_hap_config_t *bind = &s->bind;
    const duckvep_hap_row_t *r = &s->row;
    for (unsigned i = 0u; i < (bind->source_records ? 17u : 15u); i++) {
        if (i != 9u && i != 10u && (r->null_mask >> i & 1u)) {
            snprintf(error, error_size, "duckvep_haplotypes: required input column %u is NULL", i + 1u);
            return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        }
    }
    if (r->copies != 1 || r->versions != 1 || r->ploidies != 1) {
        duckvep_sql_set_error(error, error_size, bind->source_records
            ? "duckvep_haplotypes: duplicate call, inconsistent source record identity or source GT"
            : "duckvep_haplotypes: duplicate call, inconsistent event identity or changed sample/transcript ploidy");
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    }
    uint32_t chrom = r->chrom, allele_index = r->allele_index;
    uint64_t pos = r->pos;
    uint32_t ref_len = r->ref_len, alt_len = r->alt_len;
    if (chrom > UINT16_MAX || !pos || pos > UINT32_MAX || !ref_len || ref_len > UINT16_MAX ||
        (!alt_len && !(bind->source_records && allele_index == UINT32_MAX)) || alt_len > UINT16_MAX)
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    duckvep_haplotype_source_t source = {r->event_id, r->ref, r->alt,
        (uint32_t)pos, (uint16_t)chrom, (uint16_t)ref_len, (uint16_t)alt_len,
        bind->source_records ? allele_index : 0u, (uint8_t)bind->source_records,
        bind->source_records ? r->replay_order : 0u, allele_index};
    uint32_t tx = r->transcript;
    duckvep_haplotype_stream_status_t status;
    int new_event = !s->stream.have_input || source.event_id != s->stream.last_event_id ||
        source.allele_index != s->stream.last_allele_index;
    if (new_event) {
        status = duckvep_haplotype_stream_begin(&s->stream, &source);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
        s->have_call = 0;
    }
    if (!s->have_call || tx != s->last_tx) {
        status = duckvep_haplotype_stream_project(&s->stream, tx);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) return status;
    }
    if (bind->source_records) {
        duckvep_raw_gt_status_t parsed = (duckvep_raw_gt_status_t)r->raw[0];
        duckvep_raw_gt_t call = {{r->raw[1], r->raw[2]}, r->raw[3], (uint16_t)r->raw[4],
            (uint8_t)r->raw[5], (duckvep_raw_gt_disposition_t)r->raw[6]};
        if (parsed != DUCKVEP_RAW_GT_OK) {
            snprintf(error, error_size, "duckvep_haplotypes: raw GT status %u at event %llu, sample %u",
                (unsigned)parsed, (unsigned long long)source.event_id, r->sample);
            return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        }
        if (call.source_ploidy > bind->limits[DUCKVEP_HAP_LIMIT_PLOIDY] || bind->limits[DUCKVEP_HAP_LIMIT_PLOIDY] < 2u) {
            snprintf(error, error_size, "duckvep_haplotypes: max_ploidy=%zu exceeded by raw source/file ploidy",
                bind->limits[DUCKVEP_HAP_LIMIT_PLOIDY]);
            return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        }
        status = duckvep_haplotype_stream_push_raw_call(&s->stream, tx, r->sample, &call,
            (uint8_t)(r->source_selected != 0));
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            s->have_call = 1; s->last_tx = tx; s->have_row = 0;
        }
        return status;
    }
    size_t gt_length = r->gt_length;
    int have_phase = r->have_phase;
    if (!gt_length || gt_length > bind->limits[DUCKVEP_HAP_LIMIT_PLOIDY] ||
        (have_phase && r->phase_length != gt_length) || r->sets_length > bind->limits[DUCKVEP_HAP_LIMIT_PHASE_SETS]) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid GT/phase lengths or max_ploidy/max_phase_sets exceeded");
        return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
    }
    for (size_t i = 0u; i < gt_length; i++) {
        int missing = !bit_valid(r->gt_validity, r->gt_offset + i);
        s->gt[i] = missing ? -1 : r->gt[r->gt_offset + i];
        if (!missing && s->gt[i] < 0) return DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
        s->phase[i] = have_phase && bit_valid(r->phase_validity, r->phase_offset + i) &&
            r->phase[r->phase_offset + i];
    }
    for (size_t i = 0u; i < r->sets_length; i++) {
        s->sets[i].present = bit_valid(r->sets_validity, r->sets_offset + i);
        s->sets[i].value = s->sets[i].present ? r->sets[r->sets_offset + i] : 0;
    }
    duckvep_haplotype_call_t call = {s->gt, s->phase, r->sample, r->allele_index,
        (uint16_t)gt_length, {0, 0u}, bind->policy};
    call.phase_set.present = r->phase_set_present;
    if (call.phase_set.present) call.phase_set.value = r->phase_set;
    status = duckvep_haplotype_stream_push_call(&s->stream, tx, &call, s->sets, r->sets_length);
    if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        s->have_call = 1; s->last_tx = tx;
        s->have_row = 0;
    }
    return status;
}

size_t duckvep_hap_workspace_bytes(const duckvep_hap_state_t *state) {
    return state->workspace_bytes;
}

duckvep_hap_state_t *duckvep_hap_open(const duckvep_hap_config_t *config, char *error, size_t error_size) {
    duckvep_hap_state_t *s = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1u, sizeof(*s));
    duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: workspace allocation or configured limit exceeded");
    duckvep_budget_clear_failure();
    if (!s) return NULL;
    s->bind = *config;
    if (!workspace_allocate(s, &s->bind)) goto failed;
    const duckvep_owned_model_t *m = s->bind.model;
    if (s->bind.hgvs && !duckvep_reference_reader_init(&s->reference, m,
            s->reference.bases, s->reference.capacity, error, error_size)) goto failed;
    if (duckvep_haplotype_stream_init(&s->stream, &m->transcripts, &m->exons, &m->sequences,
        &s->buffers) != DUCKVEP_HAPLOTYPE_STREAM_OK) {
        duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: invalid native model/workspace"); goto failed;
    }
    return s;
failed:
    duckvep_hap_close(s);
    return NULL;
}

int duckvep_hap_scan(duckvep_hap_state_t *s, const duckvep_hap_input_t *input, void *chunk, size_t capacity,
    size_t *rows_out, char *error, size_t error_size) {
    const duckvep_hap_config_t *bind = &s->bind;
    duckvep_h_chunk output = chunk;
    size_t rows = 0u;
    error[0] = '\0';
    while (rows < capacity) {
        duckvep_haplotype_stream_status_t status;
        if (s->stream.closing) {
            duckvep_haplotype_leaf_t leaf;
            status = duckvep_haplotype_stream_next(&s->stream, &leaf);
            if (status == DUCKVEP_HAPLOTYPE_STREAM_DONE) continue;
            if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) {
                if (!append_leaf(output, rows, s, bind, &leaf, error, error_size)) {
                    if (!error[0]) duckvep_sql_set_error(error, error_size,
                        "duckvep_haplotypes: output list allocation failed");
                    return 0;
                }
                rows++; continue;
            }
        } else {
            if (!s->eof && !s->have_row) {
                int fetched = input->next(input->context, &s->row, error, error_size);
                if (fetched < 0) {
                    if (!error[0]) duckvep_sql_set_error(error, error_size, "duckvep_haplotypes: could not read the input");
                    return 0;
                }
                if (fetched == 0) s->eof = 1; else s->have_row = 1;
            }
            status = s->eof ? duckvep_haplotype_stream_finish(&s->stream) : consume_call(s, error, error_size);
            if (status == DUCKVEP_HAPLOTYPE_STREAM_DONE) break;
            if (status == DUCKVEP_HAPLOTYPE_STREAM_OK || status == DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) continue;
        }
        if (!error[0]) {
            uint64_t event = s->stream.last_event_id, position = s->stream.last_pos1;
            if (!s->eof && s->have_row) { event = s->row.event_id; position = s->row.pos; }
            int limit = exhausted_limit(status, s->stream.carrier_error);
            if (limit >= 0) snprintf(error, error_size,
                "duckvep_haplotypes: %s=%zu exhausted at event %llu, position %llu",
                duckvep_hap_limit_names[limit], bind->limits[limit], (unsigned long long)event, (unsigned long long)position);
            else snprintf(error, error_size,
                "duckvep_haplotypes: native status %u, carrier status %u at event %llu, position %llu; invalid input or candidate",
                (unsigned)status, (unsigned)s->stream.carrier_error,
                (unsigned long long)event, (unsigned long long)position);
        }
        return 0;
    }
    *rows_out = rows;
    return 1;
}
