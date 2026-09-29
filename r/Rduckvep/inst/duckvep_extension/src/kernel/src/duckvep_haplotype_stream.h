/* Model-scoped literal-event replay over sparse carrier paths (INTERNAL).
 * All storage is caller-owned. begin copies REF/ALT once; project adds one
 * candidate transcript; push_call interprets a decoded GT and push_raw_call
 * interprets retained source-record slots. Candidate and phase-domain discovery and output
 * materialization belong to the host query plan.
 *
 * Input is sorted by (chrom_id, pos1, event_id, allele_index). begin may report a transcript
 * ready before consuming the input: drain next until DONE, then retry begin with
 * the same input. finish uses the same drain protocol. Each next result borrows
 * scratch until the next next/begin/finish call; its carrier list stays valid
 * throughout that transcript's drain. The owner may pause between any two calls.
 *
 * Event/projection/allele rings retain only the oldest still-active genomic
 * window, including intervening events pinned behind a longer transcript.
 * Limits are explicit; exhaustion is an error, never silent loss or growth.
 * Errors latch until reinitialization. Per-leaf projection/edit conflicts are
 * results with complete provenance, not errors that drop an occupied path.
 */
#ifndef DUCKVEP_HAPLOTYPE_STREAM_H
#define DUCKVEP_HAPLOTYPE_STREAM_H

#include "duckvep_carriers.h"
#include "duckvep_delta.h"
#include "duckvep_phase.h"

/* Pinned Haplosaurus input-buffer planning over a complete coordinate-sorted
 * source stream and borrowed, start-sorted transcript spans. The first record
 * that overlaps a transcript closes over all connected transcript spans. Later
 * records enter that buffer by start coordinate; their ends do not extend it.
 * Records outside every buffer return buffer=ordinal=0. No allocation. */
typedef struct {
    const duckvep_transcript_model_t *model;
    size_t transcript;
    uint64_t buffer, ordinal;
    uint32_t end1, last_pos1;
    uint16_t chrom;
    uint8_t have_input;
} duckvep_haplotype_record_plan_t;

int duckvep_haplotype_record_plan_init(duckvep_haplotype_record_plan_t *plan,
    const duckvep_transcript_model_t *model);
int duckvep_haplotype_record_plan_next(duckvep_haplotype_record_plan_t *plan,
    uint16_t chrom, uint32_t start1, uint32_t end1, uint64_t *buffer, uint64_t *ordinal);
/* One-based traversal rank after sorted insertion into Set::IntervalTree 0.12.
 * Counts every record in the source buffer, including unretained REF calls.
 * Zero indicates invalid count/ordinal. Constant space and logarithmic time. */
uint64_t duckvep_haplotype_record_order(uint64_t count, uint64_t ordinal);

typedef enum {
    DUCKVEP_HAPLOTYPE_STREAM_OK,
    DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
    DUCKVEP_HAPLOTYPE_STREAM_DONE,
    DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG,
    DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER,
    DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR,
    DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR
} duckvep_haplotype_stream_status_t;

typedef struct {
    uint64_t event_id;
    const uint8_t *ref, *alt;
    uint32_t pos1;
    uint16_t chrom_id, ref_len, alt_len;
    /* With source_record set, event_id identifies a source record and this is
     * its REF (zero), positive ALT ordinal, or UINT32_MAX for the undefined
     * file-slot interpretation (empty ALT, complete source REF deletion). */
    uint32_t allele_index;
    uint8_t source_record;
    uint64_t replay_order; /* Unique positive planned source rank; all-zero uses caller order. */
    /* Original ALT ordinal supplied by the caller (0 = REF, UINT32_MAX = undefined
     * file slot). Provenance only: it never affects ordering or identity. */
    uint32_t alt_ordinal;
} duckvep_haplotype_source_t;

/* Versioned coding-v1 status of one completed leaf. ELIGIBLE means the leaf is inside
 * the supported domain but its whole-haplotype consequence needs a classifier that has
 * not landed (no path is left in this state since the start/stop classifier); PREDICTED means the
 * same-codon, frame/stop-gain or start/stop classifier decided it and leaf.haplotype_so_mask is authoritative.
 * Failures keep every contributor. */
typedef enum {
    DUCKVEP_PREDICTION_ELIGIBLE = 0,
    DUCKVEP_PREDICTION_INCOMPLETE_INPUT,
    DUCKVEP_PREDICTION_EDIT_CONFLICT,
    DUCKVEP_PREDICTION_UNSUPPORTED_OVERLAP,
    DUCKVEP_PREDICTION_UNSUPPORTED_CONTEXT,
    DUCKVEP_PREDICTION_PREDICTED
} duckvep_prediction_status_t;

typedef enum {
    DUCKVEP_REASON_SUPPORTED_DOMAIN = 0,
    DUCKVEP_REASON_MISSING_CALL,
    DUCKVEP_REASON_UNPHASED_HETEROZYGOUS,
    DUCKVEP_REASON_CROSS_PS_UNRESOLVED,
    DUCKVEP_REASON_CONTRADICTORY_EDITS,
    DUCKVEP_REASON_OVERLAPPING_EDITS,
    DUCKVEP_REASON_DUPLICATE_EDITS,
    DUCKVEP_REASON_SAME_GAP_INSERTIONS,
    DUCKVEP_REASON_NON_STRICT_PHASE_POLICY,
    DUCKVEP_REASON_NON_DIPLOID_CALL,
    DUCKVEP_REASON_PROJECTION, /* leaf.prediction_projection names the preserved projection status */
    DUCKVEP_REASON_REFERENCE_MISMATCH,
    DUCKVEP_REASON_INVALID_BASE,
    DUCKVEP_REASON_INVALID_SEQUENCE,
    DUCKVEP_REASON_TRANSCRIPT_NOT_CODING,
    DUCKVEP_REASON_NON_STANDARD_CODON_TABLE,
    DUCKVEP_REASON_CURATED_TRANSCRIPT,
    DUCKVEP_REASON_INCOMPLETE_CDS,
    DUCKVEP_REASON_NONCANONICAL_START,
    DUCKVEP_REASON_NONCANONICAL_STOP,
    DUCKVEP_REASON_INTERNAL_STOP,
    DUCKVEP_REASON_NON_LITERAL_ALLELE,
    DUCKVEP_REASON_ALLELE_OVER_50,
    /* No longer emitted since slice 5 (every eligible path is decided); kept so the numbering is stable. */
    DUCKVEP_REASON_START_STOP_CLASSIFIER_PENDING
} duckvep_prediction_reason_t;

/* Whole-haplotype NMD prediction under rule ejc50-v1 (coding-v1 slice 6). An EJC-distance heuristic on the
 * edited spliced transcript only: no reinitiation, no long-exon exception, no NMD_transcript_variant biotype
 * term, and not the single-allele VEP NMD plugin. UNKNOWN (zero) is every path that is not PREDICTED, a
 * lost start, an edited CDS with no stop (termination unavailable) and incomplete exon topology;
 * NOT_APPLICABLE is a reference lane and known termination without a newly premature stop. */
#define DUCKVEP_HAPLOTYPE_NMD_RULE "ejc50-v1"
#define DUCKVEP_HAPLOTYPE_NMD_THRESHOLD 50
typedef enum {
    DUCKVEP_HAPLOTYPE_NMD_UNKNOWN = 0,
    DUCKVEP_HAPLOTYPE_NMD_NOT_APPLICABLE,
    DUCKVEP_HAPLOTYPE_NMD_ESCAPE,
    DUCKVEP_HAPLOTYPE_NMD_TRIGGER
} duckvep_haplotype_nmd_t;

/* Contributor role in the completed leaf. Every contributor is retained. */
typedef enum {
    DUCKVEP_ROLE_OMITTED = 0, /* No projected sequence edit (outside CDS, no called allele, no net change). */
    DUCKVEP_ROLE_SHADOWED,    /* Selected away by another source at the same locus. */
    DUCKVEP_ROLE_UNAPPLIED,   /* Projected edits exist, but the leaf has no rebuilt sequence. */
    DUCKVEP_ROLE_APPLIED,
    DUCKVEP_ROLE_POST_STOP    /* Applied only after the first stop of the rebuilt protein. */
} duckvep_haplotype_role_t;

typedef struct {
    duckvep_haplotype_source_t source;
    duckvep_event_t prepared;
    uint64_t serial;
    size_t allele_consumed;
    uint32_t projection_begin, projection_count, last_end1;
} duckvep_haplotype_stored_event_t;

typedef struct {
    duckvep_haplotype_edit_t edit;
    uint32_t transcript_index;
    duckvep_cds_edit_status_t status;
    uint8_t cds_unaffected; /* Proven no literal CDS contribution, not an ignored projection failure. */
    uint8_t source_selected; /* Raw source selected for this candidate's replacement. */
    uint8_t selection_set;
    uint8_t source_exonic; /* Full raw REF span reaches an exon admitted by Haplosaurus. */
} duckvep_haplotype_projection_t;

typedef struct {
    /* Source bytes are contiguous REF then ALT in the retained allele ring;
     * neither span crosses its end. They remain valid throughout this drain. */
    duckvep_haplotype_source_t source;
    duckvep_cds_edit_status_t projection_status;
    uint8_t evidence_flags;
    const duckvep_event_t *prepared; /* Borrowed source geometry for this transcript drain. */
    const duckvep_haplotype_edit_t *projected; /* Successful physical projection, or NULL. */
    uint8_t source_replaced; /* Ordered replay changed the then-current sequence. */
    uint8_t role;            /* duckvep_haplotype_role_t. */
    uint32_t edit_count;     /* Differing edit islands contributed by this source. */
    uint32_t post_stop_edits; /* Of those, islands after the first stop of the rebuilt protein. */
} duckvep_haplotype_contributor_t;

typedef struct {
    int64_t value;
    uint8_t present;
} duckvep_haplotype_phase_set_t;

typedef struct {
    const int32_t *alleles;
    const uint8_t *phase_before; /* NULL means no phase information. */
    uint32_t sample_index, alt_index;
    uint16_t ploidy;
    duckvep_haplotype_phase_set_t phase_set;
    duckvep_phase_policy_t policy;
} duckvep_haplotype_call_t;

typedef struct {
    duckvep_carrier_buffers_t carriers;
    duckvep_haplotype_stored_event_t *events;
    duckvep_haplotype_projection_t *projections;
    uint8_t *alleles;
    uint32_t event_capacity, projection_capacity;
    size_t allele_capacity;
    /* Scratch for one distinct occupied path, reused across its carriers. */
    duckvep_carrier_event_t *leaf_events;
    duckvep_haplotype_contributor_t *contributors;
    duckvep_haplotype_edit_t *edits;
    uint64_t *edit_event_ids; /* Edit payloads, compacted to applied source IDs for ordered replay. */
    duckvep_haplotype_block_t *blocks; /* At most edit_capacity interaction blocks. */
    size_t leaf_capacity, edit_capacity;
    uint8_t *cds, *protein, *reference_protein, *reference_coding_protein;
    size_t cds_capacity, protein_capacity, reference_protein_capacity;
} duckvep_haplotype_stream_buffers_t;

typedef struct {
    duckvep_carrier_leaf_t carriers;
    const duckvep_haplotype_contributor_t *contributors;
    size_t contributor_count;
    size_t edit_count; /* Differing islands, or applied full spans when ordered_replacements is set. */
    /* One source identity per differing island or applied full-span operation,
     * in ascending reference CDS order (tied-source order is not guaranteed).
     * A source may occur more than once or in several blocks. Borrowed only for
     * known sequences; block.edit_begin/edit_count select the corresponding IDs. */
    const uint64_t *edit_event_ids;
    const duckvep_haplotype_block_t *blocks; /* Ascending reference CDS order. */
    size_t block_count;
    const uint8_t *reference_cds; /* Model-owned; block spans borrow this and cds. */
    /* Raw mutation protein is a first-stop prefix. Without a retained exonic
     * genotype, raw replay uses the curated reference peptide instead. */
    const uint8_t *cds, *protein;
    const uint8_t *reference_protein; /* Worker-owned; NULL without a complete reference codon. */
    size_t reference_protein_length;
    const uint8_t *reference_coding_protein; /* Full consensus view, without reference curation. */
    duckvep_translation_t reference_coding_translation;
    duckvep_translation_t translation; /* Full raw translation remains in buffers.protein. */
    size_t cds_length, protein_length;
    /* Sum of replayed ALT lengths minus nominal REF spans, before clipping or
     * equality checks. Valid when cds is present, including conditional replay. */
    int64_t nominal_length_diff;
    uint32_t flags;
    uint8_t evidence_flags; /* OR of contributor evidence, distinct from sequence flags. */
    uint8_t stop_in_displaced_frame; /* First translated stop intersects a frame excursion. */
    /* Source-order components have exact sequences/provenance, but do not provide
     * a disjoint physical edit history for local SO or frame-span interpretation. */
    uint8_t ordered_replacements;
    /* First failed coding projection, or edit/rebuild status when projection is OK.
     * Proven noncoding contributors keep their own OUT_OF_CDS status but do not
     * suppress a coding transcript's literal CDS. This does not predict splicing.
     * Failed paths have no CDS/protein/blocks; all contributors/carriers remain.
     * CONDITIONAL has sequence with explicitly interpreted source evidence;
     * INPUT_INCOMPLETE has no sequence. Projection failures override both. */
    duckvep_cds_edit_status_t projection_status;
    duckvep_haplotype_status_t sequence_status;
    /* Coding-v1 eligibility; see duckvep_prediction_status_t. The path fields depend only
     * on the shared edit path. The row summary is ELIGIBLE only when every carrier is;
     * otherwise it is the first ineligible carrier's result. Phase-domain completeness
     * and ploidy belong to each sample/phase/lane key: use
     * duckvep_haplotype_carrier_prediction for the keyed result. */
    duckvep_prediction_status_t prediction_status, path_status;
    duckvep_prediction_reason_t prediction_reason, path_reason;
    duckvep_cds_edit_status_t prediction_projection; /* Valid for DUCKVEP_REASON_PROJECTION. */
    /* Same-codon and frame/stop-gain classifier (coding-v1 slices 3 and 4). Set only when path_status is PREDICTED;
     * the reduced whole-protein SO set as DUCKVEP_SO() bits, empty for a lane equal to the
     * reference (which has no IMPACT). */
    uint64_t haplotype_so_mask;
    /* Slice 6: ejc50-v1. S is the final nucleotide of the first stop codon and J the final nucleotide of the
     * penultimate exon, both 1-based in edited spliced-transcript (cDNA, 5' UTR included) coordinates. Each is
     * valid only when its flag is set: S for a newly premature stop, J additionally for a multi-exon
     * transcript. The contributors that put the stop there are those with role APPLIED. */
    duckvep_haplotype_nmd_t nmd;
    uint8_t nmd_stop_valid, nmd_junction_valid;
    uint64_t nmd_stop_position1, nmd_junction_position1;
    /* Ascending-CDS edit islands with source IDs, in stream buffers edits/edit_event_ids.
     * Listed for known sequences and for failed decoded-call leaves, so conflicts and
     * omitted or post-stop sources stay attributable. */
    size_t listed_edit_count;
} duckvep_haplotype_leaf_t;

/* Keyed coding-v1 result of one carrier of a completed leaf: unresolved cross-PS phase
 * and non-diploid calls are properties of the sample/phase/lane, not of the shared path. */
void duckvep_haplotype_carrier_prediction(const duckvep_haplotype_leaf_t *leaf,
    const duckvep_carrier_call_t *call, duckvep_prediction_status_t *status,
    duckvep_prediction_reason_t *reason);

typedef struct {
    duckvep_carriers_t carriers;
    const duckvep_exon_model_t *exons;
    const duckvep_sequence_pool_t *sequences;
    duckvep_haplotype_stream_buffers_t buffers;
    uint32_t event_begin, event_count, projection_begin, projection_count;
    uint32_t current_event, closing;
    size_t allele_begin, allele_count;
    uint64_t serial, last_event_id;
    uint32_t last_pos1, last_allele_index;
    uint16_t last_chrom;
    uint8_t have_input, have_current, initialized;
    uint8_t have_phase_policy;
    duckvep_phase_policy_t phase_policy;
    uint32_t reference_transcript;
    size_t reference_protein_length;
    duckvep_translation_t reference_coding_translation;
    uint8_t have_reference_protein, reference_protein_known;
    uint8_t have_domain;
    uint32_t domain_transcript;
    duckvep_prediction_reason_t domain_reason;
    duckvep_haplotype_stream_status_t error;
    duckvep_carriers_status_t carrier_error;
    uint64_t input_events, projected_events, completed_leaves, translated_bases;
    uint32_t peak_events, peak_projections;
    size_t peak_alleles;
} duckvep_haplotype_stream_t;

/* Distinct buffers; immutable model views remain pinned until teardown. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_init(
    duckvep_haplotype_stream_t *stream,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t *exons,
    const duckvep_sequence_pool_t *sequences,
    const duckvep_haplotype_stream_buffers_t *buffers);

/* The owner may discard input alleles after OK; input must not alias workspace
 * storage. Event IDs need only
 * be unique, with increasing IDs used to order events at the same coordinate.
 * Interpretations of one source record share its ID and geometry and arrive in
 * increasing allele_index order; UINT32_MAX is last, never a real ALT. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_begin(
    duckvep_haplotype_stream_t *stream, const duckvep_haplotype_source_t *event);

/* Add candidates in strictly increasing model-local ordinal order. Calls for
 * an earlier candidate may be consumed before projecting the next candidate;
 * an event's entire cohort need not be buffered. Each candidate must overlap
 * the genomic span (an insertion may touch the transcript end). */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_project(
    duckvep_haplotype_stream_t *stream, uint32_t transcript_index);

/* Interpret the complete decoded GT for this source ALT and one candidate
 * transcript. The caller's query plan supplies ALL phase sets for this sample
 * and transcript, including those first observed in later records. Arrays are
 * borrowed for the call, sorted absent-first then by signed value, and unique.
 * An empty domain means only the absent/default set. Domains must stay logically
 * identical across the transcript; discovery/broadcast planning belongs to the
 * host, not a second first-party store of the whole query's phase-set catalogue.
 * Homozygous/haploid and wholly unphased evidence is broadcast across the domain.
 * Partial phase affects only its declared set and the unresolved slots within it.
 * Compatibility mode requires the absent-only domain, compacts called slots,
 * and retains missing-call evidence on every lane. Uncertain paths return
 * INPUT_INCOMPLETE with no CDS/protein, including in compatibility mode; this
 * does not certify VEP's conditional sequence output for missing genotypes.
 * Source ALT ordinals must be positive and refer to the current begin event.
 * Missing GT with unknown ploidy is an error; known-ploidy missing slots use -1.
 * One phase policy applies to the stream; changing it is an error.
 * A failure latches even if a broadcast already updated an earlier carrier.
 */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push_call(
    duckvep_haplotype_stream_t *stream, uint32_t transcript_index,
    const duckvep_haplotype_call_t *call,
    const duckvep_haplotype_phase_set_t *phase_sets, size_t phase_set_count);

/* Route a validated raw-parser result to this source-record interpretation in
 * one candidate transcript. The host supplies every selected ALT interpretation
 * once, plus an empty-ALT interpretation if a retained call has an undefined
 * file slot. Retained REF slots participate in ordered replacements; omitted missing
 * calls retain conditional observations without replacing sequence. Nonmutating known
 * REF observations are omitted from emitted contributors. File ploidy is two and PS is
 * ignored; the input source ploidy remains in the parser result, not the key.
 * Missing-source, undefined-slot, SOURCE_UNMAPPED and SOURCE_ALLELE_SKIPPED
 * observations make replay explicitly CONDITIONAL. SOURCE_UNMAPPED retains a
 * layout/REF-checked mapper omission without executing a partial edit.
 * SOURCE_ALLELE_SKIPPED retains a checked raw ALT excluded from mutation,
 * with complete provenance but no edit or coding block; other valid sources
 * still replay. A checked REF slot excluded from mutation is a nonmutating
 * reference observation; its missing-call evidence, if present, is preserved.
 * An exon-admitted skipped source selects normal CDS translation.
 * Other projection errors withhold sequence. Mixing raw and decoded policies
 * is an error.
 * source_selected is the candidate-wide mapping choice across all samples;
 * zero retains contributor evidence without executing the source replacement.
 * The choice must agree for every sample at this interpretation/candidate. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push_raw_call(
    duckvep_haplotype_stream_t *stream, uint32_t transcript_index,
    uint32_t sample_index, const duckvep_raw_gt_t *call, uint8_t source_selected);

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_finish(
    duckvep_haplotype_stream_t *stream);

/* Each occupied edit prefix is rebuilt/translated once. DONE releases the
 * completed transcript; retry begin/finish to accept input or close the next.
 * Iterate leaf.carriers.first_call using duckvep_carriers_call(&stream->carriers,
 * id), following next_leaf, before calling next again. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_next(
    duckvep_haplotype_stream_t *stream, duckvep_haplotype_leaf_t *leaf);

#endif
