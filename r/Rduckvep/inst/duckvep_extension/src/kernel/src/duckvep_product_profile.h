#ifndef DUCKVEP_PRODUCT_PROFILE_H
#define DUCKVEP_PRODUCT_PROFILE_H

#include <stddef.h>
#include <stdint.h>
#include "duckvep_haplotype.h"

/* Positional geometry over borrowed equal-length genomic sequences. The
 * execute_position entry point composes the original sources before use.
 * This probe reports intron occupancy and feature intersection only, not
 * splice consequences or a biological product. Coordinates are one-based. */
typedef enum {
    DUCKVEP_PROFILE_OK = 0,
    DUCKVEP_PROFILE_INVALID,
    DUCKVEP_PROFILE_REFERENCE_MISMATCH,
    DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE,
    DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,
    DUCKVEP_PROFILE_PARTIAL_MODEL,
    DUCKVEP_PROFILE_SHORT_BUFFER,
    DUCKVEP_PROFILE_CAPACITY,
    DUCKVEP_PROFILE_UNSUPPORTED_PRODUCT
} duckvep_profile_status_t;

enum {
    DUCKVEP_PROFILE_INTRON = 1u,
    DUCKVEP_PROFILE_SPLICE_REGION = 2u,
    DUCKVEP_PROFILE_DONOR = 4u,
    DUCKVEP_PROFILE_ACCEPTOR = 8u,
    DUCKVEP_PROFILE_REGULATORY_INTERVAL = 16u
};

typedef struct { uint32_t start1, end1; } duckvep_profile_interval_t;
typedef struct {
    duckvep_profile_interval_t interval;
    uint64_t feature_id; /* Nonzero caller-owned source/provenance identifier. */
} duckvep_profile_feature_t;
typedef struct {
    uint32_t reference_start1;
    const uint8_t *reference_dna, *reference_witness, *final_dna;
    size_t dna_length, witness_length, final_length;
    duckvep_profile_status_t joint_status; /* Positional-only API input; ignored by source execution. */
    duckvep_profile_interval_t transcript;
    const duckvep_profile_interval_t *exons; /* Genomic ascending, complete model. */
    size_t exon_count;
    int8_t strand;
    uint8_t model_complete;
    const duckvep_profile_feature_t *features; /* Ascending (start,end,id). */
    size_t feature_count;
} duckvep_profile_v1_input_t;
typedef struct {
    uint32_t pos1, flags;
    uint64_t feature_id; /* Zero for positional row; nonzero for feature row. */
} duckvep_profile_v1_fact_t;

/* On failure count=0 and out is untouched. On success count includes the
 * positional row for every changed base and one row per overlapping feature.
 * Nonchanged bases have no rows: a successful empty result is a known zero.
 * INTRON means occupancy strictly between modeled exons. No splice motif
 * interpretation or VEP SO predicate is represented.
 */
duckvep_profile_status_t duckvep_product_position_v1(
    const duckvep_profile_v1_input_t *input, duckvep_profile_v1_fact_t *out,
    size_t capacity, size_t *count);

typedef struct {
    duckvep_profile_interval_t genomic_source;
    uint32_t rna_start1; /* Blocks ordered in transcript orientation, consecutive on RNA. */
} duckvep_profile_rna_block_v1_t;
typedef struct {
    const char *identity, *evidence_id, *assumption_id;
    const duckvep_profile_rna_block_v1_t *blocks;
    size_t block_count;
    const uint8_t *reference_rna; size_t rna_length;
    uint32_t cds_start1, cds_stop_end1; /* 1-based RNA, inclusive. */
    uint32_t annotated_start1, alternative_start1; /* 0 if not proposed. */
    const char *trigger_id, *counterfactual_id; /* Required for alternative start. */
} duckvep_profile_programme_v1_t;
/* Opt-in fixed-length source-to-product probe. Original source coordinates are
 * genomic, one-based; edits must be distinct, disjoint and in descending order.
 * IDs identify original sources, not inferred aliases. Workspace spans are
 * caller-owned, distinct from input and from each other. No persistence. */
typedef struct {
    const duckvep_profile_v1_input_t *geometry; /* reference_dna + exon topology */
    const duckvep_haplotype_edit_t *sources; /* cds_start is genomic start1 */
    const uint64_t *source_ids;
    size_t source_count;
    const duckvep_profile_programme_v1_t *programme;
    int32_t codon_table;
} duckvep_profile_execution_t;
typedef struct {
    uint64_t source_id;
    uint32_t genomic_start1, genomic_end1;
    uint32_t rna_first1, rna_last1; /* zero if this source does not alter the named RNA */
    uint32_t changed_rna_bases;
    uint8_t changed_dna;
} duckvep_profile_source_effect_t;
typedef struct {
    uint8_t *dna, *reference_rna, *alternate_rna, *reference_peptide, *alternate_peptide;
    size_t dna_capacity, rna_capacity, peptide_capacity;
    duckvep_profile_source_effect_t *source_effects;
    size_t source_effect_capacity;
    duckvep_haplotype_edit_t *edits;
    size_t edit_capacity;
} duckvep_profile_workspace_t;
typedef struct {
    const char *identity, *evidence_id, *assumption_id;
    size_t rna_length, reference_peptide_length, alternate_peptide_length;
    size_t reference_stop1, alternate_stop1;
    uint32_t reference_start1, alternate_start1;
    uint32_t last_junction1; /* RNA position of final upstream exon base; zero if absent. */
    uint8_t annotated_start_lost, alternative_start_used;
    size_t verified_sources; /* original REF-checked sources, including silent edits */
} duckvep_profile_execution_result_t;
/* Failure leaves result untouched; workspace is scratch and may be modified.
 * Exon blocks must be exact complete source exons or a declared included
 * intron block; strict increasing transcript orientation. Alias/overlapping
 * source representations are unsupported, not evidence of alias invariance. */
duckvep_profile_status_t duckvep_product_execute_v1(
    const duckvep_profile_execution_t *input,
    duckvep_profile_workspace_t *workspace,
    duckvep_profile_execution_result_t *result);
/* Joint execution with geometry-only positional facts. No VEP/SO prediction. */
duckvep_profile_status_t duckvep_product_execute_position_v1(
    const duckvep_profile_execution_t *input,
    duckvep_profile_workspace_t *workspace,
    duckvep_profile_execution_result_t *result,
    duckvep_profile_v1_fact_t *facts, size_t capacity, size_t *count);
#endif
