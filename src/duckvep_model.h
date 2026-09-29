/* Private resident-model registry shared by the DuckDB adapter functions. */
#ifndef DUCKVEP_MODEL_H
#define DUCKVEP_MODEL_H

#include "duckdb_extension.h"
#include "cgranges.h"
#include "kernel/include/duckvep_kernel.h"
#include "duckvep_reference.h"
#include "kernel/src/duckvep_lift.h"

#include <htslib/faidx.h>

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>

#define DUCKVEP_SQL_ERROR_SIZE 512

typedef struct duckvep_reference_file_identity {
	uint64_t device;
	uint64_t inode;
	uint64_t size;
	int64_t mtime_seconds;
	uint32_t mtime_nanoseconds;
	int present;
} duckvep_reference_file_identity_t;

struct duckvep_lifted_model;

typedef struct duckvep_owned_model {
	duckvep_transcript_model_t transcripts;
	duckvep_exon_model_t exons;
	duckvep_sequence_pool_t sequences;
	duckvep_interval_feature_model_t interval_features;
	duckvep_model_t *kernel;
	uint16_t *known_seq_regions;
	uint32_t *sequence_lengths;
	uint8_t *region_circular;
	int has_wrapped_coordinates;
	/* Non-NULL when a circular region carries wrapped objects. Its objects
	 * are represented by lifted images and every annotation entry point reads
	 * this one execution view; this struct keeps the source contract. */
	struct duckvep_lifted_model *lifted;
	char **sequence_names;
	char *reference_fasta_path;
	char *reference_fai_path;
	char *reference_gzi_path;
	char *reference_fasta_open_path;
	char *reference_fai_open_path;
	char *reference_gzi_open_path;
	duckvep_reference_file_identity_t reference_fasta_identity;
	duckvep_reference_file_identity_t reference_fai_identity;
	duckvep_reference_file_identity_t reference_gzi_identity;
	int reference_fasta_descriptor;
	int reference_fai_descriptor;
	int reference_gzi_descriptor;
	int reference_descriptors_open;
	size_t known_seq_region_count;
	uint16_t *seq_regions;
	uint32_t *transcript_starts;
	uint32_t *transcript_ends;
	int8_t *strands;
	uint64_t *transcript_flags;
	uint32_t *gene_indices;
	uint32_t *exon_offsets;
	uint16_t *exon_counts;
	uint32_t *cds_starts;
	uint32_t *cds_ends;
	uint64_t *cds_sequence_offsets;
	uint32_t *cds_sequence_lengths;
	uint8_t *codon_tables;
	uint64_t *pre_cds_sequence_offsets;
	uint32_t *pre_cds_sequence_lengths;
	uint64_t *post_cds_sequence_offsets;
	uint32_t *post_cds_sequence_lengths;
	uint32_t *exon_starts;
	uint32_t *exon_ends;
	uint32_t *exon_cdna_starts;
	uint32_t *exon_cdna_ends;
	int8_t *exon_phases;
	int8_t *exon_end_phases;
	uint32_t *mature_mirna_offsets;
	uint32_t *mature_mirna_starts;
	uint32_t *mature_mirna_ends;
	size_t mature_mirna_count;
	uint32_t *peptide_edit_offsets;
	uint32_t *peptide_edit_positions;
	uint8_t *peptide_edit_alts;
	size_t peptide_edit_count;
	uint8_t *cds_sequence_bytes;
	size_t cds_sequence_length;
	uint8_t *flank_sequence_bytes;
	size_t flank_sequence_length;
	uint16_t *interval_feature_seq_regions;
	uint32_t *interval_feature_starts;
	uint32_t *interval_feature_ends;
	uint8_t *interval_feature_kinds;
	size_t interval_feature_count;
	cgranges_t *interval_index;
	int interval_index_complete;
	cgranges_t *interval_feature_index;
	int interval_feature_index_complete;
	int transcript_coverage_complete;
	int transcript_flanks_complete;
	size_t known_seq_region_capacity;
	size_t transcript_capacity;
	size_t exon_capacity;
	size_t mature_mirna_capacity;
	size_t peptide_edit_capacity;
	size_t cds_sequence_capacity;
	size_t flank_sequence_capacity;
	size_t interval_feature_capacity;
} duckvep_owned_model_t;

/* Lifted-interval execution view of a model with wrapped circular objects
 * (see kernel/src/duckvep_lift.h). `model` is a facade over the lifted linear
 * arrays: the prepared kernel views, kernel, interval indexes and per-object
 * gene ordinals that annotation consumes. Annotation rows carry lifted object
 * indices until duckvep_lift_resolve maps them to the source model. */
typedef struct duckvep_lifted_model {
	duckvep_owned_model_t model;
	duckvep_lift_t *lift;
} duckvep_lifted_model_t;

/* The model whose views annotation consumes: the lifted view when present. */
static inline const duckvep_owned_model_t *
duckvep_model_active(const duckvep_owned_model_t *model)
{
	return model->lifted != NULL ? &model->lifted->model : model;
}

/* Lift parameters of a region; returns 0 when the region runs unlifted. */
int duckvep_model_region_lift(const duckvep_owned_model_t *model,
	uint16_t seq_region, uint32_t *length, uint32_t *base,
	uint32_t *virtual_length);

typedef struct duckvep_workspace_cache {
	duckvep_workspace_t *workspace;
	/* faidx_t carries mutable seek/decompression state and therefore belongs
	 * to one checked-out worker cache, never the shared immutable model. */
	duckvep_reference_reader_t reference;
	struct duckvep_workspace_cache *next;
} duckvep_workspace_cache_t;

typedef struct duckvep_model_entry {
	char *name;
	duckvep_owned_model_t model;
	duckvep_workspace_cache_t *workspaces;
	size_t pins;
	struct duckvep_model_entry *next;
} duckvep_model_entry_t;

typedef struct duckvep_registry {
	pthread_mutex_t mutex;
	pthread_mutex_t query_mutex;
	duckdb_connection query_connection;
	duckvep_model_entry_t *models;
	void *annotation_state_pool;
	void (*annotation_state_pool_destroy)(void *);
	size_t references;
} duckvep_registry_t;

void duckvep_sql_set_error(char *, size_t, const char *);
int duckvep_sql_resize(void **, size_t, size_t);
size_t duckvep_sql_next_capacity(size_t, size_t);
int duckvep_row_is_null(duckdb_vector, idx_t);
char *duckvep_vector_string(duckdb_vector, idx_t);

duckvep_registry_t *duckvep_registry_create(duckdb_database);
void duckvep_registry_retain(duckvep_registry_t *);
void duckvep_registry_release(void *);
int duckvep_registry_query_acquire(duckvep_registry_t *, char *, size_t);
duckvep_model_entry_t *duckvep_registry_pin(duckvep_registry_t *,
	const char *);
void duckvep_registry_unpin(duckvep_registry_t *, duckvep_model_entry_t *);
duckvep_workspace_cache_t *duckvep_registry_workspace_take(
	duckvep_registry_t *, duckvep_model_entry_t *, char *, size_t);
void duckvep_registry_workspace_return(duckvep_registry_t *,
	duckvep_model_entry_t *, duckvep_workspace_cache_t *);
void duckvep_workspace_cache_destroy(duckvep_workspace_cache_t *);
int duckvep_model_reference_identity_matches(
	const duckvep_owned_model_t *);
void duckvep_register_model_functions(duckdb_connection,
	duckvep_registry_t *);
void duckvep_register_haplotypes(duckdb_connection, duckvep_registry_t *);

#endif /* DUCKVEP_MODEL_H */
