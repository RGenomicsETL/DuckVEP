/* Host-neutral resident-model types, loaders and registry: no DuckDB symbol.
 * A host supplies the six model relations as duckvep_source_t row sources. */
#ifndef DUCKVEP_CORE_MODEL_H
#define DUCKVEP_CORE_MODEL_H

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
	/* A model restored from a snapshot borrows its arrays from this storage (a read-only file
	 * mapping, or one block where the platform has no mapping) instead of owning each one. */
	void *snapshot_base;
	size_t snapshot_bytes;
	int snapshot_mapped;
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
	void *query_connection;
	void (*host_release)(struct duckvep_registry *);
	duckvep_model_entry_t *models;
	void *annotation_state_pool;
	size_t annotation_state_pool_count;
	pthread_cond_t admission;
	size_t admitted;
	void (*annotation_state_pool_destroy)(void *);
	size_t references;
} duckvep_registry_t;

void duckvep_sql_set_error(char *, size_t, const char *);
const char *duckvep_sql_final_error(char *, size_t, const char *, const char *);
int duckvep_sql_resize(void **, size_t, size_t);
size_t duckvep_sql_next_capacity(size_t, size_t);

void duckvep_registry_retain(duckvep_registry_t *);
void duckvep_registry_release(void *);
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

/* ---- Row sources -------------------------------------------------------------
 * The loaders read each model relation through a source: a typed, columnar
 * batch stream. v1 backs it with a query on a private connection; v2 with the
 * rows its COPY sink staged. Batches are flat (no selection vector). */

typedef enum {
	DUCKVEP_CT_INVALID = 0, DUCKVEP_CT_BOOLEAN, DUCKVEP_CT_TINYINT, DUCKVEP_CT_UTINYINT,
	DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UBIGINT, DUCKVEP_CT_VARCHAR, DUCKVEP_CT_BLOB,
	DUCKVEP_CT_OTHER
} duckvep_ctype_t;

/* The DuckDB string layout (the same in the v1 and v2 C APIs). */
typedef struct {
	union {
		struct { uint32_t length; char prefix[4]; char *ptr; } pointer;
		struct { uint32_t length; char inlined[12]; } inlined;
	} value;
} duckvep_string_t;

static inline uint32_t duckvep_string_t_length(duckvep_string_t s) { return s.value.inlined.length; }
static inline const char *duckvep_string_t_data(const duckvep_string_t *s) {
	return s->value.inlined.length <= 12 ? s->value.inlined.inlined : s->value.pointer.ptr;
}

typedef struct duckvep_col {
	void *data;
	const uint64_t *validity; /* NULL: every row valid */
} duckvep_col_t;

typedef struct duckvep_batch {
	size_t rows;
	duckvep_col_t *columns;
	void *context;
	void (*release)(struct duckvep_batch *);
} duckvep_batch_t;

typedef struct duckvep_source {
	void *context;
	int (*count)(void *, size_t *, char *, size_t);
	int (*open)(void *, char *, size_t);
	size_t (*column_count)(void *);
	const char *(*column_name)(void *, size_t);
	duckvep_ctype_t (*column_type)(void *, size_t);
	duckvep_batch_t *(*next)(void *);
	void (*close)(void *);
} duckvep_source_t;

static inline int duckvep_source_count(duckvep_source_t *s, size_t *n, char *e, size_t es) { return s->count(s->context, n, e, es); }
static inline int duckvep_source_open(duckvep_source_t *s, char *e, size_t es) { return s->open(s->context, e, es); }
static inline size_t duckvep_source_columns(duckvep_source_t *s) { return s->column_count(s->context); }
static inline const char *duckvep_source_column_name(duckvep_source_t *s, size_t i) { return s->column_name(s->context, i); }
static inline duckvep_ctype_t duckvep_source_column_type(duckvep_source_t *s, size_t i) { return s->column_type(s->context, i); }
static inline duckvep_batch_t *duckvep_source_next(duckvep_source_t *s) { return s->next(s->context); }
static inline void duckvep_source_close(duckvep_source_t *s) { s->close(s->context); }

static inline size_t duckvep_batch_rows(const duckvep_batch_t *b) { return b->rows; }
static inline const duckvep_col_t *duckvep_batch_column(const duckvep_batch_t *b, size_t i) { return &b->columns[i]; }
static inline void duckvep_batch_release(duckvep_batch_t **b) { if (*b) { (*b)->release(*b); *b = NULL; } }
static inline void *duckvep_col_data(const duckvep_col_t *c) { return c->data; }
int duckvep_col_is_null(const duckvep_col_t *, size_t);
char *duckvep_col_string(const duckvep_col_t *, size_t);
int duckvep_col_string_wellformed(const duckvep_col_t *, size_t);

typedef struct duckvep_model_sources {
	duckvep_source_t *regions, *transcripts, *exons;
	duckvep_source_t *mature_mirna, *peptide_edits, *interval_features; /* optional */
	const char *reference_fasta;
	int transcript_coverage_complete;
} duckvep_model_sources_t;

/* Reads the relations into *model (zeroed first). On failure the caller
 * destroys the model. The v1 host runs this inside its private transaction. */
int duckvep_core_model_load_relations(duckvep_owned_model_t *model,
	const duckvep_model_sources_t *sources, char *error, size_t error_size);
/* Publishes, checks region references, lifts circular models and opens the
 * kernel. Destroys the model on failure. */
int duckvep_core_model_finish(duckvep_owned_model_t *model, char *error, size_t error_size);
void duckvep_owned_model_destroy(duckvep_owned_model_t *model);
duckvep_registry_t *duckvep_core_registry_create(void);
int duckvep_core_model_install(duckvep_registry_t *registry, const char *name,
	const duckvep_model_sources_t *sources, const char *label,
	char *final_error, size_t final_error_size);
uint64_t duckvep_core_model_fingerprint(const duckvep_owned_model_t *model);
/* Opens and checks the model's reference FASTA against its regions, as a relation load does. */
int duckvep_core_model_validate_reference(const char *reference_fasta, duckvep_owned_model_t *model,
	char *error, size_t error_size);
void duckvep_model_entry_destroy(duckvep_model_entry_t *entry);
duckvep_model_entry_t *duckvep_registry_find_locked(duckvep_registry_t *, const char *);

#endif /* DUCKVEP_CORE_MODEL_H */
