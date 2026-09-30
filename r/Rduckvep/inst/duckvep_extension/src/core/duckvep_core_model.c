#include "duckvep_core_model.h"
#include "kernel/src/duckvep_model_internal.h"
#include "kernel/src/duckvep_budget.h"

#include <htslib/faidx.h>

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#if defined(_WIN32)
#include <io.h>
#include <share.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

static char *
duckvep_string_copy(const char *source)
{
	size_t length;
	char *copy;

	if (source == NULL)
		return NULL;
	length = strlen(source);
	copy = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
	if (copy != NULL)
		memcpy(copy, source, length + 1);
	return copy;
}

int
duckvep_col_is_null(const duckvep_col_t *column, size_t row)
{
	return column->validity != NULL &&
	    ((column->validity[row / 64] >> (row % 64)) & UINT64_C(1)) == 0;
}

/* True for a non-NULL, non-empty string without embedded NUL bytes, i.e. one
 * duckvep_col_string can only fail to copy for lack of memory. */
int
duckvep_col_string_wellformed(const duckvep_col_t *column, size_t row)
{
	const duckvep_string_t *strings;
	uint32_t length;

	if (duckvep_col_is_null(column, row))
		return 0;
	strings = column->data;
	length = duckvep_string_t_length(strings[row]);
	return length != 0 &&
	    memchr(duckvep_string_t_data(&strings[row]), '\0', length) == NULL;
}

char *
duckvep_col_string(const duckvep_col_t *column, size_t row)
{
	const duckvep_string_t *strings;
	const char *data;
	uint32_t length;
	char *copy;

	if (duckvep_col_is_null(column, row))
		return NULL;
	strings = column->data;
	length = duckvep_string_t_length(strings[row]);
	data = duckvep_string_t_data(&strings[row]);
	if (memchr(data, '\0', length) != NULL)
		return NULL;
	copy = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, (size_t)length + 1);
	if (copy == NULL)
		return NULL;
	memcpy(copy, data, length);
	copy[length] = '\0';
	return copy;
}



static void
duckvep_reference_descriptor_close(int descriptor)
{
	if (descriptor < 0)
		return;
#if defined(_WIN32)
	(void)_close(descriptor);
#else
	(void)close(descriptor);
#endif
}

void
duckvep_sql_set_error(char *error, size_t error_size, const char *message)
{
	if (error != NULL && error_size != 0)
		(void)snprintf(error, error_size, "%s",
		    message != NULL ? message : "unknown error");
}

/* If a native-budget or lease refusal is pending on this thread, the message
 * becomes an explicit capacity error that names the exhausted owner. */
const char *
duckvep_sql_final_error(char *out, size_t out_size, const char *error,
	const char *fallback)
{
	char capacity[192];
	const char *detail;

	detail = error != NULL && error[0] != '\0' ? error : fallback;
	if (!duckvep_budget_take_failure(capacity, sizeof(capacity)))
		return detail;
	(void)snprintf(out, out_size, "%s (%s)", capacity, detail);
	return out;
}

int
duckvep_sql_resize(void **pointer, size_t width, size_t count)
{
	void *resized;

	if (width != 0 && count > SIZE_MAX / width)
		return 0;
	resized = duckvep_budget_realloc(DUCKVEP_OWNER_MODEL, *pointer, width * count);
	if (resized == NULL && count != 0)
		return 0;
	*pointer = resized;
	return 1;
}

size_t
duckvep_sql_next_capacity(size_t current, size_t needed)
{
	size_t capacity;

	capacity = current != 0 ? current : 64;
	while (capacity < needed) {
		if (capacity > SIZE_MAX / 2)
			return needed;
		capacity *= 2;
	}
	return capacity;
}

static int
duckvep_model_reserve_regions(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->known_seq_region_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(
	    model->known_seq_region_capacity, needed);
	if (!duckvep_sql_resize((void **)&model->known_seq_regions,
	    sizeof(*model->known_seq_regions), capacity))
		return 0;
	if (!duckvep_sql_resize((void **)&model->sequence_lengths,
	    sizeof(*model->sequence_lengths), capacity))
		return 0;
	if (!duckvep_sql_resize((void **)&model->region_circular,
	    sizeof(*model->region_circular), capacity))
		return 0;
	if (!duckvep_sql_resize((void **)&model->sequence_names,
	    sizeof(*model->sequence_names), capacity))
		return 0;
	model->known_seq_region_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_transcripts(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->transcript_capacity)
		return 1;
	capacity = model->transcript_capacity == 0 ? needed :
	    duckvep_sql_next_capacity(model->transcript_capacity, needed);
#define DUCKVEP_RESIZE_TRANSCRIPT(member) \
	if (!duckvep_sql_resize((void **)&model->member, \
	    sizeof(*model->member), capacity)) \
		return 0
	DUCKVEP_RESIZE_TRANSCRIPT(seq_regions);
	DUCKVEP_RESIZE_TRANSCRIPT(transcript_starts);
	DUCKVEP_RESIZE_TRANSCRIPT(transcript_ends);
	DUCKVEP_RESIZE_TRANSCRIPT(strands);
	DUCKVEP_RESIZE_TRANSCRIPT(transcript_flags);
	DUCKVEP_RESIZE_TRANSCRIPT(gene_indices);
	DUCKVEP_RESIZE_TRANSCRIPT(exon_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(exon_counts);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_starts);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_ends);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_sequence_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(cds_sequence_lengths);
	DUCKVEP_RESIZE_TRANSCRIPT(codon_tables);
	DUCKVEP_RESIZE_TRANSCRIPT(pre_cds_sequence_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(pre_cds_sequence_lengths);
	DUCKVEP_RESIZE_TRANSCRIPT(post_cds_sequence_offsets);
	DUCKVEP_RESIZE_TRANSCRIPT(post_cds_sequence_lengths);
#undef DUCKVEP_RESIZE_TRANSCRIPT
	model->transcript_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_exons(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->exon_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->exon_capacity, needed);
#define DUCKVEP_RESIZE_EXON(member) \
	if (!duckvep_sql_resize((void **)&model->member, \
	    sizeof(*model->member), capacity)) \
		return 0
	DUCKVEP_RESIZE_EXON(exon_starts);
	DUCKVEP_RESIZE_EXON(exon_ends);
	DUCKVEP_RESIZE_EXON(exon_cdna_starts);
	DUCKVEP_RESIZE_EXON(exon_cdna_ends);
	DUCKVEP_RESIZE_EXON(exon_phases);
	DUCKVEP_RESIZE_EXON(exon_end_phases);
#undef DUCKVEP_RESIZE_EXON
	model->exon_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_sequence(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->cds_sequence_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->cds_sequence_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->cds_sequence_bytes,
	    sizeof(*model->cds_sequence_bytes), capacity))
		return 0;
	model->cds_sequence_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_mature_mirna(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->mature_mirna_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->mature_mirna_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->mature_mirna_starts,
	    sizeof(*model->mature_mirna_starts), capacity) ||
	    !duckvep_sql_resize((void **)&model->mature_mirna_ends,
	    sizeof(*model->mature_mirna_ends), capacity))
		return 0;
	model->mature_mirna_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_peptide_edits(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->peptide_edit_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->peptide_edit_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->peptide_edit_positions,
	    sizeof(*model->peptide_edit_positions), capacity) ||
	    !duckvep_sql_resize((void **)&model->peptide_edit_alts,
	    sizeof(*model->peptide_edit_alts), capacity))
		return 0;
	model->peptide_edit_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_flanks(duckvep_owned_model_t *model, size_t needed)
{
	size_t capacity;

	if (needed <= model->flank_sequence_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->flank_sequence_capacity,
	    needed);
	if (!duckvep_sql_resize((void **)&model->flank_sequence_bytes,
	    sizeof(*model->flank_sequence_bytes), capacity))
		return 0;
	model->flank_sequence_capacity = capacity;
	return 1;
}

static int
duckvep_model_reserve_interval_features(duckvep_owned_model_t *model,
	size_t needed)
{
	size_t capacity;

	if (needed <= model->interval_feature_capacity)
		return 1;
	capacity = duckvep_sql_next_capacity(model->interval_feature_capacity,
	    needed);
#define DUCKVEP_RESIZE_INTERVAL_FEATURE(member) \
	if (!duckvep_sql_resize((void **)&model->member, \
	    sizeof(*model->member), capacity)) \
		return 0
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_seq_regions);
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_starts);
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_ends);
	DUCKVEP_RESIZE_INTERVAL_FEATURE(interval_feature_kinds);
#undef DUCKVEP_RESIZE_INTERVAL_FEATURE
	model->interval_feature_capacity = capacity;
	return 1;
}

static void duckvep_lifted_destroy(duckvep_lifted_model_t *);

void
duckvep_owned_model_destroy(duckvep_owned_model_t *model)
{
	if (model == NULL)
		return;
	if (model->lifted != NULL) {
		duckvep_lifted_destroy(model->lifted);
		model->lifted = NULL;
	}
	if (model->kernel != NULL)
		duckvep_model_close(model->kernel);
	duckvep_budget_free(model->known_seq_regions);
	duckvep_budget_free(model->sequence_lengths);
	duckvep_budget_free(model->region_circular);
	if (model->sequence_names != NULL) {
		size_t region;

		for (region = 0; region < model->known_seq_region_count; region++)
			duckvep_budget_free(model->sequence_names[region]);
	}
	duckvep_budget_free(model->sequence_names);
	duckvep_budget_free(model->reference_fasta_path);
	duckvep_budget_free(model->reference_fai_path);
	duckvep_budget_free(model->reference_gzi_path);
	duckvep_budget_free(model->reference_fasta_open_path);
	duckvep_budget_free(model->reference_fai_open_path);
	duckvep_budget_free(model->reference_gzi_open_path);
	if (model->reference_descriptors_open) {
		duckvep_reference_descriptor_close(
		    model->reference_fasta_descriptor);
		duckvep_reference_descriptor_close(model->reference_fai_descriptor);
		duckvep_reference_descriptor_close(model->reference_gzi_descriptor);
	}
	duckvep_budget_free(model->seq_regions);
	duckvep_budget_free(model->transcript_starts);
	duckvep_budget_free(model->transcript_ends);
	duckvep_budget_free(model->strands);
	duckvep_budget_free(model->transcript_flags);
	duckvep_budget_free(model->gene_indices);
	duckvep_budget_free(model->exon_offsets);
	duckvep_budget_free(model->exon_counts);
	duckvep_budget_free(model->cds_starts);
	duckvep_budget_free(model->cds_ends);
	duckvep_budget_free(model->cds_sequence_offsets);
	duckvep_budget_free(model->cds_sequence_lengths);
	duckvep_budget_free(model->codon_tables);
	duckvep_budget_free(model->pre_cds_sequence_offsets);
	duckvep_budget_free(model->pre_cds_sequence_lengths);
	duckvep_budget_free(model->post_cds_sequence_offsets);
	duckvep_budget_free(model->post_cds_sequence_lengths);
	duckvep_budget_free(model->exon_starts);
	duckvep_budget_free(model->exon_ends);
	duckvep_budget_free(model->exon_cdna_starts);
	duckvep_budget_free(model->exon_cdna_ends);
	duckvep_budget_free(model->exon_phases);
	duckvep_budget_free(model->exon_end_phases);
	duckvep_budget_free(model->mature_mirna_offsets);
	duckvep_budget_free(model->mature_mirna_starts);
	duckvep_budget_free(model->mature_mirna_ends);
	duckvep_budget_free(model->peptide_edit_offsets);
	duckvep_budget_free(model->peptide_edit_positions);
	duckvep_budget_free(model->peptide_edit_alts);
	duckvep_budget_free(model->cds_sequence_bytes);
	duckvep_budget_free(model->flank_sequence_bytes);
	duckvep_budget_free(model->interval_feature_seq_regions);
	duckvep_budget_free(model->interval_feature_starts);
	duckvep_budget_free(model->interval_feature_ends);
	duckvep_budget_free(model->interval_feature_kinds);
	if (model->interval_index != NULL) {
		/* cgranges 0.1.1 does not release its interval array. */
		duckvep_budget_free(model->interval_index->r);
		model->interval_index->r = NULL;
		cr_destroy(model->interval_index);
	}
	if (model->interval_feature_index != NULL) {
		/* cgranges 0.1.1 does not release its interval array. */
		duckvep_budget_free(model->interval_feature_index->r);
		model->interval_feature_index->r = NULL;
		cr_destroy(model->interval_feature_index);
	}
	memset(model, 0, sizeof(*model));
}

static void
duckvep_owned_model_publish(duckvep_owned_model_t *model)
{
	model->transcripts.chrom_id = model->seq_regions;
	model->transcripts.start1 = model->transcript_starts;
	model->transcripts.end1 = model->transcript_ends;
	model->transcripts.strand = model->strands;
	model->transcripts.flags = model->transcript_flags;
	model->transcripts.exon_offset = model->exon_offsets;
	model->transcripts.exon_count = model->exon_counts;
	model->transcripts.cds_start1 = model->cds_starts;
	model->transcripts.cds_end1 = model->cds_ends;
	model->transcripts.mature_mirna_offset = model->mature_mirna_offsets;
	model->transcripts.mature_mirna_start1 = model->mature_mirna_starts;
	model->transcripts.mature_mirna_end1 = model->mature_mirna_ends;
	model->transcripts.mature_mirna_count = model->mature_mirna_count;
	model->exons.start1 = model->exon_starts;
	model->exons.end1 = model->exon_ends;
	model->exons.cdna_start1 = model->exon_cdna_starts;
	model->exons.cdna_end1 = model->exon_cdna_ends;
	model->exons.phase = model->exon_phases;
	model->exons.end_phase = model->exon_end_phases;
	model->sequences.cds_bytes = model->cds_sequence_bytes;
	model->sequences.cds_bytes_len = model->cds_sequence_length;
	model->sequences.cds_offset = model->cds_sequence_offsets;
	model->sequences.cds_length = model->cds_sequence_lengths;
	model->sequences.codon_table = model->codon_tables;
	model->sequences.transcript_count = model->transcripts.transcript_count;
	model->sequences.peptide_edit_offset = model->peptide_edit_offsets;
	model->sequences.peptide_edit_position1 = model->peptide_edit_positions;
	model->sequences.peptide_edit_alt = model->peptide_edit_alts;
	model->sequences.peptide_edit_count = model->peptide_edit_count;
	model->sequences.flank_bytes = model->flank_sequence_bytes;
	model->sequences.flank_bytes_len = model->flank_sequence_length;
	model->sequences.pre_cds_offset = model->pre_cds_sequence_offsets;
	model->sequences.pre_cds_length = model->pre_cds_sequence_lengths;
	model->sequences.post_cds_offset = model->post_cds_sequence_offsets;
	model->sequences.post_cds_length = model->post_cds_sequence_lengths;
	model->sequences.flanks_complete =
	    (uint8_t)(model->transcript_flanks_complete != 0);
	model->interval_features.chrom_id =
	    model->interval_feature_seq_regions;
	model->interval_features.start1 = model->interval_feature_starts;
	model->interval_features.end1 = model->interval_feature_ends;
	model->interval_features.kind = model->interval_feature_kinds;
	model->interval_features.feature_count = model->interval_feature_count;
}

static int
duckvep_owned_interval_index(const uint16_t *seq_regions, const uint32_t *starts,
	const uint32_t *ends, size_t count, const char *object_name,
	cgranges_t **result, int *complete, char *error, size_t error_size)
{
	size_t index;
	char message[128], region_name[16];

	*complete = 0;
	if (count > (size_t)INT32_MAX)
		return 1;
	for (index = 0; index < count; index++) {
		if (starts[index] > (uint32_t)INT32_MAX ||
		    ends[index] > (uint32_t)INT32_MAX)
			return 1;
	}
	*result = cr_init();
	if (*result == NULL) {
		(void)snprintf(message, sizeof(message),
		    "out of memory building %s interval index", object_name);
		duckvep_sql_set_error(error, error_size, message);
		return 0;
	}
	/* cr_add registers only occupied regions. Empty pre-registered contigs have
	 * no interval offset and cannot be indexed by cgranges. The model's region
	 * relation, not this lookup accelerator, owns known-region identity. */
	for (index = 0; index < count; index++) {
		(void)snprintf(region_name, sizeof(region_name), "%u",
		    seq_regions[index]);
		if (cr_add(*result, region_name,
		    (int32_t)(starts[index] - 1), (int32_t)ends[index],
		    (int32_t)index) == NULL) {
			(void)snprintf(message, sizeof(message),
			    "could not build %s interval index", object_name);
			duckvep_sql_set_error(error, error_size, message);
			return 0;
		}
	}
	cr_index(*result);
	*complete = 1;
	return 1;
}

static int
duckvep_owned_model_index(duckvep_owned_model_t *model, char *error,
	size_t error_size)
{
	return duckvep_owned_interval_index(model->seq_regions,
	    model->transcript_starts, model->transcript_ends,
	    model->transcripts.transcript_count, "transcript",
	    &model->interval_index, &model->interval_index_complete,
	    error, error_size) &&
	    duckvep_owned_interval_index(
	    model->interval_feature_seq_regions, model->interval_feature_starts,
	    model->interval_feature_ends, model->interval_feature_count,
	    "regulation feature", &model->interval_feature_index,
	    &model->interval_feature_index_complete, error, error_size);
}

static int
duckvep_result_schema(duckvep_source_t *source, const char *const *names,
	const duckvep_ctype_t *types, size_t count, size_t flexible_string_column,
	char *error, size_t error_size)
{
	size_t column;

	if (duckvep_source_columns(source) != (size_t)count) {
		duckvep_sql_set_error(error, error_size,
		    "query returned the wrong number of columns");
		return 0;
	}
	for (column = 0; column < count; column++) {
		const char *name;
		duckvep_ctype_t type;

		name = duckvep_source_column_name(source, (size_t)column);
		type = duckvep_source_column_type(source, (size_t)column);
		if (name == NULL || strcmp(name, names[column]) != 0) {
			(void)snprintf(error, error_size,
			    "query column %zu must be named %s", column + 1,
			    names[column]);
			return 0;
		}
		if (column == flexible_string_column &&
		    (type == DUCKVEP_CT_VARCHAR || type == DUCKVEP_CT_BLOB))
			continue;
		if (type != types[column]) {
			(void)snprintf(error, error_size,
			    "query column %s has the wrong type", names[column]);
			return 0;
		}
	}
	return 1;
}

typedef struct {
	uint16_t region;
	uint32_t length;
	char *name;
	bool circular;
} duckvep_region_row_t;

static int
duckvep_region_row_compare(const void *left, const void *right)
{
	const duckvep_region_row_t *a = left, *b = right;
	return (a->region > b->region) - (a->region < b->region);
}

static int
duckvep_load_regions(duckvep_source_t *source,
	int require_sequence_names, duckvep_owned_model_t *model,
	char *error, size_t error_size)
{
	static const char *const one_name[] = {"seq_region"};
	static const duckvep_ctype_t one_type[] = {DUCKVEP_CT_UINTEGER};
	static const char *const two_names[] = {"seq_region", "sequence_length"};
	static const duckvep_ctype_t two_types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UBIGINT
	};
	static const char *const three_names[] = {
		"seq_region", "sequence_length", "seq_region_name"
	};
	static const char *const four_names[] = {
		"seq_region", "sequence_length", "seq_region_name", "circular"
	};
	static const duckvep_ctype_t four_types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UBIGINT,
		DUCKVEP_CT_VARCHAR, DUCKVEP_CT_BOOLEAN
	};
	static const duckvep_ctype_t three_types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UBIGINT, DUCKVEP_CT_VARCHAR
	};
		duckvep_batch_t *chunk;
	size_t column_count;
	uint8_t *seen;
	int ok;

	seen = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, (size_t)UINT16_MAX + 1u, 1u);
	if (seen == NULL) {
		duckvep_sql_set_error(error, error_size, "out of memory loading sequence regions");
		return 0;
	}
	if (!duckvep_source_open(source, error, error_size)) {
		duckvep_budget_free(seen);
		return 0;
	}
	ok = 0;
	column_count = duckvep_source_columns(source);
	if (column_count != 1 && column_count != 2 && column_count != 3 &&
	    column_count != 4) {
		duckvep_sql_set_error(error, error_size,
		    "seq_region query must return seq_region, optional sequence_length, optional seq_region_name, and optional circular");
		goto done;
	}
	if (model->transcript_coverage_complete && column_count < 2) {
		duckvep_sql_set_error(error, error_size,
		    "complete transcript coverage requires sequence_length for every region");
		goto done;
	}
	if (require_sequence_names && column_count < 3) {
		duckvep_sql_set_error(error, error_size,
		    "reference_fasta requires seq_region, sequence_length, and seq_region_name in the region query");
		goto done;
	}
	if (!duckvep_result_schema(source,
	    column_count == 1 ? one_name :
	    (column_count == 2 ? two_names :
	    (column_count == 3 ? three_names : four_names)),
	    column_count == 1 ? one_type :
	    (column_count == 2 ? two_types :
	    (column_count == 3 ? three_types : four_types)),
	    (size_t)column_count, SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckvep_source_next(source)) != NULL) {
		const duckvep_col_t *region_vector, *length_vector, *name_vector, *circular_vector;
		uint32_t *values;
		uint64_t *lengths;
		size_t row, rows;

		rows = duckvep_batch_rows(chunk);
		region_vector = duckvep_batch_column(chunk, 0);
		length_vector = column_count >= 2
		    ? duckvep_batch_column(chunk, 1) : NULL;
		name_vector = column_count >= 3
		    ? duckvep_batch_column(chunk, 2) : NULL;
		circular_vector = column_count == 4
		    ? duckvep_batch_column(chunk, 3) : NULL;
		values = (uint32_t *)duckvep_col_data(region_vector);
		lengths = length_vector != NULL
		    ? (uint64_t *)duckvep_col_data(length_vector) : NULL;
		for (row = 0; row < rows; row++) {
			size_t index;

			if (duckvep_col_is_null(region_vector, row) ||
			    (length_vector != NULL &&
			    duckvep_col_is_null(length_vector, row)) ||
			    (name_vector != NULL &&
			    duckvep_col_is_null(name_vector, row)) ||
		    (circular_vector != NULL &&
		    duckvep_col_is_null(circular_vector, row))) {
				duckvep_sql_set_error(error, error_size,
				    "seq_region query contains NULL");
				duckvep_batch_release(&chunk);
				goto done;
			}
			index = model->known_seq_region_count;
			if (values[row] > UINT16_MAX) {
				duckvep_sql_set_error(error, error_size,
				    "seq_region exceeds the compact uint16 model limit");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (lengths != NULL &&
			    (lengths[row] == 0 || lengths[row] > UINT32_MAX)) {
				duckvep_sql_set_error(error, error_size,
				    "sequence_length must fit a positive uint32 coordinate");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (seen[values[row]]) {
				duckvep_sql_set_error(error, error_size,
				    "seq_region query contains a duplicate region");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (!duckvep_model_reserve_regions(model, index + 1)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading sequence regions");
				duckvep_batch_release(&chunk);
				goto done;
			}
			model->sequence_names[index] = NULL;
			if (name_vector != NULL) {
				model->sequence_names[index] =
				    duckvep_col_string(name_vector, row);
				if (model->sequence_names[index] == NULL &&
				    duckvep_col_string_wellformed(name_vector, row)) {
					duckvep_sql_set_error(error, error_size,
					    "out of memory copying a sequence region name");
					duckvep_batch_release(&chunk);
					goto done;
				}
				if (model->sequence_names[index] == NULL ||
				    model->sequence_names[index][0] == '\0') {
					duckvep_budget_free(model->sequence_names[index]);
					model->sequence_names[index] = NULL;
					duckvep_sql_set_error(error, error_size,
					    "seq_region_name must be a non-empty string without embedded NUL bytes");
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			model->known_seq_regions[index] = (uint16_t)values[row];
			model->sequence_lengths[index] = lengths != NULL
			    ? (uint32_t)lengths[row] : 0;
			model->region_circular[index] = circular_vector != NULL
			    ? ((bool *)duckvep_col_data(circular_vector))[row] : 0;
			seen[values[row]] = 1u;
			model->known_seq_region_count++;
		}
		duckvep_batch_release(&chunk);
	}
	if (model->known_seq_region_count == 0) {
		duckvep_sql_set_error(error, error_size,
		    "seq_region query returned no rows");
		goto done;
	}
	{
		duckvep_region_row_t *sorted;
		size_t i, n = model->known_seq_region_count;

		sorted = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, n * sizeof(*sorted));
		if (sorted == NULL) {
			duckvep_sql_set_error(error, error_size,
			    "out of memory sorting sequence regions");
			goto done;
		}
		for (i = 0; i < n; i++) {
			sorted[i].region = model->known_seq_regions[i];
			sorted[i].length = model->sequence_lengths[i];
			sorted[i].name = model->sequence_names[i];
			sorted[i].circular = model->region_circular[i];
		}
		qsort(sorted, n, sizeof(*sorted), duckvep_region_row_compare);
		for (i = 0; i < n; i++) {
			model->known_seq_regions[i] = sorted[i].region;
			model->sequence_lengths[i] = sorted[i].length;
			model->sequence_names[i] = sorted[i].name;
			model->region_circular[i] = sorted[i].circular;
		}
		duckvep_budget_free(sorted);
	}
	ok = 1;
done:
	duckvep_budget_free(seen);
	duckvep_source_close(source);
	return ok;
}

static int
duckvep_sequence_name_compare(const void *left, const void *right)
{
	const char *const *a, *const *b;

	a = left;
	b = right;
	return strcmp(*a, *b);
}

static char *
duckvep_reference_path_with_suffix(const char *path, const char *suffix)
{
	size_t path_length, suffix_length;
	char *result;

	if (path == NULL || suffix == NULL)
		return NULL;
	path_length = strlen(path);
	suffix_length = strlen(suffix);
	if (path_length > SIZE_MAX - suffix_length - 1u)
		return NULL;
	result = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, path_length + suffix_length + 1u);
	if (result == NULL)
		return NULL;
	memcpy(result, path, path_length);
	memcpy(result + path_length, suffix, suffix_length + 1u);
	return result;
}

static int
duckvep_reference_file_identity_read_descriptor(int descriptor,
	duckvep_reference_file_identity_t *identity)
{
#if defined(_WIN32)
	struct _stat64 status;
#else
	struct stat status;
#endif

	if (identity == NULL)
		return 0;
	memset(identity, 0, sizeof(*identity));
	if (descriptor < 0)
		return 1;
	errno = 0;
#if defined(_WIN32)
	if (_fstat64(descriptor, &status) != 0)
#else
	if (fstat(descriptor, &status) != 0)
#endif
		return 0;
	if (status.st_size < 0)
		return 0;
	identity->device = (uint64_t)status.st_dev;
	identity->inode = (uint64_t)status.st_ino;
	identity->size = (uint64_t)status.st_size;
	identity->mtime_seconds = (int64_t)status.st_mtime;
#if defined(__APPLE__)
	identity->mtime_nanoseconds =
	    (uint32_t)status.st_mtimespec.tv_nsec;
#elif !defined(_WIN32)
	identity->mtime_nanoseconds = (uint32_t)status.st_mtim.tv_nsec;
#endif
	identity->present = 1;
	return 1;
}

static int
duckvep_reference_descriptor_open(const char *path, int required,
	int *descriptor_out, duckvep_reference_file_identity_t *identity)
{
	int descriptor;

	if (path == NULL || descriptor_out == NULL || identity == NULL)
		return 0;
	*descriptor_out = -1;
	memset(identity, 0, sizeof(*identity));
	errno = 0;
#if defined(_WIN32)
	errno_t open_error;

	descriptor = -1;
	open_error = _sopen_s(&descriptor, path, _O_RDONLY | _O_BINARY,
	    _SH_DENYWR, 0);
	if (open_error != 0) {
		if (!required && open_error == ENOENT)
			return 1;
		return 0;
	}
#else
	{
		int flags;

		flags = O_RDONLY;
#if defined(O_CLOEXEC)
		flags |= O_CLOEXEC;
#endif
		descriptor = open(path, flags);
	}
	if (descriptor < 0) {
		if (!required && errno == ENOENT)
			return 1;
		return 0;
	}
#endif
	if (!duckvep_reference_file_identity_read_descriptor(descriptor,
	    identity)) {
		duckvep_reference_descriptor_close(descriptor);
		return 0;
	}
	*descriptor_out = descriptor;
	return 1;
}

static char *
duckvep_reference_descriptor_path(int descriptor, const char *source_path)
{
#if defined(_WIN32)
	HANDLE handle;
	DWORD needed, written;
	char *path;

	(void)source_path;
	if (descriptor < 0)
		return NULL;
	handle = (HANDLE)_get_osfhandle(descriptor);
	if (handle == INVALID_HANDLE_VALUE)
		return NULL;
	needed = GetFinalPathNameByHandleA(handle, NULL, 0,
	    FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);
	/* `needed + 1` is passed back to a DWORD-sized Win32 API. */
	if (needed == 0 || needed == (DWORD)-1)
		return NULL;
	path = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, (size_t)needed + 1u);
	if (path == NULL)
		return NULL;
	written = GetFinalPathNameByHandleA(handle, path, needed + 1u,
	    FILE_NAME_NORMALIZED | VOLUME_NAME_DOS);
	if (written == 0 || written > needed) {
		duckvep_budget_free(path);
		return NULL;
	}
	return path;
#else
	if (descriptor < 0)
		return NULL;
#if defined(__linux__)
	char buffer[64];
	int written;

	(void)source_path;
	written = snprintf(buffer, sizeof(buffer), "/proc/self/fd/%d",
	    descriptor);
	if (written < 0 || (size_t)written >= sizeof(buffer))
		return NULL;
	return duckvep_string_copy(buffer);
#elif defined(__APPLE__)
	{
		char path[PATH_MAX];

		(void)source_path;
		if (fcntl(descriptor, F_GETPATH, path) != 0)
			return NULL;
		return duckvep_string_copy(path);
	}
#else
	return duckvep_string_copy(source_path);
#endif
#endif
}

static int
duckvep_reference_file_identity_equal(
	const duckvep_reference_file_identity_t *left,
	const duckvep_reference_file_identity_t *right)
{
	if (left == NULL || right == NULL || left->present != right->present)
		return 0;
	if (!left->present)
		return 1;
	/*
	 * POSIX rename can update ctime without changing the open file or its
	 * contents.  Device/inode detect a different object; size/mtime detect
	 * ordinary in-place mutation without rejecting a safely pinned rename.
	 */
	return left->device == right->device &&
	    left->inode == right->inode && left->size == right->size &&
	    left->mtime_seconds == right->mtime_seconds &&
	    left->mtime_nanoseconds == right->mtime_nanoseconds;
}

#if !defined(__linux__) && !defined(_WIN32)
static int
duckvep_reference_open_path_identity_matches(const char *path,
	const duckvep_reference_file_identity_t *expected)
{
	duckvep_reference_file_identity_t observed;
	int descriptor;
	int matches;

	if (expected == NULL)
		return 0;
	if (!expected->present)
		return path == NULL;
	if (path == NULL || !duckvep_reference_descriptor_open(path, 1,
	    &descriptor, &observed))
		return 0;
	matches = duckvep_reference_file_identity_equal(expected, &observed);
	duckvep_reference_descriptor_close(descriptor);
	return matches;
}
#endif

int
duckvep_model_reference_identity_matches(const duckvep_owned_model_t *model)
{
	duckvep_reference_file_identity_t fasta, fai, gzi;

	if (model == NULL)
		return 0;
	if (model->reference_fasta_path == NULL)
		return 1;
	if (!model->reference_descriptors_open ||
	    !duckvep_reference_file_identity_read_descriptor(
	    model->reference_fasta_descriptor, &fasta) ||
	    !duckvep_reference_file_identity_read_descriptor(
	    model->reference_fai_descriptor, &fai) ||
	    !duckvep_reference_file_identity_read_descriptor(
	    model->reference_gzi_descriptor, &gzi))
		return 0;
#if !defined(__linux__) && !defined(_WIN32)
	if (!duckvep_reference_open_path_identity_matches(
	    model->reference_fasta_open_path,
	    &model->reference_fasta_identity) ||
	    !duckvep_reference_open_path_identity_matches(
	    model->reference_fai_open_path,
	    &model->reference_fai_identity) ||
	    !duckvep_reference_open_path_identity_matches(
	    model->reference_gzi_open_path,
	    &model->reference_gzi_identity))
		return 0;
#endif
	return duckvep_reference_file_identity_equal(
	    &model->reference_fasta_identity, &fasta) &&
	    duckvep_reference_file_identity_equal(
	    &model->reference_fai_identity, &fai) &&
	    duckvep_reference_file_identity_equal(
	    &model->reference_gzi_identity, &gzi);
}

static int
duckvep_validate_reference_fasta(const char *reference_fasta,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	faidx_t *fai;
	uint64_t fai_reserved;
	char **sorted_names;
	size_t region;
	int ok;

	if (reference_fasta == NULL)
		return 1;
	model->reference_fasta_path = duckvep_string_copy(reference_fasta);
	model->reference_fai_path = duckvep_reference_path_with_suffix(
	    reference_fasta, ".fai");
	model->reference_gzi_path = duckvep_reference_path_with_suffix(
	    reference_fasta, ".gzi");
	if (model->reference_fasta_path == NULL ||
	    model->reference_fai_path == NULL ||
	    model->reference_gzi_path == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory retaining reference FASTA identity paths");
		return 0;
	}
	model->reference_fasta_descriptor = -1;
	model->reference_fai_descriptor = -1;
	model->reference_gzi_descriptor = -1;
	model->reference_descriptors_open = 1;
	if (!duckvep_reference_descriptor_open(model->reference_fasta_path, 1,
	    &model->reference_fasta_descriptor,
	    &model->reference_fasta_identity)) {
		duckvep_sql_set_error(error, error_size,
		    "could not pin the reference FASTA");
		return 0;
	}
	if (!duckvep_reference_descriptor_open(model->reference_fai_path, 1,
	    &model->reference_fai_descriptor,
	    &model->reference_fai_identity)) {
		duckvep_sql_set_error(error, error_size,
		    "could not pin the reference FASTA .fai index");
		return 0;
	}
	if (!duckvep_reference_descriptor_open(model->reference_gzi_path, 0,
	    &model->reference_gzi_descriptor,
	    &model->reference_gzi_identity)) {
		duckvep_sql_set_error(error, error_size,
		    "could not pin the reference FASTA .gzi index");
		return 0;
	}
	model->reference_fasta_open_path = duckvep_reference_descriptor_path(
	    model->reference_fasta_descriptor, model->reference_fasta_path);
	model->reference_fai_open_path = duckvep_reference_descriptor_path(
	    model->reference_fai_descriptor, model->reference_fai_path);
	if (model->reference_gzi_identity.present)
		model->reference_gzi_open_path =
		    duckvep_reference_descriptor_path(
		    model->reference_gzi_descriptor, model->reference_gzi_path);
	if (model->reference_fasta_open_path == NULL ||
	    model->reference_fai_open_path == NULL ||
	    (model->reference_gzi_identity.present &&
	    model->reference_gzi_open_path == NULL)) {
		duckvep_sql_set_error(error, error_size,
		    "could not retain stable reference descriptor paths");
		return 0;
	}
	if (model->known_seq_region_count == 0 ||
	    model->known_seq_region_count > SIZE_MAX / sizeof(*sorted_names)) {
		duckvep_sql_set_error(error, error_size,
		    "reference FASTA region map exceeds addressable memory");
		return 0;
	}
	sorted_names = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, model->known_seq_region_count *
	    sizeof(*sorted_names));
	if (sorted_names == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory validating reference FASTA region names");
		return 0;
	}
	for (region = 0; region < model->known_seq_region_count; region++) {
		if (model->sequence_names[region] == NULL ||
		    model->sequence_lengths[region] == 0) {
			duckvep_budget_free(sorted_names);
			duckvep_sql_set_error(error, error_size,
			    "reference FASTA requires a name and length for every model region");
			return 0;
		}
		sorted_names[region] = model->sequence_names[region];
	}
	qsort(sorted_names, model->known_seq_region_count,
	    sizeof(*sorted_names), duckvep_sequence_name_compare);
	for (region = 1; region < model->known_seq_region_count; region++) {
		if (strcmp(sorted_names[region - 1], sorted_names[region]) == 0) {
			duckvep_budget_free(sorted_names);
			duckvep_sql_set_error(error, error_size,
			    "seq_region_name values must be unique");
			return 0;
		}
	}
	duckvep_budget_free(sorted_names);
	if (!duckvep_model_reference_identity_matches(model)) {
		duckvep_sql_set_error(error, error_size,
		    "reference FASTA or index changed while the model was loading");
		return 0;
	}

	/* Transient htslib index: reserve its estimated size before opening. */
	fai_reserved = 2u * model->reference_fai_identity.size + 4u * 65536u;
	if (!duckvep_budget_reserve(DUCKVEP_OWNER_REFERENCE, fai_reserved)) {
		duckvep_sql_set_error(error, error_size,
		    "could not reserve the reference index");
		return 0;
	}
	fai = fai_load3_format(model->reference_fasta_open_path,
	    model->reference_fai_open_path, model->reference_gzi_open_path,
	    0, FAI_FASTA);
	if (fai == NULL) {
		duckvep_budget_unreserve(DUCKVEP_OWNER_REFERENCE, fai_reserved);
		duckvep_sql_set_error(error, error_size,
		    "could not open the indexed reference FASTA without creating an index");
		return 0;
	}
	ok = 1;
	for (region = 0; region < model->known_seq_region_count; region++) {
		hts_pos_t length;

		length = faidx_seq_len64(fai, model->sequence_names[region]);
		if (length < 0 || (uint64_t)length !=
		    (uint64_t)model->sequence_lengths[region]) {
			(void)snprintf(error, error_size,
			    "reference FASTA contig %s is absent or has the wrong length",
			    model->sequence_names[region]);
			ok = 0;
			break;
		}
	}
	fai_destroy(fai);
	duckvep_budget_unreserve(DUCKVEP_OWNER_REFERENCE, fai_reserved);
	if (!ok)
		return 0;
	if (!duckvep_model_reference_identity_matches(model)) {
		duckvep_sql_set_error(error, error_size,
		    "reference FASTA or index changed while the model was loading");
		return 0;
	}
	return 1;
}

static int
duckvep_model_region_topology(const duckvep_owned_model_t *model,
	uint32_t seq_region, uint32_t *length, int *circular)
{
	size_t lo = 0, hi = model->known_seq_region_count;
	while (lo < hi) {
		size_t mid = lo + (hi - lo) / 2;
		if (model->known_seq_regions[mid] < seq_region)
			lo = mid + 1;
		else
			hi = mid;
	}
	if (lo == model->known_seq_region_count ||
	    model->known_seq_regions[lo] != seq_region)
		return 0;
	if (length) *length = model->sequence_lengths[lo];
	if (circular) *circular = model->region_circular[lo];
	return 1;
}

static int
duckvep_load_transcripts(duckvep_source_t *source,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "seq_region", "transcript_start",
		"transcript_end", "strand", "gene_index", "transcript_flags",
		"cds_start", "cds_end", "cds_sequence", "codon_table",
		"pre_cds_sequence", "post_cds_sequence"
	};
	static const duckvep_ctype_t types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UINTEGER,
		DUCKVEP_CT_UBIGINT, DUCKVEP_CT_UBIGINT,
		DUCKVEP_CT_TINYINT, DUCKVEP_CT_UINTEGER,
		DUCKVEP_CT_UBIGINT, DUCKVEP_CT_UBIGINT,
		DUCKVEP_CT_UBIGINT, DUCKVEP_CT_BLOB,
		DUCKVEP_CT_UTINYINT, DUCKVEP_CT_BLOB,
		DUCKVEP_CT_BLOB
	};
		duckvep_batch_t *chunk;
	size_t column_count;
	size_t expected, received;
	uint8_t *seen;
	int ok;

	if (!duckvep_source_count(source, &expected, error, error_size))
		return 0;
	if (expected == 0 || expected > UINT32_MAX ||
	    !duckvep_model_reserve_transcripts(model, expected)) {
		duckvep_sql_set_error(error, error_size,
		    "transcript count is empty, exceeds uint32, or cannot be allocated");
		return 0;
	}
	seen = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, expected, 1u);
	if (seen == NULL) {
		duckvep_sql_set_error(error, error_size, "out of memory tracking transcript indexes");
		return 0;
	}
	if (!duckvep_source_open(source, error, error_size)) {
		duckvep_budget_free(seen);
		return 0;
	}
	received = 0;
	ok = 0;
	column_count = duckvep_source_columns(source);
	if (column_count != 11 && column_count != 13) {
		duckvep_sql_set_error(error, error_size,
		    "transcript query must return 11 CDS-only columns or 13 with complete pre_cds_sequence and post_cds_sequence");
		goto done;
	}
	if (!duckvep_result_schema(source, names, types,
	    (size_t)column_count, 9,
	    error, error_size))
		goto done;
	model->transcript_flanks_complete = column_count == 13;
	while ((chunk = duckvep_source_next(source)) != NULL) {
		const duckvep_col_t *vectors[13];
		uint32_t *transcript_indices, *seq_regions, *gene_indices;
		uint64_t *starts, *ends, *flags, *cds_starts, *cds_ends;
		int8_t *strands;
		duckvep_string_t *sequences, *pre_flanks;
		duckvep_string_t *post_flanks;
		uint8_t *tables;
		size_t row, rows;
		size_t column;

		rows = duckvep_batch_rows(chunk);
		for (column = 0; column < (size_t)column_count; column++)
			vectors[column] = duckvep_batch_column(chunk,
			    (size_t)column);
		transcript_indices = duckvep_col_data(vectors[0]);
		seq_regions = duckvep_col_data(vectors[1]);
		starts = duckvep_col_data(vectors[2]);
		ends = duckvep_col_data(vectors[3]);
		strands = duckvep_col_data(vectors[4]);
		gene_indices = duckvep_col_data(vectors[5]);
		flags = duckvep_col_data(vectors[6]);
		cds_starts = duckvep_col_data(vectors[7]);
		cds_ends = duckvep_col_data(vectors[8]);
		sequences = duckvep_col_data(vectors[9]);
		tables = duckvep_col_data(vectors[10]);
		pre_flanks = column_count == 13
		    ? duckvep_col_data(vectors[11]) : NULL;
		post_flanks = column_count == 13
		    ? duckvep_col_data(vectors[12]) : NULL;
		for (row = 0; row < rows; row++) {
			size_t flank_offset, index, post_length, pre_length;
			size_t sequence_length, sequence_offset;
			int cds_nulls, sequence_nulls, circular;
			uint32_t region_length;

			for (column = 0; column < 7; column++) {
				if (duckvep_col_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "transcript query contains NULL in %s",
					    names[column]);
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			index = transcript_indices[row];
			if (index >= expected || seen[index] || received >= expected) {
				duckvep_sql_set_error(error, error_size,
				    "transcript_index is duplicated or outside the dense zero-based range");
				duckvep_batch_release(&chunk);
				goto done;
			}
			circular = 0;
			region_length = 0;
			seen[index] = 1u;
			received++;
			if (seq_regions[row] > UINT16_MAX || starts[row] == 0 ||
			    ends[row] == 0 || starts[row] > UINT32_MAX ||
			    ends[row] > UINT32_MAX ||
			    !duckvep_model_region_topology(model, seq_regions[row],
			    &region_length, &circular) ||
			    (region_length != 0 &&
			    (starts[row] > region_length || ends[row] > region_length)) ||
			    (starts[row] > ends[row] && (!circular || !region_length)) ||
			    (strands[row] != 1 && strands[row] != -1)) {
				duckvep_sql_set_error(error, error_size,
				    "transcript row has an invalid region, span, or strand");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (starts[row] > ends[row])
				model->has_wrapped_coordinates = 1;
			model->seq_regions[index] = (uint16_t)seq_regions[row];
			model->transcript_starts[index] = (uint32_t)starts[row];
			model->transcript_ends[index] = (uint32_t)ends[row];
			model->strands[index] = strands[row];
			model->transcript_flags[index] = flags[row];
			model->gene_indices[index] = gene_indices[row];
			model->exon_offsets[index] = 0;
			model->exon_counts[index] = 0;
			cds_nulls = duckvep_col_is_null(vectors[7], row) +
			    duckvep_col_is_null(vectors[8], row);
			sequence_nulls = duckvep_col_is_null(vectors[9], row) +
			    duckvep_col_is_null(vectors[10], row);
			if (cds_nulls == 1 || sequence_nulls == 1 ||
			    (cds_nulls == 2 && sequence_nulls == 0)) {
				duckvep_sql_set_error(error, error_size,
				    "CDS span and CDS sequence/table must each be both NULL or both present");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (column_count == 13) {
				int pre_null = duckvep_col_is_null(vectors[11], row);
				int post_null = duckvep_col_is_null(vectors[12], row);

				if (pre_null != post_null ||
				    ((sequence_nulls == 0) == (pre_null != 0))) {
					duckvep_sql_set_error(error, error_size,
					    "complete transcript flanks must both be present exactly when CDS sequence is present");
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			sequence_offset = model->cds_sequence_length;
			sequence_length = 0;
			if (cds_nulls == 0) {
				if (cds_starts[row] == 0 || cds_ends[row] == 0 ||
				    cds_starts[row] > UINT32_MAX ||
				    cds_ends[row] > UINT32_MAX ||
				    (region_length != 0 &&
				    (cds_starts[row] > region_length ||
				    cds_ends[row] > region_length)) ||
				    (starts[row] <= ends[row] &&
				    (cds_ends[row] < cds_starts[row] ||
				    cds_starts[row] < starts[row] ||
				    cds_ends[row] > ends[row])) ||
				    (starts[row] > ends[row] &&
				    ((cds_starts[row] < starts[row] &&
				    cds_starts[row] > ends[row]) ||
				    (cds_ends[row] < starts[row] &&
				    cds_ends[row] > ends[row])))) {
					duckvep_sql_set_error(error, error_size,
					    "coding transcript has an invalid CDS span");
					duckvep_batch_release(&chunk);
					goto done;
				}
				model->cds_starts[index] = (uint32_t)cds_starts[row];
				model->cds_ends[index] = (uint32_t)cds_ends[row];
			} else {
				model->cds_starts[index] = 0;
				model->cds_ends[index] = 0;
			}
			model->codon_tables[index] = 1;
			if (sequence_nulls == 0) {
				sequence_length = (size_t)duckvep_string_t_length(
				    sequences[row]);
				if (sequence_length == 0 || sequence_length > UINT32_MAX ||
				    !duckvep_codon_table_supported(
				    (duckvep_codon_table_t)tables[row])) {
					duckvep_sql_set_error(error, error_size,
					    "coding transcript has an invalid CDS sequence or codon table");
					duckvep_batch_release(&chunk);
					goto done;
				}
				if (sequence_length > SIZE_MAX - sequence_offset ||
				    !duckvep_model_reserve_sequence(model,
				    sequence_offset + sequence_length)) {
					duckvep_sql_set_error(error, error_size,
					    "out of memory loading coding sequences");
					duckvep_batch_release(&chunk);
					goto done;
				}
				memcpy(model->cds_sequence_bytes + sequence_offset,
				    duckvep_string_t_data(&sequences[row]),
				    sequence_length);
				model->codon_tables[index] = tables[row];
			}
			flank_offset = model->flank_sequence_length;
			pre_length = 0;
			post_length = 0;
			if (pre_flanks != NULL) {
				if (sequence_nulls == 0) {
					pre_length = (size_t)duckvep_string_t_length(
					    pre_flanks[row]);
					post_length = (size_t)duckvep_string_t_length(
					    post_flanks[row]);
				}
				if (pre_length > UINT32_MAX || post_length > UINT32_MAX) {
					duckvep_sql_set_error(error, error_size,
					    "transcript flank exceeds the uint32 model limit");
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			if (pre_length > SIZE_MAX - flank_offset ||
			    post_length > SIZE_MAX - flank_offset - pre_length ||
			    !duckvep_model_reserve_flanks(model,
			    flank_offset + pre_length + post_length)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading transcript flanks");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (pre_length != 0)
				memcpy(model->flank_sequence_bytes + flank_offset,
				    duckvep_string_t_data(&pre_flanks[row]), pre_length);
			if (post_length != 0)
				memcpy(model->flank_sequence_bytes + flank_offset +
				    pre_length, duckvep_string_t_data(&post_flanks[row]),
				    post_length);
			model->cds_sequence_offsets[index] =
			    (uint64_t)sequence_offset;
			model->cds_sequence_lengths[index] =
			    (uint32_t)sequence_length;
			model->pre_cds_sequence_offsets[index] =
			    (uint64_t)flank_offset;
			model->pre_cds_sequence_lengths[index] =
			    (uint32_t)pre_length;
			model->post_cds_sequence_offsets[index] =
			    (uint64_t)(flank_offset + pre_length);
			model->post_cds_sequence_lengths[index] =
			    (uint32_t)post_length;
			model->cds_sequence_length += sequence_length;
			model->flank_sequence_length += pre_length + post_length;
			model->transcripts.transcript_count++;
		}
		duckvep_batch_release(&chunk);
	}
	if (received != expected) {
		duckvep_sql_set_error(error, error_size,
		    "transcript_index has a gap in the dense zero-based range");
		goto done;
	}
	{
		uint8_t *cds, *flanks;
		size_t i, cds_offset = 0u, flank_offset = 0u;
		int in_order = 1;

		/* Rows normally arrive in transcript_index order, so the pools are already laid out in index order. Then
		 * only a shrink to the exact length is needed; reordering would copy every sequence byte a second time
		 * and hold both copies at once (the load high-water mark). */
		for (i = 0; i < expected && in_order; i++) {
			in_order = model->cds_sequence_offsets[i] == cds_offset &&
			    model->pre_cds_sequence_offsets[i] == flank_offset &&
			    model->post_cds_sequence_offsets[i] == flank_offset +
			    model->pre_cds_sequence_lengths[i];
			cds_offset += model->cds_sequence_lengths[i];
			flank_offset += (size_t)model->pre_cds_sequence_lengths[i] +
			    model->post_cds_sequence_lengths[i];
		}
		if (in_order && cds_offset == model->cds_sequence_length &&
		    flank_offset == model->flank_sequence_length) {
			size_t cds_exact = model->cds_sequence_length == 0 ? 1u : model->cds_sequence_length;
			size_t flank_exact = model->flank_sequence_length == 0 ? 1u : model->flank_sequence_length;

			cds = model->cds_sequence_capacity == cds_exact ? model->cds_sequence_bytes :
			    duckvep_budget_realloc(DUCKVEP_OWNER_MODEL, model->cds_sequence_bytes, cds_exact);
			if (cds != NULL)
				model->cds_sequence_bytes = cds;
			flanks = model->flank_sequence_capacity == flank_exact ? model->flank_sequence_bytes :
			    duckvep_budget_realloc(DUCKVEP_OWNER_MODEL, model->flank_sequence_bytes, flank_exact);
			if (flanks != NULL)
				model->flank_sequence_bytes = flanks;
			if (cds == NULL || flanks == NULL) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory ordering transcript sequences");
				goto done;
			}
			model->cds_sequence_capacity = model->cds_sequence_length;
			model->flank_sequence_capacity = model->flank_sequence_length;
			ok = 1;
			goto done;
		}
		cds_offset = 0u;
		flank_offset = 0u;
		cds = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, model->cds_sequence_length == 0 ? 1u :
		    model->cds_sequence_length);
		flanks = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, model->flank_sequence_length == 0 ? 1u :
		    model->flank_sequence_length);
		if (cds == NULL || flanks == NULL) {
			duckvep_budget_free(cds);
			duckvep_budget_free(flanks);
			duckvep_sql_set_error(error, error_size,
			    "out of memory ordering transcript sequences");
			goto done;
		}
		for (i = 0; i < expected; i++) {
			size_t length = model->cds_sequence_lengths[i];
			size_t pre = model->pre_cds_sequence_lengths[i];
			size_t post = model->post_cds_sequence_lengths[i];

			if (length != 0)
				memcpy(cds + cds_offset, model->cds_sequence_bytes +
				    model->cds_sequence_offsets[i], length);
			if (pre + post != 0)
				memcpy(flanks + flank_offset,
				    model->flank_sequence_bytes +
				    model->pre_cds_sequence_offsets[i], pre + post);
			model->cds_sequence_offsets[i] = cds_offset;
			model->pre_cds_sequence_offsets[i] = flank_offset;
			model->post_cds_sequence_offsets[i] = flank_offset + pre;
			cds_offset += length;
			flank_offset += pre + post;
		}
		duckvep_budget_free(model->cds_sequence_bytes);
		duckvep_budget_free(model->flank_sequence_bytes);
		model->cds_sequence_bytes = cds;
		model->flank_sequence_bytes = flanks;
		model->cds_sequence_capacity = model->cds_sequence_length;
		model->flank_sequence_capacity = model->flank_sequence_length;
	}
	ok = 1;
done:
	duckvep_budget_free(seen);
	duckvep_source_close(source);
	return ok;
}

/* Unordered side-relation rows are streamed into one compact native array,
 * counting-sorted by transcript_index, then sorted per transcript on the
 * secondary key. The scratch is native, bounded by the row count, and freed
 * before the load returns. */
static int
duckvep_rows_append(void **rows, size_t *count, size_t *capacity,
	size_t width, const void *record)
{
	if (*count == *capacity) {
		size_t next = duckvep_sql_next_capacity(*capacity, *count + 1u);

		if (!duckvep_sql_resize(rows, width, next))
			return 0;
		*capacity = next;
	}
	memcpy((char *)*rows + *count * width, record, width);
	(*count)++;
	return 1;
}

/* Counting sort of fixed-width rows by uint32 transcript (first member).
 * offsets has transcript_count + 1 entries holding per-transcript counts in
 * slots 1..n on entry (counted while streaming) and group boundaries on exit. */
static void *
duckvep_rows_group(const void *rows, size_t count, size_t width,
	size_t transcript_count, size_t *offsets)
{
	char *sorted;
	size_t *cursor, i;

	sorted = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, count == 0 ? 1u : count * width);
	cursor = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, (transcript_count + 1u) * sizeof(*cursor));
	if (sorted == NULL || cursor == NULL) {
		duckvep_budget_free(sorted);
		duckvep_budget_free(cursor);
		return NULL;
	}
	for (i = 1u; i <= transcript_count; i++)
		offsets[i] += offsets[i - 1u];
	memcpy(cursor, offsets, (transcript_count + 1u) * sizeof(*cursor));
	for (i = 0; i < count; i++) {
		const char *row = (const char *)rows + i * width;

		memcpy(sorted + cursor[*(const uint32_t *)row]++ * width, row,
		    width);
	}
	duckvep_budget_free(cursor);
	return sorted;
}

typedef struct {
	uint32_t transcript, start, end, cdna_start, cdna_end;
	int8_t phase, end_phase;
} duckvep_exon_row_t;

static int
duckvep_exon_row_compare(const void *left, const void *right)
{
	const duckvep_exon_row_t *a = left, *b = right;

	return (a->cdna_start > b->cdna_start) - (a->cdna_start < b->cdna_start);
}

static int
duckvep_load_exons(duckvep_source_t *source,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "exon_start", "exon_end",
		"exon_cdna_start", "exon_cdna_end", "phase", "end_phase"
	};
	static const duckvep_ctype_t types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UBIGINT,
		DUCKVEP_CT_UBIGINT, DUCKVEP_CT_UBIGINT,
		DUCKVEP_CT_UBIGINT, DUCKVEP_CT_TINYINT,
		DUCKVEP_CT_TINYINT
	};
		duckvep_batch_t *chunk;
	duckvep_exon_row_t *rows, *sorted;
	size_t row_count, row_capacity, *offsets, transcript_count, i;
	int ok;

	transcript_count = model->transcripts.transcript_count;
	rows = sorted = NULL;
	row_count = row_capacity = 0;
	offsets = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, transcript_count + 1u, sizeof(*offsets));
	if (offsets == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading exons");
		return 0;
	}
	if (!duckvep_source_open(source, error, error_size)) {
		duckvep_budget_free(offsets);
		return 0;
	}
	ok = 0;
	if (!duckvep_result_schema(source, names, types, 7, SIZE_MAX,
	    error, error_size))
		goto done;
	while ((chunk = duckvep_source_next(source)) != NULL) {
		const duckvep_col_t *vectors[7];
		uint32_t *transcript_indices;
		uint64_t *starts, *ends, *cdna_starts, *cdna_ends;
		int8_t *phases, *end_phases;
		size_t row, count;
		size_t column;

		count = duckvep_batch_rows(chunk);
		for (column = 0; column < 7; column++)
			vectors[column] = duckvep_batch_column(chunk,
			    (size_t)column);
		transcript_indices = duckvep_col_data(vectors[0]);
		starts = duckvep_col_data(vectors[1]);
		ends = duckvep_col_data(vectors[2]);
		cdna_starts = duckvep_col_data(vectors[3]);
		cdna_ends = duckvep_col_data(vectors[4]);
		phases = duckvep_col_data(vectors[5]);
		end_phases = duckvep_col_data(vectors[6]);
		for (row = 0; row < count; row++) {
			duckvep_exon_row_t record;
			uint32_t transcript_index, region_length;
			int circular;

			for (column = 0; column < 7; column++) {
				if (duckvep_col_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "exon query contains NULL in %s",
					    names[column]);
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			transcript_index = transcript_indices[row];
			if (transcript_index >= transcript_count) {
				duckvep_sql_set_error(error, error_size,
				    "exon transcript_index is outside the loaded transcript range");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (offsets[transcript_index + 1u] == UINT16_MAX) {
				duckvep_sql_set_error(error, error_size,
				    "transcript has too many exons");
				duckvep_batch_release(&chunk);
				goto done;
			}
			region_length = 0;
			circular = 0;
			(void)duckvep_model_region_topology(model,
			    model->seq_regions[transcript_index], &region_length,
			    &circular);
			if (row_count >= UINT32_MAX || starts[row] == 0 ||
			    ends[row] == 0 || starts[row] > UINT32_MAX ||
			    ends[row] > UINT32_MAX ||
			    (region_length != 0 && (starts[row] > region_length ||
			    ends[row] > region_length)) ||
			    (starts[row] > ends[row] && (!circular || !region_length)) ||
			    cdna_starts[row] == 0 ||
			    cdna_starts[row] > UINT32_MAX ||
			    cdna_ends[row] > UINT32_MAX ||
			    cdna_ends[row] < cdna_starts[row] ||
			    (model->transcript_starts[transcript_index] <=
			    model->transcript_ends[transcript_index] &&
			    (starts[row] < model->transcript_starts[transcript_index] ||
			    ends[row] > model->transcript_ends[transcript_index])) ||
			    (model->transcript_starts[transcript_index] >
			    model->transcript_ends[transcript_index] &&
			    ((starts[row] < model->transcript_starts[transcript_index] &&
			    starts[row] > model->transcript_ends[transcript_index]) ||
			    (ends[row] < model->transcript_starts[transcript_index] &&
			    ends[row] > model->transcript_ends[transcript_index]))) ||
			    phases[row] < -1 || phases[row] > 2 ||
			    end_phases[row] < -1 || end_phases[row] > 2) {
				duckvep_sql_set_error(error, error_size,
				    "exon row has invalid coordinates or phase");
				duckvep_batch_release(&chunk);
				goto done;
			}
			record.transcript = transcript_index;
			record.start = (uint32_t)starts[row];
			record.end = (uint32_t)ends[row];
			record.cdna_start = (uint32_t)cdna_starts[row];
			record.cdna_end = (uint32_t)cdna_ends[row];
			record.phase = phases[row];
			record.end_phase = end_phases[row];
			if (!duckvep_rows_append((void **)&rows, &row_count,
			    &row_capacity, sizeof(record), &record)) {
				duckvep_sql_set_error(error, error_size,
				    "out of memory loading exons");
				duckvep_batch_release(&chunk);
				goto done;
			}
			offsets[transcript_index + 1u]++;
		}
		duckvep_batch_release(&chunk);
	}
	sorted = duckvep_rows_group(rows, row_count, sizeof(*rows),
	    transcript_count, offsets);
	duckvep_budget_free(rows);
	rows = NULL;
	if (sorted == NULL || !duckvep_model_reserve_exons(model, row_count)) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory ordering exons");
		goto done;
	}
	for (i = 0; i < transcript_count; i++) {
		size_t first = offsets[i], count = offsets[i + 1u] - offsets[i], k;

		int wrapped_transcript, transcript_circular;
		uint32_t region_length = 0;
		int crossings = 0;

		if (count == 0) {
			duckvep_sql_set_error(error, error_size,
			    "every transcript must have at least one exon");
			goto done;
		}
		qsort(sorted + first, count, sizeof(*sorted),
		    duckvep_exon_row_compare);
		wrapped_transcript = model->transcript_starts[i] >
		    model->transcript_ends[i];
		transcript_circular = 0;
		(void)duckvep_model_region_topology(model, model->seq_regions[i],
		    &region_length, &transcript_circular);
		/* Ordering, overlap and circular traversal checks need the
		 * per-transcript exon rank order, so they run after the sort. */
		for (k = first; k < first + count; k++) {
			int ordered = 1;

			if (k != first) {
				ordered = sorted[k].cdna_start > sorted[k - 1u].cdna_end;
				if (!wrapped_transcript) {
					if (model->strands[i] > 0)
						ordered = ordered && sorted[k].start >
						    sorted[k - 1u].end;
					else
						ordered = ordered && sorted[k].end <
						    sorted[k - 1u].start;
				}
			}
			if (!ordered) {
				duckvep_sql_set_error(error, error_size,
				    "exons must be non-overlapping and ordered on the transcript");
				goto done;
			}
			if (transcript_circular && region_length) {
				uint64_t span = sorted[k].start > sorted[k].end
				    ? (uint64_t)region_length - sorted[k].start + 1u +
				    sorted[k].end
				    : (uint64_t)sorted[k].end - sorted[k].start + 1u;

				if ((uint64_t)sorted[k].cdna_end -
				    sorted[k].cdna_start + 1u != span) {
					duckvep_sql_set_error(error, error_size,
					    "circular exon cDNA length does not match its genomic span");
					goto done;
				}
				crossings += sorted[k].start > sorted[k].end;
				if (k != first)
					crossings += model->strands[i] == 1
					    ? sorted[k - 1u].end > sorted[k].start
					    : sorted[k - 1u].start < sorted[k].end;
			}
			if (sorted[k].start > sorted[k].end)
				model->has_wrapped_coordinates = 1;
			model->exon_starts[k] = sorted[k].start;
			model->exon_ends[k] = sorted[k].end;
			model->exon_cdna_starts[k] = sorted[k].cdna_start;
			model->exon_cdna_ends[k] = sorted[k].cdna_end;
			model->exon_phases[k] = sorted[k].phase;
			model->exon_end_phases[k] = sorted[k].end_phase;
		}
		if (transcript_circular && region_length &&
		    (crossings != (wrapped_transcript ? 1 : 0) ||
		    (wrapped_transcript &&
		    (model->strands[i] == 1
		    ? (sorted[first].start != model->transcript_starts[i] ||
		    sorted[first + count - 1u].end != model->transcript_ends[i])
		    : (sorted[first].end != model->transcript_ends[i] ||
		    sorted[first + count - 1u].start != model->transcript_starts[i]))))) {
			duckvep_sql_set_error(error, error_size,
			    "circular exon ranks must traverse the origin exactly once");
			goto done;
		}
		model->exon_offsets[i] = (uint32_t)first;
		model->exon_counts[i] = (uint16_t)count;
	}
	model->exons.exon_count = row_count;
	ok = 1;
done:
	duckvep_budget_free(rows);
	duckvep_budget_free(sorted);
	duckvep_budget_free(offsets);
	duckvep_source_close(source);
	return ok;
}

static int
duckvep_mature_mirna_segment_is_exonic(
	const duckvep_owned_model_t *model, uint32_t transcript_index,
	uint32_t start1, uint32_t end1)
{
	size_t offset, count, exon;

	offset = model->exon_offsets[transcript_index];
	count = model->exon_counts[transcript_index];
	for (exon = offset; exon < offset + count; exon++) {
		if (model->exon_starts[exon] <= model->exon_ends[exon]) {
			if (start1 <= end1 &&
			    start1 >= model->exon_starts[exon] &&
			    end1 <= model->exon_ends[exon])
				return 1;
		} else if ((start1 > end1 &&
		    start1 >= model->exon_starts[exon] &&
		    end1 <= model->exon_ends[exon]) ||
		    (start1 <= end1 &&
		    (start1 >= model->exon_starts[exon] ||
		    end1 <= model->exon_ends[exon]))) {
			return 1;
		}
	}
	return 0;
}

typedef struct {
	uint32_t transcript, start, end;
} duckvep_mirna_row_t;

static int
duckvep_mirna_row_compare(const void *left, const void *right)
{
	const duckvep_mirna_row_t *a = left, *b = right;

	if (a->start != b->start)
		return a->start < b->start ? -1 : 1;
	return (a->end > b->end) - (a->end < b->end);
}

static int
duckvep_load_mature_mirna(duckvep_source_t *source,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "mature_mirna_start", "mature_mirna_end"
	};
	static const duckvep_ctype_t types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UBIGINT, DUCKVEP_CT_UBIGINT
	};
		duckvep_batch_t *chunk;
	duckvep_mirna_row_t *rows, *sorted;
	size_t row_count, row_capacity, *offsets, transcript_count, i;
	int ok;

	transcript_count = model->transcripts.transcript_count;
	if (transcript_count >
	    SIZE_MAX / sizeof(*model->mature_mirna_offsets) - 1u) {
		duckvep_sql_set_error(error, error_size,
		    "mature-miRNA row-offset array exceeds addressable memory");
		return 0;
	}
	rows = sorted = NULL;
	row_count = row_capacity = 0;
	model->mature_mirna_offsets = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, transcript_count + 1u,
	    sizeof(*model->mature_mirna_offsets));
	offsets = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, transcript_count + 1u, sizeof(*offsets));
	if (model->mature_mirna_offsets == NULL || offsets == NULL) {
		duckvep_budget_free(offsets);
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading mature-miRNA row offsets");
		return 0;
	}
	if (!duckvep_source_open(source, error, error_size)) {
		duckvep_budget_free(offsets);
		return 0;
	}
	ok = 0;
	if (!duckvep_result_schema(source, names, types, 3,
	    SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckvep_source_next(source)) != NULL) {
		const duckvep_col_t *vectors[3];
		uint32_t *transcript_indices;
		uint64_t *starts, *ends;
		size_t row, count;
		size_t column;

		count = duckvep_batch_rows(chunk);
		for (column = 0; column < 3; column++)
			vectors[column] = duckvep_batch_column(chunk,
			    (size_t)column);
		transcript_indices = duckvep_col_data(vectors[0]);
		starts = duckvep_col_data(vectors[1]);
		ends = duckvep_col_data(vectors[2]);
		for (row = 0; row < count; row++) {
			duckvep_mirna_row_t record;
			uint32_t transcript_index, region_length;
			int circular;

			for (column = 0; column < 3; column++) {
				if (duckvep_col_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "mature-miRNA query contains NULL in %s",
					    names[column]);
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			transcript_index = transcript_indices[row];
			if (transcript_index >= transcript_count) {
				duckvep_sql_set_error(error, error_size,
				    "mature-miRNA transcript_index is outside the loaded transcript range");
				duckvep_batch_release(&chunk);
				goto done;
			}
			region_length = 0;
			circular = 0;
			(void)duckvep_model_region_topology(model,
			    model->seq_regions[transcript_index], &region_length,
			    &circular);
			if (starts[row] == 0 || ends[row] == 0 ||
			    starts[row] > UINT32_MAX || ends[row] > UINT32_MAX ||
			    (region_length != 0 &&
			    (starts[row] > region_length || ends[row] > region_length)) ||
			    (starts[row] > ends[row] && !circular) ||
			    (model->transcript_starts[transcript_index] <=
			    model->transcript_ends[transcript_index] &&
			    (starts[row] < model->transcript_starts[transcript_index] ||
			    ends[row] > model->transcript_ends[transcript_index])) ||
			    (model->transcript_starts[transcript_index] >
			    model->transcript_ends[transcript_index] &&
			    ((starts[row] < model->transcript_starts[transcript_index] &&
			    starts[row] > model->transcript_ends[transcript_index]) ||
			    (ends[row] < model->transcript_starts[transcript_index] &&
			    ends[row] > model->transcript_ends[transcript_index]))) ||
			    (model->transcript_flags[transcript_index] &
			    (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA) == 0u ||
			    !duckvep_mature_mirna_segment_is_exonic(model,
			    transcript_index, (uint32_t)starts[row],
			    (uint32_t)ends[row])) {
				duckvep_sql_set_error(error, error_size,
				    "mature-miRNA row is not an exonic interval of a miRNA transcript");
				duckvep_batch_release(&chunk);
				goto done;
			}
			record.transcript = transcript_index;
			record.start = (uint32_t)starts[row];
			record.end = (uint32_t)ends[row];
			if (row_count == (size_t)UINT32_MAX ||
			    !duckvep_rows_append((void **)&rows, &row_count,
			    &row_capacity, sizeof(record), &record)) {
				duckvep_sql_set_error(error, error_size,
				    "mature-miRNA side relation exceeds the uint32 model limit");
				duckvep_batch_release(&chunk);
				goto done;
			}
			offsets[transcript_index + 1u]++;
		}
		duckvep_batch_release(&chunk);
	}
	sorted = duckvep_rows_group(rows, row_count, sizeof(*rows),
	    transcript_count, offsets);
	duckvep_budget_free(rows);
	rows = NULL;
	if (sorted == NULL ||
	    !duckvep_model_reserve_mature_mirna(model, row_count)) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory ordering mature miRNAs");
		goto done;
	}
	for (i = 0; i < transcript_count; i++) {
		size_t k;

		qsort(sorted + offsets[i], offsets[i + 1u] - offsets[i],
		    sizeof(*sorted), duckvep_mirna_row_compare);
		for (k = offsets[i]; k < offsets[i + 1u]; k++) {
			model->mature_mirna_starts[k] = sorted[k].start;
			model->mature_mirna_ends[k] = sorted[k].end;
		}
	}
	for (i = 0; i <= transcript_count; i++)
		model->mature_mirna_offsets[i] = (uint32_t)offsets[i];
	model->mature_mirna_count = row_count;
	ok = 1;
done:
	duckvep_budget_free(rows);
	duckvep_budget_free(sorted);
	duckvep_budget_free(offsets);
	duckvep_source_close(source);
	return ok;
}

static int
duckvep_peptide_edit_alt_valid(uint8_t amino_acid)
{
	return amino_acid == (uint8_t)'*' ||
	    (amino_acid >= (uint8_t)'A' && amino_acid <= (uint8_t)'Z');
}

typedef struct {
	uint32_t transcript, position;
	uint8_t amino_acid;
} duckvep_peptide_row_t;

static int
duckvep_peptide_row_compare(const void *left, const void *right)
{
	const duckvep_peptide_row_t *a = left, *b = right;

	return (a->position > b->position) - (a->position < b->position);
}

static int
duckvep_load_peptide_edits(duckvep_source_t *source,
	duckvep_owned_model_t *model, char *error, size_t error_size)
{
	static const char *const names[] = {
		"transcript_index", "protein_position", "alternate_amino_acid"
	};
	static const duckvep_ctype_t types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UINTEGER, DUCKVEP_CT_VARCHAR
	};
		duckvep_batch_t *chunk;
	duckvep_peptide_row_t *rows, *sorted;
	size_t row_count, row_capacity, *offsets, transcript_count, i;
	int ok;

	transcript_count = model->transcripts.transcript_count;
	if (transcript_count >
	    SIZE_MAX / sizeof(*model->peptide_edit_offsets) - 1u) {
		duckvep_sql_set_error(error, error_size,
		    "peptide-edit row-offset array exceeds addressable memory");
		return 0;
	}
	rows = sorted = NULL;
	row_count = row_capacity = 0;
	model->peptide_edit_offsets = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, transcript_count + 1u,
	    sizeof(*model->peptide_edit_offsets));
	offsets = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, transcript_count + 1u, sizeof(*offsets));
	if (model->peptide_edit_offsets == NULL || offsets == NULL) {
		duckvep_budget_free(offsets);
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading peptide-edit row offsets");
		return 0;
	}
	if (!duckvep_source_open(source, error, error_size)) {
		duckvep_budget_free(offsets);
		return 0;
	}
	ok = 0;
	if (!duckvep_result_schema(source, names, types, 3,
	    SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckvep_source_next(source)) != NULL) {
		const duckvep_col_t *vectors[3];
		uint32_t *transcript_indices, *positions;
		duckvep_string_t *alternates;
		size_t row, count;
		size_t column;

		count = duckvep_batch_rows(chunk);
		for (column = 0; column < 3; column++)
			vectors[column] = duckvep_batch_column(chunk,
			    (size_t)column);
		transcript_indices = duckvep_col_data(vectors[0]);
		positions = duckvep_col_data(vectors[1]);
		alternates = duckvep_col_data(vectors[2]);
		for (row = 0; row < count; row++) {
			duckvep_peptide_row_t record;
			const char *alternate;
			uint32_t transcript_index, position;
			uint32_t alternate_length;
			uint8_t amino_acid;

			for (column = 0; column < 3; column++) {
				if (duckvep_col_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "peptide-edit query contains NULL in %s",
					    names[column]);
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			transcript_index = transcript_indices[row];
			position = positions[row];
			if (transcript_index >= transcript_count) {
				duckvep_sql_set_error(error, error_size,
				    "peptide-edit transcript_index is outside the loaded transcript range");
				duckvep_batch_release(&chunk);
				goto done;
			}
			alternate_length = duckvep_string_t_length(alternates[row]);
			alternate = duckvep_string_t_data(&alternates[row]);
			amino_acid = alternate_length == 1u
			    ? (uint8_t)alternate[0] : 0u;
			if (position == 0 ||
			    position > model->cds_sequence_lengths[transcript_index] / 3u ||
			    alternate_length != 1u ||
			    !duckvep_peptide_edit_alt_valid(amino_acid)) {
				duckvep_sql_set_error(error, error_size,
				    "peptide edit must replace one in-range protein position with one uppercase amino-acid code");
				duckvep_batch_release(&chunk);
				goto done;
			}
			record.transcript = transcript_index;
			record.position = position;
			record.amino_acid = amino_acid;
			if (row_count == (size_t)UINT32_MAX ||
			    !duckvep_rows_append((void **)&rows, &row_count,
			    &row_capacity, sizeof(record), &record)) {
				duckvep_sql_set_error(error, error_size,
				    "peptide-edit side relation exceeds the uint32 model limit");
				duckvep_batch_release(&chunk);
				goto done;
			}
			offsets[transcript_index + 1u]++;
		}
		duckvep_batch_release(&chunk);
	}
	sorted = duckvep_rows_group(rows, row_count, sizeof(*rows),
	    transcript_count, offsets);
	duckvep_budget_free(rows);
	rows = NULL;
	if (sorted == NULL ||
	    !duckvep_model_reserve_peptide_edits(model, row_count)) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory ordering peptide edits");
		goto done;
	}
	for (i = 0; i < transcript_count; i++) {
		size_t k;

		qsort(sorted + offsets[i], offsets[i + 1u] - offsets[i],
		    sizeof(*sorted), duckvep_peptide_row_compare);
		for (k = offsets[i]; k < offsets[i + 1u]; k++) {
			if (k != offsets[i] &&
			    sorted[k].position == sorted[k - 1u].position) {
				duckvep_sql_set_error(error, error_size,
				    "peptide-edit query must be unique by transcript_index and protein_position");
				goto done;
			}
			model->peptide_edit_positions[k] = sorted[k].position;
			model->peptide_edit_alts[k] = sorted[k].amino_acid;
		}
	}
	for (i = 0; i <= transcript_count; i++)
		model->peptide_edit_offsets[i] = (uint32_t)offsets[i];
	model->peptide_edit_count = row_count;
	ok = 1;
done:
	duckvep_budget_free(rows);
	duckvep_budget_free(sorted);
	duckvep_budget_free(offsets);
	duckvep_source_close(source);
	return ok;
}

static int
duckvep_load_interval_features(duckvep_source_t *source,
	duckvep_owned_model_t *model, char *error,
	size_t error_size)
{
	static const char *const names[] = {
		"regulation_feature_index", "seq_region", "feature_start",
		"feature_end", "feature_kind"
	};
	static const duckvep_ctype_t types[] = {
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UINTEGER,
		DUCKVEP_CT_UINTEGER, DUCKVEP_CT_UINTEGER,
		DUCKVEP_CT_UTINYINT
	};
		duckvep_batch_t *chunk;
	size_t expected, received, i;
	uint8_t *seen;
	int ok;

	if (!duckvep_source_count(source, &expected, error, error_size))
		return 0;
	if (expected > UINT32_MAX ||
	    (expected != 0 && !duckvep_model_reserve_interval_features(model,
	    expected))) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading interval features");
		return 0;
	}
	seen = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, expected == 0 ? 1u : expected, 1u);
	if (seen == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory loading interval features");
		return 0;
	}
	if (!duckvep_source_open(source, error, error_size)) {
		duckvep_budget_free(seen);
		return 0;
	}
	ok = 0;
	received = 0;
	if (!duckvep_result_schema(source, names, types, 5,
	    SIZE_MAX, error, error_size))
		goto done;
	while ((chunk = duckvep_source_next(source)) != NULL) {
		const duckvep_col_t *vectors[5];
		uint32_t *indices, *seq_regions, *starts, *ends;
		uint8_t *kinds;
		size_t row, rows;
		size_t column;

		rows = duckvep_batch_rows(chunk);
		for (column = 0u; column < 5u; column++)
			vectors[column] = duckvep_batch_column(chunk,
			    (size_t)column);
		indices = duckvep_col_data(vectors[0]);
		seq_regions = duckvep_col_data(vectors[1]);
		starts = duckvep_col_data(vectors[2]);
		ends = duckvep_col_data(vectors[3]);
		kinds = duckvep_col_data(vectors[4]);
		for (row = 0; row < rows; row++) {
			size_t index;
			uint32_t sequence_length;
			int circular;

			for (column = 0u; column < 5u; column++) {
				if (duckvep_col_is_null(vectors[column], row)) {
					(void)snprintf(error, error_size,
					    "interval-feature query contains NULL in %s",
					    names[column]);
					duckvep_batch_release(&chunk);
					goto done;
				}
			}
			index = indices[row];
			sequence_length = 0u;
			circular = 0;
			if (index >= expected || seen[index] ||
			    seq_regions[row] > UINT16_MAX || starts[row] == 0u ||
			    ends[row] == 0u ||
			    (kinds[row] !=
			    (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION &&
			    kinds[row] !=
			    (uint8_t)DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE) ||
			    !duckvep_model_region_topology(model, seq_regions[row],
			    &sequence_length, &circular) ||
			    (sequence_length != 0u && (starts[row] > sequence_length ||
			    ends[row] > sequence_length)) ||
			    (starts[row] > ends[row] && (!circular || !sequence_length))) {
				duckvep_sql_set_error(error, error_size,
				    "interval-feature row has an invalid dense index, region, coordinate, or kind");
				duckvep_batch_release(&chunk);
				goto done;
			}
			if (starts[row] > ends[row])
				model->has_wrapped_coordinates = 1;
			model->interval_feature_seq_regions[index] =
			    (uint16_t)seq_regions[row];
			model->interval_feature_starts[index] = starts[row];
			model->interval_feature_ends[index] = ends[row];
			model->interval_feature_kinds[index] = kinds[row];
			seen[index] = 1u;
			received++;
		}
		duckvep_batch_release(&chunk);
	}
	if (received != expected) {
		duckvep_sql_set_error(error, error_size,
		    "interval-feature row has an invalid dense index, region, coordinate, or kind");
		goto done;
	}
	for (i = 1u; i < expected; i++) {
		if (model->interval_feature_seq_regions[i] <
		    model->interval_feature_seq_regions[i - 1u] ||
		    (model->interval_feature_seq_regions[i] ==
		    model->interval_feature_seq_regions[i - 1u] &&
		    model->interval_feature_starts[i] <
		    model->interval_feature_starts[i - 1u])) {
			duckvep_sql_set_error(error, error_size,
			    "interval-feature query must be ordered by seq_region, feature_start, and regulation_feature_index");
			goto done;
		}
	}
	model->interval_feature_count = expected;
	ok = 1;
done:
	duckvep_budget_free(seen);
	duckvep_source_close(source);
	return ok;
}

/* ---------------------------------------------------------------- lifting --
 * Circular regions carrying wrapped objects execute on a lifted linear copy of
 * the model (kernel/src/duckvep_lift.c). The source arrays stay the contract
 * for provenance, metadata ordinals and output; the facade below is the one
 * execution authority for annotation. */

int
duckvep_model_region_lift(const duckvep_owned_model_t *model,
	uint16_t seq_region, uint32_t *length, uint32_t *base,
	uint32_t *virtual_length)
{
	if (model == NULL || model->lifted == NULL)
		return 0;
	return duckvep_lift_region(model->lifted->lift, seq_region, length, base,
	    virtual_length);
}

static void
duckvep_lifted_destroy(duckvep_lifted_model_t *lifted)
{
	duckvep_owned_model_t *lm;

	if (lifted == NULL)
		return;
	lm = &lifted->model;
	if (lm->kernel != NULL)
		duckvep_model_close(lm->kernel);
	duckvep_budget_free(lm->gene_indices);
	if (lm->interval_index != NULL) {
		duckvep_budget_free(lm->interval_index->r);
		lm->interval_index->r = NULL;
		cr_destroy(lm->interval_index);
	}
	if (lm->interval_feature_index != NULL) {
		duckvep_budget_free(lm->interval_feature_index->r);
		lm->interval_feature_index->r = NULL;
		cr_destroy(lm->interval_feature_index);
	}
	duckvep_lift_close(lifted->lift);
	duckvep_budget_free(lifted);
}

static int
duckvep_lifted_build(duckvep_owned_model_t *src, char *error,
	size_t error_size)
{
	duckvep_lifted_model_t *lifted;
	duckvep_owned_model_t *lm;
	duckvep_lift_regions_t regions;
	duckvep_error_t kernel_error;
	duckvep_status_t status;
	size_t index;

	regions.chrom_id = src->known_seq_regions;
	regions.length = src->sequence_lengths;
	regions.circular = src->region_circular;
	regions.count = src->known_seq_region_count;
	lifted = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, 1u, sizeof(*lifted));
	if (lifted == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory building the lifted circular model");
		return 0;
	}
	lm = &lifted->model;
	memset(&kernel_error, 0, sizeof(kernel_error));
	status = duckvep_lift_open(&regions, &src->transcripts, &src->exons,
	    &src->sequences, &src->interval_features, &lifted->lift,
	    &kernel_error);
	if (status != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "circular model: %s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "lifting failed");
		duckvep_budget_free(lifted);
		return 0;
	}
	lm->transcripts = lifted->lift->transcripts;
	lm->exons = lifted->lift->exons;
	lm->sequences = lifted->lift->sequences;
	lm->interval_features = lifted->lift->interval_features;
	lm->seq_regions = (uint16_t *)lifted->lift->transcripts.chrom_id;
	lm->transcript_starts = (uint32_t *)lifted->lift->transcripts.start1;
	lm->transcript_ends = (uint32_t *)lifted->lift->transcripts.end1;
	lm->interval_feature_seq_regions =
	    (uint16_t *)lifted->lift->interval_features.chrom_id;
	lm->interval_feature_starts =
	    (uint32_t *)lifted->lift->interval_features.start1;
	lm->interval_feature_ends =
	    (uint32_t *)lifted->lift->interval_features.end1;
	lm->interval_feature_count =
	    lifted->lift->interval_features.feature_count;
	lm->transcript_coverage_complete = src->transcript_coverage_complete;
	lm->transcript_flanks_complete = src->transcript_flanks_complete;
	lm->gene_indices = duckvep_budget_calloc(DUCKVEP_OWNER_MODEL, lm->transcripts.transcript_count + 1u,
	    sizeof(*lm->gene_indices));
	if (lm->gene_indices == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory building the lifted circular model");
		duckvep_lifted_destroy(lifted);
		return 0;
	}
	for (index = 0; index < lm->transcripts.transcript_count; index++)
		lm->gene_indices[index] =
		    src->gene_indices[lifted->lift->transcript_source[index]];
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (duckvep_model_open(&lm->transcripts, &lm->exons, &lm->sequences,
	    &lm->interval_features, &lm->kernel, &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size,
		    "invalid lifted circular transcript model: %s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "kernel validation failed");
		duckvep_lifted_destroy(lifted);
		return 0;
	}
	lm->transcripts = *duckvep_model_prepared_transcripts(lm->kernel);
	lm->exons = *duckvep_model_prepared_exons(lm->kernel);
	lm->sequences = *duckvep_model_prepared_sequences(lm->kernel);
	if (!duckvep_owned_model_index(lm, error, error_size)) {
		duckvep_lifted_destroy(lifted);
		return 0;
	}
	src->lifted = lifted;
	return 1;
}

void
duckvep_model_entry_destroy(duckvep_model_entry_t *entry)
{
	duckvep_workspace_cache_t *workspace, *next;

	if (entry == NULL)
		return;
	for (workspace = entry->workspaces; workspace != NULL;
	    workspace = next) {
		next = workspace->next;
		duckvep_workspace_cache_destroy(workspace);
	}
	duckvep_budget_free(entry->name);
	duckvep_owned_model_destroy(&entry->model);
	duckvep_budget_free(entry);
}

void
duckvep_workspace_cache_destroy(duckvep_workspace_cache_t *cache)
{
	if (cache == NULL)
		return;
	duckvep_reference_reader_close(&cache->reference);
	duckvep_budget_free(cache->reference.bases);
	duckvep_workspace_close(cache->workspace);
	duckvep_budget_free(cache);
}

duckvep_model_entry_t *
duckvep_registry_find_locked(duckvep_registry_t *registry, const char *name)
{
	duckvep_model_entry_t *entry;

	for (entry = registry->models; entry != NULL; entry = entry->next) {
		if (strcmp(entry->name, name) == 0)
			return entry;
	}
	return NULL;
}

duckvep_model_entry_t *
duckvep_registry_pin(duckvep_registry_t *registry, const char *name)
{
	duckvep_model_entry_t *entry;

	pthread_mutex_lock(&registry->mutex);
	entry = duckvep_registry_find_locked(registry, name);
	if (entry != NULL)
		entry->pins++;
	pthread_mutex_unlock(&registry->mutex);
	return entry;
}

void
duckvep_registry_unpin(duckvep_registry_t *registry,
	duckvep_model_entry_t *entry)
{
	if (registry == NULL || entry == NULL)
		return;
	pthread_mutex_lock(&registry->mutex);
	if (entry->pins != 0)
		entry->pins--;
	pthread_mutex_unlock(&registry->mutex);
}

duckvep_workspace_cache_t *
duckvep_registry_workspace_take(duckvep_registry_t *registry,
	duckvep_model_entry_t *entry, char *error, size_t error_size)
{
	duckvep_workspace_cache_t *cache;
	duckvep_error_t kernel_error;

	if (registry == NULL || entry == NULL)
		return NULL;
	pthread_mutex_lock(&registry->mutex);
	cache = entry->workspaces;
	if (cache != NULL) {
		entry->workspaces = cache->next;
		cache->next = NULL;
	}
	pthread_mutex_unlock(&registry->mutex);
	if (cache != NULL)
		return cache;
	cache = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1, sizeof(*cache));
	if (cache == NULL) {
		duckvep_sql_set_error(error, error_size,
		    "out of memory allocating a DuckVEP workspace");
		return NULL;
	}
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (duckvep_workspace_open(duckvep_model_active(&entry->model)->kernel,
	    &cache->workspace,
	    &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "%s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "could not open a DuckVEP workspace");
		duckvep_budget_free(cache);
		return NULL;
	}
	return cache;
}

void
duckvep_registry_workspace_return(duckvep_registry_t *registry,
	duckvep_model_entry_t *entry, duckvep_workspace_cache_t *cache)
{
	if (registry == NULL || entry == NULL || cache == NULL)
		return;
	pthread_mutex_lock(&registry->mutex);
	cache->next = entry->workspaces;
	entry->workspaces = cache;
	pthread_mutex_unlock(&registry->mutex);
}

void
duckvep_registry_retain(duckvep_registry_t *registry)
{
	pthread_mutex_lock(&registry->mutex);
	registry->references++;
	pthread_mutex_unlock(&registry->mutex);
}

void
duckvep_registry_release(void *pointer)
{
	duckvep_registry_t *registry;
	duckvep_model_entry_t *entry, *next;
	int destroy;

	registry = pointer;
	if (registry == NULL)
		return;
	pthread_mutex_lock(&registry->mutex);
	if (registry->references != 0)
		registry->references--;
	destroy = registry->references == 0;
	pthread_mutex_unlock(&registry->mutex);
	if (!destroy)
		return;
	for (entry = registry->models; entry != NULL; entry = next) {
		next = entry->next;
		duckvep_model_entry_destroy(entry);
	}
	if (registry->annotation_state_pool_destroy != NULL)
		registry->annotation_state_pool_destroy(
		    registry->annotation_state_pool);
	if (registry->host_release != NULL)
		registry->host_release(registry);
	pthread_cond_destroy(&registry->admission);
	pthread_mutex_destroy(&registry->query_mutex);
	pthread_mutex_destroy(&registry->mutex);
	duckvep_budget_free(registry);
}
int
duckvep_core_model_load_relations(duckvep_owned_model_t *model,
	const duckvep_model_sources_t *sources, char *error, size_t error_size)
{
	int ok;

	memset(model, 0, sizeof(*model));
	model->transcript_coverage_complete = sources->transcript_coverage_complete;
	ok = duckvep_load_regions(sources->regions,
	    sources->reference_fasta != NULL, model, error, error_size) &&
	    duckvep_validate_reference_fasta(sources->reference_fasta, model, error,
	    error_size) &&
	    duckvep_load_transcripts(sources->transcripts, model, error,
	    error_size) &&
	    duckvep_load_exons(sources->exons, model, error, error_size) &&
	    (sources->mature_mirna == NULL ||
	    duckvep_load_mature_mirna(sources->mature_mirna, model,
	    error, error_size)) &&
	    (sources->peptide_edits == NULL ||
	    duckvep_load_peptide_edits(sources->peptide_edits, model,
	    error, error_size)) &&
	    (sources->interval_features == NULL ||
	    duckvep_load_interval_features(sources->interval_features,
	    model, error, error_size));
	return ok;
}

int
duckvep_core_model_finish(duckvep_owned_model_t *model, char *error, size_t error_size)
{
	duckvep_error_t kernel_error;
	size_t transcript;

	duckvep_owned_model_publish(model);
	for (transcript = 0;
	    transcript < model->transcripts.transcript_count; transcript++) {
		size_t region;
		int found;

		found = 0;
		for (region = 0; region < model->known_seq_region_count; region++) {
			if (model->known_seq_regions[region] ==
			    model->seq_regions[transcript]) {
				found = 1;
				break;
			}
		}
		if (!found) {
			duckvep_sql_set_error(error, error_size,
			    "transcript seq_region is absent from the region query");
			duckvep_owned_model_destroy(model);
			return 0;
		}
	}
	/* Wrapped coordinates cannot enter the linear kernel. Annotation executes
	 * a lifted linear view of the model; the source arrays stay the contract
	 * for provenance, metadata ordinals and output. */
	if (model->has_wrapped_coordinates) {
		if (!duckvep_lifted_build(model, error, error_size)) {
			duckvep_owned_model_destroy(model);
			return 0;
		}
		return 1;
	}
	memset(&kernel_error, 0, sizeof(kernel_error));
	if (duckvep_model_open(&model->transcripts,
	    &model->exons, &model->sequences, &model->interval_features,
	    &model->kernel, &kernel_error) != DUCKVEP_OK) {
		(void)snprintf(error, error_size, "invalid transcript model: %s",
		    kernel_error.message[0] != '\0' ? kernel_error.message :
		    "kernel validation failed");
		duckvep_owned_model_destroy(model);
		return 0;
	}
	/* The kernel validates the borrowed model once and owns the canonical
	 * immutable projection caches. Publish those prepared views to the adapter
	 * instead of retaining a second transcript authority that lacks the caches. */
	model->transcripts = *duckvep_model_prepared_transcripts(model->kernel);
	model->exons = *duckvep_model_prepared_exons(model->kernel);
	model->sequences = *duckvep_model_prepared_sequences(model->kernel);
	if (!duckvep_owned_model_index(model, error, error_size)) {
		duckvep_owned_model_destroy(model);
		return 0;
	}
	return 1;
}

/* ---- Registry creation and installation (hosts without a query connection) ---- */

duckvep_registry_t *
duckvep_core_registry_create(void)
{
	duckvep_registry_t *registry;

	registry = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*registry));
	if (registry == NULL)
		return NULL;
	(void)pthread_mutex_init(&registry->mutex, NULL);
	(void)pthread_mutex_init(&registry->query_mutex, NULL);
	(void)pthread_cond_init(&registry->admission, NULL);
	registry->references = 1;
	return registry;
}

/* Loads the relations and installs the model under `name`. The error is final
 * (budget refusals named). */
int
duckvep_core_model_install(duckvep_registry_t *registry, const char *name,
	const duckvep_model_sources_t *sources, const char *label,
	char *final_error, size_t final_error_size)
{
	duckvep_model_entry_t *entry;
	char error[DUCKVEP_SQL_ERROR_SIZE];
	char message[DUCKVEP_SQL_ERROR_SIZE + 256];
	char fallback[128];
	int loaded;

	duckvep_budget_clear_failure();
	pthread_mutex_lock(&registry->mutex);
	entry = duckvep_registry_find_locked(registry, name);
	pthread_mutex_unlock(&registry->mutex);
	if (entry != NULL) {
		(void)snprintf(final_error, final_error_size, "%s: model name already exists", label);
		return 0;
	}
	entry = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*entry));
	if (entry == NULL || (entry->name = duckvep_string_copy(name)) == NULL) {
		duckvep_model_entry_destroy(entry);
		(void)snprintf(fallback, sizeof(fallback), "%s: out of memory", label);
		duckvep_sql_set_error(final_error, final_error_size,
		    duckvep_sql_final_error(message, sizeof(message), NULL, fallback));
		return 0;
	}
	memset(error, 0, sizeof(error));
	loaded = duckvep_core_model_load_relations(&entry->model, sources, error, sizeof(error)) &&
	    duckvep_core_model_finish(&entry->model, error, sizeof(error));
	if (!loaded) {
		/* finish() destroys the model itself on failure; a failed load leaves it zeroed-or-partial. */
		duckvep_model_entry_destroy(entry);
		(void)snprintf(fallback, sizeof(fallback), "%s: model load failed", label);
		duckvep_sql_set_error(final_error, final_error_size,
		    duckvep_sql_final_error(message, sizeof(message), error, fallback));
		return 0;
	}
	pthread_mutex_lock(&registry->mutex);
	if (duckvep_registry_find_locked(registry, entry->name) != NULL) {
		pthread_mutex_unlock(&registry->mutex);
		duckvep_model_entry_destroy(entry);
		(void)snprintf(final_error, final_error_size, "%s: model name was created concurrently", label);
		return 0;
	}
	entry->next = registry->models;
	registry->models = entry;
	pthread_mutex_unlock(&registry->mutex);
	return 1;
}

/* ---- Model fingerprint: FNV-1a over the loaded arrays ------------------------- */

static uint64_t
fingerprint_bytes(uint64_t hash, const void *data, size_t length)
{
	const unsigned char *bytes = data;
	size_t i;

	for (i = 0; i < length; i++) {
		hash ^= bytes[i];
		hash *= UINT64_C(1099511628211);
	}
	return hash;
}

#define FP(hash, array, count) \
	do { if ((array) != NULL && (count) != 0) \
		(hash) = fingerprint_bytes((hash), (array), (count) * sizeof(*(array))); } while (0)

uint64_t
duckvep_core_model_fingerprint(const duckvep_owned_model_t *model)
{
	uint64_t hash = UINT64_C(14695981039346656037);
	size_t regions = model->known_seq_region_count;
	size_t transcripts = model->transcripts.transcript_count;
	size_t exons = 0, index;

	for (index = 0; index < transcripts; index++)
		exons += model->exon_counts[index];
	FP(hash, &regions, 1);
	FP(hash, &transcripts, 1);
	FP(hash, &exons, 1);
	FP(hash, model->known_seq_regions, regions);
	FP(hash, model->sequence_lengths, regions);
	FP(hash, model->region_circular, regions);
	for (index = 0; model->sequence_names != NULL && index < regions; index++)
		if (model->sequence_names[index] != NULL)
			hash = fingerprint_bytes(hash, model->sequence_names[index],
			    strlen(model->sequence_names[index]) + 1);
	FP(hash, model->seq_regions, transcripts);
	FP(hash, model->transcript_starts, transcripts);
	FP(hash, model->transcript_ends, transcripts);
	FP(hash, model->strands, transcripts);
	FP(hash, model->transcript_flags, transcripts);
	FP(hash, model->gene_indices, transcripts);
	FP(hash, model->exon_offsets, transcripts);
	FP(hash, model->exon_counts, transcripts);
	FP(hash, model->cds_starts, transcripts);
	FP(hash, model->cds_ends, transcripts);
	FP(hash, model->cds_sequence_offsets, transcripts);
	FP(hash, model->cds_sequence_lengths, transcripts);
	FP(hash, model->codon_tables, transcripts);
	FP(hash, model->pre_cds_sequence_offsets, transcripts);
	FP(hash, model->pre_cds_sequence_lengths, transcripts);
	FP(hash, model->post_cds_sequence_offsets, transcripts);
	FP(hash, model->post_cds_sequence_lengths, transcripts);
	FP(hash, model->exon_starts, exons);
	FP(hash, model->exon_ends, exons);
	FP(hash, model->exon_cdna_starts, exons);
	FP(hash, model->exon_cdna_ends, exons);
	FP(hash, model->exon_phases, exons);
	FP(hash, model->exon_end_phases, exons);
	FP(hash, model->cds_sequence_bytes, model->cds_sequence_length);
	FP(hash, model->flank_sequence_bytes, model->flank_sequence_length);
	FP(hash, &model->mature_mirna_count, 1);
	if (model->mature_mirna_count != 0) {
		FP(hash, model->mature_mirna_offsets, transcripts + 1);
		FP(hash, model->mature_mirna_starts, model->mature_mirna_count);
		FP(hash, model->mature_mirna_ends, model->mature_mirna_count);
	}
	FP(hash, &model->peptide_edit_count, 1);
	if (model->peptide_edit_count != 0) {
		FP(hash, model->peptide_edit_offsets, transcripts + 1);
		FP(hash, model->peptide_edit_positions, model->peptide_edit_count);
		FP(hash, model->peptide_edit_alts, model->peptide_edit_count);
	}
	FP(hash, &model->interval_feature_count, 1);
	FP(hash, model->interval_feature_seq_regions, model->interval_feature_count);
	FP(hash, model->interval_feature_starts, model->interval_feature_count);
	FP(hash, model->interval_feature_ends, model->interval_feature_count);
	FP(hash, model->interval_feature_kinds, model->interval_feature_count);
	FP(hash, &model->transcript_coverage_complete, 1);
	FP(hash, &model->transcript_flanks_complete, 1);
	return hash;
}
