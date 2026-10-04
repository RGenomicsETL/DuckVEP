/* Model snapshots: the validated arrays of a loaded model written to one file, and that file mapped
 * read-only as the arrays of a new model.
 *
 * Loading a model from relations copies every byte several times (table scan, sort, materialized
 * result, native arrays). A snapshot is the native arrays themselves: restoring maps the file, so the
 * operating system shares one copy between every process that restores it, and the load costs the
 * validation passes only.
 *
 * A snapshot is not trusted. The header, section bounds and a checksum of every byte are verified, the
 * region and coordinate rules of the relation loaders that memory safety rests on are checked again,
 * and duckvep_core_model_finish runs the kernel's full model validation, as it does for a relation load. */
#include "duckvep_core_snapshot.h"

#include "kernel/src/duckvep_budget.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#define DUCKVEP_SNAPSHOT_MMAP 1
#else
#define DUCKVEP_SNAPSHOT_MMAP 0
#endif

#define SNAP_MAGIC "DVEPSNAP"
#define SNAP_VERSION UINT64_C(1)
#define SNAP_ENDIAN UINT64_C(0x0102030405060708)
#define SNAP_ALIGN ((uint64_t)64)
#define SNAP_ARRAYS 38u
#define SNAP_NAMES SNAP_ARRAYS
#define SNAP_REFERENCE (SNAP_ARRAYS + 1u)
#define SNAP_SECTIONS (SNAP_ARRAYS + 2u)

enum { SNAP_REGIONS, SNAP_TRANSCRIPTS, SNAP_EXONS, SNAP_MIRNA, SNAP_PEPTIDE, SNAP_CDS, SNAP_FLANK,
       SNAP_FEATURES, SNAP_COUNTS };
enum { SNAP_FLAG_COVERAGE = 1u, SNAP_FLAG_FLANKS = 2u };

/* Every field is a uint64_t, so the layout has no padding. */
typedef struct {
	char magic[8];
	uint64_t version;
	uint64_t endian;
	uint64_t header_bytes;
	uint64_t file_bytes;
	uint64_t checksum;      /* of the header with this field zero, then of every section in order */
	uint64_t fingerprint;   /* duckvep_core_model_fingerprint of the saved model; informational */
	uint64_t counts[SNAP_COUNTS];
	uint64_t flags;
	uint64_t present;       /* bit i: array i exists (an array may exist with zero elements) */
	uint64_t width[SNAP_SECTIONS];
	uint64_t offset[SNAP_SECTIONS];
	uint64_t bytes[SNAP_SECTIONS];
} snap_header_t;

typedef struct {
	void **slot;
	size_t width;
	uint64_t count;
} snap_array_t;

/* The arrays of a model, in file order, with the element count each one must have. */
static void
snap_arrays(duckvep_owned_model_t *m, const uint64_t counts[SNAP_COUNTS], snap_array_t out[SNAP_ARRAYS])
{
	uint64_t r = counts[SNAP_REGIONS], t = counts[SNAP_TRANSCRIPTS], e = counts[SNAP_EXONS];
	uint64_t mi = counts[SNAP_MIRNA], p = counts[SNAP_PEPTIDE], f = counts[SNAP_FEATURES];
	size_t i = 0;

#define A(field, n) do { out[i].slot = (void **)&m->field; out[i].width = sizeof(*m->field); \
	out[i].count = (n); i++; } while (0)
	A(known_seq_regions, r); A(sequence_lengths, r); A(region_circular, r);
	A(seq_regions, t); A(transcript_starts, t); A(transcript_ends, t); A(strands, t);
	A(transcript_flags, t); A(gene_indices, t); A(exon_offsets, t); A(exon_counts, t);
	A(cds_starts, t); A(cds_ends, t); A(cds_sequence_offsets, t); A(cds_sequence_lengths, t);
	A(codon_tables, t); A(pre_cds_sequence_offsets, t); A(pre_cds_sequence_lengths, t);
	A(post_cds_sequence_offsets, t); A(post_cds_sequence_lengths, t);
	A(exon_starts, e); A(exon_ends, e); A(exon_cdna_starts, e); A(exon_cdna_ends, e);
	A(exon_phases, e); A(exon_end_phases, e);
	A(mature_mirna_offsets, t + 1u); A(mature_mirna_starts, mi); A(mature_mirna_ends, mi);
	A(peptide_edit_offsets, t + 1u); A(peptide_edit_positions, p); A(peptide_edit_alts, p);
	A(cds_sequence_bytes, counts[SNAP_CDS]); A(flank_sequence_bytes, counts[SNAP_FLANK]);
	A(interval_feature_seq_regions, f); A(interval_feature_starts, f); A(interval_feature_ends, f);
	A(interval_feature_kinds, f);
#undef A
}

static void
snap_counts(const duckvep_owned_model_t *m, uint64_t counts[SNAP_COUNTS])
{
	counts[SNAP_REGIONS] = m->known_seq_region_count;
	counts[SNAP_TRANSCRIPTS] = m->transcripts.transcript_count;
	counts[SNAP_EXONS] = m->exons.exon_count;
	counts[SNAP_MIRNA] = m->mature_mirna_count;
	counts[SNAP_PEPTIDE] = m->peptide_edit_count;
	counts[SNAP_CDS] = m->cds_sequence_length;
	counts[SNAP_FLANK] = m->flank_sequence_length;
	counts[SNAP_FEATURES] = m->interval_feature_count;
}

static uint64_t
snap_rotl(uint64_t value, unsigned shift)
{
	return (value << shift) | (value >> (64u - shift));
}

/* A 64-bit checksum, four independent lanes of eight bytes. It detects damage; it is not a signature. */
static uint64_t
snap_hash(uint64_t seed, const void *data, size_t length)
{
	const unsigned char *bytes = data;
	uint64_t lane[4];
	size_t i;

	lane[0] = seed ^ UINT64_C(0x9E3779B97F4A7C15);
	lane[1] = seed ^ UINT64_C(0xC2B2AE3D27D4EB4F);
	lane[2] = seed ^ UINT64_C(0x165667B19E3779F9);
	lane[3] = seed ^ UINT64_C(0x27D4EB2F165667C5);
	while (length >= 32u) {
		uint64_t word[4];

		memcpy(word, bytes, 32u);
		for (i = 0; i < 4u; i++)
			lane[i] = snap_rotl(lane[i] ^ word[i], 29u) * UINT64_C(0x9FB21C651E98DF25);
		bytes += 32u;
		length -= 32u;
	}
	for (i = 0; i < length; i++)
		lane[i & 3u] = snap_rotl(lane[i & 3u] ^ bytes[i], 11u) * UINT64_C(0x9FB21C651E98DF25);
	lane[0] ^= snap_rotl(lane[1], 17u) ^ snap_rotl(lane[2], 31u) ^ snap_rotl(lane[3], 47u) ^ (uint64_t)length;
	lane[0] ^= lane[0] >> 32;
	lane[0] *= UINT64_C(0xD6E8FEB86659FD93);
	lane[0] ^= lane[0] >> 32;
	return lane[0];
}

static void
snap_error(char *error, size_t error_size, const char *label, const char *message, const char *path)
{
	if (path != NULL)
		(void)snprintf(error, error_size, "%s: %s: %s", label, message, path);
	else
		(void)snprintf(error, error_size, "%s: %s", label, message);
}

/* ---- Save -------------------------------------------------------------------------------------- */

static char *
snap_names_blob(const duckvep_owned_model_t *model, size_t *bytes)
{
	size_t region, total = 0, at = 0;
	char *blob;

	*bytes = 0;
	if (model->sequence_names == NULL)
		return NULL;
	for (region = 0; region < model->known_seq_region_count; region++)
		total += 1u + (model->sequence_names[region] != NULL
		    ? strlen(model->sequence_names[region]) + 1u : 0u);
	blob = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, total != 0 ? total : 1u);
	if (blob == NULL)
		return NULL;
	for (region = 0; region < model->known_seq_region_count; region++) {
		const char *name = model->sequence_names[region];

		blob[at++] = name != NULL;
		if (name != NULL) {
			size_t length = strlen(name) + 1u;

			memcpy(blob + at, name, length);
			at += length;
		}
	}
	*bytes = total;
	return blob;
}

int
duckvep_core_model_snapshot_save(const duckvep_owned_model_t *model, const char *path,
	char *error, size_t error_size)
{
	static const char label[] = "duckvep_model_save";
	static const unsigned char padding[SNAP_ALIGN];
	snap_header_t header;
	snap_array_t arrays[SNAP_ARRAYS];
	const void *data[SNAP_SECTIONS];
	duckvep_owned_model_t view = *model; /* snap_arrays takes slots; the copy is only read */
	char *names = NULL, *partial = NULL;
	size_t names_bytes = 0, path_length, section;
	uint64_t at;
	FILE *file = NULL;
	int ok = 0;

	memset(&header, 0, sizeof(header));
	memcpy(header.magic, SNAP_MAGIC, sizeof(header.magic));
	header.version = SNAP_VERSION;
	header.endian = SNAP_ENDIAN;
	header.header_bytes = sizeof(header);
	header.fingerprint = duckvep_core_model_fingerprint(model);
	snap_counts(model, header.counts);
	header.flags = (model->transcript_coverage_complete ? SNAP_FLAG_COVERAGE : 0u) |
	    (model->transcript_flanks_complete ? SNAP_FLAG_FLANKS : 0u);
	snap_arrays(&view, header.counts, arrays);
	for (section = 0; section < SNAP_ARRAYS; section++) {
		data[section] = *arrays[section].slot;
		header.width[section] = arrays[section].width;
		if (data[section] != NULL) {
			header.present |= UINT64_C(1) << section;
			header.bytes[section] = arrays[section].count * arrays[section].width;
		}
	}
	names = snap_names_blob(model, &names_bytes);
	if (model->sequence_names != NULL && names == NULL) {
		snap_error(error, error_size, label, "out of memory", NULL);
		goto done;
	}
	data[SNAP_NAMES] = names;
	header.width[SNAP_NAMES] = 1u;
	header.bytes[SNAP_NAMES] = names_bytes;
	data[SNAP_REFERENCE] = model->reference_fasta_path;
	header.width[SNAP_REFERENCE] = 1u;
	header.bytes[SNAP_REFERENCE] = model->reference_fasta_path != NULL
	    ? strlen(model->reference_fasta_path) + 1u : 0u;
	at = (sizeof(header) + SNAP_ALIGN - 1u) / SNAP_ALIGN * SNAP_ALIGN;
	for (section = 0; section < SNAP_SECTIONS; section++) {
		header.offset[section] = at;
		at += (header.bytes[section] + SNAP_ALIGN - 1u) / SNAP_ALIGN * SNAP_ALIGN;
	}
	header.file_bytes = at;
	header.checksum = snap_hash(SNAP_VERSION, &header, sizeof(header));
	for (section = 0; section < SNAP_SECTIONS; section++)
		header.checksum = snap_hash(header.checksum, data[section], (size_t)header.bytes[section]);

	/* Written beside the target and renamed, so a reader never maps a partial file. */
	path_length = strlen(path);
	partial = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, path_length + sizeof(".partial"));
	if (partial == NULL) {
		snap_error(error, error_size, label, "out of memory", NULL);
		goto done;
	}
	memcpy(partial, path, path_length);
	memcpy(partial + path_length, ".partial", sizeof(".partial"));
	file = fopen(partial, "wb");
	if (file == NULL) {
		snap_error(error, error_size, label, "cannot create the snapshot file", partial);
		goto done;
	}
	ok = fwrite(&header, sizeof(header), 1, file) == 1;
	at = sizeof(header);
	for (section = 0; ok && section < SNAP_SECTIONS; section++) {
		size_t gap = (size_t)(header.offset[section] - at);

		ok = (gap == 0 || fwrite(padding, 1, gap, file) == gap) &&
		    (header.bytes[section] == 0 ||
		    fwrite(data[section], 1, (size_t)header.bytes[section], file) == header.bytes[section]);
		at = header.offset[section] + header.bytes[section];
	}
	if (ok && header.file_bytes > at)
		ok = fwrite(padding, 1, (size_t)(header.file_bytes - at), file) == header.file_bytes - at;
	if (fclose(file) != 0)
		ok = 0;
	file = NULL;
#ifdef _WIN32
	if (ok)
		(void)remove(path); /* rename does not replace an existing file there */
#endif
	if (ok && rename(partial, path) != 0)
		ok = 0;
	if (!ok) {
		(void)remove(partial);
		snap_error(error, error_size, label, "cannot write the snapshot file", path);
	}
done:
	duckvep_budget_free(names);
	duckvep_budget_free(partial);
	return ok;
}

/* ---- Restore ----------------------------------------------------------------------------------- */

void
duckvep_core_snapshot_release(void *base, size_t bytes, int mapped)
{
	if (base == NULL)
		return;
#if DUCKVEP_SNAPSHOT_MMAP
	if (mapped) {
		(void)munmap(base, bytes);
		duckvep_budget_unreserve(DUCKVEP_OWNER_MODEL, bytes);
		return;
	}
#else
	(void)bytes;
	(void)mapped;
#endif
	duckvep_budget_free(base);
}

/* The whole file as read-only memory: a private mapping where the platform has one, else one block. */
static void *
snap_open(const char *path, size_t *bytes, int *mapped, const char *label, char *error, size_t error_size)
{
	void *base = NULL;

	*bytes = 0;
	*mapped = 0;
#if DUCKVEP_SNAPSHOT_MMAP
	{
		struct stat status;
		int descriptor = open(path, O_RDONLY);

		if (descriptor < 0 || fstat(descriptor, &status) != 0 || !S_ISREG(status.st_mode) ||
		    status.st_size < (off_t)sizeof(snap_header_t)) {
			if (descriptor >= 0)
				(void)close(descriptor);
			snap_error(error, error_size, label, "cannot open the snapshot file", path);
			return NULL;
		}
		/* The mapping is charged to the model owner: it is the model's memory, shared or not. */
		if (!duckvep_budget_reserve(DUCKVEP_OWNER_MODEL, (uint64_t)status.st_size)) {
			(void)close(descriptor);
			snap_error(error, error_size, label, "the snapshot does not fit the native budget", path);
			return NULL;
		}
		base = mmap(NULL, (size_t)status.st_size, PROT_READ, MAP_PRIVATE, descriptor, 0);
		(void)close(descriptor);
		if (base == MAP_FAILED) {
			duckvep_budget_unreserve(DUCKVEP_OWNER_MODEL, (uint64_t)status.st_size);
			snap_error(error, error_size, label, "cannot map the snapshot file", path);
			return NULL;
		}
		*bytes = (size_t)status.st_size;
		*mapped = 1;
	}
#else
	{
		FILE *file = fopen(path, "rb");
		long size = -1;

		if (file != NULL && fseek(file, 0, SEEK_END) == 0)
			size = ftell(file);
		if (file == NULL || size < (long)sizeof(snap_header_t) || fseek(file, 0, SEEK_SET) != 0) {
			if (file != NULL)
				(void)fclose(file);
			snap_error(error, error_size, label, "cannot open the snapshot file", path);
			return NULL;
		}
		base = duckvep_budget_malloc(DUCKVEP_OWNER_MODEL, (size_t)size);
		if (base == NULL || fread(base, 1, (size_t)size, file) != (size_t)size) {
			(void)fclose(file);
			duckvep_budget_free(base);
			snap_error(error, error_size, label, "cannot read the snapshot file", path);
			return NULL;
		}
		(void)fclose(file);
		*bytes = (size_t)size;
	}
#endif
	return base;
}

static int
snap_region(const duckvep_owned_model_t *m, uint16_t region, uint32_t *length, int *circular)
{
	size_t lo = 0, hi = m->known_seq_region_count;

	while (lo < hi) {
		size_t mid = lo + (hi - lo) / 2u;

		if (m->known_seq_regions[mid] < region)
			lo = mid + 1u;
		else
			hi = mid;
	}
	if (lo == m->known_seq_region_count || m->known_seq_regions[lo] != region)
		return 0;
	*length = m->sequence_lengths[lo];
	*circular = m->region_circular[lo] != 0;
	return 1;
}

/* A span on its region: one-based, inside a known length, and start > end only on a circular region. */
static int
snap_span(uint32_t start, uint32_t end, uint32_t length, int circular, int *wrapped)
{
	if (start == 0 || end == 0 || (length != 0 && (start > length || end > length)) ||
	    (start > end && (!circular || length == 0)))
		return 0;
	if (start > end)
		*wrapped = 1;
	return 1;
}

/* Row offsets of a per-transcript side relation: absent with no rows, or zero-based, monotone and
 * ending at the row count, with the row columns present. */
static int
snap_offsets(const uint32_t *offsets, size_t transcripts, size_t rows, int columns)
{
	size_t i;

	if (offsets == NULL)
		return rows == 0;
	if (offsets[0] != 0 || offsets[transcripts] != rows || (rows != 0 && !columns))
		return 0;
	for (i = 0; i < transcripts; i++) {
		if (offsets[i] > offsets[i + 1u])
			return 0;
	}
	return 1;
}

/* Per-transcript slices of a byte pool. */
static int
snap_slices(const uint64_t *offsets, const uint32_t *lengths, size_t transcripts, const uint8_t *pool,
	size_t pool_bytes)
{
	size_t i;

	if (offsets == NULL || lengths == NULL)
		return offsets == NULL && lengths == NULL;
	for (i = 0; i < transcripts; i++) {
		if (lengths[i] != 0 && (pool == NULL || offsets[i] > pool_bytes || lengths[i] > pool_bytes - offsets[i]))
			return 0;
	}
	return 1;
}

/* The rules of the relation loaders that the host relies on before the kernel sees the model. */
static const char *
snap_validate(duckvep_owned_model_t *m)
{
	size_t regions = m->known_seq_region_count, transcripts = m->transcripts.transcript_count;
	size_t exons = m->exons.exon_count, i, k;
	uint64_t exon_total = 0;
	int wrapped = 0;

	for (i = 0; i < regions; i++) {
		if (m->region_circular[i] > 1u || (i != 0 && m->known_seq_regions[i] <= m->known_seq_regions[i - 1u]))
			return "sequence regions are not strictly ascending or carry an invalid circular flag";
	}
	for (i = 0; i < transcripts; i++) {
		uint32_t length = 0;
		int circular = 0;

		if (!snap_region(m, m->seq_regions[i], &length, &circular) ||
		    !snap_span(m->transcript_starts[i], m->transcript_ends[i], length, circular, &wrapped) ||
		    (m->strands[i] != 1 && m->strands[i] != -1))
			return "a transcript has an invalid region, span or strand";
		if (m->exon_offsets[i] != exon_total || m->exon_counts[i] > exons - exon_total)
			return "transcript exon ranges do not tile the exon arrays";
		for (k = (size_t)exon_total; k < (size_t)exon_total + m->exon_counts[i]; k++) {
			if (!snap_span(m->exon_starts[k], m->exon_ends[k], length, circular, &wrapped))
				return "an exon has an invalid span";
		}
		exon_total += m->exon_counts[i];
	}
	if (exon_total != exons)
		return "transcript exon ranges do not tile the exon arrays";
	for (i = 0; i < m->interval_feature_count; i++) {
		uint32_t length = 0;
		int circular = 0;

		if (!snap_region(m, m->interval_feature_seq_regions[i], &length, &circular) ||
		    !snap_span(m->interval_feature_starts[i], m->interval_feature_ends[i], length, circular, &wrapped))
			return "an interval feature has an invalid region or span";
	}
	/* Side relations and sequence slices: the lifted build of a wrapped model reads them before the
	 * kernel validates the lifted arrays. */
	if (!snap_offsets(m->mature_mirna_offsets, transcripts, m->mature_mirna_count,
	    m->mature_mirna_starts != NULL && m->mature_mirna_ends != NULL) ||
	    !snap_offsets(m->peptide_edit_offsets, transcripts, m->peptide_edit_count,
	    m->peptide_edit_positions != NULL && m->peptide_edit_alts != NULL))
		return "a side relation's row offsets do not cover it";
	if (!snap_slices(m->cds_sequence_offsets, m->cds_sequence_lengths, transcripts,
	    m->cds_sequence_bytes, m->cds_sequence_length) ||
	    !snap_slices(m->pre_cds_sequence_offsets, m->pre_cds_sequence_lengths, transcripts,
	    m->flank_sequence_bytes, m->flank_sequence_length) ||
	    !snap_slices(m->post_cds_sequence_offsets, m->post_cds_sequence_lengths, transcripts,
	    m->flank_sequence_bytes, m->flank_sequence_length))
		return "a transcript sequence slice is outside its pool";
	m->has_wrapped_coordinates = wrapped;
	return NULL;
}

static const char *
snap_restore_names(duckvep_owned_model_t *m, const char *blob, size_t bytes)
{
	size_t region, at = 0, regions = m->known_seq_region_count;

	m->sequence_names = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, regions != 0 ? regions : 1u,
	    sizeof(*m->sequence_names));
	if (m->sequence_names == NULL)
		return "out of memory";
	for (region = 0; region < regions; region++) {
		const char *end;

		if (at >= bytes)
			return "the sequence-region names are truncated";
		if (blob[at++] == 0)
			continue;
		end = memchr(blob + at, '\0', bytes - at);
		if (end == NULL)
			return "the sequence-region names are truncated";
		m->sequence_names[region] = duckvep_budget_strdup(DUCKVEP_OWNER_CONTROL, blob + at);
		if (m->sequence_names[region] == NULL)
			return "out of memory";
		at = (size_t)(end - blob) + 1u;
	}
	return at == bytes ? NULL : "the sequence-region names have trailing bytes";
}

/* Fills *model from the snapshot at `path`. On failure the model is destroyed. */
static int
snap_load(duckvep_owned_model_t *model, const char *path, const char *label, char *error, size_t error_size)
{
	const snap_header_t *header;
	snap_array_t arrays[SNAP_ARRAYS];
	const unsigned char *base;
	const char *problem = NULL, *reference = NULL;
	uint64_t at, checksum;
	snap_header_t zeroed;
	size_t bytes = 0, section;
	int mapped = 0;

	memset(model, 0, sizeof(*model));
	base = snap_open(path, &bytes, &mapped, label, error, error_size);
	if (base == NULL)
		return 0;
	model->snapshot_base = (void *)base;
	model->snapshot_bytes = bytes;
	model->snapshot_mapped = mapped;
	header = (const snap_header_t *)base;
	if (memcmp(header->magic, SNAP_MAGIC, sizeof(header->magic)) != 0)
		problem = "not a DuckVEP model snapshot";
	else if (header->version != SNAP_VERSION || header->endian != SNAP_ENDIAN ||
	    header->header_bytes != sizeof(*header))
		problem = "the snapshot was written by an incompatible DuckVEP build";
	else if (header->file_bytes != bytes)
		problem = "the snapshot is truncated or has trailing bytes";
	else if (header->counts[SNAP_TRANSCRIPTS] > UINT32_MAX || header->counts[SNAP_REGIONS] > UINT16_MAX + 1u ||
	    header->counts[SNAP_EXONS] > UINT32_MAX || header->counts[SNAP_MIRNA] > UINT32_MAX ||
	    header->counts[SNAP_PEPTIDE] > UINT32_MAX || header->counts[SNAP_FEATURES] > UINT32_MAX ||
	    header->counts[SNAP_CDS] > bytes || header->counts[SNAP_FLANK] > bytes)
		problem = "the snapshot header has counts out of range";
	if (problem != NULL)
		goto fail;
	snap_arrays(model, header->counts, arrays);
	at = (sizeof(*header) + SNAP_ALIGN - 1u) / SNAP_ALIGN * SNAP_ALIGN;
	if (at > bytes) {
		problem = "the snapshot is truncated or has trailing bytes";
		goto fail;
	}
	for (section = 0; section < SNAP_SECTIONS; section++) {
		uint64_t expected = header->bytes[section];
		int present = section >= SNAP_ARRAYS || ((header->present >> section) & 1u) != 0;

		if (section < SNAP_ARRAYS) {
			if (header->width[section] != arrays[section].width)
				problem = "the snapshot was written by an incompatible DuckVEP build";
			expected = present ? arrays[section].count * arrays[section].width : 0u;
		} else if (header->width[section] != 1u) {
			problem = "the snapshot was written by an incompatible DuckVEP build";
		}
		if (problem == NULL && (header->bytes[section] != expected || header->offset[section] != at ||
		    header->bytes[section] > bytes - at))
			problem = "the snapshot sections do not match its header";
		if (problem != NULL)
			goto fail;
		at += (header->bytes[section] + SNAP_ALIGN - 1u) / SNAP_ALIGN * SNAP_ALIGN;
		if (at > bytes) {
			problem = "the snapshot sections do not match its header";
			goto fail;
		}
	}
	zeroed = *header;
	zeroed.checksum = 0;
	checksum = snap_hash(SNAP_VERSION, &zeroed, sizeof(zeroed));
	for (section = 0; section < SNAP_SECTIONS; section++)
		checksum = snap_hash(checksum, base + header->offset[section], (size_t)header->bytes[section]);
	if (checksum != header->checksum) {
		problem = "the snapshot checksum does not match its contents";
		goto fail;
	}
	for (section = 0; section < SNAP_ARRAYS; section++) {
		if ((header->present >> section) & 1u)
			*arrays[section].slot = (void *)(base + header->offset[section]);
	}
	/* The relation loaders allocate every region, transcript and exon array; the optional miRNA,
	 * peptide-edit, sequence-byte and feature arrays are checked by the kernel against their counts. */
	for (section = 0; section < 26u; section++) {
		uint64_t rows = header->counts[section < 3u ? SNAP_REGIONS : section < 20u ? SNAP_TRANSCRIPTS : SNAP_EXONS];

		if (rows != 0 && *arrays[section].slot == NULL) {
			problem = "the snapshot lacks a required array";
			goto fail;
		}
	}
	model->known_seq_region_count = (size_t)header->counts[SNAP_REGIONS];
	model->transcripts.transcript_count = (size_t)header->counts[SNAP_TRANSCRIPTS];
	model->exons.exon_count = (size_t)header->counts[SNAP_EXONS];
	model->mature_mirna_count = (size_t)header->counts[SNAP_MIRNA];
	model->peptide_edit_count = (size_t)header->counts[SNAP_PEPTIDE];
	model->cds_sequence_length = (size_t)header->counts[SNAP_CDS];
	model->flank_sequence_length = (size_t)header->counts[SNAP_FLANK];
	model->interval_feature_count = (size_t)header->counts[SNAP_FEATURES];
	model->transcript_coverage_complete = (header->flags & SNAP_FLAG_COVERAGE) != 0;
	model->transcript_flanks_complete = (header->flags & SNAP_FLAG_FLANKS) != 0;
	if ((model->interval_feature_count != 0 && (model->interval_feature_seq_regions == NULL ||
	    model->interval_feature_starts == NULL || model->interval_feature_ends == NULL)) ||
	    (problem = snap_validate(model)) != NULL) {
		if (problem == NULL)
			problem = "the snapshot lacks a required array";
		goto fail;
	}
	if (header->bytes[SNAP_NAMES] != 0 &&
	    (problem = snap_restore_names(model, (const char *)base + header->offset[SNAP_NAMES],
	    (size_t)header->bytes[SNAP_NAMES])) != NULL)
		goto fail;
	if (header->bytes[SNAP_REFERENCE] != 0) {
		reference = (const char *)base + header->offset[SNAP_REFERENCE];
		if (reference[header->bytes[SNAP_REFERENCE] - 1u] != '\0' ||
		    strlen(reference) + 1u != header->bytes[SNAP_REFERENCE]) {
			problem = "the snapshot reference path is malformed";
			goto fail;
		}
		for (section = 0; section < model->known_seq_region_count; section++) {
			if (model->sequence_names == NULL || model->sequence_names[section] == NULL) {
				problem = "the snapshot names a reference FASTA but not every sequence region";
				goto fail;
			}
		}
		/* The model's FASTA is opened and checked against the regions again, as at a relation load. */
		if (!duckvep_core_model_validate_reference(reference, model, error, error_size)) {
			duckvep_owned_model_destroy(model);
			return 0;
		}
	}
	return 1;
fail:
	snap_error(error, error_size, label, problem, path);
	duckvep_owned_model_destroy(model);
	return 0;
}

int
duckvep_core_model_snapshot_install(duckvep_registry_t *registry, const char *name, const char *path,
	char *final_error, size_t final_error_size)
{
	static const char label[] = "duckvep_model_restore";
	duckvep_model_entry_t *entry;
	char error[DUCKVEP_SQL_ERROR_SIZE];
	char message[DUCKVEP_SQL_ERROR_SIZE + 256];

	duckvep_budget_clear_failure();
	pthread_mutex_lock(&registry->mutex);
	entry = duckvep_registry_find_locked(registry, name);
	pthread_mutex_unlock(&registry->mutex);
	if (entry != NULL) {
		(void)snprintf(final_error, final_error_size, "%s: model name already exists", label);
		return 0;
	}
	entry = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1, sizeof(*entry));
	if (entry == NULL || (entry->name = duckvep_budget_strdup(DUCKVEP_OWNER_CONTROL, name)) == NULL) {
		duckvep_model_entry_destroy(entry);
		duckvep_sql_set_error(final_error, final_error_size, duckvep_sql_final_error(message,
		    sizeof(message), NULL, "duckvep_model_restore: out of memory"));
		return 0;
	}
	memset(error, 0, sizeof(error));
	/* Both steps destroy the model themselves on failure. */
	if (!snap_load(&entry->model, path, label, error, sizeof(error)) ||
	    !duckvep_core_model_finish(&entry->model, error, sizeof(error))) {
		duckvep_model_entry_destroy(entry);
		duckvep_sql_set_error(final_error, final_error_size, duckvep_sql_final_error(message,
		    sizeof(message), error, "duckvep_model_restore: model restore failed"));
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
