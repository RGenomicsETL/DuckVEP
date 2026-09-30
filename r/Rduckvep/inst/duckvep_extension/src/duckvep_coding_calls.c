/* duckvep_coding_calls(model, path): a fused native VCF/BCF reader that returns the calls relation consumed by
 * duckvep_haplotypes, decoding genotypes only for the records that touch coding sequence.
 *
 * Every record is read as text (or as a BCF record), its CHROM is mapped to the model's seq_region by name, and each ALT
 * allele goes through the same discovery as the scalar duckvep_coding_transcripts (duckvep_discovery.c). Only when an
 * allele overlaps the CDS of at least one transcript does the reader parse FORMAT and decode GT and PS; the other records
 * (99% of a genome) cost one line read, four field splits and a few interval lookups.
 *
 * Output schema and semantics are those of the calls relation of the mode B scale worker: one row per (ALT allele,
 * transcript, sample) with event_index = (record ordinal << 6) | (alt ordinal - 1), a 1-based record ordinal counting every
 * data record of the file, alt_index 1-based, alleles INTEGER[] (NULL for a missing allele), phase_before BOOLEAN[] with a
 * leading false and then one flag per separator ('|' is true), and phase_set the integer PS of the sample (labels such as
 * PATMAT, '.' and absent values are NULL). htslib (bundled, zlib only) provides the BGZF, VCF and BCF decoding. */
#include "duckdb_extension.h"
#include "kernel/src/duckvep_budget.h"
DUCKDB_EXTENSION_EXTERN
#include "duckvep_list.h"

#include "duckvep_discovery.h"
#include "duckvep_model.h"

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>

#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CC "duckvep_coding_calls: "
enum { CC_MAX_ALTS = 64 };
/* Fixed reservation for htslib's own buffers (BGZF block, block cache, line buffer), which cannot be routed through
 * the budget; the header and the growing buffers add to it while the file is read. */
#define CC_BASE_RESERVE ((uint64_t)2u << 20)

typedef struct { const char *name; uint32_t length; uint32_t region; } cc_region_t;

typedef struct {
    duckvep_registry_t *registry;
    duckvep_model_entry_t *entry;
    char *path;
    cc_region_t *regions;   /* sorted by name */
    size_t region_count;
} cc_bind_t;

typedef struct { uint32_t transcript, alt_index, alt_offset, alt_length; } cc_pair_t;

typedef struct {
    const cc_bind_t *bind;
    htsFile *fp;
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    kstring_t line;
    int binary, samples, gt_id, ps_id, ps_type;
    uint64_t reserved, line_reserved;
    uint64_t record_index;
    duckvep_discovery_t scratch;
    duckvep_u32_list_t found;
    /* The record being emitted. */
    int have_record;
    uint32_t region;
    int64_t position;
    char *text; size_t text_length, text_capacity;   /* REF at 0, then the ALT alleles that produced pairs */
    uint32_t ref_length;
    cc_pair_t *pairs; size_t pair_count, pair_capacity;
    int32_t *alleles; size_t allele_capacity;         /* samples x lanes, -1 for a missing allele */
    uint8_t *phase; size_t phase_capacity;
    uint16_t *ploidy; size_t ploidy_capacity;         /* lanes per sample */
    int64_t *phase_set; size_t phase_set_capacity; uint8_t *phase_set_present; size_t present_capacity;
    uint32_t lanes;
    size_t pair_cursor; uint32_t sample_cursor;
    /* Last CHROM text -> region (text VCF), last rid -> region (BCF). */
    char last_name[64]; size_t last_name_length; int last_found; uint32_t last_region;
    int last_rid;
    int eof;
} cc_state_t;

static void cc_bind_destroy(void *pointer) {
    cc_bind_t *b = pointer;
    if (!b) return;
    duckvep_registry_unpin(b->registry, b->entry);
    duckvep_registry_release(b->registry);
    duckdb_free(b->path);
    duckvep_budget_free(b->regions);
    duckvep_budget_free(b);
}

static int cc_region_compare(const void *left, const void *right) {
    const cc_region_t *a = left, *b = right;
    return strcmp(a->name, b->name);
}

static void cc_add_column(duckdb_bind_info info, const char *name, duckdb_type type, int list) {
    duckdb_logical_type element = duckdb_create_logical_type(type), result = element;
    if (list) result = duckdb_create_list_type(element);
    duckdb_bind_add_result_column(info, name, result);
    if (list) duckdb_destroy_logical_type(&result);
    duckdb_destroy_logical_type(&element);
}

static void cc_bind(duckdb_bind_info info) {
    cc_bind_t *b = duckvep_budget_calloc(DUCKVEP_OWNER_CONTROL, 1u, sizeof(*b));
    if (!b) { duckdb_bind_set_error(info, CC "bind allocation failed"); return; }
    b->registry = duckdb_bind_get_extra_info(info);
    duckvep_registry_retain(b->registry);
    duckdb_value value = duckdb_bind_get_parameter(info, 0u);
    char *name = value && !duckdb_is_null_value(value) ? duckdb_get_varchar(value) : NULL;
    duckdb_destroy_value(&value);
    value = duckdb_bind_get_parameter(info, 1u);
    if (value && !duckdb_is_null_value(value)) b->path = duckdb_get_varchar(value);
    duckdb_destroy_value(&value);
    if (name) b->entry = duckvep_registry_pin(b->registry, name);
    duckdb_free(name);
    if (!b->entry || !b->path || !b->path[0]) {
        duckdb_bind_set_error(info, CC "require a loaded model name and a nonempty VCF or BCF path");
        cc_bind_destroy(b); return;
    }
    const duckvep_owned_model_t *m = &b->entry->model;
    if (m->lifted) {
        duckdb_bind_set_error(info, CC "transcript discovery is not supported for models with wrapped circular objects");
        cc_bind_destroy(b); return;
    }
    if (!m->sequence_names || !m->known_seq_region_count) {
        duckdb_bind_set_error(info, CC "the model has no seq_region_name values; load it with seq_region_name in the regions query");
        cc_bind_destroy(b); return;
    }
    b->regions = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, m->known_seq_region_count * sizeof(*b->regions));
    if (!b->regions) { duckdb_bind_set_error(info, CC "bind allocation failed"); cc_bind_destroy(b); return; }
    for (size_t i = 0u; i < m->known_seq_region_count; i++) {
        if (!m->sequence_names[i]) continue;
        b->regions[b->region_count].name = m->sequence_names[i];
        b->regions[b->region_count].length = (uint32_t)strlen(m->sequence_names[i]);
        b->regions[b->region_count].region = m->known_seq_regions[i];
        b->region_count++;
    }
    if (!b->region_count) {
        duckdb_bind_set_error(info, CC "the model has no seq_region_name values; load it with seq_region_name in the regions query");
        cc_bind_destroy(b); return;
    }
    qsort(b->regions, b->region_count, sizeof(*b->regions), cc_region_compare);
    for (size_t i = 1u; i < b->region_count; i++)
        if (!strcmp(b->regions[i - 1u].name, b->regions[i].name)) {
            duckdb_bind_set_error(info, CC "seq_region_name values must be unique");
            cc_bind_destroy(b); return;
        }
    cc_add_column(info, "event_index", DUCKDB_TYPE_BIGINT, 0);
    cc_add_column(info, "seq_region", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "position", DUCKDB_TYPE_BIGINT, 0);
    cc_add_column(info, "reference", DUCKDB_TYPE_VARCHAR, 0);
    cc_add_column(info, "alternate", DUCKDB_TYPE_VARCHAR, 0);
    cc_add_column(info, "alt_index", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "transcript_index", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "sample_index", DUCKDB_TYPE_INTEGER, 0);
    cc_add_column(info, "alleles", DUCKDB_TYPE_INTEGER, 1);
    cc_add_column(info, "phase_before", DUCKDB_TYPE_BOOLEAN, 1);
    cc_add_column(info, "phase_set", DUCKDB_TYPE_BIGINT, 0);
    duckdb_bind_set_bind_data(info, b, cc_bind_destroy);
}

static void cc_state_destroy(void *pointer) {
    cc_state_t *s = pointer;
    if (!s) return;
    if (s->rec) bcf_destroy(s->rec);
    if (s->hdr) bcf_hdr_destroy(s->hdr);
    if (s->fp) hts_close(s->fp);
    free(s->line.s);
    if (s->reserved) duckvep_budget_unreserve(DUCKVEP_OWNER_WORKSPACE, s->reserved);
    duckvep_discovery_release(&s->scratch);
    duckvep_u32_list_release(&s->found);
    duckvep_budget_free(s->text); duckvep_budget_free(s->pairs); duckvep_budget_free(s->alleles);
    duckvep_budget_free(s->phase); duckvep_budget_free(s->ploidy); duckvep_budget_free(s->phase_set);
    duckvep_budget_free(s->phase_set_present);
    duckvep_budget_free(s);
}

static int cc_reserve(cc_state_t *s, uint64_t bytes) {
    if (!duckvep_budget_reserve(DUCKVEP_OWNER_WORKSPACE, bytes)) return 0;
    s->reserved += bytes;
    return 1;
}

/* Grow a budget-owned array to at least `need` elements. */
static int cc_grow(void **array, size_t *capacity, size_t need, size_t width) {
    if (need <= *capacity) return 1;
    size_t next = *capacity ? *capacity : 64u;
    while (next < need) {
        if (next > SIZE_MAX / 2u) return 0;
        next *= 2u;
    }
    if (next > SIZE_MAX / width) return 0;
    void *grown = duckvep_budget_realloc(DUCKVEP_OWNER_WORKSPACE, *array, next * width);
    if (!grown) return 0;
    *array = grown; *capacity = next;
    return 1;
}

static int cc_open(cc_state_t *s, char *error, size_t error_size) {
    if (!cc_reserve(s, CC_BASE_RESERVE)) goto nomem;
    int level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_OFF);   /* failures get explicit messages below */
    s->fp = hts_open(s->bind->path, "r");
    hts_set_log_level(level);
    if (!s->fp) {
        snprintf(error, error_size, CC "cannot open '%s'%s%s", s->bind->path, errno ? ": " : "", errno ? strerror(errno) : "");
        return 0;
    }
    const htsFormat *format = hts_get_format(s->fp);
    if (!format || (format->format != vcf && format->format != bcf)) {
        snprintf(error, error_size, CC "'%s' is not a VCF or BCF file", s->bind->path);
        return 0;
    }
    s->binary = format->format == bcf;
    {   /* the header of phased callsets often declares PS as a String: a label, not an error */
        hts_set_log_level(HTS_LOG_ERROR);
        s->hdr = bcf_hdr_read(s->fp);
        hts_set_log_level(level);
    }
    if (!s->hdr) { snprintf(error, error_size, CC "cannot read the header of '%s'", s->bind->path); return 0; }
    s->samples = bcf_hdr_nsamples(s->hdr);
    if (s->samples < 1) {
        snprintf(error, error_size, CC "'%s' has no samples, so there are no genotypes to call", s->bind->path);
        return 0;
    }
    s->gt_id = bcf_hdr_id2int(s->hdr, BCF_DT_ID, "GT");
    if (!bcf_hdr_idinfo_exists(s->hdr, BCF_HL_FMT, s->gt_id)) {
        snprintf(error, error_size, CC "the header of '%s' declares no FORMAT/GT field", s->bind->path);
        return 0;
    }
    s->ps_id = bcf_hdr_id2int(s->hdr, BCF_DT_ID, "PS");
    s->ps_type = -1;
    if (bcf_hdr_idinfo_exists(s->hdr, BCF_HL_FMT, s->ps_id)) s->ps_type = (int)bcf_hdr_id2type(s->hdr, BCF_HL_FMT, s->ps_id);
    else s->ps_id = -1;
    if (!cc_reserve(s, 2048u * (uint64_t)s->hdr->nhrec + 256u * (uint64_t)s->samples)) goto nomem;
    s->rec = bcf_init();
    if (!s->rec) goto nomem;
    return 1;
nomem:
    if (!duckvep_budget_take_failure(error, error_size)) snprintf(error, error_size, CC "out of memory or native budget exceeded");
    return 0;
}

static void cc_init(duckdb_init_info info) {
    const cc_bind_t *bind = duckdb_init_get_bind_data(info);
    cc_state_t *s = duckvep_budget_calloc(DUCKVEP_OWNER_WORKSPACE, 1u, sizeof(*s));
    char error[DUCKVEP_SQL_ERROR_SIZE] = CC "allocation failed or native budget exceeded";
    duckdb_init_set_max_threads(info, 1u);
    duckvep_budget_clear_failure();
    if (!s) goto failed;
    s->bind = bind;
    s->last_rid = -1;
    duckvep_discovery_init(&s->scratch);
    if (!cc_open(s, error, sizeof error)) goto failed;
    duckdb_init_set_init_data(info, s, cc_state_destroy);
    return;
failed:
    duckdb_init_set_error(info, error);
    cc_state_destroy(s);
}

/* Name -> seq_region, by binary search over the sorted names. */
static int cc_lookup(const cc_bind_t *b, const char *name, size_t length, uint32_t *region) {
    size_t lo = 0u, hi = b->region_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        const cc_region_t *r = &b->regions[mid];
        size_t common = length < r->length ? length : r->length;
        int c = memcmp(name, r->name, common);
        if (!c) c = length < r->length ? -1 : length > r->length;
        if (!c) { *region = r->region; return 1; }
        if (c < 0) hi = mid; else lo = mid + 1u;
    }
    return 0;
}

/* Integer lane `k` of sample `i` of a BCF-packed FORMAT field. Returns 0 for a value, 1 missing, 2 past the vector end. */
static int cc_fmt_int(const bcf_fmt_t *f, int sample, int k, int32_t *out) {
    const uint8_t *p = f->p + (size_t)sample * (size_t)f->size;
    switch (f->type) {
    case BCF_BT_INT8: { int8_t v = ((const int8_t *)p)[k];
        if (v == bcf_int8_vector_end) return 2;
        if (v == bcf_int8_missing) return 1;
        *out = v; return 0; }
    case BCF_BT_INT16: { int16_t v; memcpy(&v, p + 2 * k, 2);
        if (v == bcf_int16_vector_end) return 2;
        if (v == bcf_int16_missing) return 1;
        *out = v; return 0; }
    case BCF_BT_INT32: { int32_t v; memcpy(&v, p + 4 * k, 4);
        if (v == bcf_int32_vector_end) return 2;
        if (v == bcf_int32_missing) return 1;
        *out = v; return 0; }
    default: return -1;
    }
}

/* try_cast(text AS BIGINT) for a PS label: a whole-string decimal integer, else none. */
static int cc_parse_integer(const char *text, size_t length, int64_t *out) {
    size_t i = 0u;
    int negative = 0;
    if (i < length && (text[i] == '+' || text[i] == '-')) negative = text[i++] == '-';
    if (i == length) return 0;
    uint64_t value = 0u;
    for (; i < length; i++) {
        if (text[i] < '0' || text[i] > '9') return 0;
        unsigned digit = (unsigned)(text[i] - '0');
        if (value > (UINT64_C(9223372036854775807) - digit) / 10u) return 0;
        value = value * 10u + digit;
    }
    *out = negative ? -(int64_t)value : (int64_t)value;
    return 1;
}

/* Decode GT and PS of every sample of the parsed record. */
static int cc_decode_samples(cc_state_t *s, char *error, size_t error_size) {
    bcf1_t *rec = s->rec;
    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t *gt = bcf_get_fmt(s->hdr, rec, "GT");
    if (!gt || gt->n < 1 || gt->type == BCF_BT_NULL) {
        snprintf(error, error_size, CC "record %llu at position %lld has no GT; genotypes are required for every record that touches coding sequence",
            (unsigned long long)s->record_index, (long long)s->position);
        return 0;
    }
    uint32_t lanes = (uint32_t)gt->n;
    size_t cells = (size_t)s->samples * lanes;
    if (lanes > 1024u || cells / lanes != (size_t)s->samples) {
        snprintf(error, error_size, CC "record %llu has unsupported ploidy %u", (unsigned long long)s->record_index, lanes);
        return 0;
    }
    if (!cc_grow((void **)&s->alleles, &s->allele_capacity, cells, sizeof *s->alleles) ||
        !cc_grow((void **)&s->phase, &s->phase_capacity, cells, sizeof *s->phase) ||
        !cc_grow((void **)&s->ploidy, &s->ploidy_capacity, (size_t)s->samples, sizeof *s->ploidy) ||
        !cc_grow((void **)&s->phase_set, &s->phase_set_capacity, (size_t)s->samples, sizeof *s->phase_set) ||
        !cc_grow((void **)&s->phase_set_present, &s->present_capacity, (size_t)s->samples, sizeof *s->phase_set_present)) {
        if (!duckvep_budget_take_failure(error, error_size)) snprintf(error, error_size, CC "out of memory decoding genotypes");
        return 0;
    }
    s->lanes = lanes;
    for (int i = 0; i < s->samples; i++) {
        uint32_t n = 0u;
        for (uint32_t k = 0u; k < lanes; k++) {
            /* The packed GT value is ((allele + 1) << 1) | phased; 0 is a missing allele. Read it raw so that the
             * phase bit of a missing allele survives. */
            const uint8_t *p = gt->p + (size_t)i * (size_t)gt->size;
            int32_t raw;
            if (gt->type == BCF_BT_INT8) { int8_t x = ((const int8_t *)p)[k]; if (x == bcf_int8_vector_end) break; raw = x; }
            else if (gt->type == BCF_BT_INT16) { int16_t x; memcpy(&x, p + 2u * k, 2u); if (x == bcf_int16_vector_end) break; raw = x; }
            else if (gt->type == BCF_BT_INT32) { int32_t x; memcpy(&x, p + 4u * k, 4u); if (x == bcf_int32_vector_end) break; raw = x; }
            else {
                snprintf(error, error_size, CC "record %llu has a GT field that is not integer-coded", (unsigned long long)s->record_index);
                return 0;
            }
            s->alleles[(size_t)i * lanes + k] = (raw >> 1) ? (raw >> 1) - 1 : -1;
            s->phase[(size_t)i * lanes + k] = k ? (uint8_t)(raw & 1) : 0u;
            n++;
        }
        if (!n) {   /* a sample with no GT values at all: one missing allele, as the text '.' would give */
            s->alleles[(size_t)i * lanes] = -1; s->phase[(size_t)i * lanes] = 0u; n = 1u;
        }
        s->ploidy[i] = (uint16_t)n;
        s->phase_set_present[i] = 0u;
    }
    if (s->ps_id >= 0) {
        bcf_fmt_t *ps = bcf_get_fmt(s->hdr, rec, "PS");
        if (ps && ps->n >= 1) {
            for (int i = 0; i < s->samples; i++) {
                if (s->ps_type == BCF_HT_STR) {
                    const char *text = (const char *)ps->p + (size_t)i * (size_t)ps->size;
                    size_t length = 0u;
                    while (length < (size_t)ps->n && text[length]) length++;
                    int64_t value;
                    if (cc_parse_integer(text, length, &value)) { s->phase_set[i] = value; s->phase_set_present[i] = 1u; }
                } else if (s->ps_type == BCF_HT_INT) {
                    int32_t value;
                    if (cc_fmt_int(ps, i, 0, &value) == 0) { s->phase_set[i] = value; s->phase_set_present[i] = 1u; }
                }
            }
        }
    }
    return 1;
}

/* Copies `length` bytes of text into the record buffer and returns the offset, or UINT32_MAX. */
static uint32_t cc_store(cc_state_t *s, const char *text, size_t length) {
    if (length > UINT32_MAX / 2u || s->text_length > UINT32_MAX / 2u - length ||
        !cc_grow((void **)&s->text, &s->text_capacity, s->text_length + length + 1u, 1u)) return UINT32_MAX;
    uint32_t offset = (uint32_t)s->text_length;
    memcpy(s->text + offset, text, length);
    s->text_length += length;
    return offset;
}

/* Discovery for one allele; appends pairs. Returns 1 on success. */
static int cc_allele(cc_state_t *s, const duckvep_owned_model_t *m, uint32_t alt_index, const char *alt, size_t alt_length,
    char *error, size_t error_size) {
    size_t added = 0u;
    s->found.count = 0u;
    duckvep_discovery_status_t status = duckvep_discover_coding(m, &s->scratch, (int64_t)s->region, s->position,
        (const uint8_t *)s->text, s->ref_length, (const uint8_t *)alt, alt_length, &s->found, &added);
    if (status == DUCKVEP_DISCOVERY_BAD_SPAN) {
        snprintf(error, error_size, CC "record %llu: position must be from 1 through 2147483647 and alleles at most 65535 bases",
            (unsigned long long)s->record_index);
        return 0;
    }
    if (status == DUCKVEP_DISCOVERY_NOMEM) {
        if (!duckvep_budget_take_failure(error, error_size)) snprintf(error, error_size, CC "out of memory collecting transcripts");
        return 0;
    }
    if (!added) return 1;
    if (alt_index > CC_MAX_ALTS) {
        snprintf(error, error_size, CC "more than 64 ALT alleles in one record (record %llu)", (unsigned long long)s->record_index);
        return 0;
    }
    uint32_t offset = cc_store(s, alt, alt_length);
    if (offset == UINT32_MAX || !cc_grow((void **)&s->pairs, &s->pair_capacity, s->pair_count + added, sizeof *s->pairs)) {
        if (!duckvep_budget_take_failure(error, error_size)) snprintf(error, error_size, CC "out of memory collecting calls");
        return 0;
    }
    for (size_t i = 0u; i < added; i++)
        s->pairs[s->pair_count++] = (cc_pair_t){s->found.items[i], alt_index, offset, (uint32_t)alt_length};
    return 1;
}

/* Reads records until one touches coding sequence (then decodes its samples and returns 1), or the file ends (0), or an
 * error occurs (-1). */
static int cc_next_record(cc_state_t *s, char *error, size_t error_size) {
    const duckvep_owned_model_t *m = &s->bind->entry->model;
    for (;;) {
        s->text_length = 0u; s->pair_count = 0u;
        const char *chrom = NULL, *alts = NULL, *ref = NULL;
        size_t chrom_length = 0u, alts_length = 0u, ref_length = 0u;
        uint32_t region = 0u;
        int found = 0;
        if (s->binary) {
            int r = bcf_read(s->fp, s->hdr, s->rec);
            if (r == -1) return 0;
            if (r < 0) { snprintf(error, error_size, CC "malformed BCF record after record %llu", (unsigned long long)s->record_index); return -1; }
            s->record_index++;
            if (s->rec->rid != s->last_rid) {
                const char *name = bcf_hdr_id2name(s->hdr, s->rec->rid);
                s->last_found = name && cc_lookup(s->bind, name, strlen(name), &s->last_region);
                s->last_rid = s->rec->rid;
            }
            found = s->last_found; region = s->last_region;
            if (!found) continue;
            bcf_unpack(s->rec, BCF_UN_STR);
            if (s->rec->n_allele < 2) continue;   /* no ALT */
            s->position = (int64_t)s->rec->pos + 1;
            ref = s->rec->d.allele[0]; ref_length = strlen(ref);
        } else {
            int r = hts_getline(s->fp, '\n', &s->line);
            if (r == -1) return 0;
            if (r < -1) { snprintf(error, error_size, CC "cannot read the file after record %llu", (unsigned long long)s->record_index); return -1; }
            if (s->line.l == 0u) continue;
            if ((size_t)s->line.m > s->line_reserved) {
                uint64_t more = (uint64_t)s->line.m - s->line_reserved;
                if (!cc_reserve(s, more)) {
                    if (!duckvep_budget_take_failure(error, error_size)) snprintf(error, error_size, CC "native budget exceeded by a record line");
                    return -1;
                }
                s->line_reserved = (uint64_t)s->line.m;
            }
            s->record_index++;
            /* CHROM \t POS \t ID \t REF \t ALT */
            const char *f[5] = {s->line.s, NULL, NULL, NULL, NULL}, *end = s->line.s + s->line.l;
            size_t n = 0u;
            for (const char *p = s->line.s; n < 4u && p < end; p++) {
                p = memchr(p, '\t', (size_t)(end - p));
                if (!p) break;
                f[++n] = p + 1;
            }
            if (n < 4u) {
                snprintf(error, error_size, CC "record %llu has fewer than 5 tab-separated fields", (unsigned long long)s->record_index);
                return -1;
            }
            chrom = f[0]; chrom_length = (size_t)(f[1] - f[0]) - 1u;
            if (chrom_length == s->last_name_length && chrom_length <= sizeof s->last_name && s->last_name_length &&
                !memcmp(chrom, s->last_name, chrom_length)) { found = s->last_found; region = s->last_region; }
            else {
                found = cc_lookup(s->bind, chrom, chrom_length, &region);
                if (chrom_length <= sizeof s->last_name) {
                    memcpy(s->last_name, chrom, chrom_length);
                    s->last_name_length = chrom_length; s->last_found = found; s->last_region = region;
                }
            }
            if (!found) continue;
            ref = f[3]; ref_length = (size_t)(f[4] - f[3]) - 1u;
            alts = f[4];
            const char *alts_end = memchr(alts, '\t', (size_t)(end - alts));
            alts_length = (size_t)((alts_end ? alts_end : end) - alts);
            if (alts_length == 1u && alts[0] == '.') continue;
            /* POS */
            const char *pos = f[1], *pos_end = f[2] - 1;
            int64_t position = 0;
            if (pos == pos_end || pos_end - pos > 18) position = -1;
            for (const char *p = pos; p < pos_end && position >= 0; p++) {
                if (*p < '0' || *p > '9') position = -1; else position = position * 10 + (*p - '0');
            }
            if (position < 0) {
                snprintf(error, error_size, CC "record %llu has an invalid POS", (unsigned long long)s->record_index);
                return -1;
            }
            s->position = position;
        }
        s->region = region;
        if (ref_length > UINT32_MAX / 4u) { snprintf(error, error_size, CC "REF too long at record %llu", (unsigned long long)s->record_index); return -1; }
        s->ref_length = (uint32_t)ref_length;
        if (cc_store(s, ref, ref_length) == UINT32_MAX) goto nomem;
        if (s->binary) {
            for (int a = 1; a < s->rec->n_allele; a++) {
                const char *alt = s->rec->d.allele[a];
                if (!cc_allele(s, m, (uint32_t)a, alt, strlen(alt), error, error_size)) return -1;
            }
        } else {
            uint32_t alt_index = 0u;
            for (size_t at = 0u; at <= alts_length;) {
                const char *comma = memchr(alts + at, ',', alts_length - at);
                size_t length = comma ? (size_t)(comma - (alts + at)) : alts_length - at;
                if (!cc_allele(s, m, ++alt_index, alts + at, length, error, error_size)) return -1;
                if (!comma) break;
                at += length + 1u;
            }
        }
        if (!s->pair_count) continue;
        /* The record touches coding sequence: only now parse FORMAT and decode the genotypes. */
        if (!s->binary) {
            if (vcf_parse(&s->line, s->hdr, s->rec) < 0) {
                snprintf(error, error_size, CC "malformed VCF record %llu at position %lld", (unsigned long long)s->record_index,
                    (long long)s->position);
                return -1;
            }
        }
        if (!cc_decode_samples(s, error, error_size)) return -1;
        s->pair_cursor = 0u; s->sample_cursor = 0u;
        return 1;
    }
nomem:
    if (!duckvep_budget_take_failure(error, error_size)) snprintf(error, error_size, CC "out of memory reading a record");
    return -1;
}

static void cc_scan_chunk(duckdb_function_info info, duckdb_data_chunk output) {
    cc_state_t *s = duckdb_function_get_init_data(info);
    char error[DUCKVEP_SQL_ERROR_SIZE] = {0};
    idx_t capacity = duckdb_vector_size(), rows = 0u;
    duckdb_vector v[11];
    for (unsigned i = 0u; i < 11u; i++) {
        v[i] = duckdb_data_chunk_get_vector(output, i);
        duckdb_vector_ensure_validity_writable(v[i]);
    }
    if (duckdb_list_vector_set_size(v[8], 0u) != DuckDBSuccess || duckdb_list_vector_set_size(v[9], 0u) != DuckDBSuccess) {
        duckdb_function_set_error(info, CC "cannot reset output list"); return;
    }
    while (rows < capacity) {
        if (!s->have_record || s->pair_cursor >= s->pair_count) {
            s->have_record = 0;
            if (s->eof) break;
            int r = cc_next_record(s, error, sizeof error);
            if (r < 0) { duckdb_function_set_error(info, error); return; }
            if (r == 0) { s->eof = 1; break; }
            s->have_record = 1;
        }
        const cc_pair_t *pair = &s->pairs[s->pair_cursor];
        uint32_t sample = s->sample_cursor;
        uint32_t lanes = s->ploidy[sample];
        ((int64_t *)duckdb_vector_get_data(v[0]))[rows] = (int64_t)((s->record_index << 6) | (pair->alt_index - 1u));
        ((int32_t *)duckdb_vector_get_data(v[1]))[rows] = (int32_t)s->region;
        ((int64_t *)duckdb_vector_get_data(v[2]))[rows] = s->position;
        duckdb_vector_assign_string_element_len(v[3], rows, s->text, s->ref_length);
        duckdb_vector_assign_string_element_len(v[4], rows, s->text + pair->alt_offset, pair->alt_length);
        ((int32_t *)duckdb_vector_get_data(v[5]))[rows] = (int32_t)pair->alt_index;
        ((int32_t *)duckdb_vector_get_data(v[6]))[rows] = (int32_t)pair->transcript;
        ((int32_t *)duckdb_vector_get_data(v[7]))[rows] = (int32_t)sample;
        duckdb_list_entry *alleles_entry = (duckdb_list_entry *)duckdb_vector_get_data(v[8]) + rows;
        duckdb_list_entry *phase_entry = (duckdb_list_entry *)duckdb_vector_get_data(v[9]) + rows;
        if (!duckvep_list_extend(v[8], lanes, alleles_entry) || !duckvep_list_extend(v[9], lanes, phase_entry)) {
            duckdb_function_set_error(info, CC "output list allocation failed"); return;
        }
        duckdb_vector allele_child = duckdb_list_vector_get_child(v[8]), phase_child = duckdb_list_vector_get_child(v[9]);
        duckdb_vector_ensure_validity_writable(allele_child);
        uint64_t *allele_validity = duckdb_vector_get_validity(allele_child);
        int32_t *allele_data = duckdb_vector_get_data(allele_child);
        bool *phase_data = duckdb_vector_get_data(phase_child);
        for (uint32_t k = 0u; k < lanes; k++) {
            int32_t a = s->alleles[(size_t)sample * s->lanes + k];
            if (a < 0) { allele_data[alleles_entry->offset + k] = 0; duckdb_validity_set_row_invalid(allele_validity, alleles_entry->offset + k); }
            else { allele_data[alleles_entry->offset + k] = a; duckdb_validity_set_row_valid(allele_validity, alleles_entry->offset + k); }
            phase_data[phase_entry->offset + k] = s->phase[(size_t)sample * s->lanes + k] != 0u;
        }
        if (s->phase_set_present[sample]) {
            ((int64_t *)duckdb_vector_get_data(v[10]))[rows] = s->phase_set[sample];
            duckdb_validity_set_row_valid(duckdb_vector_get_validity(v[10]), rows);
        } else duckdb_validity_set_row_invalid(duckdb_vector_get_validity(v[10]), rows);
        rows++;
        if (++s->sample_cursor >= (uint32_t)s->samples) { s->sample_cursor = 0u; s->pair_cursor++; }
    }
    duckdb_data_chunk_set_size(output, rows);
}

static void cc_scan(duckdb_function_info info, duckdb_data_chunk output) {
    int level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_ERROR);   /* htslib's own notices would repeat on stderr; errors get explicit messages */
    cc_scan_chunk(info, output);
    hts_set_log_level(level);
}

void duckvep_register_coding_calls(duckdb_connection connection, duckvep_registry_t *registry) {
    duckdb_table_function function = duckdb_create_table_function();
    duckdb_logical_type string = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_set_name(function, "duckvep_coding_calls");
    duckdb_table_function_add_parameter(function, string);
    duckdb_table_function_add_parameter(function, string);
    duckvep_registry_retain(registry);
    duckdb_table_function_set_extra_info(function, registry, duckvep_registry_release);
    duckdb_table_function_set_bind(function, cc_bind);
    duckdb_table_function_set_init(function, cc_init);
    duckdb_table_function_set_function(function, cc_scan);
    (void)duckdb_register_table_function(connection, function);
    duckdb_destroy_table_function(&function);
    duckdb_destroy_logical_type(&string);
}
