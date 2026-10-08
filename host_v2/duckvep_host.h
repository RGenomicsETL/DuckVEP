/* The v2 host layer for src/core/duckvep_core_annotate_run.c: the vector operations of the
 * annotation natives over the v2 C API. A call (v2_call) owns small wrappers for the flat
 * argument vectors and for the result vectors it opens; output validity masks start all valid.
 * Failures of an API call are recorded on the call and reported after the core returns. */
#ifndef DUCKVEP_HOST_H
#define DUCKVEP_HOST_H

#include "host_v2_common.h"
#include "core/duckvep_core_model.h"
#include "core/duckvep_core_haplotypes.h"

#define V2_VEC_KIDS 16

typedef struct v2_vec {
    duckdb_v2_vector_handle handle;
    void *data;
    uint64_t *validity;
    size_t size;       /* elements the vector is opened (and its mask sized) for */
    size_t final_size; /* the logical size once the call is done (<= size) */
    size_t child_size; /* a list's final element count */
    size_t child_capacity; /* a list's reserved element count */
    duckdb_v2_arena_handle arena;
    struct v2_call *call;
    bool leaf;                        /* has a data pointer */
    struct v2_vec *kids[V2_VEC_KIDS]; /* the children already opened, by child index */
} v2_vec;

#define V2_CALL_INPUTS 32
/* Flat wrappers serve scalar arguments or table-function result columns. */
#define V2_CALL_COLUMNS ((V2_CALL_INPUTS > DUCKVEP_HAP_OUTPUT_COLUMNS) ? V2_CALL_INPUTS : DUCKVEP_HAP_OUTPUT_COLUMNS)
#define V2_CALL_POOL 192

typedef struct v2_call {
    duckdb_v2_scalar_function_exec_info_handle info;
    duckdb_v2_error_info_handle *error;
    void *user_data;
    size_t rows;
    size_t argc;
    v2_vec inputs[V2_CALL_COLUMNS];
    v2_vec output;
    v2_vec pool[V2_CALL_POOL];
    size_t pool_used;
    bool failed;
} v2_call;

typedef v2_vec *duckvep_h_vector;
typedef v2_call *duckvep_h_chunk;
typedef v2_call *duckvep_h_info;
typedef duckdb_v2_list_entry duckvep_h_list_entry;

/* The registry behind a function's user data (defined in host_v2_model.c). */
duckvep_registry_t *host_v2_registry_of(void *user_data);

static inline void v2_fail(v2_call *call, DUCKDB_V2_ERROR status, duckdb_v2_error_info_handle detail) {
    if (!call->failed) {
        call->failed = true;
        copy_duckdb_error(*call->error, status, detail);
    }
}

/* Sets up `v` as a flat, writable vector of `size` elements whose mask starts all valid. */
static inline bool v2_open_writable(v2_call *call, v2_vec *v, duckdb_v2_vector_handle handle, size_t size,
                                    bool leaf) {
    duckdb_v2_error_info_handle detail = NULL;
    DUCKDB_V2_ERROR status;
    memset(v, 0, sizeof(*v));
    v->handle = handle;
    v->size = size;
    v->final_size = size;
    v->call = call;
    v->leaf = leaf;
    status = duckdb_v2_vector_set_size(handle, size, &detail);
    if (status == DUCKDB_V2_ERROR_NONE && leaf) {
        status = duckdb_v2_vector_get_data_mutable(handle, &v->data, &detail);
    }
    if (status == DUCKDB_V2_ERROR_NONE) {
        status = duckdb_v2_vector_flat_get_validity_mutable(handle, &v->validity, &detail);
    }
    if (status != DUCKDB_V2_ERROR_NONE) {
        v2_fail(call, status, detail);
    } else {
        for (size_t word = 0; word < (size + 63) / 64; ++word) {
            v->validity[word] = ~UINT64_C(0);
        }
    }
    (void)duckdb_v2_error_info_destroy(&detail);
    return !call->failed;
}

/* Grows an opened vector to `size` elements. The elements already written stay; the data and mask
 * pointers are fetched again, and the new mask words start all valid. */
static inline void v2_grow(v2_call *call, v2_vec *v, size_t size) {
    duckdb_v2_error_info_handle detail = NULL;
    size_t old_words = (v->size + 63) / 64;
    DUCKDB_V2_ERROR status = duckdb_v2_vector_set_size(v->handle, size, &detail);
    if (status == DUCKDB_V2_ERROR_NONE && v->leaf) {
        status = duckdb_v2_vector_get_data_mutable(v->handle, &v->data, &detail);
    }
    if (status == DUCKDB_V2_ERROR_NONE) {
        status = duckdb_v2_vector_flat_get_validity_mutable(v->handle, &v->validity, &detail);
    }
    if (status != DUCKDB_V2_ERROR_NONE) {
        v2_fail(call, status, detail);
    } else {
        for (size_t word = old_words; word < (size + 63) / 64; ++word) {
            v->validity[word] = ~UINT64_C(0);
        }
        v->size = size;
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* The child `index` of `parent`, opened for at least `size` elements. A child is opened once per call and
 * grown when a later row needs more room, so several rows of one chunk can share a list. */
static inline v2_vec *v2_child_of(v2_vec *parent, size_t index, size_t size, bool leaf) {
    v2_call *call = parent->call;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle handle = NULL;
    v2_vec *child;
    DUCKDB_V2_ERROR status;
    if (index < V2_VEC_KIDS && parent->kids[index]) {
        child = parent->kids[index];
        if (size > child->size && !call->failed) {
            v2_grow(call, child, size);
        }
        return child;
    }
    if (call->pool_used >= V2_CALL_POOL) {
        v2_fail(call, DUCKDB_V2_ERROR_INPUT_INVALID, NULL);
        return &call->pool[0];
    }
    child = &call->pool[call->pool_used++];
    status = duckdb_v2_vector_get_child(parent->handle, index, &handle, &detail);
    if (status != DUCKDB_V2_ERROR_NONE) {
        v2_fail(call, status, detail);
        memset(child, 0, sizeof(*child));
        (void)duckdb_v2_error_info_destroy(&detail);
        return child;
    }
    if (v2_open_writable(call, child, handle, size, leaf) && index < V2_VEC_KIDS) {
        parent->kids[index] = child;
    }
    (void)duckdb_v2_error_info_destroy(&detail);
    return child;
}

/* The element capacity to open a list's child for: the reserved or extended count, doubled when an
 * opened child has to grow so that a chunk of rows costs a logarithmic number of resizes. */
static inline size_t v2_list_capacity(const v2_vec *list) {
    size_t capacity = list->child_capacity > list->child_size ? list->child_capacity : list->child_size;
    const v2_vec *child = list->kids[0];
    if (child && capacity > child->size && capacity < 2 * child->size) {
        capacity = 2 * child->size;
    }
    return capacity;
}

#define duckvep_h_data(v) ((v)->data)
#define duckvep_h_validity(v) ((v)->validity)
#define duckvep_h_ensure_validity(v) ((void)0)
#define duckvep_h_chunk_vector(c, i) (&(c)->inputs[i])
#define duckvep_h_chunk_rows(c) ((c)->rows)
#define duckvep_h_chunk_columns(c) ((c)->argc)
#define duckvep_h_set_valid(v, row) mark_valid((v)->validity, (row))
#define duckvep_h_set_null(v, row) mark_null((v)->validity, (row))
#define duckvep_h_set_error(info, message) set_error(*(info)->error, DUCKDB_V2_ERROR_INPUT_INVALID, (message))
#define duckvep_h_extra_info(info) host_v2_registry_of((info)->user_data)
#define duckvep_h_vector_size() ((size_t)2048)

/* The record vector of a result list: its fields are opened by duckvep_h_struct_child. The list is
 * opened for its reserved capacity (writers size masks by an upper bound) and shrunk to its final
 * element count when the call ends (v2_finish). */
static inline v2_vec *duckvep_h_list_child(v2_vec *list) {
    v2_vec *child = v2_child_of(list, 0, v2_list_capacity(list), false);
    child->final_size = list->child_size;
    return child;
}

/* The element vector of a list of primitives or text: a leaf, so it has data. */
static inline v2_vec *duckvep_h_list_values(v2_vec *list) {
    v2_vec *child = v2_child_of(list, 0, v2_list_capacity(list), true);
    child->final_size = list->child_size;
    return child;
}

static inline v2_vec *duckvep_h_struct_child(v2_vec *record, size_t index) {
    v2_vec *field = v2_child_of(record, index, record->size, true);
    field->final_size = record->final_size;
    return field;
}

/* Shrinks the vectors opened for a reserved capacity to their final sizes. */
static inline void v2_finish(v2_call *call) {
    for (size_t i = call->pool_used; i-- > 0;) {
        v2_vec *v = &call->pool[i];
        if (v->handle && v->final_size != v->size) {
            duckdb_v2_error_info_handle detail = NULL;
            DUCKDB_V2_ERROR status = duckdb_v2_vector_set_size(v->handle, v->final_size, &detail);
            if (status != DUCKDB_V2_ERROR_NONE) {
                v2_fail(call, status, detail);
            }
            (void)duckdb_v2_error_info_destroy(&detail);
        }
    }
}

/* A list's elements are sized once, by the final count (set_size reserves). */
static inline int duckvep_h_list_reserve(v2_vec *list, size_t count) {
    list->child_capacity = count;
    return 1;
}

static inline int duckvep_h_list_set_size(v2_vec *list, size_t count) {
    list->child_size = count;
    return 1;
}

/* Appends `count` elements to a result list and returns where they start. The child is opened, or grown,
 * by the next duckvep_h_list_child / duckvep_h_list_values call. */
static inline int duckvep_h_list_extend(v2_vec *list, size_t count, duckdb_v2_list_entry *entry) {
    if (count > SIZE_MAX - list->child_size) return 0;
    entry->offset = list->child_size;
    entry->length = count;
    list->child_size += count;
    if (list->child_capacity < list->child_size) list->child_capacity = list->child_size;
    return 1;
}

/* Text is written through the vector's arena; invalid UTF-8 becomes NULL, as in v1. */
static inline void duckvep_h_assign_string(v2_vec *v, size_t row, const char *text, size_t length) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_str view = {text, (idx_t)length};
    DUCKDB_V2_ERROR status = duckdb_v2_validate_utf8(view, NULL);
    if (status != DUCKDB_V2_ERROR_NONE) {
        v->validity[row >> 6] &= ~(UINT64_C(1) << (row & 63));
        return;
    }
    if (!v->arena) {
        status = duckdb_v2_vector_get_arena(v->handle, &v->arena, &detail);
        if (status != DUCKDB_V2_ERROR_NONE) {
            v2_fail(v->call, status, detail);
            (void)duckdb_v2_error_info_destroy(&detail);
            return;
        }
    }
    status = write_string(v->arena, &((duckdb_v2_bytes *)v->data)[row], text, length, &detail);
    if (status != DUCKDB_V2_ERROR_NONE) {
        v2_fail(v->call, status, detail);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

static inline char *duckvep_h_string(v2_vec *v, size_t row) {
    duckvep_col_t column;
    column.data = v->data;
    column.validity = v->validity;
    return duckvep_col_string(&column, row);
}

static inline int duckvep_h_string_wellformed(v2_vec *v, size_t row) {
    duckvep_col_t column;
    column.data = v->data;
    column.validity = v->validity;
    return duckvep_col_string_wellformed(&column, row);
}

#endif
