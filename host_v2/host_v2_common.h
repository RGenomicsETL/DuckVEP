/* Shared DuckDB-facing helpers of the v2 host: errors, vector views, validity,
 * strings, struct results, function registration. Included by each host_v2
 * translation unit; everything is static inline. */
#ifndef DUCKVEP_HOST_V2_COMMON_H
#define DUCKVEP_HOST_V2_COMMON_H

#define DUCKDB_EXTENSION_NAME duckvep
#define DUCKDB_V2_API_ALLOW_UNSTABLE 0
#define DUCKDB_V2_API_ALLOW_DEPRECATED 0
#include "duckdb_extension_v2.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if DUCKDB_EXTENSION_HEADER_VERSION != 2 || DUCKDB_V2_API_VERSION_MAJOR != 2
#error "host_v2 requires the DuckDB C API v2 SDK"
#endif
#if DUCKDB_V2_API_ALLOW_UNSTABLE || DUCKDB_V2_API_ALLOW_DEPRECATED
#error "host_v2 must compile with the unstable and deprecated surfaces disabled"
#endif

DUCKDB_EXTENSION_EXTERN

/* ---------------------------------------------------------------------------
 * Adapter helpers: errors, vector views, validity, strings.
 * ------------------------------------------------------------------------- */

static inline duckdb_v2_str string_view(const char *text) {
    duckdb_v2_str view = {text, (idx_t)strlen(text)};
    return view;
}

/* Callback error handles are borrowed; each fallible call gets its own owned
 * `detail` slot, copied into the callback's handle and destroyed. */
static inline void set_error(duckdb_v2_error_info_handle target, DUCKDB_V2_ERROR code, const char *text) {
    (void)duckdb_v2_error_info_set_code(target, code);
    (void)duckdb_v2_error_info_set_text(target, string_view(text));
}

static inline void copy_duckdb_error(duckdb_v2_error_info_handle target, DUCKDB_V2_ERROR code,
                              duckdb_v2_error_info_handle source) {
    duckdb_v2_str text = string_view("duckvep: DuckDB C API v2 call failed");
    if (source) {
        (void)duckdb_v2_error_info_get_text(source, &text);
    }
    (void)duckdb_v2_error_info_set_code(target, code);
    (void)duckdb_v2_error_info_set_text(target, text);
}

/* Requires locals `error` (duckdb_v2_error_info_handle *), `detail`
 * (duckdb_v2_error_info_handle) and a `cleanup:` label. */
#define DUCKDB_CALL(expression)                                                                    \
    do {                                                                                           \
        DUCKDB_V2_ERROR call_status = (expression);                                                \
        if (call_status != DUCKDB_V2_ERROR_NONE) {                                                 \
            copy_duckdb_error(*error, call_status, detail);                                        \
            goto cleanup;                                                                          \
        }                                                                                          \
    } while (0)

#define INPUT_ERROR(text)                                                                          \
    do {                                                                                           \
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, text);                                    \
        goto cleanup;                                                                              \
    } while (0)

/* Selection vector: a view addresses logical row i at data[sel[i]] when sel is
 * set. Validity is indexed by the same physical position. */
static inline idx_t physical_row(const duckdb_v2_vector_view *view, idx_t row) {
    return view->sel ? view->sel[row] : row;
}

static inline bool row_is_valid(const duckdb_v2_vector_view *view, idx_t row) {
    idx_t position = physical_row(view, row);
    return !view->validity ||
           (view->validity[position >> 6] & (UINT64_C(1) << (position & 63))) != 0;
}

static inline void mark_valid(uint64_t *mask, idx_t row) {
    mask[row >> 6] |= UINT64_C(1) << (row & 63);
}

static inline void mark_null(uint64_t *mask, idx_t row) {
    mask[row >> 6] &= ~(UINT64_C(1) << (row & 63));
}

static inline uint64_t integer_at(const duckdb_v2_vector_view *view, idx_t row) {
    return ((const uint64_t *)view->data)[physical_row(view, row)];
}

/* VARCHAR payloads are borrowed from the input chunk. */
static inline duckdb_v2_str string_at(const duckdb_v2_vector_view *view, idx_t row) {
    const duckdb_v2_bytes *bytes = &((const duckdb_v2_bytes *)view->data)[physical_row(view, row)];
    duckdb_v2_str text;
    text.len = bytes->value.inlined.length;
    text.ptr = text.len <= sizeof(bytes->value.inlined.inlined) ? bytes->value.inlined.inlined
                                                                : bytes->value.pointer.ptr;
    return text;
}

static inline DUCKDB_V2_ERROR write_string(duckdb_v2_arena_handle arena, duckdb_v2_bytes *output,
                                    const char *bytes, size_t length,
                                    duckdb_v2_error_info_handle *detail) {
    memset(output, 0, sizeof(*output));
    if (length > UINT32_MAX) {
        return DUCKDB_V2_ERROR_INPUT_INVALID;
    }
    output->value.inlined.length = (uint32_t)length;
    if (length <= sizeof(output->value.inlined.inlined)) {
        if (length > 0) {
            memcpy(output->value.inlined.inlined, bytes, length);
        }
        return DUCKDB_V2_ERROR_NONE;
    }
    uint8_t *target = NULL;
    DUCKDB_V2_ERROR status = duckdb_v2_arena_allocate(arena, (idx_t)length, &target, detail);
    if (status != DUCKDB_V2_ERROR_NONE) {
        return status;
    }
    memcpy(target, bytes, length);
    memcpy(output->value.pointer.prefix, bytes, sizeof(output->value.pointer.prefix));
    output->value.pointer.ptr = (char *)target;
    return DUCKDB_V2_ERROR_NONE;
}

/* Borrows the views of `count` arguments. Flat, constant and dictionary
 * vectors are read in place through their selection; representations the view
 * getter rejects are flattened first. Reading a dictionary flattens its child
 * in place, which can invalidate a view taken earlier over an aliased buffer,
 * so every view is taken twice and only the second, stable set is returned. */
static inline bool load_views(duckdb_v2_scalar_function_exec_info_handle info, idx_t count,
                       duckdb_v2_vector_view *views, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle vectors[4] = {0};
    bool success = false;
    for (idx_t argument = 0; argument < count; ++argument) {
        DUCKDB_V2_VECTOR_TYPE type = DUCKDB_V2_VECTOR_TYPE_OTHER;
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, argument, &vectors[argument], &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_vector_type(vectors[argument], &type, &detail));
        if (type == DUCKDB_V2_VECTOR_TYPE_OTHER) {
            DUCKDB_CALL(duckdb_v2_vector_flatten(vectors[argument], &detail));
        }
    }
    for (int pass = 0; pass < 2; ++pass) {
        for (idx_t argument = 0; argument < count; ++argument) {
            DUCKDB_CALL(duckdb_v2_vector_get_view(vectors[argument], &views[argument], &detail));
        }
    }
    success = true;
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

/* Prepares a flat STRUCT result: sized, every row and field initially valid
 * (callers null rows through duckdb_v2_vector_set_null, which also nulls the
 * children). `fields`, `data` and `validity` receive the borrowed children. */
static inline bool open_struct_result(duckdb_v2_scalar_function_exec_info_handle info, idx_t row_count,
                               idx_t field_count, duckdb_v2_vector_handle *result,
                               duckdb_v2_vector_handle *fields, void **data, uint64_t **validity,
                               duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    uint64_t *struct_validity = NULL;
    bool success = false;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, result, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(*result, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(*result, row_count, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(*result, &struct_validity, &detail));
    for (idx_t row = 0; row < row_count; ++row) {
        mark_valid(struct_validity, row);
    }
    for (idx_t field = 0; field < field_count; ++field) {
        DUCKDB_CALL(duckdb_v2_vector_get_child(*result, field, &fields[field], &detail));
        DUCKDB_CALL(duckdb_v2_vector_set_size(fields[field], row_count, &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(fields[field], &data[field], &detail));
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(fields[field], &validity[field], &detail));
        for (idx_t row = 0; row < row_count; ++row) {
            mark_valid(validity[field], row);
        }
    }
    success = true;
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

/* Type text is parsed by the running context, so types need no enum table. */
static inline DUCKDB_V2_ERROR make_type(duckdb_v2_context_handle context, const char *text,
                                 duckdb_v2_logical_type_handle *type,
                                 duckdb_v2_error_info_handle *detail) {
    /* ANY is a signature wildcard; it has no text form. */
    if (strcmp(text, "ANY") == 0) {
        return duckdb_v2_context_create_type_from_id(context, DUCKDB_V2_LOGICAL_TYPE_ID_ANY, NULL,
                                                     NULL, 0, type, detail);
    }
    return duckdb_v2_context_create_type_from_text(context, string_view(text), type, detail);
}

static inline bool register_scalar(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                            const char *name_text, const char *const *parameter_types,
                            const char *const *parameter_names, idx_t parameter_count,
                            const char *return_type_text,
                            duckdb_v2_scalar_function_exec_callback_fn callback,
                            duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_scalar_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    bool success = false;
    duckdb_v2_str name = string_view(name_text);
    DUCKDB_CALL(duckdb_v2_scalar_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_name(function, &name, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_get_signature(function, &signature, &detail));
    for (idx_t parameter = 0; parameter < parameter_count; ++parameter) {
        DUCKDB_CALL(make_type(context, parameter_types[parameter], &type, &detail));
        DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(
            signature, string_view(parameter_names[parameter]), type, NULL, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    DUCKDB_CALL(make_type(context, return_type_text, &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_set_return_type(signature, type, &detail));
    /* The callbacks handle NULL inputs themselves (NULL in, NULL out). */
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(
        function, DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING,
        DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING_SPECIAL, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_exec_callback(function, callback, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_register(function, &detail));
    success = true;
cleanup:
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_scalar_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}


#endif
