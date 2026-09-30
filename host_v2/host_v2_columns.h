/* Typed columns, cells and option STRUCTs of the v2 host: reads vectors (flat,
 * constant, dictionary, nested) into the neutral cells of src/core. Included by
 * each host_v2 translation unit that reads arguments. */
#ifndef DUCKVEP_HOST_V2_COLUMNS_H
#define DUCKVEP_HOST_V2_COLUMNS_H

#include "host_v2_common.h"

#include "core/duckvep_core_cells.h"
#include "kernel/src/duckvep_budget.h"

/* ---------------------------------------------------------------------------
 * Typed columns
 * ------------------------------------------------------------------------- */

typedef struct {
    DUCKDB_V2_LOGICAL_TYPE_ID id;
    duckvep_cell_kind_t kind; /* a DECIMAL's storage kind */
    uint8_t scale;
} type_info;

typedef struct {
    duckdb_v2_vector_handle vector;
    duckdb_v2_vector_view view;
    type_info type;
    bool has_view; /* false for a SQLNULL-typed vector: every cell is NULL */
} column;

static inline duckvep_cell_kind_t kind_of(DUCKDB_V2_LOGICAL_TYPE_ID id) {
    switch (id) {
    case DUCKDB_V2_LOGICAL_TYPE_ID_TINYINT: return DUCKVEP_CELL_TINYINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_SMALLINT: return DUCKVEP_CELL_SMALLINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_INTEGER: return DUCKVEP_CELL_INTEGER;
    case DUCKDB_V2_LOGICAL_TYPE_ID_BIGINT: return DUCKVEP_CELL_BIGINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UTINYINT: return DUCKVEP_CELL_UTINYINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_USMALLINT: return DUCKVEP_CELL_USMALLINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UINTEGER: return DUCKVEP_CELL_UINTEGER;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UBIGINT: return DUCKVEP_CELL_UBIGINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_FLOAT: return DUCKVEP_CELL_FLOAT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_DOUBLE: return DUCKVEP_CELL_DOUBLE;
    case DUCKDB_V2_LOGICAL_TYPE_ID_HUGEINT: return DUCKVEP_CELL_HUGEINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_VARCHAR: return DUCKVEP_CELL_VARCHAR;
    case DUCKDB_V2_LOGICAL_TYPE_ID_BOOLEAN: return DUCKVEP_CELL_BOOLEAN;
    default: return DUCKVEP_CELL_UNSUPPORTED;
    }
}

/* The logical type of a vector; a DECIMAL reports its storage kind (width <= 4
 * SMALLINT, <= 9 INTEGER, <= 18 BIGINT, else HUGEINT) and scale. */
static inline bool read_type(duckdb_v2_vector_handle vector, type_info *info,
                      duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_value_handle width_value = NULL, scale_value = NULL;
    duckdb_v2_identifier_t name;
    uint8_t width = 0, scale = 0;
    bool success = false;
    DUCKDB_CALL(duckdb_v2_vector_get_logical_type(vector, &type, &detail));
    DUCKDB_CALL(duckdb_v2_logical_type_get_id(type, &info->id, &detail));
    info->kind = kind_of(info->id);
    info->scale = 0;
    if (info->id == DUCKDB_V2_LOGICAL_TYPE_ID_DECIMAL) {
        DUCKDB_CALL(duckdb_v2_logical_type_get_param(type, 0, &name, &width_value, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_get_param(type, 1, &name, &scale_value, &detail));
        DUCKDB_CALL(duckdb_v2_value_get_utinyint(width_value, &width, &detail));
        DUCKDB_CALL(duckdb_v2_value_get_utinyint(scale_value, &scale, &detail));
        info->kind = width <= 4 ? DUCKVEP_CELL_SMALLINT : width <= 9 ? DUCKVEP_CELL_INTEGER
            : width <= 18 ? DUCKVEP_CELL_BIGINT : DUCKVEP_CELL_HUGEINT;
        info->scale = scale;
    }
    success = true;
cleanup:
    (void)duckdb_v2_value_destroy(&width_value);
    (void)duckdb_v2_value_destroy(&scale_value);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

/* Opens a vector for reading. Nested vectors (and, to keep it simple, any
 * vector the caller asks for) are flattened first so that their children can be
 * addressed by list offset; flat, constant and dictionary leaves are read in
 * place through their selection. Views are taken twice (see load_views). */
static inline bool column_open(duckdb_v2_vector_handle vector, bool force_flatten, column *out,
                        duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    DUCKDB_V2_VECTOR_TYPE vector_type = DUCKDB_V2_VECTOR_TYPE_OTHER;
    bool success = false;
    memset(out, 0, sizeof(*out));
    out->vector = vector;
    if (!read_type(vector, &out->type, error)) {
        return false;
    }
    DUCKDB_CALL(duckdb_v2_vector_get_vector_type(vector, &vector_type, &detail));
    if (force_flatten || vector_type == DUCKDB_V2_VECTOR_TYPE_OTHER) {
        DUCKDB_CALL(duckdb_v2_vector_flatten(vector, &detail));
    }
    if (out->type.id != DUCKDB_V2_LOGICAL_TYPE_ID_SQLNULL) {
        DUCKDB_CALL(duckdb_v2_vector_get_view(vector, &out->view, &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_view(vector, &out->view, &detail));
        out->has_view = true;
    }
    success = true;
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

static inline bool column_child(const column *parent, idx_t index, column *child,
                         duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle vector = NULL;
    bool success = false;
    DUCKDB_CALL(duckdb_v2_vector_get_child(parent->vector, index, &vector, &detail));
    success = column_open(vector, true, child, error);
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

static inline bool column_valid(const column *c, idx_t index) {
    return c->has_view && row_is_valid(&c->view, index);
}

/* One element of a leaf column as a neutral cell. */
static inline void fill_cell(const column *c, idx_t index, duckvep_cell_t *cell) {
    *cell = (duckvep_cell_t){0};
    cell->kind = c->type.kind;
    cell->scale = c->type.scale;
    cell->decimal = c->type.id == DUCKDB_V2_LOGICAL_TYPE_ID_DECIMAL;
    cell->valid = column_valid(c, index);
    if (!c->has_view) {
        cell->kind = DUCKVEP_CELL_UNSUPPORTED;
        return;
    }
    if (!cell->valid) {
        return; /* a NULL element's storage is uninitialized: do not read it */
    }
    idx_t at = physical_row(&c->view, index);
    const void *data = c->view.data;
    switch (cell->kind) {
    case DUCKVEP_CELL_TINYINT: cell->i = ((const int8_t *)data)[at]; break;
    case DUCKVEP_CELL_SMALLINT: cell->i = ((const int16_t *)data)[at]; break;
    case DUCKVEP_CELL_INTEGER: cell->i = ((const int32_t *)data)[at]; break;
    case DUCKVEP_CELL_BIGINT: cell->i = ((const int64_t *)data)[at]; break;
    case DUCKVEP_CELL_UTINYINT: cell->u = ((const uint8_t *)data)[at]; break;
    case DUCKVEP_CELL_USMALLINT: cell->u = ((const uint16_t *)data)[at]; break;
    case DUCKVEP_CELL_UINTEGER: cell->u = ((const uint32_t *)data)[at]; break;
    case DUCKVEP_CELL_UBIGINT: cell->u = ((const uint64_t *)data)[at]; break;
    case DUCKVEP_CELL_FLOAT: cell->d = ((const float *)data)[at]; break;
    case DUCKVEP_CELL_DOUBLE: cell->d = ((const double *)data)[at]; break;
    case DUCKVEP_CELL_HUGEINT:
        cell->hugeint_lower = ((const duckdb_v2_hugeint_t *)data)[at].lower;
        cell->hugeint_upper = ((const duckdb_v2_hugeint_t *)data)[at].upper;
        break;
    case DUCKVEP_CELL_VARCHAR: {
        duckdb_v2_str text = string_at(&c->view, index);
        cell->text = text.ptr;
        cell->text_length = text.len;
        break;
    }
    case DUCKVEP_CELL_BOOLEAN: cell->boolean = ((const uint8_t *)data)[at] != 0; break;
    default: break;
    }
}

static inline duckdb_v2_list_entry list_at(const column *c, idx_t index) {
    return ((const duckdb_v2_list_entry *)c->view.data)[physical_row(&c->view, index)];
}

/* A native-budget refusal becomes an explicit capacity error. */
static inline void report(duckdb_v2_error_info_handle target, const char *message) {
    char capacity[192], text[768];
    if (duckvep_budget_take_failure(capacity, sizeof capacity)) {
        (void)snprintf(text, sizeof text, "%s (%s)", capacity, message);
        message = text;
    }
    set_error(target, DUCKDB_V2_ERROR_INPUT_INVALID, message);
}

/* ---------------------------------------------------------------------------
 * Option STRUCT arguments
 * ------------------------------------------------------------------------- */

#define MAX_OPTIONS 6

typedef struct {
    column record;
    bool present;
    size_t count;
    column fields[MAX_OPTIONS];
    bool has[MAX_OPTIONS];
    bool not_struct;
    bool rejected;
    char message[256];
} options;

static inline duckvep_core_field_type_t field_category(DUCKDB_V2_LOGICAL_TYPE_ID id) {
    switch (id) {
    case DUCKDB_V2_LOGICAL_TYPE_ID_SQLNULL: return DUCKVEP_CORE_FIELD_SQLNULL;
    case DUCKDB_V2_LOGICAL_TYPE_ID_VARCHAR: return DUCKVEP_CORE_FIELD_VARCHAR;
    case DUCKDB_V2_LOGICAL_TYPE_ID_BOOLEAN: return DUCKVEP_CORE_FIELD_BOOLEAN;
    case DUCKDB_V2_LOGICAL_TYPE_ID_TINYINT: case DUCKDB_V2_LOGICAL_TYPE_ID_SMALLINT:
    case DUCKDB_V2_LOGICAL_TYPE_ID_INTEGER: case DUCKDB_V2_LOGICAL_TYPE_ID_BIGINT:
    case DUCKDB_V2_LOGICAL_TYPE_ID_UTINYINT: case DUCKDB_V2_LOGICAL_TYPE_ID_USMALLINT:
    case DUCKDB_V2_LOGICAL_TYPE_ID_UINTEGER: case DUCKDB_V2_LOGICAL_TYPE_ID_UBIGINT:
        return DUCKVEP_CORE_FIELD_INTEGER;
    case DUCKDB_V2_LOGICAL_TYPE_ID_HUGEINT: case DUCKDB_V2_LOGICAL_TYPE_ID_DECIMAL:
    case DUCKDB_V2_LOGICAL_TYPE_ID_FLOAT: case DUCKDB_V2_LOGICAL_TYPE_ID_DOUBLE:
        return DUCKVEP_CORE_FIELD_NUMERIC;
    default: return DUCKVEP_CORE_FIELD_OTHER;
    }
}

/* Indexes the struct's fields against `names`; a field the function does not
 * know, or of the wrong type, is an error the first valid row reports (as the
 * v1 host does: a NULL options row is never checked). */
static inline bool options_open(duckdb_v2_vector_handle vector, const char *const *names,
                         const duckvep_core_option_kind_t *kinds, size_t count, options *out,
                         duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_value_handle value = NULL;
    idx_t fields = 0;
    bool success = false;
    memset(out, 0, sizeof(*out));
    out->count = count;
    if (count > MAX_OPTIONS) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, "duckvep host_v2: too many options");
        return false;
    }
    if (!vector) {
        return true;
    }
    out->present = true;
    if (!column_open(vector, true, &out->record, error)) {
        return false;
    }
    if (out->record.type.id != DUCKDB_V2_LOGICAL_TYPE_ID_STRUCT) {
        out->not_struct = true;
        return true;
    }
    DUCKDB_CALL(duckdb_v2_vector_get_logical_type(vector, &type, &detail));
    DUCKDB_CALL(duckdb_v2_logical_type_get_param_count(type, &fields, &detail));
    for (idx_t i = 0; i < fields && !out->rejected; ++i) {
        duckdb_v2_identifier_t name;
        char key[128];
        column field;
        size_t at;
        DUCKDB_CALL(duckdb_v2_logical_type_get_param(type, i, &name, &value, &detail));
        DUCKDB_CALL(duckdb_v2_value_destroy(&value));
        (void)snprintf(key, sizeof key, "%.*s", (int)(name.len < sizeof key - 1 ? name.len : sizeof key - 1),
                       name.ptr);
        if (!column_child(&out->record, i, &field, error)) {
            goto cleanup;
        }
        at = duckvep_core_option_index(key, names, count);
        if (!duckvep_core_option_permitted(key, at, count, kinds, field_category(field.type.id),
                                           out->message, sizeof out->message)) {
            out->rejected = true;
            break;
        }
        out->fields[at] = field;
        out->has[at] = true;
    }
    success = true;
cleanup:
    (void)duckdb_v2_value_destroy(&value);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

/* Reports a bad options argument on the first row where it is non-NULL. */
static inline bool options_check(const options *o, idx_t row, duckdb_v2_error_info_handle *error) {
    if (!o->present || !column_valid(&o->record, row)) {
        return true;
    }
    if (o->not_struct) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, duckvep_core_option_not_struct);
        return false;
    }
    if (o->rejected) {
        set_error(*error, DUCKDB_V2_ERROR_INPUT_INVALID, o->message);
        return false;
    }
    return true;
}

/* The cell of option `at` for `row`; false when the option is absent for it. */
static inline bool option_cell(const options *o, size_t at, idx_t row, duckvep_cell_t *cell) {
    if (!o->present || !column_valid(&o->record, row) || !o->has[at]) {
        return false;
    }
    fill_cell(&o->fields[at], row, cell);
    return true;
}


#endif
