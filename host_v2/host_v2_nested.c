/* The nested-vector family of the v2 host: duckvep_repeat_alleles and
 * duckvep_phase_call (LIST and STRUCT inputs, numeric cells, LIST<STRUCT>
 * output) and the small kernels _duckvep_revcomp, _duckvep_raw_gt and
 * _duckvep_record_order. The row logic is in src/core; this file reads v2
 * vectors into neutral cells and writes the results back. */
#include "host_v2_common.h"

#include "core/duckvep_core_cells.h"
#include "core/duckvep_core_phase.h"
#include "core/duckvep_core_repeat.h"
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

static duckvep_cell_kind_t kind_of(DUCKDB_V2_LOGICAL_TYPE_ID id) {
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
static bool read_type(duckdb_v2_vector_handle vector, type_info *info,
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
static bool column_open(duckdb_v2_vector_handle vector, bool force_flatten, column *out,
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

static bool column_child(const column *parent, idx_t index, column *child,
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

static bool column_valid(const column *c, idx_t index) {
    return c->has_view && row_is_valid(&c->view, index);
}

/* One element of a leaf column as a neutral cell. */
static void fill_cell(const column *c, idx_t index, duckvep_cell_t *cell) {
    *cell = (duckvep_cell_t){0};
    cell->kind = c->type.kind;
    cell->scale = c->type.scale;
    cell->valid = column_valid(c, index);
    if (!c->has_view) {
        cell->kind = DUCKVEP_CELL_UNSUPPORTED;
        return;
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
    case DUCKVEP_CELL_BOOLEAN: cell->boolean = ((const bool *)data)[at]; break;
    default: break;
    }
}

static duckdb_v2_list_entry list_at(const column *c, idx_t index) {
    return ((const duckdb_v2_list_entry *)c->view.data)[physical_row(&c->view, index)];
}

/* A native-budget refusal becomes an explicit capacity error. */
static void report(duckdb_v2_error_info_handle target, const char *message) {
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

#define MAX_OPTIONS 2

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

static duckvep_core_field_type_t field_category(DUCKDB_V2_LOGICAL_TYPE_ID id) {
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
static bool options_open(duckdb_v2_vector_handle vector, const char *const *names,
                         const duckvep_core_option_kind_t *kinds, size_t count, options *out,
                         duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_value_handle value = NULL;
    idx_t fields = 0;
    bool success = false;
    memset(out, 0, sizeof(*out));
    out->count = count;
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
static bool options_check(const options *o, idx_t row, duckdb_v2_error_info_handle *error) {
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
static bool option_cell(const options *o, size_t at, idx_t row, duckvep_cell_t *cell) {
    if (!o->present || !column_valid(&o->record, row) || !o->has[at]) {
        return false;
    }
    fill_cell(&o->fields[at], row, cell);
    return true;
}

/* ---------------------------------------------------------------------------
 * duckvep_repeat_alleles(reference, alternate, sequence_exact [, options])
 * ------------------------------------------------------------------------- */

typedef struct {
    bool is_list;
    bool usable;           /* elements are {unit VARCHAR, count numeric} records */
    column list, records, units, counts;
} axis;

/* An axis is a LIST of {unit, count} records (fields in either order), a list of
 * NULLs, or a bare NULL. Anything else is the "expected lists" error. */
static bool axis_open(duckdb_v2_vector_handle vector, idx_t rows, axis *out, bool *well_formed,
                      duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_value_handle value = NULL;
    idx_t fields = 0, unit_index = 0;
    bool success = false;
    memset(out, 0, sizeof(*out));
    *well_formed = true;
    if (!column_open(vector, true, &out->list, error)) {
        return false;
    }
    if (out->list.type.id != DUCKDB_V2_LOGICAL_TYPE_ID_LIST) {
        /* Only an all-NULL argument (an untyped NULL) is accepted. */
        *well_formed = rows == 0 || !column_valid(&out->list, 0);
        return true;
    }
    out->is_list = true;
    if (!column_child(&out->list, 0, &out->records, error)) {
        return false;
    }
    if (out->records.type.id == DUCKDB_V2_LOGICAL_TYPE_ID_SQLNULL) {
        return true;
    }
    if (out->records.type.id != DUCKDB_V2_LOGICAL_TYPE_ID_STRUCT) {
        *well_formed = false;
        return true;
    }
    DUCKDB_CALL(duckdb_v2_vector_get_logical_type(out->records.vector, &type, &detail));
    DUCKDB_CALL(duckdb_v2_logical_type_get_param_count(type, &fields, &detail));
    if (fields == 2) {
        duckdb_v2_identifier_t first, second;
        DUCKDB_CALL(duckdb_v2_logical_type_get_param(type, 0, &first, &value, &detail));
        DUCKDB_CALL(duckdb_v2_value_destroy(&value));
        DUCKDB_CALL(duckdb_v2_logical_type_get_param(type, 1, &second, &value, &detail));
        DUCKDB_CALL(duckdb_v2_value_destroy(&value));
        unit_index = first.len == 4 && memcmp(first.ptr, "unit", 4) == 0 ? 0 : 1;
        bool named = unit_index == 0
            ? (second.len == 5 && memcmp(second.ptr, "count", 5) == 0)
            : (first.len == 5 && memcmp(first.ptr, "count", 5) == 0 &&
               second.len == 4 && memcmp(second.ptr, "unit", 4) == 0);
        if (!named) {
            *well_formed = false;
            success = true;
            goto cleanup;
        }
    } else {
        *well_formed = false;
        success = true;
        goto cleanup;
    }
    if (!column_child(&out->records, unit_index, &out->units, error) ||
        !column_child(&out->records, 1 - unit_index, &out->counts, error)) {
        goto cleanup;
    }
    if (out->units.type.id != DUCKDB_V2_LOGICAL_TYPE_ID_VARCHAR) {
        *well_formed = false;
        success = true;
        goto cleanup;
    }
    out->usable = true;
    success = true;
cleanup:
    (void)duckdb_v2_value_destroy(&value);
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

typedef struct {
    duckvep_repeat_element_t *elements;
    size_t capacity;
} element_buffer;

/* One list row of an axis as neutral elements; false on allocation failure. */
static bool axis_row(const axis *a, idx_t row, element_buffer *buffer, duckvep_repeat_axis_t *out) {
    *out = (duckvep_repeat_axis_t){0};
    if (!a->is_list || !a->usable || !column_valid(&a->list, row)) {
        return true;
    }
    duckdb_v2_list_entry list = list_at(&a->list, row);
    if (list.length > buffer->capacity) {
        duckvep_repeat_element_t *grown = duckvep_budget_realloc(
            DUCKVEP_OWNER_CONTROL, buffer->elements, list.length * sizeof(*buffer->elements));
        if (!grown) {
            return false;
        }
        buffer->elements = grown;
        buffer->capacity = list.length;
    }
    for (idx_t i = 0; i < list.length; ++i) {
        idx_t at = list.offset + i;
        duckvep_repeat_element_t *element = &buffer->elements[i];
        duckdb_v2_str unit;
        *element = (duckvep_repeat_element_t){0};
        element->present = column_valid(&a->records, at) && column_valid(&a->units, at) &&
                           column_valid(&a->counts, at);
        unit = string_at(&a->units.view, at);
        element->unit = unit.ptr;
        element->unit_length = (uint32_t)unit.len;
        fill_cell(&a->counts, at, &element->count);
    }
    out->usable = true;
    out->elements = buffer->elements;
    out->count = list.length;
    return true;
}

enum { REPEAT_REFERENCE, REPEAT_ALTERNATE, REPEAT_REFERENCE_LENGTH, REPEAT_ALTERNATE_LENGTH,
       REPEAT_LENGTH_CHANGE, REPEAT_DIRECTION, REPEAT_STATUS, REPEAT_FIELDS };

static void repeat_exec(duckdb_v2_scalar_function_exec_info_handle info,
                        duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    static const char *const keys[] = {"max_allele_bases"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_NUMERIC};
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle arguments[4] = {0};
    duckdb_v2_vector_handle result = NULL;
    duckdb_v2_vector_handle fields[REPEAT_FIELDS] = {0};
    void *data[REPEAT_FIELDS] = {0};
    uint64_t *validity[REPEAT_FIELDS] = {0};
    duckdb_v2_arena_handle arenas[REPEAT_FIELDS] = {0};
    axis axes[2];
    column exact;
    options option;
    element_buffer buffers[2] = {{NULL, 0}, {NULL, 0}};
    bool well_formed[2] = {true, true};
    idx_t rows = 0;
    uint32_t argc = 0;
    bool opened = false;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    for (uint32_t i = 0; i < argc; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
    }
    memset(axes, 0, sizeof axes);
    if (!axis_open(arguments[0], rows, &axes[0], &well_formed[0], error) ||
        !axis_open(arguments[1], rows, &axes[1], &well_formed[1], error) ||
        !column_open(arguments[2], false, &exact, error) ||
        !options_open(argc == 4 ? arguments[3] : NULL, keys, kinds, 1, &option, error)) {
        goto cleanup;
    }
    opened = true;
    if (!well_formed[0] || !well_formed[1]) {
        report(*error, duckvep_core_repeat_expected_lists);
        goto cleanup;
    }
    if (!open_struct_result(info, rows, REPEAT_FIELDS, &result, fields, data, validity, error)) {
        goto cleanup;
    }
    for (idx_t field = 0; field < REPEAT_FIELDS; ++field) {
        if (field == REPEAT_REFERENCE || field == REPEAT_ALTERNATE || field == REPEAT_DIRECTION ||
            field == REPEAT_STATUS) {
            DUCKDB_CALL(duckdb_v2_vector_get_arena(fields[field], &arenas[field], &detail));
        }
    }
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_cell_t cap_cell;
        duckvep_repeat_axis_t parts[2];
        duckvep_repeat_plan_t plan;
        char message[240];
        long double cap = 5000;
        bool have_cap;
        if (!options_check(&option, row, error)) {
            goto cleanup;
        }
        if (!column_valid(&exact, row)) {
            report(*error, duckvep_core_repeat_exact_required);
            goto cleanup;
        }
        have_cap = option_cell(&option, 0, row, &cap_cell);
        if (!duckvep_core_repeat_cap(have_cap ? &cap_cell : NULL, &cap)) {
            report(*error, duckvep_core_repeat_cap_invalid);
            goto cleanup;
        }
        if (!axis_row(&axes[0], row, &buffers[0], &parts[0]) ||
            !axis_row(&axes[1], row, &buffers[1], &parts[1])) {
            report(*error, "duckvep_repeat_alleles: allocation failed");
            goto cleanup;
        }
        duckvep_core_repeat_plan(parts, ((const bool *)exact.view.data)[physical_row(&exact.view, row)],
                                 cap, &plan, message, sizeof message);
        if (plan.error) {
            report(*error, plan.error);
            goto cleanup;
        }
        DUCKDB_CALL(write_string(arenas[REPEAT_STATUS], &((duckdb_v2_bytes *)data[REPEAT_STATUS])[row],
                                 plan.status, strlen(plan.status), &detail));
        if (strcmp(plan.status, "ok") != 0) {
            for (idx_t field = 0; field < REPEAT_STATUS; ++field) {
                mark_null(validity[field], row);
            }
            continue;
        }
        for (idx_t part = 0; part < 2; ++part) {
            size_t length = (size_t)plan.required[part];
            char *text = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, length + 1);
            size_t written;
            DUCKDB_V2_ERROR status;
            if (!text) {
                report(*error, "duckvep_repeat_alleles: allocation failed");
                goto cleanup;
            }
            written = duckvep_core_repeat_render(&parts[part], text);
            status = write_string(arenas[part], &((duckdb_v2_bytes *)data[part])[row], text, written,
                                  &detail);
            duckvep_budget_free(text);
            DUCKDB_CALL(status);
        }
        ((uint64_t *)data[REPEAT_REFERENCE_LENGTH])[row] = (uint64_t)plan.required[0];
        ((uint64_t *)data[REPEAT_ALTERNATE_LENGTH])[row] = (uint64_t)plan.required[1];
        ((int64_t *)data[REPEAT_LENGTH_CHANGE])[row] =
            (int64_t)plan.required[1] - (int64_t)plan.required[0];
        {
            const char *direction = duckvep_core_repeat_direction(&plan);
            DUCKDB_CALL(write_string(arenas[REPEAT_DIRECTION],
                                     &((duckdb_v2_bytes *)data[REPEAT_DIRECTION])[row], direction,
                                     strlen(direction), &detail));
        }
    }
cleanup:
    (void)opened;
    duckvep_budget_free(buffers[0].elements);
    duckvep_budget_free(buffers[1].elements);
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_phase_call(alleles, phase_before [, options]) -> LIST<STRUCT>
 * ------------------------------------------------------------------------- */

enum { PHASE_INPUT_SLOT, PHASE_ALLELE_INDEX, PHASE_LANE, PHASE_PLOIDY, PHASE_SET, PHASE_SCOPE,
       PHASE_STATUS, PHASE_FIELDS };

typedef struct {
    duckdb_v2_error_info_handle *error;
    duckdb_v2_error_info_handle detail;
    const column *alleles;   /* the allele list's element column (NULL: an untyped NULL list) */
    const column *phases;
    idx_t allele_base, phase_base, at;
    void *data[PHASE_FIELDS];
    uint64_t *validity[PHASE_FIELDS];
    duckdb_v2_arena_handle scope_arena, status_arena;
    bool has_phase_set;
    int64_t phase_set;
} phase_row;

static void phase_allele_cell(void *pointer, size_t slot, duckvep_cell_t *cell) {
    phase_row *row = pointer;
    fill_cell(row->alleles, row->allele_base + slot, cell);
}

static void phase_flag_cell(void *pointer, size_t slot, duckvep_cell_t *cell) {
    phase_row *row = pointer;
    fill_cell(row->phases, row->phase_base + slot, cell);
}

static bool phase_emit(void *pointer, size_t slot, const duckvep_core_phase_slot_t *out) {
    phase_row *row = pointer;
    idx_t at = row->at + slot;
    duckdb_v2_error_info_handle *error = row->error;
    duckdb_v2_error_info_handle detail = NULL;
    bool success = false;
    ((uint16_t *)row->data[PHASE_INPUT_SLOT])[at] = out->input_slot;
    ((int32_t *)row->data[PHASE_ALLELE_INDEX])[at] = out->allele_index;
    if (!out->allele_called) {
        mark_null(row->validity[PHASE_ALLELE_INDEX], at);
    }
    ((uint16_t *)row->data[PHASE_LANE])[at] = out->lane;
    if (!out->lane) {
        mark_null(row->validity[PHASE_LANE], at);
    }
    ((uint16_t *)row->data[PHASE_PLOIDY])[at] = out->ploidy;
    if (out->phase_set_applies && row->has_phase_set) {
        ((int64_t *)row->data[PHASE_SET])[at] = row->phase_set;
    } else {
        mark_null(row->validity[PHASE_SET], at);
    }
    DUCKDB_CALL(write_string(row->scope_arena, &((duckdb_v2_bytes *)row->data[PHASE_SCOPE])[at],
                             out->scope, strlen(out->scope), &detail));
    DUCKDB_CALL(write_string(row->status_arena, &((duckdb_v2_bytes *)row->data[PHASE_STATUS])[at],
                             out->status, strlen(out->status), &detail));
    success = true;
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

static void phase_exec(duckdb_v2_scalar_function_exec_info_handle info,
                       duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    static const char *const keys[] = {"phase_set", "phase_policy"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_INTEGER,
                                                       DUCKVEP_CORE_OPTION_TEXT};
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle arguments[3] = {0};
    duckdb_v2_vector_handle result = NULL, records = NULL;
    duckdb_v2_vector_handle fields[PHASE_FIELDS] = {0};
    uint64_t *list_validity = NULL;
    void *list_data = NULL;
    column alleles, phases, allele_values, phase_values;
    options option;
    bool allele_list, phase_list;
    idx_t rows = 0, at = 0;
    uint32_t argc = 0;
    size_t total = 0;
    phase_row row_state;
    memset(&row_state, 0, sizeof row_state);
    row_state.error = error;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    for (uint32_t i = 0; i < argc; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
    }
    if (!column_open(arguments[0], true, &alleles, error) ||
        !column_open(arguments[1], true, &phases, error) ||
        !options_open(argc == 3 ? arguments[2] : NULL, keys, kinds, 2, &option, error)) {
        goto cleanup;
    }
    allele_list = alleles.type.id == DUCKDB_V2_LOGICAL_TYPE_ID_LIST;
    phase_list = phases.type.id == DUCKDB_V2_LOGICAL_TYPE_ID_LIST;
    memset(&allele_values, 0, sizeof allele_values);
    memset(&phase_values, 0, sizeof phase_values);
    if ((allele_list && !column_child(&alleles, 0, &allele_values, error)) ||
        (phase_list && !column_child(&phases, 0, &phase_values, error))) {
        goto cleanup;
    }
    /* A phase flag's DECIMAL is unsupported (only an allele index is scaled). */
    if (phase_values.type.id == DUCKDB_V2_LOGICAL_TYPE_ID_DECIMAL) {
        phase_values.type.kind = DUCKVEP_CELL_UNSUPPORTED;
        phase_values.type.scale = 0;
    }
    row_state.alleles = &allele_values;
    row_state.phases = &phase_values;

    for (idx_t row = 0; row < rows; ++row) {
        bool have_gt = allele_list && column_valid(&alleles, row);
        bool have_phase = phase_list && column_valid(&phases, row);
        const char *message = duckvep_core_phase_check_row(
            have_gt, have_phase, have_gt ? list_at(&alleles, row).length : 0,
            have_phase ? list_at(&phases, row).length : 0, &total);
        if (message) {
            INPUT_ERROR(message);
        }
    }
    /* The output list's child is sized once for every slot of the chunk: a
     * chunk of 2,048 diploid rows needs 4,096 records, more than one vector. */
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &result, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(result, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(result, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(result, &list_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(result, &list_validity, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_child(result, 0, &records, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(records, total, &detail));
    for (idx_t field = 0; field < PHASE_FIELDS; ++field) {
        DUCKDB_CALL(duckdb_v2_vector_get_child(records, field, &fields[field], &detail));
        DUCKDB_CALL(duckdb_v2_vector_set_size(fields[field], total, &detail));
        DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(fields[field], &row_state.data[field], &detail));
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(fields[field],
                                                               &row_state.validity[field], &detail));
        for (idx_t slot = 0; slot < total; ++slot) {
            mark_valid(row_state.validity[field], slot);
        }
    }
    DUCKDB_CALL(duckdb_v2_vector_get_arena(fields[PHASE_SCOPE], &row_state.scope_arena, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(fields[PHASE_STATUS], &row_state.status_arena, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_cell_t cell;
        duckvep_phase_policy_t policy = DUCKVEP_PHASE_STRICT;
        const char *message = NULL;
        bool have_phase_set = false;
        int64_t phase_set = 0;
        idx_t count;
        bool have_phase;
        duckvep_core_phase_result_t outcome;
        if (!options_check(&option, row, error)) {
            goto cleanup;
        }
        if (option_cell(&option, 0, row, &cell) && cell.valid) {
            have_phase_set = true;
            if (!duckvep_core_phase_set(&cell, &phase_set)) {
                INPUT_ERROR(duckvep_core_phase_set_error);
            }
        }
        if (option_cell(&option, 1, row, &cell) && cell.valid &&
            !duckvep_core_phase_policy(cell.text, cell.text_length, &policy)) {
            INPUT_ERROR(duckvep_core_phase_policy_error);
        }
        ((duckdb_v2_list_entry *)list_data)[row] = (duckdb_v2_list_entry){at, 0};
        if (!allele_list || !column_valid(&alleles, row)) {
            mark_null(list_validity, row);
            continue;
        }
        mark_valid(list_validity, row);
        count = list_at(&alleles, row).length;
        have_phase = phase_list && column_valid(&phases, row);
        row_state.allele_base = list_at(&alleles, row).offset;
        row_state.phase_base = have_phase ? list_at(&phases, row).offset : 0;
        row_state.at = at;
        row_state.has_phase_set = have_phase_set;
        row_state.phase_set = phase_set;
        outcome = duckvep_core_phase_row(
            &(duckvep_core_phase_reader_t){&row_state, phase_allele_cell, phase_flag_cell, phase_emit},
            count, have_phase, policy, &message);
        if (outcome == DUCKVEP_CORE_PHASE_ERROR) {
            INPUT_ERROR(message);
        }
        if (outcome == DUCKVEP_CORE_PHASE_HOST_FAILED) {
            goto cleanup;
        }
        ((duckdb_v2_list_entry *)list_data)[row].length = count;
        at += count;
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * _duckvep_revcomp(VARCHAR), _duckvep_raw_gt(VARCHAR, UINTEGER),
 * _duckvep_record_order(UBIGINT, UBIGINT)
 * ------------------------------------------------------------------------- */

static void revcomp_exec(duckdb_v2_scalar_function_exec_info_handle info,
                         duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[1];
    duckdb_v2_vector_handle output = NULL;
    duckdb_v2_arena_handle arena = NULL;
    void *data = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 1, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(output, &arena, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        duckdb_v2_str sequence;
        char *reversed;
        DUCKDB_V2_ERROR status;
        if (!row_is_valid(&views[0], row)) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(output, row, &detail));
            continue;
        }
        sequence = string_at(&views[0], row);
        reversed = duckvep_budget_malloc(DUCKVEP_OWNER_CONTROL, sequence.len ? sequence.len : 1);
        if (!reversed) {
            set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "_duckvep_revcomp: allocation failed");
            goto cleanup;
        }
        duckvep_core_revcomp(sequence.ptr, sequence.len, reversed);
        status = write_string(arena, &((duckdb_v2_bytes *)data)[row], reversed, sequence.len, &detail);
        duckvep_budget_free(reversed);
        DUCKDB_CALL(status);
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void raw_gt_exec(duckdb_v2_scalar_function_exec_info_handle info,
                        duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[2];
    duckdb_v2_vector_handle result = NULL;
    duckdb_v2_vector_handle fields[7] = {0};
    void *data[7] = {0};
    uint64_t *validity[7] = {0};
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 2, views, error) ||
        !open_struct_result(info, rows, 7, &result, fields, data, validity, error)) {
        goto cleanup;
    }
    for (idx_t row = 0; row < rows; ++row) {
        uint32_t values[7];
        duckdb_v2_str gt;
        if (!row_is_valid(&views[0], row) || !row_is_valid(&views[1], row)) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(result, row, &detail));
            continue;
        }
        gt = string_at(&views[0], row);
        duckvep_core_raw_gt(gt.ptr, gt.len, ((const uint32_t *)views[1].data)[physical_row(&views[1], row)],
                            values);
        for (idx_t field = 0; field < 7; ++field) {
            ((uint32_t *)data[field])[row] = values[field];
        }
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void record_order_exec(duckdb_v2_scalar_function_exec_info_handle info,
                              duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_view views[2];
    duckdb_v2_vector_handle output = NULL;
    uint64_t *ranks = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 2, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, (void **)&ranks, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        if (!row_is_valid(&views[0], row) || !row_is_valid(&views[1], row)) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(output, row, &detail));
            continue;
        }
        ranks[row] = duckvep_core_record_order(integer_at(&views[0], row), integer_at(&views[1], row));
        if (!ranks[row]) {
            INPUT_ERROR(duckvep_core_record_order_error);
        }
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * Registration
 * ------------------------------------------------------------------------- */

bool host_v2_register_nested(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error) {
    static const char *const repeat_types[] = {"ANY", "ANY", "BOOLEAN", "ANY"};
    static const char *const repeat_names[] = {"reference_components", "alternate_components",
                                               "sequence_exact", "options"};
    static const char *const phase_types[] = {"ANY", "ANY", "ANY"};
    static const char *const phase_names[] = {"alleles", "phase_before", "options"};
    static const char *const text_types[] = {"VARCHAR", "UINTEGER"};
    static const char *const text_names[] = {"genotype", "source_alt_count"};
    static const char *const order_types[] = {"UBIGINT", "UBIGINT"};
    static const char *const order_names[] = {"record_count", "ordinal"};
    const char *repeat_result =
        "STRUCT(reference VARCHAR, alternate VARCHAR, reference_length UBIGINT, "
        "alternate_length UBIGINT, length_change BIGINT, length_direction VARCHAR, status VARCHAR)";
    const char *phase_result =
        "STRUCT(input_slot USMALLINT, allele_index INTEGER, haplotype_lane USMALLINT, "
        "ploidy USMALLINT, phase_set BIGINT, phase_scope VARCHAR, status VARCHAR)[]";
    for (idx_t arity = 3; arity <= 4; ++arity) {
        if (!register_scalar(extension, context, "duckvep_repeat_alleles", repeat_types, repeat_names,
                             arity, repeat_result, repeat_exec, error)) {
            return false;
        }
    }
    for (idx_t arity = 2; arity <= 3; ++arity) {
        if (!register_scalar(extension, context, "duckvep_phase_call", phase_types, phase_names, arity,
                             phase_result, phase_exec, error)) {
            return false;
        }
    }
    static const char *const one_text[] = {"VARCHAR"};
    static const char *const one_name[] = {"sequence"};
    return register_scalar(extension, context, "_duckvep_revcomp", one_text, one_name, 1, "VARCHAR",
                           revcomp_exec, error) &&
           register_scalar(extension, context, "_duckvep_raw_gt", text_types, text_names, 2,
                           "STRUCT(status UINTEGER, allele0 UINTEGER, allele1 UINTEGER, "
                           "parsed_slots UINTEGER, source_ploidy UINTEGER, source_has_missing UINTEGER, "
                           "disposition UINTEGER)",
                           raw_gt_exec, error) &&
           register_scalar(extension, context, "_duckvep_record_order", order_types, order_names, 2,
                           "UBIGINT", record_order_exec, error);
}
