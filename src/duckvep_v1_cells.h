/* v1-host reader: one element of a v1 vector into a host-neutral cell. The
 * conversions themselves live in src/core/duckvep_core_cells.c. */
#ifndef DUCKVEP_V1_CELLS_H
#define DUCKVEP_V1_CELLS_H

#include "duckdb_extension.h"
#include "core/duckvep_core_cells.h"

static inline duckvep_cell_kind_t duckvep_v1_cell_kind(duckdb_type id) {
    switch (id) {
    case DUCKDB_TYPE_TINYINT: return DUCKVEP_CELL_TINYINT;
    case DUCKDB_TYPE_SMALLINT: return DUCKVEP_CELL_SMALLINT;
    case DUCKDB_TYPE_INTEGER: return DUCKVEP_CELL_INTEGER;
    case DUCKDB_TYPE_BIGINT: return DUCKVEP_CELL_BIGINT;
    case DUCKDB_TYPE_UTINYINT: return DUCKVEP_CELL_UTINYINT;
    case DUCKDB_TYPE_USMALLINT: return DUCKVEP_CELL_USMALLINT;
    case DUCKDB_TYPE_UINTEGER: return DUCKVEP_CELL_UINTEGER;
    case DUCKDB_TYPE_UBIGINT: return DUCKVEP_CELL_UBIGINT;
    case DUCKDB_TYPE_FLOAT: return DUCKVEP_CELL_FLOAT;
    case DUCKDB_TYPE_DOUBLE: return DUCKVEP_CELL_DOUBLE;
    case DUCKDB_TYPE_HUGEINT: return DUCKVEP_CELL_HUGEINT;
    case DUCKDB_TYPE_VARCHAR: return DUCKVEP_CELL_VARCHAR;
    case DUCKDB_TYPE_BOOLEAN: return DUCKVEP_CELL_BOOLEAN;
    default: return DUCKVEP_CELL_UNSUPPORTED;
    }
}

/* `id` is the storage type (a DECIMAL's internal type) and `scale` its scale. */
static inline void duckvep_v1_fill_cell(duckdb_vector vector, duckdb_type id, uint8_t scale,
                                        idx_t at, duckvep_cell_t *cell) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    void *data = duckdb_vector_get_data(vector);
    *cell = (duckvep_cell_t){0};
    cell->kind = duckvep_v1_cell_kind(id);
    cell->scale = scale;
    cell->valid = !validity || duckdb_validity_row_is_valid(validity, at);
    switch (cell->kind) {
    case DUCKVEP_CELL_TINYINT: cell->i = ((int8_t *)data)[at]; break;
    case DUCKVEP_CELL_SMALLINT: cell->i = ((int16_t *)data)[at]; break;
    case DUCKVEP_CELL_INTEGER: cell->i = ((int32_t *)data)[at]; break;
    case DUCKVEP_CELL_BIGINT: cell->i = ((int64_t *)data)[at]; break;
    case DUCKVEP_CELL_UTINYINT: cell->u = ((uint8_t *)data)[at]; break;
    case DUCKVEP_CELL_USMALLINT: cell->u = ((uint16_t *)data)[at]; break;
    case DUCKVEP_CELL_UINTEGER: cell->u = ((uint32_t *)data)[at]; break;
    case DUCKVEP_CELL_UBIGINT: cell->u = ((uint64_t *)data)[at]; break;
    case DUCKVEP_CELL_FLOAT: cell->d = ((float *)data)[at]; break;
    case DUCKVEP_CELL_DOUBLE: cell->d = ((double *)data)[at]; break;
    case DUCKVEP_CELL_HUGEINT:
        cell->hugeint_lower = ((duckdb_hugeint *)data)[at].lower;
        cell->hugeint_upper = ((duckdb_hugeint *)data)[at].upper;
        break;
    case DUCKVEP_CELL_VARCHAR: {
        duckdb_string_t *string = &((duckdb_string_t *)data)[at];
        cell->text = duckdb_string_t_data(string);
        cell->text_length = duckdb_string_t_length(*string);
        break;
    }
    case DUCKVEP_CELL_BOOLEAN: cell->boolean = ((bool *)data)[at]; break;
    default: break;
    }
}

#endif
