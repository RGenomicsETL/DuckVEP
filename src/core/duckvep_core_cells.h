/* Host-neutral typed cells and option checks. A host reads one element of a
 * vector into a cell (its type switch is the only DuckDB-facing part); the row
 * logic in src/core converts and validates cells identically for every host.
 * No DuckDB symbol may appear in src/core. */
#ifndef DUCKVEP_CORE_CELLS_H
#define DUCKVEP_CORE_CELLS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
    DUCKVEP_CELL_UNSUPPORTED = 0, /* any type the native logic does not read, including SQLNULL */
    DUCKVEP_CELL_TINYINT, DUCKVEP_CELL_SMALLINT, DUCKVEP_CELL_INTEGER, DUCKVEP_CELL_BIGINT,
    DUCKVEP_CELL_UTINYINT, DUCKVEP_CELL_USMALLINT, DUCKVEP_CELL_UINTEGER, DUCKVEP_CELL_UBIGINT,
    DUCKVEP_CELL_FLOAT, DUCKVEP_CELL_DOUBLE, DUCKVEP_CELL_HUGEINT,
    DUCKVEP_CELL_VARCHAR, DUCKVEP_CELL_BOOLEAN
} duckvep_cell_kind_t;

/* A DECIMAL is its storage kind (TINYINT is never used; SMALLINT, INTEGER, BIGINT
 * or HUGEINT) plus a nonzero scale. Only the member matching `kind` is read. */
typedef struct {
    duckvep_cell_kind_t kind;
    uint8_t scale;
    bool decimal; /* a DECIMAL (kind is its storage kind); callers that accept only plain numbers refuse it */
    bool valid;
    int64_t i;
    uint64_t u;
    double d;
    uint64_t hugeint_lower;
    int64_t hugeint_upper;
    const char *text;
    size_t text_length;
    bool boolean;
} duckvep_cell_t;

/* The numeric value of a valid numeric cell, and whether it has a fractional
 * part. False for an invalid, VARCHAR, BOOLEAN or unsupported cell. */
bool duckvep_core_cell_number(const duckvep_cell_t *cell, long double *number, bool *fractional);

/* duckvep_phase_call inputs. An allele accepts any numeric kind (DECIMAL scaled,
 * rounded to an INTEGER) and a base-10 VARCHAR; a flag accepts BOOLEAN, the
 * strings true/false (any case) and 1/0, and any allele-convertible number. */
bool duckvep_core_phase_allele(const duckvep_cell_t *cell, int32_t *out);
bool duckvep_core_phase_flag(const duckvep_cell_t *cell, bool *out);

/* Option STRUCT fields: what a host must check before reading a field. */
typedef enum {
    DUCKVEP_CORE_OPTION_TEXT, DUCKVEP_CORE_OPTION_INTEGER,
    DUCKVEP_CORE_OPTION_NUMERIC, DUCKVEP_CORE_OPTION_BOOLEAN
} duckvep_core_option_kind_t;

typedef enum {
    DUCKVEP_CORE_FIELD_OTHER, DUCKVEP_CORE_FIELD_SQLNULL, DUCKVEP_CORE_FIELD_VARCHAR,
    DUCKVEP_CORE_FIELD_BOOLEAN,
    DUCKVEP_CORE_FIELD_INTEGER, /* TINYINT .. BIGINT and UTINYINT .. UBIGINT */
    DUCKVEP_CORE_FIELD_NUMERIC  /* HUGEINT, DECIMAL, FLOAT, DOUBLE */
} duckvep_core_field_type_t;

/* The index of `key` in names[count], or count when unknown. */
size_t duckvep_core_option_index(const char *key, const char *const *names, size_t count);

/* Writes the builder's error text into message and returns false when the field
 * `key` (index `at`, count meaning unknown) of type `field` is not permitted. */
bool duckvep_core_option_permitted(const char *key, size_t at, size_t count,
    const duckvep_core_option_kind_t *kinds, duckvep_core_field_type_t field,
    char *message, size_t message_size);

extern const char duckvep_core_option_not_struct[];

#endif
