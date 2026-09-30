#include "duckvep_core_cells.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char duckvep_core_option_not_struct[] = "DuckVEP builder: options must be a STRUCT";

static bool huge_fractional(uint64_t lower, int64_t upper, uint8_t scale)
{
    uint64_t high = (uint64_t)upper, low = lower;
    if (upper < 0) {
        low = ~low + 1;
        high = ~high + (low == 0);
    }
    uint32_t limbs[4] = {(uint32_t)(high >> 32), (uint32_t)high,
        (uint32_t)(low >> 32), (uint32_t)low};
    for (uint8_t digit = 0; digit < scale; digit++) {
        uint64_t remainder = 0;
        for (size_t i = 0; i < 4; i++) {
            uint64_t value = (remainder << 32) | limbs[i];
            limbs[i] = (uint32_t)(value / 10);
            remainder = value % 10;
        }
        if (remainder) return true;
    }
    return false;
}

/* Every numeric kind as a long double; false for the others. */
static bool raw_number(const duckvep_cell_t *cell, long double *number)
{
    switch (cell->kind) {
    case DUCKVEP_CELL_TINYINT: case DUCKVEP_CELL_SMALLINT:
    case DUCKVEP_CELL_INTEGER: case DUCKVEP_CELL_BIGINT:
        *number = cell->i; return true;
    case DUCKVEP_CELL_UTINYINT: case DUCKVEP_CELL_USMALLINT:
    case DUCKVEP_CELL_UINTEGER: case DUCKVEP_CELL_UBIGINT:
        *number = cell->u; return true;
    case DUCKVEP_CELL_FLOAT: case DUCKVEP_CELL_DOUBLE:
        *number = cell->d; return true;
    case DUCKVEP_CELL_HUGEINT:
        *number = (long double)cell->hugeint_upper * 18446744073709551616.0L + cell->hugeint_lower;
        return true;
    default: return false;
    }
}

bool duckvep_core_cell_number(const duckvep_cell_t *cell, long double *number, bool *fractional)
{
    if (!cell || !cell->valid || !raw_number(cell, number)) return false;
    uint8_t scale = cell->scale;
    if (cell->kind == DUCKVEP_CELL_HUGEINT && scale)
        *fractional = huge_fractional(cell->hugeint_lower, cell->hugeint_upper, scale);
    if (scale) {
        long double base = 1;
        for (uint8_t i = 0; i < scale; i++) base *= 10;
        if (cell->kind != DUCKVEP_CELL_HUGEINT)
            *fractional = fmodl(*number, base) != 0;
        *number /= base;
    } else *fractional = isfinite(*number) && truncl(*number) != *number;
    return true;
}

bool duckvep_core_phase_allele(const duckvep_cell_t *cell, int32_t *out)
{
    long double number;
    if (cell->kind == DUCKVEP_CELL_VARCHAR) {
        size_t length = cell->text_length;
        if (!length || length >= 32) return false;
        char text[32];
        memcpy(text, cell->text, length);
        text[length] = '\0';
        char *end;
        long value = strtol(text, &end, 10);
        if (end != text + length || value < INT32_MIN || value > INT32_MAX) return false;
        *out = (int32_t)value;
        return true;
    }
    if (!raw_number(cell, &number)) return false;
    for (uint8_t i = 0; i < cell->scale; i++) number /= 10;
    if (!isfinite(number) || number < INT32_MIN - 0.5L || number >= INT32_MAX + 0.5L)
        return false;
    *out = (int32_t)roundl(number);
    return true;
}

bool duckvep_core_phase_flag(const duckvep_cell_t *cell, bool *result)
{
    if (cell->kind == DUCKVEP_CELL_BOOLEAN) {
        *result = cell->boolean;
        return true;
    }
    if (cell->kind == DUCKVEP_CELL_VARCHAR) {
        const char *data = cell->text;
        size_t length = cell->text_length;
        if (length == 4 || length == 5) {
            const char *expected = length == 4 ? "true" : "false";
            bool matches = true;
            for (size_t i = 0; i < length; i++) {
                char c = data[i];
                if (c >= 'A' && c <= 'Z') c = (char)(c + ('a' - 'A'));
                if (c != expected[i]) { matches = false; break; }
            }
            if (matches) { *result = length == 4; return true; }
        }
        if (length == 1 && (*data == '1' || *data == '0')) {
            *result = *data == '1'; return true;
        }
        return false;
    }
    duckvep_cell_t plain = *cell;
    int32_t value = 0;
    plain.scale = 0; /* a flag is never scaled */
    if (!duckvep_core_phase_allele(&plain, &value)) return false;
    *result = value != 0;
    return true;
}

size_t duckvep_core_option_index(const char *key, const char *const *names, size_t count)
{
    size_t at = 0;
    while (at < count && strcmp(key, names[at]) != 0) at++;
    return at;
}

bool duckvep_core_option_permitted(const char *key, size_t at, size_t count,
    const duckvep_core_option_kind_t *kinds, duckvep_core_field_type_t field,
    char *message, size_t message_size)
{
    bool permitted = at < count && (field == DUCKVEP_CORE_FIELD_SQLNULL ||
        (kinds[at] == DUCKVEP_CORE_OPTION_TEXT && field == DUCKVEP_CORE_FIELD_VARCHAR) ||
        (kinds[at] == DUCKVEP_CORE_OPTION_INTEGER && field == DUCKVEP_CORE_FIELD_INTEGER) ||
        (kinds[at] == DUCKVEP_CORE_OPTION_BOOLEAN && field == DUCKVEP_CORE_FIELD_BOOLEAN) ||
        (kinds[at] == DUCKVEP_CORE_OPTION_NUMERIC &&
         (field == DUCKVEP_CORE_FIELD_INTEGER || field == DUCKVEP_CORE_FIELD_NUMERIC)));
    if (permitted) return true;
    const char *type_name = at == count ? "" : kinds[at] == DUCKVEP_CORE_OPTION_TEXT ?
        "VARCHAR" : kinds[at] == DUCKVEP_CORE_OPTION_INTEGER ? "INTEGER" :
        kinds[at] == DUCKVEP_CORE_OPTION_BOOLEAN ? "BOOLEAN" : "numeric";
    if (at == count || kinds[at] == DUCKVEP_CORE_OPTION_TEXT)
        snprintf(message, message_size, "DuckVEP builder: %s option '%s'",
                 at == count ? "unknown" : "expected VARCHAR for", key);
    else
        snprintf(message, message_size, "DuckVEP builder: expected %s for option '%s'",
                 type_name, key);
    return false;
}
