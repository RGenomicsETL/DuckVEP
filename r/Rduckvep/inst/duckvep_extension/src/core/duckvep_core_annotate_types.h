/* The result columns of the _duckvep_annotate_* natives (LIST of STRUCT), shared by the hosts'
 * registrations so that the v1 and v2 result types cannot drift. */
#ifndef DUCKVEP_CORE_ANNOTATE_TYPES_H
#define DUCKVEP_CORE_ANNOTATE_TYPES_H

#include <stddef.h>

typedef enum {
	DUCKVEP_COLUMN_UINTEGER, DUCKVEP_COLUMN_VARCHAR, DUCKVEP_COLUMN_BOOLEAN,
	DUCKVEP_COLUMN_UBIGINT, DUCKVEP_COLUMN_UTINYINT
} duckvep_column_type_t;

typedef struct {
	const char *name;
	duckvep_column_type_t type;
} duckvep_result_column_t;

#define DUCKVEP_RESULT_COLUMNS_MAX 57

typedef enum {
	DUCKVEP_RESULT_RICH,              /* 29 columns; with_hgvs adds 7, with_projection 21 */
	DUCKVEP_RESULT_COMPACT,           /* 16 */
	DUCKVEP_RESULT_COMPACT_HGVS       /* 23 */
} duckvep_result_kind_t;

size_t duckvep_core_result_columns(duckvep_result_kind_t kind, int with_hgvs, int with_projection,
	duckvep_result_column_t out[DUCKVEP_RESULT_COLUMNS_MAX]);
const char *duckvep_core_column_type_name(duckvep_column_type_t type);

#endif
