/* Relation-oriented DuckVEP SQL surface and static annotation metadata. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckvep_so.h"
#include "duckvep_codon.h"
#include "duckvep_sql.h"
#include "duckvep_builder.h"

typedef struct {
	idx_t offset;
} duckvep_so_scan_t;

static void
duckvep_so_terms_bind(duckdb_bind_info info)
{
	duckdb_logical_type utinyint_type, ubigint_type, varchar_type;

	utinyint_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
	varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_bind_add_result_column(info, "bit_index", utinyint_type);
	duckdb_bind_add_result_column(info, "consequence_mask", ubigint_type);
	duckdb_bind_add_result_column(info, "consequence", varchar_type);
	duckdb_bind_add_result_column(info, "impact_code", utinyint_type);
	duckdb_bind_add_result_column(info, "impact", varchar_type);
	duckdb_bind_add_result_column(info, "severity_rank", utinyint_type);
	duckdb_bind_add_result_column(info, "evaluator_tier", utinyint_type);
	duckdb_destroy_logical_type(&utinyint_type);
	duckdb_destroy_logical_type(&ubigint_type);
	duckdb_destroy_logical_type(&varchar_type);
}

static void
duckvep_so_scan_destroy(void *data)
{
	duckdb_free(data);
}

static void
duckvep_so_terms_init(duckdb_init_info info)
{
	duckvep_so_scan_t *scan;

	scan = duckdb_malloc(sizeof(*scan));
	if (scan == NULL) {
		duckdb_init_set_error(info,
		    "duckvep_so_terms: failed to allocate scan state");
		return;
	}
	scan->offset = 0;
	duckdb_init_set_init_data(info, scan, duckvep_so_scan_destroy);
}

static void
duckvep_so_terms_scan(duckdb_function_info info, duckdb_data_chunk output)
{
	duckvep_so_scan_t *scan;
	duckdb_vector bit_vector, mask_vector, consequence_vector;
	duckdb_vector impact_code_vector, impact_vector, rank_vector, tier_vector;
	uint8_t *bits, *impact_codes, *ranks, *tiers;
	uint64_t *masks;
	idx_t count, vector_size;

	scan = duckdb_function_get_init_data(info);
	if (scan == NULL) {
		duckdb_data_chunk_set_size(output, 0);
		return;
	}
	bit_vector = duckdb_data_chunk_get_vector(output, 0);
	mask_vector = duckdb_data_chunk_get_vector(output, 1);
	consequence_vector = duckdb_data_chunk_get_vector(output, 2);
	impact_code_vector = duckdb_data_chunk_get_vector(output, 3);
	impact_vector = duckdb_data_chunk_get_vector(output, 4);
	rank_vector = duckdb_data_chunk_get_vector(output, 5);
	tier_vector = duckdb_data_chunk_get_vector(output, 6);
	bits = duckdb_vector_get_data(bit_vector);
	masks = duckdb_vector_get_data(mask_vector);
	impact_codes = duckdb_vector_get_data(impact_code_vector);
	ranks = duckdb_vector_get_data(rank_vector);
	tiers = duckdb_vector_get_data(tier_vector);
	vector_size = duckdb_vector_size();
	count = 0;
	while (count < vector_size && scan->offset < DUCKVEP_SO_BIT_COUNT) {
		duckvep_so_bit_t bit;
		duckvep_impact_t impact;
		const char *name;

		bit = (duckvep_so_bit_t)scan->offset;
		impact = duckvep_so_bit_impact(bit);
		name = duckvep_so_name(bit);
		bits[count] = (uint8_t)bit;
		masks[count] = DUCKVEP_SO(bit);
		impact_codes[count] = (uint8_t)impact;
		ranks[count] = duckvep_so_rank(bit);
		tiers[count] = duckvep_so_tier(bit);
		duckdb_vector_assign_string_element(consequence_vector, count,
		    name != NULL ? name : "");
		duckdb_vector_assign_string_element(impact_vector, count,
		    duckvep_impact_name(impact));
		count++;
		scan->offset++;
	}
	duckdb_data_chunk_set_size(output, count);
}

static bool
duckvep_register_so_terms(duckdb_connection connection)
{
	duckdb_table_function function;
	duckdb_state state;

	function = duckdb_create_table_function();
	duckdb_table_function_set_name(function, "duckvep_so_terms");
	duckdb_table_function_set_bind(function, duckvep_so_terms_bind);
	duckdb_table_function_set_init(function, duckvep_so_terms_init);
	duckdb_table_function_set_function(function, duckvep_so_terms_scan);
	state = duckdb_register_table_function(connection, function);
	duckdb_destroy_table_function(&function);
	return state == DuckDBSuccess;
}

static bool
duckvep_register_annotate_relation(duckhts_registration_t *registration)
{
	static const char *const sql[] = {
		"CREATE OR REPLACE MACRO duckvep_annotate(events_table, model_name, ",
		"hgvs := false, upstream_distance := 5000, downstream_distance := 5000, ",
		"rich := false) AS TABLE ",
		"WITH parameters AS MATERIALIZED (SELECT CAST(model_name AS VARCHAR) AS model_name, ",
		"coalesce(CAST(hgvs AS BOOLEAN), false) AS include_hgvs, ",
		"coalesce(CAST(rich AS BOOLEAN), false) AS include_rich, ",
		"CAST(upstream_distance AS UBIGINT) AS upstream_distance, ",
		"CAST(downstream_distance AS UBIGINT) AS downstream_distance), ",
		"source AS (SELECT CAST(e.event_index AS UBIGINT) AS __duckvep_event_index, ",
		"CAST(e.seq_region AS UINTEGER) AS __duckvep_seq_region, ",
		"CAST(e.position AS UBIGINT) AS __duckvep_position, ",
		"CAST(e.reference AS VARCHAR) AS __duckvep_reference, ",
		"CAST(e.alternate AS VARCHAR) AS __duckvep_alternate, ",
		"CAST(e.end_position AS UBIGINT) AS __duckvep_end_position, ",
		"upper(CAST(e.structural_type AS VARCHAR)) AS __duckvep_explicit_structural_type, ",
		"upper(CAST(e.copy_change AS VARCHAR)) AS __duckvep_explicit_copy_change, ",
		"CAST(e.mate_seq_region AS UINTEGER) AS __duckvep_mate_seq_region, ",
		"CAST(e.mate_position AS UBIGINT) AS __duckvep_mate_position ",
		"FROM query_table(events_table) e), ",
		"typed_base AS MATERIALIZED (SELECT *, ",
		"CASE __duckvep_alternate ",
		"WHEN '<DEL>' THEN 'DEL' WHEN '<DUP>' THEN 'DUP' ",
		"WHEN '<TDUP>' THEN 'TDUP' WHEN '<STR>' THEN 'STR' ",
		"WHEN '<INV>' THEN 'INV' WHEN '<INS>' THEN 'INS' ",
		"WHEN '<CNV>' THEN 'CNV' WHEN '<UNKNOWN>' THEN 'UNKNOWN' ",
		"ELSE NULL END AS __duckvep_symbolic_structural_type, ",
		"__duckvep_alternate = '<*>' AS __duckvep_unspecified_alt, ",
		"__duckvep_alternate IN ('<NON_REF>', '*', '.') ",
		"AS __duckvep_non_variant FROM source), ",
		"typed AS MATERIALIZED (SELECT *, ",
		"coalesce(__duckvep_explicit_structural_type, ",
		"__duckvep_symbolic_structural_type) AS __duckvep_structural_type ",
		"FROM typed_base), ",
		"classified AS MATERIALIZED (SELECT *, CASE ",
		"WHEN __duckvep_non_variant THEN 'non_variant' ",
		"WHEN __duckvep_unspecified_alt THEN 'small_variant' ",
		"WHEN __duckvep_mate_seq_region IS NOT NULL OR __duckvep_mate_position IS NOT NULL ",
		"THEN 'breakend' WHEN __duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL ",
		"OR (__duckvep_alternate IS NOT NULL AND ",
		"(starts_with(__duckvep_alternate, '<') OR contains(__duckvep_alternate, '[') ",
		"OR contains(__duckvep_alternate, ']'))) ",
		"OR (__duckvep_end_position IS NOT NULL AND ",
		"(__duckvep_reference IS NULL OR __duckvep_alternate IS NULL)) ",
		"THEN 'structural_variant' ",
		"ELSE 'small_variant' END AS __duckvep_event_kind, ",
		"coalesce(__duckvep_explicit_copy_change, CASE __duckvep_structural_type ",
		"WHEN 'DEL' THEN 'LOSS' WHEN 'DUP' THEN 'GAIN' WHEN 'TDUP' THEN 'GAIN' ",
		"ELSE 'UNKNOWN' END) AS __duckvep_copy_change FROM typed), ",
		"validation AS MATERIALIZED (SELECT CASE ",
		"WHEN (SELECT model_name IS NULL OR model_name = '' FROM parameters) ",
		"THEN error('duckvep_annotate: model_name must be non-empty') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_index IS NULL) ",
		"THEN error('duckvep_annotate: event_index is required') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_seq_region IS NULL ",
		"OR __duckvep_position IS NULL OR __duckvep_position = 0) ",
		"THEN error('duckvep_annotate: seq_region and positive one-based position are required') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'non_variant' ",
		"AND (__duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL ",
		"OR __duckvep_mate_seq_region IS NOT NULL ",
		"OR __duckvep_mate_position IS NOT NULL)) ",
		"THEN error('duckvep_annotate: a non-variant gVCF allele cannot carry structural or breakend fields') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_unspecified_alt ",
		"AND (__duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL ",
		"OR __duckvep_mate_seq_region IS NOT NULL ",
		"OR __duckvep_mate_position IS NOT NULL)) ",
		"THEN error('duckvep_annotate: <*> cannot carry structural or breakend fields') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE ",
		"(__duckvep_mate_seq_region IS NULL) != (__duckvep_mate_position IS NULL)) ",
		"THEN error('duckvep_annotate: a breakend requires both mate_seq_region and mate_position') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'breakend' AND ",
		"(__duckvep_end_position IS NOT NULL OR __duckvep_explicit_structural_type IS NOT NULL ",
		"OR __duckvep_explicit_copy_change IS NOT NULL)) ",
		"THEN error('duckvep_annotate: breakend mate coordinates cannot be mixed with single-locus structural fields') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'breakend' ",
		"AND __duckvep_mate_position = 0) ",
		"THEN error('duckvep_annotate: mate_position must be positive and one-based') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'small_variant' ",
		"AND (__duckvep_reference IS NULL OR __duckvep_reference = '' ",
		"OR __duckvep_alternate IS NULL OR __duckvep_alternate = '')) ",
		"THEN error('duckvep_annotate: literal small variants require reference and alternate') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE ",
		"__duckvep_explicit_structural_type IS NOT NULL ",
		"AND __duckvep_symbolic_structural_type IS NOT NULL ",
		"AND __duckvep_explicit_structural_type <> __duckvep_symbolic_structural_type) ",
		"THEN error('duckvep_annotate: symbolic alternate and structural_type disagree') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'structural_variant' ",
		"AND (__duckvep_structural_type IS NULL OR __duckvep_end_position IS NULL)) ",
		"THEN error('duckvep_annotate: a structural event requires end_position and a supported structural_type') ",
		"WHEN EXISTS (SELECT 1 FROM classified WHERE __duckvep_event_kind = 'structural_variant' ",
		"AND __duckvep_structural_type NOT IN ('DEL','DUP','TDUP','STR','INV','INS','CNV','UNKNOWN')) ",
		"THEN error('duckvep_annotate: unsupported structural_type') ",
		"ELSE true END AS valid), ",
		"validated AS (SELECT classified.* FROM classified CROSS JOIN validation ",
		"WHERE validation.valid), ",
		"small_events AS (SELECT * FROM validated ",
		"WHERE __duckvep_event_kind = 'small_variant'), ",
		"structural_events AS (SELECT * FROM validated ",
		"WHERE __duckvep_event_kind = 'structural_variant'), ",
		"breakend_events AS (SELECT * FROM validated ",
		"WHERE __duckvep_event_kind = 'breakend'), ",
		"small_compact_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_compact(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE (NOT p.include_hgvs OR __duckvep_unspecified_alt) ",
		"AND NOT p.include_rich), ",
		"small_hgvs_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_hgvs(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE p.include_hgvs AND NOT __duckvep_unspecified_alt ",
		"AND NOT p.include_rich), ",
		"small_rich_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_rich(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE (NOT p.include_hgvs OR __duckvep_unspecified_alt) ",
		"AND p.include_rich), ",
		"small_rich_hgvs_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_small_rich_hgvs(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_reference, __duckvep_alternate, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM small_events e CROSS JOIN parameters p ",
		"WHERE p.include_hgvs AND NOT __duckvep_unspecified_alt ",
		"AND p.include_rich), ",
		"structural_compact_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_structural_compact(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_end_position, __duckvep_structural_type, ",
		"__duckvep_copy_change, p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM structural_events e CROSS JOIN parameters p WHERE NOT p.include_rich), ",
		"structural_rich_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_structural_rich(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_end_position, __duckvep_structural_type, ",
		"__duckvep_copy_change, p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM structural_events e CROSS JOIN parameters p WHERE p.include_rich), ",
		"breakend_compact_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_breakend_compact(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_mate_seq_region, __duckvep_mate_position, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM breakend_events e CROSS JOIN parameters p WHERE NOT p.include_rich), ",
		"breakend_rich_raw AS (SELECT e.*, ",
		"unnest(_duckvep_annotate_breakend_rich(p.model_name, __duckvep_seq_region, ",
		"__duckvep_position, __duckvep_mate_seq_region, __duckvep_mate_position, ",
		"p.upstream_distance, p.downstream_distance)) AS annotation ",
		"FROM breakend_events e CROSS JOIN parameters p WHERE p.include_rich), ",
		"small_compact AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* ",
		"FROM small_compact_raw), ",
		"small_hgvs AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* FROM small_hgvs_raw), ",
		"small_rich AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason FROM small_rich_raw), ",
		"small_rich_hgvs AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason ",
		"FROM small_rich_hgvs_raw), ",
		"structural_compact AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* ",
		"FROM structural_compact_raw), ",
		"structural_rich AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason ",
		"FROM structural_rich_raw), ",
		"breakend_compact AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, annotation.* ",
		"FROM breakend_compact_raw), ",
		"breakend_rich AS (SELECT __duckvep_event_index AS event_index, ",
		"__duckvep_event_kind AS duckvep_event_kind, ",
		"annotation.* EXCLUDE (status, reason), annotation.status AS duckvep_status, ",
		"annotation.reason AS duckvep_reason ",
		"FROM breakend_rich_raw) ",
		"SELECT * FROM small_compact UNION ALL BY NAME SELECT * FROM small_hgvs ",
		"UNION ALL BY NAME SELECT * FROM small_rich ",
		"UNION ALL BY NAME SELECT * FROM small_rich_hgvs ",
		"UNION ALL BY NAME SELECT * FROM structural_compact ",
		"UNION ALL BY NAME SELECT * FROM structural_rich ",
		"UNION ALL BY NAME SELECT * FROM breakend_compact ",
		"UNION ALL BY NAME SELECT * FROM breakend_rich"
	};

	return duckhts_register_sql_parts(registration, sql,
	    sizeof(sql) / sizeof(sql[0]));
}

/* SQL presentation uses the kernel's immutable genetic-code authority. The
 * string is indexed by T/C/A/G at each of the three codon positions. */
static void
duckvep_projection_code(duckdb_function_info info, duckdb_data_chunk input,
	duckdb_vector output)
{
	duckdb_vector vector = duckdb_data_chunk_get_vector(input, 0);
	const uint8_t *codes = duckdb_vector_get_data(vector);
	uint64_t *validity = duckdb_vector_get_validity(vector);
	idx_t count = duckdb_data_chunk_get_size(input);

	duckdb_vector_ensure_validity_writable(output);
	for (idx_t row = 0; row < count; row++) {
		const char *amino_acids;

		if (validity != NULL && !duckdb_validity_row_is_valid(validity, row)) {
			duckdb_validity_set_row_invalid(duckdb_vector_get_validity(output), row);
			continue;
		}
		amino_acids = duckvep_codon_table_amino_acids((duckvep_codon_table_t)codes[row]);
		if (amino_acids == NULL) {
			duckdb_scalar_function_set_error(info,
			    "duckvep_transcript_projection: unsupported genetic code");
			return;
		}
		duckdb_vector_assign_string_element_len(output, row, amino_acids, 64);
	}
}

static bool
duckvep_register_projection_code(duckdb_connection connection)
{
	duckdb_scalar_function function = duckdb_create_scalar_function();
	duckdb_logical_type code_type = duckdb_create_logical_type(DUCKDB_TYPE_UTINYINT);
	duckdb_logical_type text_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
	duckdb_state state;

	duckdb_scalar_function_set_name(function, "__duckvep_projection_code");
	duckdb_scalar_function_add_parameter(function, code_type);
	duckdb_scalar_function_set_return_type(function, text_type);
	duckdb_scalar_function_set_function(function, duckvep_projection_code);
	state = duckdb_register_scalar_function(connection, function);
	duckdb_destroy_scalar_function(&function);
	duckdb_destroy_logical_type(&code_type);
	duckdb_destroy_logical_type(&text_type);
	return state == DuckDBSuccess;
}

bool
register_duckvep_sql_kernels(duckhts_registration_t *registration)
{
	if (!duckvep_register_so_terms(registration->connection)) {
		return duckhts_registration_error(registration,
		    "DuckHTS could not register duckvep_so_terms");
	}
	if (!duckvep_register_projection_code(registration->connection)) {
		return duckhts_registration_error(registration,
		    "DuckHTS could not register __duckvep_projection_code");
	}
	if (!duckvep_register_phase_kernels(registration->connection)) {
		return duckhts_registration_error(registration,
		    "DuckHTS could not register DuckVEP phase kernels");
	}
	return true;
}

static bool
duckvep_projection_table(duckvep_sql_text *sql, const char *name)
{
    const char *dot = strchr(name, '.');
    if (dot) {
        size_t length = (size_t)(dot - name);
        char *schema = malloc(length + 1);
        if (!schema) return false;
        memcpy(schema, name, length);
        schema[length] = 0;
        bool ok = duckvep_sql_identifier(sql, schema) && duckvep_sql_append(sql, ".") &&
            duckvep_sql_identifier(sql, dot + 1);
        free(schema);
        return ok;
    }
    return duckvep_sql_identifier(sql, name);
}

#include "duckvep_projection_template.h"

static void
duckvep_projection_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output)
{
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[4];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *names[3] = {0};
        bool ok = true;
        for (idx_t i = 0; i < 3; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) { ok = false; break; }
            names[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!names[i] || !*names[i]) { ok = false; break; }
        }
        if (ok && argc == 4) {
            duckdb_vector field = NULL;
            const char *const keys[] = {"unused"};
            const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_TEXT};
            ok = duckvep_builder_option_vectors(info, args[3], row, keys, kinds, 0, &field);
            if (!ok) {
                for (idx_t i = 0; i < 3; i++) free(names[i]);
                return;
            }
        }
        duckvep_sql_text sql = {0};
        for (size_t i = 0; ok && i < sizeof(duckvep_projection_parts) / sizeof(*duckvep_projection_parts); i++) {
            const char *part = duckvep_projection_parts[i];
            while (ok && *part) {
                const char *mark = strstr(part, "__DUCKVEP_");
                if (!mark) { ok = duckvep_sql_append(&sql, part); break; }
                size_t length = (size_t)(mark - part);
                char *prefix = malloc(length + 1);
                if (!prefix) { ok = false; break; }
                memcpy(prefix, part, length); prefix[length] = 0;
                ok = duckvep_sql_append(&sql, prefix); free(prefix);
                const char *token = NULL; idx_t which = 0;
                if (!strncmp(mark, "__DUCKVEP_EVENTS_TABLE__", strlen("__DUCKVEP_EVENTS_TABLE__"))) token = "__DUCKVEP_EVENTS_TABLE__";
                else if (!strncmp(mark, "__DUCKVEP_ANNOTATIONS_TABLE__", strlen("__DUCKVEP_ANNOTATIONS_TABLE__"))) { token = "__DUCKVEP_ANNOTATIONS_TABLE__"; which = 1; }
                else if (!strncmp(mark, "__DUCKVEP_TRANSCRIPTS_TABLE__", strlen("__DUCKVEP_TRANSCRIPTS_TABLE__"))) { token = "__DUCKVEP_TRANSCRIPTS_TABLE__"; which = 2; }
                else ok = false;
                if (ok) ok = duckvep_projection_table(&sql, names[which]);
                if (token) part = mark + strlen(token);
            }
        }
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckdb_scalar_function_set_error(info, "duckvep_transcript_projection_sql: invalid table name, options, or allocation failure");
        duckvep_sql_free(&sql);
        for (idx_t i = 0; i < 3; i++) free(names[i]);
        if (!ok) return;
    }
}

static bool
duckvep_register_projection_relation(duckhts_registration_t *registration)
{
    return duckvep_register_builder(registration->connection,
        "duckvep_transcript_projection_sql", 3, duckvep_projection_builder);
}

bool
register_duckvep_sql_functions(duckhts_registration_t *registration)
{
	return duckvep_register_repeat_alleles(registration->connection) &&
	    duckvep_register_annotate_relation(registration) &&
	    duckvep_register_projection_relation(registration);
}
