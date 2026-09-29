#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_registration.h"
#include <stdlib.h>
#include <string.h>
DUCKDB_EXTENSION_EXTERN

static const char sv_sql[] =
"), parts AS (SELECT *, string_split(info, ';') AS tokens, "
"try_cast(try_cast(pos AS HUGEINT) + length(ref) AS BIGINT) AS nominal_start, "
"regexp_full_match(alt, '<[^>]+>') AS symbolic, "
"regexp_full_match(ref, '[ACGT]+') AND regexp_full_match(alt, '[ACGT]+') "
"AND length(alt)>length(ref) AND starts_with(alt,ref) AS literal FROM source), "
"fields AS (SELECT *, list_transform(tokens, lambda x: split_part(x,'=',1)) AS keys, "
"nullif(regexp_extract(info, '(?:^|;)END=([^;]*)', 1), '') AS end_text, "
"nullif(regexp_extract(info, '(?:^|;)SEQ=([^;]*)', 1), '') AS source_sequence, "
"nullif(regexp_extract(info, '(?:^|;)CIPOS=([^;]*)', 1), '') AS cipos_text, "
"nullif(regexp_extract(info, '(?:^|;)CIEND=([^;]*)', 1), '') AS ciend_text "
"FROM parts), bounds AS (SELECT *, CASE WHEN symbolic THEN try_cast(end_text AS BIGINT) "
"ELSE try_cast(nominal_start::HUGEINT-1 AS BIGINT) END AS nominal_end, "
"CASE WHEN cipos_text IS NOT NULL AND regexp_full_match(cipos_text,'[+-]?[0-9]+,[+-]?[0-9]+') "
"THEN try_cast(split_part(cipos_text,',',1) AS BIGINT) END AS cipos_lo, "
"CASE WHEN cipos_text IS NOT NULL AND regexp_full_match(cipos_text,'[+-]?[0-9]+,[+-]?[0-9]+') "
"THEN try_cast(split_part(cipos_text,',',2) AS BIGINT) END AS cipos_hi, "
"CASE WHEN ciend_text IS NOT NULL AND regexp_full_match(ciend_text,'[+-]?[0-9]+,[+-]?[0-9]+') "
"THEN try_cast(split_part(ciend_text,',',1) AS BIGINT) END AS ciend_lo, "
"CASE WHEN ciend_text IS NOT NULL AND regexp_full_match(ciend_text,'[+-]?[0-9]+,[+-]?[0-9]+') "
"THEN try_cast(split_part(ciend_text,',',2) AS BIGINT) END AS ciend_hi FROM fields), "
"geometry AS (SELECT *, try_cast(nominal_start::HUGEINT+cipos_lo AS BIGINT) AS outer_start, "
"try_cast(nominal_start::HUGEINT+cipos_hi AS BIGINT) AS inner_start, "
"try_cast(nominal_end::HUGEINT+ciend_lo AS BIGINT) AS inner_end, "
"try_cast(nominal_end::HUGEINT+ciend_hi AS BIGINT) AS outer_end FROM bounds) "
"SELECT event_index, CASE WHEN "
"(symbolic OR literal) AND pos BETWEEN 1 AND 2147483647 AND nominal_start BETWEEN 0 AND 2147483647 "
"AND nominal_end BETWEEN nominal_start::HUGEINT-1 AND 2147483647 "
"AND (NOT symbolic OR (regexp_full_match(ref,'[ACGTN]') AND regexp_full_match(end_text,'[1-9][0-9]*'))) "
"AND (cipos_text IS NULL OR (outer_start BETWEEN 0 AND 2147483647 AND inner_start BETWEEN 0 AND 2147483647 AND outer_start<=inner_start)) "
"AND (ciend_text IS NULL OR (inner_end BETWEEN 0 AND 2147483647 AND outer_end BETWEEN 0 AND 2147483647 AND inner_end<=outer_end)) "
"AND length(keys)=length(list_distinct(keys)) "
"THEN 'ok' ELSE 'unsupported_geometry' END AS status, "
"CASE WHEN symbolic THEN 'structural' WHEN literal THEN 'literal_insertion' ELSE 'unsupported' END AS mode, "
"nominal_start, nominal_end, outer_start, inner_start, inner_end, outer_end, "
"CASE WHEN literal THEN substring(alt,length(ref)+1) END AS inserted_sequence, source_sequence FROM geometry";

static void prepare_sv(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector name_vector = duckdb_data_chunk_get_vector(input, 0);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        duckdb_vector option = NULL;
        const char *const names[] = {"unused"};
        const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_TEXT};
        if (argc == 2 && !duckvep_builder_option_vectors(info,
            duckdb_data_chunk_get_vector(input, 1), row, names, kinds, 0, &option)) return;
        uint64_t *valid = duckdb_vector_get_validity(name_vector);
        char *name = valid && !duckdb_validity_row_is_valid(valid, row) ? NULL :
            duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(name_vector))[row]);
        duckvep_sql_text sql = {0};
        bool ok = name && duckvep_sql_append(&sql, "WITH source AS (SELECT event_index, pos, ref, alt, info FROM ") &&
            duckvep_sql_relation(&sql, name) && duckvep_sql_append(&sql, sv_sql);
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckvep_builder_set_error(info, "duckvep_prepare_sv_geometry_sql: invalid relation name or allocation failure");
        duckvep_sql_free(&sql);
        duckvep_budget_free(name);
        if (!ok) return;
    }
}

static const char str_sql[] =
"), joined AS (SELECT s.*, r.reference_sequence, r.ref_rows FROM source s LEFT JOIN "
"(SELECT event_index, count(*) AS ref_rows, CASE WHEN count(*) > 1 THEN NULL ELSE min(reference_sequence) END AS reference_sequence FROM " ;

static const char str_tail[] =
" GROUP BY event_index) r USING (event_index)), parsed AS (SELECT *, string_split(alt,',') alts, "
"string_split(info,';') tokens, string_split(format,':') fields, string_split(\"sample\",':') entries "
"FROM joined), metadata AS (SELECT *, list_transform(tokens, lambda x: split_part(x,'=',1)) keys, "
"nullif(regexp_extract(info,'(?:^|;)REF=([^;]*)',1),'') ref_text, "
"nullif(regexp_extract(info,'(?:^|;)RL=([^;]*)',1),'') rl_text, "
"nullif(regexp_extract(info,'(?:^|;)RU=([^;]*)',1),'') unit, "
"nullif(regexp_extract(info,'(?:^|;)END=([^;]*)',1),'') end_text, "
"list_extract(alts, try_cast(alt_index AS INTEGER)) selected FROM parsed), counts AS (SELECT *, "
"try_cast(ref_text AS DOUBLE) ref_count, try_cast(rl_text AS DOUBLE) ref_length, "
"try_cast(regexp_extract(selected,'^<STR([0-9]+)>$',1) AS DOUBLE) alt_count, "
"CASE WHEN list_has_all(fields,['CN','CI']) THEN 'CN' ELSE 'REPCN' END count_key, "
"CASE WHEN list_has_all(fields,['CN','CI']) THEN 'CI' ELSE 'REPCI' END ci_key "
"FROM metadata), samples AS (SELECT *, "
"string_split(replace(coalesce(list_extract(entries,list_position(fields,'GT')),''),'|','/'),'/') gt, "
"string_split(coalesce(list_extract(entries,list_position(fields,count_key)),''),'/') copies, "
"string_split(coalesce(list_extract(entries,list_position(fields,ci_key)),''),'/') ranges, "
"string_split(coalesce(list_extract(entries,list_position(fields,'SO')),''),'/') support "
"FROM counts), calls AS (SELECT *, "
"list_filter(range(1,length(gt)+1),lambda i: gt[i]=cast(try_cast(alt_index AS BIGINT) AS VARCHAR)) selected_calls "
"FROM samples), ";

static const char str_tail2[] =
"decisions AS (SELECT *, CASE "
"WHEN ref_rows > 1 THEN 'ambiguous_reference_sequence' "
"WHEN info IS NULL OR format IS NULL OR \"sample\" IS NULL OR ref IS NULL OR alt IS NULL OR reference_sequence IS NULL THEN 'missing_field' "
"WHEN alt_index IS NULL OR NOT isfinite(try_cast(alt_index AS DOUBLE)) OR "
"try_cast(alt_index AS DOUBLE) != floor(try_cast(alt_index AS DOUBLE)) OR "
"try_cast(alt_index AS DOUBLE) NOT BETWEEN 1 AND length(alts) THEN 'alt_index' "
"WHEN NOT regexp_full_match(ref,'[ACGT]') OR NOT regexp_full_match(reference_sequence,'[ACGT]+') THEN 'literal_reference' "
"WHEN length(keys)!=length(list_distinct(keys)) THEN 'duplicate_info' "
"WHEN NOT list_has_all(keys,['REF','RL','RU','END']) THEN 'missing_info' "
"WHEN NOT regexp_full_match(unit,'[ACGT]+') OR "
"NOT regexp_full_match(ref_text,'(0|[1-9][0-9]*)') OR "
"NOT regexp_full_match(rl_text,'(0|[1-9][0-9]*)') OR "
"NOT regexp_full_match(end_text,'(0|[1-9][0-9]*)') OR "
"NOT isfinite(ref_count) OR NOT isfinite(ref_length) OR "
"NOT isfinite(try_cast(end_text AS DOUBLE)) THEN 'metadata_syntax' "
"WHEN ref_length>5000 THEN 'allele_capacity' "
"WHEN ref_count*length(unit)!=ref_length OR ref_length!=length(reference_sequence) OR "
"repeat(unit,try_cast(ref_count AS INTEGER))!=reference_sequence THEN 'reference_mismatch' "
"WHEN NOT regexp_full_match(selected,'<STR(0|[1-9][0-9]*)>') THEN 'unsupported_alt' "
"WHEN NOT isfinite(alt_count) THEN 'metadata_syntax' "
"WHEN length(fields)!=length(entries) OR length(fields)!=length(list_distinct(fields)) THEN 'format_shape' "
"WHEN NOT list_has_all(fields,['GT','SO']) OR NOT (list_has_all(fields,['CN','CI']) OR "
"list_has_all(fields,['REPCN','REPCI'])) THEN 'missing_format' "
"WHEN list_has_all(fields,['CN','CI','REPCN','REPCI']) THEN 'ambiguous_count_fields' "
"WHEN length(gt)!=length(copies) OR length(gt)!=length(ranges) OR length(gt)!=length(support) THEN 'format_shape' "
"WHEN NOT coalesce(list_bool_and(list_transform(gt,lambda x: x='.' OR "
"(regexp_full_match(x,'(0|[1-9][0-9]*)') AND try_cast(x AS DOUBLE) BETWEEN 0 AND length(alts)))),false) THEN 'genotype_index' "
"WHEN length(selected_calls)=0 THEN 'uncalled_alt' "
"WHEN NOT coalesce(list_bool_and(list_transform(selected_calls,lambda i: "
"regexp_full_match(copies[i],'(0|[1-9][0-9]*)') AND try_cast(copies[i] AS DOUBLE)=alt_count)),false) THEN 'count_mismatch' "
"WHEN NOT coalesce(list_bool_or(list_transform(selected_calls,lambda i: support[i]='SPANNING' AND "
"ranges[i]=cast(try_cast(alt_count AS BIGINT) AS VARCHAR)||'-'||cast(try_cast(alt_count AS BIGINT) AS VARCHAR))),false) "
"THEN 'estimated_count' "
"WHEN alt_count*length(unit)>5000 THEN 'allele_capacity' ELSE 'exact' END reason FROM calls) "
"SELECT event_index, CASE WHEN reason='exact' THEN 'ok' "
"WHEN reason IN ('missing_field','missing_info','missing_format') THEN 'incomplete' "
"WHEN reason IN ('unsupported_alt','uncalled_alt','estimated_count','allele_capacity') THEN 'summary_only' "
"ELSE 'invalid' END status, reason, "
"struct_pack(info := info, format := format, sample := \"sample\", ref := ref, alt := alt, "
"reference_sequence := reference_sequence, alt_index := alt_index) AS \"source\", "
"CASE WHEN reason NOT IN ('ambiguous_reference_sequence','missing_field','alt_index','literal_reference','duplicate_info',"
"'missing_info','metadata_syntax','reference_mismatch','allele_capacity') "
"OR (reason='allele_capacity' AND ref_length<=5000) "
"THEN struct_pack(unit := unit, count := ref_count) END reference_components, "
"CASE WHEN reason='exact' THEN struct_pack(unit := unit, count := alt_count) END alternate_components, "
"reason='exact' sequence_exact FROM decisions";

static void prepare_str(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    const char *const names[] = {"unused"};
    const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_TEXT};
    duckdb_vector args[3];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        duckdb_vector option = NULL;
        if (argc == 3 && !duckvep_builder_option_vectors(info, args[2], row, names, kinds, 0, &option)) return;
        char *tables[2] = {0};
        bool ok = true;
        for (idx_t i = 0; i < 2; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) { ok = false; break; }
            tables[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!tables[i]) { ok = false; break; }
        }
        duckvep_sql_text sql = {0};
        if (ok) ok = duckvep_sql_append(&sql, "WITH source AS (SELECT event_index, info, format, \"sample\", ref, alt, alt_index FROM ") &&
            duckvep_sql_relation(&sql, tables[0]) && duckvep_sql_append(&sql, str_sql) &&
            duckvep_sql_relation(&sql, tables[1]) && duckvep_sql_append(&sql, str_tail) &&
            duckvep_sql_append(&sql, str_tail2);
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else duckvep_builder_set_error(info, "duckvep_prepare_expansionhunter_sql: invalid relation name or allocation failure");
        duckvep_sql_free(&sql);
        for (size_t i = 0; i < 2; i++) duckvep_budget_free(tables[i]);
        if (!ok) return;
    }
}

bool register_duckvep_prepare_sql(duckdb_connection connection) {
    return duckvep_register_builder(connection, "duckvep_prepare_sv_geometry_sql", 1, prepare_sv) &&
        duckvep_register_builder(connection, "duckvep_prepare_expansionhunter_sql", 2, prepare_str);
}
