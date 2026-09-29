#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_registration.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
DUCKDB_EXTENSION_EXTERN

/*
 * BND record identity, breakend-gene evidence and structural HGVS builders.
 *
 * Every builder returns SQL text. The SQL keeps one output row per physical
 * source record: a mate is looked up, never merged. Reasons are stable
 * lowercase identifiers documented in docs/structural-identity-hgvs.md.
 */

typedef enum { PIECE_TEXT, PIECE_RELATION, PIECE_INTEGER } piece_kind;
typedef struct {
    piece_kind kind;
    const char *text; /* PIECE_TEXT: SQL; PIECE_RELATION: relation ordinal as "0"/"1" */
} piece;

static bool assemble(duckvep_sql_text *sql, const piece *pieces, size_t count,
                     char *const *relations, int64_t integer) {
    for (size_t i = 0; i < count; i++) {
        bool ok;
        if (pieces[i].kind == PIECE_TEXT) ok = duckvep_sql_append(sql, pieces[i].text);
        else if (pieces[i].kind == PIECE_RELATION)
            ok = duckvep_sql_relation(sql, relations[pieces[i].text[0] - '0']);
        else {
            char digits[32];
            snprintf(digits, sizeof(digits), "%lld", (long long)integer);
            ok = duckvep_sql_append(sql, digits);
        }
        if (!ok) return false;
    }
    return true;
}

/* Reads a bounded integer option; absent or NULL keeps the default. */
static bool integer_option(duckdb_function_info info, duckdb_vector value, idx_t row,
                           const char *name, int64_t low, int64_t high, int64_t *out) {
    if (!value) return true;
    uint64_t *validity = duckdb_vector_get_validity(value);
    if (validity && !duckdb_validity_row_is_valid(validity, row)) return true;
    duckdb_logical_type type = duckdb_vector_get_column_type(value);
    duckdb_type id = duckdb_get_type_id(type);
    duckdb_destroy_logical_type(&type);
    void *data = duckdb_vector_get_data(value);
    int64_t parsed = 0;
    bool overflow = false;
    switch (id) {
    case DUCKDB_TYPE_TINYINT: parsed = ((int8_t *)data)[row]; break;
    case DUCKDB_TYPE_SMALLINT: parsed = ((int16_t *)data)[row]; break;
    case DUCKDB_TYPE_INTEGER: parsed = ((int32_t *)data)[row]; break;
    case DUCKDB_TYPE_BIGINT: parsed = ((int64_t *)data)[row]; break;
    case DUCKDB_TYPE_UTINYINT: parsed = ((uint8_t *)data)[row]; break;
    case DUCKDB_TYPE_USMALLINT: parsed = ((uint16_t *)data)[row]; break;
    case DUCKDB_TYPE_UINTEGER: parsed = ((uint32_t *)data)[row]; break;
    case DUCKDB_TYPE_UBIGINT: {
        uint64_t big = ((uint64_t *)data)[row];
        if (big > (uint64_t)INT64_MAX) overflow = true; else parsed = (int64_t)big;
        break;
    }
    default: overflow = true; break;
    }
    if (overflow || parsed < low || parsed > high) {
        char message[160];
        snprintf(message, sizeof(message), "DuckVEP builder: option '%s' must be between %lld and %lld",
                 name, (long long)low, (long long)high);
        duckdb_scalar_function_set_error(info, message);
        return false;
    }
    *out = parsed;
    return true;
}

typedef struct {
    const char *name;
    size_t relations;
    const piece *pieces;
    size_t piece_count;
    bool has_max_span;
} builder_spec;

static void run_builder(const builder_spec *spec, duckdb_function_info info,
                        duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[3];
    for (idx_t i = 0; i < argc && i < 3; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        const char *const names[] = {"max_span"};
        const duckvep_option_kind kinds[] = {DUCKVEP_OPTION_INTEGER};
        duckdb_vector option_fields[1] = {NULL};
        int64_t max_span = 5000;
        char *relations[2] = {0};
        bool ok = true;
        if (argc == spec->relations + 1) {
            ok = duckvep_builder_option_vectors(info, args[spec->relations], row, names, kinds,
                                                spec->has_max_span ? 1u : 0u, option_fields);
            if (ok && spec->has_max_span)
                ok = integer_option(info, option_fields[0], row, "max_span", 1, 60000, &max_span);
            if (!ok) return;
        }
        for (size_t i = 0; ok && i < spec->relations; i++) {
            uint64_t *valid = duckdb_vector_get_validity(args[i]);
            if (valid && !duckdb_validity_row_is_valid(valid, row)) { ok = false; break; }
            relations[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            if (!relations[i]) ok = false;
        }
        duckvep_sql_text sql = {0};
        if (ok) ok = assemble(&sql, spec->pieces, spec->piece_count, relations, max_span);
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        else {
            char message[160];
            snprintf(message, sizeof(message),
                     "%s: invalid relation name or allocation failure", spec->name);
            duckdb_scalar_function_set_error(info, message);
        }
        duckvep_sql_free(&sql);
        for (size_t i = 0; i < 2; i++) duckvep_budget_free(relations[i]);
        if (!ok) return;
    }
}

/* ---- BND record identity -------------------------------------------------- */

static const char pairs_head[] =
"WITH source AS (SELECT CAST(event_index AS BIGINT) AS event_index, chrom, "
"CASE WHEN try_cast(pos AS DOUBLE) = floor(try_cast(pos AS DOUBLE)) THEN try_cast(pos AS BIGINT) END AS pos, "
"id, ref, alt, info FROM ";

static const char pairs_body_0[] =
"), parsed AS (SELECT *, TRY(duckvep_breakend_geometry(alt)) AS g, "
"string_split(coalesce(info, ''), ';') AS tokens, "
"CASE WHEN id IS NULL OR id IN ('', '.') THEN NULL ELSE id END AS rid FROM source), "
"tokened AS (SELECT *, "
"list_filter(tokens, lambda x: split_part(x, '=', 1) = 'MATEID') AS mate_tokens, "
"list_filter(tokens, lambda x: split_part(x, '=', 1) = 'EVENT') AS event_tokens "
"FROM parsed), "
"records AS (SELECT event_index, CAST(chrom AS VARCHAR) AS chrom, pos, rid AS id, g, "
"CASE WHEN g IS NOT NULL THEN CASE WHEN g.mate_chrom IS NULL THEN 'single_breakend' "
"ELSE 'paired_breakend' END "
"WHEN alt IS NOT NULL AND (regexp_matches(alt, '[\\[\\]]') OR "
"(length(alt) > 1 AND (starts_with(alt, '.') OR ends_with(alt, '.')))) THEN 'malformed_alt' "
"ELSE 'not_breakend' END AS record_kind, "
"length(mate_tokens) AS mate_token_count, length(event_tokens) AS event_token_count, "
"CASE WHEN length(mate_tokens) = 1 THEN nullif(substring(mate_tokens[1], 8), '') END AS mate_id, "
"CASE WHEN length(event_tokens) = 1 THEN nullif(substring(event_tokens[1], 7), '') END AS event_id, "
"g.mate_chrom AS declared_mate_chrom, "
"CAST(g.mate_position AS BIGINT) AS declared_mate_position, "
"g.local_join_after AS local_join_after, g.mate_extends_right AS mate_extends_right, "
"CASE WHEN g IS NOT NULL THEN greatest(length(g.replacement_sequence) - 1, 0) END AS inserted_length, "
"CASE WHEN g IS NOT NULL AND length(g.replacement_sequence) > 1 THEN "
"CASE WHEN g.local_join_after THEN substring(g.replacement_sequence, 2) "
"ELSE substring(g.replacement_sequence, 1, length(g.replacement_sequence) - 1) END END AS inserted_sequence "
"FROM tokened), "
"ids AS (SELECT id, count(*) AS n, min(event_index) AS first_index FROM records "
"WHERE id IS NOT NULL GROUP BY id), "
"joined AS (SELECT r.*, own.n AS id_count, mi.n AS mate_id_count, m.event_index AS m_index, "
"m.chrom AS m_chrom, m.pos AS m_pos, m.mate_id AS m_mate_id, m.event_id AS m_event_id, "
"m.declared_mate_chrom AS m_dchrom, m.declared_mate_position AS m_dpos, "
"m.local_join_after AS m_join, m.mate_extends_right AS m_right, m.inserted_length AS m_ins, "
"m.record_kind AS m_kind, m.mate_token_count AS m_mate_count, m.event_token_count AS m_event_count FROM records r "
"LEFT JOIN ids own ON own.id = r.id LEFT JOIN ids mi ON mi.id = r.mate_id "
"LEFT JOIN records m ON m.event_index = CASE WHEN mi.n = 1 THEN mi.first_index END), "
"checks AS (SELECT *, "
"coalesce(m_kind = 'paired_breakend' AND m_mate_id = id AND m_index <> event_index, false) AS id_reciprocal, "
"coalesce(m_kind = 'paired_breakend' AND declared_mate_chrom = m_chrom AND "
"declared_mate_position = m_pos AND m_dchrom = chrom AND m_dpos = pos, false) AS coordinate_reciprocal, "
"coalesce(m_kind = 'paired_breakend' AND local_join_after = NOT m_right AND "
"mate_extends_right = NOT m_join, false) AS orientation_reciprocal, "
"coalesce(m_kind = 'paired_breakend' AND inserted_length = m_ins, false) AS insert_agree, "
"event_id IS NOT DISTINCT FROM m_event_id AS event_agree, "
"CASE WHEN event_id IS NULL THEN NULL ELSE count(*) OVER (PARTITION BY event_id) END AS event_record_count "
"FROM joined), "
"decided AS (SELECT *, CASE record_kind "
"WHEN 'not_breakend' THEN 'not_breakend' "
"WHEN 'malformed_alt' THEN 'malformed_alt' "
"WHEN 'single_breakend' THEN CASE WHEN mate_token_count > 0 THEN 'single_breakend_with_mateid' "
"ELSE 'single_breakend' END "
"ELSE CASE "
"WHEN id IS NULL THEN 'missing_id' "
"WHEN id_count > 1 THEN 'duplicate_id' ";

static const char pairs_body_1[] =
"WHEN mate_token_count > 1 OR contains(coalesce(mate_id, ''), ',') THEN 'duplicate_mateid' "
"WHEN mate_id IS NULL THEN 'missing_mateid' "
"WHEN mate_id = id THEN 'self_mate' "
"WHEN event_token_count > 1 OR coalesce(m_event_count, 0) > 1 THEN 'duplicate_event' "
"WHEN mate_id_count IS NULL THEN 'mate_not_found' "
"WHEN mate_id_count > 1 THEN 'ambiguous_mate_id' "
"WHEN m_kind <> 'paired_breakend' THEN 'mate_not_paired_breakend' "
"WHEN m_mate_count > 1 OR contains(coalesce(m_mate_id, ''), ',') THEN 'mate_duplicate_mateid' "
"WHEN m_mate_id IS NULL THEN 'mate_missing_mateid' "
"WHEN NOT id_reciprocal THEN 'mate_not_reciprocal' "
"WHEN NOT coordinate_reciprocal THEN 'mate_coordinate_conflict' "
"WHEN NOT orientation_reciprocal THEN 'mate_orientation_conflict' "
"WHEN NOT insert_agree THEN 'mate_insert_conflict' "
"WHEN NOT event_agree THEN 'event_conflict' "
"ELSE 'reciprocal' END END AS reason FROM checks) "
"SELECT event_index, id, chrom, pos, record_kind, "
"CASE WHEN reason = 'reciprocal' THEN 'reciprocal' "
"WHEN reason IN ('not_breakend', 'single_breakend') THEN 'not_applicable' "
"WHEN reason = 'malformed_alt' THEN 'invalid' "
"WHEN reason IN ('missing_id', 'missing_mateid', 'mate_not_found', 'mate_missing_mateid') THEN 'unproven' "
"ELSE 'conflict' END AS status, reason, mate_id, event_id, "
"CASE WHEN record_kind = 'paired_breakend' AND mate_id_count = 1 AND mate_id <> id THEN m_index END AS mate_event_index, "
"declared_mate_chrom, declared_mate_position, "
"CASE WHEN record_kind = 'paired_breakend' THEN id_reciprocal END AS id_reciprocal, "
"CASE WHEN record_kind = 'paired_breakend' THEN coordinate_reciprocal END AS coordinate_reciprocal, "
"CASE WHEN record_kind = 'paired_breakend' THEN orientation_reciprocal END AS orientation_reciprocal, "
"CASE WHEN record_kind = 'paired_breakend' THEN insert_agree END AS insert_agree, "
"CASE WHEN record_kind = 'paired_breakend' THEN event_agree END AS event_agree, "
"event_record_count, inserted_length, inserted_sequence, "
"CASE WHEN reason = 'reciprocal' THEN CAST(least(event_index, m_index) AS VARCHAR) || '|' || "
"CAST(greatest(event_index, m_index) AS VARCHAR) END AS pair_key, "
"'not_evaluated' AS phase_status, 'not_asserted' AS fusion_status FROM decided ORDER BY event_index";

static const piece pairs_pieces[] = {
    {PIECE_TEXT, pairs_head},
    {PIECE_RELATION, "0"},
    {PIECE_TEXT, pairs_body_0},
    {PIECE_TEXT, pairs_body_1},
};
static const builder_spec pairs_spec = {"duckvep_prepare_breakend_pairs_sql", 1, pairs_pieces,
    sizeof(pairs_pieces) / sizeof(pairs_pieces[0]), false};
static void prepare_pairs(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    run_builder(&pairs_spec, info, input, output);
}

/* ---- Breakend gene evidence ---------------------------------------------- */

static const char fusion_head[] =
"WITH pairs AS (SELECT event_index, mate_event_index, pair_key, status, reason, "
"id_reciprocal, coordinate_reciprocal, orientation_reciprocal, event_agree FROM ";

static const char fusion_mid[] =
"), gene_sets AS (SELECT event_index, "
"list_sort(list_distinct(list(gene_id) FILTER (WHERE gene_id IS NOT NULL))) AS gene_ids FROM ";

static const char fusion_tail[] =
" GROUP BY event_index), "
"joined AS (SELECT p.*, coalesce(g.gene_ids, CAST([] AS VARCHAR[])) AS endpoint_genes, "
"coalesce(mg.gene_ids, CAST([] AS VARCHAR[])) AS mate_endpoint_genes "
"FROM pairs p LEFT JOIN gene_sets g ON g.event_index = p.event_index "
"LEFT JOIN gene_sets mg ON mg.event_index = p.mate_event_index) "
"SELECT event_index, mate_event_index, pair_key, "
"CASE WHEN NOT coalesce(id_reciprocal AND coordinate_reciprocal AND event_agree, false) "
"THEN 'identity_unproven' "
"WHEN length(endpoint_genes) = 0 OR length(mate_endpoint_genes) = 0 THEN 'endpoint_without_gene' "
"WHEN list_has_any(endpoint_genes, mate_endpoint_genes) THEN 'shared_gene_endpoints' "
"WHEN orientation_reciprocal THEN 'candidate_partner_genes' "
"ELSE 'candidate_orientation_conflict' END AS status, "
"CASE WHEN NOT coalesce(id_reciprocal AND coordinate_reciprocal AND event_agree, false) "
"THEN reason END AS reason, "
"endpoint_genes, mate_endpoint_genes, "
"false AS fusion_asserted, 'unproven' AS phase_status FROM joined ORDER BY event_index";

static const piece fusion_pieces[] = {
    {PIECE_TEXT, fusion_head},
    {PIECE_RELATION, "0"},
    {PIECE_TEXT, fusion_mid},
    {PIECE_RELATION, "1"},
    {PIECE_TEXT, fusion_tail},
};
static const builder_spec fusion_spec = {"duckvep_prepare_breakend_fusion_sql", 2, fusion_pieces,
    sizeof(fusion_pieces) / sizeof(fusion_pieces[0]), false};
static void prepare_fusion(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    run_builder(&fusion_spec, info, input, output);
}

/* ---- Structural HGVS ------------------------------------------------------- */

static const char hgvs_head[] =
"WITH source AS (SELECT CAST(event_index AS BIGINT) AS event_index, chrom, pos AS raw_pos, "
"CASE WHEN try_cast(pos AS DOUBLE) = floor(try_cast(pos AS DOUBLE)) THEN try_cast(pos AS BIGINT) END AS pos, "
"ref, alt, info FROM ";

static const char hgvs_mid[] =
"), refs AS (SELECT CAST(event_index AS BIGINT) AS event_index, count(*) AS n, min(upper(reference_sequence)) AS refseq FROM ";

static const char hgvs_a[] =
" GROUP BY 1), "
"parsed AS (SELECT s.*, r.n AS ref_rows, r.refseq, string_split(coalesce(s.info, ''), ';') AS tokens, "
"TRY(duckvep_breakend_geometry(s.alt)) AS g FROM source s LEFT JOIN refs r USING (event_index)), "
"fields AS (SELECT *, "
"list_filter(tokens, lambda x: split_part(x, '=', 1) = 'END') AS end_tokens, "
"list_filter(tokens, lambda x: split_part(x, '=', 1) IN ('CIPOS', 'CIEND', 'IMPRECISE')) AS imprecise_tokens "
"FROM parsed), "
"spans AS (SELECT *, "
"CASE WHEN length(end_tokens) = 1 THEN nullif(substring(end_tokens[1], 5), '') END AS end_text, "
"CASE WHEN alt IN ('<DEL>', '<DUP>', '<DUP:TANDEM>', '<INV>') THEN alt END AS symbolic_kind FROM fields), "
"bounds AS (SELECT *, "
"CASE WHEN regexp_full_match(end_text, '[1-9][0-9]{0,15}') THEN CAST(end_text AS BIGINT) END AS end_pos "
"FROM spans), "
"sequences AS (SELECT *, "
"CASE WHEN refseq IS NOT NULL AND end_pos > pos THEN substring(refseq, 1, end_pos - pos + 1) END AS core, "
"CASE WHEN refseq IS NOT NULL AND end_pos > pos THEN substring(refseq, 2, end_pos - pos) END AS span, "
"CASE WHEN refseq IS NOT NULL AND end_pos > pos AND length(refseq) = 2 * (end_pos - pos) + 1 THEN "
"substring(refseq, end_pos - pos + 2, end_pos - pos) END AS flank FROM bounds), "
"inverted AS (SELECT *, reverse(translate(span, 'ACGT', 'TGCA')) AS inverse FROM sequences), "
"mismatches AS (SELECT *, "
"CASE WHEN span IS NOT NULL AND regexp_full_match(span, '[ACGT]+') THEN "
"list_filter(range(1, length(span) + 1), lambda i: span[i] != inverse[i]) END AS differing "
"FROM inverted), "
"cores AS (SELECT *, "
"CASE WHEN differing IS NOT NULL AND length(differing) > 0 THEN differing[1] - 1 END AS lead_equal, "
"CASE WHEN differing IS NOT NULL AND length(differing) > 0 THEN length(span) - differing[-1] END AS trail_equal "
"FROM mismatches), "
"inversions AS (SELECT *, "
"CASE WHEN lead_equal IS NOT NULL THEN substring(span, lead_equal + 1, length(span) - lead_equal - trail_equal) END AS core_ref, "
"CASE WHEN lead_equal IS NOT NULL THEN substring(inverse, lead_equal + 1, length(span) - lead_equal - trail_equal) END AS core_alt "
"FROM cores), "
"decided AS (SELECT *, CASE "
"WHEN alt IS NULL OR raw_pos IS NULL OR ref IS NULL THEN 'unavailable:missing_field' "
"WHEN pos IS NULL THEN 'unavailable:invalid_pos' "
"WHEN g IS NOT NULL THEN 'unsupported:breakend' "
"WHEN regexp_full_match(alt, '<INS(:[^>]*)?>') THEN 'unsupported:symbolic_insertion' "
"WHEN regexp_full_match(alt, '<(CNV|CN[0-9]+)>') THEN 'unsupported:copy_number' "
"WHEN regexp_full_match(alt, '<STR[0-9]*>') THEN 'unsupported:repeat_expansion' "
"WHEN regexp_full_match(alt, '<[^>]+>') AND symbolic_kind IS NULL THEN 'unsupported:symbolic_allele' "
"WHEN symbolic_kind IS NULL AND regexp_full_match(alt, '[ACGTN]+') THEN 'unsupported:literal_allele' "
"WHEN symbolic_kind IS NULL THEN 'unsupported:alt_syntax' "
"WHEN chrom IS NULL OR chrom = '' THEN 'unavailable:missing_chrom' "
"WHEN length(imprecise_tokens) > 0 THEN 'unsupported:imprecise' "
"WHEN length(end_tokens) = 0 THEN 'unavailable:missing_end' "
"WHEN length(end_tokens) > 1 THEN 'unavailable:duplicate_end' "
"WHEN end_pos IS NULL OR pos < 1 OR end_pos <= pos THEN 'unavailable:invalid_end' "
"WHEN end_pos - pos > ";

static const char hgvs_b[] =
" THEN 'unsupported:span_capacity' "
"WHEN ref_rows IS NULL OR refseq IS NULL THEN 'unavailable:missing_reference_sequence' "
"WHEN ref_rows > 1 THEN 'unavailable:ambiguous_reference_sequence' "
"WHEN length(refseq) NOT IN (end_pos - pos + 1, 2 * (end_pos - pos) + 1) THEN 'unavailable:reference_length_mismatch' "
"WHEN NOT regexp_full_match(refseq, '[ACGT]+') THEN 'unavailable:reference_alphabet' "
"WHEN NOT regexp_full_match(upper(ref), '[ACGTN]') THEN 'unavailable:invalid_ref' "
"WHEN upper(ref) != 'N' AND upper(ref) != substring(refseq, 1, 1) THEN 'unavailable:reference_anchor_mismatch' "
"WHEN symbolic_kind IN ('<DUP>', '<DUP:TANDEM>') AND flank IS NULL THEN 'unavailable:missing_flank_sequence' "
"WHEN symbolic_kind IN ('<DUP>', '<DUP:TANDEM>') AND flank = span THEN 'unsupported:duplication_adjacent_repeat' "
"WHEN symbolic_kind = '<INV>' AND lead_equal IS NULL THEN 'unsupported:inversion_identity' "
"WHEN symbolic_kind = '<INV>' AND (length(core_ref) < 2 OR "
"core_alt != reverse(translate(core_ref, 'ACGT', 'TGCA'))) THEN 'unsupported:inversion_not_reducible' "
"ELSE 'supported:' || CASE symbolic_kind WHEN '<DEL>' THEN 'deletion' WHEN '<INV>' THEN 'inversion' "
"ELSE 'duplication' END END AS decision FROM inversions) "
"SELECT event_index, split_part(decision, ':', 1) AS hgvs_status, "
"CASE WHEN decision LIKE 'supported:%' THEN NULL ELSE split_part(decision, ':', 2) END AS hgvs_reason, "
"CASE WHEN decision LIKE 'supported:%' THEN split_part(decision, ':', 2) END AS edit, "
"CASE WHEN decision LIKE 'supported:%' THEN CAST(chrom AS VARCHAR) || ':g.' || CAST(pos + 1 AS VARCHAR) || "
"CASE WHEN end_pos > pos + 1 THEN '_' || CAST(end_pos AS VARCHAR) ELSE '' END || "
"CASE decision WHEN 'supported:deletion' THEN 'del' WHEN 'supported:inversion' THEN 'inv' "
"ELSE 'dup' END END AS hgvs_g, "
"CASE WHEN decision LIKE 'supported:%' THEN 'none' END AS normalization, "
"CASE split_part(decision, ':', 1) WHEN 'supported' THEN 'literal_equivalent' "
"ELSE split_part(decision, ':', 1) END AS transcript_hgvs_route, "
"CASE decision WHEN 'supported:deletion' THEN CAST(pos AS BIGINT) "
"WHEN 'supported:duplication' THEN CAST(end_pos AS BIGINT) "
"WHEN 'supported:inversion' THEN CAST(pos + 1 AS BIGINT) END AS literal_position, "
"CASE decision WHEN 'supported:deletion' THEN core "
"WHEN 'supported:duplication' THEN substring(core, length(core), 1) "
"WHEN 'supported:inversion' THEN span END AS literal_reference, "
"CASE decision WHEN 'supported:deletion' THEN substring(core, 1, 1) "
"WHEN 'supported:duplication' THEN substring(core, length(core), 1) || span "
"WHEN 'supported:inversion' THEN inverse END AS literal_alternate FROM decided ORDER BY event_index";

static const piece hgvs_pieces[] = {
    {PIECE_TEXT, hgvs_head},
    {PIECE_RELATION, "0"},
    {PIECE_TEXT, hgvs_mid},
    {PIECE_RELATION, "1"},
    {PIECE_TEXT, hgvs_a},
    {PIECE_INTEGER, NULL},
    {PIECE_TEXT, hgvs_b},
};
static const builder_spec hgvs_spec = {"duckvep_prepare_structural_hgvs_sql", 2, hgvs_pieces,
    sizeof(hgvs_pieces) / sizeof(hgvs_pieces[0]), true};
static void prepare_hgvs(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    run_builder(&hgvs_spec, info, input, output);
}

bool register_duckvep_structural_sql(duckdb_connection connection) {
    return duckvep_register_builder(connection, "duckvep_prepare_breakend_pairs_sql", 1, prepare_pairs) &&
        duckvep_register_builder(connection, "duckvep_prepare_breakend_fusion_sql", 2, prepare_fusion) &&
        duckvep_register_builder(connection, "duckvep_prepare_structural_hgvs_sql", 2, prepare_hgvs);
}
