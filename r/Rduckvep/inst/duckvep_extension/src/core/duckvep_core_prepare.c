#include "duckvep_core_prepare.h"

#include <stdio.h>
#include <string.h>

const char duckvep_core_prepare_sv_failed[] =
    "duckvep_prepare_sv_geometry_sql: invalid relation name or allocation failure";
const char duckvep_core_prepare_expansionhunter_failed[] =
    "duckvep_prepare_expansionhunter_sql: invalid relation name or allocation failure";

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

bool duckvep_core_prepare_sv_sql(const char *relation, duckvep_sql_text *out)
{
    return relation && duckvep_sql_append(out, "WITH source AS (SELECT event_index, pos, ref, alt, info FROM ") &&
        duckvep_sql_relation(out, relation) && duckvep_sql_append(out, sv_sql);
}

bool duckvep_core_prepare_expansionhunter_sql(const char *events, const char *reference,
    duckvep_sql_text *out)
{
    return duckvep_sql_append(out, "WITH source AS (SELECT event_index, info, format, \"sample\", ref, alt, alt_index FROM ") &&
        duckvep_sql_relation(out, events) && duckvep_sql_append(out, str_sql) &&
        duckvep_sql_relation(out, reference) && duckvep_sql_append(out, str_tail) &&
        duckvep_sql_append(out, str_tail2);
}

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

typedef struct {
    const char *name;
    size_t relations;
    const piece *pieces;
    size_t piece_count;
    bool has_max_span;
} builder_spec;

static const builder_spec specs[] = {
    {"duckvep_prepare_breakend_pairs_sql", 1, pairs_pieces, sizeof(pairs_pieces) / sizeof(pairs_pieces[0]), false},
    {"duckvep_prepare_breakend_fusion_sql", 2, fusion_pieces, sizeof(fusion_pieces) / sizeof(fusion_pieces[0]), false},
    {"duckvep_prepare_structural_hgvs_sql", 2, hgvs_pieces, sizeof(hgvs_pieces) / sizeof(hgvs_pieces[0]), true},
};

size_t duckvep_core_structural_relations(duckvep_structural_kind_t kind) { return specs[kind].relations; }
bool duckvep_core_structural_has_max_span(duckvep_structural_kind_t kind) { return specs[kind].has_max_span; }
const char *duckvep_core_structural_name(duckvep_structural_kind_t kind) { return specs[kind].name; }

bool duckvep_core_structural_max_span(const duckvep_cell_t *cell, int64_t *max_span,
    char *message, size_t message_size)
{
    int64_t parsed = 0;
    bool overflow = false;
    *max_span = 5000;
    if (!cell || !cell->valid) return true;
    switch (cell->kind) {
    case DUCKVEP_CELL_TINYINT: case DUCKVEP_CELL_SMALLINT:
    case DUCKVEP_CELL_INTEGER: case DUCKVEP_CELL_BIGINT:
        parsed = cell->i; break;
    case DUCKVEP_CELL_UTINYINT: case DUCKVEP_CELL_USMALLINT: case DUCKVEP_CELL_UINTEGER:
        parsed = (int64_t)cell->u; break;
    case DUCKVEP_CELL_UBIGINT:
        if (cell->u > (uint64_t)INT64_MAX) overflow = true; else parsed = (int64_t)cell->u;
        break;
    default: overflow = true; break;
    }
    if (overflow || parsed < 1 || parsed > 60000) {
        snprintf(message, message_size, "DuckVEP builder: option '%s' must be between %lld and %lld",
                 "max_span", (long long)1, (long long)60000);
        return false;
    }
    *max_span = parsed;
    return true;
}

bool duckvep_core_structural_sql(duckvep_structural_kind_t kind, char *const relations[2],
    int64_t max_span, duckvep_sql_text *out)
{
    return assemble(out, specs[kind].pieces, specs[kind].piece_count, relations, max_span);
}
