#include "duckvep_core_prepare.h"

#include <stdio.h>
#include <string.h>

/* ---- Reading read_bcf-shaped columns -------------------------------------
 * The VCF header decides an INFO or FORMAT column's shape: read_bcf gives a list
 * for Number=A, R or ., a scalar for Number=1, and VARCHAR for Type=String.
 * DUCKVEP_VCF_VALUES reads any of these as VARCHAR[] without an extension: a list
 * column's elements, or a scalar's text split at the commas VCF writes between list
 * values; NULL stays NULL. TRY_CAST of a scalar to a list is NULL at run time, not
 * a bind error (measured on DuckDB 1.5.1, 1.5.6 and 2.0.0), so one SQL text serves
 * every shape. A scalar string written like a bracketed list reads as that list. A
 * whole number written with a zero fraction, as a DOUBLE column renders 40 ("40.0"),
 * reads as that whole number.
 * Optional columns come from a template of untyped NULLs joined UNION ALL BY NAME,
 * which keeps each present column's own type. */
#define DUCKVEP_VCF_VALUES(x) \
    "list_transform(CASE WHEN try_cast(" x " AS VARCHAR[]) IS NOT NULL THEN try_cast(" x " AS VARCHAR[]) " \
    "ELSE string_split(try_cast(" x " AS VARCHAR), ',') END, " \
    "lambda v: regexp_replace(v, '^([+-]?[0-9]+)\\.0+$', '\\1'))"
/* The first value, trimmed, or NULL when empty. VEP reads the first comma value. */
#define DUCKVEP_VCF_FIRST(x) "nullif(trim(" DUCKVEP_VCF_VALUES(x) "[1]), '')"
/* An integer read from such text: only whole-number text counts, because TRY_CAST
 * rounds '20.5' to 21 and reads '1e3' as 1000. */
#define DUCKVEP_VCF_INTEGER(text) \
    "CASE WHEN regexp_full_match(" text ", '[+-]?[0-9]{1,18}') THEN try_cast(" text " AS BIGINT) END"
/* The allele a row describes: element alt_index (1-based, default 1) of ALT. */
#define DUCKVEP_VCF_ALT \
    DUCKVEP_VCF_VALUES("\"ALT\"") "[coalesce(try_cast(alt_index AS BIGINT), 1)]"
/* VEP 116 ensembl-io BaseVCF4 get_start/get_end for a structural record: it starts
 * at POS + 1 and ends at POS + |SVLEN| when SVLEN is set, else at END, else at its
 * start. The arguments are BIGINT expressions. */
#define DUCKVEP_VEP_SV_START(pos) "try_cast(" pos "::HUGEINT + 1 AS BIGINT)"
#define DUCKVEP_VEP_SV_END(pos, svlen, end) \
    "CASE WHEN " svlen " IS NOT NULL THEN try_cast(" pos "::HUGEINT + abs(" svlen "::HUGEINT) AS BIGINT) " \
    "WHEN " end " IS NOT NULL THEN " end " ELSE try_cast(" pos "::HUGEINT + 1 AS BIGINT) END"

const char duckvep_core_prepare_sv_failed[] =
    "duckvep_prepare_sv_geometry_sql: invalid relation name or allocation failure";
const char duckvep_core_prepare_expansionhunter_failed[] =
    "duckvep_prepare_expansionhunter_sql: invalid relation name or allocation failure";
const char duckvep_core_vcf_events_failed[] =
    "duckvep_vcf_events_sql: invalid relation name or allocation failure";

/* Records are read_bcf-shaped: event_index, POS, REF and ALT (a list, with the
 * 1-based alt_index choosing the allele, or one allele), and the optional INFO_END,
 * INFO_SVLEN, INFO_CIPOS, INFO_CIEND and INFO_SEQ. A symbolic allele's nominal span
 * follows VEP's get_start/get_end, the same rule duckvep_vcf_events_sql applies. */
static const char sv_head[] =
"WITH __duckvep_sv_source AS (FROM (SELECT NULL AS alt_index, NULL AS \"INFO_END\", NULL AS \"INFO_SVLEN\", "
"NULL AS \"INFO_CIPOS\", NULL AS \"INFO_CIEND\", NULL AS \"INFO_SEQ\" WHERE false) UNION ALL BY NAME "
"SELECT COLUMNS('(?i)^(event_index|pos|ref|alt|alt_index|info_end|info_svlen|info_cipos|info_ciend|info_seq)$') FROM ";

static const char sv_sql[] =
"), __duckvep_sv_rows AS (SELECT event_index, pos, ref, alt, alt_index, \"INFO_END\", \"INFO_SVLEN\", "
"\"INFO_CIPOS\", \"INFO_CIEND\", \"INFO_SEQ\" FROM __duckvep_sv_source), "
"source AS (SELECT event_index, try_cast(pos AS BIGINT) AS pos, CAST(ref AS VARCHAR) AS ref, "
DUCKVEP_VCF_ALT " AS alt, "
DUCKVEP_VCF_FIRST("\"INFO_END\"") " AS end_text, " DUCKVEP_VCF_FIRST("\"INFO_SVLEN\"") " AS svlen_text, "
DUCKVEP_VCF_FIRST("\"INFO_SEQ\"") " AS source_sequence, "
DUCKVEP_VCF_VALUES("\"INFO_CIPOS\"") " AS cipos, " DUCKVEP_VCF_VALUES("\"INFO_CIEND\"") " AS ciend "
"FROM __duckvep_sv_rows), parts AS (SELECT *, "
"regexp_full_match(alt, '<[^>]+>') AS symbolic, "
"regexp_full_match(ref, '[ACGT]+') AND regexp_full_match(alt, '[ACGT]+') "
"AND length(alt)>length(ref) AND starts_with(alt,ref) AS literal, "
DUCKVEP_VCF_INTEGER("end_text") " AS info_end, " DUCKVEP_VCF_INTEGER("svlen_text") " AS info_svlen FROM source), "
"";

/* ISO C99 caps a string literal at 4,095 characters, so the text continues here. */
static const char sv_sql2[] =
"bounds AS (SELECT *, CASE WHEN symbolic THEN " DUCKVEP_VEP_SV_START("pos") " "
"ELSE try_cast(pos::HUGEINT + length(ref) AS BIGINT) END AS nominal_start, "
"CASE WHEN symbolic THEN " DUCKVEP_VEP_SV_END("pos", "info_svlen", "info_end") " "
"ELSE try_cast(pos::HUGEINT + length(ref) - 1 AS BIGINT) END AS nominal_end, "
"CASE WHEN len(cipos) = 2 THEN " DUCKVEP_VCF_INTEGER("trim(cipos[1])") " END AS cipos_lo, "
"CASE WHEN len(cipos) = 2 THEN " DUCKVEP_VCF_INTEGER("trim(cipos[2])") " END AS cipos_hi, "
"CASE WHEN len(ciend) = 2 THEN " DUCKVEP_VCF_INTEGER("trim(ciend[1])") " END AS ciend_lo, "
"CASE WHEN len(ciend) = 2 THEN " DUCKVEP_VCF_INTEGER("trim(ciend[2])") " END AS ciend_hi FROM parts), "
"geometry AS (SELECT *, try_cast(nominal_start::HUGEINT+cipos_lo AS BIGINT) AS outer_start, "
"try_cast(nominal_start::HUGEINT+cipos_hi AS BIGINT) AS inner_start, "
"try_cast(nominal_end::HUGEINT+ciend_lo AS BIGINT) AS inner_end, "
"try_cast(nominal_end::HUGEINT+ciend_hi AS BIGINT) AS outer_end FROM bounds) "
"SELECT event_index, CASE WHEN "
"(symbolic OR literal) AND pos BETWEEN 1 AND 2147483647 AND nominal_start BETWEEN 0 AND 2147483647 "
"AND nominal_end BETWEEN nominal_start::HUGEINT-1 AND 2147483647 "
"AND (end_text IS NULL OR info_end > 0) AND (svlen_text IS NULL OR info_svlen IS NOT NULL) "
"AND (NOT symbolic OR (regexp_full_match(ref,'[ACGTN]') AND (info_end IS NOT NULL OR info_svlen IS NOT NULL))) "
"AND (cipos IS NULL OR (outer_start BETWEEN 0 AND 2147483647 AND inner_start BETWEEN 0 AND 2147483647 AND outer_start<=inner_start)) "
"AND (ciend IS NULL OR (inner_end BETWEEN 0 AND 2147483647 AND outer_end BETWEEN 0 AND 2147483647 AND inner_end<=outer_end)) "
"THEN 'ok' ELSE 'unsupported_geometry' END AS status, "
"CASE WHEN symbolic THEN 'structural' WHEN literal THEN 'literal_insertion' ELSE 'unsupported' END AS mode, "
"nominal_start, nominal_end, outer_start, inner_start, inner_end, outer_end, "
"CASE WHEN literal THEN substring(alt,length(ref)+1) END AS inserted_sequence, source_sequence FROM geometry";

/* Rows are read_bcf-shaped, one per event and sample as read_bcf(..., tidy_format :=
 * true) gives them: event_index, REF, ALT (a list or one comma-separated text) and the
 * 1-based alt_index, the optional INFO_END, INFO_REF, INFO_RL and INFO_RU, and the
 * optional FORMAT_GT, FORMAT_SO and either FORMAT_REPCN with FORMAT_REPCI or FORMAT_CN
 * with FORMAT_CI. The '/'-separated values inside the FORMAT fields are
 * ExpansionHunter's own notation and are read here. */
static const char str_head[] =
"WITH __duckvep_str_source AS (FROM (SELECT NULL AS \"INFO_END\", NULL AS \"INFO_REF\", NULL AS \"INFO_RL\", "
"NULL AS \"INFO_RU\", NULL AS \"FORMAT_GT\", NULL AS \"FORMAT_SO\", NULL AS \"FORMAT_REPCN\", "
"NULL AS \"FORMAT_REPCI\", NULL AS \"FORMAT_CN\", NULL AS \"FORMAT_CI\" WHERE false) UNION ALL BY NAME "
"SELECT COLUMNS('(?i)^(event_index|ref|alt|alt_index|info_end|info_ref|info_rl|info_ru|"
"format_gt|format_so|format_repcn|format_repci|format_cn|format_ci)$') FROM ";

static const char str_source[] =
"), __duckvep_str_rows AS (SELECT event_index, ref, alt, alt_index, \"INFO_END\", \"INFO_REF\", \"INFO_RL\", \"INFO_RU\", "
"\"FORMAT_GT\", \"FORMAT_SO\", \"FORMAT_REPCN\", \"FORMAT_REPCI\", \"FORMAT_CN\", \"FORMAT_CI\" FROM __duckvep_str_source), "
"source AS (SELECT event_index, CAST(ref AS VARCHAR) AS ref, " DUCKVEP_VCF_VALUES("\"ALT\"") " AS alts, alt_index, "
DUCKVEP_VCF_FIRST("\"INFO_REF\"") " AS ref_text, " DUCKVEP_VCF_FIRST("\"INFO_RL\"") " AS rl_text, "
DUCKVEP_VCF_FIRST("\"INFO_RU\"") " AS unit, " DUCKVEP_VCF_FIRST("\"INFO_END\"") " AS end_text, "
DUCKVEP_VCF_FIRST("\"FORMAT_GT\"") " AS gt_text, " DUCKVEP_VCF_FIRST("\"FORMAT_SO\"") " AS so_text, "
DUCKVEP_VCF_FIRST("\"FORMAT_REPCN\"") " AS repcn_text, " DUCKVEP_VCF_FIRST("\"FORMAT_REPCI\"") " AS repci_text, "
DUCKVEP_VCF_FIRST("\"FORMAT_CN\"") " AS cn_text, " DUCKVEP_VCF_FIRST("\"FORMAT_CI\"") " AS ci_text "
"FROM __duckvep_str_rows";

static const char str_sql[] =
"), joined AS (SELECT s.*, r.reference_sequence, r.ref_rows FROM source s LEFT JOIN "
"(SELECT event_index, count(*) AS ref_rows, CASE WHEN count(*) > 1 THEN NULL ELSE min(reference_sequence) END AS reference_sequence FROM " ;

static const char str_tail[] =
" GROUP BY event_index) r USING (event_index)), metadata AS (SELECT *, "
"list_extract(alts, try_cast(alt_index AS INTEGER)) selected, "
"cn_text IS NOT NULL AND ci_text IS NOT NULL has_cn, repcn_text IS NOT NULL AND repci_text IS NOT NULL has_repcn "
"FROM joined), counts AS (SELECT *, "
"try_cast(ref_text AS DOUBLE) ref_count, try_cast(rl_text AS DOUBLE) ref_length, "
"try_cast(regexp_extract(selected,'^<STR([0-9]+)>$',1) AS DOUBLE) alt_count, "
"CASE WHEN has_cn THEN cn_text ELSE repcn_text END count_text, "
"CASE WHEN has_cn THEN ci_text ELSE repci_text END ci_value "
"FROM metadata), samples AS (SELECT *, "
"string_split(replace(coalesce(gt_text,''),'|','/'),'/') gt, "
"string_split(coalesce(count_text,''),'/') copies, "
"string_split(coalesce(ci_value,''),'/') ranges, "
"string_split(coalesce(so_text,''),'/') support "
"FROM counts), calls AS (SELECT *, "
"list_filter(range(1,length(gt)+1),lambda i: gt[i]=cast(try_cast(alt_index AS BIGINT) AS VARCHAR)) selected_calls "
"FROM samples), ";

static const char str_tail2[] =
"decisions AS (SELECT *, CASE "
"WHEN ref_rows > 1 THEN 'ambiguous_reference_sequence' "
"WHEN ref IS NULL OR alts IS NULL OR reference_sequence IS NULL THEN 'missing_field' "
"WHEN alt_index IS NULL OR NOT isfinite(try_cast(alt_index AS DOUBLE)) OR "
"try_cast(alt_index AS DOUBLE) != floor(try_cast(alt_index AS DOUBLE)) OR "
"try_cast(alt_index AS DOUBLE) NOT BETWEEN 1 AND length(alts) THEN 'alt_index' "
"WHEN NOT regexp_full_match(ref,'[ACGT]') OR NOT regexp_full_match(reference_sequence,'[ACGT]+') THEN 'literal_reference' "
"WHEN ref_text IS NULL OR rl_text IS NULL OR unit IS NULL OR end_text IS NULL THEN 'missing_info' "
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
"WHEN gt_text IS NULL OR so_text IS NULL OR NOT (has_cn OR has_repcn) THEN 'missing_format' "
"WHEN has_cn AND has_repcn THEN 'ambiguous_count_fields' "
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
"struct_pack(ref := ref, alts := alts, alt_index := alt_index, \"end\" := end_text, reference_copies := ref_text, "
"reference_length := rl_text, repeat_unit := unit, genotype := gt_text, read_support := so_text, "
"copy_numbers := count_text, copy_number_intervals := ci_value, "
"reference_sequence := reference_sequence) AS \"source\", "
"CASE WHEN reason NOT IN ('ambiguous_reference_sequence','missing_field','alt_index','literal_reference',"
"'missing_info','metadata_syntax','reference_mismatch','allele_capacity') "
"OR (reason='allele_capacity' AND ref_length<=5000) "
"THEN struct_pack(unit := unit, count := ref_count) END reference_components, "
"CASE WHEN reason='exact' THEN struct_pack(unit := unit, count := alt_count) END alternate_components, "
"reason='exact' sequence_exact FROM decisions";

bool duckvep_core_prepare_sv_sql(const char *relation, duckvep_sql_text *out)
{
    return relation && duckvep_sql_append(out, sv_head) &&
        duckvep_sql_relation(out, relation) && duckvep_sql_append(out, sv_sql) &&
        duckvep_sql_append(out, sv_sql2);
}

bool duckvep_core_prepare_expansionhunter_sql(const char *events, const char *reference,
    duckvep_sql_text *out)
{
    return duckvep_sql_append(out, str_head) &&
        duckvep_sql_relation(out, events) && duckvep_sql_append(out, str_source) && duckvep_sql_append(out, str_sql) &&
        duckvep_sql_relation(out, reference) && duckvep_sql_append(out, str_tail) &&
        duckvep_sql_append(out, str_tail2);
}

/* ---- VCF records to the events relation ------------------------------------ */

/* Records are read_bcf-shaped (CHROM, POS, ID, REF, ALT[], optional INFO_END,
 * INFO_SVTYPE and INFO_SVLEN). Dispatch and coordinates follow VEP 116
 * Parser::VCF::create_VariationFeatures and BaseVCF4 get_start/get_end: a
 * structural record starts at POS + 1 and ends at POS + |SVLEN|, else END, else
 * its start. Records outside the supported families keep a row with a reason. */
static const char vcf_events_head[] =
"WITH __duckvep_ve_source AS (FROM (SELECT NULL::VARCHAR AS \"CHROM\", NULL::BIGINT AS \"POS\", NULL::VARCHAR AS \"ID\", NULL::VARCHAR AS \"REF\", NULL::VARCHAR[] AS \"ALT\", NULL AS \"INFO_END\", NULL AS \"INFO_SVTYPE\", NULL AS \"INFO_SVLEN\" WHERE false) UNION ALL BY NAME SELECT COLUMNS('(?i)^(chrom|pos|id|ref|alt|info_end|info_svtype|info_svlen)$') FROM ";

static const char vcf_events_middle[] =
"),\n"
"__duckvep_ve_regions AS (SELECT CAST(seq_region_name AS VARCHAR) AS seq_region_name, min(CAST(seq_region AS UINTEGER)) AS seq_region, count(*) AS name_rows FROM ";

static const char vcf_events_records[] =
" GROUP BY 1),\n"
"__duckvep_ve_records AS (SELECT CAST(\"CHROM\" AS VARCHAR) AS source_chrom, CAST(\"POS\" AS BIGINT) AS source_pos, CAST(\"ID\" AS VARCHAR) AS source_id, CAST(\"REF\" AS VARCHAR) AS source_ref, CAST(\"ALT\" AS VARCHAR[]) AS alts,\n"
"  " DUCKVEP_VCF_FIRST("\"INFO_END\"") " AS end_text,\n"
"  " DUCKVEP_VCF_FIRST("\"INFO_SVLEN\"") " AS svlen_text,\n"
"  upper(" DUCKVEP_VCF_FIRST("\"INFO_SVTYPE\"") ") AS svtype FROM __duckvep_ve_source),\n"
"__duckvep_ve_classified AS (SELECT *, len(coalesce(alts, [])) AS alt_count,\n"
"  NOT regexp_full_match(coalesce(source_ref, '') || coalesce(array_to_string(alts, ''), ''), '[ACGT]+')\n"
"  AND (svtype IS NOT NULL OR regexp_matches(coalesce(array_to_string(alts, ','), ''), '[<\\[\\]][^*]+[>\\]\\[]|^\\.\\w+|\\w+\\.$')) AS structural_record,\n"
"  " DUCKVEP_VCF_INTEGER("end_text") " AS info_end, " DUCKVEP_VCF_INTEGER("svlen_text") " AS info_svlen FROM __duckvep_ve_records),\n"
"__duckvep_ve_alleles AS (SELECT c.*, a.source_alt, a.alt_index FROM __duckvep_ve_classified c CROSS JOIN UNNEST(CASE WHEN c.alt_count = 0 THEN [NULL::VARCHAR] ELSE c.alts END) WITH ORDINALITY AS a(source_alt, alt_index)),\n";

static const char vcf_events_types[] =
"__duckvep_ve_typed AS (SELECT l.*, r.seq_region AS local_region, r.name_rows AS local_rows,\n"
"  CASE WHEN structural_record AND alt_count = 1 THEN TRY(duckvep_breakend_geometry(source_alt)) END AS bnd,\n"
"  CASE WHEN regexp_matches(source_alt, '^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)') THEN CASE source_alt WHEN '<DEL>' THEN 'DEL' WHEN '<DUP>' THEN 'DUP' WHEN '<DUP:TANDEM>' THEN 'TDUP' WHEN '<INV>' THEN 'INV' WHEN '<INS>' THEN 'INS' WHEN '<CNV>' THEN 'CNV' ELSE 'unsupported_symbolic_allele' END\n"
"       WHEN regexp_full_match(source_alt, '<[^>]+>') THEN CASE WHEN svtype IS NULL THEN 'unsupported_symbolic_allele' WHEN svtype IN ('DEL', 'DUP', 'INV', 'INS', 'CNV') THEN svtype ELSE 'unsupported_structural_type' END\n"
"       ELSE NULL END AS sv_type\n"
"  FROM __duckvep_ve_alleles l LEFT JOIN __duckvep_ve_regions r ON r.seq_region_name = l.source_chrom),\n"
"__duckvep_ve_spans AS (SELECT t.*, m.seq_region AS mate_region, m.name_rows AS mate_rows,\n"
"  " DUCKVEP_VEP_SV_START("source_pos") " AS vep_start,\n"
"  " DUCKVEP_VEP_SV_END("source_pos", "info_svlen", "info_end") " AS vep_end\n"
"  FROM __duckvep_ve_typed t LEFT JOIN __duckvep_ve_regions m ON m.seq_region_name = t.bnd.mate_chrom),\n";

static const char vcf_events_decisions[] =
"__duckvep_ve_decided AS (SELECT *, CASE\n"
"  WHEN source_chrom IS NULL OR local_region IS NULL THEN 'contig_not_in_model'\n"
"  WHEN local_rows > 1 THEN 'ambiguous_contig'\n"
"  WHEN source_pos IS NULL OR source_pos < 1 OR source_pos > 4294967294 THEN 'invalid_position'\n"
"  WHEN source_alt IS NULL THEN 'no_alternate_allele'\n"
"  WHEN source_ref IS NULL OR source_ref = '' THEN 'missing_reference'\n"
"  WHEN NOT structural_record THEN NULL\n"
"  WHEN alt_count <> 1 THEN 'multiallelic_structural_record'\n"
"  WHEN bnd IS NOT NULL AND bnd.mate_chrom IS NULL THEN 'single_breakend'\n"
"  WHEN bnd IS NOT NULL AND (mate_region IS NULL OR mate_rows > 1) THEN 'mate_contig_not_in_model'\n"
"  WHEN bnd IS NOT NULL AND (bnd.mate_position < 1 OR bnd.mate_position > 4294967294) THEN 'invalid_position'\n"
"  WHEN bnd IS NOT NULL THEN NULL\n"
"  WHEN regexp_matches(source_alt, '[\\[\\]]') OR regexp_matches(source_alt, '^\\.\\w+|\\w+\\.$') THEN 'malformed_breakend'\n"
"  WHEN sv_type IS NULL THEN 'literal_alleles_with_svtype'\n"
"  WHEN sv_type IN ('unsupported_symbolic_allele', 'unsupported_structural_type') THEN sv_type\n"
"  WHEN end_text IS NOT NULL AND info_end IS NULL THEN 'invalid_end'\n"
"  WHEN svlen_text IS NOT NULL AND info_svlen IS NULL THEN 'invalid_svlen'\n"
"  WHEN sv_type = 'DEL' AND (vep_start > vep_end OR (info_end IS NULL AND info_svlen IS NULL)) THEN 'incomplete_deletion'\n"
"  WHEN sv_type = 'INS' AND vep_end = vep_start THEN 'single_base_insertion_span'\n"
"  WHEN sv_type = 'INS' AND vep_end < source_pos THEN 'invalid_span'\n"
"  WHEN sv_type <> 'INS' AND vep_end < vep_start THEN 'empty_span'\n"
"  WHEN vep_end > 4294967294 THEN 'invalid_position'\n"
"  ELSE NULL END AS event_reason FROM __duckvep_ve_spans),\n";

static const char vcf_events_select[] =
"__duckvep_ve_events AS (SELECT *, event_reason IS NULL AS ok,\n"
"  CASE WHEN NOT structural_record THEN 'small_variant' WHEN bnd IS NOT NULL THEN 'breakend' WHEN sv_type IS NOT NULL AND NOT starts_with(sv_type, 'unsupported') THEN 'structural_variant' END AS event_kind,\n"
"  CASE WHEN NOT structural_record OR bnd IS NOT NULL THEN source_pos WHEN sv_type = 'INS' AND vep_end < vep_start THEN source_pos ELSE vep_start END AS event_position,\n"
"  CASE WHEN structural_record AND bnd IS NULL THEN CASE WHEN sv_type = 'INS' AND vep_end < vep_start THEN source_pos ELSE vep_end END END AS event_end\n"
"  FROM __duckvep_ve_decided)\n"
"SELECT CAST(row_number() OVER (ORDER BY CASE WHEN ok THEN local_region END NULLS LAST, CASE WHEN ok THEN event_position END NULLS LAST, CASE WHEN ok THEN event_end END NULLS FIRST, source_chrom, source_pos, source_id, alt_index) AS UBIGINT) AS event_index,\n"
"  CASE WHEN ok THEN local_region END AS seq_region,\n"
"  CASE WHEN ok THEN CAST(event_position AS UBIGINT) END AS position,\n"
"  CASE WHEN ok THEN source_ref END AS reference,\n"
"  CASE WHEN ok THEN source_alt END AS alternate,\n"
"  CASE WHEN ok THEN CAST(event_end AS UBIGINT) END AS end_position,\n"
"  CASE WHEN ok AND event_kind = 'structural_variant' THEN sv_type END AS structural_type,\n"
"  CASE WHEN ok AND event_kind = 'structural_variant' THEN CASE sv_type WHEN 'DEL' THEN 'LOSS' WHEN 'DUP' THEN 'GAIN' WHEN 'TDUP' THEN 'GAIN' WHEN 'INV' THEN 'NEUTRAL' ELSE 'UNKNOWN' END END AS copy_change,\n"
"  CASE WHEN ok AND event_kind = 'breakend' THEN mate_region END AS mate_seq_region,\n"
"  CASE WHEN ok AND event_kind = 'breakend' THEN CAST(bnd.mate_position AS UBIGINT) END AS mate_position,\n"
"  event_kind, CASE WHEN ok THEN 'ok' ELSE 'unsupported' END AS event_status, event_reason,\n"
"  source_chrom, source_pos, source_id, source_ref, source_alt, CAST(alt_index AS UINTEGER) AS alt_index\n"
"FROM __duckvep_ve_events ORDER BY event_index";

/* The tail is split by common table expression: ISO C99 caps a string literal at
 * 4,095 characters. */
static const char *const vcf_events_tail[] = {
    vcf_events_records, vcf_events_types, vcf_events_decisions, vcf_events_select
};

bool duckvep_core_vcf_events_sql(const char *records, const char *regions, duckvep_sql_text *out)
{
    if (!records || !regions || !duckvep_sql_append(out, vcf_events_head) ||
        !duckvep_sql_relation(out, records) || !duckvep_sql_append(out, vcf_events_middle) ||
        !duckvep_sql_relation(out, regions))
        return false;
    for (size_t i = 0; i < sizeof(vcf_events_tail) / sizeof(vcf_events_tail[0]); i++)
        if (!duckvep_sql_append(out, vcf_events_tail[i])) return false;
    return true;
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

/* Records are read_bcf-shaped: event_index, CHROM, POS, ID and ALT (a list, with the
 * 1-based alt_index choosing the allele, or one allele), and the optional INFO_MATEID
 * and INFO_EVENT read as value lists. */
static const char pairs_head[] =
"WITH __duckvep_bnd_source AS (FROM (SELECT NULL AS alt_index, NULL AS \"INFO_MATEID\", NULL AS \"INFO_EVENT\" "
"WHERE false) UNION ALL BY NAME "
"SELECT COLUMNS('(?i)^(event_index|chrom|pos|id|alt|alt_index|info_mateid|info_event)$') FROM ";

static const char pairs_body_0[] =
"), __duckvep_bnd_rows AS (SELECT event_index, chrom, pos, id, alt, alt_index, \"INFO_MATEID\", \"INFO_EVENT\" "
"FROM __duckvep_bnd_source), "
"source AS (SELECT CAST(event_index AS BIGINT) AS event_index, chrom, "
"CASE WHEN try_cast(pos AS DOUBLE) = floor(try_cast(pos AS DOUBLE)) THEN try_cast(pos AS BIGINT) END AS pos, "
"CAST(id AS VARCHAR) AS id, " DUCKVEP_VCF_ALT " AS alt, "
"list_filter(coalesce(" DUCKVEP_VCF_VALUES("\"INFO_MATEID\"") ", []), lambda x: nullif(trim(x), '') IS NOT NULL) AS mate_tokens, "
"list_filter(coalesce(" DUCKVEP_VCF_VALUES("\"INFO_EVENT\"") ", []), lambda x: nullif(trim(x), '') IS NOT NULL) AS event_tokens "
"FROM __duckvep_bnd_rows), "
"tokened AS (SELECT *, TRY(duckvep_breakend_geometry(alt)) AS g, "
"CASE WHEN id IS NULL OR id IN ('', '.') THEN NULL ELSE id END AS rid FROM source), "
"";

static const char pairs_body_0b[] =
"records AS (SELECT event_index, CAST(chrom AS VARCHAR) AS chrom, pos, rid AS id, g, "
"CASE WHEN g IS NOT NULL THEN CASE WHEN g.mate_chrom IS NULL THEN 'single_breakend' "
"ELSE 'paired_breakend' END "
"WHEN alt IS NOT NULL AND (regexp_matches(alt, '[\\[\\]]') OR "
"(length(alt) > 1 AND (starts_with(alt, '.') OR ends_with(alt, '.')))) THEN 'malformed_alt' "
"ELSE 'not_breakend' END AS record_kind, "
"length(mate_tokens) AS mate_token_count, length(event_tokens) AS event_token_count, "
"CASE WHEN length(mate_tokens) > 0 THEN array_to_string(list_transform(mate_tokens, lambda x: trim(x)), ',') END AS mate_id, "
"CASE WHEN length(event_tokens) = 1 THEN trim(event_tokens[1]) END AS event_id, "
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
    {PIECE_TEXT, pairs_body_0b},
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

/* Records are read_bcf-shaped: event_index, CHROM, POS, REF and ALT (a list, with the
 * 1-based alt_index choosing the allele, or one allele), and the optional INFO_END,
 * INFO_SVLEN, INFO_CIPOS, INFO_CIEND and INFO_IMPRECISE. The span follows VEP's
 * get_end; an END and SVLEN that disagree are refused, not reconciled. */
static const char hgvs_head[] =
"WITH __duckvep_shgvs_source AS (FROM (SELECT NULL AS alt_index, NULL AS \"INFO_END\", NULL AS \"INFO_SVLEN\", "
"NULL AS \"INFO_CIPOS\", NULL AS \"INFO_CIEND\", NULL AS \"INFO_IMPRECISE\" WHERE false) UNION ALL BY NAME "
"SELECT COLUMNS('(?i)^(event_index|chrom|pos|ref|alt|alt_index|info_end|info_svlen|info_cipos|info_ciend|info_imprecise)$') FROM ";

static const char hgvs_source[] =
"), __duckvep_shgvs_rows AS (SELECT event_index, chrom, pos, ref, alt, alt_index, \"INFO_END\", \"INFO_SVLEN\", "
"\"INFO_CIPOS\", \"INFO_CIEND\", \"INFO_IMPRECISE\" FROM __duckvep_shgvs_source), "
"source AS (SELECT CAST(event_index AS BIGINT) AS event_index, chrom, pos AS raw_pos, "
"CASE WHEN try_cast(pos AS DOUBLE) = floor(try_cast(pos AS DOUBLE)) THEN try_cast(pos AS BIGINT) END AS pos, "
"ref, " DUCKVEP_VCF_ALT " AS alt, "
DUCKVEP_VCF_FIRST("\"INFO_END\"") " AS end_text, " DUCKVEP_VCF_FIRST("\"INFO_SVLEN\"") " AS svlen_text, "
DUCKVEP_VCF_VALUES("\"INFO_CIPOS\"") " AS cipos, " DUCKVEP_VCF_VALUES("\"INFO_CIEND\"") " AS ciend, "
"coalesce(try_cast(" DUCKVEP_VCF_FIRST("\"INFO_IMPRECISE\"") " AS BOOLEAN), false) AS imprecise_flag "
"FROM __duckvep_shgvs_rows";

static const char hgvs_mid[] =
"), refs AS (SELECT CAST(event_index AS BIGINT) AS event_index, count(*) AS n, min(upper(reference_sequence)) AS refseq FROM ";

static const char hgvs_a[] =
" GROUP BY 1), "
"parsed AS (SELECT s.*, r.n AS ref_rows, r.refseq, "
"TRY(duckvep_breakend_geometry(s.alt)) AS g FROM source s LEFT JOIN refs r USING (event_index)), "
"spans AS (SELECT *, "
"CASE WHEN regexp_full_match(end_text, '[1-9][0-9]{0,15}') THEN CAST(end_text AS BIGINT) END AS info_end, "
DUCKVEP_VCF_INTEGER("svlen_text") " AS info_svlen, "
"CASE WHEN alt IN ('<DEL>', '<DUP>', '<DUP:TANDEM>', '<INV>') THEN alt END AS symbolic_kind FROM parsed), "
"bounds AS (SELECT *, CASE WHEN info_end IS NOT NULL OR info_svlen IS NOT NULL THEN "
DUCKVEP_VEP_SV_END("pos", "info_svlen", "info_end") " END AS end_pos FROM spans), "
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
"WHEN imprecise_flag OR cipos IS NOT NULL OR ciend IS NOT NULL THEN 'unsupported:imprecise' "
"WHEN end_text IS NULL AND svlen_text IS NULL THEN 'unavailable:missing_end' "
"WHEN (end_text IS NOT NULL AND info_end IS NULL) OR (svlen_text IS NOT NULL AND info_svlen IS NULL) "
"THEN 'unavailable:invalid_end' "
"WHEN info_end IS NOT NULL AND info_svlen IS NOT NULL AND info_end <> pos + abs(info_svlen) "
"THEN 'unavailable:inconsistent_end' "
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
    {PIECE_TEXT, hgvs_source},
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
