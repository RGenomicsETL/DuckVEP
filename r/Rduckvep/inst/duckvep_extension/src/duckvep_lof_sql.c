#include "duckvep_builder.h"
#include "kernel/src/duckvep_budget.h"
#include "duckvep_registration.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
DUCKDB_EXTENSION_EXTERN

/*
 * duckvep_lof_sql: the LOFTEE loss-of-function relation as SQL over rows DuckVEP
 * already produces. The builder only names the relations and the constants; the
 * rules are ordinary joins and CASE expressions (see docs/functions.md and
 * benchmarks/duckvep_lof.md). The rules follow konradjk/loftee at a46b502.
 */

/* BEGIN LOF TAIL (generated from the SQL source; edit the SQL, not this block) */
static const char lof_tail_0[] =
"cand AS (\n"
" SELECT a.event_index, a.transcript_index,\n"
"  CAST(a.geometry.feature_start0 AS BIGINT) + 1 AS vs, CAST(a.geometry.feature_end0 AS BIGINT) AS ve,\n"
"  substring(a.reference, a.geometry.reference_difference_offset + 1, a.geometry.reference_difference_length) AS ref_edit,\n"
"  substring(a.alternate, a.geometry.alternate_difference_offset + 1, a.geometry.alternate_difference_length) AS alt_edit,\n"
"  CAST(a.exon_first AS BIGINT) AS exon_first, CAST(a.exon_last AS BIGINT) AS exon_last, CAST(a.exon_total AS BIGINT) AS exon_total, CAST(a.intron_first AS BIGINT) AS intron_first, CAST(a.intron_total AS BIGINT) AS intron_total, CAST(CASE WHEN a.interbase THEN a.cds_start ELSE a.cds_end END AS BIGINT) AS cds_end, a.cds_start_nf, a.cds_end_nf,\n"
"  t.chrom, t.strand, t.coding_start, t.coding_end, t.transcript_id, t.exons,\n"
"  list_contains(a.cs, 'stop_gained') OR list_contains(a.cs, 'frameshift_variant') AS other_lof,\n"
"  list_contains(a.cs, 'splice_acceptor_variant') AS is_acc,\n"
"  list_contains(a.cs, 'splice_donor_variant') AS is_don\n"
" FROM (SELECT *, string_split(consequence, '&') AS cs FROM ann WHERE transcript_index IS NOT NULL) a\n"
" JOIN tx t USING (transcript_index)\n"
" WHERE t.transcript_biotype = 'protein_coding'\n"
"  AND (list_contains(a.cs, 'stop_gained') OR list_contains(a.cs, 'frameshift_variant')\n"
"   OR list_contains(a.cs, 'splice_acceptor_variant') OR list_contains(a.cs, 'splice_donor_variant'))),\n"
"ex AS (\n"
" SELECT c.event_index, c.transcript_index, c.chrom, c.strand, c.vs, c.coding_start, c.coding_end, c.exon_first,\n"
"  c.other_lof AND c.cds_end IS NOT NULL AS walk,\n"
"  CASE WHEN c.strand > 0 THEN c.coding_end ELSE c.coding_start END AS stop,\n"
"  row_number() OVER (PARTITION BY c.event_index, c.transcript_index\n"
"   ORDER BY CASE WHEN c.strand > 0 THEN u.e.exon_start ELSE -u.e.exon_start END) AS idx,\n"
"  CAST(u.e.exon_start AS BIGINT) AS es, CAST(u.e.exon_end AS BIGINT) AS ee\n"
" FROM cand c, unnest(c.exons) AS u(e)),\n"
"seg AS (\n"
" SELECT *, (strand > 0 AND vs > ee) OR (strand < 0 AND vs < es) AS skipped,\n"
"  (strand > 0 AND es < stop AND ee >= stop) OR (strand < 0 AND ee > stop AND es <= stop) AS is_last,\n"
"  vs >= es AND vs <= ee AS affected\n"
" FROM ex WHERE walk),\n"
"segs AS (\n"
" SELECT *, CASE WHEN is_last THEN (CASE WHEN affected THEN vs WHEN strand > 0 THEN es ELSE ee END)\n"
"   WHEN affected THEN vs ELSE es END AS s0,\n"
"  CASE WHEN is_last THEN stop WHEN affected THEN (CASE WHEN strand > 0 THEN ee ELSE es END) ELSE ee END AS s1\n"
" FROM seg),\n"
"segw AS (\n"
" SELECT g.event_index, g.transcript_index, g.idx,\n"
"  sum(greatest(0, least(greatest(g.s0, g.s1), p.\"end\") - greatest(least(g.s0, g.s1), p.\"start\" + 1) + 1) * p.score) AS w\n"
" FROM segs g JOIN gerp p ON p.chrom = g.chrom AND p.\"start\" < greatest(g.s0, g.s1) AND p.\"end\" >= least(g.s0, g.s1)\n"
" GROUP BY g.event_index, g.transcript_index, g.idx),\n"
"dist AS (\n"
" SELECT g.event_index, g.transcript_index,\n"
"  coalesce(sum(abs(g.s1 - g.s0)) FILTER (WHERE NOT g.skipped), 0) AS bp_dist,\n"
"  max(CASE WHEN g.idx >= g.exon_first AND g.es <= g.stop AND g.ee >= g.stop\n"
"   THEN CASE WHEN g.strand > 0 THEN g.stop - g.es ELSE g.ee - g.stop END END) AS last_len,\n"
"  coalesce(sum(coalesce(w.w, 0)) FILTER (WHERE NOT g.skipped), 0) AS gerp_dist\n"
" FROM segs g LEFT JOIN segw w USING (event_index, transcript_index, idx)\n"
" GROUP BY g.event_index, g.transcript_index),\n"
"cdslen AS (\n"
" SELECT event_index, transcript_index,\n";

static const char lof_tail_1[] =
"  sum(greatest(0, least(ee, coding_end) - greatest(es, coding_start) + 1)) AS cds_len\n"
" FROM ex WHERE walk GROUP BY event_index, transcript_index),\n"
"intr AS (\n"
" SELECT c.event_index, c.transcript_index, c.chrom,\n"
"  CASE WHEN c.strand > 0 THEN a.ee + 1 ELSE b.ee + 1 END AS istart,\n"
"  CASE WHEN c.strand > 0 THEN b.es - 1 ELSE a.es - 1 END AS iend\n"
" FROM cand c\n"
"  JOIN ex a ON a.event_index = c.event_index AND a.transcript_index = c.transcript_index AND a.idx = c.intron_first\n"
"  JOIN ex b ON b.event_index = c.event_index AND b.transcript_index = c.transcript_index AND b.idx = c.intron_first + 1),\n"
"req AS (\n"
" SELECT event_index, transcript_index, 'L' AS tag, chrom, istart AS a, istart + 1 AS b FROM intr\n"
" UNION ALL SELECT event_index, transcript_index, 'R', chrom, iend - 1, iend FROM intr\n"
" UNION ALL SELECT event_index, transcript_index, 'W', chrom, vs - 4, ve + 4 FROM cand\n"
"  WHERE is_acc AND intron_first IS NOT NULL AND ve = vs),\n"
"got AS (\n"
" SELECT r.event_index, r.transcript_index, r.tag,\n"
"  CASE WHEN sum(least(r.b, c.\"end\") - greatest(r.a, c.\"start\" + 1) + 1) = r.b - r.a + 1\n"
"   THEN upper(string_agg(substring(c.seq, greatest(r.a, c.\"start\" + 1) - c.\"start\",\n"
"    least(r.b, c.\"end\") - greatest(r.a, c.\"start\" + 1) + 1), '' ORDER BY c.\"start\")) END AS s\n"
" FROM req r JOIN refs c ON c.chrom = r.chrom AND c.\"start\" < r.b AND c.\"end\" >= r.a\n"
" GROUP BY r.event_index, r.transcript_index, r.tag, r.a, r.b),\n"
"motif AS (\n"
" SELECT event_index, transcript_index,\n"
"  max(CASE WHEN tag = 'L' THEN s END) AS left2, max(CASE WHEN tag = 'R' THEN s END) AS right2,\n"
"  max(CASE WHEN tag = 'W' THEN s END) AS win\n"
" FROM got GROUP BY event_index, transcript_index),\n"
"isz AS (SELECT event_index, transcript_index, iend - istart + 1 AS isize FROM intr),\n"
"anc_got AS (\n"
" SELECT c.event_index, c.transcript_index, upper(substring(n.seq, c.vs - n.\"start\", 1)) AS base\n"
" FROM cand c JOIN anc n ON n.chrom = c.chrom AND n.\"start\" < c.vs AND n.\"end\" >= c.vs\n"
" WHERE length(c.ref_edit) = 1 AND length(c.alt_edit) = 1),\n"
"pc AS (SELECT transcript, exon, min(orf) AS orf, min(maxs) AS maxs FROM pcsf GROUP BY transcript, exon),\n"
"facts AS (\n"
" SELECT c.*, cfg.*, d.bp_dist, d.last_len, d.gerp_dist, l.cds_len, z.isize, ag.base AS anc_base,\n"
"  p.orf, p.maxs, p.transcript IS NOT NULL AS pcsf_row,\n"
"  CASE WHEN c.strand > 0 THEN m.left2 ELSE reverse(translate(m.right2, 'ACGT', 'TGCA')) END AS first2,\n"
"  CASE WHEN c.strand > 0 THEN m.right2 ELSE reverse(translate(m.left2, 'ACGT', 'TGCA')) END AS last2,\n"
"  CASE WHEN c.strand > 0 THEN c.ref_edit ELSE reverse(translate(c.ref_edit, 'ACGT', 'TGCA')) END AS ref_t,\n"
"  CASE WHEN c.strand > 0 THEN c.alt_edit ELSE reverse(translate(c.alt_edit, 'ACGT', 'TGCA')) END AS alt_t,\n"
"  CASE WHEN c.strand > 0 THEN m.win ELSE reverse(translate(m.win, 'ACGT', 'TGCA')) END AS win_t,\n"
"  c.other_lof AND c.cds_end IS NOT NULL AS end_block,\n"
"  d.bp_dist IS NOT NULL AS have_geom,\n"
"  c.other_lof AND c.exon_first IS NOT NULL AS exon_block,\n"
"  c.intron_first IS NOT NULL AND c.intron_total IS NOT NULL AS intron_block,\n"
"  (c.strand > 0 AND c.ve < c.coding_start) OR (c.strand < 0 AND c.vs > c.coding_end) AS utr5,\n"
"  (c.strand > 0 AND c.vs > c.coding_end) OR (c.strand < 0 AND c.ve < c.coding_start) AS utr3,\n"
"  c.is_acc OR c.is_don AS splice_lof\n"
" FROM cand c CROSS JOIN cfg\n"
"  LEFT JOIN dist d USING (event_index, transcript_index)\n"
"  LEFT JOIN cdslen l USING (event_index, transcript_index)\n";

static const char lof_tail_2[] =
"  LEFT JOIN motif m USING (event_index, transcript_index)\n"
"  LEFT JOIN isz z USING (event_index, transcript_index)\n"
"  LEFT JOIN anc_got ag USING (event_index, transcript_index)\n"
"  LEFT JOIN pc p ON p.transcript = c.transcript_id AND p.exon = c.exon_first AND c.exon_first = c.exon_last),\n"
"ev AS (\n"
" SELECT *,\n"
"  CASE WHEN end_block AND have_geom AND exon_first IS NOT NULL AND bp_dist - coalesce(last_len, -1000) <= 50\n"
"    AND (NOT has_gerp OR coalesce(gerp_dist, 0) <= gerp_cut) THEN 'END_TRUNC' END AS f_end,\n"
"  CASE WHEN exon_block AND exon_total IS NULL THEN 'EXON_INTRON_UNDEF'\n"
"   WHEN exon_block AND exon_total <> 1 AND check_cds AND (cds_start_nf OR cds_end_nf) THEN 'INCOMPLETE_CDS' END AS f_exon,\n"
"  CASE WHEN intron_first IS NOT NULL AND intron_total IS NULL THEN 'EXON_INTRON_UNDEF'\n"
"   WHEN intron_first IS NOT NULL AND isize < min_intron THEN 'SMALL_INTRON' END AS f_intron,\n"
"  CASE WHEN intron_block AND is_don AND first2 = 'GC' AND ref_t = 'C' AND alt_t = 'T' THEN 'GC_TO_GT_DONOR' END AS f_gc,\n"
"  CASE WHEN intron_block AND splice_lof AND utr5 THEN '5UTR_SPLICE' END AS f_5utr,\n"
"  CASE WHEN intron_block AND splice_lof AND utr3 THEN '3UTR_SPLICE' END AS f_3utr,\n"
"  CASE WHEN has_anc AND anc_base = alt_edit THEN 'ANC_ALLELE' END AS f_anc,\n"
"  CASE WHEN end_block AND exon_first IS NULL THEN 'NO_EXON_NUMBER' END AS g_noexon,\n"
"  CASE WHEN exon_block AND exon_total = 1 THEN 'SINGLE_EXON' END AS g_single,\n"
"  CASE WHEN exon_block AND has_pcsf AND pcsf_row AND orf < 0\n"
"   THEN (CASE WHEN maxs > 0 THEN 'PHYLOCSF_UNLIKELY_ORF' ELSE 'PHYLOCSF_WEAK' END) END AS g_pcsf,\n"
"  CASE WHEN intron_block AND splice_lof AND (first2 <> 'GT' OR last2 <> 'AG') THEN 'NON_CAN_SPLICE' END AS g_noncan,\n"
"  CASE WHEN intron_block AND is_acc AND regexp_matches(win_t, 'AG.AG') THEN 'NAGNAG_SITE' END AS g_nagnag\n"
" FROM facts),\n"
"res AS (\n"
" SELECT event_index, transcript_index,\n"
"  nullif(concat_ws(',', f_end, f_exon, f_intron, f_gc, f_5utr, f_3utr, f_anc), '') AS lof_filter,\n"
"  nullif(concat_ws(',', g_noexon, g_single, g_pcsf, g_noncan, g_nagnag), '') AS lof_flags,\n"
"  nullif(concat_ws(',',\n"
"   CASE WHEN end_block AND have_geom THEN 'PERCENTILE:' || printf('%.15g', cds_end / cds_len::DOUBLE) END,\n"
"   CASE WHEN end_block AND have_geom AND has_gerp THEN 'GERP_DIST:' || printf('%.15g', gerp_dist) END,\n"
"   CASE WHEN end_block AND have_geom THEN 'BP_DIST:' || bp_dist::VARCHAR END,\n"
"   CASE WHEN end_block AND have_geom AND exon_first IS NOT NULL THEN 'DIST_FROM_LAST_EXON:' || (bp_dist - coalesce(last_len, -1000))::VARCHAR END,\n"
"   CASE WHEN end_block AND have_geom AND exon_first IS NOT NULL THEN '50_BP_RULE:' || CASE WHEN bp_dist - coalesce(last_len, -1000) <= 50 THEN 'FAIL' ELSE 'PASS' END END,\n"
"   CASE WHEN exon_block AND has_pcsf AND pcsf_row THEN 'ANN_ORF:' || printf('%.15g', orf) END,\n"
"   CASE WHEN exon_block AND has_pcsf AND pcsf_row THEN 'MAX_ORF:' || printf('%.15g', maxs) END,\n"
"   CASE WHEN exon_block AND has_pcsf AND NOT pcsf_row THEN 'PHYLOCSF_TOO_SHORT' END,\n"
"   CASE WHEN intron_block THEN 'INTRON_SIZE:' || isize::VARCHAR END), '') AS lof_info,\n"
"  nullif(concat_ws(',',\n"
"   CASE WHEN end_block AND exon_first IS NOT NULL AND NOT have_geom THEN 'END_TRUNC' END,\n"
"   CASE WHEN end_block AND exon_first IS NOT NULL AND have_geom AND NOT has_gerp THEN 'GERP_END_TRUNC' END,\n"
"   CASE WHEN intron_block AND isize IS NULL THEN 'SMALL_INTRON' END,\n";

static const char lof_tail_3[] =
"   CASE WHEN intron_block AND is_don AND first2 IS NULL AND ref_t = 'C' AND alt_t = 'T' THEN 'GC_TO_GT_DONOR' END,\n"
"   CASE WHEN exon_block AND NOT has_pcsf THEN 'PHYLOCSF' END,\n"
"   CASE WHEN NOT has_anc AND length(ref_edit) = 1 AND length(alt_edit) = 1 THEN 'ANC_ALLELE' END,\n"
"   CASE WHEN intron_block AND splice_lof AND (first2 IS NULL OR last2 IS NULL) THEN 'NON_CAN_SPLICE' END,\n"
"   CASE WHEN intron_block AND is_acc AND ve = vs AND win_t IS NULL THEN 'NAGNAG_SITE' END,\n"
"   CASE WHEN has_anc AND length(ref_edit) = 1 AND length(alt_edit) = 1 AND anc_base IS NULL THEN 'ANC_ALLELE' END), '') AS lof_unchecked,\n"
"  CASE WHEN concat_ws(',', f_end, f_exon, f_intron, f_gc, f_5utr, f_3utr, f_anc) <> '' THEN 'LC' ELSE 'HC' END AS lof\n"
" FROM ev)\n"
"SELECT a.event_index, a.transcript_index, r.lof, r.lof_filter, r.lof_flags, r.lof_info, r.lof_unchecked\n"
"FROM (SELECT event_index, transcript_index FROM ann WHERE transcript_index IS NOT NULL) a\n"
" LEFT JOIN res r USING (event_index, transcript_index)";

static const char *const lof_tail[] = {lof_tail_0, lof_tail_1, lof_tail_2, lof_tail_3};
/* END LOF TAIL */

enum { LOF_GERP, LOF_ANCESTOR, LOF_PHYLOCSF, LOF_MIN_INTRON, LOF_CUTOFF, LOF_CHECK_CDS, LOF_OPTIONS };

typedef struct {
    char *relation[3];       /* gerp, ancestor, phylocsf; NULL when omitted */
    int64_t min_intron_size;
    double gerp_cutoff;
    bool check_complete_cds;
} lof_options;

static bool option_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return !validity || duckdb_validity_row_is_valid(validity, row);
}

static bool integer_value(duckdb_vector vector, idx_t row, int64_t *out) {
    duckdb_logical_type type = duckdb_vector_get_column_type(vector);
    duckdb_type id = duckdb_get_type_id(type);
    duckdb_destroy_logical_type(&type);
    void *data = duckdb_vector_get_data(vector);
    switch (id) {
    case DUCKDB_TYPE_TINYINT: *out = ((int8_t *)data)[row]; return true;
    case DUCKDB_TYPE_SMALLINT: *out = ((int16_t *)data)[row]; return true;
    case DUCKDB_TYPE_INTEGER: *out = ((int32_t *)data)[row]; return true;
    case DUCKDB_TYPE_BIGINT: *out = ((int64_t *)data)[row]; return true;
    case DUCKDB_TYPE_UTINYINT: *out = ((uint8_t *)data)[row]; return true;
    case DUCKDB_TYPE_USMALLINT: *out = ((uint16_t *)data)[row]; return true;
    case DUCKDB_TYPE_UINTEGER: *out = ((uint32_t *)data)[row]; return true;
    case DUCKDB_TYPE_UBIGINT: {
        uint64_t big = ((uint64_t *)data)[row];
        if (big > (uint64_t)INT64_MAX) return false;
        *out = (int64_t)big;
        return true;
    }
    default: return false;
    }
}

static bool double_value(duckdb_vector vector, idx_t row, double *out) {
    duckdb_logical_type type = duckdb_vector_get_column_type(vector);
    duckdb_type id = duckdb_get_type_id(type);
    duckdb_destroy_logical_type(&type);
    if (id == DUCKDB_TYPE_DOUBLE) { *out = ((double *)duckdb_vector_get_data(vector))[row]; return true; }
    if (id == DUCKDB_TYPE_FLOAT) { *out = ((float *)duckdb_vector_get_data(vector))[row]; return true; }
    int64_t whole;
    if (!integer_value(vector, row, &whole)) return false;
    *out = (double)whole;
    return true;
}

static void free_options(lof_options *options) {
    for (size_t i = 0; i < 3; i++) duckvep_budget_free(options->relation[i]);
}

static bool read_options(duckdb_function_info info, duckdb_vector vector, idx_t row, lof_options *options) {
    const char *const names[LOF_OPTIONS] = {"gerp", "ancestor", "phylocsf", "min_intron_size",
        "gerp_end_trunc_cutoff", "check_complete_cds"};
    const duckvep_option_kind kinds[LOF_OPTIONS] = {DUCKVEP_OPTION_TEXT, DUCKVEP_OPTION_TEXT,
        DUCKVEP_OPTION_TEXT, DUCKVEP_OPTION_INTEGER, DUCKVEP_OPTION_NUMERIC, DUCKVEP_OPTION_BOOLEAN};
    duckdb_vector fields[LOF_OPTIONS] = {0};
    options->min_intron_size = 15;
    options->gerp_cutoff = -58;
    if (!vector) return true;
    if (!duckvep_builder_option_vectors(info, vector, row, names, kinds, LOF_OPTIONS, fields)) return false;
    for (size_t i = 0; i < 3; i++) {
        if (!fields[i] || !option_valid(fields[i], row)) continue;
        options->relation[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(fields[i]))[row]);
        if (!options->relation[i]) {
            duckvep_builder_set_error(info, "duckvep_lof_sql: invalid option string or allocation failure");
            return false;
        }
    }
    if (fields[LOF_MIN_INTRON] && option_valid(fields[LOF_MIN_INTRON], row) &&
        (!integer_value(fields[LOF_MIN_INTRON], row, &options->min_intron_size) ||
         options->min_intron_size < 0 || options->min_intron_size > 1000000000)) {
        duckvep_builder_set_error(info, "duckvep_lof_sql: min_intron_size must be an integer from 0 through 1000000000");
        return false;
    }
    if (fields[LOF_CUTOFF] && option_valid(fields[LOF_CUTOFF], row) &&
        (!double_value(fields[LOF_CUTOFF], row, &options->gerp_cutoff) || !isfinite(options->gerp_cutoff))) {
        duckvep_builder_set_error(info, "duckvep_lof_sql: gerp_end_trunc_cutoff must be a finite number");
        return false;
    }
    if (fields[LOF_CHECK_CDS] && option_valid(fields[LOF_CHECK_CDS], row))
        options->check_complete_cds = ((uint8_t *)duckdb_vector_get_data(fields[LOF_CHECK_CDS]))[row] != 0;
    return true;
}

/* One optional relation: the named relation, or an empty relation of the same shape. */
static bool optional_relation(duckvep_sql_text *sql, const char *name, const char *relation,
                              const char *columns, const char *empty) {
    if (!duckvep_sql_append(sql, name) || !duckvep_sql_append(sql, " AS (SELECT ")) return false;
    if (!relation) return duckvep_sql_append(sql, empty) && duckvep_sql_append(sql, " WHERE false),\n");
    return duckvep_sql_append(sql, columns) && duckvep_sql_append(sql, " FROM ") &&
        duckvep_sql_relation(sql, relation) && duckvep_sql_append(sql, "),\n");
}

static bool assemble(duckvep_sql_text *sql, char *const *names, const lof_options *options) {
    char number[64];
    bool ok = duckvep_sql_append(sql, "WITH ann AS (SELECT * FROM ") && duckvep_sql_relation(sql, names[0]) &&
        duckvep_sql_append(sql, "),\ntx AS (SELECT transcript_index, seq_region_name AS chrom, strand, "
            "CAST(cds_start AS BIGINT) AS coding_start, CAST(cds_end AS BIGINT) AS coding_end, "
            "transcript_biotype, ") &&
        duckvep_sql_append(sql, options->relation[LOF_PHYLOCSF] ? "CAST(transcript_stable_id AS VARCHAR)" : "NULL::VARCHAR") &&
        duckvep_sql_append(sql, " AS transcript_id, exons FROM ") && duckvep_sql_relation(sql, names[1]) &&
        duckvep_sql_append(sql, "),\nrefs AS (SELECT CAST(chrom AS VARCHAR) AS chrom, CAST(\"start\" AS BIGINT) AS \"start\", "
            "CAST(\"end\" AS BIGINT) AS \"end\", CAST(seq AS VARCHAR) AS seq FROM ") &&
        duckvep_sql_relation(sql, names[2]) && duckvep_sql_append(sql, "),\n") &&
        optional_relation(sql, "gerp", options->relation[LOF_GERP],
            "CAST(chrom AS VARCHAR) AS chrom, CAST(\"start\" AS BIGINT) AS \"start\", CAST(\"end\" AS BIGINT) AS \"end\", CAST(score AS DOUBLE) AS score",
            "NULL::VARCHAR AS chrom, NULL::BIGINT AS \"start\", NULL::BIGINT AS \"end\", NULL::DOUBLE AS score") &&
        optional_relation(sql, "anc", options->relation[LOF_ANCESTOR],
            "CAST(chrom AS VARCHAR) AS chrom, CAST(\"start\" AS BIGINT) AS \"start\", CAST(\"end\" AS BIGINT) AS \"end\", CAST(seq AS VARCHAR) AS seq",
            "NULL::VARCHAR AS chrom, NULL::BIGINT AS \"start\", NULL::BIGINT AS \"end\", NULL::VARCHAR AS seq") &&
        optional_relation(sql, "pcsf", options->relation[LOF_PHYLOCSF],
            "CAST(transcript AS VARCHAR) AS transcript, CAST(exon AS BIGINT) AS exon, CAST(corresponding_orf_score AS DOUBLE) AS orf, CAST(max_score AS DOUBLE) AS maxs",
            "NULL::VARCHAR AS transcript, NULL::BIGINT AS exon, NULL::DOUBLE AS orf, NULL::DOUBLE AS maxs");
    if (!ok) return false;
    snprintf(number, sizeof(number), "%lld", (long long)options->min_intron_size);
    ok = duckvep_sql_append(sql, "cfg AS (SELECT CAST(") && duckvep_sql_append(sql, number) &&
        duckvep_sql_append(sql, " AS BIGINT) AS min_intron, CAST(");
    snprintf(number, sizeof(number), "%.17g", options->gerp_cutoff);
    ok = ok && duckvep_sql_append(sql, number) && duckvep_sql_append(sql, " AS DOUBLE) AS gerp_cut, ") &&
        duckvep_sql_append(sql, options->check_complete_cds ? "true" : "false") &&
        duckvep_sql_append(sql, " AS check_cds, ") &&
        duckvep_sql_append(sql, options->relation[LOF_GERP] ? "true" : "false") &&
        duckvep_sql_append(sql, " AS has_gerp, ") &&
        duckvep_sql_append(sql, options->relation[LOF_ANCESTOR] ? "true" : "false") &&
        duckvep_sql_append(sql, " AS has_anc, ") &&
        duckvep_sql_append(sql, options->relation[LOF_PHYLOCSF] ? "true" : "false") &&
        duckvep_sql_append(sql, " AS has_pcsf),\n");
    for (size_t i = 0; ok && i < sizeof(lof_tail) / sizeof(*lof_tail); i++)
        ok = duckvep_sql_append(sql, lof_tail[i]);
    return ok;
}

static void lof_builder(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    idx_t argc = duckdb_data_chunk_get_column_count(input);
    duckdb_vector args[4];
    for (idx_t i = 0; i < argc; i++) args[i] = duckdb_data_chunk_get_vector(input, i);
    for (idx_t row = 0; row < duckdb_data_chunk_get_size(input); row++) {
        char *names[3] = {0};
        lof_options options;
        memset(&options, 0, sizeof(options));
        bool ok = true;
        for (size_t i = 0; ok && i < 3; i++) {
            ok = option_valid(args[i], row);
            if (ok) names[i] = duckvep_builder_string(((duckdb_string_t *)duckdb_vector_get_data(args[i]))[row]);
            ok = ok && names[i] != NULL;
        }
        if (!ok) duckvep_builder_set_error(info, "duckvep_lof_sql: invalid relation name or allocation failure");
        if (ok) ok = read_options(info, argc == 4 ? args[3] : NULL, row, &options);
        duckvep_sql_text sql = {0};
        if (ok && !assemble(&sql, names, &options)) {
            duckvep_builder_set_error(info, "duckvep_lof_sql: invalid relation name or allocation failure");
            ok = false;
        }
        if (ok) duckdb_vector_assign_string_element_len(output, row, sql.data, sql.length);
        duckvep_sql_free(&sql);
        free_options(&options);
        for (size_t i = 0; i < 3; i++) duckvep_budget_free(names[i]);
        if (!ok) return;
    }
}

bool register_duckvep_lof_sql(duckdb_connection connection) {
    return duckvep_register_builder(connection, "duckvep_lof_sql", 3, lof_builder);
}
