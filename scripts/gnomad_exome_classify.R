# Exome panel classification, shared by gnomad_exome_panel.R and its equivalence test.
#
# exome_prepare(execute) builds `intervals` from the attached model (`model`) and the disjoint,
# sorted unions `cds_union` and `splice_union` (splice windows are exon boundary +/- 8 bp).
# exome_classify_sql classifies `distinct_alleles` with one ASOF lookup per allele (the last
# union interval starting at or before the allele's last base) instead of scanning every
# interval on its contig; the bins are those of the original EXISTS formulation.
exome_prepare <- function(execute) {
  execute("CREATE TEMP TABLE intervals AS SELECT t.seq_region,
  CASE WHEN t.cds_start IS NOT NULL AND t.cds_end IS NOT NULL
    THEN greatest(e.exon_start, t.cds_start) END AS cds_first,
  CASE WHEN t.cds_start IS NOT NULL AND t.cds_end IS NOT NULL
    THEN least(e.exon_end, t.cds_end) END AS cds_last,
  CASE WHEN e.exon_cdna_start > 1 THEN e.exon_start END AS exon_first,
  CASE WHEN e.exon_cdna_end < max(e.exon_cdna_end) OVER (PARTITION BY e.transcript_index)
    THEN e.exon_end END AS exon_last
  FROM model.duckvep_exons e JOIN model.model_transcripts t USING(transcript_index)
  WHERE t.seq_region IN (SELECT DISTINCT seq_region FROM source_alleles)")
  merge <- function(name, source) execute(paste0("CREATE TEMP TABLE ", name, " AS
  WITH s AS (SELECT seq_region, first, last, max(last) OVER (PARTITION BY seq_region
      ORDER BY first, last ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_max FROM (", source, ")),
  g AS (SELECT *, sum(CASE WHEN prev_max IS NULL OR first > prev_max THEN 1 ELSE 0 END)
      OVER (PARTITION BY seq_region ORDER BY first, last) AS grp FROM s)
  SELECT seq_region, min(first) AS first, max(last) AS last FROM g GROUP BY seq_region, grp"))
  merge("cds_union", "SELECT seq_region, cds_first AS first, cds_last AS last FROM intervals
  WHERE cds_first IS NOT NULL AND cds_first <= cds_last")
  merge("splice_union", "SELECT seq_region, greatest(1, exon_first - 8) AS first, exon_first + 8 AS last
  FROM intervals WHERE exon_first IS NOT NULL UNION ALL
  SELECT seq_region, greatest(1, exon_last - 8), exon_last + 8 FROM intervals WHERE exon_last IS NOT NULL")
}
exome_classify_sql <- "SELECT * FROM (WITH literal AS (
 SELECT *, position + length(reference) - 1 AS last_base FROM distinct_alleles WHERE status = 'literal'
 ), annotated AS (
 SELECT v.* EXCLUDE(last_base), coalesce(c.last >= v.position, false) AS in_cds,
   coalesce(sp.last >= v.position, false) AS near_splice
 FROM literal v
 ASOF LEFT JOIN cds_union c ON v.seq_region = c.seq_region AND v.last_base >= c.first
 ASOF LEFT JOIN splice_union sp ON v.seq_region = sp.seq_region AND v.position >= sp.first
 ) SELECT *, CASE
 WHEN in_cds AND length(reference) = 1 AND length(alternate) = 1 THEN 'cds_snv'
 WHEN in_cds AND length(reference) != length(alternate) THEN 'cds_indel'
 WHEN in_cds AND length(reference) = length(alternate) AND length(reference) > 1 THEN 'cds_mnv'
 WHEN near_splice THEN 'splice_remainder'
 ELSE 'other' END AS bin FROM annotated) WHERE bin <> 'other'"
