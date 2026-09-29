#!/usr/bin/env Rscript
# Compare the GFF3 transcript envelope, exon boundaries and each CDS segment
# with the Ensembl model. The latter derives CDS segments by clipping exons
# to the model's genomic coding interval.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
stopifnot(requireNamespace('DBI', quietly = TRUE), requireNamespace('duckdb', quietly = TRUE))
con <- DBI::dbConnect(duckdb::duckdb(dbdir = args[1], read_only = TRUE, shared_home = FALSE))
q <- function(x) as.character(DBI::dbQuoteString(con, normalizePath(x, mustWork = FALSE)))
sql <- sprintf("COPY (
WITH ids AS (
  SELECT column0 AS tx FROM read_csv(%s, header=false, columns={'column0':'VARCHAR'})
), features AS (
  SELECT g.* FROM read_csv(%s, header=false, delim='\\t',
    columns={'tx':'VARCHAR','chr':'VARCHAR','kind':'VARCHAR','s':'BIGINT',
             'e':'BIGINT','strand':'VARCHAR','phase':'VARCHAR'}) g JOIN ids USING (tx)
), gg AS (
  SELECT tx, any_value(chr) FILTER (WHERE kind='mRNA') AS chr,
    min(s) FILTER (WHERE kind='mRNA') AS tx_start,
    max(e) FILTER (WHERE kind='mRNA') AS tx_end,
    any_value(strand) FILTER (WHERE kind='mRNA') AS strand,
    string_agg(CAST(s AS VARCHAR)||'-'||CAST(e AS VARCHAR), ',' ORDER BY s,e)
      FILTER (WHERE kind='exon') AS exons,
    string_agg(CAST(s AS VARCHAR)||'-'||CAST(e AS VARCHAR)||':'||phase, ',' ORDER BY s,e)
      FILTER (WHERE kind='CDS') AS cds
  FROM features GROUP BY tx
), mm AS (
  SELECT m.transcript_stable_id AS tx, m.seq_region_name AS chr,
    m.transcript_start AS tx_start, m.transcript_end AS tx_end,
    CASE m.strand WHEN 1 THEN '+' ELSE '-' END AS strand,
    array_to_string(list_transform(list_sort(list_transform(m.exons,
      lambda x: struct_pack(s:=x.exon_start, e:=x.exon_end))),
      lambda x: CAST(x.s AS VARCHAR)||'-'||CAST(x.e AS VARCHAR)), ',') AS exons,
    array_to_string(list_transform(list_sort(list_filter(list_transform(m.exons, lambda x:
      CASE WHEN m.cds_start IS NOT NULL AND x.exon_end >= m.cds_start AND x.exon_start <= m.cds_end
      THEN struct_pack(s:=greatest(x.exon_start,m.cds_start), e:=least(x.exon_end,m.cds_end),
        phase:=CASE WHEN x.phase < 0 THEN 0 ELSE (3 - x.phase) %% 3 END)
      ELSE NULL END), lambda x: x IS NOT NULL)),
      lambda x: CAST(x.s AS VARCHAR)||'-'||CAST(x.e AS VARCHAR)||':'||CAST(x.phase AS VARCHAR)), ',') AS cds
  FROM model_transcripts m JOIN ids ON ids.tx = m.transcript_stable_id
)
SELECT ids.tx,
  CASE WHEN gg.tx IS NULL OR gg.chr IS NULL OR gg.exons IS NULL OR gg.cds IS NULL THEN 'absent_gff'
       WHEN mm.tx IS NULL THEN 'absent_model'
       WHEN gg.chr=mm.chr AND gg.tx_start=mm.tx_start AND gg.tx_end=mm.tx_end
         AND gg.strand=mm.strand AND gg.exons=mm.exons AND gg.cds=mm.cds THEN 'identical'
       ELSE 'geometry_mismatch' END AS parity
FROM ids LEFT JOIN gg USING (tx) LEFT JOIN mm USING (tx)
ORDER BY ids.tx
) TO %s (HEADER, DELIMITER '\\t')", q(args[3]), q(args[2]), q(args[4]))
DBI::dbExecute(con, sql)
DBI::dbDisconnect(con, shutdown = TRUE)
