#!/usr/bin/env Rscript
# Whole-input check that duckvep_coding_transcripts returns exactly the (event, transcript) pairs that the public
# annotation builder reports as overlapping a CDS (region_mask bit 16), no more and no fewer, and reports how many of the records the committed csq
# accounting compares are discovered (untimed; needs the regenerated ledger).
#
#   Rscript benchmarks/haplotype_scale/verify_discovery.R WORK/hg002_domain.counts.tsv.records.tsv.gz [EXTENSION]
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(trailingOnly = TRUE)
ledger <- args[[1L]]
here <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
root <- normalizePath(file.path(dirname(here), "..", ".."))
extension <- if (length(args) >= 2L) args[[2L]] else file.path(root, "build", "release", "duckvep.duckdb_extension")
vcf <- "/root/duckvep/data/hg002-csq/hg002.ens.vcf.gz"
model <- "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb"
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true", threads = "4")))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
q <- function(x) as.character(dbQuoteString(con, x))
sql <- function(s) invisible(dbExecute(con, s))
get <- function(s) dbGetQuery(con, s)
sql(paste("LOAD", q(extension))); sql(paste0("ATTACH ", q(model), " AS m (READ_ONLY)"))
sql("CREATE TABLE regions AS SELECT seq_region::BIGINT AS seq_region, seq_region_name, sequence_length FROM m.model_regions")
sql("SELECT loaded FROM duckvep_model_load('hap', 'SELECT seq_region::UINTEGER AS seq_region, sequence_length FROM regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM m.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM m.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM m.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid FROM m.model_transcripts, unnest(peptide_edits) u(x)',
 transcript_coverage_complete := TRUE)")
n <- as.integer(system2("sh", c("-c", shQuote(paste0("zcat ", shQuote(vcf), " | head -n 5000 | grep -c '^#'"))), stdout = TRUE))
sql(sprintf("CREATE TABLE ev AS SELECT (row_number() OVER ())::UBIGINT AS event_index, (v.record_index)::UBIGINT AS record_index, r.seq_region::UINTEGER AS seq_region, v.pos::UBIGINT AS position,
  v.ref AS reference, a.alt AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
  FROM (SELECT row_number() OVER () AS record_index, chrom, pos, ref, alt FROM read_csv(%s, delim='\\t', header=false, skip=%d, auto_detect=false, quote='', escape='', strict_mode=false,
  columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR','fmt':'VARCHAR','sample':'VARCHAR'}, compression='gzip')) v
  JOIN regions r ON r.seq_region_name = v.chrom, unnest(string_split(v.alt, ',')) a(alt) WHERE v.alt <> '.'", q(vcf), n))
sql("CREATE VIEW ordered_events AS SELECT event_index, seq_region, position, reference, alternate, end_position, structural_type, copy_change, mate_seq_region, mate_position
  FROM ev ORDER BY seq_region, position, event_index")
sql("CREATE TABLE annotate_pairs AS SELECT DISTINCT event_index, transcript_index FROM query(duckvep_annotate_sql('ordered_events', 'hap')) WHERE region_mask & 16 <> 0")
sql("CREATE TABLE scalar_pairs AS SELECT DISTINCT event_index, t AS transcript_index FROM (SELECT event_index, unnest(duckvep_coding_transcripts('hap', seq_region, position, reference, alternate)) AS t FROM ev)")
res <- get("SELECT (SELECT count(*) FROM annotate_pairs) AS annotate_pairs, (SELECT count(*) FROM scalar_pairs) AS scalar_pairs,
  (SELECT count(*) FROM annotate_pairs a ANTI JOIN scalar_pairs s USING (event_index, transcript_index)) AS only_annotate,
  (SELECT count(*) FROM scalar_pairs s ANTI JOIN annotate_pairs a USING (event_index, transcript_index)) AS only_scalar")
print(res)
sql(sprintf("CREATE TABLE ledger AS SELECT record_index, category FROM read_csv(%s, delim='\\t', header=true, columns={'record_index':'BIGINT','chrom':'VARCHAR','pos':'BIGINT','ref':'VARCHAR','alt':'VARCHAR','category':'VARCHAR','reason':'VARCHAR'})", q(ledger)))
cmp <- get("SELECT count(*) AS compared, count(*) FILTER (WHERE d.record_index IS NOT NULL) AS discovered
  FROM ledger l LEFT JOIN (SELECT DISTINCT record_index FROM ev JOIN scalar_pairs USING (event_index)) d ON d.record_index = l.record_index WHERE l.category = 'compared'")
print(cmp)
stopifnot(res$only_annotate == 0, res$only_scalar == 0)
cat(sprintf("discovery verified: identical (event, transcript) pairs (0 missing, 0 extra); %d of %d compared records are discovered (the others touch a CDS only through an indel anchor base, which the annotation builder does not count)\n", cmp$discovered, cmp$compared))
