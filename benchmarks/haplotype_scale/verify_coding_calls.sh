#!/usr/bin/env bash
# Proves the fused reader's calls relation equal to mode B's: both are built in one DuckDB CLI session from the same
# VCF and compared as multisets (EXCEPT ALL in both directions, every column), plus the order-independent hash of each.
#
#   verify_coding_calls.sh EXT VCF [MODEL.duckdb]
# Prints "calls_reader=N calls_modeB=M only_reader=a only_modeB=b" and "IDENTICAL" when a = b = 0 and N = M.
set -euo pipefail
ext=$1; vcf=$2; model=${3:-/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb}
skip=$(zcat "$vcf" 2>/dev/null | head -n 5000 | grep -c '^#' || true)
sql=$(mktemp)
cat > "$sql" <<SQL
SET threads = 1;
LOAD '$ext';
ATTACH '$model' AS m (READ_ONLY);
SELECT loaded FROM duckvep_model_load('hap', 'SELECT seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name FROM m.model_regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM m.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM m.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM m.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid FROM m.model_transcripts, unnest(peptide_edits) u(x)',
 transcript_coverage_complete := TRUE);
CREATE TABLE regions AS SELECT seq_region::BIGINT AS seq_region, seq_region_name, sequence_length FROM m.model_regions;
CREATE TABLE reader AS SELECT * FROM duckvep_coding_calls('hap', '$vcf');
CREATE TABLE modeb AS SELECT ((c.record_index << 6) | (c.i - 1))::BIGINT AS event_index, c.seq_region::INTEGER AS seq_region,
 c.pos::BIGINT AS position, c.ref AS reference, c.alt AS alternate, c.i::INTEGER AS alt_index, t.transcript_index::INTEGER AS transcript_index,
 0::INTEGER AS sample_index, list_transform(string_split_regex(split_part(c.sample, ':', 1), '[/|]'), lambda x: try_cast(x AS INTEGER)) AS alleles,
 [false] || list_transform(regexp_extract_all(split_part(c.sample, ':', 1), '[/|]'), lambda s: s = '|') AS phase_before,
 try_cast(list_last(string_split(c.sample, ':')) AS BIGINT) AS phase_set
 FROM (SELECT * FROM (SELECT v.record_index, a.i, r.seq_region, v.pos, v.ref, a.alt AS alt, v.sample,
   duckvep_coding_transcripts('hap', r.seq_region, v.pos, v.ref, a.alt) AS tx
   FROM (SELECT row_number() OVER () AS record_index, chrom, pos, ref, alt, sample FROM read_csv('$vcf', delim='\t', header=false, skip=$skip,
    auto_detect=false, quote='', escape='', strict_mode=false, compression='gzip',
    columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR','fmt':'VARCHAR','sample':'VARCHAR'})) v
   JOIN regions r ON r.seq_region_name = v.chrom, unnest(string_split(v.alt, ',')) WITH ORDINALITY a(alt, i) WHERE v.alt <> '.') WHERE len(tx) > 0) c, unnest(c.tx) t(transcript_index);
SELECT 'calls_reader=' || (SELECT count(*) FROM reader) || ' calls_modeB=' || (SELECT count(*) FROM modeb)
 || ' only_reader=' || (SELECT count(*) FROM (SELECT * FROM reader EXCEPT ALL SELECT * FROM modeb))
 || ' only_modeB=' || (SELECT count(*) FROM (SELECT * FROM modeb EXCEPT ALL SELECT * FROM reader))
 || ' hash_reader=' || (SELECT sum(hash(r)::HUGEINT) FROM reader r) || ' hash_modeB=' || (SELECT sum(hash(r)::HUGEINT) FROM modeb r);
SELECT CASE WHEN (SELECT count(*) FROM reader) = (SELECT count(*) FROM modeb)
 AND NOT EXISTS (SELECT * FROM reader EXCEPT ALL SELECT * FROM modeb) AND NOT EXISTS (SELECT * FROM modeb EXCEPT ALL SELECT * FROM reader)
 THEN 'IDENTICAL' ELSE 'DIFFERENT' END;
SQL
duckdb -unsigned -csv -noheader -f "$sql" | tail -n 2
rm -f "$sql"
