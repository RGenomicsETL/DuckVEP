#!/usr/bin/env bash
# DuckVEP against vep-rs (https://github.com/natera-open-source/vep-rs) on one sites-only VCF, with
# executable Ensembl VEP adjudicating every disagreement.
#
# Both tools read the same file: one ALT allele per record, the ID column a unique key, so their output rows
# join on (ID, transcript). vep-rs writes VEP tab text; DuckVEP writes its rich relation as Parquet.
#
#   benchmark_duckvep_vep_rs.sh --vcf GENOTYPES.vcf.gz --work DIR [--threads 16] [--runs 3]
#
# Environment (all required unless a default is given):
#   DUCKVEP_EXT          the duckvep extension            (default build/release/duckvep.duckdb_extension)
#   DUCKVEP_READER_EXT   a DuckHTS extension, for read_bcf while the input is prepared
#   DUCKVEP_MODEL        DuckVEP model database (homo_sapiens_<release>_GRCh38_final.duckdb)
#   VEP_RS_DIR           directory holding the vep-rs `vep` and `vep-cache-builder` binaries
#   VEP_RS_CACHE         vep-rs JSON cache built by vep-cache-builder from the same Ensembl release's GTF
#   VEP_GTF              that GTF, sorted, bgzipped and tabix-indexed, for the adjudicating VEP run
#   VEP_FASTA            indexed reference FASTA
#   VEP_PREFIX           conda prefix of executable Ensembl VEP (default /root/miniconda3/envs/vep)
#   BENCH_OUT            where timed outputs go; memory-backed storage keeps disk writeback out of the timing
#                        (default /dev/shm/duckvep-vep-rs)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
vcf=; work=; threads=16; runs=3
while [[ $# -gt 0 ]]; do case $1 in
    --vcf) vcf=$2; shift 2;; --work) work=$2; shift 2;; --threads) threads=$2; shift 2;; --runs) runs=$2; shift 2;;
    *) echo "unknown option $1" >&2; exit 2;; esac; done
[[ -n $vcf && -n $work ]] || { echo "need --vcf and --work" >&2; exit 2; }
ext=${DUCKVEP_EXT:-$root/build/release/duckvep.duckdb_extension}
: "${DUCKVEP_READER_EXT:?}" "${DUCKVEP_MODEL:?}" "${VEP_RS_DIR:?}" "${VEP_RS_CACHE:?}" "${VEP_GTF:?}" "${VEP_FASTA:?}"
vep_prefix=${VEP_PREFIX:-/root/miniconda3/envs/vep}
out=${BENCH_OUT:-/dev/shm/duckvep-vep-rs}
mkdir -p "$work" "$out"
cd "$work"

# ---- Untimed preparation: the shared input, a model snapshot, and identifier tables -----------------------
if [[ ! -s sites.vcf ]]; then
duckdb -unsigned <<SQL
LOAD '$DUCKVEP_READER_EXT';
CREATE TABLE sites AS SELECT CHROM, POS, REF, unnest(ALT) AS ALT FROM read_bcf('$vcf');
COPY (SELECT line FROM (
  SELECT 0 o, 0 r, '##fileformat=VCFv4.2' line UNION ALL SELECT 1, 0, '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
  UNION ALL SELECT 2, row_number() OVER (), CHROM || chr(9) || POS || chr(9) || 'e' || row_number() OVER () || chr(9) ||
    REF || chr(9) || ALT || chr(9) || '.	.	.' FROM sites) ORDER BY o, r)
 TO 'sites.vcf' (FORMAT csv, HEADER false, QUOTE '', ESCAPE '', DELIMITER '\x01');
LOAD '$ext';
ATTACH '$DUCKVEP_MODEL' AS m (READ_ONLY);
SELECT loaded FROM duckvep_model_load('g', 'SELECT seq_region, sequence_length FROM m.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM m.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM m.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM m.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM m.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 transcript_coverage_complete := TRUE);
SELECT duckvep_model_save('g', 'model.dvsnap');
COPY (SELECT transcript_index, transcript_stable_id FROM m.model_transcripts) TO 'transcript_ids.parquet' (FORMAT parquet);
COPY (SELECT seq_region, seq_region_name FROM m.model_regions) TO 'regions.parquet' (FORMAT parquet);
SQL
fi

# ---- The timed DuckVEP run: extension load, snapshot restore, VCF read, coordinate sort, annotation, Parquet ----
duckvep_sql() {
cat <<SQL
SET threads=$1; SET memory_limit='8GB';
LOAD '$ext';
SELECT duckvep_model_restore('g', 'model.dvsnap');
CREATE TEMP TABLE ev AS
 SELECT row_number() OVER (ORDER BY r.seq_region, v.pos, v.id)::UBIGINT AS event_index, v.id AS vcf_id,
        r.seq_region::UINTEGER AS seq_region, v.pos::UBIGINT AS "position", v.ref AS reference, v.alt AS alternate
 FROM read_csv('sites.vcf', delim='\t', header=false, skip=2, quote='', escape='', auto_detect=false,
   columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR'}) v
 JOIN read_parquet('regions.parquet') r ON r.seq_region_name = v.chrom;
CREATE TEMP TABLE alleles AS SELECT event_index, seq_region, "position", reference, alternate, NULL::UBIGINT AS end_position,
  NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
 FROM ev ORDER BY event_index;
COPY (SELECT e.vcf_id, t.transcript_stable_id, a.*
      FROM query(duckvep_annotate_sql('alleles', 'g', struct_pack(rich := true))) a
      JOIN ev e USING (event_index) LEFT JOIN read_parquet('transcript_ids.parquet') t USING (transcript_index))
 TO '$out/duckvep.parquet' (FORMAT parquet);
SQL
}
: > timings.tsv
for t in "$threads" 1; do
  n=$runs; [[ $t == 1 ]] && n=1
  for run in $(seq 1 "$n"); do
    /usr/bin/time -f "vep-rs\t$t\t$run\t%e\t%U\t%S\t%M" -a -o timings.tsv "$VEP_RS_DIR/vep" --json_cache "$VEP_RS_CACHE" \
      --assembly GRCh38 -i sites.vcf -o "$out/veprs.tab" --tab --fork "$t" --force_overwrite > /dev/null 2>&1
    duckvep_sql "$t" > run.sql
    /usr/bin/time -f "duckvep\t$t\t$run\t%e\t%U\t%S\t%M" -a -o timings.tsv duckdb -unsigned < run.sql > /dev/null
  done
done
echo "tool	threads	run	wall_s	user_s	sys_s	maxrss_kib"; cat timings.tsv

# ---- Tuples: (ID, transcript, sorted consequence set), on the transcripts both tools emit ------------------
tab_columns="columns={'c00':'VARCHAR','c01':'VARCHAR','c02':'VARCHAR','c03':'VARCHAR','c04':'VARCHAR','c05':'VARCHAR','c06':'VARCHAR','c07':'VARCHAR','c08':'VARCHAR','c09':'VARCHAR','c10':'VARCHAR','c11':'VARCHAR','c12':'VARCHAR','c13':'VARCHAR','c14':'VARCHAR','c15':'VARCHAR','c16':'VARCHAR','c17':'VARCHAR','c18':'VARCHAR','c19':'VARCHAR','c20':'VARCHAR'}, null_padding=true, ignore_errors=true"
rm -f compare.duckdb
duckdb compare.duckdb <<SQL
SET threads=$threads; SET preserve_insertion_order=false;
CREATE TABLE d AS SELECT vcf_id AS id, coalesce(transcript_stable_id, '-') AS feature,
  array_to_string(list_sort(string_split_regex(consequence, '[,&]')), ',') AS cons FROM '$out/duckvep.parquet';
CREATE TABLE v AS SELECT c00 AS id, c04 AS feature, array_to_string(list_sort(string_split(c06, ',')), ',') AS cons
  FROM read_csv('$out/veprs.tab', delim='\t', header=false, comment='#', quote='', escape='', auto_detect=false, $tab_columns);
CREATE TABLE shared AS SELECT feature FROM (SELECT DISTINCT feature FROM d) JOIN (SELECT DISTINCT feature FROM v) USING (feature) WHERE feature <> '-';
CREATE TABLE only_d AS SELECT d.* FROM d SEMI JOIN shared USING (feature) ANTI JOIN v USING (id, feature, cons);
CREATE TABLE only_v AS SELECT v.* FROM v SEMI JOIN shared USING (feature) ANTI JOIN d USING (id, feature, cons);
SELECT (SELECT count(*) FROM d) AS duckvep_rows, (SELECT count(*) FROM v) AS veprs_rows,
  (SELECT count(DISTINCT feature) FROM d WHERE feature <> '-') AS duckvep_transcripts,
  (SELECT count(DISTINCT feature) FROM v WHERE feature <> '-') AS veprs_transcripts, (SELECT count(*) FROM shared) AS shared_transcripts;
SELECT (SELECT count(*) FROM d SEMI JOIN shared USING (feature)) AS duckvep_tuples,
  (SELECT count(*) FROM v SEMI JOIN shared USING (feature)) AS veprs_tuples,
  (SELECT count(*) FROM d SEMI JOIN shared USING (feature) SEMI JOIN v USING (id, feature, cons)) AS matching,
  (SELECT count(*) FROM only_d) AS only_duckvep, (SELECT count(*) FROM only_v) AS only_veprs;
COPY (SELECT line FROM (SELECT 0 o, 0::BIGINT k, '##fileformat=VCFv4.2' line UNION ALL SELECT 1, 0, '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
  UNION ALL SELECT 2, row_number() OVER (ORDER BY chrom, pos::BIGINT), chrom || chr(9) || pos || chr(9) || id || chr(9) || ref || chr(9) || alt || chr(9) || '.	.	.'
  FROM (SELECT column2 AS id, column0 AS chrom, column1 AS pos, column3 AS ref, column4 AS alt
        FROM read_csv('sites.vcf', delim='\t', header=false, skip=2, quote='', all_varchar=true)) s SEMI JOIN only_v USING (id)) ORDER BY o, k)
 TO 'disagree.vcf' (FORMAT csv, HEADER false, QUOTE '', DELIMITER '\x01');
SQL

# ---- Executable VEP on every site where the two tools disagree --------------------------------------------
PATH=$vep_prefix/bin:$PATH vep --gtf "$VEP_GTF" --fasta "$VEP_FASTA" -i disagree.vcf -o disagree.vep.tab --tab \
  --force_overwrite --no_stats --fork 8 > vep.log 2>&1
duckdb compare.duckdb <<SQL
CREATE OR REPLACE TABLE p AS SELECT c00 AS id, c04 AS feature, array_to_string(list_sort(string_split(c06, ',')), ',') AS cons
  FROM read_csv('disagree.vep.tab', delim='\t', header=false, comment='#', quote='', escape='', auto_detect=false, $tab_columns);
SELECT CASE WHEN p.cons IS NULL THEN 'VEP has no row' WHEN p.cons = d.cons THEN 'VEP agrees with DuckVEP'
            WHEN p.cons = v.cons THEN 'VEP agrees with vep-rs' ELSE 'VEP differs from both' END AS verdict, count(*) AS n
FROM only_d d JOIN only_v v USING (id, feature) LEFT JOIN p USING (id, feature) GROUP BY ALL ORDER BY n DESC;
CREATE OR REPLACE TEMP TABLE pid AS SELECT DISTINCT id FROM p;
SELECT (SELECT count(*) FROM p) AS vep_tuples,
  (SELECT count(*) FROM p SEMI JOIN d USING (id, feature, cons)) AS duckvep_matches,
  (SELECT count(*) FROM p SEMI JOIN v USING (id, feature, cons)) AS veprs_matches,
  (SELECT count(*) FROM p ANTI JOIN d USING (id, feature)) AS vep_rows_absent_from_duckvep,
  (SELECT count(*) FROM p ANTI JOIN v USING (id, feature)) AS vep_rows_absent_from_veprs,
  (SELECT count(*) FROM (SELECT * FROM d SEMI JOIN pid USING (id)) x ANTI JOIN p USING (id, feature)) AS duckvep_rows_absent_from_vep,
  (SELECT count(*) FROM (SELECT * FROM v SEMI JOIN pid USING (id)) x ANTI JOIN p USING (id, feature)) AS veprs_rows_absent_from_vep;
SQL
