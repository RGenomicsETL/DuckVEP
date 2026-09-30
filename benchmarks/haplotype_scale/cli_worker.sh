#!/usr/bin/env bash
# One fresh DuckDB CLI process of the #34 slice 1 gate: model load, stage the coding calls with the fused native reader
# (duckvep_coding_calls) or with mode B's SQL (read_csv + duckvep_coding_transcripts), then duckvep_haplotypes into Parquet.
# No R: the gate does not depend on R start-up (report that separately, worker.R).
#
#   cli_worker.sh --extension EXT --mode F|I|B --vcf VCF --out OUT.parquet [--model-file M.duckdb] [--native-budget BYTES]
#                 [--spill DIR] [--audit]
#
# Modes   F  fused reader staged into a calls table, then duckvep_haplotypes over the table (mode B's structure).
#         I  inline: duckvep_haplotypes over duckvep_coding_calls directly (nothing staged).
#         B  mode B's SQL (DuckDB read_csv, scalar discovery, SQL genotype parsing), the #2 slice 7 reference.
# Prints one line per stage, "STAGE <name> <seconds>" (the CLI's own wall timer), and with --audit the output's
# row count and full-output checksum (sum(hash(row)) over the Parquet, as in worker.R), "AUDIT rows=.. checksum=..".
set -euo pipefail
ext=; mode=F; vcf=; out=; model=/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb
budget=4294967296; spill=${TMPDIR:-/tmp}/duckvep-cli-spill; audit=0
while [[ $# -gt 0 ]]; do case $1 in
    --extension) ext=$2; shift 2;; --mode) mode=$2; shift 2;; --vcf) vcf=$2; shift 2;; --out) out=$2; shift 2;;
    --model-file) model=$2; shift 2;; --native-budget) budget=$2; shift 2;; --spill) spill=$2; shift 2;;
    --audit) audit=1; shift;; *) echo "unknown option $1" >&2; exit 2;; esac; done
[[ -n $ext && -n $vcf && -n $out ]] || { echo "need --extension, --vcf and --out" >&2; exit 2; }
mkdir -p "$spill"; rm -f "$out" "$out.partial"
duckdb_bin=${DUCKDB:-duckdb}
skip=0; [[ $mode == B ]] && skip=$(zcat "$vcf" 2>/dev/null | head -n 5000 | grep -c '^#' || true)
sql=$(mktemp "$spill/cli-XXXXXX.sql")
{
cat <<SQL
.timer on
SET threads = 1;
SET memory_limit = '8GB';
SET temp_directory = '$spill';
SET max_temp_directory_size = '32GiB';
LOAD '$ext';
ATTACH '$model' AS m (READ_ONLY);
SELECT duckvep_native_budget_set($budget);
.print STAGE-BEGIN load
SELECT loaded FROM duckvep_model_load('hap', 'SELECT seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name FROM m.model_regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM m.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM m.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM m.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid FROM m.model_transcripts, unnest(peptide_edits) u(x)',
 transcript_coverage_complete := TRUE);
.print STAGE-END load
.print STAGE-BEGIN stage
SQL
case $mode in
F) echo "CREATE TABLE calls AS SELECT * FROM duckvep_coding_calls('hap', '$vcf');" ;;
I) echo "SELECT 1;" ;;
B) cat <<SQL
CREATE TABLE regions AS SELECT seq_region::BIGINT AS seq_region, seq_region_name, sequence_length FROM m.model_regions;
CREATE TABLE cand AS SELECT * FROM (SELECT v.record_index, a.i, r.seq_region, v.pos, v.ref, a.alt AS alt, v.sample,
 duckvep_coding_transcripts('hap', r.seq_region, v.pos, v.ref, a.alt) AS tx
 FROM (SELECT row_number() OVER () AS record_index, chrom, pos, ref, alt, sample FROM read_csv('$vcf', delim='\t', header=false, skip=$skip, auto_detect=false, quote='', escape='', strict_mode=false,
  columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR','fmt':'VARCHAR','sample':'VARCHAR'}, compression='gzip')) v
 JOIN regions r ON r.seq_region_name = v.chrom, unnest(string_split(v.alt, ',')) WITH ORDINALITY a(alt, i) WHERE v.alt <> '.') WHERE len(tx) > 0;
CREATE TABLE calls AS SELECT ((c.record_index << 6) | (c.i - 1))::BIGINT AS event_index, c.seq_region::INTEGER AS seq_region,
 c.pos::BIGINT AS position, c.ref AS reference, c.alt AS alternate, c.i::INTEGER AS alt_index, t.transcript_index::INTEGER AS transcript_index,
 0::INTEGER AS sample_index, list_transform(string_split_regex(split_part(c.sample, ':', 1), '[/|]'), lambda x: try_cast(x AS INTEGER)) AS alleles,
 [false] || list_transform(regexp_extract_all(split_part(c.sample, ':', 1), '[/|]'), lambda s: s = '|') AS phase_before,
 try_cast(list_last(string_split(c.sample, ':')) AS BIGINT) AS phase_set
 FROM cand c, unnest(c.tx) t(transcript_index)
 WHERE CASE WHEN c.i > 64 THEN error('more than 64 ALT alleles in one record') ELSE true END
 ORDER BY seq_region, position, event_index, transcript_index;
SQL
;;
esac
echo ".print STAGE-END stage"
echo ".print STAGE-BEGIN predict"
if [[ $mode == I ]]; then src="SELECT * FROM duckvep_coding_calls(''hap'', ''$vcf'')"; else src="SELECT * FROM calls"; fi
cat <<SQL
COPY (SELECT * FROM duckvep_haplotypes('$src', 'hap', max_alignment_cells := 268435456, workspace_limit := 1073741824)) TO '$out.partial' (FORMAT parquet);
.print STAGE-END predict
SQL
if [[ $audit == 1 ]]; then
cat <<SQL
.timer off
SELECT 'AUDIT rows=' || count(*) || ' checksum=' || sum(hash(t)::HUGEINT)::VARCHAR FROM read_parquet('$out.partial') t;
SELECT 'AUDITCALLS calls=' || count(*) FROM (SELECT * FROM duckvep_coding_calls('hap', '$vcf'));
SQL
fi
} > "$sql"
start=$(date +%s.%N)
"$duckdb_bin" -unsigned -csv -noheader -f "$sql" 2>&1 | awk '
  /^STAGE-BEGIN/ { name = $2; next }
  /^STAGE-END/ { printf "STAGE %s %.3f\n", name, sum[name]; name = ""; next }
  /^Run Time/ { if (name != "") sum[name] += $5; next }
  /^AUDIT/ || /Error/ { print; next }'
end=$(date +%s.%N)
mv "$out.partial" "$out"
printf 'WALL %.3f\n' "$(echo "$end - $start" | bc)"
rm -f "$sql"
