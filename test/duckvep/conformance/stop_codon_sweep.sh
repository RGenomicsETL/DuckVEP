#!/usr/bin/env bash
# Compare DuckVEP with executable VEP on every 1 to 4 bp deletion around each
# stop codon of one chromosome (cache path, --hgvs). Consequence sets, HGVSc and
# HGVSp are compared separately for every (allele, transcript) pair.
#
# Usage: stop_codon_sweep.sh CHROM GTF FASTA VEP_CACHE_DIR MODEL_DB OUT_DIR
# Environment: VEP (default vep), DUCKDB (default duckdb), SAMTOOLS, FORK (8),
#   CACHE_VERSION (116), ASSEMBLY (GRCh38), DUCKVEP_EXTENSION.
# OUT_DIR must not exist. Exit status 1 means at least one pair differs.
set -euo pipefail
[[ $# == 6 ]] || { sed -n 2,9p "$0" >&2; exit 2; }
chrom=$1; gtf=$2; fasta=$3; cache=$4; model=$5; out=$6
root=$(cd "$(dirname "$0")/../../.." && pwd)
vep=${VEP:-vep}; duckdb_bin=${DUCKDB:-duckdb}; samtools=${SAMTOOLS:-samtools}
extension=${DUCKVEP_EXTENSION:-$root/build/release/duckvep.duckdb_extension}
[[ ! -e $out ]] || { echo "OUT_DIR exists: $out" >&2; exit 2; }
mkdir -p "$out"

python3 - "$chrom" "$gtf" "$fasta" "$samtools" "$out/events.tsv" <<'PY'
import gzip, subprocess, sys
chrom, gtf, fasta, samtools, out = sys.argv[1:]
stops = set()
opener = gzip.open if gtf.endswith(".gz") else open
with opener(gtf, "rt") as handle:
    for line in handle:
        if not line.startswith(chrom + "\t"):
            continue
        fields = line.split("\t")
        if fields[2] == "stop_codon":
            stops.add((int(fields[3]), int(fields[4])))
rows = set()
for start, end in sorted(stops):
    low, high = start - 9, end + 12
    text = subprocess.run([samtools, "faidx", fasta, f"{chrom}:{low}-{high}"],
                          check=True, capture_output=True, text=True).stdout
    sequence = "".join(text.split("\n")[1:]).upper()
    for position in range(low, end + 7):
        for length in (1, 2, 3, 4):
            reference = sequence[position - low:position - low + length + 1]
            if len(reference) == length + 1 and "N" not in reference:
                rows.add((position, reference))
with open(out, "w") as handle:
    for position, reference in sorted(rows):
        handle.write(f"{chrom}\t{position}\ts{position}_{len(reference) - 1}\t{reference}\t{reference[0]}\t.\t.\t.\n")
print(f"{len(stops)} stop codons, {len(rows)} alleles", file=sys.stderr)
PY
{ printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'; cat "$out/events.tsv"; } > "$out/events.vcf"

"$vep" -i "$out/events.vcf" --cache --offline --dir_cache "$cache" \
  --cache_version "${CACHE_VERSION:-116}" --assembly "${ASSEMBLY:-GRCh38}" --fasta "$fasta" \
  --hgvs --numbers --distance 5000 --tab -o "$out/vep.raw.tsv" --no_stats --fork "${FORK:-8}"
grep -v '^##' "$out/vep.raw.tsv" > "$out/vep.tsv"

"$duckdb_bin" -unsigned "$out/sweep.duckdb" <<SQL
LOAD '$extension';
ATTACH '$model' AS src (READ_ONLY);
SELECT loaded FROM duckvep_model_load('sweep',
 'SELECT seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name FROM src.model_regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM src.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM src.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM src.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid, x.edit_code FROM src.model_transcripts, unnest(peptide_edits) u(x) ORDER BY 1, 2',
 transcript_coverage_complete := TRUE, reference_fasta := '$fasta');
CREATE TABLE events AS
SELECT row_number() OVER (ORDER BY pos, ref)::UBIGINT AS event_index,
       (SELECT seq_region FROM src.model_regions WHERE seq_region_name = '$chrom')::UINTEGER AS seq_region,
       pos::UBIGINT AS position, ref AS reference, alt AS alternate,
       NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
       NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position, id
FROM read_csv('$out/events.tsv', delim = '\t', header = false, quote = '',
  columns = {'chrom': 'VARCHAR', 'pos': 'BIGINT', 'id': 'VARCHAR', 'ref': 'VARCHAR', 'alt': 'VARCHAR', 'q': 'VARCHAR', 'f': 'VARCHAR', 'i': 'VARCHAR'});
CREATE TABLE duckvep AS
SELECT e.id, t.transcript_stable_id AS transcript, a.consequence, a.transcript_hgvs AS hgvsc, a.protein_hgvs AS hgvsp
FROM query(duckvep_annotate_sql('events', 'sweep', struct_pack(rich := true, hgvs := true))) a
JOIN events e USING (event_index) JOIN src.model_transcripts t USING (transcript_index);
CREATE TABLE vep AS
SELECT "#Uploaded_variation" AS id, Feature AS transcript, Consequence AS consequence,
       nullif(regexp_replace(HGVSc, '^[^:]+:', ''), '-') AS hgvsc,
       nullif(replace(regexp_replace(HGVSp, '^[^:]+:', ''), '%3D', '='), '-') AS hgvsp
FROM read_csv('$out/vep.tsv', delim = '\t', header = true, all_varchar = true, quote = '')
WHERE Feature_type = 'Transcript';
CREATE TABLE pairs AS
SELECT id, transcript, v.consequence AS vep_consequence, d.consequence AS duckvep_consequence,
       v.hgvsc AS vep_hgvsc, d.hgvsc AS duckvep_hgvsc, v.hgvsp AS vep_hgvsp, d.hgvsp AS duckvep_hgvsp,
       v.id IS NULL AS duckvep_only, d.id IS NULL AS vep_only,
       list_sort(string_split(v.consequence, ',')) IS DISTINCT FROM list_sort(string_split(d.consequence, '&')) AS consequence_differs,
       v.hgvsc IS DISTINCT FROM d.hgvsc AS hgvsc_differs, v.hgvsp IS DISTINCT FROM d.hgvsp AS hgvsp_differs
FROM vep v FULL JOIN duckvep d USING (id, transcript);
COPY (SELECT count(*) AS pairs, count(*) FILTER (WHERE vep_only) AS vep_only,
        count(*) FILTER (WHERE duckvep_only) AS duckvep_only,
        count(*) FILTER (WHERE consequence_differs) AS consequence_differences,
        count(*) FILTER (WHERE hgvsc_differs) AS hgvsc_differences,
        count(*) FILTER (WHERE hgvsp_differs) AS hgvsp_differences,
        count(*) FILTER (WHERE vep_hgvsc IS NOT NULL) AS vep_hgvsc, count(*) FILTER (WHERE vep_hgvsp IS NOT NULL) AS vep_hgvsp
      FROM pairs) TO '$out/summary.csv' (HEADER);
COPY (SELECT * FROM pairs WHERE consequence_differs OR hgvsc_differs OR hgvsp_differs ORDER BY id, transcript)
  TO '$out/differences.tsv' (HEADER, DELIMITER '\t');
SQL
cat "$out/summary.csv"
[[ $(wc -l < "$out/differences.tsv") == 1 ]]
