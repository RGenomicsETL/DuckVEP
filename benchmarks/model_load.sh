#!/usr/bin/env bash
# Model load from relations versus snapshot restore, fresh DuckDB CLI processes.
#
#   bash benchmarks/model_load.sh [MODEL.duckdb] [RUNS]
#
# Each run is one process: LOAD the extension, then either duckvep_model_load from the model database
# (threads 1 and 8) or duckvep_model_restore from a snapshot saved once beforehand. Appends one row per run
# to benchmarks/data/duckvep_model_load.csv: date, revision, extension SHA-256, model SHA-256, method,
# threads, seconds (the CLI's timer around the load statement), peak RSS (KiB, GNU time) and CPU.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
model=${1:-/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb}; runs=${2:-3}
ext=$root/build/release/duckvep.duckdb_extension
out=$root/benchmarks/data/duckvep_model_load.csv
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$ext" "$work/duckvep.duckdb_extension"; ext=$work/duckvep.duckdb_extension
load="SELECT loaded FROM duckvep_model_load('m',
 'SELECT seq_region::UINTEGER AS seq_region, sequence_length, seq_region_name FROM src.model_regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM src.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM src.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM src.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid, x.edit_code FROM src.model_transcripts, unnest(peptide_edits) u(x) ORDER BY 1, 2',
 transcript_coverage_complete := TRUE);"
duckdb -unsigned :memory: -c "LOAD '$ext'; ATTACH '$model' AS src (READ_ONLY); $load SELECT duckvep_model_save('m', '$work/model.dvsnap');" > /dev/null
[ -s "$out" ] || echo "run_date,source_revision,extension_sha256,model_sha256,method,threads,seconds,max_rss_kib,cpu" > "$out"
rev=$(git -C "$root" rev-parse HEAD); esha=$(sha256sum "$ext" | cut -d' ' -f1); msha=$(sha256sum "$model" | cut -d' ' -f1)
cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')
one() { # method threads statement
    local seconds
    seconds=$(/usr/bin/time -f "%M" -o "$work/rss" duckdb -unsigned :memory: -c "SET threads = $2; LOAD '$ext'; ATTACH '$model' AS src (READ_ONLY);" \
        -c ".timer on" -c "$3" 2>&1 | awk '/^Run Time/ {s = $5} END {print s}')
    printf '%s,%s,%s,%s,%s,%s,%s,%s,"%s"\n' "$(date +%F)" "$rev" "$esha" "$msha" "$1" "$2" "$seconds" "$(cat "$work/rss")" "$cpu" >> "$out"
}
for run in $(seq "$runs"); do
    one relation_load 1 "$load"
    one relation_load 8 "$load"
    one snapshot_restore 1 "SELECT duckvep_model_restore('m', '$work/model.dvsnap');"
done
