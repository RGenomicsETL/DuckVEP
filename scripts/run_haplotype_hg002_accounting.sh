#!/usr/bin/env bash
set -euo pipefail
# Usage: run_haplotype_hg002_accounting.sh phased.vcf.gz reference.fa gff3.gz model.duckdb output-prefix
if (( $# != 5 )); then
    echo 'usage: run_haplotype_hg002_accounting.sh phased.vcf.gz reference.fa gff3.gz model.duckdb output-prefix' >&2
    exit 2
fi
input=$1; reference=$2; gff=$3; model=$4; prefix=$5
[[ -f $model && ! -e $prefix.counts.tsv ]] || { echo 'missing model or existing counts' >&2; exit 1; }
[[ $(/usr/local/bin/bcftools query -l "$input" | wc -l) -eq 1 ]] || {
    echo 'one sample required for record-level phase accounting' >&2; exit 1;
}
scripts/run_haplotype_csq_oracle.sh "$input" "$reference" "$gff" "$prefix.strict" s
scripts/run_haplotype_csq_oracle.sh "$input" "$reference" "$gff" "$prefix.domain" a
# Extract only the GFF fields used in the exon/CDS geometry join.
gzip -cd "$gff" | awk -F '\t' 'BEGIN {OFS="\t"} $3=="mRNA" || $3=="exon" || $3=="CDS" {
  id=$9; if ($3=="mRNA") sub(/^.*ID=transcript:/,"",id); else sub(/^.*Parent=transcript:/,"",id)
  sub(/;.*/,"",id)
  if (id ~ /^ENST[0-9]+$/) print id,$1,$3,$4,$5,$7,$8
}' > "$prefix.geometry.tsv"
/usr/local/bin/bcftools query -f '%INFO/BCSQ\n' "$prefix.domain.bcf" |
    awk -F '[,|]' '{for (i=1;i<=NF-3;i++) if ($(i+3)=="protein_coding" && $(i+2) ~ /^ENST[0-9]+$/) print $(i+2)}' |
    sort -u > "$prefix.coding-transcripts.txt"
Rscript scripts/haplotype_geometry_parity.R "$model" "$prefix.geometry.tsv" \
    "$prefix.coding-transcripts.txt" "$prefix.parity.tsv"
for mode in strict domain; do
    /usr/local/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\t[%GT]\t[%PS]\t[%TBCSQ]\n' \
        "$prefix.$mode.bcf" > "$prefix.$mode.tsv"
done
Rscript scripts/account_haplotype_csq.R "$prefix.strict.tsv" "$prefix.domain.tsv" \
    "$prefix.parity.tsv" "$prefix.counts.tsv"
{
    printf 'field\tvalue\ncomparison_scope\tduckvep-coding\n'
    printf 'domain_annotated_records\t%s\n' "$(awk -F '\t' '$5!="."{n++} END{print n+0}' "$prefix.domain.tsv")"
    printf 'domain_unannotated_records\t%s\n' "$(awk -F '\t' '$5=="."{n++} END{print n+0}' "$prefix.domain.tsv")"
    printf 'model_sha256\t%s\n' "$(sha256sum "$model" | cut -d' ' -f1)"
    printf 'parity_sha256\t%s\n' "$(sha256sum "$prefix.parity.tsv" | cut -d' ' -f1)"
    printf 'strict_stream_sha256\t%s\n' "$(sha256sum "$prefix.strict.tsv" | cut -d' ' -f1)"
    printf 'domain_stream_sha256\t%s\n' "$(sha256sum "$prefix.domain.tsv" | cut -d' ' -f1)"
    printf 'counts_sha256\t%s\n' "$(sha256sum "$prefix.counts.tsv" | cut -d' ' -f1)"
    printf 'record_ledger_content_sha256\t%s\n' "$(gzip -cd "$prefix.counts.tsv.records.tsv.gz" | sha256sum | cut -d' ' -f1)"
} > "$prefix.accounting.tsv"
rm "$prefix.geometry.tsv" "$prefix.coding-transcripts.txt" "$prefix.strict.tsv" "$prefix.domain.tsv"
