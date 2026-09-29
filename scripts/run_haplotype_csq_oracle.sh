#!/usr/bin/env bash
set -euo pipefail
# Usage: run_haplotype_csq_oracle.sh input.vcf[.gz] reference.fa annotation.gff3[.gz] output-prefix [s|a]
# Mode a supplies domain accounting only; strict comparisons use mode s.
if (( $# < 4 || $# > 5 )); then
    echo 'usage: run_haplotype_csq_oracle.sh input.vcf[.gz] reference.fa annotation.gff3[.gz] output-prefix [s|a]' >&2
    exit 2
fi
mode=${5:-s}
[[ $mode == s || $mode == a ]] || { echo 'phase mode must be s or a' >&2; exit 2; }
input=$(realpath "$1")
reference=$(realpath "$2")
gff=$(realpath "$3")
prefix=$4
bcftools=/usr/local/bin/bcftools
expected_bcftools=6dbd8fef51e529755a4a81544075dc6a43ff3cd3
expected_htslib=c1f35d67dd5ff1e226d94abe7b850f28c60f5910
expected_gff=08e881d96ab6385a2c31f063a018be4b2c36860b323f2724be07022deeef21ce
expected_reference=1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b
version=$($bcftools --version-only)
[[ $version == '1.23.1-70-g6dbd8fef+htslib-1.22.1' ]] || { echo "unexpected bcftools/htslib: $version" >&2; exit 1; }
[[ $(sha256sum "$gff" | cut -d' ' -f1) == "$expected_gff" ]] || { echo 'GFF3 sha256 mismatch' >&2; exit 1; }
[[ $(sha256sum "$reference" | cut -d' ' -f1) == "$expected_reference" ]] || { echo 'FASTA sha256 mismatch' >&2; exit 1; }
[[ -f $reference.fai ]] || { echo 'reference .fai required' >&2; exit 1; }
[[ ! -e $prefix.bcf && ! -e $prefix.oracle.tsv ]] || { echo 'output already exists' >&2; exit 1; }
mkdir -p "$(dirname "$prefix")"
# -p s skips unphased heterozygotes; these must be counted as phase disagreements.
# -p a (used in the historical throughput run) silently assigns their phase.
"$bcftools" csq -f "$reference" -g "$gff" -p "$mode" -n 1024 -Ob -o "$prefix.bcf" "$input" 2>"$prefix.log"
if grep -q 'Too many consequences for sample' "$prefix.log"; then
    echo 'csq consequence cap reached; oracle output is incomplete' >&2
    rm "$prefix.bcf"
    exit 1
fi
{
    printf 'field\tvalue\n'
    printf 'bcftools_commit\t%s\nhtslib_commit\t%s\n' "$expected_bcftools" "$expected_htslib"
    printf 'bcftools_version\t%s\nhtslib_version\t1.22.1\nphase_mode\t%s\nncsq\t1024\n' "$version" "$mode"
    for entry in "input:$input" "reference:$reference" "gff3:$gff" "output:$prefix.bcf"; do
        key=${entry%%:*}
        file=${entry#*:}
        printf '%s_sha256\t%s\n' "$key" "$(sha256sum "$file" | cut -d' ' -f1)"
        printf '%s_path\t%s\n' "$key" "$file"
    done
} >"$prefix.oracle.tsv"
