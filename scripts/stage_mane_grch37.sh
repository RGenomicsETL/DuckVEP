#!/bin/sh
# Stage immutable MANE and GRCh37.p13 inputs outside the repository.
set -eu
if [ "$#" -ne 1 ]; then
  echo "usage: scripts/stage_mane_grch37.sh STAGING_DIR" >&2
  exit 2
fi
stage=$1
mkdir -p "$stage"
ncbi=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13
mane=https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5
ensembl=https://ftp.ensembl.org/pub/grch37/release-116/fasta/homo_sapiens/dna
fetch() {
  expected=$1
  url=$2
  name=${url##*/}
  file=$stage/$name
  if [ -f "$file" ] && printf '%s  %s\n' "$expected" "$file" | sha256sum --check --status; then
    return
  fi
  curl --fail --location --retry 5 --continue-at - --output "$file" "$url"
  printf '%s  %s\n' "$expected" "$file" | sha256sum --check
}
fetch d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c "$mane/MANE.GRCh38.v1.5.summary.txt.gz"
fetch a6cf8300aa2cef9188590bad2d9d54a5909f6d4d2da3f22aa5c5ba2fda1adab3 "$ncbi/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
fetch 5fcadac26be5d82a1f1c52e33cc5047247f719a60b50182c6ba298bda77cc80f "$ncbi/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
fetch 0f48eeebe6e5631ff8b346e41db8a9e0d7a3506e1acc9df3fcdcdac8139beb2b "$ncbi/GCF_000001405.25_GRCh37.p13_rna.fna.gz"
fetch 7697d816c1a1639ed304e0fc0fff6867d55df91a5e50f6b6f7287055f389f514 "$ncbi/GCF_000001405.25_GRCh37.p13_protein.faa.gz"
fetch 0a43b56dec40debae976d6e70cac68ea6ed874f9fb7c8c814363702ff1d47865 "$ensembl/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
fasta=$stage/Homo_sapiens.GRCh37.dna.primary_assembly.fa
if [ ! -f "$fasta" ]; then
  gzip -dc "$fasta.gz" > "$fasta"
fi
printf '%s  %s\n' 3a3872e7bdd1532fdbfbc3afd70c04ff8a12a109319ae050c91d962f32e16a4f "$fasta" | sha256sum --check
if [ ! -f "$fasta.fai" ]; then
  samtools faidx "$fasta"
fi
