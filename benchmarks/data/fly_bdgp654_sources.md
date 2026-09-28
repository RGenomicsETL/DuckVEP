# Drosophila melanogaster BDGP6.54 source pins (Ensembl 116)

| Resource | Source | SHA-256 |
| --- | --- | --- |
| Public core | `mysql://anonymous@ensembldb.ensembl.org:3306/drosophila_melanogaster_core_116_11` | Local eleven-table `core.duckdb` snapshot: `f27482b51e42d056b86bcfe86d4a8f6d288c0d6756b9675d65b7c58c90e493a8` |
| Toplevel FASTA gzip | `https://ftp.ensembl.org/pub/release-116/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz` | `f4e13ad38156784b7bc2ed2ddab0e09f0da0964acd90c9a75348862429e136c1` |
| Decompressed FASTA | Same archive | `ce9277317eca7467e8143c422ec46335f929eb1363136ba3fe5f56a88ef93879` |
| FASTA index | `samtools faidx` on decompressed FASTA | `8a28f7d34a3c8f5315abc5fb495fb4f2b0285e2e6a99c14dceab1020d0fde769` |
| Reference chunks | 2,005 contiguous chunks across 1,870 contigs; 143,726,002 bases | Staged Parquet: `ee06d208916ab8728cad1c915a06972037a37dffed7e2180cf97e8f3376679f8` |
| Indexed VEP cache | `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/drosophila_melanogaster_vep_116_BDGP6.54.tar.gz` | Archive: `59f20ab8896a11b5d63aab33dc4defdc886ff24f0c05ce637262754d8252353f`; `drosophila_melanogaster/116_BDGP6.54/info.txt`: `bc88d54faa86abbac91ffc9291301f1a29e454796dfec4772d5281f253cfa4d1` |
| Oracle image | `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` | Image digest in reference |

The core hash identifies the copied tables, not a hash published by MySQL; a fresh download must be checked against this snapshot before claiming the same receipt. The toplevel FASTA contains `mitochondrion_genome` and all 1,870 modeled regions.

```sh
mkdir -p /tmp/duckvep-fly/cache
# Download the FASTA and cache archives above as fly.fa.gz and fly-cache.tar.gz.
duckdb -bail -unsigned /tmp/duckvep-fly/core.duckdb < scripts/stage_fly_core.sql
gzip -dc /tmp/duckvep-fly/fly.fa.gz > /tmp/duckvep-fly/fly.fa
samtools faidx /tmp/duckvep-fly/fly.fa
Rscript scripts/stage_species_reference.R /tmp/duckvep-fly/fly.fa \
  /tmp/duckvep-fly/reference.duckdb /tmp/duckvep-fly/reference_chunks.parquet
tar -xzf /tmp/duckvep-fly/fly-cache.tar.gz -C /tmp/duckvep-fly/cache
```
