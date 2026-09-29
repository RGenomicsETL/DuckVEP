# Arabidopsis thaliana TAIR10 source pins (Ensembl Plants 63 / VEP 116)

| Resource | Source | SHA-256 |
| --- | --- | --- |
| Public core | `mysql://anonymous@mysql-eg-publicsql.ebi.ac.uk:4157/arabidopsis_thaliana_core_63_116_11` | Local eleven-table `core.duckdb` snapshot: `7e804577dc733a44f8cbb805f171c5a5e212606d07ce31a20d1768cc34486b2b` |
| Toplevel FASTA gzip | `https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz` | `80a13166003333ba4982b57abdb3d2b62037bc831f4da29ea378c1b776fff8e2` |
| Decompressed FASTA | Same archive | `85c83b6dd6820769dae5df190adda7ca0762872d61e5634a91f1d14c8f4fd8ce` |
| FASTA index | `samtools faidx` on decompressed FASTA | `cc29e77ff371b947f8d8f1045769c5026cb35b20119c4128112fbe8d5859c28e` |
| Reference chunks | 123 contiguous chunks across 7 contigs; 119,667,750 bases | Staged Parquet: `f46a71495ab8c6ceb6a94b8773c80c098c3535fac6a94b60bc90bed30348d5ae` |
| Indexed VEP cache | `https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/variation/indexed_vep_cache/arabidopsis_thaliana_vep_63_TAIR10.tar.gz` | Archive: `a62b4d9dd5e70fc0342a3ce1756cb927edad4d5440a2292da634363b24f729b2`; `arabidopsis_thaliana/63_TAIR10/info.txt`: `7f740ccf1055d2f3a356dfa80cdaa7cc523c70b18eca15830496bf36ae37b322` |
| Oracle image | `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` | Image digest in reference |

The snapshot hash identifies the copied MySQL tables, not a published upstream hash; verify a new snapshot against it before attributing the same receipt. The toplevel FASTA includes chromosomes `1`–`5`, `Mt` (366,924 bases) and `Pt` (154,478 bases). The pinned core annotates `Mt` with codon table 1 and `Pt` with table 11.

```sh
mkdir -p /tmp/duckvep-plant/cache
# Download the FASTA and cache archives above as plant.fa.gz and plant-cache.tar.gz.
duckdb -bail -unsigned /tmp/duckvep-plant/core.duckdb < scripts/stage_plant_core.sql
gzip -dc /tmp/duckvep-plant/plant.fa.gz > /tmp/duckvep-plant/plant.fa
samtools faidx /tmp/duckvep-plant/plant.fa
Rscript scripts/stage_species_reference.R /tmp/duckvep-plant/plant.fa \
  /tmp/duckvep-plant/reference.duckdb /tmp/duckvep-plant/reference_chunks.parquet
tar -xzf /tmp/duckvep-plant/plant-cache.tar.gz -C /tmp/duckvep-plant/cache
```
