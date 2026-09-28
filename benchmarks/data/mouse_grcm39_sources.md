# Mus musculus GRCm39 source pins (Ensembl 116)

| Resource | Source | SHA-256 |
| --- | --- | --- |
| Public core | `mysql://anonymous@ensembldb.ensembl.org:3306/mus_musculus_core_116_39` | Snapshot `9a14c954620c71ffe388f020209088d952be8ab85c393dbc161b5f6137cb276a` (`core.duckdb`) |
| Primary-assembly FASTA gzip | `http://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz` | `c661d19cfdbbee7ffbafa9bffb44581c6306480b9fef7b70e1d9c173782d370f` |
| Decompressed FASTA | Same archive | `14571f7559e292baf0a40f9d155c41ede19a04d80fdeb59a0c2dfe566db90552` |
| FASTA index | `samtools faidx` on decompressed FASTA | `17e430fb8bba1dd9ca5c1c3aa4c9a5a5af51a3f9a565af5aec33100fbc5fba94` |
| Staged reference chunks | 2,774 contiguous 1 Mb chunks, 2,728,222,451 bases across 61 contigs | `438627a7deab2b8b5e7c22274d153c5eb31f704f1c574489ce214f91551ee588` (`reference_chunks.parquet`; byte hash is staging-run-specific) |
| Indexed VEP cache | `http://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/mus_musculus_vep_116_GRCm39.tar.gz` | `d82d158f22cbb99de5ae82780f8ad83b8fb5081631a204897652b7d91b8bdcf7` |
| Oracle image | `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` | Image digest in reference |

The core hash pins the local snapshot of the eleven tables copied by `scripts/stage_mouse_core.sql`, not a hash asserted by the public MySQL server. The FASTA contains 61 contigs including MT; every indexed contig and length resolves to a release-116 core sequence region. The VEP cache contains `mus_musculus/116_GRCm39/info.txt`. The model builder loads DuckVEP alone in its DuckDB connection.

```sh
mkdir -p /tmp/duckvep-mouse
# Download the FASTA and cache from the URLs above into this directory.
duckdb -bail -unsigned /tmp/duckvep-mouse/core.duckdb < scripts/stage_mouse_core.sql
gzip -dc /tmp/duckvep-mouse/mouse.fa.gz > /tmp/duckvep-mouse/mouse.fa
samtools faidx /tmp/duckvep-mouse/mouse.fa
Rscript scripts/stage_species_reference.R /tmp/duckvep-mouse/mouse.fa /tmp/duckvep-mouse/reference.duckdb /tmp/duckvep-mouse/reference_chunks.parquet
mkdir -p /tmp/duckvep-mouse/cache
tar -xzf /tmp/duckvep-mouse/mouse-cache.tar.gz -C /tmp/duckvep-mouse/cache
```

The public MySQL endpoint can change its contents without changing its database name; verify the snapshot digest before claiming the model receipt below is the same build input.
