# Tetrahymena thermophila JCVI-TTA1-2.2 VEP 116 differential

The pinned Docker image `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` annotated the [seed-11606 corpus](tetrahymena_seed11606_corpus.md) against the [Protists 63 indexed cache](tetrahymena_jcvi_tta1_sources.md). Oracle JSON SHA-256: `a6de6345e2202817e07b250e23c2539365a2876c0ed53d46b3319c2712df9e24`. The comparator normalizes consequence-term sets by `(variant ID, transcript ID)` and full-outer-joins the VEP and [model](tetrahymena_jcvi_tta1_model.md) outputs.

| Comparison | VEP pairs | DuckVEP pairs | Exact | Disagreements |
| --- | ---: | ---: | ---: | ---: |
| Seed 11606 | 100,201 | 100,201 | 100,201 | 0 |

All 12 codon witnesses have exact synonymous consequences in table 6. Both TAA and TAG glutamine witnesses differ from their table-1 stop-lost interpretations in the committed offline SQL and R tests. The comparator outputs `verdicts.csv` (SHA-256 `2d6cb573b84f53f23780cd16674f59e6e865d80aebbe7f8a75687aa457f106ab`) and `codon_witness_pairs.csv` (SHA-256 `f7dc9d21f4a9383b856847f7e4d930b4caefebbca0958421aeb48fd5bab7318d`).

```sh
scripts/run_species_vep116_docker.sh tetrahymena_thermophila JCVI-TTA1-2.2 63 \
  /tmp/duckvep-tetrahymena-species-final/cache \
  /tmp/duckvep-tetrahymena-species-final/tetra.fa \
  /tmp/duckvep-tetrahymena-species-final/tetra-seed11606.vcf \
  /tmp/duckvep-tetrahymena-species-final/oracle-seed11606.json
DUCKVEP_EXTENSION_FILE=/path/to/immutable/duckvep.duckdb_extension \
  Rscript scripts/compare_species_corpus.R MODEL_DB CORPUS.vcf CORPUS.vcf.provenance.tsv \
  ORACLE.json DIFFERENTIAL_DIR tetrahymena_thermophila
```
