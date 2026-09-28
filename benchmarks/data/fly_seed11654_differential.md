# Drosophila melanogaster BDGP6.54 VEP 116 differential

The pinned Docker VEP 116 image `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` uses the [indexed cache and FASTA](fly_bdgp654_sources.md) with 5,000-bp distance. The [seed-11654 corpus](fly_seed11654_corpus.md) yields 3,527 JSON records (SHA-256 `d8c7a0bea9f4d387d621c9a8271243c32458abe2eae27871c2cda58c1b483b16`). The comparator independently normalizes consequence-term sets for each `(variant ID, transcript ID)` and full-outer-joins VEP with the [DuckVEP model](fly_bdgp654_model.md) SHA-256 `e6deea1ac2b0097589df9da501bd4a1ddf47ec5f8e0d4291ee49e727352533c8`.

```sh
scripts/run_species_vep116_docker.sh drosophila_melanogaster BDGP6.54 116 \
  /tmp/duckvep-fly/cache /tmp/duckvep-fly/fly.fa \
  /tmp/duckvep-fly/fly-seed11654.vcf /tmp/duckvep-fly/oracle-seed11654.json
Rscript scripts/compare_species_corpus.R \
  /tmp/duckvep-fly/fly.duckdb /tmp/duckvep-fly/fly-seed11654.vcf \
  /tmp/duckvep-fly/fly-seed11654.vcf.provenance.tsv \
  /tmp/duckvep-fly/oracle-seed11654.json \
  /tmp/duckvep-fly/differential-seed11654 drosophila_melanogaster
```

| Stratification | Exact | Oracle-only | DuckVEP-only | Terms differ |
| --- | ---: | ---: | ---: | ---: |
| **All transcript pairs** | **36,136** | **0** | **0** | **0** |
| Nuclear | 33,392 | 0 | 0 | 0 |
| Mitochondrial | 2,744 | 0 | 0 | 0 |
| Codon table 1 | 17,279 | 0 | 0 | 0 |
| Codon table 5 | 1,040 | 0 | 0 | 0 |
| No codon table | 17,817 | 0 | 0 | 0 |
| Forward strand | 17,910 | 0 | 0 | 0 |
| Reverse strand | 18,226 | 0 | 0 | 0 |
| SNV | 9,142 | 0 | 0 | 0 |
| MNV | 8,998 | 0 | 0 | 0 |
| Insertion | 8,998 | 0 | 0 | 0 |
| Deletion | 8,998 | 0 | 0 | 0 |

All ten biotypes have exact pairs; [stratified tables](fly_seed11654_differential/) include biotype and SO-term counts. Every DuckVEP pair has `supported` status, and the [disagreement table](fly_seed11654_differential/disagreements.csv) has no data rows. All 12 designated table-1 TGG→TGA witnesses are `stop_gained`, while all three designated table-5 mitochondrial TGG→TGA witnesses are `synonymous_variant` in both engines.

**Admission: yes.** The `drosophila_melanogaster` [current-definition receipt](duckvep_model_receipts.csv) is admitted on this exact seeded differential. The comparison covers these model-derived variants and transcript consequences, not every possible allele or regulatory/HGVS annotation.
