# Arabidopsis thaliana TAIR10 VEP 116 differential

The pinned Docker VEP 116 image `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` uses the [Plants 63 indexed cache and FASTA](plant_tair10_sources.md) with 5,000-bp distance. The [seed-11663 corpus](plant_seed11663_corpus.md) yields 2,667 JSON records (SHA-256 `c71375512e6fb66b43e13843008eba0580e49a617230f48e0c7d5023edd1051e`). The comparator loads the 325 mature-miRNA segments from the [model](plant_tair10_model.md) SHA-256 `5b3fac7c92db9c2dfd13cc83df457a993ece29b5b703343194c5363136d0a12f`, normalizes consequence-term sets for each `(variant ID, transcript ID)` and full-outer-joins VEP and DuckVEP.

```sh
scripts/run_species_vep116_docker.sh arabidopsis_thaliana TAIR10 63 \
  /tmp/duckvep-plant/cache /tmp/duckvep-plant/plant.fa \
  /tmp/duckvep-plant/plant-seed11663.vcf /tmp/duckvep-plant/oracle-seed11663.json
Rscript scripts/compare_species_corpus.R \
  /tmp/duckvep-plant/plant.duckdb /tmp/duckvep-plant/plant-seed11663.vcf \
  /tmp/duckvep-plant/plant-seed11663.vcf.provenance.tsv \
  /tmp/duckvep-plant/oracle-seed11663.json \
  /tmp/duckvep-plant/differential-seed11663 arabidopsis_thaliana
```

| Stratification | Exact | Oracle-only | DuckVEP-only | Terms differ |
| --- | ---: | ---: | ---: | ---: |
| **All transcript pairs** | **21,772** | **0** | **0** | **0** |
| Nuclear | 19,148 | 0 | 0 | 0 |
| Mitochondrial (`Mt`) | 804 | 0 | 0 | 0 |
| Chloroplast (`Pt`) | 1,820 | 0 | 0 | 0 |
| Codon table 1 | 14,180 | 0 | 0 | 0 |
| Codon table 11 | 1,098 | 0 | 0 | 0 |
| No codon table | 6,494 | 0 | 0 | 0 |
| Forward strand | 10,638 | 0 | 0 | 0 |
| Reverse strand | 11,134 | 0 | 0 | 0 |
| SNV | 5,716 | 0 | 0 | 0 |
| MNV | 5,352 | 0 | 0 | 0 |
| Insertion | 5,352 | 0 | 0 | 0 |
| Deletion | 5,352 | 0 | 0 | 0 |

All nine biotypes have exact pairs; [stratified tables](plant_seed11663_differential/) include biotype and SO-term counts. All DuckVEP pairs have `supported` status; the disagreement table has only its header. Four designated table-11 GTG→GCG initiation witnesses and 26 ATG→GTG witnesses match VEP `start_lost`; nine table-1 TGG→TGA witnesses match `stop_gained`.

**Admission: yes.** The `arabidopsis_thaliana` [current-definition receipt](duckvep_model_receipts.csv) is admitted on this exact seeded differential. The comparison covers these model-derived variants and transcript consequences, not every possible allele or regulatory/HGVS annotation.
