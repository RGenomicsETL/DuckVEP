# Mus musculus GRCm39 VEP 116 differential

The oracle runs `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` with the [pinned indexed cache and FASTA](mouse_grcm39_sources.md). It consumes the [seed-11639 VCF](mouse_seed11639_corpus.md) and produces 22,204 JSON records (SHA-256 `bd74474a730ffd7aa157d14929a83c9b0b6f5d9c67824e4d89c8fd93259abf4e`). The generalized comparator matches normalized consequence-term **sets for each variant/transcript pair** to the DuckVEP model SHA-256 `ae39ffc9e647d0a096a13b737599d480fae2b94d3a202938a8eeaa44208c62be`.

```sh
scripts/run_species_vep116_docker.sh mus_musculus GRCm39 116 \
  /tmp/duckvep-mouse/cache /tmp/duckvep-mouse/mouse.fa \
  /tmp/duckvep-mouse/mouse-seed11639.vcf /tmp/duckvep-mouse/oracle-seed11639.json
Rscript scripts/compare_species_corpus.R \
  /tmp/duckvep-mouse/mouse.duckdb /tmp/duckvep-mouse/mouse-seed11639.vcf \
  /tmp/duckvep-mouse/mouse-seed11639.vcf.provenance.tsv \
  /tmp/duckvep-mouse/oracle-seed11639.json \
  /tmp/duckvep-mouse/differential-seed11639 mus_musculus
```

| Stratification | Exact | Oracle-only | DuckVEP-only | Terms differ |
| --- | ---: | ---: | ---: | ---: |
| **All transcript pairs** | **400,145** | **0** | **0** | **0** |
| Nuclear | 397,608 | 0 | 0 | 0 |
| Mitochondrial | 2,537 | 0 | 0 | 0 |
| Codon table 1 | 254,309 | 0 | 0 | 0 |
| Codon table 2 | 950 | 0 | 0 | 0 |
| No codon table | 144,886 | 0 | 0 | 0 |
| Forward strand | 202,615 | 0 | 0 | 0 |
| Reverse strand | 197,530 | 0 | 0 | 0 |
| SNV | 100,208 | 0 | 0 | 0 |
| MNV | 99,980 | 0 | 0 | 0 |
| Insertion | 99,978 | 0 | 0 | 0 |
| Deletion | 99,979 | 0 | 0 | 0 |

All 43 modeled biotypes represented in the corpus have exact pairs; their counts and the sequence-ontology term breakdown are in [mouse_seed11639_differential/](mouse_seed11639_differential/). All 400,145 DuckVEP pairs have `supported` status. Twelve table-1 witnesses are `stop_gained`, and all four mitochondrial table-2 TGG→TGA witnesses are `synonymous_variant` in both engines. The [disagreement table](mouse_seed11639_differential/disagreements.csv) has no data rows.

**Admission: yes.** The mouse current-definition receipt in [duckvep_model_receipts.csv](duckvep_model_receipts.csv) is admitted on this exact seeded differential. This tests these 22,204 model-derived variants, not all possible GRCm39 alleles or every noncoding feature.
