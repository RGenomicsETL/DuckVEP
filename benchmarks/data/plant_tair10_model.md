# Arabidopsis thaliana TAIR10 DuckVEP model

The Ensembl Plants 63 core, FASTA and indexed cache are pinned in [plant_tair10_sources.md](plant_tair10_sources.md). The core snapshot is attached read-only to the generalized builder:

```sh
Rscript scripts/build_species_model.R \
  /tmp/duckvep-plant/core.duckdb /tmp/duckvep-plant/reference_chunks.parquet \
  /tmp/duckvep-plant/plant.fa /tmp/duckvep-plant/plant.duckdb TAIR10 63/116 \
  "$(pwd)/build/release/extension/duckvep/duckvep.duckdb_extension"
```

The [current-definition receipt](duckvep_model_receipts.csv) pins model SHA-256 `5b3fac7c92db9c2dfd13cc83df457a993ece29b5b703343194c5363136d0a12f`: 7 regions, 119,667,750 reference bases, 54,013 transcripts, 32,833 genes, 48,321 coding transcripts, 48,316 sequence-backed transcripts, 5 sequence-withheld transcripts, 313,952 exon memberships and 325 mature-miRNA segments. The sequence-backed coding transcripts use table 1 (48,228; chromosomes `1`–`5` and `Mt`) or table 11 (88; `Pt`). Local DuckDB artifact SHA-256: `35b751d1e80777b310b0da4b78b4b9718683e06382804d30d2fcf2d3a1b56cde`.

`source_revision` is `013ce4fdc0134b7247614a8515f1c908165d85ea`, already on `origin/main`. The extension was built with `make release`; `git diff origin/main --exit-code -- src duckdb_capi third_party cmake CMakeLists.txt` confirmed identical extension sources. An independent build from the same staged inputs reproduced the model SHA-256.
