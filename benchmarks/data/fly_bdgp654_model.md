# Drosophila melanogaster BDGP6.54 DuckVEP model

The Ensembl 116 model uses the pinned core snapshot and FASTA in [fly_bdgp654_sources.md](fly_bdgp654_sources.md):

```sh
Rscript scripts/build_species_model.R \
  /tmp/duckvep-fly/core.duckdb /tmp/duckvep-fly/reference_chunks.parquet \
  /tmp/duckvep-fly/fly.fa /tmp/duckvep-fly/fly.duckdb BDGP6.54 116 \
  "$(pwd)/build/release/extension/duckvep/duckvep.duckdb_extension"
```

The builder exposes eleven read-only core tables as local views to its SQL generators. The [current-definition receipt](duckvep_model_receipts.csv) pins model SHA-256 `e6deea1ac2b0097589df9da501bd4a1ddf47ec5f8e0d4291ee49e727352533c8`: 1,870 regions, 143,726,002 reference bases, 41,600 transcripts, 24,254 genes, 30,802 coding transcripts, 30,710 sequence-backed transcripts and 92 sequence-withheld transcripts. There are 30,697 table-1 coding transcripts and 13 table-5 mitochondrial coding transcripts. Local DuckDB artifact SHA-256: `a2049b827b6a59199be39e4248b882fc2b11ef618d83cf329bca79bd4f692172`.

The receipt's `source_revision` is the full hash `0ceeeba4adeece46a44f333f71de77922ab10447` on `origin/main`. The release extension was built with `make release`; `git diff origin/main --exit-code -- src duckdb_capi third_party cmake CMakeLists.txt` confirmed that its extension sources match that revision. The model SHA-256 was reproduced by an independent build from the same staged inputs.
