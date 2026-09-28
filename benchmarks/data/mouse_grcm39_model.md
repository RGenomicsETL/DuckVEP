# Mus musculus GRCm39 DuckVEP model

The Ensembl 116 mouse model is built from the snapshot and reference pinned in [mouse_grcm39_sources.md](mouse_grcm39_sources.md):

```sh
Rscript scripts/build_species_model.R \
  /tmp/duckvep-mouse/core.duckdb /tmp/duckvep-mouse/reference_chunks.parquet \
  /tmp/duckvep-mouse/mouse.fa /tmp/duckvep-mouse/mouse.duckdb GRCm39 116 \
  "$(pwd)/build/release/duckvep.duckdb_extension" scripts/mouse_core_input.sql
```

The build connection loads DuckVEP alone. The snapshot retains every source row. One miRNA attribute on transcript `ENSMUST00000175144` (`transcript_id=144081`, `attrib_type_id=15`, value `35-56`) extends three bases beyond its 53-nt spliced cDNA; the core-input view omits that attribute while retaining the transcript and its valid `1-23` mature segment. This prevents an invalid cDNA range from entering the model without inventing a corrected mature sequence.

The current-definition receipt is the `mus_musculus` row in [duckvep_model_receipts.csv](duckvep_model_receipts.csv): model SHA-256 `ae39ffc9e647d0a096a13b737599d480fae2b94d3a202938a8eeaa44208c62be`, 61 regions, 2,728,222,451 reference bases, 481,483 transcripts, 78,077 genes, and 269,905 coding transcripts. Codon tables are 1 (269,892 coding transcripts) and 2 (13 mitochondrial coding transcripts). The local DuckDB artifact SHA-256 is `da6e3eb7e86dc39f0ea9499f8b4b2dd239c3310e9d00ff8cfc60721c2542807e`.
