# DuckDB v2 host tests

Run the build, binary, static and SQL gates with a DuckDB CLI built at the SDK revision pinned in `duckvep-package.json`:

```sh
DUCKVEP_V2_DUCKDB=/path/to/pinned/duckdb make test_v2
```

The runtime suite checks `LOAD` twice on writable and read-only databases, no catalog or database-file changes, v2-only assertions, model-sink and native-budget behavior, and the 585 equality cases in `equality_golden.json`. The golden contains v1-host results; if a v1 extension exists at `build/release/duckvep.duckdb_extension`, the runner compares it with the live v1 host too. A case named `... error` must produce a DuckVEP error; other cases must succeed. The suite and API pin are detailed in [the v2 host guide](../../docs/v2-host.md).

## Whole-HG002 parity

`scripts/check_v2_hg002.sh` is a separate full-data check. It runs v1 and v2 with separate DuckDB CLIs, records stage timings, peak RSS, native-budget snapshots, SQL and Parquet outputs, and compares fingerprints, schema and the complete nested row multiset with two-way `EXCEPT ALL`. `--v1-reference` additionally compares the fresh v1 output with the supplied saved Parquet file. The measured parity receipt is [`benchmarks/data/v2_hg002/bbe2ec4/`](../../benchmarks/data/v2_hg002/bbe2ec4/README.md); its wall times are descriptive measurements, not a throughput qualification.

Build matching extensions before the check. The v2 CLI must match the pinned SDK revision. Supply local data paths through flags or environment variables; the output directory must be empty:

```sh
scripts/check_v2_hg002.sh \
  --v1-cli "$(command -v duckdb)" \
  --v2-cli "$DUCKVEP_V2_DUCKDB" \
  --v1-extension build/release/duckvep.duckdb_extension \
  --v2-extension build/release_v2/duckvep.duckdb_extension \
  --vcf "$HG002_VCF" --model "$ENSEMBL_MODEL" \
  --v1-reference "$HG002_V1_PARQUET" --out-dir "$HG002_RUN_DIR" \
  --cpu-list 2 --threads 1
```
