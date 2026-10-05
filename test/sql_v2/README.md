# v2 host tests

- `equality_cases.sql`: queries run on both hosts; results are compared with
  `equality_golden.json`, the recorded output of the v1 host (`run_v2_tests.py --record`).
  A case named `... error` must fail with a `duckvep_*:` message; every other case must succeed.
- `v2_native.sql`: v2-only assertions (registration kinds, no catalog objects, NULL
  invariants of struct results, table-function rescans).
- Run with `DUCKVEP_V2_DUCKDB=<DuckDB CLI built at the pinned revision> make test_v2`.
  The PyPI `duckdb>=2.0.0.dev0` wheels are different snapshots of a preview ABI that is not
  frozen and cannot load this extension; see `docs/v2-host.md`.

## HG002 host-scale check

`scripts/check_v2_hg002.sh` runs the v1 haplotype query and the v2 builder/COPY/scan
path with separate DuckDB CLIs, pins each process with `taskset`, and records CLI stage
timers, peak RSS, native-budget snapshots, SQL scripts and Parquet outputs. It checks
model fingerprints, schema and the complete nested row multiset with two-way
`EXCEPT ALL`; it does not use hashes as an equality proof. `--v1-reference` additionally
compares the freshly generated v1 output with a saved v1 Parquet file.

Build compatible extensions first and use a v2 CLI from the revision pinned in
`duckvep-package.json`. Example inputs are supplied through environment variables or
flags, not repository-specific data paths:

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

The output directory must be empty. The full-data check is separate from `make test_v2`
and must be run serially by the coordinator; no full HG002 timing run is implied by
this SQL test suite.
