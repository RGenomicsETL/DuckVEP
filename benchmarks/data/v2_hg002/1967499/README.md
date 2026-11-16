# Whole-HG002 v1/v2 parity

Measured 2026-11-15 with the biological implementation at `1967499`, the
`check_v2_hg002.sh` runner, DuckDB v1 source `7dbb2e646f`, and the released DuckDB
v2.0.0 source `d54cbbe3df17e8e5ce68312f483b4155836112aa`. `inputs.sha256` identifies
the binaries, model, VCF and runner used. Paths in the receipts are measurement
locators; supply local paths to the runner.

- 158,137 rows on each host, including 157,635 predicted paths.
- Exact schema and bidirectional `EXCEPT ALL` equality, including nested fields.
- Model fingerprint `8650913055114199712` on both hosts.
- Native high-water charges below the 4 GiB ceiling, with zero refusals.
- Peak RSS: v1 4,039,124 KiB (3.85 GiB), v2 4,692,000 KiB (4.47 GiB).

Each host ran once on CPU 4 with one thread and an 8 GB DuckDB memory limit.
Wall times were 9.76 s and 9.57 s. These timings are descriptive, not a throughput
qualification. v2's capture collections are managed by DuckDB and do not count
against the native budget.

The output has 36 columns. The saved reference `hg002_bbe2ec4.parquet` has 32 and
predicts 157,343 paths, so this run makes no comparison with it.

Reproduce with `scripts/check_v2_hg002.sh --help`.
