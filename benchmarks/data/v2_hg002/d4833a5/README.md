# Whole-HG002 v1/v2 parity at the 2.0.0 release candidate

Measured 2026-11-19 with the implementation at `d4833a5`, the `check_v2_hg002.sh`
runner, DuckDB v1 source `7dbb2e646f`, and the released DuckDB v2.0.0 source
`d54cbbe3df17e8e5ce68312f483b4155836112aa`. `inputs.sha256` identifies the binaries,
model, VCF and runner used. Paths in the receipts are measurement locators; supply
local paths to the runner.

- 158,137 rows on each host.
- Exact schema and bidirectional `EXCEPT ALL` equality, including nested fields.
- Model fingerprint `8650913055114199712` on both hosts.
- The v1 Parquet output is byte-identical to the output of the run at
  [`1967499`](../1967499/README.md), which therefore describes the same rows.
- Native high-water charges at most 1,758,619,872 bytes under the 4 GiB ceiling,
  with zero refusals.
- Peak RSS: v1 4,038,968 KiB (3.85 GiB), v2 4,689,272 KiB (4.47 GiB).

Each host ran once on CPU 4 with one thread and an 8 GB DuckDB memory limit.
Wall times were 10.07 s and 9.72 s. These timings are descriptive, not a throughput
qualification. v2's capture collections are managed by DuckDB and do not count
against the native budget.

Reproduce with `scripts/check_v2_hg002.sh --help`.
