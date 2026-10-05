# Whole-HG002 v1/v2 parity

Measured 2026-10-05 with the biological implementation at `bbe2ec4`, the integrated
`check_v2_hg002.sh` runner, DuckDB v1 source `7dbb2e646f`, and v2 source
`fece4143738e2b1d05a851d5c5dc036838aff8ec`. `inputs.sha256` identifies the binaries,
model, VCF, runner and saved reference used. Paths in the receipts are measurement
locators; supply local paths to the runner.

- 158,137 rows on each host, including 157,343 predicted paths.
- Exact schema and bidirectional `EXCEPT ALL` equality, including nested fields.
- Exact equality between the live v1 output and `hg002_bbe2ec4.parquet`.
- Model fingerprint `8650913055114199712` on both hosts.
- Native high-water charges below the 4 GiB ceiling, with zero refusals.
- Peak RSS: v1 4,025,860 KiB (3.84 GiB), v2 4,666,540 KiB (4.45 GiB).

Each host ran once on CPU 19 with one thread and an 8 GB DuckDB memory limit.
Wall times were 15.03 s and 15.28 s. The host had other active workloads;
these timings are descriptive, not a throughput qualification. v2's capture
collections are managed by DuckDB and do not count against the native budget.

Reproduce with `scripts/check_v2_hg002.sh --help`. The full output Parquets are
retained outside Git under `work/luna-six/artifacts/integration/v2-hg002/`.
