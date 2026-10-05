# DuckDB C API v2 SDK (pinned preview)

`duckdb_v2.h` and `duckdb_extension_v2.h` match the revision and SHA-256 values recorded under `v2_host` in `../duckvep-package.json`; `REVISION` repeats that commit. Only the v2 host in `../host_v2/` uses these headers. The distributed v1 host uses `../duckdb_capi/`.

From the repository root, `scripts/fetch-v2-sdk.py` verifies the checksums and `make check-v2-sdk` runs verification without refreshing the files.

The v2 host compiles with `DUCKDB_V2_API_ALLOW_UNSTABLE=0` and `DUCKDB_V2_API_ALLOW_DEPRECATED=0`; its footer uses ABI `C_STRUCT` and extension API `v2.0.0`. The preview ABI is tied to the DuckDB source revision. Validate the extension against a DuckDB build at the pinned revision; the `v2.0.0` label alone does not establish cross-snapshot compatibility.
