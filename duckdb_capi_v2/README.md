# DuckDB C API v2 SDK (pinned preview)

`duckdb_v2.h` and `duckdb_extension_v2.h` are vendored at the DuckDB revision
and SHA-256 checksums recorded under `v2_host` in `../duckvep-package.json`
(`REVISION` repeats the revision). They are used only by the v2 host
(`host_v2/`); the v1 host keeps using `../duckdb_capi/`.

`scripts/fetch-v2-sdk.py` verifies the checksums, and refreshes the files from
`raw.githubusercontent.com/duckdb/duckdb/<revision>/src/include/` when they
differ. `make check-v2-sdk` runs it in verify-only mode.

The v2 host compiles with `DUCKDB_V2_API_ALLOW_UNSTABLE=0` and
`DUCKDB_V2_API_ALLOW_DEPRECATED=0`, and its footer is `C_STRUCT` with
extension API `v2.0.0`. The inspected preview SDK still warns that the
extension ABI is not frozen: validate against a DuckDB built at this revision,
and do not infer a cross-release compatibility guarantee from the version string.
