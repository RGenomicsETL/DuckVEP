# DuckDB C API v1 headers

This directory supplies the headers used by the distributed v1 host. The root `Makefile` pins the extension API footer to `TARGET_DUCKDB_VERSION=v1.2.0` and the header source version to `DUCKDB_HEADER_VERSION=v1.5.3`; these pins serve different purposes.

Refresh the headers with:

```sh
make update_duckdb_headers
```

Keep `duckdb.h` and `duckdb_extension.h` from the same header release. The extension API target remains the version in `TARGET_DUCKDB_VERSION`.
