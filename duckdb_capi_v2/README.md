# DuckDB v2.0.0 C API SDK

`duckdb_v2.h` and `duckdb_extension_v2.h` come from the DuckDB `v2.0.0` tag at commit `d54cbbe3df17e8e5ce68312f483b4155836112aa`. The `v2_host` section of `../duckvep-package.json` records their SHA-256 values. `REVISION` contains the commit.

Only `../host_v2/` uses these headers. The default stable-v1 host uses `../duckdb_capi/`. The v2 host requires `DUCKDB_V2_API_ALLOW_UNSTABLE=1` and keeps `DUCKDB_V2_API_ALLOW_DEPRECATED=0`. DuckDB v2.0.0 has no stable v2 extension API surface.

## Verify the SDK

From the repository root, run:

```sh
python3 scripts/fetch-v2-sdk.py
```

The command verifies the headers and revision without a download. Its `--refresh` option retrieves the same files from the pinned commit.

## Build and package

Build the separate CMake project in `host_v2/`. Append a footer with the repository tool `extension-ci-tools/scripts/append_extension_metadata.py`. Use these values:

- `--abi-type C_STRUCT_UNSTABLE`
- `-dv v2.0.0` (the exact DuckDB engine version)
- The target platform and extension version

The extension entrypoint requests the v2 API family. The unstable footer restricts LOAD to DuckDB v2.0.0. A `C_STRUCT` footer does not satisfy this contract. Do not use preview snapshots or another DuckDB version with this artifact.

The root Makefile, footer checks, static API checks, and CI must use this contract before a release. Those integration files are separate from this SDK. The native port does not establish full-suite equality, performance, or release readiness.
