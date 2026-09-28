# Rduckvep 0.1.0 build diagnosis

[r-universe run 36446091352](https://github.com/r-universe/rgenomicsetl/actions/runs/36446091352) contains these failed jobs:

| Job | Failure in log |
| --- | --- |
| Windows devel | `test_duckvep.R: Error: DuckVEP artifact for installed DuckDB is not available` |
| Windows release | `test_duckvep.R: Error: DuckVEP artifact for installed DuckDB is not available` |
| Windows oldrel | `test_duckvep.R: Error: DuckVEP artifact for installed DuckDB is not available` |
| macOS arm64 release | HTSlib external-project configure exited 77; CMake reports `Command failed: 77` and `ERROR: configuration failed for package ‘Rduckvep’`. Its private configure log is not included in the job output. |
| macOS arm64 oldrel | Same HTSlib configure exit 77 and package configuration failure; the private configure log is not included. |
| wasm release | `Unsupported DuckDB extension platform` in the package configure script, then `Building wasm binary for package 'Rduckvep' failed.` |

The three Windows binary installations have no `configure.win`, so no extension artifact is built. The macOS command passes all detected libraries to HTSlib's configure, including curl and crypto; exit 77 alone does not establish which configure probe failed. The wasm package build needs a separate toolchain and runtime design; it cannot use the native extension artifact.
