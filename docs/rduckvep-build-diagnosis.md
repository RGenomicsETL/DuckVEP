# Rduckvep 0.1.0 build diagnosis

The [r-universe run 36446091352](https://github.com/r-universe/rgenomicsetl/actions/runs/36446091352), recorded 2026-09-29, reported these failures:

| Job | Recorded failure |
| --- | --- |
| Windows devel, release and oldrel | `test_duckvep.R: Error: DuckVEP artifact for installed DuckDB is not available` |
| macOS arm64 release and oldrel | HTSlib external-project configure exited 77; CMake reported `Command failed: 77` and `ERROR: configuration failed for package ‘Rduckvep’`. The job output omitted the private configure log. |
| wasm release | `Unsupported DuckDB extension platform`, followed by `Building wasm binary for package 'Rduckvep' failed.` |

The Windows binary-installation jobs in that run had no `configure.win` and produced no extension artifact. This checkout contains `r/Rduckvep/configure.win`; the logged failures do not establish the outcome of this checkout's Windows configuration.

The macOS configure command in those jobs supplied all detected libraries, including curl and crypto. Exit 77 does not identify the failing configure probe; the missing private logs leave that cause unresolved. The wasm error is specific to the logged package configuration. A native extension artifact is not a wasm artifact.
