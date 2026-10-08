# DuckDB v2 release preflight

`test/scripts/check_v2_release_preflight.R` checks released-host qualification inputs against concrete artifacts. `ADMITTED` applies to this preflight; biological conformance, platform checks and signed distribution remain separate release gates. The checker fails closed and does not install software, fetch headers, or access the network unless `--live-status` is explicit.

An admission requires all of the following:

- the checked-out `duckvep-package.json`, `duckdb_capi_v2/REVISION`, and pinned header hashes pass `scripts/fetch-v2-sdk.py`;
- the repository manifest marks the v2 host as `release`;
- `--live-status` finds DuckDB `v2.0.0` through `gh api` and R `duckdb` `2.0.0` in explicit CRAN metadata;
- `--extension` passes the existing v2 footer checker as `v2.0.0` with `C_STRUCT`, then loads in `--v2-cli`;
- the supplied CLI and the R `duckdb` package loaded from `--r-library` report `2.0.0` and a `source_id` matching the pinned SDK revision; the extension also loads in that R runtime;
- `--hg002-dir` is a completed `scripts/check_v2_hg002.sh` output directory whose receipt, concrete model/output paths and hashes for both CLIs, both extensions, VCF and model validate. Its exact comparison must report equal row counts and zero differences; fingerprints are compared as exact strings, and nonnegative integer budget measurements remain within their limits. Its v2 extension path and source revision must match the supplied artifact and CLI.

An offline `--offline-status` JSON file is diagnostic input only and can never admit a release. It is useful for a disconnected review, but live upstream status and local artifact qualification remain required. The current preview manifest therefore reports `BLOCKED`; it is not release evidence.

Run the offline negative controls:

```sh
make check-v2-release-selftest
```

They reject an unstable footer ABI, a mismatched SDK source, missing qualification, unequal or invalid row counts, differing 64-bit fingerprints, negative budgets, and incomplete input-hash coverage. Supply all real paths to run a preflight:

```sh
make check-v2-release-preflight V2_RELEASE_PREFLIGHT_ARGS='\
  --live-status \
  --extension /absolute/path/duckvep.duckdb_extension \
  --v2-cli /absolute/path/duckdb \
  --r-library /absolute/path/R/library \
  --hg002-dir /absolute/path/hg002-run'
```

`ADMITTED` is emitted only after every check succeeds. `READY_TO_QUALIFY` means that live upstream versions are visible but artifact or HG002 qualification is missing. `BLOCKED` identifies missing, invalid, preview, or mismatched evidence.
