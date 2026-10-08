# DuckDB hosts and release profile

DuckVEP has two independently built DuckDB C API hosts. **Stable C API v1 is the default native build and the community-submission target.** It works on released DuckDB 1.5.x. The v2 host is a separate preview workflow with a pinned SDK; its ABI is not frozen across DuckDB snapshots.

## Stable-v1 submission profile

The merged implementation passes the complete 35-file SQL suite on released DuckDB 1.5.6 (`069cc9f9b5`), alongside the existing DuckDB 1.5.1 CLI and R `duckdb` 1.5.5 qualification. Native distribution CI covers Linux and macOS x86_64/ARM64, Windows MinGW x86_64, and wasm. R package Windows ARM checks are separate from native MSVC extension qualification.

The profile includes documented independent-event annotation and the bounded haplotype domains: observed-phase coding replay, nonoverlapping compound literal HGVS, explicit diploid phase alternatives, equal-length exonic noncoding/UTR replay, and conditional typed-curation peptide pairs. The [function reference](functions.md) and [haplotype contract](../design/duckvep_haplotype_contract.md) define the admission rules, fields and exclusions. Hypotheses are not observed cis; conditional recoding is not an unconditional biological prediction.

Community submission requires the native, stable-v1 SQL, memory/lifecycle, package/platform and generated-source checks to pass for the pinned source revision, together with executable documentation and accurate scope. Broader work in [#11](https://github.com/RGenomicsETL/DuckVEP/issues/11) and [#50](https://github.com/RGenomicsETL/DuckVEP/issues/50) stays open. A defect within the declared profile remains a blocker; completing every future domain or qualifying an unreleased v2 runtime does not gate this stable-v1 submission.

The community entry is being prepared, not yet published. Use the [source-build instructions](../README.md#install) for local evaluation. Signed publication is coordinated with CRAN: a community PR is not a signed release or a CRAN submission.

## Host and distribution boundaries

| | v1 host | v2 host |
| --- | --- | --- |
| Extension API | `v1.2.0` | `v2.0.0`, footer ABI `C_STRUCT` |
| Build | `make release` | `make release_v2` |
| Artifact | `build/release/duckvep.duckdb_extension` | `build/release_v2/duckvep.duckdb_extension` |
| Headers | `duckdb_capi/` | pinned `duckdb_capi_v2/` |

The R package builds and loads the v1 host. Its `v1.2.0` extension API marker is a C API compatibility version, not a claim that every SQL workflow was tested on DuckDB 1.2. [#48](https://github.com/RGenomicsETL/DuckVEP/issues/48) governs released-v2 qualification and the default-host switch; it requires released DuckDB 2.0.0 and R `duckdb` 2.0 on CRAN. Preview parity does not satisfy that gate, and that future switch is separate from the stable-v1 community submission.

`host_v2/` is a separate CMake project. A v1 build does not compile or link its sources. Both hosts compile the host-neutral implementation in `src/kernel/` and `src/core/`; the v1 adapter is in `src/`, and the v2 adapter is in `host_v2/`. The adapters translate DuckDB vectors and callbacks to the shared row logic and write host-specific results.

## v2 SDK pin and ABI

`duckvep-package.json` pins revision `fece4143738e2b1d05a851d5c5dc036838aff8ec` and the SHA-256 values of `duckdb_v2.h` and `duckdb_extension_v2.h`. `duckdb_capi_v2/REVISION` records the same commit. `make check-v2-sdk` verifies these files. The host compiles with `DUCKDB_V2_API_ALLOW_UNSTABLE=0` and `DUCKDB_V2_API_ALLOW_DEPRECATED=0`.

The v2 ABI is tied to the SDK snapshot, not just the `v2.0.0` API label. The test runner checks that the DuckDB CLI's `source_id` matches the pinned revision. Measured PyPI preview builds had snapshot-specific outcomes: source `ca15f79c32` (`2.0.0.dev2609250715`) segfaulted while loading the extension; `dev2609221243` and `dev2609222040` failed initialization with out-of-memory errors. The v2 test host is built at the pinned revision (`scripts/stage-v2-duckdb.py`); those observations do not characterize other builds.

## Build and gates

```sh
make release_v2
make test-extension-symbols-v2
make check-v2-footer
python3 scripts/stage-v2-duckdb.py
DUCKVEP_V2_DUCKDB=.deps-v2/duckdb-build/duckdb make test_v2
```

`make test_v2` checks the SDK pin, builds the separate host, verifies that the binary exports only `duckvep_init_c_api_v2` and imports no `duckdb_*` symbols, checks the `C_STRUCT` / `v2.0.0` footer, and runs `test/scripts/check_v2_host.py static`. The static gate confirms that v2 calls are stable functions from the pinned header and that `src/core/` and `src/kernel/` contain no DuckDB API calls.

The runtime suite tests repeated `LOAD` on writable and read-only databases, with no catalog DDL or database-file changes; v2-only SQL assertions; model and budget behavior; and the 601 cases in `test/sql_v2/equality_golden.json`. The golden records v1 results. When `build/release/duckvep.duckdb_extension` exists, the runner also checks those results against the live v1 extension. `test/sql_v2/README.md` documents the SQL and the separate whole-HG002 check.

## Release preflight

`make check-v2-release-selftest` runs offline negative controls. A release admission uses `make check-v2-release-preflight` with concrete arguments for the extension, pinned CLI, R library, and completed HG002 output directory. The preflight checks the SDK manifest and headers, the actual extension footer and load, runtime and source revisions, and the complete hashed HG002 comparison receipt. It contacts GitHub and CRAN only with explicit `--live-status`. The current preview manifest is blocked; the command can emit `READY_TO_QUALIFY` after upstream releases are visible but before local artifacts and qualification outputs are inspected. Details and the required invocation are in [`test/sql_v2/release-readiness.md`](../test/sql_v2/release-readiness.md).

`LOAD` registers the native functions without executing SQL or creating persistent catalog objects. The read-only database test covers this contract.

Geometry, repeat/phase, SQL-builder, annotation, coding-discovery and resource-control functions retain their v1 SQL names. The model and haplotype query workflows use v2-specific load builders; the workflow differences are described below. `duckvep_model_drop`, `duckvep_model_save` and `duckvep_model_restore` retain their v1 names.

## Ownership, borrowed data and concurrency

The v2 callback's error-info handle is borrowed. Each fallible API call returns its own detail handle; the adapter copies the error into the callback handle and destroys the detail. Input vector views and VARCHAR payloads are borrowed from the execution chunk and are valid only during that call. Core haplotype row views remain valid until the host fetches the next input row.

Vector resizing can invalidate data and validity pointers. The adapter reacquires them after `vector_set_size`; nested child pointers are taken after sizing their vectors. Dictionary reads can flatten an aliased child, so `load_views` takes the returned views twice and uses the second set. Returned owned logical-type, value and error-detail handles are destroyed by their owner.

Registration user data carries destructors and reference counts. The model registry is mutex-protected; scans pin model entries for the duration of use, so a pinned entry remains alive while a scan reads it. Staged model relations and haplotype jobs are shared state protected by the v2 host's mutex. Cursor-bearing scans (`duckvep_so_terms`, haplotype scan and record-plan scan, and `duckvep_coding_calls`) set their global table-function state to one thread. Annotation workers use the shared process-wide native budget and worker limits.

## Model loading

The v1 `duckvep_model_load(name, regions_query, transcripts_query, exons_query, ...)` executes query strings on its private connection in a transaction. It sees committed permanent relations, not the caller's TEMP tables or uncommitted rows.

The v2 host cannot execute SQL on a private connection inside a callback. `duckvep_model_load_sql` returns ordered caller-side statements: a `COPY ... TO 'duckvep_stage'` for each supplied relation, followed by `duckvep_model_publish`. The caller executes the statements in order and in its own session. The three required relations are regions, transcripts and exons; mature miRNA, peptide edits and interval features are optional. Queries must return the expected typed columns in the required order, and the COPY preserves row order.

```sql
SELECT unnest(duckvep_model_load_sql(
  'my_model', 'SELECT ...', 'SELECT ...', 'SELECT ...'));
-- Execute each returned statement in order, including the COPY statements and publish call.
```

For v2, TEMP tables and uncommitted rows visible to the caller's COPY can supply the model. Staged rows live in DuckDB-managed column collections, which can spill under DuckDB's memory limit; model arrays and indexes are charged to DuckVEP's native-memory budget. A failed or cancelled COPY stages no rows. A completed COPY replaces the staged relation of the same name. `duckvep_model_publish` consumes the staged relations on success or failure, validates and installs the model through the same core loaders as v1. Model installation is extension state, not SQL transaction state: a surrounding transaction rollback does not remove a published model on either host. Use `duckvep_model_drop` to release it.

The visible workflow differs by host: v1 performs loading in one table-function call; v2 returns SQL for the caller to run and then publishes the staged rows. Validation uses the shared core loaders. `duckvep_model_save` and `duckvep_model_restore` use the shared snapshot implementation on both hosts.

## Haplotype input

The v1 `duckvep_haplotypes(calls_query, model, ...)` executes and normalizes its query on a private connection. On v2, `duckvep_haplotype_load_sql` returns caller-side normalization and staging statements; run them in order, then read `duckvep_haplotype_scan(job)`. The scan replays the captured rows through `src/core/duckvep_core_haplotypes.c`, the same core used by v1. With `input_mode = 'source_records'`, set `phase_policy = 'vep_compat'`; the builder returns four ordered statements for raw input, record-plan staging, normalized-call staging and cleanup. A job scan consumes the job once and releases it at completion or error; use `duckvep_haplotype_drop` to discard a staged job that will not be scanned.

The staging collections are DuckDB-managed and governed by DuckDB's memory limit, not the native-memory budget. The scan's workspace and core allocations use the shared native budget. `duckvep_coding_calls` uses the same budget for its HTSlib reader and discovery workspace.

## Haplotype arrangements

`duckvep_haplotype_arrangements(calls_query, model, ...)` on v1 normalizes and replays strict diploid calls on its private query connection. It therefore sees committed permanent relations, but not caller TEMP relations or uncommitted rows. On v2, call `duckvep_haplotype_arrangements_load_sql(calls_query, model, job [, options])`, execute its single returned caller-side `COPY`, then query `duckvep_haplotype_arrangements(job)`. The staged query uses the same shared strict normalization recipe as v1 and preserves event IDs, REF/ALT values, both original alleles, `phase_before`, and `phase_set`.

The v2 arrangement job is consumed once. The scan fully validates and charges the native arrangement workspace before returning an output row, then releases its cursor and staged collection. A failed or cancelled `COPY` publishes no job. Arrangement limits are accepted in the loader's optional options struct: `max_sites`, `max_calls`, `max_arrangements`, and `max_replays`; zero, invalid, and over-capacity values are rejected before a job is published. `CALLED` is not evidence of cis phase: unphased calls produce bounded hypotheses, while only observed phase metadata constrains the lanes.

## Whole-HG002 parity measurement

The receipt at [`benchmarks/data/v2_hg002/bbe2ec4/`](../benchmarks/data/v2_hg002/bbe2ec4/README.md) records a run measured on 2026-10-05 with biological implementation `bbe2ec4`, v1 DuckDB source `7dbb2e646f`, and v2 DuckDB source `fece4143738e2b1d05a851d5c5dc036838aff8ec`. The run used one thread on CPU 19 and an 8 GB DuckDB memory limit.

Both hosts produced 158,137 rows, including 157,343 predicted paths. Schema and the complete nested row multiset matched by bidirectional `EXCEPT ALL`; the v1 output also matched the saved reference Parquet. The model fingerprint was `8650913055114199712` on each host. Native-memory high-water stayed below 4 GiB, with zero capacity refusals. Peak RSS was 4,025,860 KiB (3.84 GiB) on v1 and 4,666,540 KiB (4.45 GiB) on v2. The v2 capture collections are DuckDB-managed and outside the native budget. Wall times were 15.03 s on v1 and 15.28 s on v2; the recorded host had other active workloads, so these values are descriptive rather than a throughput qualification.

The receipt includes hashes, schemas, budget snapshots, inputs and outputs' provenance. `scripts/check_v2_hg002.sh --help` describes the reproduction inputs and options. The full-data run is separate from `make test_v2`.
