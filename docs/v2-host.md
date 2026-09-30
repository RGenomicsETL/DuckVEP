# The DuckDB C API v2 host (preview)

DuckVEP is being ported to DuckDB's v2 C API as a second *host* next to the
existing v1 build (issue #8). This page is the plan and the porting map: how the
two hosts are laid out, what is ported, and which slice moves each of the 26
public functions.

## Status

| | host_v1 | host_v2 |
| --- | --- | --- |
| Extension API | stable C API v1.2.0 | stable C API v2.0.0, footer `C_STRUCT` |
| Build | `make release` (the root `CMakeLists.txt`) | `make release_v2` (`host_v2/CMakeLists.txt`) |
| Artifact | `build/release/duckvep.duckdb_extension` | `build/release_v2/duckvep.duckdb_extension` |
| Ships | yes: CRAN and the community repository | no: preview, waits for DuckDB 2.0.0 |
| Functions | all 26 public functions | slices 1 to 4a: `duckvep_so_terms`, `duckvep_allele_geometry`, `duckvep_breakend_geometry`, `duckvep_repeat_alleles`, `duckvep_phase_call`, the twelve `duckvep_*_sql` builders, `duckvep_model_load` / `duckvep_model_drop` (as the COPY sink below), and the internal `_duckvep_revcomp`, `_duckvep_raw_gt`, `_duckvep_record_order` |

The v2 host is a separate CMake project, so a v1 build does not compile or link
any v2 file, and `duckdb_capi/` and MainDistributionPipeline are untouched. Since
slice 2 the v1 host calls the same `src/core` row logic as the v2 host for the
functions ported so far (its callbacks keep their own DuckDB vector reads and
writes). The v2 host uses **stable v2 functions only**, compiled with
`DUCKDB_V2_API_ALLOW_UNSTABLE=0` and `DUCKDB_V2_API_ALLOW_DEPRECATED=0`, against a
**pinned** preview SDK. It is C only.

## Layout

```text
src/kernel/          host-neutral kernel (no DuckDB symbol)
src/core/            host-neutral row logic shared by both hosts (no DuckDB symbol)
src/*.c              DuckDB-facing v1 adapter (+ src/duckvep_v1_cells.h, its cell reader)
host_v2/             DuckDB-facing v2 adapter: host_v2.c, host_v2_nested.c,
                     host_v2_common.h (helpers), and its own CMakeLists.txt
duckdb_capi/         v1 headers (R bundle and CI depend on this path)
duckdb_capi_v2/      pinned v2 headers, REVISION, LICENSE, README
duckvep-package.json the pins for both hosts
```

`src/core/` holds what used to be inside v1 callbacks, so that the two hosts
cannot drift:

| File | Logic |
| --- | --- |
| `duckvep_core_geometry.{c,h}` | `duckvep_so_terms` rows, `duckvep_allele_geometry`, the breakend messages |
| `duckvep_core_cells.{c,h}` | typed cells (`duckvep_cell_t`: a host reads one vector element into one), numeric conversion of any integer, float, HUGEINT or scaled DECIMAL cell, the phase allele and flag conversions, and the option-STRUCT field checks |
| `duckvep_core_repeat.{c,h}` | `duckvep_repeat_alleles`: option cap, row validation and status, required lengths, rendering, direction |
| `duckvep_core_sql.{c,h}` | the SQL text buffer, identifier and literal quoting, relation-name splitting, the NUL-checked string copy |
| `duckvep_core_annotate.{c,h}` | `duckvep_annotate_sql`, `duckvep_annotate_projected_sql`, `duckvep_transcript_projection_sql`: templates and option rendering |
| `duckvep_core_ensembl.{c,h}` | the three Ensembl builders and `duckvep_model_receipt_sql`: templates, species option, receipt body |
| `duckvep_core_prepare.{c,h}` | `duckvep_prepare_sv_geometry_sql`, `duckvep_prepare_expansionhunter_sql` and the BND identity, breakend gene and structural HGVS builders, with the `max_span` option |
| `duckvep_core_model.{c,h}` | the resident model: arrays, the six relation loaders and their validation, FASTA check, publication, lifting of circular regions, kernel open, the model registry, and the fingerprint. Loaders read *row sources* (`duckvep_source_t`: typed, flat column batches) and never a DuckDB result |
| `duckvep_core_model_script.{c,h}` | the statements of the v2 load (below) |
| `duckvep_core_phase.{c,h}` | `duckvep_phase_call`: list checks, policy and phase-set options, the two-pass slot reducer with a reader interface (`allele`, `phase`, `emit` callbacks), and the `_duckvep_revcomp`, `_duckvep_raw_gt` and `_duckvep_record_order` kernels |

A host's job is to turn vectors into cells (its type switch and its selection
and validity handling), call the core, and write the results: v1 does this in
`src/duckvep_repeat.c`, `src/duckvep_phase_sql.c`, `src/duckvep.c`,
`src/duckvep_annotate.c` and `src/duckvep_sql.c`; v2 in `host_v2/`. `src/core` is
under `src/`, so it is part of the R bundle (`Rscript r/Rduckvep/bootstrap.R .`).
What is still embedded in v1 callbacks (the builders, the model, the haplotype
scanner, the annotation natives) moves to the core as its slice reaches it.
The v1-vs-v2 equality test is the parity guard for the ported functions.

## Pins

The pins live in `duckvep-package.json` (not `description.yml`, which is the
community-repository manifest and has no place for a second host):

- `v2_host.duckdb_sdk_revision`: `fece4143738e2b1d05a851d5c5dc036838aff8ec`,
  the same DuckDB revision ducksassy (the maintainer's reference) pins.
- `v2_host.duckdb_sdk_sha256`: checksums of `duckdb_v2.h` and
  `duckdb_extension_v2.h`; `make check-v2-sdk` verifies them
  (`scripts/fetch-v2-sdk.py --refresh` re-fetches from the revision).
- `v1_host`: points at `duckdb_capi/` and the Makefile's `TARGET_DUCKDB_VERSION`.

The v2 ABI is not frozen. Two preview binaries that both say `v2.0.0` can have
different function tables, and a v2 extension loaded by the wrong snapshot
crashes or fails to initialize. Measured on this slice: the host loads and passes
on a DuckDB built at `fece414373`, segfaults on the PyPI wheel
`2.0.0.dev2609250715` (`ca15f79c32`), and fails with an out-of-memory initialization
error on `dev2609221243` and `dev2609222040`. So the v2 tests run on a DuckDB CLI
built at the pinned revision (`scripts/stage-v2-duckdb.py`), and the runner refuses
any other engine. `pip install --pre duckdb` is therefore not a usable v2 test host
until the extension ABI is frozen; the cross-release matrix waits for that.

## Build and test

```sh
make release_v2                 # build/release_v2/duckvep.duckdb_extension
make test-extension-symbols-v2  # exports exactly duckvep_init_c_api_v2; no duckdb_* import
make check-v2-footer            # ABI C_STRUCT, extension API v2.0.0
python3 scripts/stage-v2-duckdb.py                # DuckDB CLI at the pinned revision
DUCKVEP_V2_DUCKDB=.deps-v2/duckdb-build/duckdb make test_v2
DUCKVEP_V2_DUCKDB=... make test_v2_asan           # the same tests under ASan and UBSan
```

`make test_v2` runs `test/scripts/check_v2_host.py static` (only stable
`duckdb_v2_*` calls, opt-ins off, `host_v2/core` and `src/kernel` free of DuckDB),
then `test/scripts/run_v2_tests.py`: LOAD on a writable and a read-only primary
(twice, no DDL, database bytes unchanged), `test/sql_v2/v2_native.sql`, and the 326
cases of `test/sql_v2/equality_cases.sql` compared with the recorded v1 host
(`equality_golden.json`, from `run_v2_tests.py --record`, also re-checked live
against the v1 build when `build/release` exists). CI is `.github/workflows/v2-host.yml`.

A "no unstable symbol" check for a v2 binary cannot be an `nm` import scan like
`test-extension-symbols`, because a v2 extension imports nothing from DuckDB: it
reaches the engine through the function table passed to its entry point. The gate
is therefore three checks together: the source compiles with the unstable surface
compiled out, every called function has `stable: v2.0.0` and no
`unstable`/`deprecated` history in the pinned header, and the binary imports no
`duckdb_*` symbol and exports only the entry point.

## Porting map

Every public function is in exactly one family and one planned slice. The six
slices are Astra's (memo section 4); slice 1 is this one.

| Slice | Family | What it proves | Functions |
| --- | --- | --- | --- |
| 1 | ABI, vectors and types | entry point, footer, types, NULL/selection/constant/dictionary vectors, errors, table-function state | `duckvep_so_terms`, `duckvep_allele_geometry`, `duckvep_breakend_geometry` (done) |
| 2 | Vectors and types | nested LIST/STRUCT input and output, list-child capacity, >2,048-element nested output, DECIMAL/HUGEINT, ANY parameters and overloads, option STRUCTs (done) | `duckvep_repeat_alleles`, `duckvep_phase_call` (and internal `_duckvep_revcomp`, `_duckvep_raw_gt`, `_duckvep_record_order`) |
| 3 | SQL-builder surface | builders only: VARCHAR and STRUCT option inputs, SQL text out, no macros, no DDL at LOAD (done) | `duckvep_ensembl_regions_sql`, `duckvep_ensembl_transcripts_sql`, `duckvep_ensembl_regulation_features_sql`, `duckvep_model_receipt_sql`, `duckvep_annotate_sql`, `duckvep_annotate_projected_sql`, `duckvep_transcript_projection_sql`, `duckvep_prepare_sv_geometry_sql`, `duckvep_prepare_breakend_pairs_sql`, `duckvep_prepare_breakend_fusion_sql`, `duckvep_prepare_structural_hgvs_sql`, `duckvep_prepare_expansionhunter_sql` |
| 4a | Model sink (COPY staging) | stable-v2 COPY callbacks take the caller's relations; exact-size native arrays; explicit publish and drop (done) | `duckvep_model_load`, `duckvep_model_drop` |
| 4b | Annotation natives | the internal natives that read a pinned model; then `duckvep_annotate_sql` and the projection builders execute on v2 | internal `_duckvep_annotate_*`, `__duckvep_projection_code` |
| 5 | Haplotype capture | caller-side normalization into spillable column collections, serial scanner, no appender | `duckvep_haplotypes`, `duckvep_coding_transcripts`, `duckvep_coding_calls` (all three share `src/duckvep_discovery.c`, which reads the model) |
| 6 | Bounded parallelism | partitioned workers, spill, cancellation and cleanup, quota accounting; 5M variants and a ten-job stress | `duckvep_native_budget`, `duckvep_native_budget_set`, `duckvep_native_budget_reset_high_water`, `duckvep_worker_limits_set` |

That is 3 + 2 + 12 + 2 + 3 + 4 = 26. The budget and worker-limit functions are
cheap to register (a BIGINT scalar and a table function), and may be ported
earlier as extra type exercises; their family is slice 6 because it is where the
limits are enforced across workers. Slice 3's builders return SQL that calls the
internal natives of slice 4, so they are verified by text equality in slice 3 and
end to end once slice 4 lands.

The 14 internal helpers (names starting `_` or `__`) are not public and not in
`docs/functions.md`: `_duckvep_revcomp`, `_duckvep_raw_gt`, `_duckvep_record_order`,
`__duckvep_projection_code` and the ten `_duckvep_annotate_{small,structural,breakend}_*`
natives. They follow the family of the function that needs them, as listed.

## Host seams

What is DuckDB-facing today, and what each file needs from a v2 host adapter.
Files not listed (`src/kernel/**`, `duckvep_reference.c`, `third_party/`,
`duckvep_model.h`, `duckvep_list.h`) are host-neutral or contain no DuckDB call.

| File | DuckDB-facing surface | What the v2 adapter must supply |
| --- | --- | --- |
| `duckvep.c` | entry point, version gate, `_duckvep_revcomp` (done), registration order | the `duckdb_v2` entry point (`duckvep_init_c_api_v2`); no version check needed (the footer is the contract); no DDL at LOAD. Done; `_duckvep_revcomp` is in `host_v2_nested.c` |
| `duckvep_registration.c` | `duckdb_query` of registration SQL, error reporting | nothing: LOAD runs no SQL in v2, so this file is not ported (the v1 build keeps it until its own cleanup) |
| `duckvep_builder.c` / `.h` | scalar-function sets, option STRUCT reading, VARCHAR out (done) | `host_v2_columns.h` option-STRUCT reader, overload registration, SQL text into an arena-backed string: `host_v2_builders.c` |
| `duckvep_sql.c` | `duckvep_so_terms`, `duckvep_annotate_sql`, `duckvep_annotate_projected_sql`, `duckvep_transcript_projection_sql` | table-function bind/global state/exec (`duckvep_so_terms` done); the three builders (done) |
| `duckvep_prepare_sql.c`, `duckvep_structural_sql.c` | five preparation builders (done) | builder adapter: `host_v2_builders.c` |
| `duckvep_ensembl.c` | four Ensembl/receipt builders (done) | builder adapter: `host_v2_builders.c` |
| `duckvep_annotate.c` | `duckvep_allele_geometry`, `duckvep_breakend_geometry` (ported), ten internal annotation natives (registry pin, `volatile`, extra info, STRUCT/LIST output) | scalar user data for the registry, STRUCT output, LIST output with child capacity, per-worker scratch |
| `duckvep_repeat.c` | `duckvep_repeat_alleles` (done) | LIST and STRUCT input, DECIMAL/HUGEINT cells, NULL-in-child propagation, option STRUCT. Done in `host_v2_nested.c` |
| `duckvep_phase_sql.c` | `duckvep_phase_call`, `_duckvep_raw_gt`, `_duckvep_record_order` (all done) | the same as repeat, plus LIST<STRUCT> output sized once per chunk. Done in `host_v2_nested.c` |
| `duckvep_model.c` | `duckvep_model_load`/`_drop` and the registry with its private connection (v1 only; the model itself is `src/core`) | done for v2 in `host_v2_model.c`: COPY-to bind/init/batch/flush/finalize, ordered merge into a column-data collection, scan sources, publish/drop/fingerprint, the statement builder |
| `duckvep_haplotype_sql.c`, `duckvep_coding_calls.c`, `duckvep_discovery.c` | `duckvep_haplotypes`, `duckvep_coding_transcripts`, `duckvep_coding_calls`: private query/fetch, statement parsing, an appender, LIST-of-STRUCT output | column-data-collection capture in place of the appender, the native scan source over the captured input, LIST-of-STRUCT output, cancellation |
| `duckvep_budget_sql.c` | four resource-control functions, `set_max_threads` | BIGINT scalar, table function with `max_threads`, volatile stability |

The shared pieces of the v2 adapter, in `host_v2/host_v2_common.h` and
`host_v2_nested.c`: an error helper, vector views with selection and validity,
string read and arena write, struct output, `register_scalar` (overloads are
repeated registrations of one name; `ANY` parameters), typed columns with
logical-type inspection (DECIMAL storage kind and scale), nested column opening,
cell reading, option-STRUCT reading and LIST<STRUCT> output. Still to grow: scalar
user data and init data, table-function named parameters and local state, COPY
callbacks, column-data-collection capture.

## The model sink on v2 (slice 4a)

v1's `duckvep_model_load(name, regions_query, transcripts_query, exons_query, ...)` runs the
caller's queries on a **private connection** inside its own transaction and stays exactly as it
is. v2 cannot run them there (a private connection inside a callback is not allowed), so the load
becomes caller-side statements plus one native call, all in the caller's own transaction:

| Step | v1 | v2 |
| --- | --- | --- |
| build | `SELECT loaded FROM duckvep_model_load(name, regions_q, transcripts_q, exons_q, mature_mirna_query := ..., peptide_edit_query := ..., interval_feature_query := ..., reference_fasta := ..., transcript_coverage_complete := ...)` | `SELECT unnest(duckvep_model_load_sql(name, regions_q, transcripts_q, exons_q [, {mature_mirna_query: ..., peptide_edit_query: ..., interval_feature_query: ..., reference_fasta: ..., transcript_coverage_complete: ...}]))` returns the statements to run, in order |
| stage | (inside the call) | one `COPY (<query>) TO 'duckvep_stage' (FORMAT duckvep_stage, MODEL '<name>', RELATION '<relation>', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)` per relation |
| publish | (inside the call) | `SELECT duckvep_model_publish('<name>' [, {reference_fasta: ..., transcript_coverage_complete: ...}])`, which returns `true` |
| drop | `SELECT duckvep_model_drop(name)` | the same; it also discards rows staged for that name |

`rduckvep_load_model(con, name, regions_query, ...)` hides the difference: it detects the host and
runs the right form, so R code is the same on both.

Differences a user can see:

- **Visibility.** On v2 the queries run in the caller's session and transaction, so TEMP tables and
  uncommitted rows are visible to them; on v1 they see only committed, permanent relations.
  (Publishing itself is not transactional: a model stays installed if the surrounding transaction rolls back,
  on both hosts. Drop it explicitly.)
- **Shape.** v1 is one table-function call with named parameters; v2 is a list of statements (from a
  scalar with an options STRUCT, like the other builders) plus the publish call. The `COPY`s are ordinary SQL and
  cannot run inside `query()`, which is why `duckvep_model_load_sql` returns text to execute instead of executing it.
- **Errors.** Validation errors are the same text (the loaders are the same code); the publish
  call reports them under `duckvep_model_publish:`.

How it works: `COPY ... (FORMAT duckvep_stage)` runs the caller's typed, ordered SELECT and merges its
ordered batches into a DuckDB column-data collection (buffer-manager backed, so it spills under the
memory limit and is counted against DuckDB's, not the native budget). A COPY that fails or is cancelled
destroys its rows and stages nothing; a finished COPY replaces earlier rows of the same relation. Publish
takes the staged relations out (whatever the outcome: a failed publish publishes nothing and consumes the
staging), scans them through the same `src/core` loaders v1 uses, allocates the exact-size native arrays
charged to the native budget (an exhausted budget is the same named capacity error), validates, builds the
kernel, and installs the model under the registry lock. The SDK exposes neither bytes nor a size hint for a
collection, so staging is bounded by its row counts (the loaders' `uint32` limits) and DuckDB's memory limit.

`_duckvep_model_fingerprint(name)` (internal, on both hosts) is an FNV-1a hash of the loaded arrays. The equality
suite loads the README model, a richer model (strands, a non-coding transcript, flanks, mature miRNA, peptide edits,
regulatory features, complete coverage) and a model with a reference FASTA on both hosts, and requires equal
fingerprints. `test/sql_v2/v2_model.sql` covers what only v2 has: TEMP and uncommitted visibility, a failed COPY,
invalid staged rows, a wrong column type, discard by drop, and option errors.

v2 does not yet have `duckvep_native_budget_set` (slice 6), so the budget refusal of a model load is tested on v1
(`make test_fault_injection` fails every allocation of a load in turn) and shared by v2 through the same loaders.

## Builders: text equality, and what runs on v2

The builders' product is SQL text, so the gate is byte equality: the equality suite compares
the md5 and length of the text each builder returns on v1 and v2 over every option key,
the defaults, NULL and invalid options, identifier quoting (quotes, schema dots, empty,
unicode, embedded NUL), the error message and its order, and multi-row, constant and
selected inputs (192 cases). The recorded results come from the v1 build *before* the text moved
to `src/core`, so the same check proves the v1 text did not change.

Which built queries execute on v2 today (slice 3), by the natives their text calls:

| Builder | Text calls | Runs on v2 now |
| --- | --- | --- |
| `duckvep_prepare_sv_geometry_sql`, `_prepare_expansionhunter_sql`, `_prepare_breakend_fusion_sql`, `_prepare_structural_hgvs_sql`, `_prepare_breakend_pairs_sql` | SQL only; `duckvep_breakend_geometry` (pairs, hgvs) | yes: run against fixtures in the suite and equal to v1 (the expansionhunter rows also feed `duckvep_repeat_alleles`) |
| `duckvep_ensembl_regions_sql`, `_ensembl_regulation_features_sql`, `duckvep_model_receipt_sql` | SQL only | yes in principle (plain SQL over caller tables); only the text is tested |
| `duckvep_ensembl_transcripts_sql` | `_duckvep_revcomp` | yes in principle; only the text is tested |
| `duckvep_annotate_sql`, `duckvep_annotate_projected_sql` | `_duckvep_annotate_*` natives | no: slice 4 (they read a pinned model) |
| `duckvep_transcript_projection_sql` | `__duckvep_projection_code` (and `duckvep_so_terms`, `_duckvep_revcomp`, `duckvep_allele_geometry`) | no: `__duckvep_projection_code` is slice 4 |

`duckvep_coding_transcripts` has no SQL-text part: it is a model lookup and moves with slice 4.

## What the v2 SDK lacks, for later slices

Read from the pinned headers and the slice-1 and slice-2 experience:

- **No appender.** `column_data_collection_{create_with_context,append}` replaces
  buffering, not insertion; haplotype capture (slice 5) removes its appender use
  and any insertion goes through caller-side `INSERT ... SELECT` or COPY.
- **No separate list-child reserve, and none is needed on the pinned revision.**
  `vector_set_size` reserves space for the new size (the memo saw it not allocate on
  another snapshot). Slice 2 sizes a list's child once per chunk and writes through
  it: 18,000 `duckvep_phase_call` records from 9,000 rows, a single row of 5,000
  slots, and lists of 5,000 elements as input all pass on the pinned build and under
  ASan. Borrowed child pointers must still be taken after the `set_size` calls.
- **No cost hint** for scalar and aggregate functions (only a cast cost exists).
  Planner behavior for the annotation natives is therefore volatile-or-not and nothing
  finer.
- **Ordered aggregates** crash upstream on the memo's snapshot; nothing in DuckVEP
  should be built on one.
- **No safe re-entry.** `parse_sql`/`statement_bind` parse and bind but do not make
  callback re-entry safe; a private connection inside a callback is still not allowed,
  which is why the model sink uses COPY.
- **Vectors are richer and stricter.** A view can be flat, constant or dictionary
  (selection); the view getter rejects FSST, sequence and shredded vectors until
  they are flattened; reading a dictionary flattens its child in place, which can
  invalidate a view taken earlier over an aliased buffer (the adapter takes views twice
  and keeps the second set). A NULL struct row needs NULL children: write it through
  `vector_set_null`, not a raw mask.
- **Type inspection costs a value per parameter.** A DECIMAL's width and scale, and a
  STRUCT's field names, come from `logical_type_get_param` as owned values that must
  be destroyed; `ANY` parameters exist only in signatures (`create_type_from_id`), so
  a function that takes `ANY` reads each argument's type at execution, as v1 does.
  Nested arguments are flattened (children too) before their children are addressed
  by list offset; there is no per-chunk bind-time type hook that is cheaper.
- **Error text is per call.** Every fallible call can return its own error handle;
  the adapter copies it into the callback's borrowed handle.
- **The CLI of the pinned snapshot swallows statement errors** in its streaming output
  modes (`-list`, `-csv`, `-json`, `-line`), reporting exit 0. The test runner uses
  `-column` for that reason.
- **The preview ABI moves between snapshots** (see Pins), so the pinned revision
  must be bumped deliberately, with the headers, checksums and CI cache key together.
