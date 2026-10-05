# DuckVEP scale runner

`scripts/duckvep_scale_run.sh` runs N concurrent DuckVEP annotation jobs over one input panel, each in its own cgroup v2 scope with the configured ceilings, and writes one receipt (`receipt.csv`) and one summary (`summary.md`). It annotates through the public builders (`FROM query(duckvep_annotate_sql(...))` for compact output, the projected builder plus the complete-17 field query for complete output), streams each result to Parquet, and never collects output in R.

## Ceilings (enforced, never advisory)

| Layer | Ceiling |
|---|---|
| Process | cgroup v2 `memory.max` 16 GiB, `memory.swap.max` 0, one scope per job |
| DuckDB | `memory_limit='8GB'`, job-local `temp_directory`, `max_temp_directory_size='32GiB'` |
| DuckVEP native | `duckvep_native_budget_set` 4 GiB; `duckvep_worker_limits_set` 6 workers x 128 MiB scratch + 256 MiB emit (the extension default); 6 threads (see the emit note below) |

Every job ends `ok`, `capacity_error` (an explicit DuckVEP or DuckDB capacity error such as `capacity error: native memory (model) budget exceeded` or `Out of Memory Error`) or `failed`, with the reason. Output is written to `*.partial` and renamed only after success, so a capacity error publishes nothing. A cgroup OOM kill (`memory.events` `oom_kill` > 0) is always `failed`, even if the job had reported success. The receipt's `ceilings_met` is true only when the job was `ok` and the ceilings were enforced.

Emit allowance: the extension default for the per-worker emitted-output allowance is 256 MiB. In the gene-dense 1M panel run (13.5M output rows), complete-17 output used about 110 to 128 MiB of emitted buffers and a 64 MiB cap produced `capacity error: per-worker emitted-output lease budget exceeded`; receipt: `benchmarks/data/scale_contracts/runs/smoke-2x1M-emit64/`. `--emit-mib 64` reproduces that measured failure. The allowance is charged as it grows, so the cap reserves nothing; six workers x (128 + 256) MiB is a ceiling, not a reservation. Measured native high-water was about 1.5 GiB of the 4 GiB budget.

Regulation: jobs load the resident RegulatoryFeature and MotifFeature intervals from the model's `duckvep_regulation_features` (the production configuration) by default and record their count in the receipt and summary. `--no-regulation` runs without them.

## Certification and frozen artifacts

A run certifies (exit 0) only when all of the following hold; anything else exits 3 with the reason in `summary.md` (for example `checksums disagree in complete17`, `native_hw_total missing for job 2`, `swap not disabled`, `model changed during the run`):

- every job and mode is `ok`;
- within each mode, the checksums (rows, `hash_sum`, `hash_xor`) agree across all ok jobs, whenever there are two or more;
- every mandatory counter is present in every row: `memory.peak`, `memory.max`, `memory.swap.max` (which must be 0), `oom_kill`, native high-water total, spill peak and mode seconds;
- the ceilings were enforced, and the model, panel and extension were unchanged from start to end.

A capacity error during certification means the qualification failed; it is never retried into a pass.

Frozen artifacts: at start the runner copies the extension into `OUT/artifacts/` (read-only) and loads only that copy, so rebuilding or replacing the original mid-campaign cannot crash running jobs. A path under a build directory is refused unless `--allow-live-extension` is given. The SHA-256 of the extension, model and panel are recorded in `run.tsv` at start and re-hashed once at the end (`*_sha256_end`; not per job); a difference fails the run. All three appear in `summary.md`.

Disk pre-flight: jobs x spill quota, plus an output budget of jobs x (estimated output bytes per mode, summed over modes) x 2 for the `.partial` copy x 1.5 margin, plus `--min-free-gib`. The estimate is the panel row count (or `--limit-rows`) times 12 B/row (compact) and 48 B/row (complete17), rounded up from the committed smoke receipts (`runs/smoke-2x1M-regulation`: 11.2 and 46.2 B/row). Outputs are deleted after hashing without `--keep-output`, but they coexist while jobs run, so they are budgeted either way. Spill lives in a dedicated `$SPILL_PARENT/<run-id>` per campaign; only that directory is removed.

## Run

A certified run requires Linux with cgroup v2 and either usable `systemd-run` or permission to write to the cgroup. R with `DBI` and `duckdb` (DuckDB 1.5 or later), `sha256sum`, `awk` and `du` are also required. `--cgroup none` requires `--allow-unenforced` and cannot produce a certified run. Supply readable paths to the panel, model and extension; the output directory must be new or empty.

```sh
scripts/duckvep_scale_run.sh --jobs 2 --panel /data/panel.parquet \
  --model /data/model.duckdb --extension /data/duckvep.duckdb_extension \
  --threads-per-job 6 --out /scratch/duckvep-scale
```

`--panel` accepts a Parquet path or a name (`1M`, `5M`, `25M`, `100M`, `exome2M`, `structural`) resolved through `$DUCKVEP_GNOMAD_ROOT/panels/`. The runner records SHA-256 values for the panel, model and extension at the start and end; a change during the run fails certification. It copies the extension into `OUT/artifacts/` and loads that copy. Extensions under a build directory require `--allow-live-extension`; the loaded copy remains fixed for the campaign.

Pre-flight checks available RAM against jobs x process memory ceiling plus 4 GiB, total job threads against available cores, free space against spill quotas and the output budget plus `--min-free-gib` (default 10), and cgroup support. `--oversubscribe` overrides host-capacity checks only; the configured ceilings remain enforced. `--cgroup manual` uses direct cgroup writes. The receipt and summary mark `--cgroup none` runs as unenforced.

`run.tsv` records the host, kernel, core and RAM counts, load average, input hashes, git revision and applied ceilings. Keep `receipt.csv`, `summary.md`, `run.tsv`, and `job-*/job.log` for any job that is not `ok` when reviewing a run.

## What is recorded

Per job and mode (`compact`, `complete17`): model load seconds, annotation-and-write seconds, alleles/s, output rows and bytes, a full-output checksum (row count, sum and XOR of DuckDB `hash()` over every output row, independent of row order), cgroup `memory.peak` and `memory.events`, the peak size of the job's spill directory (sampled every 0.25 s), the process peak RSS, the native budget high-water per owner after load and during each mode (`duckvep_native_budget()`), and the outcome.

Aggregate per mode: alleles/s (sum of panel rows over the slowest job's annotation time, excluding load), p50 and max per-job latency, whether the checksums of all jobs agree, and whether every job met every ceiling.

Exit status: 0 all jobs ok with ceilings enforced; 2 usage error or pre-flight refusal; 3 some job was not ok (or ceilings not enforced).

## Self-test

```sh
scripts/duckvep_scale_run.sh --self-test --panel 1M --out /tmp/duckvep-scale-selftest
```

First runs the aggregator guards (`scripts/duckvep_scale_selftest_aggregate.sh`, also `make test_scale_aggregate`), which feed synthetic receipts to the aggregator and need no DuckVEP run: a clean set exits 0; disagreeing checksums, a missing native high-water or other counter, `memory.swap.max` not 0, `oom_kill` above 0, a changed model hash, missing end hashes and unenforced ceilings each exit 3 with their reason. Then runs four small jobs and checks the classification: native budget 64 MiB gives a capacity error, publishes nothing, leaves no model or index bytes charged and the same connection annotates again after the budget is restored; DuckDB `memory_limit` 100MB with a 1 MiB spill quota gives DuckDB's explicit out-of-memory error; a cgroup `MemoryMax` of 700M (far above the bare process, below the model load) is classified `failed` with `oom_kill` counted, never `ok`.

## Panels

`scripts/gnomad_build_panels.sh` rebuilds every panel in one command from the staged gnomAD v4.1 lean Parquet: distinct 1M/5M/25M/100M genome panels, the 2M coding/splice exome panel (cds_snv 1M, cds_indel 600k, splice_remainder 400k; gnomAD sites contain no MNVs, so MNV coverage comes from the existing ClinVar and indel/MNV controls, not from gnomAD) and the 12,000-row structural controls. Panels are data under `$DUCKVEP_GNOMAD_ROOT/panels/` and are never committed; `benchmarks/data/scale_contracts/panels/panel-receipts.tsv` (rows, file SHA-256, content checksum, chromosomes, sets, source-part count and source-manifest SHA-256) and `panel-shards.tsv` (per source shard) are. A panel whose quota exceeds the distinct alleles staged so far is recorded as `insufficient_distinct_alleles` until more genome shards land; run the same command again then. DuckDB spill for the builders goes to `$DUCKVEP_SCALE_TMP` (default `/root/duckvep/data/scale-tmp`), is removed at exit and is capped so that 50 GiB stay free.
