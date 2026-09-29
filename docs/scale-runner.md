# DuckVEP scale runner

`scripts/duckvep_scale_run.sh` runs N concurrent DuckVEP annotation jobs over one gnomAD panel, each in its own cgroup v2 scope with the agreed ceilings, and writes one receipt (`receipt.csv`) and one summary (`summary.md`). It annotates through the public builders (`FROM query(duckvep_annotate_sql(...))` for compact output, the projected builder plus the complete-17 field query for complete output), streams each result to Parquet, and never collects output in R.

## Ceilings (enforced, never advisory)

| Layer | Ceiling |
|---|---|
| Process | cgroup v2 `memory.max` 16 GiB, `memory.swap.max` 0, one scope per job |
| DuckDB | `memory_limit='8GB'`, job-local `temp_directory`, `max_temp_directory_size='32GiB'` |
| DuckVEP native | `duckvep_native_budget_set` 4 GiB; `duckvep_worker_limits_set` 6 workers x 128 MiB scratch + 256 MiB emit (the extension default); 6 threads (see the emit note below) |

Every job ends `ok`, `capacity_error` (an explicit DuckVEP or DuckDB capacity error such as `capacity error: native memory (model) budget exceeded` or `Out of Memory Error`) or `failed`, with the reason. Output is written to `*.partial` and renamed only after success, so a capacity error publishes nothing. A cgroup OOM kill (`memory.events` `oom_kill` > 0) is always `failed`, even if the job had reported success. The receipt's `ceilings_met` is true only when the job was `ok` and the ceilings were enforced.

Emit allowance: the extension default for the per-worker emitted-output allowance is 256 MiB (it was 64 MiB). On the gene-dense 1M panel (13.5M output rows) complete-17 output peaks at about 110 to 128 MiB of emitted buffers and fails with an explicit `capacity error: per-worker emitted-output lease budget exceeded` at 64 MiB; that run's receipt is `benchmarks/data/scale_contracts/runs/smoke-2x1M-emit64/`. The allowance is charged as it grows, so the cap reserves nothing; six workers x (128 + 256) MiB is a ceiling, not a reservation, and the measured native high-water is about 1.5 GiB of the 4 GiB budget. `--emit-mib 64` reproduces the old ceiling.

Regulation: jobs load the resident RegulatoryFeature and MotifFeature intervals from the model's `duckvep_regulation_features` (the production configuration) by default and record their count in the receipt and summary. `--no-regulation` runs without them.

## The consumer's run (10 jobs x 5M alleles, 64 cores, 256 GB)

Prerequisites: Linux with cgroup v2 mounted at `/sys/fs/cgroup`; root (or a systemd that lets you `systemd-run --scope -p MemoryMax=...`); R with the `DBI` and `duckdb` packages (DuckDB 1.5 or later); `sha256sum`, `awk`, `du`; a built extension (`make release`, giving `build/release/extension/duckvep/duckvep.duckdb_extension`) or a released one; the GRCh38 model `homo_sapiens_116_GRCh38_final.duckdb`; the 5M panel Parquet; at least 10 x 32 GiB + 10 GiB free disk on the spill filesystem; about 176 GiB of RAM for the ceilings (10 x 16 GiB) plus the OS, and 60 cores (10 x 6 threads).

```sh
scripts/duckvep_scale_run.sh --jobs 10 --panel /data/genomes-5M.parquet \
  --model /data/homo_sapiens_116_GRCh38_final.duckdb \
  --extension build/release/extension/duckvep/duckvep.duckdb_extension \
  --threads-per-job 6 --out /scratch/duckvep-scale-10x5M
```

`--panel` takes a file path or a name (`1M`, `5M`, `25M`, `100M`, `exome2M`, `structural`) that resolves through `$DUCKVEP_GNOMAD_ROOT/panels/`. The consumer has no staged gnomAD tree, so pass the path of the `panel-5000000.parquet` we send. Check its SHA-256 against `benchmarks/data/scale_contracts/panels/panel-receipts.tsv` (`genomes-5M`, `file_sha256`); the runner records the SHA-256 of the panel, model and extension it used.

The pre-flight refuses, with the reason, when the host cannot hold the jobs: RAM below jobs x 16 GiB + 4 GiB, jobs x threads above the core count, free disk below jobs x spill quota + `--min-free-gib` (10), no cgroup v2, or no usable `systemd-run`. `--oversubscribe` overrides the capacity checks only; it never disables enforcement. `--cgroup manual` writes `cgroup.procs` directly on hosts without `systemd-run`; `--cgroup none` needs `--allow-unenforced` and the receipt and summary then say "ceilings not enforced".

To send back: the whole `--out` directory except spill leftovers, or at minimum `receipt.csv`, `summary.md` and `run.tsv`, plus `job-*/job.log` if any job is not `ok`. `run.tsv` carries the host, kernel, core and RAM counts, load average, panel/model/extension SHA-256, git revision and every ceiling that was applied.

## What is recorded

Per job and mode (`compact`, `complete17`): model load seconds, annotation-and-write seconds, alleles/s, output rows and bytes, a full-output checksum (row count, sum and XOR of DuckDB `hash()` over every output row, independent of row order), cgroup `memory.peak` and `memory.events`, the peak size of the job's spill directory (sampled every 0.25 s), the process peak RSS, the native budget high-water per owner after load and during each mode (`duckvep_native_budget()`), and the outcome.

Aggregate per mode: alleles/s (sum of panel rows over the slowest job's annotation time, excluding load), p50 and max per-job latency, whether the checksums of all jobs agree, and whether every job met every ceiling.

Exit status: 0 all jobs ok with ceilings enforced; 2 usage error or pre-flight refusal; 3 some job was not ok (or ceilings not enforced).

## Self-test

```sh
scripts/duckvep_scale_run.sh --self-test --panel 1M --out /tmp/duckvep-scale-selftest
```

Runs four small jobs and checks the classification: native budget 64 MiB gives a capacity error, publishes nothing, leaves no model or index bytes charged and the same connection annotates again after the budget is restored; DuckDB `memory_limit` 100MB with a 1 MiB spill quota gives DuckDB's explicit out-of-memory error; a cgroup `MemoryMax` of 700M (far above the bare process, below the model load) is classified `failed` with `oom_kill` counted, never `ok`.

## Panels

`scripts/gnomad_build_panels.sh` rebuilds every panel in one command from the staged gnomAD v4.1 lean Parquet: distinct 1M/5M/25M/100M genome panels, the 2M coding/splice exome panel (cds_snv 1M, cds_indel 600k, splice_remainder 400k; gnomAD sites contain no MNVs, so MNV coverage comes from the existing ClinVar and indel/MNV controls, not from gnomAD) and the 12,000-row structural controls. Panels are data under `$DUCKVEP_GNOMAD_ROOT/panels/` and are never committed; `benchmarks/data/scale_contracts/panels/panel-receipts.tsv` (rows, file SHA-256, content checksum, chromosomes, sets, source-part count and source-manifest SHA-256) and `panel-shards.tsv` (per source shard) are. A panel whose quota exceeds the distinct alleles staged so far is recorded as `insufficient_distinct_alleles` until more genome shards land; run the same command again then. DuckDB spill for the builders goes to `$DUCKVEP_SCALE_TMP` (default `/root/duckvep/data/scale-tmp`), is removed at exit and is capped so that 50 GiB stay free.
