# DuckVEP scale run scale-20260929T203851Z-218970

Every job met every ceiling: **no**. Outcomes: ok 2, capacity_error 2, failed 0.

## Run

| item | value |
|---|---|
| host | Ubuntu-2404-noble-amd64-base |
| cores | 12 |
| RAM (GiB) | 62.58 |
| kernel | 6.8.0-78-generic |
| jobs x threads | 2 x 6 |
| panel | 1M (all) |
| panel rows | 1000000 |
| panel sha256 | 51af4833b5cee778e2f8bcf68ea8fc714d2f813b541829a436bf9fb4c2bd82f4 |
| model sha256 | 0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b |
| extension sha256 | ab23ed5c0f6066991db7f97c6e8fbf18fd454c6d4d5d9dc3e11ebfd843d2b5ff |
| git revision (dirty files) | 290d13a28781e065a2fe24c1e1a2a27672e5a020 (5) |
| cgroup mode | systemd |
| wall (s) | 13.750 |
| load average at start | 6.05 7.14 6.75 |

## Ceilings

Process memory.max 16.00 GiB with memory.swap.max 0; DuckDB memory_limit 8GB, max_temp_directory_size 4.00 GiB; native budget 4096 MiB, 6 workers x 128 MiB scratch + 64 MiB emit.

## Aggregate

| mode | jobs_ok | jobs | aggregate_alleles_per_s | mode_s_p50 | mode_s_max | job_wall_s_p50 | job_wall_s_max | checksums_agree |
|---|---|---|---|---|---|---|---|---|
| compact | 2 | 2 | 822707 | 2.203 | 2.431 | 13.468 | 13.709 | TRUE |
| complete17 | 0 | 2 |  |  |  | 13.468 | 13.709 | FALSE |

Aggregate alleles/s is the sum of the jobs' panel rows over the slowest job's annotation-and-write time; it excludes model load.

## Per job

| job | mode | outcome | load_s | mode_s | alleles_per_s | out_rows | out_MiB | memory_peak_GiB | oom_kill | spill_peak_MiB | native_hw_total_MiB | native_hw_model_MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | compact | ok | 3.897 | 2.431 | 411353.4 | 13518845 | 10.0 | 4.98 | 0 | 0.0 | 1377.4 | 1304.1 |
| 1 | complete17 | capacity_error | 3.897 | FALSE |  |  | NA | 4.98 | 0 | 0.0 | 1430.5 | 1304.1 |
| 2 | compact | ok | 3.624 | 2.203 | 453926.5 | 13518845 | 10.0 | 4.98 | 0 | 0.0 | 1385.5 | 1304.1 |
| 2 | complete17 | capacity_error | 3.624 | FALSE |  |  | NA | 4.98 | 0 | 0.0 | 1428.7 | 1304.1 |

## Checksums

| job | mode | out_rows | hash_sum | hash_xor |
|---|---|---|---|---|
| 1 | compact | 13518845 | 124686083948546364800189331 | 17690839170935794317 |
| 1 | complete17 |  |  |  |
| 2 | compact | 13518845 | 124686083948546364800189331 | 17690839170935794317 |
| 2 | complete17 |  |  |  |

The checksum is the row count with the sum and XOR of DuckDB `hash()` over every output row, independent of row order.

## Non-ok jobs

- job 1 complete17: **capacity_error**: Invalid Error: Invalid Input Error: capacity error: per-worker emitted-output lease budget exceeded (requested 27262976 bytes, 57683968 in use, limit 67108864) (duckvep_annotate: fused projection stream diverged) ℹ Context: rapi_execute ℹ Error type: INVALID
- job 2 complete17: **capacity_error**: Invalid Error: Invalid Input Error: capacity error: per-worker emitted-output lease budget exceeded (requested 27262976 bytes, 57683968 in use, limit 67108864) (duckvep_annotate: fused projection stream diverged) ℹ Context: rapi_execute ℹ Error type: INVALID
