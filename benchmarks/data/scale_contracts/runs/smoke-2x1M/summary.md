# DuckVEP scale run scale-20260929T204055Z-242740

Every job met every ceiling: **yes**. Outcomes: ok 4, capacity_error 0, failed 0.

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
| wall (s) | 14.591 |
| load average at start | 6.97 7.43 6.92 |

## Ceilings

Process memory.max 16.00 GiB with memory.swap.max 0; DuckDB memory_limit 8GB, max_temp_directory_size 4.00 GiB; native budget 4096 MiB, 6 workers x 128 MiB scratch + 128 MiB emit.

## Aggregate

| mode | jobs_ok | jobs | aggregate_alleles_per_s | mode_s_p50 | mode_s_max | job_wall_s_p50 | job_wall_s_max | checksums_agree |
|---|---|---|---|---|---|---|---|---|
| compact | 2 | 2 | 1058201 | 1.872 | 1.890 | 14.252 | 14.552 | TRUE |
| complete17 | 2 | 2 |  270599 | 7.072 | 7.391 | 14.252 | 14.552 | TRUE |

Aggregate alleles/s is the sum of the jobs' panel rows over the slowest job's annotation-and-write time; it excludes model load.

## Per job

| job | mode | outcome | load_s | mode_s | alleles_per_s | out_rows | out_MiB | memory_peak_GiB | oom_kill | spill_peak_MiB | native_hw_total_MiB | native_hw_model_MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | compact | ok | 3.224 | 1.89 | 529100.5 | 13518845 | 10.0 | 5.10 | 0 | 0.0 | 1380.9 | 1304.1 |
| 1 | complete17 | ok | 3.224 | 7.072 | 141402.7 | 13518845 | 42.2 | 5.10 | 0 | 0.0 | 1507.0 | 1304.1 |
| 2 | compact | ok | 3.249 | 1.872 | 534188 | 13518845 | 9.9 | 5.09 | 0 | 0.0 | 1383.6 | 1304.1 |
| 2 | complete17 | ok | 3.249 | 7.391 | 135299.7 | 13518845 | 42.3 | 5.09 | 0 | 0.0 | 1476.8 | 1304.1 |

## Checksums

| job | mode | out_rows | hash_sum | hash_xor |
|---|---|---|---|---|
| 1 | compact | 13518845 | 124686083948546364800189331 | 17690839170935794317 |
| 1 | complete17 | 13518845 | 124697337108154286232220601 | 7186442940544995373 |
| 2 | compact | 13518845 | 124686083948546364800189331 | 17690839170935794317 |
| 2 | complete17 | 13518845 | 124697337108154286232220601 | 7186442940544995373 |

The checksum is the row count with the sum and XOR of DuckDB `hash()` over every output row, independent of row order.
