# DuckVEP scale run scale-20260929T205000Z-275080

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
| regulation (resident features) | on (1383580) |
| panel sha256 | 51af4833b5cee778e2f8bcf68ea8fc714d2f813b541829a436bf9fb4c2bd82f4 |
| model sha256 | 0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b |
| extension sha256 | fb567cbca680237ad5402a31357b9e3ee2566e80adc5d917413783c8c2704480 |
| git revision (dirty files) | 0691d2d0b5b02cc86d6ae4f19b71adbbf74a936a (0) |
| cgroup mode | systemd |
| wall (s) | 17.359 |
| load average at start | 6.19 5.91 6.15 |

## Ceilings

Process memory.max 16.00 GiB with memory.swap.max 0; DuckDB memory_limit 8GB, max_temp_directory_size 4.00 GiB; native budget 4096 MiB, 6 workers x 128 MiB scratch + 256 MiB emit.

## Aggregate

| mode | jobs_ok | jobs | aggregate_alleles_per_s | mode_s_p50 | mode_s_max | job_wall_s_p50 | job_wall_s_max | checksums_agree |
|---|---|---|---|---|---|---|---|---|
| compact | 2 | 2 | 842460 | 2.243 | 2.374 | 16.423 | 17.318 | TRUE |
| complete17 | 2 | 2 | 228781 | 7.885 | 8.742 | 16.423 | 17.318 | TRUE |

Aggregate alleles/s is the sum of the jobs' panel rows over the slowest job's annotation-and-write time; it excludes model load.

## Per job

| job | mode | outcome | load_s | mode_s | alleles_per_s | out_rows | out_MiB | memory_peak_GiB | oom_kill | spill_peak_MiB | native_hw_total_MiB | native_hw_model_MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | compact | ok | 3.976 | 2.374 | 421230 | 13588086 | 10.7 | 4.98 | 0 | 0.0 | 1433.6 | 1326.1 |
| 1 | complete17 | ok | 3.976 | 8.742 | 114390.3 | 13588086 | 44.1 | 4.98 | 0 | 0.0 | 1525.3 | 1326.1 |
| 2 | compact | ok | 3.904 | 2.243 | 445831.5 | 13588086 | 10.7 | 5.03 | 0 | 0.0 | 1435.2 | 1326.1 |
| 2 | complete17 | ok | 3.904 | 7.885 | 126823.1 | 13588086 | 44.0 | 5.03 | 0 | 0.0 | 1505.2 | 1326.1 |

## Checksums

| job | mode | out_rows | hash_sum | hash_xor |
|---|---|---|---|---|
| 1 | compact | 13588086 | 125326712289080813184522349 | 14348256641605856663 |
| 1 | complete17 | 13588086 | 125339144276500816740216863 | 10632275200237930663 |
| 2 | compact | 13588086 | 125326712289080813184522349 | 14348256641605856663 |
| 2 | complete17 | 13588086 | 125339144276500816740216863 | 10632275200237930663 |

The checksum is the row count with the sum and XOR of DuckDB `hash()` over every output row, independent of row order.
