# DuckVEP scale run scale-20260930T141326Z-1334504

Certified: **yes**. Outcomes: ok 6, capacity_error 0, failed 0.

## Run

| item | value |
|---|---|
| host | Ubuntu-2404-noble-amd64-base |
| cores | 18 |
| RAM (GiB) | 62.58 |
| kernel | 6.8.0-78-generic |
| jobs x threads | 3 x 6 |
| panel | 5M (all) |
| panel rows | 5000000 |
| regulation (resident features) | on (1383580) |
| panel sha256 (start / end) | 2d68aefab5e1d425e192c97f28bc02114fca4d19a6ee30a4e3a216aa722b7c6c / 2d68aefab5e1d425e192c97f28bc02114fca4d19a6ee30a4e3a216aa722b7c6c |
| model sha256 (start / end) | 0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b / 0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b |
| extension sha256 (start / end) | c914483a5ef7c44f2b8878fa3e79c8c2ab89827bc2522a33fc11fa961eeead9c / c914483a5ef7c44f2b8878fa3e79c8c2ab89827bc2522a33fc11fa961eeead9c |
| git revision (dirty files) | d03b97aeb79b05feb77e57a839123b0e05f82273 (0) |
| cgroup mode | systemd |
| wall (s) | 54.256 |
| load average at start | 3.04 1.76 1.53 |

## Ceilings

Process memory.max 16.00 GiB with memory.swap.max 0; DuckDB memory_limit 8GB, max_temp_directory_size 4.00 GiB; native budget 4096 MiB, 6 workers x 128 MiB scratch + 256 MiB emit.

## Aggregate

| mode | jobs_ok | jobs | aggregate_alleles_per_s | mode_s_p50 | mode_s_max | job_wall_s_p50 | job_wall_s_max | checksums_agree |
|---|---|---|---|---|---|---|---|---|
| compact | 3 | 3 | 1406206 | 10.303 | 10.667 | 53.357 | 54.214 | TRUE |
| complete17 | 3 | 3 |  429923 | 33.862 | 34.890 | 53.357 | 54.214 | TRUE |

Aggregate alleles/s is the sum of the jobs' panel rows over the slowest job's annotation-and-write time; it excludes model load.

## Per job

| job | mode | outcome | load_s | mode_s | alleles_per_s | out_rows | out_MiB | memory_peak_GiB | oom_kill | spill_peak_MiB | native_hw_total_MiB | native_hw_model_MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | compact | ok | 2.805 | 10.667 | 468735.4 | 68104887 | 38.0 | 6.67 | 0 | 0.0 | 1482.3 | 1326.1 |
| 1 | complete17 | ok | 2.805 | 34.89 | 143307.5 | 68104887 | 179.3 | 6.67 | 0 | 0.0 | 1550.2 | 1326.1 |
| 2 | compact | ok | 2.684 | 10.293 | 485767 | 68104887 | 37.4 | 6.65 | 0 | 0.0 | 1473.7 | 1326.1 |
| 2 | complete17 | ok | 2.684 | 33.746 | 148165.7 | 68104887 | 178.4 | 6.65 | 0 | 0.0 | 1565.9 | 1326.1 |
| 3 | compact | ok | 2.606 | 10.303 | 485295.5 | 68104887 | 38.1 | 6.67 | 0 | 0.0 | 1501.9 | 1326.1 |
| 3 | complete17 | ok | 2.606 | 33.862 | 147658.1 | 68104887 | 178.9 | 6.67 | 0 | 0.0 | 1545.5 | 1326.1 |

## Checksums

| job | mode | out_rows | hash_sum | hash_xor |
|---|---|---|---|---|
| 1 | compact | 68104887 | 628201490641970538314146085 | 575190312610025655 |
| 1 | complete17 | 68104887 | 628184535783468984172792541 | 15426324378123767487 |
| 2 | compact | 68104887 | 628201490641970538314146085 | 575190312610025655 |
| 2 | complete17 | 68104887 | 628184535783468984172792541 | 15426324378123767487 |
| 3 | compact | 68104887 | 628201490641970538314146085 | 575190312610025655 |
| 3 | complete17 | 68104887 | 628184535783468984172792541 | 15426324378123767487 |

The checksum is the row count with the sum and XOR of DuckDB `hash()` over every output row, independent of row order.
