# GRCh38 model load trials

Input: `/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb` (SHA-256 `0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b`), attached read-only. Each trial used a fresh R process, DuckDB R 1.5.5, one DuckDB thread and an 8 GB DuckDB memory limit. `benchmarks/scale_model_load.R` loaded the model through `duckvep_model_load`; `/proc/self/status` VmHWM and `/usr/bin/time` maximum RSS agreed. Both binaries were loaded from immutable copies named `duckvep.duckdb_extension`.

| Build | Binary SHA-256 | Source | Median load (s) | Peak RSS range (KiB) |
|---|---|---|---:|---:|
| origin-main | `1c29eed6407c7e83472dcff772ef47010aeb72b54d5b1f3284ee64f4ef3bda73` | `eb7c6861c00c741896feeb9e5d4aa71b32e512c6` | 3.553 | 5,284,020–5,284,416 |
| exact-capacity | `7b114ce738735aaa94fb0df7bde61bbacd786f6ef78d52f569b6917b573bd86e` | local commit `9747bf3b1b550e91894c48e5e2a2b81ca5daa45e` | 4.656 | 5,461,964–5,462,700 |

Individual observations are in `model-load-trials.tsv`. The exact-capacity loader performs additional source scans; observed median load time increased by 1.103 s and process peak RSS increased by at least 177,548 KiB. VmHWM includes DuckDB and R allocations; it cannot substitute for model/index reservation or native owner counters. No owner-specific reserved/used measurements, 4 GiB budget admission, 32 MiB page cap or 2 GiB model/index certification is available from these trials.
