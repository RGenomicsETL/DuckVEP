# GRCh38 model load trials

Input: `/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb` (SHA-256 `0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b`), attached read-only. Each trial used a fresh R process, DuckDB R 1.5.5, one DuckDB thread and an 8 GB DuckDB memory limit. `benchmarks/scale_model_load.R` loaded the model through `duckvep_model_load`; `/proc/self/status` VmHWM and `/usr/bin/time` maximum RSS agreed. Both binaries were loaded from immutable copies named `duckvep.duckdb_extension`.

| Build | Binary SHA-256 | Source | Median load (s) | Peak RSS range (KiB) |
|---|---|---|---:|---:|
| origin-main | `1c29eed6407c7e83472dcff772ef47010aeb72b54d5b1f3284ee64f4ef3bda73` | `eb7c6861c00c741896feeb9e5d4aa71b32e512c6` | 3.553 | 5,284,020–5,284,416 |
| exact-capacity | `7b114ce738735aaa94fb0df7bde61bbacd786f6ef78d52f569b6917b573bd86e` | local commit `9747bf3b1b550e91894c48e5e2a2b81ca5daa45e` | 4.656 | 5,461,964–5,462,700 |

Individual observations are in `model-load-trials.tsv`. The exact-capacity loader performs additional source scans; observed median load time increased by 1.103 s and process peak RSS increased by at least 177,548 KiB. VmHWM includes DuckDB and R allocations; it cannot substitute for model/index reservation or native owner counters. No owner-specific reserved/used measurements, 4 GiB budget admission, 32 MiB page cap or 2 GiB model/index certification is available from these trials.

## Unordered scatter loader (`scale-slice3`, commit `096b945`)

Same setup (fresh R process, one thread, 8 GB limit, immutable binary copies, a copy of the read-only model file). `origin-main-rerun` is the origin/main-equivalent binary `5e77dc2c175901afc635121b0958729e8b6a59ae880e82a4ffe5e77fb191b61e` loading the **ordered** queries; `scatter-unordered` is `b831bd915f138d2b3fb0d267663b8ddb02c3ab72c4916b56d6ed1394def0fc11` loading the same tables **without ORDER BY** (`benchmarks/scale_model_load.R EXTENSION MODEL unordered`).

| Build | Median load (s) | Peak RSS (KiB), 3 runs |
|---|---:|---:|
| origin-main-rerun | 3.438 | 5,285,120 / 5,285,404 / 5,285,564 |
| scatter-unordered | 3.299 | 5,232,852 / 5,232,888 / 5,233,200 |

Load time improved by 0.139 s (4.0%); peak RSS fell by about 52,500 KiB (1.0%, 51 MiB). The peak is not materially reduced: the stable v1 C API has no streaming result (`duckdb_execute_prepared_streaming` is unstable-only), so every query result is still fully materialized by DuckDB before rows are scanned; only the sort buffers are gone. This is consistent with the earlier 2.8 GiB pre-array figure being result materialization and the attached model buffer cache rather than the sort; that attribution is inferred, not separately measured.

## Scatter loader rebased onto circular-topology main

Same setup as above (fresh R process, one thread, 8 GB limit, immutable binary copies, copy of the read-only model). `main` is origin/main including the #19 circular validations (binary `f7a19dfca143e936b755fccc36a54ac4451b2a8c50e0eda538049c389d904101`, ordered queries). `scatter` is the rebased branch (binary `4907665be5961850d2e70c1f18c035a3a20bab38b70c9b61762d5ce87b42f813`) loading the same tables without ORDER BY. For reference, the rebased binary was also run on the ordered queries.

| Build | Load times (s), 3 runs | Median load (s) | Peak RSS (KiB), 3 runs |
|---|---|---:|---:|
| main (ordered) | 3.509 / 3.532 / 3.493 | 3.509 | 5,285,640 / 5,284,912 / 5,285,212 |
| scatter (unordered) | 3.312 / 3.377 / 3.364 | 3.364 | 5,233,072 / 5,233,328 / 5,232,960 |
| scatter (ordered queries, reference) | 3.969 / 3.966 / 3.955 | 3.966 | 6,433,940 / 6,433,404 / 6,433,416 |

The ordered-query row only shows the cost of DuckDB sorting inside the query (about 1.2 GB extra RSS); it is not a loader comparison.
