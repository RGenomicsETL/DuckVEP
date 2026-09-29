# GRCh38 model-load attribution

Instrumented extension built from `7291608acaf61ea2094e4adf146180ea59da2666` plus the trace-only hooks in `src/duckvep_model.c`. Binary SHA-256: `9c666b9c05c185f0224f28ad16456ce50d389eb98bcc9fd4593f3ba2b3986b01`. Each of three fresh R processes loaded an immutable copy with `DUCKVEP_LOAD_TRACE=1 Rscript benchmarks/scale_model_load.R /tmp/scale3-attrib/duckvep.duckdb_extension /root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb`. Times include the entire load; RSS is `/proc/self/status` VmRSS/VmHWM, KiB. The trace rows are ordered, with UTC timestamp, RSS, high-water RSS, estimated primary used bytes, primary reserved bytes, and interval-array reserved bytes.

| Stage/owner | Median elapsed from start (s) | Median RSS (GiB) | Median high-water (GiB) | Ownership estimate |
| --- | ---: | ---: | ---: | --- |
| Start (R + DuckDB + attached model) | 0 | 0.107 | 0.107 | baseline |
| Pre-count (including materializing ordered transcript subquery) | 1.54 | 2.80 | 2.80 | query execution / DuckDB |
| Arrays reserved, before touching pages | 1.54 | 2.80 | 2.80 | native reserved 1.20 GiB, not yet resident |
| Region and transcript query materialized | 2.29 | 4.13 | 4.13 | DuckDB query materialization overlaps arrays |
| Exon query materialized | 3.14 | 4.04 | 5.21 | high-water overlaps materialization with resident arrays |
| All primary arrays populated | 3.30 | 4.12 | 5.21 | native used/reserved 1.20/1.20 GiB |
| Kernel constructed | 4.33 | 4.12 | 5.21 | projection caches not separately accounted |
| Index built | 4.36 | 4.12 | 5.21 | cgranges interval array 13.7 MiB plus small metadata |
| Load completed | 4.39 | 4.12 | 5.21 | reference reader cache 0 (no FASTA requested) |

Medians are rounded; the three complete load times were 4.424, 4.390, 4.289 seconds and peaks 5,460,604, 5,460,008, 5,461,400 KiB. The difference between RSS and attributed native allocations is **not** a precise DuckDB allocation census: it includes allocator retention, R, the DuckDB buffer manager and transient query results. The 2.8 GiB before reserving native arrays is strong evidence that ordered query execution, especially the count over the transcript query, dominates array slack; the peak is not attributable to the 13.7 MiB index or reference cache. The counters describe only primary arrays and cgranges interval backing arrays; they exclude cgranges metadata and kernel projection caches. Snapshot of the uninstrumented parent remains in `../model-load-trials.md`.
