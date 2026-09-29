# Fused presentation timing

From `/root/duckvep/wt-scale`, run `taskset -c 2 Rscript benchmarks/scale_fusion_timing.R <checkout> <bench|gnomad1m>` for each checkout and corpus. Reference checkout: `/root/duckvep/DuckVEP` at `013ce4fdc0134b7247614a8515f1c908165d85ea` (locally rebuilt release binary). Fused checkout: `scale-slice2`. The model and canonical snapshot are the frozen GRCh38 scale-contract inputs. The panel is `panels/genomes-v1-77d2cdff65171780/panel-1000000.parquet` under `/root/duckvep/data/gnomad-v4.1` (SHA-256 `da912d30ed038cf38f6b84b8f817611e3b008e307d8a61d19f0aa6f4a2073b01`). DuckDB R 1.5.5 ran with one DuckDB thread, one pinned CPU core, an 8 GB memory limit and five replicates per lane. The source events, metadata snapshot and input relation were staged before measurement.

The load scope times `duckvep_model_load` from the attached read-only model, excluding event staging and canonical-dimension preparation. Annotation times a `count(*)` over the annotation cursor and `UNNEST` (including its input validation); this is not the full presentation plan. Total times a `COPY` of the formatted result to a fresh local uncompressed Parquet file. Formatting-plus-sink is the paired difference, total minus annotation, not an isolated operator timer; plan pruning, input sorting, materialization, caching and output-order-dependent Parquet encoding make this an approximate attribution rather than a separately runnable presentation cost. Each output file was deleted between repetitions. The production complete-17 query streams directly from the annotation cursor to formatting and `COPY` without materializing a rich-result table.

Median wall seconds (five repetitions per checkout and corpus; raw timings in `fusion-timings.tsv`):

| Corpus (inputs → outputs) | Contract | Checkout | Load | Annotation | Formatting + sink (difference) | Total COPY |
|---|---|---|---:|---:|---:|---:|
| bench (100,957 → 1,174,245) | compact | reference | 2.618 | 0.144 | 0.142 | 0.289 |
| bench (100,957 → 1,174,245) | compact | fused | 2.610 | 0.146 | 0.143 | 0.289 |
| bench (100,957 → 1,174,245) | complete-17 | reference | 2.618 | 0.476 | 1.240 | 1.716 |
| bench (100,957 → 1,174,245) | complete-17 | fused | 2.610 | 0.483 | 1.177 | 1.659 |
| gnomAD 1M (1,000,000 → 4,224,817) | compact | reference | 2.622 | 0.602 | 0.403 | 1.003 |
| gnomAD 1M (1,000,000 → 4,224,817) | compact | fused | 2.541 | 0.599 | 0.400 | 0.999 |
| gnomAD 1M (1,000,000 → 4,224,817) | complete-17 | reference | 2.622 | 1.458 | 3.312 | 4.779 |
| gnomAD 1M (1,000,000 → 4,224,817) | complete-17 | fused | 2.541 | 1.475 | 3.237 | 4.708 |

For the bench corpus, reference and fused full-row DuckDB hash sums and XORs match in both contracts: compact `10832993039152099474950971` / `15844975574530270781`, complete-17 `10829798798206598272015142` / `6196530389359083174`. For gnomAD 1M they also match: compact `38967796526644455110324607` / `17777311707198045539`, complete-17 `38979161771609265969412183` / `14456585455263769313`. The frozen bench gate separately checks full-row partition SHA-256 digests.

On the million panel, the materialized projected plan added a second million-row `ORDER_BY` (0.246 s), a source `CTE` (0.046 s), a second source scan and a validation `FILTER` (0.031 s) in a detailed single-run DuckDB operator profile. The reference plan had one input `ORDER_BY` (0.092 s); its `UNNEST` took 1.551 s against 1.757 s in the materialized projected plan. The ordered input view supplies that ordering without another sort. The scalar validates literal REF/ALT and coordinates; the builder rejects empty model names and guards missing event identity and the special `<*>` allele. The guarded plan has one `ORDER_BY` (0.088 s), no materialized source or validation filter, and `UNNEST` at 1.407 s in a separate single-run profile. The instrumented-binary comparison found one `ORDER_BY` (0.085 s reference; 0.087 s fused) and `UNNEST` (1.427 s; 1.422 s).

Temporary native `CLOCK_MONOTONIC` timers on one instrumented extension found 3,643,164 projected pairs per path, 3,553,217 outside exons/UTRs/CDS. Cursor fill took 0.798 s reference and 0.784 s fused; sampling every 128th projected pair measured 28,462 projection stores at 0.003082 s and 0.003051 s respectively. Both paths request the same projected facts, including noncoding hits, so projection itself does not explain the SQL-plan regression. The complete-17 wall medians are below the reference on both corpora; operator timings and wall times vary with caching and scheduling.
