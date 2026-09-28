# Fused presentation timing

From `/root/duckvep/wt-scale`, run `taskset -c 2 Rscript benchmarks/scale_fusion_timing.R <checkout> <bench|gnomad1m>` for each checkout and corpus. Reference checkout: `/root/duckvep/DuckVEP` at `0ceeeba4adeece46a44f333f71de77922ab10447` (locally rebuilt release binary). Fused checkout: `scale-slice2`. The model and canonical snapshot are the frozen GRCh38 scale-contract inputs. The panel is `panels/genomes-v1-77d2cdff65171780/panel-1000000.parquet` under `/root/duckvep/data/gnomad-v4.1` (SHA-256 `da912d30ed038cf38f6b84b8f817611e3b008e307d8a61d19f0aa6f4a2073b01`). DuckDB R 1.5.5 ran with one DuckDB thread, one pinned CPU core, an 8 GB memory limit and five replicates per lane. The source events, metadata snapshot and input relation were staged before measurement.

The load scope times `duckvep_model_load` from the attached read-only model, excluding event staging and canonical-dimension preparation. Annotation times a `count(*)` over the annotation cursor and `UNNEST` (including its input validation); this is not the full presentation plan. Total times a `COPY` of the formatted result to a fresh local uncompressed Parquet file. Formatting-plus-sink is the paired difference, total minus annotation, not an isolated operator timer; plan pruning, input sorting, materialization, caching and output-order-dependent Parquet encoding make this an approximate attribution rather than a separately runnable presentation cost. Each output file was deleted between repetitions. The production complete-17 query streams directly from the annotation cursor to formatting and `COPY` without materializing a rich-result table.

Median wall seconds (five paired repetitions; raw timings in `fusion-timings.tsv`):

| Corpus (inputs → outputs) | Contract | Checkout | Load | Annotation | Formatting + sink (difference) | Total COPY |
|---|---|---|---:|---:|---:|---:|
| bench (100,957 → 1,174,245) | compact | reference | 2.767 | 0.143 | 0.146 | 0.287 |
| bench (100,957 → 1,174,245) | compact | fused | 2.648 | 0.144 | 0.149 | 0.293 |
| bench (100,957 → 1,174,245) | complete-17 | reference | 2.767 | 0.472 | 1.244 | 1.716 |
| bench (100,957 → 1,174,245) | complete-17 | fused | 2.648 | 0.571 | 0.989 | 1.560 |
| gnomAD 1M (1,000,000 → 4,224,817) | compact | reference | 2.627 | 0.596 | 0.407 | 0.997 |
| gnomAD 1M (1,000,000 → 4,224,817) | compact | fused | 2.612 | 0.607 | 0.405 | 1.012 |
| gnomAD 1M (1,000,000 → 4,224,817) | complete-17 | reference | 2.627 | 1.456 | 3.347 | 4.789 |
| gnomAD 1M (1,000,000 → 4,224,817) | complete-17 | fused | 2.612 | 2.058 | 3.135 | 5.256 |

For the bench corpus, reference and fused full-row DuckDB hash sums and XORs match in both contracts: compact `10832993039152099474950971` / `15844975574530270781`, complete-17 `10829798798206598272015142` / `6196530389359083174`. For gnomAD 1M they also match: compact `38967796526644455110324607` / `17777311707198045539`, complete-17 `38979161771609265969412183` / `14456585455263769313`. The frozen bench gate separately checks full-row partition SHA-256 digests. The complete-17 total improves on bench but regresses on the million-panel workload; the added validation/source materialization and ordering in the projected builder are plausible contributors, not isolated causal estimates. The panel regression remains an optimization target rather than a claim of scale speedup.
