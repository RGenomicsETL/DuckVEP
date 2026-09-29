# Native budget and bounded execution: evidence

Measured with `benchmarks/scale_budget.R` (one pinned core, `taskset -c 2`; DuckDB R 1.5.5; one DuckDB thread; 8 GB `memory_limit`; immutable copies of each binary named `duckvep.duckdb_extension`; the read-only GRCh38 model and the 1,000,000-allele gnomAD panel `genomes-v1-77d2cdff65171780/panel-1000000.parquet`). `base` is `origin/main` at 42cc36f (circular lifted execution, same-codon classifier, species SV builders) built from a `git archive` of that commit; `budget` is this branch rebased onto it. Runs alternate base, budget; each row is the median of five.

| Build | Binary SHA-256 |
|---|---|
| base | `e697a290cfec3782b2d8791e717df20b14aa345d1fc4280a1915a3f1c07a2acd` |
| budget | `a45c354382b0901c465ae55806887566f1238833060d84f8fbddf25856d74f14` |

## Time and process memory

| Measure | base | budget | Difference |
|---|---:|---:|---:|
| Model load (s), runs 3.419 / 3.653 / 3.360 / 3.359 / 3.424 vs 3.506 / 3.822 / 3.399 / 3.444 / 3.426 | 3.419 | 3.444 | +0.025 (+0.7%) |
| Compact annotation, 1M alleles -> 4,224,817 rows (s) | 0.683 | 0.681 | -0.002 (-0.3%) |
| Complete-17 annotation, same rows (s) | 1.459 | 1.476 | +0.017 (+1.2%) |
| Peak RSS after load + annotation, incl. staged 1M-row tables (KiB) | 5,528,068 | 5,528,004 | -64 |

The budget build is not faster than the base. Load is 0.7% and complete-17 1.2% slower in these medians; both are under the 3% bar. Three further checks say this is code layout, not allocator cost: `perf record` shows no `duckvep_budget_*` function among the hot symbols and the same instruction profile in both builds; a whole GRCh38 load makes only 264 budgeted allocations; and when both trees are built with `-falign-loops=32 -falign-functions=32` (a scratch experiment, not in the tree) load medians are 3.285 (budget) and 3.284 s (base). In a seven-round alternating comparison on the 100,957-site bench corpus the complete-17 difference was +0.005 s and, in a later round, zero. The load difference of about 20-40 ms reproduced in every Release comparison.

## Native bytes by owner (MiB)

| Owner | after load | load high-water | annotation high-water (1M alleles, 1 thread) |
|---|---:|---:|---:|
| model | 1304.1 | 2680.1 | 1304.1 |
| index | 13.7 | 13.7 | 13.7 |
| reference | 0 | 0 | 0 (no FASTA) |
| workspace | 0 | 0 | 4.6 |
| scratch (worker lease) | 0 | 0 | 0.1 |
| emit (result/text arenas) | 0 | 0 | 20.3 |
| control | 0.0015 | 0.0015 | 0.0 |
| **total** | **1317.8** | **2680.1** | **1342.8** |

The load high-water is the reordering step that copies coding sequence and flank bytes into exact-capacity arrays while the doubled-capacity originals are still charged (realloc overlap counted). It fits the 4 GiB budget but not a 2.6 GiB one. The steady model is 1.30 GiB. The target of at most 2 GiB reserved for model plus index is met at rest, not at the load peak.

## Enforcement on GRCh38 (`duckvep_native_budget_set`)

| Budget | Result of `duckvep_model_load` | Bytes charged afterwards |
|---|---|---:|
| 1 GiB | capacity error while loading coding sequences (requested 256 MiB, 814 MiB in use) | 512 (registry) |
| 2 GiB | capacity error while ordering transcript sequences (requested 662 MiB, 2018 MiB in use) | 512 |
| 2.6 GiB | capacity error while ordering transcript sequences | 512 |
| 4 GiB (default) | loaded | 1317.8 MiB |

## Frozen receipts

`Rscript benchmarks/scale_contracts.R` with `DUCKVEP_SCALE_EXTENSION` set to the immutable budget binary printed `verified: compact, compact_regulation, complete17, 41942 retained cells; model 0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b`. The script fails on any row-count, hash or SHA-256 drift, and `git status` shows no receipt file changed.

## Fault injection

`make test_fault_injection` (AddressSanitizer + LeakSanitizer, `-DDUCKVEP_FAULT_INJECTION`). The fixture loads a 12-transcript model with peptide edits, a mature miRNA, interval features and a reference FASTA, then annotates small variants with HGVS, a deletion, a duplication and a breakend, and the projected builder, and phased calls through `duckvep_haplotypes` (a same-codon pair, a frame-opening insertion restored by a deletion, a stop gain: the frame classifier with per-carrier consequence and impact lists). A second model, on a circular region with a transcript, an exon and a regulatory feature that cross the origin, is executed on a lifted copy and annotated with HGVS, with regulation, and with the projected builder. Allocation counts: load 125, circular load 92, HGVS annotation 131, regulation annotation 106, projected annotation 60, haplotype classifier 29, lifted HGVS 84, lifted regulation 60, lifted projected 57. Each allocation is failed once, each in its own DuckDB process: 744 failures, 217 in load and 527 in annotation, 143 distinct call sites, among them every allocation in `duckvep_lift_open` (MODEL owner, built at load), `duckvep_lift_resolve` (WORKSPACE owner) and the lifted reference-window and resolver scratch of `duckvep_annotate.c` (worker scratch lease). Every one returned an explicit error, published no model, left the pinned model answering with the exact golden hash, returned the budget to its baseline and produced no sanitizer report; the harness also fails if any of those lifted sites was not reached.

The per-carrier consequence and impact lists of `duckvep_haplotypes` are written into DuckDB output vectors, not native memory owned by DuckVEP; the frame classifier allocates only from the already-charged haplotype workspace (29 budget allocations in the fixture query, all failed cleanly).

## Emitted-output allowance: 64 MiB to 256 MiB

Scale-runner smoke (2 jobs x 6 threads, ceilings enforced, 1M-allele genome panel from chr20/21/22/Y, 13,518,845 output rows) with the original 64 MiB per-worker emitted-output allowance: compact succeeded, complete-17 failed in both jobs with `capacity error: per-worker emitted-output lease budget exceeded (requested 27262976 bytes, 57683968 in use, limit 67108864)`. The error was explicit and nothing was published (receipt: `runs/smoke-2x1M-emit64/`). With 128 MiB complete-17 succeeded and the emit high-water was 118-128 MiB. The default is now 256 MiB, about twice the measured peak; the allowance is charged as it grows, so the higher cap reserves nothing. The 4 GiB budget, 6 workers, 128 MiB scratch and 64 MiB idle retention are unchanged.
