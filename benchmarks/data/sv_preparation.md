# R STR and structural preparation costs

Run `Rscript benchmarks/scripts/benchmark_sv_preparation.R` at this revision. Workload: 10,000 calls per path to observe R's peak Ncells/Vcells after warmup, then 2,000 requested iterations per path with `bench::mark` (GC excluded from per-call timing). Runtime: Linux x86_64 R 4.6.0, `bench`, one process, pinned ExpansionHunter v5 documented ALS allele and chr21 VEP 116 insertion witness. The timed preparation is R-only; model loading, DuckDB annotation, FASTA extraction and VEP Docker are not timed.

| Path | Median | Timed iterations | GC | Peak Ncells MB | Peak Vcells MB |
|---|---:|---:|---:|---:|---:|
| STR exact `<STR2>` | 217 µs | 1,973 | 27 | 38.4 | 6.6 |
| STR summary `<STR349>` | 157 µs | 1,991 | 9 | 38.4 | 6.6 |
| Symbolic INS + CIPOS/CIEND/SEQ | 263 µs | 1,981 | 19 | 38.4 | 6.8 |
| Literal G>GATG | 242 µs | 1,982 | 18 | 38.4 | 6.9 |

The peak numbers include the R session and package state, **not** per-event allocation or resident-set size; identical Ncells peaks reflect R's GC trigger. They bound this local short run, not production batch size. The STR preparer bounds exact components at the native repeat fact's default 5,000 bases. No C allocation path is added. The executable VEP differentials in `test/duckvep/conformance` establish correctness separately from these preparation measurements.

The native builders' time and peak resident memory on the public Sniffles2 chr22 call set and on 100,000- and 1,000,000-event corpora are in [sv_builder_performance.md](sv_builder_performance.md); the species structural differentials are in [sv_species_evidence.md](sv_species_evidence.md).
