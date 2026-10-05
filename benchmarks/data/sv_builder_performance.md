# STR, SV and BND builder time and memory

Command (one R process per builder and corpus, one DuckDB thread, immutable extension copy outside `build/release`):

```sh
Rscript benchmarks/benchmark_sv_builders.R /path/to/copy/duckvep.duckdb_extension \
  benchmarks/data/duckvep_sv_builders_100k.csv 100000
Rscript benchmarks/benchmark_sv_builders.R /path/to/copy/duckvep.duckdb_extension \
  benchmarks/data/duckvep_sv_builders_1m.csv 1000000
```

Extension SHA-256 `bbc690cc3fac46f1651a8c6854030fe8bf9e938572eb6288a153dca6970cedf3` (source `7a741c0` with the INT64 overflow guard and hash-join mate lookup); Linux x86_64,
Intel i5-13500, R `duckdb` package running DuckDB v1.5.5 with the v1.2.0 stable C API extension; raw rows are in
[duckvep_sv_builders_100k.csv](duckvep_sv_builders_100k.csv) and [duckvep_sv_builders_1m.csv](duckvep_sv_builders_1m.csv).

The timed step is `CREATE TEMP TABLE out AS SELECT * FROM query(<builder SQL>)` over an already materialized input relation, median of three runs
(min-max in the CSV). Building the SQL string takes about 1 ms for every builder. **Peak RSS is the process high-water mark** (`VmHWM`) after the three runs, including R, DuckDB
and the loaded extension (about 105-115 MB idle); `RSS before` is the high-water mark after the input tables were generated and before the builder ran, so
`peak - before` bounds the builder's own working set from above only when the input generation was smaller than the builder, which is not the case for
`breakend_fusion` (its input includes a run of the pairs builder, so its `before` is already the pairs peak).

## Real call set: Sniffles2 1KGP joint calls, GRCh38 chr22 (all 10,176 records)

Source `https://1kgp-sv-imputation.s3.eu-west-1.amazonaws.com/sv_calls/sniffles2_joint_sv_calls.vcf.gz` (Sniffles2 2.0.7; the chr22 sites extract is
recorded in the local corpus provenance): 3,452 DEL, 3,471 INS, 3,115 BND, 114 INV, 24 DUP.

| Builder | Rows | Median s | Rows/s | RSS before MB | Peak RSS MB | Statuses |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `duckvep_prepare_sv_geometry_sql` | 10,176 | 0.125 | 81,000 | 156 | 163 | ok 211; unsupported_geometry 9,965 |
| `duckvep_prepare_breakend_pairs_sql` | 10,176 | 0.055 | 185,000 | 156 | 178 | not_applicable 7,061; unproven 3,115 |

Sniffles2 writes sequence-resolved DEL/INS with `REF`/`ALT` of `N`, which the geometry builder does not treat as a symbolic or anchored-literal
allele, so 9,965 records are `unsupported_geometry` (retained, not repaired), and its BND records are single breakends without `MATEID`
(`unproven`). The timing is real; the acceptance rate is a fact about this producer. ExpansionHunter, fusion and structural-HGVS inputs
have no equally large public call set with the required literal reference, so they are measured on synthetic corpora.

## Synthetic corpora, 100,000 and 1,000,000 events

Deterministic hash-generated events: geometry mixes DEL/DUP/INV/INS/TDUP and literal insertions with CIPOS/CIEND and `SEQ`; pairs are reciprocal
BND pairs in four bracket orientations with 2% event conflicts and 1% malformed ALT; fusion joins the pairs-builder output with 0-2 genes per record;
structural HGVS uses 1-40 bp DEL/DUP/TDUP/INV/CNV with random reference sequence; ExpansionHunter mixes `REPCN` and `CN` formats, exact and summary
alleles, and 10% reference mismatches.

| Builder | Rows | Median s | Rows/s | RSS before MB | Peak RSS MB |
| --- | ---: | ---: | ---: | ---: | ---: |
| `duckvep_prepare_sv_geometry_sql` | 100,000 | 0.139 | 719,000 | 116 | 132 |
| `duckvep_prepare_breakend_pairs_sql` | 100,000 | 0.207 | 483,000 | 123 | 234 |
| `duckvep_prepare_breakend_fusion_sql` | 100,000 | 0.126 | 794,000 | 226 | 231 |
| `duckvep_prepare_structural_hgvs_sql` | 100,000 | 0.345 | 290,000 | 126 | 170 |
| `duckvep_prepare_expansionhunter_sql` | 100,000 | 0.442 | 226,000 | 127 | 179 |
| `duckvep_prepare_sv_geometry_sql` | 1,000,000 | 1.343 | 745,000 | 176 | 299 |
| `duckvep_prepare_breakend_pairs_sql` | 1,000,000 | 2.107 | 475,000 | 248 | 1,240 |
| `duckvep_prepare_breakend_fusion_sql` | 1,000,000 | 1.360 | 735,000 | 1,200 | 1,212 |
| `duckvep_prepare_structural_hgvs_sql` | 1,000,000 | 3.450 | 290,000 | 244 | 646 |
| `duckvep_prepare_expansionhunter_sql` | 1,000,000 | 4.358 | 229,000 | 297 | 665 |

Ten times the rows takes 9.7-10.8 times as long for every builder: all five are linear at these sizes. The heaviest working set is the pairs
builder at one million records (about 1 GB over the idle process).

## Performance comparisons and guard behavior

- **BND mate lookup.** The nested-loop plan joined on `m.event_index = mi.first_index AND mi.n = 1` and took 20.2 s for 100,000 records (16 s for 50,000 pairs in the profile). The hash-join formulation folds uniqueness into the join key, `m.event_index = CASE WHEN mi.n = 1 THEN mi.first_index END`; it processes 100,000 records in 0.207 s, about 100 times faster, with byte-identical results on all tests. The 40,000-record SQL test in `test/sql/duckvep_structural.test` pins correctness at this scale. The 0.055 s real-call-set result covers 10,176 unpaired records and is not a comparable pair-building workload.
- **INT64 geometry overflow.** INT64 arithmetic for `POS`, `END` and confidence-interval bounds caused a row with `POS = 9223372036854775807` to abort the statement. The geometry builder uses HUGEINT sums with `TRY_CAST`; this row returns `unsupported_geometry`. `test/sql/duckvep_builder_errors.test` covers this with the other error-recovery cases.
