# DuckVEP and vep-rs

[vep-rs](https://github.com/natera-open-source/vep-rs) is a Rust reimplementation of Ensembl VEP that reports the
highest concordance and speed of the VEP ports we know of. This page puts DuckVEP next to it on one input, on one
machine, and lets executable Ensembl VEP decide every disagreement.

[`benchmark_duckvep_vep_rs.sh`](benchmark_duckvep_vep_rs.sh) produces paired timing and adjudication receipts.
The measurements below are from 2026-10-04 at DuckVEP `125e34a`. The 33.53 s and 102.86 s figures were
measured in separate invocations and do not establish a controlled comparison.

## Setup

| | |
|---|---|
| Input | GIAB HG002 v4.2.1 GRCh38, split to one ALT allele per record, sites only: 4,070,522 alleles |
| Annotation | Ensembl 116 GRCh38. DuckVEP: its model built from the Ensembl core database. vep-rs 0.3.1: a JSON cache built by its own `vep-cache-builder` from `Homo_sapiens.GRCh38.116.gtf.gz` |
| Output | DuckVEP: the rich relation (consequences, impact, cDNA, CDS and protein positions, amino acids, NMD) as Parquet. vep-rs: VEP tab text |
| Machine | Intel Core i5-13500 (14 cores, 20 threads), 62 GiB; outputs written to memory-backed storage |
| Oracle | Ensembl VEP 116, `--gtf` on the same GTF, for the sites where the two tools disagree |

DuckVEP's timed run is one DuckDB process: extension load, model snapshot restore, VCF read, coordinate sort,
annotation and the Parquet write. vep-rs is its release binary with `--fork`. A `target-cpu=native` build of
vep-rs used the same CPU time as the release binary (which already requires AVX2), so the release binary is
reported. vep-rs has no hand-written SIMD; its AVX2 use comes from the compiler flag and its dependencies.

`VEP_RS_CACHE` must name the species/version directory containing `transcripts/`, not its parent.
vep-rs can exit successfully with zero consequences for an invalid cache root; the runner rejects that result.

For a paired one-thread, one-CPU campaign, run the script with `--threads 1 --runs 3 --cpu-affinity 2`,
replacing `2` with one logical CPU allowed by the host. Each repetition runs vep-rs and DuckVEP sequentially
with the same affinity and thread count. `--threads 16` additionally runs that setting before the one-thread comparison; the same affinity
is applied to both tools and the stage probe. Without `--cpu-affinity`, neither tool is pinned. `--fork 1` and
`SET threads=1` limit worker counts but do not bind a process to one CPU.

The work directory contains `timings.tsv` for end-to-end process measurements, `stage_timings.tsv` for DuckVEP
query-only stages, and `run_receipt.tsv` for the campaign settings. Stage capture runs after preparation and
model restore and its rows are not included in `timings.tsv`. Use `DUCKDB_BIN`, `TIME_BIN`, `TASKSET_BIN`,
`VEP_RS_BIN` and `VEP_BIN` to select executables on another installation; `TIME_BIN` must support GNU time's
`-f`, `-a` and `-o` options. Preparation is reused only when its source paths and SHA-256 content hashes match the saved receipt
and all generated files pass their recorded SHA-256 checks.

## Time and memory

| Threads | Tool | Wall (s), median of 3 | CPU (s) | Peak RSS (GiB) | Output |
|---:|---|---:|---:|---:|---|
| 16 | DuckVEP | 6.84 | 71 | 4.9 | 47,278,065 rows, 142 MB Parquet |
| 16 | vep-rs | 9.39 | 105 | 3.5 | 36,258,238 rows, 4.17 GB text |
| 1 | DuckVEP | 33.5 | 33.5 | 2.6 | |
| 1 | vep-rs, `--fork 1` | 30.0 | 53.3 | 3.1 | |
| 1 | vep-rs, `--fork 1` pinned to one core | 34.4 | 34.4 | | |

At 16 threads DuckVEP finishes 1.37× sooner on two thirds of the CPU time, while annotating 30% more rows (see
below). `--fork 1` is not one core: vep-rs used 53 CPU-seconds in 30 s. The pinned vep-rs result and the
single-threaded DuckVEP result above came from separate timing invocations, so that pair does not establish a
controlled one-core comparison. DuckVEP's resident memory includes the 1.3 GiB model, which is a file mapping
shared between processes.

The original single-threaded DuckVEP probe separates query stages from setup and process startup; it was not
CPU-affinity pinned:

| Stage | Historical time (s) |
|---|---:|
| Annotation, compact output, rows counted | 3.4 |
| Annotation, rich output, rows counted | 8.0 |
| Rich output, every column read | 11.4 |
| The benchmark query: rich output, two identifier joins, Parquet write | 30.3 |
| Compact output written to Parquet | 15.5 |

Current probe timings are written separately in `stage_timings.tsv` with stage names `compact_count`, `rich_count`,
`rich_touch`, `rich_joined_parquet` and `compact_parquet`. The stage timer excludes the setup and model restore.

Writing 47 million rows to Parquet and attaching the VCF and transcript identifiers take about 19 s of the 30 s;
both parallelize, which is why DuckVEP pulls ahead at 16 threads.

The two outputs are not the same artifact. vep-rs writes VEP's own text format, which is what a drop-in
replacement must do; DuckVEP writes typed columns. vep-rs 0.3.1 did not accept `--parquet`, so a like-for-like
format was not measured.

### Paired single-core receipt

The [2026-10-05 paired run](data/vep_rs/bbe2ec4_paired/README.md) uses three cold-process pairs on
CPU 19 (an E-core), disk-backed outputs and a shared host. Median wall times are 58.09 s for DuckVEP
and 73.26 s for vep-rs. It retains five separate query probes and reproduces the agreement counts below.
Those hardware and storage conditions differ from the historical measurements above.

## Agreement

A tuple is (allele, transcript, set of consequence terms), the comparison vep-rs itself uses.

| | DuckVEP | vep-rs |
|---|---:|---:|
| Transcripts appearing in the output | 598,993 | 421,878 |
| Tuples on the 420,436 transcripts both emit | 34,148,222 | 34,148,222 |
| Matching the other tool | 34,146,531 | 34,146,531 |
| Different | 1,691 | 1,691 |

**Transcript coverage.** The Ensembl 116 GTF holds 646,577 transcripts; `vep-cache-builder` kept 447,179. The
transcripts DuckVEP annotates and vep-rs does not are mostly lncRNA (171,530 of them in this run). We did not
investigate why the cache builder drops them, and vep-rs's published concordance is on release 115, not 116.

**The 1,691 disagreements, adjudicated.** They sit on 1,291 sites. Executable VEP 116 on those sites:

| VEP's answer | Tuples |
|---|---:|
| agrees with DuckVEP | 1,689 |
| agrees with vep-rs | 1 |
| differs from both | 1 |

On those same 1,291 sites VEP emits 41,249 tuples in all. DuckVEP matches 41,085 and emits nothing VEP does not;
162 VEP rows are for transcripts in the GTF that DuckVEP's database-built model does not hold. vep-rs matches
38,109 and lacks 1,443 VEP rows.

The dominant disagreement is the reading frame: 1,317 tuples where DuckVEP says synonymous and vep-rs says
missense, and others where a stop is gained or retained on one side only. The two remaining mismatches are retained as [two-record witnesses](../test/duckvep/conformance/veprs_residual_witness.vcf):
- `e284879 / ENST00000696609`: a UTR/start deletion leaves an upstream `ATG` before the original CDS suffix.
  DuckVEP emits both start terms; VEP and vep-rs emit only `start_lost`.
- `e391762 / ENST00001011173`: deleting the middle base of a terminal `TGA` leaves local `TA`;
  sequence replay uses the next base to form `TAG`. DuckVEP calls the stop retained, VEP lost,
  and vep-rs also adds a frameshift term.

The [diagnostic SQL](../test/duckvep/conformance/veprs_residual_witness.sql) reproduces DuckVEP's rows from
this runner's snapshot and mappings. Load the extension and set the SQL variable `comparison_dir` before
running it from the repository root. Direct VEP 116 reruns reproduce both oracle rows. The differing
sequence contexts are identified, but VEP's internal shifted allele operands remain unresolved;
these observations do not justify changing a classifier.

## What this does and does not show

- It is one genome. vep-rs publishes concordance on 1.18 billion variants; DuckVEP's evidence is different in kind:
  smaller corpora plus generated rare states (boundaries, frames, long alleles) that population callsets seldom
  contain. This run is a reminder of why: 99.995% of tuples agree, and nearly all of the rest are states where
  vep-rs departs from VEP.
- The adjudication covers only sites where the tools disagree. Tuples on which both agree were not checked
  against VEP here.
- HGVS, regulatory features and plugins are outside this comparison.
