# Haplotype scale qualification (issue #2, slice 7)

The last slice of the signed `duckvep-coding-v1` contract (`design/duckvep_haplotype_contract.md`, section 4): DuckVEP must be at least 2x faster than
the pinned `bcftools csq` on identical input, transcript model and core, under the same caps, and must qualify a 5M-physical-variant single-sample
GRCh38 job on the MANE-selected model.

## 2026-09-30, issue #34 slice 1: the fused native reader (gate met)

`duckvep_coding_calls(model, path)` reads the (bgzipped) VCF or BCF with the bundled zlib-only htslib, maps CHROM by name, runs the shared discovery (`src/core/duckvep_core_discovery.c`, the same code as `duckvep_coding_transcripts`) on every ALT allele
and decodes FORMAT, GT and PS only for the records that touch coding sequence (25,969 of 4,023,088). It emits the calls relation of mode B. The measurement below is the unchanged gate: fresh process, one core (`taskset -c 6`, cgroup 16 GiB, native budget 4 GiB,
DuckDB threads 1), median of three, csq (`-p a -Ou -o /dev/null`) alternating with DuckVEP, load at most 3 at each start (waited once for one minute), sibling thread 0% busy in every run. Files: `perf34_slice1/` (receipt `process.tsv`, per-run stage output, `identity.txt`).
The harness is `benchmarks/haplotype_scale/cli_worker.sh` (DuckDB v1.5.1 CLI, no R); the R worker has the same pipeline as `--mode F`. Extension: release build of source tree `8713686b` (sha256 `ab4dc8e7...`), copied to an immutable file first.

| | wall s (runs) | median | gate |
|---|---|---|---|
| pinned `bcftools csq -p a -Ou -o /dev/null` | 18.07 / 17.25 / 17.08 | **17.25** | target 0.5 x = 8.63 |
| DuckVEP fused, staged (`--mode F`: calls table, then `duckvep_haplotypes` to Parquet), CLI | 8.23 / 8.21 / 8.15 | **8.21** (8.15 to 8.23) | **2.10x: met** |
| DuckVEP fused, inline (`--mode I`: `duckvep_haplotypes` over the reader, nothing staged), CLI | 8.24 / 8.20 / 8.16 | 8.20 | 2.10x |
| same pipeline from R (`worker.R --mode F`, R and DuckDB R package start-up included, reported separately) | 8.59 / 8.60 / 8.64 | 8.60 | 2.01x |

Host note: csq ran 1.6 to 2.6 s faster than in the 2026-09-29 round (19.66 s) on a quieter host, so the ratio is judged inside this round. Against the old csq median the ratio would be 2.39x. Mode B on this round's host was not re-run for the gate; its last measure is 12.58 s.

Stages, medians of the three staged CLI runs (the CLI's own timers; the wall clock adds about 0.1 s of start-up and shutdown):

| stage | mode B (2026-09-29, R) | fused |
|---|---|---|
| model load | 2.60 | 2.40 |
| decode, discovery, call construction (stage) | 6.38 | **3.08** |
| discovery, prediction and complete Parquet output | 2.88 | 2.57 |
| process wall | 12.58 | 8.21 |

What remains in the 3.08 s stage is the zlib inflate of 2.9 GB (a diagnostic `gzip -dc` is 6.8 s of single-threaded zlib-ng-less work, `bgzip -dc` with libdeflate 1.0 s; htslib's BGZF inflate with zlib blocks lies between them) plus about 4 million line splits and interval lookups. Load (2.4 s) is the next lever and is outside this slice. Sys time is 1.6 s of the 8.2 (page faults of the model).

Caps observed: memory.peak 3.83 GiB (4,107,853,824 bytes; cap 16 GiB), native pass total 1,651 MiB of 4,096 (model 1,304, index 14, workspace 333), DuckDB spill 0, all 15 processes exit 0.

Output identity. The output of the fused pipeline has 157,986 rows and full-output checksum `1456007180270799092358516`, equal to slice 7's mode B, in all nine audited-by-checksum outputs (F, I and R, three each; the checksum is `sum(hash(row))` over the Parquet, as in `worker.R`). (The brief quoted 158,004 rows; the slice 7 receipts and this run both have 157,986, and the checksum is the same.) Stronger: `benchmarks/haplotype_scale/verify_coding_calls.sh` builds mode B's calls relation and the reader's in one session and compares them as multisets, every column: 247,374 calls each, 0 only in the reader, 0 only in mode B, equal order-independent hashes.

Differences from the mode B SQL that cannot show in HG002's output, by design: the PS of a sample is read by FORMAT key (mode B takes the last `:` field, which is PS in HG002), the record ordinal counts every data record, a coding record without GT or a file without FORMAT/GT is an error (mode B would parse garbage), and multi-sample files give one row per sample with `sample_index` from 0 (mode B reads the first sample only).

Gates: `make release`, `test_release` (31 SQL files including `duckvep_coding_calls.test`), `test_release_asan`, `test_properties`, `test_haplotype_contract`, `test_fault_injection` (861 failed allocations, all clean, including new call sites in the reader and the discovery module), `check-function-docs` (26 functions, 29 examples), R tinytests (10,711 expectations, all pass), `scripts/check-rduckvep-bundle.sh`, `git diff --check`, clang `-Werror=string-concatenation` on the new files.

## Verdict

**(Superseded by the 2026-09-30 section above: the gate is met with the fused reader.) As measured on 2026-09-29, the gate was not met.** The maintainer closed #2 on its correctness contract, and the gate moved unchanged to #34. Judged on mode B (identical unsorted VCF, decode through Parquet), cold process, model load included, whole-process wall clock
(R and DuckDB start-up included), one core, three fresh processes each:

| | wall s (median of 3) |
|---|---|
| pinned `bcftools csq -p a -Ou -o /dev/null`, re-measured here | **19.66** (20.90 / 19.66 / 18.07) |
| DuckVEP mode B, cold, model load included | **12.58** (12.83 / 12.58 / 12.04) |
| gate: 0.5 x csq | 9.83 |

DuckVEP takes 0.640 of csq's time: 1.56x faster, 2.75 s (28% of the target) short of the gate. Other readings:

| reading | ratio to csq `-Ou` (19.66 s) |
|---|---|
| mode B cold, whole process (the gate) | **0.640** |
| mode B cold, without R/DuckDB start-up (load + stage + predict) | 0.608 |
| mode B cold with the VCF inflated by `bgzip -dc` instead of DuckDB's miniz (diagnostic, not the gate) | 0.561 (11.03 s) |
| mode B warm (model resident, second pass in the same process) | **0.429** (8.44 s): met |
| mode A cold / warm | 0.337 / 0.134 |
| mode B cold against csq writing BCF (`-Ob`, 32.48 s) | 0.387 |

This is the second measurement of the qualification. The first (before discovery was made exact, see "Pipeline") gave 10.87 s against 18.27 s (0.595); the host was quieter then, and the exact
discovery costs about 1.2 s of stage time in mode B (per-ALT event normalization, see below). The ratio is a comparison inside one round of alternating runs, but the round to round noise of this shared host is
visible in csq itself (18.1 s to 20.9 s). The contract's recorded baseline (22.2 s on 2026-09-28, target 11.1 s) was measured on core 2.

The floor is two costs that DuckVEP cannot remove from inside the extension: DuckDB inflates the 2.9 GB VCF text with its bundled miniz (the same job with `bgzip -dc` feeding it has a stage 1.6 s shorter)
and loading the 1.8 GB transcript model takes 2.5 s (csq parses its 108 MB GFF3 on every run). Load + stage + predict is 11.8 s of the 12.6 s process; the stage (6.4 s) is dominated by the decode.

## What was measured

Host: Intel i5-13500, kernel 6.8.0-78, other agents' jobs running elsewhere. Every timed process was a fresh process pinned with `taskset -c 6`
(SMT sibling 7) inside `systemd-run --scope -p MemoryMax=16G -p MemorySwapMax=0` (`capped_run.sh`), started only when the one-minute load average was at most 4 (it
polls every 60 s for up to 30 minutes; the receipt records the load at the start, the time waited and how busy the sibling thread was: loads
2.12 to 3.94, sibling thread at most 4% busy except in four runs (8, 9, 12 and 18%, none of them a gate run), five runs waited a minute or more). Rounds alternate csq and DuckVEP so both see the same conditions. 78 processes, all exit 0 (the receipts have 90 rows because a process reports its cold and warm pass).

Caps on every DuckVEP run: 16 GiB per process (cgroup), native budget 4 GiB (`duckvep_native_budget_set`), DuckDB `memory_limit` 8GB, temp directory capped at 32GiB, one
DuckDB thread. Observed maxima are under "Caps observed". csq ran under the same cgroup cap.

Pins (`inputs.tsv`, by content hash): input `hg002.ens.vcf.gz` sha256 `54e5d09d...` (4,023,088 records, 4,070,522 ALT alleles; it is the file that
`test/data/haplotype/hg002_strict.oracle.tsv` names), Ensembl 116 GRCh38 GFF3 `08e881d9...`, FASTA `1e74081a...`, model `homo_sapiens_116_GRCh38_final.duckdb`
(`model_sha256` `03019059...`, the receipt in `hg002_domain.accounting.tsv`), bcftools `1.23.1-70-g6dbd8fef` with htslib 1.22.1, R 4.6.0 with DuckDB 1.5.5. The
DuckVEP extension is the release build of the source tree `3d23d8a6` (`git rev-parse HEAD:src`, unchanged by the later documentation and bundle commits); the binary's sha256 is in `inputs.tsv`.

### The csq baseline

`csq.log` (2026-09-28) records only stderr. Its three warnings ("keeping the first 31", the N_REF_PAD note, the `--ncsq 345` note) and its 811,296 KiB peak RSS match
`bcftools csq -f FASTA -g GFF3 -p a` at the default `-n`, and its 20.09 s user time matches an uncompressed output: writing BCF (`-Ob`) costs another 12 s of
CPU on this host (user 29.4 s). The gate baseline is therefore `-p a -Ou -o /dev/null`; csq writing a BCF file is reported beside it. The oracle invocation (`-p s -n 1024`) is for
correctness, not timing. Median of three fresh runs each (table below).

### The pipeline

Both modes go through the public builder `duckvep_haplotypes` and end in a Parquet file (written to `<out>.partial` and renamed only after success, so a failed statement
publishes nothing at the final name).

- **Mode B**, from the VCF text: `read_csv` of the bgzipped VCF (DuckDB's own decoder), one row per record; the ALT alleles are unnested and, for each, `duckvep_coding_transcripts('hap', seq_region, pos, ref, alt)`
  returns the transcripts whose coding sequence the normalized event touches (see Discovery); 25,969 of 4,023,088 records (0.65%) reach a call. One call per ALT allele and transcript is built with genotype (`GT` to
  alleles and per-lane phase flags; an integer `PS` is the phase set, the string labels `PATMAT`, `HOMVAR` and `.` are none); then `duckvep_haplotypes` sorts, replays and writes
  the complete output (all 31 columns) to Parquet.
- **Mode A**, from the preordered calls relation (the calls mode B builds, sorted by region, position, event and transcript, read back from Parquet): stage the calls, `duckvep_haplotypes` (its own
  internal sort included), Parquet.
- **Cold** is a fresh process: R start-up, model load, first execution. **Warm** is the second execution in the same process. Cold-only runs are the gate runs (their wall clock is the cold wall
  clock, and they skip the untimed audit); separate runs add the warm pass and the audit (counts and a full-output checksum, `sum(hash(row))` over the Parquet).
- **Model load** runs `duckvep_model_load` on the model file's own tables without `ORDER BY` (the file is stored in load order and the loader rejects rows that arrive out of order; see "Pipeline choices").
  The MANE Select model is the same builder with the transcripts of `mane_select_refseq IS NOT NULL` (19,297) re-indexed densely in (region, start, index) order.
- **Capacities**: `max_alignment_cells` 268,435,456 and `workspace_limit` 1 GiB (the dense control uses 1,073,741,824 and 2 GiB); with the builder's own defaults the HG002 job also succeeds (`failures.tsv`, `builder_defaults`).
- **Discovery** (`duckvep_coding_transcripts(model, seq_region, position, reference, alternate)`) returns exactly the (event, transcript) pairs that `duckvep_annotate_sql` finds with the CDS bit of `region_mask`, no more and no fewer.
  It takes the VCF record as written and trims it inside the scalar, with the annotation builder's own normalization (shared prefix and suffix removed, so a shared anchor base never creates an overlap by itself; an insertion is an
  interbase point; introns of at most 13 bases inside the CDS count as coding, as in VEP). `verify_discovery.R` on the whole HG002 input: the builder finds 247,374 pairs, the scalar 247,374, 0 missing and 0 extra. Of the 23,866 records the committed csq accounting
  compares, 23,865 are discovered; record 2402008 (`TCTCC>T`, chr10:15719771) touches the CDS only through its anchor base, and the annotation builder does not count it either. Mode B with the scalar, mode B with the
  annotation builder's discovery and mode A write the same output (checksum `1456007180270799092358516`).

Everything is in `benchmarks/haplotype_scale/`: `run_qualification.sh` (the one command), `prepare_inputs.R` (untimed inputs), `worker.R` (one fresh process), `capped_run.sh`, `failures.R`, `summarize.R`,
`tables.R`, `verify_discovery.R`.

## Results

### csq baseline and the HG002 gate

| invocation | wall s (runs) | median wall s | median user s | peak RSS KiB | memory.peak MiB | load at start |
|---|---|---|---|---|---|---|
| `-p a -Ob -o file.bcf` | 32.48 / 30.35 / 33.63 | 32.48 | 30.85 | 811,312 | 859 | 3.4 / 3.8 / 3.9 |
| `-p a -Ou -o /dev/null` | 20.90 / 19.66 / 18.07 | 19.66 | 18.51 | 811,312 | 792 | 3.8 / 3.7 / 3.9 |

| mode | cold process wall s (runs) | cold wall median | load | stage | predict + Parquet | warm pass median | warm stage | warm predict |
|---|---|---|---|---|---|---|---|---|
| B | 12.58 / 12.83 / 12.04 | 12.58 | 2.60 | 6.38 | 2.88 | 8.44 | 5.75 | 2.69 |
| A | 6.63 / 6.72 / 6.63 | 6.63 | 2.90 | 0.04 | 3.04 | 2.63 | 0.03 | 2.60 |

Load, stage and predict are the worker's internal timers (medians); the wall clock adds R and DuckDB start-up and shutdown. In mode B the stage is decode, discovery and call construction (of which the
DuckDB decode of the 2.9 GB text is the larger part); in mode A it is reading the ordered calls. Mode A and mode B write the identical output: 157,986 rows, checksum `1456007180270799092358516` in every audited run
(the same checksum across the four optimization variants, including the unmodified base, and across both discovery routes). 149,896 rows are `predicted`.

Accounting of the HG002 input: 4,023,088 records = 25,969 with a call in coding scope (247,374 calls, emitted as 157,986 output rows, of which 149,896 are `predicted` and the rest carry their status and reason) + 3,997,119 outside coding
scope (no CDS-exon overlap in the model, so no call) + 0 unavailable. The 5M job: 5,000,000 = 294,000 in coding scope (36,442 rows) + 4,706,000 outside.

### Inside and outside the csq domain

The committed accounting's partition (regenerated here and equal to the committed receipt: 23,866 compared, 3,998,802 outside, ledger hash `81c9b983...`) applied to the calls, mode A, cold:

| partition | records with calls | calls | output rows | predicted | predict + Parquet s (median) |
|---|---|---|---|---|---|
| inside the csq domain (23,866 compared records, 23,865 with a call) | 23,865 | 239,300 | 153,654 | 146,987 | 2.69 |
| outside (records that overlap a CDS but are not compared) | 2,104 | 8,074 | 6,461 | 5,043 | 0.10 |

The other 3,997,119 records of the 4,023,088 reach no call (no coding overlap): their whole cost is the shared stage, decode plus discovery, 6.4 s in mode B. So prediction inside
the csq domain is 2.69 s, outside 0.10 s, and the stage (dominated by decoding the records outside) is 6.4 s.

### Caps observed (maximum over the three cold processes of each job)

| job | memory.peak GiB (16 cap) | max RSS GiB | native high-water MiB, load (total) | native high-water MiB per owner, pass | native pass total MiB (4096 cap) | DuckDB spill peak bytes |
|---|---|---|---|---|---|---|
| hg002_B | 4.14 | 4.16 | 1,582 | model 1,304, index 14, workspace 333, scratch 0, emit 0, control 0 | 1,651 | 0 |
| hg002_A | 3.89 | 3.92 | 1,582 | model 1,304, index 14, workspace 333, scratch 0, emit 0, control 0 | 1,651 | 0 |
| qual5m_B | 3.14 | 3.11 | 97 | model 74, index 0, workspace 333, scratch 0, emit 0, control 0 | 407 | 0 |
| qual5m_A | 2.55 | 2.55 | 97 | model 74, index 0, workspace 333, scratch 0, emit 0, control 0 | 407 | 0 |
| dense_B | 2.34 | 2.38 | 97 | model 74, index 0, workspace 1,101, scratch 0, emit 0, control 0 | 1,175 | 0 |
| low_sharing_B | 1.99 | 2.02 | 97 | model 74, index 0, workspace 333, scratch 0, emit 0, control 0 | 407 | 0 |

The native model high-water at load fell from 2,680 MiB (`benchmarks/data/scale_contracts/budget-evidence.md`) to 1,582 MiB (the reorder copy is skipped when rows arrive in order). DuckDB never spilled: the staged relations are small (the largest is the 294,251-row 5M calls).

### 5M-physical-variant job on the MANE-selected model

Input: HG002's 4,023,088 records plus 976,912 records sampled by rule from the staged gnomAD v4.1 exomes = 5,000,000 physical records (5,047,434 ALT alleles), single sample, all phased.
Top-up rule (`prepare_inputs.R`, deterministic): candidate sites are the gnomAD exome ALT rows with status `literal` and FILTER `PASS` on chromosomes 1-22, X and Y whose (contig, position) is absent
from HG002; one ALT per site is kept by the smallest `hash(20260929, contig, position, ref, alt)`; the first 976,912 sites by that hash are taken; the genotype is the same hash again: `0|1` 45%, `1|0` 45%,
`1|1` 10%, with `PS` `PATMAT` for heterozygotes and `HOMVAR` for homozygotes (HG002's own string labels, so all phased records of the sample form one phase domain per transcript and
compose with HG002's). The exomes are coding-dense, which is why it is a harder job than HG002 alone. The merged file is sorted by contig and position and bgzipped (`qual5m.vcf.gz`, sha256 in `inputs.tsv`).

| mode | cold wall s (runs) | median | load | stage | predict | warm pass | output checksum |
|---|---|---|---|---|---|---|---|
| B | 16.39 / 16.65 / 16.37 | 16.39 | 1.02 | 6.08 | 8.79 | 14.80 | 335591545666261315447574 |
| A | 10.47 / 11.04 / 10.39 | 10.47 | 0.99 | 0.04 | 8.95 | 8.86 | 335591545666261315447574 |

What the job contained (audit of the warm runs, identical for modes A and B):

- physical records 5,000,000, ALT alleles 5,047,434; physical sources in coding scope 294,000 (5.9%), projected ALT events 294,017, projections (event x transcript) 294,251, calls 294,251; every one of the other
  4.7M physical records is explicitly outside coding scope (no CDS-exon overlap in the MANE model);
- carrier states 36,514, unique paths (output rows) 36,442, translated bases 63,442,127, output 70,515,196 bytes; 35,476 rows `predicted`, 966 rows not (`outside_cds` 308, `overlapping_edits` 285,
  `unphased_heterozygous` 157, `curated_transcript` 116, `allele_over_50_bases` 61, `contradictory_edits` 37, `incomplete_cds` 2): nothing is dropped silently;
- peak active window: 18,700 transcripts carry calls, at most 916 events on one transcript (the peak of events held at once; capacities are 16,384 events and 4,096 transcripts by default); native workspace high-water 333 MiB;
- caps: memory.peak 3.14 GiB, native 407 MiB, DuckDB spill 0 bytes (see the table above);
- the output checksum is the same in every audited run.

Controls (mode B, three fresh processes; `dense` on the longest MANE coding sequence, `low_sharing` one variant per MANE transcript):

| control | median wall s | predict s | calls | output rows | carriers | translated bases | max events per transcript | transcripts | output bytes |
|---|---|---|---|---|---|---|---|---|---|
| dense_B | 3.36 | 1.86 | 4,000 | 2 | 2 | 215,914 | 4,000 | 1 | 552,711 |
| low_sharing_B | 2.06 | 0.45 | 19,116 | 19,099 | 19,099 | 33,382,784 | 2 | 19,087 | 27,861,760 |

The dense control is 4,000 non-overlapping variants (75 indels) on one haplotype of the 107,976-base coding sequence of ENST00000589042 (TTN): 4,000 events on a single transcript, two paths. It needs `max_alignment_cells`
of about 525 million (the exact band for thousands of edits), hence its larger capacity. The low-sharing control is 19,065 records on 19,087 transcripts: 19,099 paths, no sharing.

### Explicit failures (`failures.tsv`, one process, each control followed by a check that the connection and the model still work)

| control | result |
|---|---|
| native budget 64 MiB before `duckvep_model_load` | `capacity error: native memory (model) budget exceeded (requested 8388608 bytes, 65164640 in use, limit 67108864)`; 400 bytes charged afterwards (the registry); no model published; after restoring the budget the load succeeds |
| budget = resident model + 64 MiB before `duckvep_haplotypes` | `capacity error: native memory (workspace) budget exceeded (requested 16781312 bytes, ...)`; no table, no output file at the final name (`COPY` leaves a 2,581-byte `.partial` file, removed); native bytes back to the resident value; the next query runs |
| `max_alignment_cells := 100000` on the HG002 calls | `CDS difference status 3, max_alignment_cells=100000, required=166049 at transcript 669`; nothing published; bytes unchanged; reusable |
| the builder's default capacities on the HG002 calls | succeeds (157,986 rows) since capacity follows the band an alignment needs |

Before the capacity change below, the unmodified build could not run the 5M MANE job at any admissible cap: `max_alignment_cells=268435456` failed with `required=451073810 at transcript 3762` and
`3221225472` with `required=11659140506 at transcript 6435` (`base_capacity_refusal.tsv`); the final build runs it at 268,435,456.

## Optimizations, each with its A/B effect

Every change leaves output byte-identical (the same full-output checksum in every variant) and the SQL, contract and goldens suites unchanged. Profiling (`perf`) of the 22 s prediction of the first
measurement showed 67% in `duckvep_sequence_differences`; the model load and the coding transcripts lookup came next.

| variant (commit) | model load s | predict + Parquet s | process wall s (mode A cold) |
|---|---|---|---|
| base (`01f37f4`) | 3.43 | 22.21 | 26.76 |
| + bound alignment work by the optimum | 3.45 | 3.40 | 7.95 |
| + table validation and no reorder copy in model load | 2.77 | 3.40 | 7.28 |
| + table fast paths for translation, first stop and normalization | 2.51 | 2.80 | 6.44 |
| final | 2.50 | 2.72 | 6.35 |

1. **Alignment work** (`duckvep_sequence_differences`): the CDS and protein difference alignment filled the band of a feasible positional path, the whole matrix after any frameshift. The band now starts at the length change and
   widens until every path leaving it provably costs more (so the traceback, ties included, is that of the full matrix); rows inside the common prefix are closed-form and not computed; equal-length comparisons
   scan for differing runs; the counting pass is skipped when capacity cannot overflow; validation is word-wise. 22.21 s to 3.40 s. Checked against the previous implementation on 8 million random and repeat-rich cases
   (also under ASan/UBSan) and, committed, against the independent full-matrix oracle on 40,000 long repeat-rich pairs.
2. **Model load**: model open validated each flank and CDS byte with a branchy test; a table over eight-byte groups does the same test. Rows that arrive in transcript order no longer have their sequence pools copied a
   second time (load high-water 2,680 to 1,582 MiB). 3.45 s to 2.77 s.
3. **Codon fast paths**: three unambiguous bases index the amino-acid table directly, in translation and in the first-stop scan, and whole-sequence normalization is a table; ambiguous and invalid codons take the unchanged
   exact path. Model load 2.77 to 2.51 s (the first-stop cache), prediction 3.40 to 2.80 s.
4. **Capacity follows the band an alignment needs** (`7e49a27`): the trace capacity was compared with the cells of the feasible bound before any work. A deletion and a substitution 8,000 bases apart in a 20 kb sequence needed 322
   million cells for a result within a few columns of the diagonal; TTN needed 11.7 billion, more than the native budget, so the 5M MANE job could not run at all. Capacity is now checked against the band of each attempt; a band still to be tried
   that does not fit is refused explicitly with its cell count (the SQL tests' pinned limits are unchanged). Successful output is identical.
5. **`duckvep_coding_transcripts`**: transcript discovery from the resident interval index, with the annotation builder's event normalization, exact pair for pair. It is public (SQL docs in the README, `rduckvep_coding_transcripts()` in R) rather than
   a runner helper; see Discovery. It is not folded into `duckvep_haplotypes` because that builder consumes a relation of calls the caller stages (the calls carry genotype and phase that only the caller can build); the scalar is the
   fast way to produce the (event, transcript) rows of that relation.

Two pipeline choices, alternating, cold, median of 3 (`abx_*` rows in `runs.csv`):

| choice (process wall includes the untimed audit) | process wall s | load s | stage s | max RSS GiB |
|---|---|---|---|---|
| discovery by the annotation builder (twelve transcripts per intronic event) | 17.18 | 2.48 | 11.04 | 6.23 |
| discovery by `duckvep_coding_transcripts` | 11.95 | 2.49 | 5.85 | 4.28 |
| model load with `ORDER BY` on the model queries | 3.65 | 3.10 | | 5.06 |
| model load without (stored order; the loader rejects out-of-order rows) | 3.01 | 2.50 | | 3.92 |

Both discovery routes give the same 247,374 calls, 157,986 output rows and checksum. Making the scalar exact (normalizing each ALT allele and applying the frameshift-intron rule) raised its stage from about 5.1 s to about 5.9 s in this comparison;
performance work on it is left to the maintainer's decision on the gate.

## What would close the gap

The 2.75 s is mostly not in DuckVEP's own path: the miniz inflate of the VCF (a libdeflate decoder saves about 1.5 s and gets the ratio to 0.561), model load 2.5 s, R start-up and shutdown about 0.8 s, and prediction and
Parquet 2.9 s. Reaching 0.5 needs another 1.2 s beyond a libdeflate decoder; the candidates are a libdeflate-backed decoder (outside this extension), a cheaper exact discovery (per-ALT normalization is about 0.8 s of the stage),
streaming the model-load queries instead of materializing 1.1 GB of results, and avoiding the flank sequence load for coding-only jobs (the loader currently requires flanks whenever a CDS is present). No performance change was made after the
discovery fix.

## Files

`runs.csv` (one row per process and pass, with process wall, user, system, RSS, cgroup peak, load, sibling busy, native budget per owner, audit), `summary.csv` (medians and ratios), `inputs.tsv`,
`input_counts.tsv`, `failures.tsv`, `base_capacity_refusal.tsv`. Regenerate: `Rscript benchmarks/haplotype_scale/prepare_inputs.R WORK && WORK=... bash benchmarks/haplotype_scale/run_qualification.sh &&
Rscript benchmarks/haplotype_scale/summarize.R WORK/results && Rscript benchmarks/haplotype_scale/tables.R`. Large inputs and outputs are not committed.
