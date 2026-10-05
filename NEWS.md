# DuckVEP

## Strict handling of a single unphased heterozygote (issue #12)

- Under `phase_policy := 'strict'`, an unphased heterozygous call that is the sample's only heterozygous or missing call on a transcript is read in slot order: one haplotype carries the allele and the other does not. Two or more such sites, or an unphased site beside a phased one, produce `incomplete_input`. In the HG002 measurement, 782 of 1,110 such paths were resolved and 157,343 of 158,137 paths were predicted.

## Lost-stop read-through and NMD exceptions (issue #13)

- When a stop is lost, `protein` continues through the stored 3′ transcript flank to the next stop or the end of that flank; `sequence_flags` bit 16 marks read-through. In the HG002 measurement, all 363 single-edit paths with a numbered new stop in per-variant protein HGVS had the corresponding protein length.
- `nmd_exceptions` reports `start_proximal` when a premature stop lies within the first 100 coding bases and `long_exon` when its exon exceeds 407 bases. The `ejc50` prediction is independent of these exception labels. In HG002, 234 of the 795 stop-gained paths predicted to `trigger` carried an exception.

## Haplotype annotation of incomplete CDS models (issue #12)

- Whole-haplotype prediction supports `cds_start_NF` and `cds_end_NF` transcripts, including CDS starts inside a codon. In the recorded HG002 comparison, predicted paths numbered 149,897 in the baseline run and 156,413 in the expanded-domain run, among 157,986 paths; paths decided by the baseline retained their results.
- Without an annotated CDS end, an edited stop is `stop_gained`, a frame still displaced at the annotation end is a frameshift, and an edit confined to the trailing partial codon is `incomplete_terminal_codon_variant`. Without an annotated start, no start-codon test is applied.
- Hand-derived cases and an HG002 comparison of single-edit haplotypes with per-variant annotation found 91,212 identical results among 91,708 paths; the others fall into named policy differences.

## vep-rs comparison (version 0.3.1)

- On GIAB HG002 (4,070,522 alleles, Ensembl 116), the recorded 16-thread end-to-end run measured 6.8 s for DuckVEP and 9.4 s for vep-rs. Separate one-core runs measured 33.5 s and 34.4 s; vep-rs `--fork 1` used more than one core. The tools agree on 34,146,531 of 34,148,222 consequence tuples on shared transcripts. Executable VEP sides with DuckVEP on 1,689 of the 1,691 disagreements. Methods and runner: [vep-rs comparison](benchmarks/benchmark_duckvep_vep_rs.md).

## Whole-haplotype input domain (issue #12)

- `duckvep_haplotypes` accepts complete phased calls of any ploidy, literal alleles of any length, and transcripts using supported genetic codes. Start and stop tests use the transcript's genetic code: table 1 uses the ATG start rule; other tables use their start-codon sets, and a change between start codons is not a peptide change.
- In the recorded HG002 comparison, one of 157,986 paths differed: a 94-base insertion is classified as a frameshift with a gained stop. Hand-derived tests cover a 51-base insertion and four edits under NCBI table 2; a metamorphic test compares haploid, triploid and tetraploid calls with diploid lanes carrying the same edits.

## Haplotype policy identifiers

- Output identifies the prediction policy as `duckvep-coding`, the NMD rule as `ejc50`, and the VEP-compatible phase policy as `vep_compat`.

## VEP conformance campaigns for long alleles

- The random-allele generator `generate_witnesses.R --max-random-length` defaults to 100 bases; the differential runner accepts alleles of 101 bases, including the VCF anchor. Two campaigns of 30,268 alleles each (seeds 173 and 29, approximately 15,100 alleles over 50 bases per campaign) matched executable VEP 116 on consequences, with each HGVS string either matching or absent on both sides.
- The runner uses `DUCKVEP_READER_EXT` to load a DuckHTS build for VCF, GFF and FASTA readers. The field projection and replay inputs used by the scale runner are `benchmarks/duckvep_field_projection.R` and `benchmarks/data/scale_contracts/field_replay_9bf888e`.

## Model snapshots

- `duckvep_model_save(name, path)` writes a loaded model's native arrays to one file. `duckvep_model_restore(name, path)` validates and maps that file read-only. For Ensembl 116 GRCh38 (644,427 transcripts, one thread), the recorded restore took 0.54 s with 1.3 GiB peak RSS; relation loading took 3.0 s and 3.4–5.2 GiB. Processes restoring the same file share mapped pages.
- Restore verifies the header, section bounds, full-file checksum, model coordinates and reference FASTA. Wrapped circular models are lifted on restore, and mapped bytes count against the native budget. Platforms without file mapping read the file into one block.

## DuckDB C API v2 preview host (issue #8)

- `make release_v2` builds a v2 preview extension for DuckDB's stable C API v2, alongside the v1 release build. It uses the pinned preview SDK in `duckdb_capi_v2/` and `duckvep-package.json`; its footer is `C_STRUCT`, its extension API version is `v2.0.0`, and loading it executes no SQL.
- The recorded v2 preview-host campaign covered 27 public functions, with behavior and outputs compared against v1 in [`docs/v2-host.md`](docs/v2-host.md). A v1 comparison of `duckvep_so_terms`, `duckvep_allele_geometry` and `duckvep_breakend_geometry` covered 52 test cases.

## Native resource control on the v2 host (issue #8)

- `duckvep_native_budget`, `duckvep_native_budget_set`, `duckvep_native_budget_reset_high_water` and `duckvep_worker_limits_set` use the process-wide budget on the v2 preview host. Model publication, annotation workers, haplotype scans and the coding-calls reader share the budget. A capacity refusal publishes no model, leaves the connection usable and reports an explicit error (`test/sql_v2/v2_budget.sql`).

## Haplotype capture on the v2 host (issue #8)

- `duckvep_haplotypes`, `duckvep_coding_transcripts` and `duckvep_coding_calls` run on the v2 preview host. Because v2 has no private connection, `duckvep_haplotype_load_sql` returns caller-side statements. The `duckvep_stage` COPY format captures rows for one `JOB` into spillable column collections, and `duckvep_haplotype_scan('<job>')` replays them once.
- Full-row hashes match v1 across vertical, same-codon, frame, start/stop and NMD suites, 47 policy, limit and error cases, and discovery and `coding_calls` fixtures. HG002 through `duckvep_coding_calls` produced 157,986 rows with checksum `1456007180270799092358516`. Details: [`docs/v2-host.md`](docs/v2-host.md).

## Loss-of-function relation (issue #41)

- `duckvep_lof_sql(annotations, transcripts, reference [, options])` returns per-variant, per-transcript LOFTEE-style `lof`, `lof_filter`, `lof_flags`, `lof_info` and `lof_unchecked` fields from annotation, transcript and reference relations. Implemented rules follow `konradjk/loftee` at `a46b502`; missing resources are named in `lof_unchecked`, and MaxEntScan splice predictions are not implemented.
- With GRCh38 chr21, ancestor, PhyloCSF and GERP weighting off, comparison with VEP 116 and LOFTEE matched `LoF`, `LoF_filter`, `LoF_flags` and `LoF_info` for 59,668 pairs from 6,845 HG002 and ClinVar variants, and 2,499 pairs from a 177-allele two-strand set. See [LoF evidence](benchmarks/duckvep_lof.md).

## Annotation and projection on the v2 host (issue #8)

- `duckvep_annotate_sql`, `duckvep_annotate_projected_sql` and `duckvep_transcript_projection_sql` execute on the v2 preview host. Their output matches v1 for compact, rich, HGVS with a reference FASTA, projected and mixed structural/breakend events across three models, plus 11,100-row multi-vector runs.
- A separate v1-host measurement used 100,000 GIAB sites, compact output, one pinned core and five alternating runs: median throughput was 1,040 ns per variant against 1,080 ns per variant for main, with identical rows.

## Model loading on the v2 host (issue #8)

- `duckvep_model_load_sql()` returns statements for caller-side execution. The v2 preview host stages relation rows with `COPY ... (FORMAT duckvep_stage)` in the caller's transaction, so temporary tables and uncommitted rows are visible; `duckvep_model_publish()` builds and installs the model. A failed or cancelled COPY stages nothing, and a failed publish installs nothing. Byte-equal model fingerprints were checked for three fixture models across hosts. The host-specific usage is documented in [`docs/v2-host.md`](docs/v2-host.md).

## SQL builders on the v2 host (issue #8)

- The v2 preview host serves the twelve SQL builders: `duckvep_ensembl_*_sql`, `duckvep_model_receipt_sql`, `duckvep_annotate_sql`, `duckvep_annotate_projected_sql`, `duckvep_transcript_projection_sql` and the five `duckvep_prepare_*_sql` functions. Their generated SQL matches v1 in a 192-case matrix covering option keys, defaults, invalid options, identifier quoting and NULL arguments. Preparation queries return the same rows on both hosts.

## Additional row functions on the v2 host (issue #8)

- `duckvep_repeat_alleles`, `duckvep_phase_call`, `_duckvep_revcomp`, `_duckvep_raw_gt` and `_duckvep_record_order` produce v1-matching results in 126 cases, including nested output beyond one vector, NULLs inside lists and structs, and errors. `_duckvep_revcomp` and `_duckvep_raw_gt` return NULL children for NULL inputs.

## Haplotype scale qualification (issue #2)

- `duckvep_coding_transcripts(model, seq_region, position, reference, alternate)` returns resident-model transcript ordinals touched by a normalized coding event. It matches the annotation builder's CDS-overlap pairs, including its shared-anchor, interbase-insertion and short-intron rules. On whole HG002 it matched 247,374 pairs with no missing or extra pair.
- CDS and protein difference alignment measured 23.3 s to 3.4 s (6.8×) on the HG002 genome. Its widening band proves the full-matrix optimum, including traceback ties. `max_alignment_cells` applies to the band required by the alignment; an unfit band is rejected with its cell count. A 20 kb coding sequence with a deletion and substitution used 322 million cells under the positional bound and 180 thousand for the required band; the longest human transcripts required 11.7 billion full-matrix cells, above the native budget.
- In Ensembl 116 GRCh38 model-load measurements, load time was 3.5 s to 2.5 s and the high-water mark was 2.7 GiB to 1.6 GiB. In the cold mode-B HG002 qualification with identical input, model, core and caps, DuckVEP took 12.58 s and `bcftools csq` 19.66 s (ratio 0.640; the 0.5 gate was not met); warm execution measured 0.429. A 5,000,000-record single-sample MANE job took 16.4 s within the 16 GiB, 4 GiB native and 8 GB DuckDB caps. See the [qualification record](benchmarks/data/haplotype_scale/README.md).

## Circular sequence regions (issue #7)

- Transcripts, exons, regulatory features and motifs that cross a circular region's origin are annotated on lifted intervals. Rotation tests compare consequences, projected edits, peptides, NMD and HGVS; circular models without wrapped objects, including human MT, match the linear path byte for byte. Ensembl Genomes 63 includes 19 public origin-crossing transcripts, 18 bacterial or archaeal.
- Structural, breakend and phased edit-set entry points reject models containing wrapped circular objects. VEP 116 models origin-crossing transcripts as reversed-bound intervals and is not an oracle for those cases. In three cached genomes, VEP agrees on the other transcripts except flank rows available only through the origin; rotation and linear-model checks provide the crossing-transcript evidence. The survey is in [`benchmarks/data/circular_source_survey.md`](benchmarks/data/circular_source_survey.md).
- An origin-focused throughput test measured a 17–20% one-thread cost for interval lifting; output checksums matched across thread counts and rotations.

## Enforced native memory budget (issue #3)

- Native owners allocate through one process-wide atomic budget, default **4 GiB**. Blocks are charged before allocation, large blocks use page-rounded capacity, and a growing `realloc` keeps the old block charged. Refusal reports requested, used and permitted bytes; model publication fails without disturbing loaded models.
- `duckvep_native_budget()` reports current and high-water bytes per owner. `duckvep_native_budget_set(bytes)` changes the ceiling, and `duckvep_native_budget_reset_high_water()` resets high-water marks.
- Annotation allows at most **6** concurrent workers. Each worker has a **128 MiB** scratch lease and a **256 MiB** emitted-output allowance; idle workers retain at most **64 MiB**. The measured emitted-output peak was 110–128 MiB. `duckvep_worker_limits_set(workers, scratch_bytes, emit_bytes, idle_bytes)` configures these limits. An allele or vector exceeding its lease returns a capacity error.
- `make test_fault_injection` tests allocation failures under AddressSanitizer and LeakSanitizer. htslib buffers are outside this budget; each open FASTA index uses a fixed reservation.

## Native SQL builder interface

DuckVEP registers functions without modifying the database catalog during `LOAD`. Preparation, annotation and projection run SQL emitted by native scalar builders in the caller's connection. Existing macro calls must use the builder interface; `LOAD` does not install compatibility macros.

| Before | After |
| --- | --- |
| `FROM duckvep_annotate('t', 'model', hgvs := true)` | `FROM query(duckvep_annotate_sql('t', 'model', {hgvs: true}))` |
| `FROM duckvep_ensembl_regions('core', 'ref', 'GRCh38')` | `FROM query(duckvep_ensembl_regions_sql('core', 'ref', 'GRCh38'))` |
| `FROM duckvep_ensembl_transcripts('core', 'ref', 'GRCh38')` | `FROM query(duckvep_ensembl_transcripts_sql('core', 'ref', 'GRCh38'))` |
| `FROM duckvep_ensembl_regulation_features('funcgen', 'regions')` | `FROM query(duckvep_ensembl_regulation_features_sql('funcgen', 'regions'))` |
| `FROM duckvep_transcript_projection('events', 'annotations', 'transcripts')` | `FROM query(duckvep_transcript_projection_sql('events', 'annotations', 'transcripts'))` |
| `FROM duckvep_model_receipt('regions', 'transcripts', ...)` | `FROM query(duckvep_model_receipt_sql('regions', 'transcripts', ...))` |

`duckvep_annotate_sql()` requires a model name. Optional builder settings are fields of a trailing STRUCT. `Rduckvep` provides matching `rduckvep_*_sql()` wrappers and `rduckvep_annotate()` for a data-frame result.

## Species structural evidence and builder behavior (issue #4)

- `test/duckvep/conformance/species_sv_vep116_differential.R` compares DEL, DUP, DUP:TANDEM, INV, symbolic INS, paired BND and structural HGVS with digest-pinned VEP 116 for mouse GRCm39, fly BDGP6.54 and Arabidopsis TAIR10. Results, pins and retained BND examples are in [`benchmarks/data/sv_species_evidence.md`](benchmarks/data/sv_species_evidence.md).
- `duckvep_prepare_breakend_pairs_sql()` resolves mates with a hash join. On 100,000 records the measured run took 20 s with the nested-loop implementation and was about 100 times faster with the hash join, with identical output.
- `duckvep_prepare_sv_geometry_sql()` reports extreme POS/END/CIPOS/CIEND values as `unsupported_geometry` rather than failing the batch with an INT64 overflow.
- `test/sql/duckvep_builder_errors.test` and `Rduckvep`'s `test_builder_errors.R` cover NULL, malformed and extreme inputs for all five STR/SV/BND/HGVS builders. `benchmarks/benchmark_sv_builders.R` records time and peak RSS.

## Breakend identity and structural HGVS (issue #4)

- `duckvep_prepare_breakend_pairs_sql()` validates MATEID/EVENT identity, mate coordinates, orientation and inserted-sequence length. It returns one row per physical VCF record with a stable reason; mate records are not merged. Fusion, phase and inserted-only sequence are not inferred from ALT syntax.
- `duckvep_prepare_breakend_fusion_sql()` joins identity evidence to caller-supplied endpoint genes and reports partner-gene candidates per physical record. It does not assert a fusion or frame.
- `duckvep_prepare_structural_hgvs_sql()` emits unshifted genomic HGVS and an equivalent literal edit for exact-span symbolic DEL, DUP and INV. Other alleles return `unavailable` or `unsupported` with a reason. VEP 116 emits no HGVS for symbolic structural alleles or breakends; see [`docs/structural-identity-hgvs.md`](docs/structural-identity-hgvs.md).
- `duckvep_breakend_geometry()` returns NULL child fields for non-breakend input.
