# DuckVEP

## Haplotype scale qualification (issue #2, slice 7)

- `duckvep_coding_transcripts(model, seq_region, position, reference, alternate)` returns the ascending model ordinals of the transcripts the VCF event touches in coding sequence, from the resident interval index. It finds exactly the (event, transcript) pairs of the annotation builder's CDS overlap: the event is normalized inside the scalar as the builder does (a shared anchor base creates no overlap, an insertion is an interbase point, introns of at most 13 bases inside the CDS count as coding), so a whole-genome VCF reduces to its coding records before any per-record work. Records in an intron, UTR, flank or intergenic sequence, alleles that are not literal bases and events without a difference return an empty list. Integer arguments of any width and signedness bind. Lifted circular models are refused, like `duckvep_haplotypes`. It stays a separate scalar next to `duckvep_haplotypes` (whose input is a caller-staged relation of calls); `rduckvep_coding_transcripts()` is its R wrapper. On whole HG002 it matches the builder on 247,374 pairs with none missing or extra.
- The CDS and protein difference alignment of `duckvep_haplotypes` is 6.8 times faster on the HG002 genome (23.3 s to 3.4 s of prediction). Its band starts at the length change and widens until the optimum is proved to lie inside it, so the traceback, ties included, is that of the full matrix; rows inside the common prefix are closed-form; equal-length comparisons scan for differing runs. Output is byte-identical. `max_alignment_cells` is now checked against the band an alignment needs, not against the band of a feasible positional bound: a 20 kb coding sequence with a deletion and a substitution needed 322 million cells and needs 180 thousand, and the longest human transcripts needed 11.7 billion (more than the native budget), so long transcripts with several edits no longer fail. A band that is still to be tried and does not fit is refused with an explicit `max_alignment_cells` error naming its cell count.
- Model load validates stored sequences with a table, computes the first-stop cache with a table lookup per codon, and no longer copies the sequence pools a second time when rows arrive in transcript order: 3.5 s to 2.5 s for Ensembl 116 GRCh38, and a load high-water mark of 1.6 GiB instead of 2.7 GiB.
- `benchmarks/haplotype_scale/` and `benchmarks/data/haplotype_scale/` hold the scale qualification: on identical HG002 input, transcript model, core and caps, cold DuckVEP mode B takes 12.58 s against 19.66 s for `bcftools csq` (ratio 0.640, the 0.5 gate is not met; warm execution is 0.429), and a 5,000,000-record single-sample MANE job runs in 16.4 s inside the 16 GiB, 4 GiB native and 8 GB DuckDB caps.

## Circular sequence regions (issue #7, slice 2)

- Models with origin-crossing transcripts, exons, regulatory or motif features can now be pinned and annotated. A circular region that carries such an object executes on a lifted linear interval: positions shift by a multiple of the region length, every object is admitted at three images, wrapped spans become `[s, e + L]`, and one row is kept per event and source object. Consequences, projected edits, peptides, NMD and HGVS are identical under rotation of reference, model and events, and equal an ordinary linear model. Reference windows that wrap are fetched from the existing worker scratch. `model_sha256` is unchanged; circular regions without a wrapped object, human MT included, produce byte-identical output.
- Structural, breakend and phased edit-set entry points still refuse a model with a wrapped circular object, with an explicit error.
- Circular coordinates are independent of the mitochondrial codon table: `codon_table` selects a translation rule and says nothing about topology.
- Public origin-crossing transcripts exist (19 in Ensembl Genomes 63, 18 of them bacterial or archaeal) but VEP 116 models them as intervals with reversed bounds, so it is not an oracle for them. On three genomes with a VEP cache, DuckVEP equals VEP on every other transcript except flank rows that exist only through the origin; the crossing object is property-proved and linear-model-proved. See `benchmarks/data/circular_source_survey.md`, `ERRATA.md` and `scripts/circular_vep_differential.py`.
- `benchmarks/duckvep_circular_origin.py` records one- and multi-thread throughput, output equality and peak memory on an origin-focused workload: lifting costs about 17-20% of one-thread throughput and the output checksum is identical across threads and rotations.
- `make test_properties` builds and runs the native theft/greatest properties (`test_properties_sanitized` adds ASan and UBSan).

## Enforced native budget and bounded execution (issue #3)

- Every native owner (model arrays, interval indexes, reference readers, workspaces, per-worker result and text arenas, SQL builders) allocates through one process-wide atomic budget, default **4 GiB**. Blocks are charged before allocation, with page-rounded capacity for large blocks and the old block still charged while a `realloc` grows. A refusal is an explicit `capacity error: ... budget exceeded (requested N bytes, M in use, limit L)`; nothing is truncated. A model load that exceeds the budget publishes nothing and leaves every loaded model usable.
- `duckvep_native_budget()` reports current and high-water bytes per owner; `duckvep_native_budget_set(bytes)` sets the ceiling (it cannot go below the bytes already charged); `duckvep_native_budget_reset_high_water()` restarts the high-water marks.
- Annotation admits at most **6** concurrent workers; each holds a **128 MiB** native scratch lease and a **256 MiB** emitted-output allowance (raised from 64 MiB, which complete-17 exceeded on gene-dense input; the peak measured is 110-128 MiB), charged as it grows, and an idle worker keeps at most **64 MiB**. `duckvep_worker_limits_set(workers, scratch_bytes, emit_bytes, idle_bytes)` changes them. An allele or vector over its lease is a capacity error.
- `make test_fault_injection` builds an AddressSanitizer + LeakSanitizer extension whose allocator can fail the Nth allocation and fails every allocation of model load and annotation in turn. htslib's own buffers are not routed through the budget; a fixed reservation stands in for each open FASTA index.

## Breaking change: native SQL builders

DuckVEP registers functions without modifying the database catalog during `LOAD`. Preparation, annotation and projection run SQL emitted by native scalar builders in the caller's connection. Existing macro invocations must be migrated; `LOAD` does not install compatibility macros.

| Before | After |
| --- | --- |
| `FROM duckvep_annotate('t', 'model', hgvs := true)` | `FROM query(duckvep_annotate_sql('t', 'model', {hgvs: true}))` |
| `FROM duckvep_ensembl_regions('core', 'ref', 'GRCh38')` | `FROM query(duckvep_ensembl_regions_sql('core', 'ref', 'GRCh38'))` |
| `FROM duckvep_ensembl_transcripts('core', 'ref', 'GRCh38')` | `FROM query(duckvep_ensembl_transcripts_sql('core', 'ref', 'GRCh38'))` |
| `FROM duckvep_ensembl_regulation_features('funcgen', 'regions')` | `FROM query(duckvep_ensembl_regulation_features_sql('funcgen', 'regions'))` |
| `FROM duckvep_transcript_projection('events', 'annotations', 'transcripts')` | `FROM query(duckvep_transcript_projection_sql('events', 'annotations', 'transcripts'))` |
| `FROM duckvep_model_receipt('regions', 'transcripts', ...)` | `FROM query(duckvep_model_receipt_sql('regions', 'transcripts', ...))` |

`duckvep_annotate_sql()` requires the model name. Optional builder settings are fields of a STRUCT rather than named function parameters. `Rduckvep` exposes matching `rduckvep_*_sql()` wrappers and `rduckvep_annotate()` for a data-frame result.

## Species structural evidence and builder robustness (issue #4 item 6)

- New executable differential `test/duckvep/conformance/species_sv_vep116_differential.R` compares DEL, DUP, DUP:TANDEM, INV, symbolic INS, paired BND and structural HGVS against digest-pinned VEP 116 for mouse GRCm39, fly BDGP6.54 and Arabidopsis TAIR10. Results, pins and the retained BND counterexamples are in `benchmarks/data/sv_species_evidence.md`.
- `duckvep_prepare_breakend_pairs_sql()` resolved mates with a nested-loop join and was quadratic (100,000 records took 20 s); it is now a hash join and about 100 times faster there, with identical output.
- `duckvep_prepare_sv_geometry_sql()` no longer fails a whole batch with an INT64 overflow error on extreme POS/END/CIPOS/CIEND values; such rows return `unsupported_geometry`.
- `test/sql/duckvep_builder_errors.test` and `Rduckvep`'s `test_builder_errors.R` cover NULL, malformed and extreme input for all five STR/SV/BND/HGVS builders; `benchmarks/benchmark_sv_builders.R` measures their time and peak RSS.
- The human and species differential scripts read model, FASTA, cache and extension paths from `DUCKVEP_GRCH38_MODEL`, `DUCKVEP_GRCH38_FASTA`, `DUCKVEP_GRCH38_VEP_CACHE`, `DUCKVEP_GRCH37_FASTA` and `DUCKVEP_EXTENSION_FILE`, keeping the previous values as defaults. `scripts/run_species_vep116_docker.sh` honors `VEP_BUFFER_SIZE`.

## BND record identity and structural HGVS (issue #4 items 3 and 4)

- `duckvep_prepare_breakend_pairs_sql()` validates MATEID/EVENT identity, mate coordinates, mate orientation and inserted-sequence length above `duckvep_breakend_geometry()`. It returns one row per physical VCF record with a stable `reason`; a mate is looked up, never merged. Fusion, phase and inserted-only sequence are never asserted from ALT syntax.
- `duckvep_prepare_breakend_fusion_sql()` joins that identity with caller-supplied endpoint genes and reports partner-gene evidence per physical record (`fusion_asserted` is always false).
- `duckvep_prepare_structural_hgvs_sql()` emits unshifted genomic HGVS and an equivalent literal edit for exact-span symbolic DEL, DUP and INV; every other allele is `unavailable` or `unsupported` with a reason. VEP 116 emits no HGVS for symbolic structural alleles or breakends. See `docs/structural-identity-hgvs.md`.
- `duckvep_breakend_geometry()` returned NULL STRUCT rows with uninitialised children for non-breakend input, which crashed when the vector was copied (for example under `TRY`); children are now NULL.
