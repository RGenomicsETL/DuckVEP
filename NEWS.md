# DuckVEP

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
