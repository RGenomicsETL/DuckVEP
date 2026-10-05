# DuckVEP extraction checklist for DuckHTS

This checklist scopes the migration of DuckVEP-owned assets from DuckHTS and identifies the DuckHTS interfaces to retain. Source of truth: DuckHTS `origin/develop` at the extraction date. The [path manifest](duckhts-removal-paths.txt) enumerates 2,230 paths to remove or repoint. The [history filter](extraction-paths.txt) identifies four shared test fixtures retained in DuckHTS.

## Native extension

- Remove `src/duckvep/` (all kernels and adapters) and `src/include/duckvep_sql.h` from the source tree and the corresponding rows in `src/duckhts_sources.tsv` (lines 43–66). Retain `src/vep_parser.c` (line 42): it parses the `CSQ` field in `read_bcf`.
- In `src/duckhts.c`, remove the four DuckVEP extern declarations (lines 91–94), `register_duckvep_functions` and `register_duckvep_sql_kernels` (263–267), and `register_duckvep_sql_functions` and `register_duckvep_ensembl_functions` (1170–1171); preserve the non-VEP registration chain.
- In `CMakeLists.txt`, remove the DuckVEP kernel include directories (150–151), but retain the `third_party/cgranges` source and include path (144, 152): `src/cgranges_api.c` uses it independently. Retain htslib for BCF/FASTA readers and `vep_parser.c`.

## SQL and evidence

- Remove/repoint the 14 `duckvep_*` entries in `functions.yaml`: `duckvep_ensembl_regions`, `duckvep_ensembl_transcripts`, `duckvep_ensembl_regulation_features`, `duckvep_model_receipt`, `duckvep_model_load`, `duckvep_model_drop`, `duckvep_allele_geometry`, `duckvep_transcript_projection`, `duckvep_repeat_alleles`, `duckvep_breakend_geometry`, `duckvep_haplotypes`, `duckvep_phase_call`, `duckvep_annotate`, `duckvep_so_terms`. Regenerate `r/Rduckhts/inst/function_catalog/functions.yaml`, `functions.tsv`, `functions.md` and `reference.md` with the DuckVEP entries removed; keep the BCF/CSQ reader documentation.
- Move or link `test/sql/duckvep_*.test`, `test/duckvep/`, DuckVEP-only `test/data/` and `test/scripts/`, `pipelines/duckvep/`, `benchmarks/*duckvep*`, `benchmarks/fastvep_*`, `benchmarks/data/duckvep_*`, `design/duckvep*.md`, and DuckVEP-only scripts. The enumerated paths are in the manifest. Replace DuckVEP Makefile recipes and DuckVEP branches in `scripts/run_sqllogictest.py` and `scripts/test_sanitized_extension.sh` with sibling-repository invocations when those workflows include cross-extension tests.
- Keep DuckVEP compatibility errata in the DuckVEP repository; replace DuckHTS's top-level `ERRATA.md` with a link to that source.
- Repoint the root `README.Rmd` DuckVEP section (lines 199–1446) and regenerate `README.md`; point readers to the sibling repository rather than presenting the model as a DuckHTS-owned extension. Repoint `r/Rduckhts/README.Rmd` lines 457–605 and regenerate its README.

## R and CI

- Drop `r/Rduckhts/R/haplotypes.R` and `r/Rduckhts/man/rduckhts_haplotypes.Rd` in favor of `Rduckvep::rduckvep_haplotypes`; repoint users of `rduckhts_haplotypes`. Remove `test_duckvep.R`, `test_duckvep_haplotypes.R`, and `test_duckvep_phase.R` from `r/Rduckhts/inst/tinytest/`, and remove DuckVEP assertions in `test_basic.R` and `test_connection.R` while retaining their other assertions. Remove the DuckVEP source-bundling branches in `r/Rduckhts/R/bootstrap.R` and `duckhts_duckvep_kernel_source_files` in `r/Rduckhts/R/source_manifest.R`; preserve core source packaging.
- Repoint `.github/workflows/duckvep-provenance.yml` to DuckVEP or delete the DuckHTS copy. Retain `.github/workflows/MainDistributionPipeline.yml` for DuckHTS core, adjusting VEP-only sanitizer or benchmark jobs if present. Keep the DuckHTS `read_bcf` function and `vep_parser.c` CSQ reader; an optional DuckVEP integration test can then load both extensions without overlapping function registrations.
- Keep `test/data/ce.fa`, `test/data/ce.fa.fai`, `test/data/geno_phase_partial.vcf` and `test/data/geno_vcf44.vcf`: core FASTA and genotype tests use them. Keep `third_party/cgranges`, htslib, `src/cgranges_api.c` and the DuckHTS sequence helpers.
