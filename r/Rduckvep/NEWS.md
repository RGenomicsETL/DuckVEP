# Rduckvep

- `rduckvep_prepare_breakend_pairs()`, `rduckvep_prepare_breakend_fusion()` and
  `rduckvep_prepare_structural_hgvs()` wrap the native BND identity, endpoint-gene
  and structural HGVS builders. No parsing or HGVS logic lives in R.
- `rduckvep_haplotypes()` adds versioned `prediction_policy`, `prediction_status` and
  `prediction_reason` (`duckvep-coding-v1` eligibility only; no SO/IMPACT/NMD yet),
  keyed `carrier_predictions`, `contributor_provenance` and `normalized_edits`. Existing
  columns are unchanged and `nominal_length_diff` stays last.
- The extension's SQL preparation and annotation entry points are native builders:
  `query(duckvep_annotate_sql('events', 'model', {hgvs: true}))` replaces
  direct annotation macro calls. `rduckvep_annotate_sql()` and the other
  `rduckvep_*_sql()` functions invoke those builders; `rduckvep_annotate()`
  returns annotated events as a data frame. See the project NEWS.md for the
  migration table.

# Rduckvep 0.1.0

- First release as a separate package. DuckVEP and its R front end were extracted
  from DuckHTS/Rduckhts with their history; `rduckhts_haplotypes()` is now
  `rduckvep_haplotypes()`, and `rduckvep_connect()`/`rduckvep_load()` replace the
  Rduckhts connection helpers for DuckVEP work.
