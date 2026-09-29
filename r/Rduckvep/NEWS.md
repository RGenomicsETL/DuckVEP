# Rduckvep

- The bundled extension annotates models that contain origin-crossing transcripts, exons or
  regulation features on circular sequence regions (lifted-interval execution); see the project
  NEWS.md. Structural, breakend and phased entry points still refuse such models.

- `rduckvep_prepare_breakend_pairs()`, `rduckvep_prepare_breakend_fusion()` and
  `rduckvep_prepare_structural_hgvs()` wrap the native BND identity, endpoint-gene
  and structural HGVS builders. No parsing or HGVS logic lives in R.
- `rduckvep_haplotypes()` adds versioned `prediction_policy`, `prediction_status` and
  `prediction_reason` (`duckvep-coding-v1` eligibility only; no SO/IMPACT/NMD yet),
  keyed `carrier_predictions`, `contributor_provenance` and `normalized_edits`. Existing
  columns are unchanged and `nominal_length_diff` stays last.
- The same-codon classifier (coding-v1 slice 3) adds `haplotype_consequences` and `haplotype_impact`
  before `nominal_length_diff`, and a `predicted` status. It classifies the combined haplotype for
  eligible paths with no frame, start or stop effect: synonymous (LOW), missense, inframe insertion or
  deletion, or protein-altering (MODERATE). Other eligible paths stay `eligible_classifier_pending`,
  now with reason `frame_classifier_pending` or `start_stop_classifier_pending`.
- The frame opening/restoration classifier (coding-v1 slice 4) decides frame-shifting and stop-gain paths from the
  translated haplotype: `stop_gained` for a new first stop before the reference terminator, `frameshift_variant`
  when that stop intersects a displaced-frame interval or the frame is still displaced when the CDS runs out,
  and, for a frame restored before termination, `synonymous_variant` or `protein_altering_variant`. Edits after
  the first stop stay contributors. `carrier_predictions` gains `haplotype_impact` and `haplotype_consequences`
  per carrier, so an ineligible carrier no longer hides the set of eligible carriers of the same row. The
  reason `frame_classifier_pending` no longer occurs; only start/terminal-codon paths stay
  `start_stop_classifier_pending`.
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
