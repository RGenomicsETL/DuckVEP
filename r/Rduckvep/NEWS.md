# Rduckvep

- `rduckvep_coding_calls()` wraps the new SQL table function `duckvep_coding_calls(model, path)`: a fused native reader of a
  (bgzipped) VCF or BCF that discards records outside coding sequence before decoding genotypes and returns the calls relation of
  `rduckvep_haplotypes()` (same discovery as `rduckvep_coding_transcripts()`, output identical to building the calls in SQL). The
  model must have been loaded with `seq_region_name` in its regions query.
- `rduckvep_coding_transcripts()` wraps the SQL scalar `duckvep_coding_transcripts(model, seq_region, position, reference, alternate)`:
  the (event, transcript) pairs of the annotation builder's CDS overlap, found from the resident interval index (the fast
  discovery route for haplotype calls on whole-genome input). The bundled extension also makes `rduckvep_haplotypes()`
  aligns CDS and protein differences much faster with unchanged output. `max_alignment_cells` is now checked against the
  band an alignment needs rather than a feasible bound, so long transcripts with several edits fit the default; see the
  project NEWS.md and `benchmarks/data/haplotype_scale/README.md`.

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
- The start/stop classifier (coding-v1 slice 5) decides the remaining paths: `start_lost` alone when the
  edited CDS does not begin with ATG (initiation and NMD unknown, other predictions suppressed), `stop_lost`
  when no stop is left (plus `frameshift_variant` for a frame still displaced at the CDS end, which previously
  gave `frameshift_variant` alone; no downstream extension is invented), and `stop_retained_variant` (LOW) for
  a changed terminal codon that is still a stop with an unchanged peptide. No eligible path stays
  `eligible_classifier_pending` and the reason `start_stop_classifier_pending` no longer occurs. Contributor
  `role` is now assigned per edit: an edit that starts after the first stop is `post_stop` even when it shares an
  interaction block with an earlier edit (a restoring deletion after an early stop was `applied`).
- NMD attribution (coding-v1 slice 6, rule `ejc50-v1`) adds `nmd_rule`, `nmd_prediction`, `nmd_stop_position`,
  `nmd_junction_position` and `nmd_contributors` to `rduckvep_haplotypes()` (before `nominal_length_diff`, which stays
  last), and `nmd_prediction`, `nmd_stop_position` and `nmd_junction_position` to each `carrier_predictions` row. For a newly
  premature stop (`stop_gained`) the prediction is `trigger` when J - S > 50 and `escape` otherwise, with S the final
  nucleotide of the first stop codon and J the final nucleotide of the penultimate exon, both in edited spliced-transcript
  coordinates (indels upstream or inside the penultimate exon move J; a single-exon transcript always escapes and has no
  J). Known termination without a new premature stop is `not_applicable`; `start_lost`, an edited CDS with no stop (a
  frame that runs off the CDS, `stop_lost`), unresolved exon topology and every failed or ineligible path are `unknown`.
  `nmd_contributors` lists the applied contributors (the edits up to and including the stop) plus the `post_stop`
  indels that moved J (changed length at or before the penultimate exon's last base; their `role` stays `post_stop`), and is
  NULL unless the prediction is `trigger` or `escape`. The
  prediction is decided on the whole haplotype, not per allele, and is an EJC-distance heuristic only: no reinitiation,
  no long-exon exception, no `NMD_transcript_variant` biotype term. Existing columns and their values are unchanged.
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
