# BND record identity, fusion evidence and structural HGVS

Native SQL builders (safe quoting, trailing typed options STRUCT, deterministic text). `Rduckvep` only wraps them.

## `duckvep_prepare_breakend_pairs_sql(events[, {}])`

Input columns: `event_index, chrom, pos, id, ref, alt, info` (raw VCF records). Output: exactly one row per input record, ordered by `event_index`. A record and its mate are separate rows; the shared `pair_key` (`min|max` event_index) exists only when both agree.

`record_kind`: `paired_breakend`, `single_breakend`, `malformed_alt`, `not_breakend`, parsed by `duckvep_breakend_geometry()`.

`status` and stable `reason` (first failing check in this order):

| status | reason |
| --- | --- |
| `not_applicable` | `not_breakend`, `single_breakend` |
| `invalid` | `malformed_alt` |
| `unproven` | `missing_id`, `missing_mateid`, `mate_not_found`, `mate_missing_mateid` |
| `conflict` | `single_breakend_with_mateid`, `duplicate_id`, `duplicate_mateid`, `self_mate`, `duplicate_event`, `ambiguous_mate_id`, `mate_not_paired_breakend`, `mate_duplicate_mateid`, `mate_not_reciprocal`, `mate_coordinate_conflict`, `mate_orientation_conflict`, `mate_insert_conflict`, `event_conflict` |
| `reciprocal` | `reciprocal` |

Evidence columns stay separate so partial agreement is visible: `id_reciprocal` (MATEID names a unique record that names this one), `coordinate_reciprocal` (each ALT mate chrom:pos equals the other record's chrom:pos), `orientation_reciprocal` (VCF 4.3 section 5.4.9 rule: the mate's kept side is the opposite of this record's, i.e. `t[p[` pairs with `]p]t`, `t]p]` with `t]p]`, `[p[t` with `[p[t`), `insert_agree` (equal inserted-only length), `event_agree`. `inserted_sequence`/`inserted_length` exclude the anchor base. `phase_status` is always `not_evaluated` and `fusion_status` `not_asserted`.

Inputs: `event_index` is cast to BIGINT and a non-integral `pos` becomes NULL.

## `duckvep_prepare_breakend_fusion_sql(pairs, genes[, {}])`

`pairs` is the relation produced above; `genes` has `event_index, gene_id` (one row per gene overlapped by that physical endpoint, produced by the caller from `duckvep_annotate`). One row per physical record, with `endpoint_genes` and `mate_endpoint_genes`. `status`: `identity_unproven` (MATEID/coordinate/EVENT evidence incomplete; `reason` carries the identity reason), `endpoint_without_gene`, `shared_gene_endpoints`, `candidate_orientation_conflict` (identity proven, orientation contradicts the VCF rule) and `candidate_partner_genes`. `fusion_asserted` is always false and `phase_status` `unproven`: partner genes are evidence, not a fusion, frame or transcript.

## `duckvep_prepare_structural_hgvs_sql(events, reference[, {max_span}])`

`events`: `event_index, chrom, pos, ref, alt, info`. `reference`: `event_index, reference_sequence` = FASTA[pos..END], plus (for DUP) an equal-length following flank, FASTA[pos..END+(END-pos)]; DEL and INV accept either. `max_span` defaults to 5000 (1 to 60000).

Supported domain (`hgvs_status = 'supported'`): precise symbolic `<DEL>`, `<DUP>`, `<DUP:TANDEM>` and `<INV>` with integral `INFO/END`, no `CIPOS`/`CIEND`/`IMPRECISE`, and a matching reference. Output: `hgvs_g` (`chrom:g.S_Edel|dup|inv`, `normalization = 'none'`, no 3' shifting), and the equivalent literal edit (`literal_position`, `literal_reference`, `literal_alternate`) that the existing small-variant path turns into `c.`/`n.` HGVS (`transcript_hgvs_route = 'literal_equivalent'`). Nothing is derived from consequence labels.

`unavailable` (an input needed for the description is absent or unusable): `missing_field`, `invalid_pos`, `missing_chrom`, `missing_end`, `duplicate_end`, `invalid_end`, `missing_reference_sequence`, `ambiguous_reference_sequence`, `reference_length_mismatch`, `reference_alphabet`, `invalid_ref`, `reference_anchor_mismatch`, `missing_flank_sequence`.

`unsupported` (well formed but outside the domain): `breakend` (paired and single), `symbolic_insertion`, `copy_number`, `repeat_expansion`, `symbolic_allele`, `literal_allele` (use the small-variant HGVS path), `alt_syntax`, `imprecise`, `span_capacity`, `duplication_adjacent_repeat`, `inversion_identity`, `inversion_not_reducible`.

## Evidence

VEP 116 (digest-pinned Docker) emits no `hgvsg`, `hgvsc` or `hgvsp` for symbolic DEL/DUP/INV/INS/CNV alleles or BND pairs, so those states are `unsupported` for BND/INS/CNV and the supported domain is validated on the equivalent literal edit:

- `test/duckvep/conformance/structural_hgvs_vep116_differential.R`: 651 exact-span events on GRCh38 chr21 exons; all 651 builder `hgvs_g` equal VEP `--hgvsg --shift_hgvs 0`; all 28,724 DuckVEP transcript HGVS strings for the literal edits equal VEP `--hgvs`; 40 symbolic alleles and one BND pair produce no HGVS in VEP.
- Retained counterexamples: `data/structural_hgvs_vep116_counterexamples.tsv` (VEP's default 3' shift changes the genomic string, so `normalization = 'none'` is not a claim about shifted output) and `data/structural_hgvs_vep116_refused.tsv` (duplications immediately followed by their own repeat are placed by VEP at the later copy; single-base and palindromic inversions are substitutions or no change in VEP, never `inv`).
- Protein HGVS is not part of this domain: the literal route can yield `protein_hgvs`, but it is not validated against structural events here.

Identity and gene evidence: `test/duckvep/conformance/breakend_fusion_corpus_differential.R` over `data/breakend_fusion/` (pins and sha256 in `PROVENANCE.tsv`): FusionCatcher v1.20 test call set (44 physical BND records, 17 published fusions, GRCh38) and a GRIDSS DO52605T excerpt.
