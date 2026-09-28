# Haplotype prediction contract: `duckvep-coding-v1`

Status: **signed by the maintainer on 2026-09-28** for RGenomicsETL/DuckVEP#2, with two amendments. The NMD rule is `ejc50-v1` as written. The scale targets in section 4 are **replaced by tighter targets that the maintainer will set**; until then, section 4's figures are placeholders, not gates. Deferred work is tracked in #11 (compound HGVS), #12 (extended domains and Haplosaurus grouping), #13 (extended consequences and NMD models) and #4 item 5 (phased structural composition).

**Scope:** close #2 on a versioned, coding-only vertical, not “general haplotype prediction.” Maintainer approval is required before classifiers change; #8 and #3 land first. Keep independent-event and haplotype classifiers separate, sharing projection/edit/translation facts only where their semantics agree. A framework rewrite is unnecessary. Neither a `bcftools csq` replacement nor a runtime hybrid establishes the chosen output authority; extend the existing native path.

## 1. Authority per output

**Haplosaurus 116:** exact reconstructed CDS/protein bytes (including terminal-stop conventions), carrier sample/count multisets, `has_indel`, CDS flags (`indel`, `frameshift`, `resolved_frameshift`) and protein flags (`indel`, `stop_change`). The documentation’s `stop_changed` spelling differs from the pinned API. Lane assignments/provenance require input-genotype truth tables: grouped public output aggregates carriers. Pin executable/API commits, reference/model, options and input representation. Compare keyed fields, not JSON ordering or generated group IDs. Never replace nominal indel/frame flags with alignment-derived guesses. Compare on unambiguous phased inputs; strict DuckVEP phase semantics are not Haplosaurus’s PS-ignoring file semantics.

Sequence-group flags can inherit the first encountered member’s edit history. Heterogeneous grouped flags therefore remain explicitly outside this vertical, with keyed disagreements retained—not accepted mismatches. Path-level flag comparisons use separately executed, unambiguous Haplosaurus fixtures. Optional pathogenicity flags require their own enabled data/resources; they are not inferred here.

**DuckVEP policy, not VEP compatibility:** Haplosaurus supplies no whole-haplotype SO set, IMPACT or NMD prediction. Publish `duckvep-coding-v1` separately. Its independent oracle must reconstruct complete edited transcripts and translate both CDSs under one stop convention with standalone R, using reviewed goldens—not DuckVEP helpers, local block masks or independent-event label unions. Haplosaurus sequence equality alone cannot certify these outputs.

## 2. First supported vertical

Support strict, complete diploid calls; literal A/C/G/T SNPs/MNVs/indels with normalized REF/ALT lengths ≤50 bases; nonoverlapping edits contained within coding exons; complete, phase-zero-start, table-1 CDSs with canonical start/terminal stop and no curated RNA/peptide edits or recoding. Support both strands, multiple coding exons and their internal phases, reference lanes, multiple transcripts and small multisample fixtures. Splicing is fixed to the selected model. Each transcript/sample must have one unambiguous heterozygous phase domain. Phased calls without PS use the documented implicit domain, not inferred cross-block phase.

Every result retains existing sequences, carriers, contributors and statuses, plus versioned `haplotype_consequences`, `haplotype_impact`, `nmd_prediction`, `prediction_status` and reason. Keys include model/transcript identity and sample/phase/lane. Preserve source-record identity, ALT ordinal, original operands, genotype evidence and normalized edit identity through every block; sequence deduplication must not erase contributors.

The proposed **reduced SO set** reports reached HIGH start/stop/frame effects; otherwise it reports one whole-protein category below. This deliberately omits lower-severity subeffects. Start loss suppresses other biological predictions; otherwise frame and stop terms may coexist. Among LOW classes, terminal-retained precedes synonymous. IMPACT is the maximum severity of emitted terms. Frame SO is normalized-edit-path-sensitive: pin/test that convention rather than claiming representation-independent labels.

| Situation | Exact v1 policy output | Oracle |
|---|---|---|
| Same-codon cis substitutions; corresponding trans control | Separate trans lanes. With no start/stop/frame effect: identical peptide → `synonymous_variant`/LOW; changed peptide → `missense_variant`/MODERATE. Classify the combined sequence, not either SNP. | Haplosaurus sequences/carriers/flags; independent R policy goldens. |
| Open/restored frames | Premature first-stop codon intersects a displaced-frame interval, or frame remains displaced at CDS exhaustion → `frameshift_variant`/HIGH. Restored **before termination** → no frameshift SO; identical peptide → synonymous/LOW, otherwise `protein_altering_variant`/MODERATE. Pure frame-preserving insertion/deletion uses `inframe_insertion`/`inframe_deletion`/MODERATE. | Same two oracles; track cumulative frame offsets and first stop. |
| Stop created/removed | New first stop before the homologous reference terminator → `stop_gained`/HIGH. Reference termination abolished without an earlier gained stop → `stop_lost`/HIGH. Downstream edits remain contributors, not expressed effects. | Same two oracles; include compensating edits both before and after the first stop. |
| Start/terminal codons | Canonical start abolished → `start_lost`/HIGH; initiation and NMD unknown. Changed terminal codon remaining a stop, without a higher effect → `stop_retained_variant`/LOW. No invented downstream extension sequence. | Same two oracles; terminal synonym, loss and restoration controls. |

Reference-only lanes have an empty SO set and NULL IMPACT. Other mixed frame-preserving replacements use `protein_altering_variant`; unchanged peptide uses synonymous.

**NMD:** propose `ejc50-v1`: for a newly premature stop, return `trigger` when J−S>50, otherwise `escape`. S is the stop’s final nucleotide; J the penultimate exon’s final nucleotide, both 1-based in edited spliced-transcript coordinates. Intronless transcripts escape. Known termination without a new premature stop → `not_applicable`; incomplete phase/exon topology, lost initiation or unavailable termination → `unknown`. Return stop/junction coordinates and joint contributor attribution. Test 49/50/51 bases, intronless/last-exon cases and indel-shifted junctions independently. This omits reinitiation/long-exon exceptions: it is an EJC-distance heuristic, not clinical NMD truth, VEP’s allele-position plugin, or the `NMD_transcript_variant` biotype term.

**Failure contract:** missing alleles, unphased heterozygosity or unresolved cross-PS phase → `incomplete_input`; contradictory edits → `edit_conflict`; other overlaps/ambiguous same-gap insertions → `unsupported_overlap`. Preserve reference-mismatch/projection reasons; excluded contexts → `unsupported_context`. Whole SO/IMPACT are NULL and NMD unknown on these paths; any compatibility replay stays explicitly conditional. Retain every contributor, including omitted, shadowed and post-stop sources. Malformed identities/budget overflow fail explicitly rather than truncate.

## 3. Named follow-ups, not closure blockers

- **#11, “Compound haplotype HGVS and overlap nomenclature”**: move checklist item 3 and its compound shifted-HGVS work from item 6 here. Sequence equality cannot validate nomenclature. Existing supported protein HGVS remains a regression gate.
- **#4 item 5 — “Phased structural composition”**: move item 4 there; typed SV/BND/STR composition requires its own event/topology contract.
- **#12, “Extended haplotype domains and Haplosaurus grouping”**: defer item 5’s arbitrary ploidy, uncertain-phase enumeration, general overlaps/raw-parser contexts, nonstandard/partial transcripts, larger alleles and heterogeneous grouped flags. Retain current replay tests and disagreement ledgers.
- **#13, “Extended haplotype consequences and translation”**: defer splice-changing/regulatory/UTR composition, alternative initiation, downstream stop-loss extension, richer SO sets and additional NMD models.

Items 1–2, supported-domain validation, original-operand ownership and item 6’s reuse/scale measurements remain required.

## 4. Scale contract

Qualify a **5M-physical-variant, single-sample GRCh38/MANE-selected-model** job, not 5M expanded call rows. Account for every input as emitted, explicitly outside coding scope, or unavailable. The existing 2,621,440-row fixture contains only 2,560 physical events; its favorable sharing is not genome-scale evidence.

Enforce **16 GiB per-job process memory, including R**, and **4 GiB aggregate DuckVEP-native memory**, including resident models, indexes and all workers—not merely `workspace_limit`. Use native admission/accounting, bounded windows, spillable DuckDB staging, a temporary-disk quota and an external job limit. Table-backed/chunked output is the scale interface; eager R collection is not the throughput benchmark. Exercise capacity failure, cleanup and connection reuse.

Acceptance targets: **to be set by the maintainer, tighter than** the originally proposed ≤15 minutes sorted execution (~5,556 source variants/s) and ≤30 minutes sort-plus-execution (~2,778/s), on a recorded deployment-class host with ≤6 allocated cores and local SSD. Do not import #3’s independent-event million-ALT/s floor.

Measure separately: (A) preordered input through prediction and complete output materialization; (B) identical unsorted input including decoding, discovery, staging and sorting. Include unavoidable internal sorts in A. Report cold model-load and warm execution separately; all phases obey memory caps. Record physical sources/ALTs, projections, calls, carrier states, unique paths, translated bases, output rows/bytes, peak active window, model/native/DuckDB/RSS peaks, spill bytes and full-output checksums. Use real phased input plus dense/long-transcript and low-sharing stress controls. Run three fresh processes per mode: median-time gates, caps on every run, complete failures reported.

## 5. Builder-sized slices and gates

Each slice is one focused builder run; split again if it changes more than its named mechanism. Every code slice requires native properties/sanitizers, SQL and installed-R tests, unchanged independent-event regressions, relevant hot-path measurements, CI and ready-head review.

1. **Contract/oracle fixtures:** maintainer approves field authority, reduced SO rules, phase failures, NMD heuristic and numerical budgets. Gate: pinned goldens and corruption controls; no classifiers yet.
2. **Eligibility/provenance:** add prediction statuses and phase-domain completeness checks. Gate: cis/trans, missing-call, overlap/conflict, reference and contributor-conservation tests, including vector boundaries.
3. **Same-codon classifier:** gate exact Haplosaurus replay fields and independent whole-policy expectations on both strands/exon layouts.
4. **Frame/restoration classifier:** gate open/restored intervals and early-stop-before-restoration counterexamples; preserve every contributor.
5. **Start/stop classifier:** gate created, abolished, retained and rescued start/termination combinations, including post-stop edits.
6. **NMD attribution:** gate edited-transcript junction geometry, threshold boundaries and unknown/not-applicable distinctions; no per-allele NMD union.
7. **Scale qualification:** stage 5M-source workload, then measure both execution modes. Gate caps, signed throughput targets, complete fingerprints and explicit failures at the closure revision.

## 6. Proposed rewritten “Done when”

> The maintainer has signed `duckvep-coding-v1`, its Haplosaurus field authority, EJC50 heuristic, supported/excluded domains and resource/throughput gates. Slices 1–7 pass at the closure revision with zero unexplained missing, extra or discordant supported outputs; expected unsupported cases retain statuses and all contributor identities. Independent-event and existing replay contracts remain regression-clean. Public SQL/R documentation states these limits, and deferred checklist work is linked to the named follow-ups. Passing these gates closes #2; it does not assert general haplotype, compound-HGVS or structural-composition compatibility.
