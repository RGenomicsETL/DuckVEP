# Completing compound HGVS and whole-haplotype prediction

**Status: proposal for owner review.** This document recommends implementation choices and acceptance experiments for [#11](https://github.com/RGenomicsETL/DuckVEP/issues/11) and [#50](https://github.com/RGenomicsETL/DuckVEP/issues/50). It does not change the [active prediction contract](duckvep_haplotype_contract.md), authorize new biological assumptions, or close either issue. The implementation baseline is `9290d8b0670b49414a9f8455654f62f7ca2920cc` (PR #54).

## Decision in brief

Extend the existing pipeline in small, independently proved slices: **source evidence → joint edits and ownership → RNA/product qualification → translation → consequences/HGVS/NMD**. Keep compatibility replay as a separate output contract. Reuse the phase, projection, translation and nomenclature kernels; do not build another inference engine or a generic plugin system.

Start with coincident identical general delins and complete diploid single-heterozygote raw evidence. The probes establish an exact-duplicate admission gap and unique unordered diploid sequences for the single-heterozygote case; raw-slot identifiability is a separate unproved acceptance obligation. Then add reached subeffects and reference-backed positional domains, while building the typed RNA/programme inputs needed for the remaining biological domains. This order reduces uncertainty; later rows remain required work, not optional follow-ups.

A single output DNA string is not proof of source attribution, phase, HGVS or a particular splice product. A lack of information is also not an implementation excuse: distinguish a reproducibly underdetermined input from a determinate case the implementation cannot yet handle.

## Original scope and current evidence

#11 requires compound phased DNA HGVS, protein HGVS beyond independent-event coverage, original/shifted operand ownership, independent nomenclature validation and retention of supported protein behavior. #50 requires consequences and NMD for the rebuilt haplotype, covering the unfinished overlap, phase/source and transcript-programme domains. Its historical HG002 counts describe `bbe2ec4`, not the current implementation.

| Obligation | Delivered profile at the baseline | Remaining acceptance work |
|---|---|---|
| Compound HGVS and literal overlaps (#11/#50) | Nonoverlapping compound edits; agreeing substitutions; duplicate insertions/deletions and guarded placed deletions; source-aware blocks, post-stop roles and joint geometry | General unequal-length delins and distinct interacting indel/gap constraints; complete source-sensitive nomenclature inventory |
| Non-strict/raw source prediction (#50) | Complete original diploid homozygous copy proof; omitted non-diploid REF evidence; sample-independent qualification views | Other provably identifiable paths; missing, mixed-ploidy and phase-block contexts with explicit observed/conditional boundaries |
| Noncoding/UTR indels (#50) | Complete noncoding RNA or one UTR, guarded nonoverlapping literal indels/compensation on both strands | Overlapping RNA edits, boundary cases and mixed domains; public whole-transcript ablation admission distinct from a native empty-RNA classifier |
| Intronic, splice, regulatory and mixed domains (#50) | Existing independent-event classification is available; it is not a compound-product oracle | Joint positional facts, explicit product assumptions, feature inputs and dependencies; no invented splice product |
| Recoding and curated translation (#50) | Typed Sec, readthrough, initial-M and amino-acid-substitution programmes for guarded substitutions, with separate conditional peptides | Coordinate-preserving indel/frame interactions, wider programmes and programme-disruption evidence |
| RNA editing (#50) | RNA flags fail closed; no typed transcript RNA prediction programme | Producer contract, raw substrate, typed operations, RNA/CDS mapping, persistence and explicit assumptions |
| Start loss and subeffects (#50) | Annotated start loss and highest-level compound effects; local coding-block facts | Evidence-backed alternative initiation and all justified, reached lower-severity effects |
| Reference completeness (#50) | Supported flagged open/partial ends have explicit rules; inconsistent reference inputs are refused | Fresh model audit, producer repairs where evidence exists, and an owner-approved disposition for genuinely unavailable information |

The issue checkboxes describing all agreeing overlaps or all noncoding indels as refused lag the source. Existing coverage does not finish the wider rows. Sequence-resolved literal structural *nomenclature* remains in #11's inventory; a new symbolic-SV/BND/STR or phased structural-composition engine is not implied. [#28](https://github.com/RGenomicsETL/DuckVEP/issues/28) is closed as not planned.

## Executed design probes

The [study script](../test/duckvep/conformance/completion_options_probe.R) runs through the public R wrapper and independently enumerates small unordered diploid haplotype pairs in base R. It asserts baseline behavior, not the desired future admission. Runtime: Linux, R 4.6.0, Rduckdb 1.5.5; extension SHA256 `3e7b867252adc19eb7b4feee80eb011672426720777ca1012cfbd0d8db21dcd6`, whose native sources match the baseline. Its footer predates the publication commit; source equality, not the footer, connects this artifact to the baseline. Do not use an unrebuilt main-checkout binary as evidence for these sources.

```sh
DUCKVEP_EXTENSION=/absolute/path/to/qualified/duckvep.duckdb_extension \
  Rscript --vanilla test/duckvep/conformance/completion_options_probe.R
```

All three studies use reference CDS `ATGGCTGCTGAAGGTTAA` (`MAAEG*`), a complete single-exon standard-code model and retained original input records.

1. **Phase identifiability:** changes c.4G>A and c.13G>C.

   | Original GTs | Independent unordered pairs | Raw source status | Decoded strict status | Public arrangement API |
   |---|---:|---|---|---|
   | `1/1`, `1/1` | 1 | `predicted` | `predicted` | Rejects homozygous rows |
   | `0/1`, `1/1` | 1 | `unsupported_context` | `predicted` | Rejects homozygous background |
   | `0/1`, `0/1` | 2 | `unsupported_context` | `incomplete_input` | Two explicitly hypothetical arrangements |

   A unique unordered pair does not establish physical lane names. The single-heterozygote case has a potential observed-copy mapping, but that mapping still needs proof against retained raw slots. The public arrangement function is not a drop-in validator for complete genotypes; filtering out homozygous background would change the biological input.

2. **A determinate #11 gap:** one c.4_9 `GCTGCT>AAA` replacement yields `ATGAAAGAAGGTTAA / MKEG*` and predicts. Two distinct IDs asserting that exact same delins yield `unsupported_overlap`. Opaque duplicate-once semantics do not require guessing which ALT base descended from which deleted REF base.

3. **A reached-subeffect gap:** c.4G>A, c.10G>T and c.13G>C yield `MTA*`. Current whole-haplotype terms are only `stop_gained`; roles are `applied/applied/post_stop`. The Ala2Thr difference is before the reached stop at residue 4. The Gly5Arg DNA change is after it and cannot become a translated missense subeffect.

These probes establish specific gaps, not whole-domain correctness. The proposed admission, nomenclature extensions and experiments below are **NOT_RUN** implementations unless explicitly identified as existing baseline tests.

## Options and recommendations

### 1. Overlap: explicit atoms and ownership, not guessed alignments

| Option | Why consider it | Strongest rejection reason / trade-off |
|---|---|---|
| Add pairwise exceptions | Small code changes; existing guarded predicates are useful | Pairwise compatibility alone does not prove a global allele or attribution. An exception matrix grows difficult to explain and exhaustively test. |
| Extend bounded reference-base/gap/opaque-replacement atoms | Fits owned scratch and preserves the distinction between source islands and molecular operations | Moderate implementation and oracle cost; a generic alignment/search step can explode or invent ALT ancestry. Completeness of the admitted interpretation space must be proved. |
| Import an aligner/normalizer/composer | Mature sequence algorithms may be reusable | No evaluated dependency establishes phase, source identity, HGVS, post-stop attribution and charged allocation together. An optimal alignment is not a biological ancestry proof. Adapter/license/maintenance cost remains unevaluated. |

**Recommend the second option incrementally, retaining existing proved predicates.** Begin with coincident identical minimal delins, independently validate every original REF, and attach every original source to one molecular atom. Unchanged padding is context, not an assignment. Use exact immutable-reference span and payload identity, not equality of two final proteins or CDS strings. Keep raw ordered replay separate.

For the first slice, permit only the alias group and disjoint already-supported edits under existing observed-cis and ordinary-CDS guards; refuse distinct edits inside the opaque replacement. Do not infer a base-level map inside its ALT. Joint molecular sequence, joint HGVS and molecular geometry must match the single-source control. Preserve and independently validate each source's original and shifted nomenclature operands; source membership and source-labelled block/provenance materialization may differ. Do not equate source HGVS merely because final DNA or protein agrees. The duplicate-assertion convention requires owner approval; it is not a claim that two records prove two independent biological events. Start/frame/stop-sensitive controls are mandatory before widening that admission; a temporary guard must remain an explicit unfinished case.

Reuse source/edit IDs and joint builders. The current CDS-sized `joint_alt` scratch does **not** establish space for expanding replacements. Core must own or pin payload storage for the complete borrow lifetime and charge any expansion; kernels allocate nothing and perform no I/O. Account for sequence, atom, membership and HGVS capacities before publishing a row. Exact/one-short bounds must fail atomically and allow connection reuse.

For distinct overlapping replacements, first enumerate tiny independent witness spaces, then justify an admissible representation. Distinguish: proved inconsistency; multiple valid alleles; one allele with incompatible attribution/nomenclature interpretations; and an exhausted work budget. Only the first is `edit_conflict`. A cap hit is not evidence of contradiction, ambiguity or uniqueness. The alias slice needs no search; future placement search needs explicit state/comparison/payload limits and a capacity result approved before implementation.

**Decision-changing evidence:** a useful distinct-source class with a complete, bounded ancestry/membership proof supports widening; a counterexample with the same final string but different post-stop or HGVS ownership rejects that widening. Keep external tools as pinned oracles unless a concrete reuse candidate passes this full contract; no new dependency is proposed here.

### 2. Phase: prove original calls, do not promote replay slots

Prefer a shared proof path rooted in the existing phase/carrier code over either unconditional raw admission or another enumeration engine. Capture required original genotype evidence before compatibility slot expansion/folding. The first target is complete diploid input with one unphased heterozygous site and homozygous background, still using disjoint literal coding replay.

Prove both the unordered pair and the mapping of each retained raw output slot to a real copy. Map every retained slot and contributor to complete original calls, including REF-only evidence; preserve raw operands byte-for-byte and do not invent observed lane labels. If that mapping fails, the raw path stays unresolved even when decoded prediction succeeds. Keep global lane-swap equivalence separate from observed phase labels. Missing alleles, unresolved PS relationships, non-diploid normalization and absent original evidence cannot inherit this proof. Genotypes omitted because they push no ALT remain evidence. Extend homogeneous qualification views so unrelated samples cannot change a carrier's eligibility.

Do not merely invoke the current public arrangement API: its heterozygous-only admission omits required background. Reuse or extend the underlying phase machinery only with a complete-call contract and source-preserving projection. Keep hypothetical arrangements hypothetical, even when their consequence strings agree. Initial cost is evidence storage/mapping and oracle work, not a new public schema; broader ploidy contexts need their own proof and retained-original-ploidy checks.

### 3. Consequences: collect proved facts before choosing severity

Retain the current most-severe result, but accumulate additional facts from the *joint reached product and coordinate map*, not from per-event masks. The executed Ala2Thr/stop example is the first target. Partition facts by translation reachability and interactions: same-codon changes are classified jointly; restored frame/sequence can remove component effects; post-stop DNA cannot add translated effects.

Blind OR of `coding_blocks.local_consequence_mask` is rejected: local masks are not proof that a product reaches those blocks under start loss, frame changes, recoding or termination. Reuse the whole-haplotype classifier with explicit preconditions for each added bit. Most-severe term and impact remain derived from the final justified mask. Each changed term needs an independent field-level explanation; unchanged raw/source operands are mandatory.

### 4. Positional domains and product selection are different claims

For intronic/splice/regulatory/mixed work, build reference-backed joint DNA first, using explicit reference resources and model/feature intervals. Positional donor/acceptor/region facts can be determinate without identifying the transcript produced. A changed donor does not decide exon skipping versus intron retention; a motif or upstream ATG does not establish expression or translation initiation.

Compare two useful, separate modes: (a) positional facts with downstream product unknown, and (b) an explicitly supplied, named transcript/translation programme producing a conditional product. Reject heuristic splice-product selection or a default downstream-ATG scan as unconditional prediction. Existing independent-event code can supply tested geometry predicates, not a union of event consequences or a new product authority.

**Output decision required:** the 36-field result has unconditional consequence fields and conditional peptide operands, not a general parallel conditional consequence/NMD relation. Do not overload `predicted` or silently mix unconditional positional terms with conditional coding terms. Before mixed-domain implementation, approve a field-by-field table of unconditional facts, conditional peptide operands, unavailable conditional consequence/NMD claims and retained raw operands. Demonstrate two supplied products for one DNA allele without promoting either to observed biology. A status/reason cannot restore facts the fields cannot represent: if the 36-field shape cannot express the required distinctions, the slice remains incomplete pending an owner-approved versioned representation decision. That representation gap is unfinished #50 work, not grounds to relabel the domain complete. No unversioned sidecar or hidden SQL object is proposed.

Regulatory composition additionally needs declared feature intervals and model provenance; transcript/exon geometry alone is not a regulatory model. Downstream product and NMD remain unknown until the relevant RNA, initiation, termination and exon-junction geometry are justified. EJC50 is a declared geometric prediction, not proof of caused decay.

### 5. RNA/curated programmes: fix the input contract before broadening inference

The pinned Ensembl 116 API assembles exon RNA, applies transcript `_rna_edit` operations in cDNA coordinates, then derives CDS with edit-dependent boundary shifts and possible `_transl_start`/`_transl_end` overrides. SeqEdit's `start end alt` format does not carry an expected REF string. The inspected human core has zero transcript `_rna_edit` rows and five translation-attached values; those translation values must not be reinterpreted as transcript SeqEdits. Zero rows in one species/release does not remove the requirement.

Relevant upstream authorities are [Transcript.pm](https://github.com/Ensembl/ensembl/blob/release/116/modules/Bio/EnsEMBL/Transcript.pm), [SeqEdit.pm](https://github.com/Ensembl/ensembl/blob/release/116/modules/Bio/EnsEMBL/SeqEdit.pm), and [Translation.pm](https://github.com/Ensembl/ensembl/blob/release/116/modules/Bio/EnsEMBL/Translation.pm). Pin actual package/source/model hashes in execution receipts; the links alone are not binary verification.

**Recommend typed, substrate-checked intervention records and separate prediction preparation.** Capture actual unedited genomic/cDNA substrate and its provenance before applying programmes. Compatibility operands are not automatically unedited biological truth: preserve them byte-for-byte, rather than rebuilding them in place. RNA edits are not amino-acid substitutions.

A producer/model contract must identify operation kind, transcript/cDNA coordinate system, start/end convention including insertions, expected substrate derived from the declared raw sequence, ALT payload, programme identity/evidence, and CDS-boundary rules. Distinguish omitted/unknown metadata from a known empty programme. First capture typed transcript SeqEdits, raw genome-spliced substrate, boundary overrides and provenance while retaining the prediction veto. Specify original cDNA coordinates, insertion conventions, edit ordering/ties, variant-to-source mapping and boundary transformation before replay admission. The pinned API applies transcript edits highest-start first, then derives edited CDS boundaries; that does not establish variant ordering or editing persistence under an allele. Derive expected REF from the pinned substrate; never claim SeqEdit supplied it. At a changed substrate, refuse or apply an explicitly supported allele-conditioned rule—never blindly overwrite it with the annotated ALT.

The proposed first prediction slice admits only complete, equal-length programmes and genomic substitutions with a fixed RNA/CDS map. Reject overlaps, ambiguous adjacent insertions, shared edited codons, upstream length changes, boundary overrides and genomic ALT-equals-edit cases; each relaxation needs its own witness and approved assumption. Validate a synthetic or nonhuman programme against pinned Perl transcript-edit execution. Even then, persistence under the allele is a named assumption, not proved by matching edited bases; conditional peptides do not imply unconditional SO or NMD. Human-116 absence cannot supply this positive oracle.

Store typed records in owned, fingerprinted model data. Introduce an explicit snapshot wire version when the representation requires it; authenticate older layouts before mapping, retain legacy fingerprints when optional data is absent, and test mmap/copy readers and old/new producer combinations. INSTALL/LOAD must not fetch genomes or construct hidden tables.

Paired prediction peptides can use existing conditional fields for the bounded programme; unconditional SO and NMD cannot be filled by implication. RNA edits that change length or span junctions require justified product geometry, not merely a copied exon table. Wider peptide programmes need a map of surviving original codons through joint edits; an unchanged ALT substring or shifted peptide index is not sufficient to inherit Sec/readthrough/initial-M/AA edits.

This has higher producer/model/persistence cost than removing a guard, but makes assumptions inspectable and reversible. Preserve refusals where programme identity, substrate or mapping is missing. A new metadata contract and migration fixture are approval gates, not implementation details to hide in a patch.

### 6. Alternative initiation and problematic references

Annotated-start disruption is already a supported fact. Predicting rescue requires a declared initiation context/programme and an explicit counterfactual: whether the reference retains its annotated start, what changes trigger an alternative, and which mapped RNA site is used for the alternate. Enumerated candidate ORFs may be conditional possibilities, not observed translation. A genome-wide ATG scan, Kozak score or per-event VEP result cannot select the true product by itself.

For truncated/internal-stop reference residuals, first rerun a named-model audit with source hashes and original completeness/curation flags. Classify producer errors, supported flagged incompleteness, genuine biological programmes and unavailable information separately. Repair producer defects against evidence; never trim CDS, insert a start, erase a stop or synthesize missing metadata to obtain `predicted`. The historical 79/76 transcript counts are not current measurements. Closure needs a fresh disposition of these original residuals and owner agreement for information that truly cannot be recovered.

## Implementation owners

| Concern | Existing owner to extend |
|---|---|
| Molecular composition, opaque aliases, post-stop/geometry facts | [`duckvep_haplotype_stream.c`](../src/kernel/src/duckvep_haplotype_stream.c): `strict_literal_union`, `island_relation`, `finish_prediction`; existing joint buffers/builders |
| Original-call evidence and qualification views | Stream `duckvep_haplotype_stream_push_raw_call` admission; [`duckvep_carriers.c`](../src/kernel/src/duckvep_carriers.c): `share_copy_unproven`, `partition_copy_proof`; the existing [phase kernel](../src/kernel/src/duckvep_phase.c) |
| Whole-product classification and conditional peptide preparation | `classify_haplotype`, `noncoding_haplotype_mask`, `conditional_prediction` in the stream; [haplotype kernel](../src/kernel/src/duckvep_haplotype.c) and [effect predicates](../src/kernel/src/duckvep_effect.c) |
| Charged workspace, output types and source/block materialization | [`duckvep_core_haplotypes.c`](../src/core/duckvep_core_haplotypes.c), then both host adapters; no kernel allocation or schema drift |
| Programme extraction, validation and persistence | [`duckvep_core_ensembl.c`](../src/core/duckvep_core_ensembl.c), [`duckvep_core_model.c`](../src/core/duckvep_core_model.c), [`duckvep_core_snapshot.c`](../src/core/duckvep_core_snapshot.c) |
| Compound nomenclature | Existing DNA/protein builders in [`duckvep_hgvs.c`](../src/kernel/src/duckvep_hgvs.c); consume proved joint edits and separate source ownership |

## Assumptions to approve and ways to falsify them

| Proposed assumption | Evidence required / falsifier |
|---|---|
| Identical reference-relative delins assert one molecular change | Every original REF validates; minimal atoms and source membership agree. Different padding/source count must not change geometry or nomenclature. A distinct interior interaction invalidates the alias-only proof. |
| A complete diploid single-heterozygote path is identifiable | Independent unordered-pair enumeration plus a real-copy mapping. A missing/ref-only/PS/ploidy counterexample or dependence on another sample rejects admission. |
| A subeffect is part of the translated joint outcome | The reached product and frame map establish it. Cancellation, alternate initiation or a prior stop falsifies a naive event-mask union. |
| Positional splice/regulatory facts do not identify a product | Distinct products consistent with the same DNA are negative controls. Only additional declared evidence/programmes resolve them. |
| A named programme persists under the allele | Programme identity, substrate and mapping preconditions hold. Site disruption or unknown dependency keeps the result conditional/unknown, not silently normal translation. |
| A reference can be repaired | A pinned source/producer witness proves the replacement data. An attractive translated peptide is not repair evidence. |

Assumptions about ambiguous source conventions, product selection or new output semantics require owner approval. None is inferred from a green host-equality test.

## Experiments and closure gates

The executed studies above are the seed cases. These extensions are **NOT_RUN** and become required acceptance work for the corresponding slice:

- **Opaque delins:** duplicate `GCT>TTTCCC` must match the single-source `ATGTTTCCCGAATAA / MFPE*` control on reference `ATGGCTGAATAA`; compare a second source c.5C>T inside the opaque replacement, whose ancestry is not supplied. Test positive/negative strands, tied-source order, redundant sources, padding, expanding payloads, frame/stop boundaries and separate source/joint HGVS expectations.
- **Placement/gap ambiguity:** repeat-slidable deletion plus a substitution; two unequal insertions in one gap; an insertion whose attachment is deleted. Enumerate consistent bounded interpretations independently. Never choose an order by event ID, source order or smallest HGVS string.
- **Phase:** one versus two unphased heterozygotes with homozygous background; missing and omitted non-diploid REF calls before/after ALT; split PS and sample addition/removal; map every retained slot to observed or hypothetical evidence. Limits must not silently drop an arrangement or call.
- **Subeffects:** reached upstream missense plus stop gain, downstream-only changes, same-codon cancellation and frame restoration. Verify the final mask, severity, impact, source roles and unchanged compatibility peptide separately.
- **Splice/mixed:** the same donor-damaging allele is compatible with at least two transcript products; without a supplied programme, assert no selected protein/NMD. With each explicitly supplied product, independently reconstruct RNA and mapping and keep conditional results labelled. Include pure intronic, UTR/CDS and regulatory-feature controls.
- **RNA/recoding:** a substrate-preserving edit, the same site disrupted, insertion/deletion shifting CDS boundaries, boundary overrides and omitted versus empty metadata; an upstream indel preserving a curated codon's ancestry versus a frameshift destroying it. Cross-junction RNA changes must not fabricate EJC geometry.
- **Initiation/completeness:** two viable downstream ORFs without selection evidence; one declared conditional programme; a correctly flagged partial CDS versus an unflagged truncated counterpart. No unqualified rescue or model rewriting.

For each row, establish separate authorities for allele composition (base-R/hand derivation), phase/source identity (original calls plus independent enumeration), nomenclature (HGVS standards and pinned independent observations), compound SO (the declared csq-compatible contract where applicable), and edited stop/junction geometry. Oracle disagreement must be adjudicated at the disputed fields. Nomenclature is not cis evidence; independent-event VEP is not a compound-product oracle.

Before merging an implementation slice, require registered native properties, public SQL and installed-R checks, actual sanitizer flags, boundary/fault-injection and cleanup tests, sample/order/strand metamorphic checks, full v1/v2 rows, unchanged decoded/raw regression operands and individually explained golden deltas. Use the existing HG002 and pinned nomenclature gates as regressions, not as proof of unrepresented overlap or programme classes. Each slice must state worst-case work and charged storage in source count, sequence length, payload bytes and programme count, including sorting/comparisons, memberships, maps, retained evidence and simultaneous buffers. Include model-owned storage and borrow lifetimes, not only query scratch. Test overflow as well as exact/one-short capacities. Alias admission avoids placement search; future search exhaustion needs a distinct approved result, never `edit_conflict` or biological ambiguity. Measure comparison counts, memory charges and representative complete-query cost before selecting wider search defaults.

**#11 closes only** when the original compound DNA/protein/ownership inventory has positive and negative witnesses with independent source-specific nomenclature evidence; determinate general cases cannot be deferred under an ambiguity label. **#50 closes only** when every original domain has implemented, tested semantics and every residual reference/context has an evidenced disposition. Scientifically unidentifiable cases require an explicit accepted uncertainty policy; refusal alone does not complete a promised domain. Release, community submission and the released DuckDB/R 2.0 gate remain separate decisions.

## Proposed delivery and review boundaries

1. **Alias and copy-proof slices:** exact duplicate general delins; complete diploid single-heterozygote raw paths. Separate commits/receipts, with no producer or schema change. Keep each reversible by restoring its admission guard.
2. **Reached effects and positional profile:** add joint subeffects, then approve the positional/conditional output contract and implement reference-backed intronic/splice/regulatory/mixed facts. No inferred product.
3. **Typed programme substrate:** producer capture, model validation/ownership, fingerprints/snapshots, RNA/CDS maps and migration tests; then conditional RNA and wider peptide programmes.
4. **Product-dependent domains and full inventory:** explicit splice/translation scenarios, alternative initiation, remaining distinct overlaps/nomenclature and reference residuals. Their witness design can proceed in parallel; none is removed from the finish line.
5. **Acceptance audit:** rerun the original class matrix and current-model residual inventory, reconcile issue checkboxes with exact receipts, and seek closure approval. Do not use issue wording changes as implementation evidence.

The owner's immediate decisions are the alias/source interpretation, the strict original-call proof boundary, the positional versus conditional output representation, and the typed programme contract. Approving this plan authorizes neither arbitrary biological defaults nor automatic issue closure.
