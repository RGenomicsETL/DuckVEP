# DuckVEP compatibility and errata

Status: current compatibility guidance for pinned Ensembl VEP 116. DuckVEP is alpha;
this document is not a claim of complete conformance or clinical validation.

## What the compatibility target means

Within its supported surface, DuckVEP must reproduce the pinned VEP 116 executable
under the same reference, transcript model and settings, including absent HGVS values
and input-representation-dependent results. The [compatibility contract](design/duckvep.md)
pins the VEP, Ensembl core and variation revisions. The
[upstream source registry](test/duckvep/upstream/sources.tsv),
[dependency lock](test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt)
and individual comparison receipts supply source anchors, dependencies and run identities.

Every physical input record and source ALT remains a separate comparison unit, including
duplicate POS/REF/ALT records. Equal reconstructed sequence does not permit replacing their
different VEP outputs with a canonical answer. A disagreement stays in the denominator.

Three verdicts must remain separate:

- **Observed upstream convention:** reproducible VEP behavior that DuckVEP must follow.
- **DuckVEP gap:** unsupported behavior or a retained compatibility disagreement.
- **Potential upstream erratum:** a result contradicting an independently established
  contract. Reproduction, interpretation and upstream acknowledgement are separate facts.

A potential upstream erratum requires a minimized input, pinned model and executable,
complete outputs and diagnostics, an independent sequence or versioned-rule argument,
and controls excluding our parser, projection and comparator. Agreement with another
tool is insufficient. No upstream acknowledgement is recorded for the observations
below. None authorizes silently changing compatibility output.

[HGVS 21.1.4 recommendations](https://hgvs-nomenclature.org/21.1.4/) are a separate
nomenclature audit reference. Matching VEP does not certify HGVS correctness; an
`ok` result is a computation status, not that certificate. The internal strict
compatibility control is neither a second public consequence standard nor a certified
HGVS implementation. Genotype `phase_policy := 'strict'` addresses a different question.

## Remaining limitations

Haplotype sequence replay is a foundation, not complete compound consequence prediction.
Whole-haplotype SO, IMPACT and NMD require their own contract: Haplosaurus reports
sequences, differences, carriers and flags, while the pinned NMD plugin evaluates one
transcript-variation allele. Unioning independent labels or applying the plugin to a
rebuilt stop position does not reproduce a whole-haplotype upstream method.
See the [haplotype follow-up](https://github.com/RGenomicsETL/duckhts/issues/92).

Compound HGVS also has unresolved presentation comparisons. In the
[internal-codon differential](test/duckvep/conformance/hgvs_cis_codon_differential.R),
two cis SNVs changing `AGA` to `CGC` reconstruct the same protein as one MNV.
The compound path emits `p.(=)`; VEP's original MNV emits `p.Arg2=`.
Whole-protein equality and equality at a named residue are different assertions.
The 320 disagreements among 26,352 cis-SNV/MNV comparisons remain failures, not evidence
that VEP is wrong or permission to merge source identities.
These counts come from an unsigned, build-unbound local diagnostic outside release
conformance history. Its untracked receipt is
`test/duckvep/conformance/results/hgvs_cis_codon_3e62f43d481358/receipt.json`;
the linked R file is the comparison driver, not a published result.

The [retained indel comparison](test/duckvep/conformance/data/ambiguous_indel_translation/summary.csv)
contains 168,000 original records and 1,344,000 HGVSp comparisons. Independent and
decoded-singleton routes agree, but raw-source routes retain 41,712 disagreements
with independent-event VEP. Their distinct input semantics are explained below.
All 336,000 independent SO comparisons agree in that finite matrix. These are recorded
diagnostics, not current-build certificates, population error rates, or proof for
compound events, reverse strands or phase-padded indels.

The [conformance guide](test/duckvep/conformance/README.md) and
[rendered report](benchmarks/duckvep_conformance.md) describe executable checks and their
scope. The [evidence policy](design/duckvep_corpus_workflow.md) governs retained failures,
comparison keys and controls; passing native properties cannot replace an executable
differential.

<a id="n-containing-codons-in-independent-and-singleton-protein-annotation"></a>

## Ambiguous bases: translation and allele validity are different

**Classification: corrected DuckVEP compatibility defects, with scoped evidence.**
An unambiguous source allele can lie in an N-containing codon whose amino acid is
determinate. For table-1 CDS `ATGGCNTAA`, CDS position 4 `G>A` changes Ala to Thr:
VEP reports `missense_variant` and `p.Ala2Thr`. Conversely, an N in the changed
source allele can make its peptide unavailable. N is never a wildcard for REF validation.

DuckVEP uses consensus translation while preserving separate checks for uploaded
alleles, raw sequence ambiguity and curated reference proteins. Unknown local residues
can coexist with missense or frame-change facts; they do not justify discarding all
sequence predicates. The [SNV comparison](test/duckvep/conformance/data/ambiguous_codon_consensus/summary.csv)
retains all 28,800 original SNVs over 125 ACGTN codons and 24 tables: 230,400 HGVSp
and 57,600 SO comparisons agree. The
[failing baseline](test/duckvep/conformance/data/ambiguous_codon_baseline) remains available.
These forward, phase-zero internal-codon tests do not establish all transcript contexts.

For length-changing records, a literally matching N in an erased anchor or shared
prefix/suffix differs from N in the changed payload. For CDS `ATGGCNGCCTAA`,
position 6 `N>NGCC` yields independent-event `inframe_insertion` and `p.Ala2dup`.
A removed N can still support a deletion while leaving the reference peptide unavailable.
Raw predicate flags and emitted SO are also distinct: a length-decreasing event may
have true raw missense and frameshift predicates but emit only `frameshift_variant`.

The [indel witnesses](test/duckvep/conformance/data/indel_translation_witnesses.jsonl.gz)
and [predicate witnesses](test/duckvep/conformance/data/indel_predicate_witnesses)
retain original records, actual CLI/direct-API observations and controls. The
[diagnostic driver](test/duckvep/conformance/ambiguous_codon_differential.R) preserves
source identities and absent outputs; direct-API evidence does not automatically
certify VCF parser behavior.

## Record geometry and consequence terms

**Classification: observed VEP-116 conventions.** The uploaded feature, minimized
physical edit, displayed coordinates and altered sequence serve different purposes.
DuckVEP keeps them distinct; normalizing an input before annotation can change the
question VEP answers.

| Topic | Observed behavior and consequence for users |
| --- | --- |
| Retained REF bases | Equal-length uploaded spans determine local peptide windows. They can change start/stop terms even when the differing base is identical. Complete uploaded REF still requires validation. |
| CDS phase | Displayed CDS/protein coordinates use the first transcript exon phase; stored CDS padding uses the first coding exon phase. These may differ when CDS begins in a later exon. |
| UTR and mapper gaps | A feature crossing the 3′ CDS end may lose peptide annotation, while one crossing the 5′ CDS start can retain start predicates. Empty annotated UTR intervals can still produce UTR terms for spanning features. |
| Partial codons and stops | Partial-codon status depends on sequence length and first affected peptide position, not only an attribute. Terminal coordinate tests, local peptides and raw-CDS fallback translation can disagree, including after reference peptide edits. |
| Predicate combinations | Start-lost/start-retained and stop-retained/protein-altering can coexist. In-frame insertion is not determined by length modulo three; deletion and insertion use different predicates. |
| Splice and noncoding features | ALT-only differing runs can reach an intron selected by VEP's expanded cache. Short-intron exon stretching is candidate selection, not exon membership. Mature-miRNA overlap replaces generic noncoding exon terms. |
| Empty consequence sets | A real transcript overlap with no successful predicate receives transcript-associated `intergenic_variant`. It must not be confused with absence of a transcript from an incomplete model. |
| NMD | The plugin projects the full uploaded feature, including reversed insertion intervals, rather than the minimized edit. It is distinct from a transcript's curated `NMD_transcript_variant` biotype term. |

Executable witnesses are in the [projection differential](test/duckvep/conformance/projection_differential.R)
and [corpus differential](test/duckvep/conformance/corpus_differential.R).
The [annotation](test/duckvep/property/duckvep_prop_annotation.c),
[coding](test/duckvep/property/duckvep_prop_coding.c) and
[classification](test/duckvep/property/duckvep_prop_classify.c) properties retain the
specific strand, phase, start/stop, UTR, splice and NMD counterexamples. Pinned
[VariationEffect](https://github.com/Ensembl/ensembl-variation/blob/2fb834b987ede3824e200197a838ce11e91aeb4b/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm)
and [NMD](https://github.com/Ensembl/VEP_plugins/blob/0082591268417af618e03850c5ffdc7c09998a5d/NMD.pm)
remain the source authorities, not biological simplifications of their term names.

## Model inputs and structural events

**Classification: upstream conventions and explicit DuckVEP scope.**

- **Assembly paths:** X/Y PAR, patches and alternate haplotypes retain path-specific
  coordinates and transcript identities. DuckVEP models the exact paths supplied by
  the reference relation; it does not implicitly project absent paths or merge
  equivalent genes. The [PAR fixture](test/duckvep/conformance/data/par_path_witnesses.vcf)
  tests sequence-dependent annotation on both X and Y.
- **Release annotations:** Ensembl release `VE` contains stored variation-effect rows;
  release `CSQ` can overwrite terms for the same allele/feature. Neither replaces
  executable VEP as the oracle. The [release audit](test/duckvep/conformance/release_vcf_differential.R)
  preserves this distinction. Regulatory model preparation excludes EMAR because
  VEP excludes those source rows before overlap evaluation.
- **Breakends:** VEP adds one to local POS but retains the mate coordinate. Ordinary
  transcript terms use the local feature; truncation and regulatory/motif overlap
  can observe the mate. Candidate discovery and the fixed 5,000-base attachment
  distance are separate rules. DuckVEP emits each event/feature once, unioning endpoint
  terms; raw VEP may emit multiple rows. It does not infer fusion sequence. The
  [BND fixture](test/duckvep/conformance/data/breakend_default_witnesses.vcf)
  retains endpoint/default-term controls.
- **BND batching:** VEP's
  [input-buffer interval tree](https://github.com/Ensembl/ensembl-vep/blob/57ea5c52340acc1f156267f810ad162e26597082/modules/Bio/EnsEMBL/VEP/InputBuffer.pm#L345)
  can mix mate coordinates without
  a chromosome key, making a record's transcript set depend on other buffered records.
  BND differentials use `buffer_size=1` to state an isolated-event oracle contract;
  this is not conformance to arbitrary batched output. The
  [multichromosome BND report](benchmarks/duckvep_conformance.md#paired-breakend-differential)
  records that scope.
- **Symbolic and repeat alleles:** VEP accepts a finite symbolic vocabulary, not arbitrary
  `<...>` geometry. A bounded tandem repeat expanded to literal sequence is a small
  variant; an unexpanded repeat uses structural gain/insertion predicates. Nominal SV
  coordinates drive consequences; confidence intervals and inserted sequence remain
  provenance, not inferred exact geometry or compound HGVS. The
  [confidence fixture](test/duckvep/conformance/data/structural_confidence_grch38.vcf)
  and [nominal-coordinate comparison](benchmarks/duckvep_conformance.md#declared-conformance-closure)
  test uncertainty metadata; VEP's
  [structural insertion predicate](https://github.com/Ensembl/ensembl-variation/blob/2fb834b987ede3824e200197a838ce11e91aeb4b/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm#L1100)
  does not inspect inserted sequence.
- **gVCF alleles:** `<*>`, `*`, `<NON_REF>` and `.` are different. The catch-all
  `<*>` has no known alternate coding sequence or HGVS, yet VEP's length predicate can
  produce ablation when a long REF contains a complete feature. The
  [gVCF fixture](test/duckvep/conformance/data/gvcf_semantics.vcf) preserves mixed ALT
  order and long-REF controls. Literal deletions can also ablate transcripts without
  becoming symbolic structural records.

These rules do not establish a general structural-variant or pangenome annotation
engine. The [structural kernel](src/duckvep/kernel/src/duckvep_sv.c) and
[classification tests](test/duckvep/property/duckvep_prop_classify.c) define the typed
subset; original ALT, confidence, orientation and source identity remain necessary.

## HGVS follows VEP's state machine

**Classification: observed executable conventions, not a certified HGVS standard.**
Protein HGVS is not simply a difference between two complete proteins, and transcript
HGVS is not simply a walk along spliced CDS.

- Transcript HGVS shifts against a genomic reference window in transcript direction.
  VEP's clipped 1,000-base flank construction, allele-length limits and cached pre-shift
  transcript slice affect endpoint cases. Duplication tests the copied source before
  ordinary insertion flanks. These operations must not rewrite source coordinates or
  consequence events; `shift_hgvs` is not `shift_3prime` or `shift_genomic`.
- Literal exonic SNPs have a phase-aware HGVS coordinate path that MNVs and indels do
  not share. A transcript-flank consequence can legitimately have no transcript HGVS;
  missing reference or failed REF validation is instead unresolved.
- Protein notation uses cached unshifted start/stop predicates and independently
  reconstructed frameshift state. Shared peptide clipping, terminal insertion flanks
  and stop-loss precedence affect both strings and applicability.
- Some HGVS paths translate with BioPerl's default table 1: frameshift alternate CDS,
  late stop search and raw-reference duplication checks. Ordinary consequences still
  use the transcript's declared table. Mitochondrial HGVSp can therefore use a stop
  incompatible with its consequence peptide's translation table.
- Other executable details include incomplete-codon assignment behavior, distinct
  terminal-partial insertion views, negative substring positions, Xaa-to-Ter conversion
  and a delins extension guard after three-letter conversion. A rendered `Ter` is
  therefore not independent proof that the reconstructed sequence contains a stop.
- Default independent-event VEP protein strings omit prediction parentheses.
  Parentheses alone do not change the asserted protein edit.

The [compatibility policy](src/duckvep/kernel/src/duckvep_compat.h) is the single inventory
of explicitly gated runtime behavior. The
[HGVS properties](test/duckvep/property/duckvep_prop_projection_hgvs.c),
[original-record witnesses](test/duckvep/conformance/data/hgvs_compatibility_witnesses.tsv)
and [SQL tests](test/sql/duckvep_hgvs.test) pin the strings and absent-output cases.
Position-zero output is preserved as a VEP convention, not endorsed as valid HGVS.

<a id="equivalent-dna-anchors-can-produce-contradictory-vep-116-protein-descriptions"></a>

### A stop-loss witness with a sequence contradiction

**Classification: reproduced synthetic sequence contradiction and HGVS-rule inference;
unpublished complete evidence, no upstream acknowledgement.**

The forward single-exon transcript `ANCHOR1` spans `chrA1:11–45`, with phase-zero
CDS `11–22`, standard table 1 and complete transcript sequence
`ATGGGTCCTTAAAAAGAACAATAATAACTAGCTGA`. Its CDS `ATGGGTCCTTAA` translates to
`MGP*`. These four physical records are separate inputs, including the duplicate allele:

| Record ID | POS | REF → ALT | VEP HGVSc | VEP HGVSp suffix |
| --- | --- | --- | --- | --- |
| `ANCHOR1_10_T_0` | 19 | T → TT | `c.10dup` | `p.Ter4LeufsTer9` |
| `ANCHOR1_10_T_1` | 20 | T → TT | `c.10dup` | `p.Ter4delinsLeuTer` |
| `ANCHOR1_11_T_0` | 20 | T → TT | `c.10dup` | `p.Ter4delinsLeuTer` |
| `ANCHOR1_11_T_1` | 21 | A → TA | `c.10dup` | `p.Ter4delinsLeuTer` |

Every record reconstructs `ATGGGTCCTTTAAAAAGAACAATAATAACTAGCTGA`.
Independent base-R and BioPerl translation agree on `MGPLKRTIITS*`, first stop at
position 12. The delins description instead asserts a stop at position 5.
The traced mechanism is an incomplete local codon rendered as Xaa and then Ter;
unmodified and observed VEP output agreed at buffer sizes 1 and 5,000.

Our [HGVS 21.1.4 extension-rule](https://hgvs-nomenclature.org/21.1.4/recommendations/protein/extension/)
inference is `p.(Ter4LeuextTer9)`: extension has priority when the reference protein
is extended. The frameshift string has the reconstructed termination distance but
a different operation; the delins string asserts a stop absent from that translation.
This does not establish prevalence, clinical impact or a general stop-loss defect.

The complete four-record outputs and trace are in the untracked receipts
`hgvs_anchor_contract_NTMkInFP/receipt.json` and
`hgvs_anchor_trace_p2DHeEUD/mechanism_receipt.json`, both under
`test/duckvep/conformance/results/`. These are local, build-unbound diagnostics,
not a published conformance pack or release-build certificate. The shipped
[terminal-anchor differential](test/duckvep/conformance/hgvs_anchor_differential.R)
provides the broader executable comparison. DuckVEP's compatibility target remains
each original record's VEP result; this inference does not replace it.

<a id="haplosaurus-reference-and-alternate-proteins-use-different-stop-rules"></a>

## Haplotype sequences, reference curation and flags

**Classification: observed Haplosaurus conventions and unresolved grouped metadata.**
Core reference translation removes the last complete translated stop, applies legitimate
start methionine and reference peptide edits, and retains internal stops. Haplosaurus
then appends `*` only for an exact uppercase raw-CDS suffix `TAA`, `TAG` or `TGA`,
regardless of table or frame. Alternate translation has no reference edits or start
override and displays only the first-stop prefix.

For example, table-1 `CTGGCCTAA` has reference `MA*` but no-edit alternate `LA*`;
table-2 `ATGGCCTGA` has reference `MAW*` but alternate `MAW`.
Such protein differences do not prove a causal genomic edit at the differing residue.
The [reference-translation evidence](test/duckvep/conformance/data/reference_translation_consensus)
covers these conventions; single-source HGVS remains separate from curated-reference
protein differences.

Raw `source_records` replay follows literal parser alleles, not decoded singleton
normalization. Haplosaurus skips non-ACGT mutation alleles, so the retained `N>NGCC`
example does not insert GCC on that route. DuckVEP preserves the skipped contributor
and conditional evidence with zero physical edits; HGVSp is NULL and input incomplete.
A validated REF slot is not a skipped alternate. Within the supported raw alphabet,
N, U and lowercase alternate bases skip; unsupported symbols remain explicit limitations.
Raw `GT=1` is not decoded haploidy: Haplosaurus retains two file lanes, and its
undefined second slot replaces the complete REF with an empty ALT. DuckVEP keeps
this input-route distinction explicit rather than inferring a second called allele.
The [raw observations](test/duckvep/conformance/raw_indel_observations.R) and
[SQL haplotype tests](test/sql/duckvep_haplotypes.test) retain these distinctions.

Haplosaurus groups lanes by final sequence but copies indel/frame flags from the first
lane encountered. Later members do not combine those flags. The
[grouped-flag experiment](test/duckvep/conformance/haplotype_grouped_flags.R), recorded in
its [results ledger](test/duckvep/conformance/data/haplotype_grouped_flags_history.csv), observes the
same 180-base sequence group with `has_indel=0` for 11 hash seeds and `1` for 21,
with identical memberships and agreement between repeats. This is order-dependent
metadata, not a biological consensus rule. DuckVEP path flags and upstream group flags
remain distinct; the four full-output disagreements in the
[publication audit](benchmarks/duckvep_conformance.md#repeated-model-publication-audit)
are not waived by matching sequences and counts. The
[pinned container implementation](https://github.com/Ensembl/ensembl-variation/blob/2fb834b987ede3824e200197a838ce11e91aeb4b/modules/Bio/EnsEMBL/Variation/TranscriptHaplotypeContainer.pm)
is the upstream authority; full grouped-metadata conformance remains unresolved.
