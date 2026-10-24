# Haplotype prediction contract: `duckvep-coding`

`duckvep_haplotypes` rebuilds each haplotype's transcript and classifies the whole product. This page states the rules it follows, the authority each rule is checked against, and where the evidence lives. Per-column details are in the [function reference](../docs/functions.md#duckvep_haplotypes).

## Authorities

- **Whole-haplotype consequences: pinned `bcftools csq`.** VEP annotates each record on its own, and its haplotype mode (Haplosaurus) reports sequences, differences and flags but no SO/IMPACT contract. csq classifies the combined haplotype directly (same-codon compounds, frame opening and restoration, stop gain and loss, start loss), so it is the executable check wherever both tools can express the answer.
- **Everything outside csq's reach: VEP itself.** Where csq has no answer (non-diploid calls, other genetic codes, open-ended CDSs, curated translations, per-record terms), single-record paths must equal DuckVEP's per-variant annotation, which is exact against VEP. Composition rules are checked by independent base-R or Python appliers and hand-derived cases.
- **NMD: the `ejc50` rule below.** It is DuckVEP's own heuristic with its own tests; csq and VEP do not predict it.
- **Compound HGVS: HGVS rules**, checked against offline Mutalyzer and against executable VEP on the equivalent merged record.

### Comparing with csq

csq reads the Ensembl 116 GRCh38 GFF3 (sha256 `08e881d96ab6385a2c31f063a018be4b2c36860b323f2724be07022deeef21ce`) and DuckVEP reads its Ensembl model; bcftools/htslib commits, the FASTA and csq's `-p a` phasing are pinned in the [HG002 oracle receipt](../test/data/haplotype/hg002_domain.oracle.tsv).

- **Vocabulary.** The [csq→SO table](../test/data/haplotype/csq_so_map_v1.tsv) maps each csq term to one SO term or marks it as having no equivalent. Comparisons use mapped term sets and amino-acid changes keyed by transcript, sample and haplotype lane.
- **Domain.** A record is compared only when it overlaps a transcript whose exon and CDS geometry is identical in both models, and only on complete calls with unambiguous phase. Every input record is accounted for as compared, outside csq's domain, a model disagreement or a phase disagreement; nothing is dropped and an out-of-domain record is never a consequence disagreement.

## Input

- **Alleles:** literal A/C/G/T SNVs, MNVs and indels of any length, bounded only by the workspace limits.
- **Phase, `strict`:** GT and PS as written. Each carrier lane (one chromosome copy) of any ploidy is a haplotype. Phased calls without PS share an implicit phase set. An unphased heterozygous call that is the sample's only heterozygous or missing call on the transcript needs no phase: the two haplotypes are the one with the allele and the one without. A missing allele, two or more unphased heterozygous sites, or phase sets that cannot be related are `incomplete_input`; [`duckvep_haplotype_arrangements`](../docs/functions.md#duckvep_haplotype_arrangements) enumerates the alternatives as labelled hypotheses.
- **Phase, `vep_compat`:** genotype slots are read in order, as VEP's haplotype mode reads every genotype as phased; paths are then classified exactly as under `strict`. Compound HGVSc is emitted only for proven cis (`strict`).
- **Raw source records** (`input_mode := 'source_records'`) replay VCF records the way VEP's haplotype mode does, byte for byte. Prediction there requires complete diploid calls with at most one heterozygous record, so that each allele slot maps to a real copy.

## Composing records

Every record applies at its own coordinates. A reference base is kept, substituted or removed, and an interbase site receives at most one inserted sequence.

- Records that agree compose once, and every record keeps its own REF, ALT, position and carrier identity. Unchanged bases inside a substitution are context, not a claim.
- Two records that claim one base or site differently are an `edit_conflict`: two different substitutions of one base, a base one record removes and another changes, or two different insertions at one site.
- An opaque delins (unequal, nonzero REF and ALT) is replaced as a unit. An identical copy composes, an insertion at its first site precedes it, and any other record that changes or removes one of its bases conflicts.
- Insertions at one spliced site from different exons lie in different genomic gaps and apply in transcript order.
- An equivalent placement of an indel elsewhere in a repeat is not a second reading of the record.

Evidence: `duckvep_compatible_overlap.test`, `test_compatible_overlap_indels.R`, `test_opaque_delins_aliases.R`, and the [compound fuzz](../test/duckvep/conformance/compound_hgvs_fuzz.py), which checks composition against an independent literal applier on random overlapping records (9,000 cases, both strands, no mismatch).

## The product

- **CDS.** The edited CDS is translated with the transcript's genetic code. A CDS needs a start codon and a terminal stop unless its start or end is not annotated (`cds_start_NF`, `cds_end_NF`). A CDS whose sequence lacks a start codon, a terminal stop or a whole last codon is read as if that end were not annotated, as VEP translates it; a start that begins inside a codon is stored padded with N, as Ensembl translates it.
- **Lost stop.** Translation continues through the transcript's stored 3′ flank to the next stop or the end of the flank; nothing beyond the stored sequence is invented.
- **Curated translations** (selenocysteine, stop-codon readthrough, initial methionine, amino-acid substitution) and **RNA editing programmes** are part of the reference, as in Ensembl's translation. An unchanged curated codon keeps its programme; a changed one translates ordinarily. The paired peptides are in `prediction_reference_protein` and `prediction_protein`; `cds` and `protein` always hold the raw replay.
- **Splicing is the model's**, as in VEP. A record that touches no exon leaves the product alone and stays a contributor with role `omitted`. A record in a frameshift intron (an intron of at most 13 bases) is coding sequence, as VEP reads it.
- **Noncoding transcripts and UTRs** are replayed on the final RNA (substitutions, indels and compensating edits within one exon, away from splice junctions) and classified from the reference and alternate RNA.

## Classification

| Situation | Terms |
|---|---|
| Edited CDS no longer starts with a start codon | `start_lost` alone: initiation is unknown, so nothing downstream is predicted |
| First stop before the reference terminator, peptide before it changed | `stop_gained`, plus `frameshift_variant` when the stop lies in a displaced frame |
| No stop where the reference terminates | `stop_lost`, plus `frameshift_variant` when the frame is still displaced |
| Frame displaced and not restored before termination | `frameshift_variant` |
| Frame restored before termination | no frameshift term: `protein_altering_variant`, or `synonymous_variant` if the peptide is unchanged |
| Same peptide, terminal codon changed | `stop_retained_variant` |
| Same peptide otherwise | `synonymous_variant` |
| Substitutions only, peptide changed | `missense_variant` |
| Pure in-frame insertions or deletions | `inframe_insertion` / `inframe_deletion` |
| Annotated start changed to another start codon | adds `start_retained_variant` |

Beside the dominant term, effects reached in the final product are added: a missense or synonymous change in a codon translated before the first stop, and an in-frame insertion or deletion whose whole tract is translated. Changes after the first stop stay contributors (role `post_stop`) and add no term. Same-codon changes are classified as the codon they make together. Reference lanes have no terms. IMPACT is the highest impact of the emitted terms.

Evidence: pinned csq through the mapping table, and independent base-R goldens for same-codon, frame, start/stop and restoration cases (`duckvep_haplotype_same_codon.test`, `_frame.test`, `_startstop.test`).

## Per-record terms

Some paths have no whole product to rebuild: an edit that crosses the CDS start or end or an exon boundary, a record that changes an essential splice site (VEP `splice_donor_variant` or `splice_acceptor_variant`), an in-frame stop in the reference CDS, a curated programme the lift cannot follow, or overlapping RNA edits outside the RNA replay. Those paths are `predicted` with reason `per_record_terms`. Their consequences are the union of the records' terms on the transcript from the per-variant engine (`duckvep_annotate`); raw replay fields stay, and NMD is not predicted.

Evidence: on HG002 every single-record path of this kind equals the per-variant annotation term for term and in IMPACT; `duckvep_haplotype_per_record.test` covers a donor-site record beside a missense change, a deep intronic record, and agreeing and conflicting RNA edits.

## NMD (`ejc50`)

For a new premature stop: `trigger` when J − S > 50, otherwise `escape`. S is the stop's last base and J the last base of the penultimate exon, in edited spliced-transcript coordinates. Intronless transcripts escape. A termination that does not move is `not_applicable`; lost initiation, unknown termination or unavailable exon geometry is `unknown`. UTR length changes move S and J. The row reports the stop and junction positions and the contributing records.

Two exceptions are reported beside the rule in `nmd_exceptions`, without changing it: `start_proximal` (the stop codon is within the first 100 coding bases) and `long_exon` (the stop lies in an exon longer than 407 bases). Whether either is an escape is the caller's decision.

Evidence: hand-derived cases at 49/50/51 bases, intronless, last-exon and indel-shifted junctions (`duckvep_haplotype_nmd.test`, `duckvep_haplotype_nmd_exceptions.test`).

## Compound HGVS

With `hgvs := true`, a path with several records under `strict` phase gets `hgvsc` and `hgvsp` for the joint allele. Changes in one codon, adjacent changes, and a change whose 3′ shift reaches the next one are one delins; other changes are listed as `c.[a;b]`. Shifts may run into `c.*`. The protein name is built from the joint allele in its minimal form; an introduced stop is a nonsense substitution. Evidence: the [runtime cases](../test/duckvep/conformance/compound_hgvs_runtime.R) and the [compound fuzz](../test/duckvep/conformance/compound_hgvs_fuzz.py): DNA names agree with Mutalyzer's normalization on 9,000 random cases, and where the joint allele is one change the protein name equals VEP's for the merged record except for the nonsense-substitution rule and stop-codon frameshift notation.

## Statuses

| `prediction_status` | Meaning |
|---|---|
| `predicted` | Terms, IMPACT and NMD are set (reason `supported_domain`), or per-record terms are set (reason `per_record_terms`) |
| `incomplete_input` | A missing allele, or phase the calls do not define |
| `edit_conflict` | Two records claim one base or site differently |
| `unsupported_overlap` | Overlapping records on a path the composition rule does not cover (raw source replay, a record spanning two exons) |
| `unsupported_context` | The reference or input cannot support a prediction: reference mismatch, ambiguous base, missing transcript flank, untyped curated metadata, capacity |

Every row keeps all its contributors (applied, omitted, shadowed and post-stop) with their original operands, carriers and phase evidence. Budget or capacity exhaustion is an explicit error, never a truncated result.

## Scale

Requirements for the whole-genome job: a fresh process per run, 16 GiB process memory including R, a 4 GiB native budget covering models, indexes and workers, spillable DuckDB staging under a disk quota, and complete output with a full checksum. The speed requirement is at least twice pinned csq on the same input, model and core. The [HG002 runner](../benchmarks/haplotype_scale/run_hg002_csq.sh) measures csq and DuckVEP alternately on one pinned core; its receipts and the [5M-record qualification](../benchmarks/data/haplotype_scale/README.md) record the results and conditions.
