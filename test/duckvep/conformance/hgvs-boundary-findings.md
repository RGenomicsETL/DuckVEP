# Protein HGVS boundary fixtures

## Scope

These fixtures test VEP 116 compatibility for the protein-HGVS boundaries in issues [#69](https://github.com/RGenomicsETL/DuckVEP/issues/69) and [#70](https://github.com/RGenomicsETL/DuckVEP/issues/70).
The oracle uses restored transcript objects through VEP's GFF path.
It is not a full-cache comparison or experimental validation of an altered protein.
Consequence sets, HGVSc, HGVSp, and their support statuses are separate results.

## Source identities

| Resource | Identity |
| --- | --- |
| ensembl-vep | `57ea5c52340acc1f156267f810ad162e26597082` |
| ensembl-variation | `2fb834b987ede3824e200197a838ce11e91aeb4b` |
| ensembl-core | `c0cf13daa961d80584bad797b2eb0ff3a7500ef3` |
| ensembl-io | `6afb5dc27a5ae6881d75959153fbe6e9a4a7e788` |
| ensembl-funcgen | `90049ea7ee4d8ae3a6d298dca46d6c6ab20538c4` |
| HF model repository | `RgenomicsETL/duckvep-models` |
| HF revision | `a56590bdd03f19dd31ca0e7f5512b9a9aadd9072` |
| Model path | `ensembl-116/homo_sapiens_116_GRCh37.duckdb` |
| Model SHA-256 | `f397a471ebd2a3b36844b0e94bfa995ad3fb0cfb042d8d3e77d9fcbea7304bd6` |
| Local GRCh37 FASTA SHA-256 | `3a3872e7bdd1532fdbfbc3afd70c04ff8a12a109319ae050c91d962f32e16a4f` |

The original DuckVEP outputs were reproduced with DuckDB R 1.5.5.
Every input REF matched the local FASTA.
Independent exon extraction matched the model's full spliced and CDS sequences for all nine transcripts.
The FASTA checksum identifies this experiment's reference, not an asserted byte match to the issue author's reference.

The restored VEP objects preserve transcript and translation identity, biotype, exon phases, end phases, and relevant attributes.
The oracle checks spliced sequence and phase-aware CDS sequence before annotation.
The restoration does not replace consequence predicates, HGVS predicates, or translation routines.

The nine expected HGVS pairs pass on stable-v1 with DuckDB R 1.5.5 and native-v2 with DuckDB R 2.0.0.
Across all 92 annotations for these input alleles, only nine protein strings and eight protein support statuses change.
All other columns, including HGVSc and consequence masks, remain identical.

## Transcript-start deletion

The GRCh37 allele is `14:102786284 TGCTGCGGGACGCGGCGCGGCG>T`, on `ENST00000559098.1`.
The genomic deletion has 21 bases.
Its original interval includes ten bases before the transcript and eleven CDS bases.
The initial coding exon phase is zero, and the CDS start is incomplete.

VEP moves the complete genomic feature 15 bases toward the transcript's 3-prime end for protein replay.
That placement deletes CDS positions 6 through 26.
HGVSc independently clips the original feature at the transcript boundary and places its eleven-base span at `c.16_26del`.
Thus the DNA notation and protein replay do not use the same edit span.

The reference CDS starts `GCGGCGCGGCGGCTGCGGGACGCGGCCCAACAG...`, with peptide `AARRLRDAAQQ...`.
The shifted 21-base deletion gives `GCGGCCCAACAG...`, with peptide `AAQQ...`.
The peptide deletion is `p.Arg3_Ala9del`.
A restored-object run without diagnostic tracing also gives this result.

The implementation applies this distinction to incomplete-CDS transcript-start deletions whose complete shifted span fits within one coding exon and the reference window.
Complete-CDS and multi-exon extensions need separate evidence.
The test uses a reference-derived CDS segment, synthetic flanks, and a reverse-complement control.

Fixture: [`duckvep_hgvs_external_feature.test`](../../sql/duckvep_hgvs_external_feature.test).
It asserts both strands, HGVS strings, support statuses, shift, named consequences, SNV controls, and the complete-CDS scope boundary.

## Short introns and the unshifted coding predicate

VEP marks a transcript with a short-intron flag when a positive intron has at most 13 bases.
It then extends exon-overlap intervals by 12 bases for that transcript.
An intronic CDS mapper gap can consequently satisfy the original coding predicate.
The complete shifted-feature translation-bounds check still applies.

This rule accounts for the eight missing proteins in issue #69.
Their baseline transcript HGVS status is `supported`, but their protein HGVS status is `not_applicable` with a NULL reason.
The scoped fix produces the eight reported VEP proteins without changes to their consequence masks or HGVSc values.

The tested controls require the coding precondition.
An independent linearized oracle keeps all six affected circular-control proteins absent.
The linearization checks the original alleles, exon sequences, phase-aware CDS sequences, and reference peptides.
It does not establish equivalence for arbitrary circular shifts or mitochondrial genetic codes.

Fixture: [`duckvep_hgvs_short_intron.test`](../../sql/duckvep_hgvs_short_intron.test).
It covers 2-, 13-, and 14-base introns, both strands, and a UTR negative control.
The [`duckvep_circular_lifted.test`](../../sql/duckvep_circular_lifted.test) retains its protein count.

## First-residue frameshift termination

VEP's frameshift formatter calls `_stop_loss_extra_AA(start - 1)`.
That function rejects zero, so a first-residue frameshift uses `Ter?`.
DuckVEP applies the same convention in ordinary and replay paths.
This changes notation, not the translated alternate sequence.

The executable differential covers both strands, shifted deletions, and later-position controls.
The mitochondrial fixture has two affected first-residue protein strings among 248,673 rows.
Its other columns are unchanged.
The residue distribution claimed for the unavailable cloud fixture was not independently verified.
An independent 28-row paired-strand probe matched VEP and DuckVEP HGVS values.
Literal and normalized alleles gave equal alternate transcripts in the tested pairs.
This agreement does not independently validate every shared HGVS value against an alternate peptide.
Its comparable minus-strand stop deletion gave Lys in both tools.
It did not reproduce the cloud report's Ser/Lys discrepancy.
The original cloud FASTA, GFF, VCF, invocation, and raw output are still needed for that observation.
No upstream erratum or additional production change is justified by this bounded probe.

Fixtures: [`hgvs_first_frameshift_differential.R`](hgvs_first_frameshift_differential.R), the `hgvs_frameshift_first_residue_termination` property, and [`duckvep_circular_mt.test`](../../sql/duckvep_circular_mt.test).

Run the differential from the repository root with DBI, DuckDB R, micromamba, samtools, bgzip, and tabix available.
Select a VEP 116 Conda environment and an output directory that does not exist.

```sh
DUCKVEP_VEP_PREFIX=/path/to/vep-environment \
Rscript test/duckvep/conformance/hgvs_first_frameshift_differential.R \
  build/release/duckvep.duckdb_extension /path/to/new-results fixed
```

## Stop-codon boundary deletions

[`stop_codon_sweep.sh`](stop_codon_sweep.sh) compares both engines on the VEP cache path.
It takes every 1 to 4 bp deletion that starts near a stop codon of one chromosome.
On chromosome 21 of GRCh38 this is 1,109 stop codons, 77,900 alleles and 2,846,960 transcript pairs.
Consequence sets, HGVSc and HGVSp agree on every pair.
The first run of this sweep found ten differing rows ([#75](https://github.com/RGenomicsETL/DuckVEP/issues/75)); four rules account for them.

1. **No splice term outside the transcript.**
   VEP sets its intron and splice pre-predicates only when the feature overlaps the transcript.
   A terminal exon shorter than three bases otherwise leaves part of the exonic splice-region window outside it.
2. **Stop overlap is genomic and unspliced.**
   For a deletion with an `X` peptide allele, `stop_retained` requires `_overlaps_stop_codon_cil`.
   That helper tests the three genomic bases that end at the coding-region end.
   When an intron splits the stop codon, the stop bases before the intron are outside that window.
   VEP then reports `coding_sequence_variant`, although the edited sequence keeps a stop.
3. **Unchanged through the stop.**
   A frameshift whose new peptide equals the reference through the stop is named at the stop (`p.Ter308=`).
4. **Frame test after a shift across a mapper gap.**
   An original feature that crosses an exon junction has no cached peptide alleles.
   When the 3-prime shift moves it wholly into the CDS, VEP's frameshift test has no leading-stop exclusion and decides by length.
   The protein is then a frameshift (`p.Ter81AsnfsTer13`), not the stop equality VEP writes for the same change without a shift.

Rules 2 and 3 were read from an instrumented VEP 116.0 run and its `VariationEffect.pm`.
Rule 4 is inferred from the installed source and confirmed only by the two chromosome 21 rows; no synthetic fixture reproduces the shift on VEP's GFF path.

Fixtures: [`duckvep_stop_boundary.test`](../../sql/duckvep_stop_boundary.test) (six rows equal to VEP 116.0 on its GFF path, rules 2 and 3) and the `splice_region_requires_within_feature` property (rule 1).
