DuckVEP loss-of-function parity against VEP 116 with LOFTEE
================

This report records the parity check of `duckvep_lof_sql` against the
real LOFTEE plugin. The builder states LOFTEE's rules as joins and a
`CASE` over rows DuckVEP already produces; the plugin is the oracle. The
gate is exact agreement of `LoF`, `LoF_filter` and `LoF_flags` for every
variant x transcript pair LOFTEE reports, with missing and extra pairs
counted. `LoF_info` is compared as well and is reported, not gated.

Everything here is reproduced by `scripts/lof_parity.py`; the receipts are
in `benchmarks/data/duckvep_lof/`.

## Setup

| Item | Value |
| --- | --- |
| VEP | 116, image `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11`, offline, `--distance 5000`, JSON output (`scripts/run_vep116_loftee_docker.sh`) |
| Cache | Ensembl 116 GRCh38 chr21 indexed cache (`cache-grch38-chr21`) |
| Reference | `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl 116) |
| LOFTEE | konradjk/loftee at `a46b502a68c812c8ae0c5a5721c0603fe81cae8d`, mounted read-only |
| DuckVEP model | `homo_sapiens_116_GRCh38_final.duckdb` (Ensembl 116 GRCh38), annotation by `duckvep_annotate_projected_sql` |
| Profile | LOFTEE defaults. Human ancestor off, PhyloCSF (`conservation_file`) off, `run_splice_predictions` off (its default). GERP is the constant-negative bigwig `gerp_const_neg.bw` (-1000 over chr21), so LOFTEE's GERP-weighted rule reduces to the unweighted 50 bp rule. |

LOFTEE at this commit always takes its GERP branch (`use_gerp_end_trunc`
is the truthy string `'true'`) and opens a bigwig, so a GERP file is
required to run it at all; the constant makes `GERP_DIST <= -58` hold
whenever the weighted distance is positive. The "GERP relation" runs below
give `duckvep_lof_sql` the same constant as a `gerp` relation, which
exercises the weighted path: there `GERP_DIST` in `LoF_info` is compared
too.

## Corpus

- **real**: 6,845 chr21 alleles of GIAB HG002 v4.2.1 and ClinVar
  (20260706) whose DuckVEP annotation, on a protein-coding transcript,
  includes `stop_gained`, `frameshift_variant`, `splice_acceptor_variant`,
  `splice_donor_variant`, `splice_region_variant`, `start_lost`,
  `stop_lost` or `protein_altering_variant`. The net is wider than the
  four LoF terms so that a pair LOFTEE calls and DuckVEP does not would
  show up as missing. Multi-allelic records are split into one record per
  allele, because LOFTEE's GC-to-GT test reads the first two alleles of the
  variant. `duckvep_lof_corpus_real.vcf`.
- **synthetic**: 177 alleles on real chr21 transcripts, both strands,
  chosen by `scripts/lof_parity.py` (the `ID` column of
  `duckvep_lof_corpus_synthetic.vcf` names the case):
  END_TRUNC boundary (1-bp deletions and insertions in the penultimate
  exon at 49, 50 and 51 bases from its 3' end, where the stop codon lies
  in the last exon), last coding exon, a stop exon followed by a 3' UTR
  exon, single-exon genes, splice donor and acceptor SNVs of the smallest
  introns, NAGNAG hits (SNV, deletion, insertion) and misses, `GC` donors
  (C>T at intron position 2, other alternates, position 1), non-canonical
  introns, and 2-bp deletions across exon-intron boundaries.

chr21 has no protein-coding intron shorter than 60 nt, so the 14/15/16-nt
boundary of SMALL_INTRON cannot be hit on real transcripts. The synthetic
runs instead move the threshold: for the smallest introns of each strand
(60, 74, 75 and 78 nt) LOFTEE runs with `min_intron_size` S-1, S and S+1,
and `duckvep_lof_sql` with the same `{min_intron_size: ...}`. The test is
the same strict `size < threshold`. The exact 14, 15 and 16 nt introns are
in `test/sql/duckvep_lof.test`, whose fixture is hand-built.

## Results

Pairs are variant x transcript pairs for which LOFTEE reports a result.
"Exact" is agreement of `LoF`, `LoF_filter` and `LoF_flags`; "info exact"
is agreement of `LoF_info` (without `GERP_DIST`, except in the GERP
relation runs).

| Run | Alleles | LOFTEE pairs | DuckVEP pairs | Exact | Missing | Extra | Different | Info exact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| real | 6,845 | 59,668 | 59,668 | 59,668 | 0 | 0 | 0 | 59,668 |
| synthetic | 177 | 2,499 | 2,499 | 2,499 | 0 | 0 | 0 | 2,499 |
| synthetic, `min_intron_size` 59, 60, 61 (60-nt introns) | 177 | 2,499 each | 2,499 | 2,499 each | 0 | 0 | 0 | 2,499 |
| synthetic, `min_intron_size` 73, 74, 75, 76 (74- and 75-nt introns) | 177 | 2,499 each | 2,499 | 2,499 each | 0 | 0 | 0 | 2,499 |
| synthetic, `min_intron_size` 77, 78, 79 (78-nt introns) | 177 | 2,499 each | 2,499 | 2,499 each | 0 | 0 | 0 | 2,499 |
| real, `gerp` relation (constant -1000) | 6,845 | 59,668 | 59,668 | 59,668 | 0 | 0 | 0 | 59,668 (with `GERP_DIST`) |
| synthetic, `gerp` relation (constant -1000) | 177 | 2,499 | 2,499 | 2,499 | 0 | 0 | 0 | 2,499 (with `GERP_DIST`) |

The per-run rows, including the HC and LC counts, are in
`duckvep_lof_parity_summary.csv`. Calls in the default runs: real 53,697
HC and 5,971 LC; synthetic 2,122 HC and 377 LC. The synthetic threshold
runs move between 2,122 and 2,058 HC, so the boundary is exercised.
`duckvep_lof_parity_tokens.tsv` counts the pairs per filter and flag
(real: END_TRUNC 5,525, 5UTR_SPLICE 410, 3UTR_SPLICE 29, GC_TO_GT_DONOR 7,
NAGNAG_SITE 372, NON_CAN_SPLICE 57, NO_EXON_NUMBER 203, SINGLE_EXON 9;
synthetic: END_TRUNC 238, 5UTR_SPLICE 57, GC_TO_GT_DONOR 83, NAGNAG_SITE
396, NON_CAN_SPLICE 503, NO_EXON_NUMBER 100, SINGLE_EXON 12), and
SMALL_INTRON appears in the threshold runs. `duckvep_lof_parity_digests.tsv`
holds the SHA-256 of each run's sorted DuckVEP pair table and of LOFTEE's
JSON output.

**Gated disagreements: none.** `duckvep_lof_parity_disagreements.tsv` has
no rows.

### Informational run: `check_complete_cds`

| Run | LOFTEE pairs | Exact | Different |
| --- | ---: | ---: | ---: |
| real, `check_complete_cds:true` and `{check_complete_cds: true}` | 59,668 | 23,610 | 36,058 |

All 36,058 differences are one class (`duckvep_lof_parity_classes.tsv`):
LOFTEE adds `INCOMPLETE_CDS` and DuckVEP does not. The reason is in
LOFTEE's `check_incomplete_cds`: it tests
`defined($transcript->get_all_Attributes('cds_start_NF'))`, and
`get_all_Attributes` returns an array reference (empty when the attribute
is absent), which is always defined. With the option on, LOFTEE therefore
filters every non-single-exon stop or frameshift, complete or not. The
builder implements the intended test, `cds_start_nf OR cds_end_nf`; no
pair differs in the other direction. The option is off by default in
both, and this run is not part of the gate.

## Decisions and known differences from LOFTEE

- **5UTR_SPLICE and 3UTR_SPLICE** are implemented although they are not
  in the requested rule list: LOFTEE applies them to splice donor and
  acceptor variants lying entirely before the CDS start or after the CDS
  end, and they account for 439 real pairs.
- **END_TRUNC distance** follows `get_gerp_weighted_dist` literally: a
  per-exon length of `end - start` (one less than the exon length), exons
  after the stop exon counted whole, the stop exon counted to the stop
  codon, and `-1000` for the last coding exon length when the stop exon is
  not found. For an insertion, LOFTEE's `cds_end` is the lower of the two
  projected CDS positions (only `PERCENTILE` in `LoF_info` uses it).
- **EXON_INTRON_UNDEF** cannot fire in LOFTEE: the test for an undefined
  exon or intron number sits inside a branch that already requires it to
  be defined. The builder keeps the filter for a malformed row (a rank
  without a count) and it does not occur in the corpora.
- **Not compared with LOFTEE**: ANC_ALLELE, PHYLOCSF_WEAK,
  PHYLOCSF_UNLIKELY_ORF, `ANN_ORF`/`MAX_ORF`, and GERP weighting on real
  scores, because the profile leaves them off and no ancestor, PhyloCSF
  database or GERP bigwig is staged. They are covered by the hand-derived
  values of `test/sql/duckvep_lof.test`. The weighted distance is an exact
  per-base sum over the `gerp` relation; LOFTEE's bigwig summary can round
  at zoom levels.
- **Missing resources** never change a call. They are named in
  `lof_unchecked`, for example `GERP_END_TRUNC` when only the unweighted
  rule ran.
- **Not implemented**: MaxEntScan splice-prediction extensions, which
  LOFTEE leaves off by default (`OS` confidence, `DE_NOVO_DONOR`).
- The builder reads the model's transcript relation (exon list, strand,
  CDS bounds, biotype, `seq_region_name`) and a reference-chunk relation,
  because the resident model cannot be read back from SQL.

## Reproduce

```sh
git clone https://github.com/konradjk/loftee.git data/loftee
git -C data/loftee checkout a46b502a68c812c8ae0c5a5721c0603fe81cae8d
make release
configure/venv/bin/python3 scripts/lof_parity.py \
  --extension build/release/duckvep.duckdb_extension \
  --model-db data/models/homo_sapiens_116_GRCh38_final.duckdb \
  --fasta data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --cache ~/.cache/duckhts/vep/cache-grch38-chr21 --loftee data/loftee \
  --gerp benchmarks/data/duckvep_lof/gerp_const_neg.bw \
  --clinvar clinvar_20260706_grch38.vcf.gz \
  --hg002 HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --out benchmarks/data/duckvep_lof
```

The script exits non-zero unless every gated run has no missing, extra or
different pair and equal `LoF_info`. `gerp_const_neg.bw` was written with
pyBigWig (one interval, `21:0-46709983`, value -1000).
