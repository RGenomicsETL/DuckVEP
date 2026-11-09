# VEP 116 structural-variant corpus audit

## Boundary records

`vep_rs_sv_boundary.R` selects the lengthless `<INS:ME...>` records and the
cross-contig breakends of both corpora. Run on 2026-11-09 with the stable-v1
host, DuckDB R 1.5.5 and the Ensembl 116 models:

| Assembly | Records | VEP tuples | DuckVEP tuples | Differences | Elapsed | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| GRCh38 | 6 | 233 | 233 | 0 | 12.3 s | 4.06 GiB |
| GRCh37 | 5 | 58 | 58 | 0 | 23.0 s | 0.87 GiB |

A lengthless symbolic insertion is VEP's one-base feature at `POS + 1`.
`duckvep_vcf_events_sql` marks it with `structural_geometry = 'POINT'`, and the
annotation reads that feature as an overlap rather than an interbase site.
The elapsed time and memory include the model load.

## Complete corpora

In the complete-corpus run of the same date, every projected difference is a
breakend. VEP emits each breakend transcript twice, once for the
single-breakend allele and once for the input ALT; DuckVEP emits one event row.
For 324 GRCh38 and 192 GRCh37 keys the two VEP rows have equal terms and DuckVEP
has that row once. For 19 and 17 keys the two VEP rows differ, and the DuckVEP
row equals the VEP row that carries `feature_truncation`. No transcript is
missing from either engine, and no DuckVEP row lacks an equal VEP row.
The table and text below are from the earlier audit.

This diagnostic accounts for the GRCh37 and GRCh38 VCF inputs from vep-rs
revision `70df135effe286bc96c8dc93381bb06f2e2979cf`. It uses the captured
Ensembl VEP 116.2 JSON and the Ensembl 116 models. The input hashes, model
hashes, extension hash and comparison limits are in each generated
`receipt.json`. The audited stable-v1 extension had SHA-256
`97817d081cc1993f54c43ad355aa1942bfc3efba538b1832e52cad8aa56e9e60` and
ran on DuckDB 2.0.0. It contained the working-tree insertion-geometry changes;
it was not a released artifact. The diagnostic script's SHA-256 was
`6879ac26afcf7fb5084a5a00fc6b58fe810a3802bad2f326d2985c06d4f7638b`.

| Measure | GRCh38 | GRCh37 |
|---|---:|---:|
| Input records | 817 | 1,008 |
| Input alleles | 817 | 1,028 |
| Unsupported input alleles | 38 | 60 |
| Records without captured VEP JSON | 0 | 6 |
| VEP transcript rows on unsupported or unmapped inputs | 2,257 | 3,400 |
| Projected VEP transcript-term tuples | 18,364 | 14,133 |
| Projected DuckVEP transcript-term tuples | 18,021 | 13,924 |
| Projected multiset intersection | 18,021 | 13,924 |
| Projected VEP-only / DuckVEP-only | 343 / 0 | 209 / 0 |

The tuples aggregate by input record, input ALT, transcript, and sorted terms;
they retain multiplicity. They collapse VEP's separate breakend endpoint
alleles onto the input ALT. Therefore the projected differences are **not a
count of confirmed biological errors** and the intersection is not full
allele-level parity. In the supported breakend rows, VEP labels 373 GRCh38
transcript rows with a single-breakend allele and 397 with the original ALT;
the corresponding GRCh37 counts are 240 and 237. The generated
`breakend-alleles.tsv` retains both labels for review.

For the GRCh38 record `21:25880900 G>[21:25881361[G`, VEP emits 80 transcript
rows: 40 for `.G` and 40 for the input ALT, with 40 distinct transcript IDs.
The two labels can have different consequence terms. DuckVEP reports the
breakend as an event-level annotation, so assigning both VEP labels to the one
DuckVEP row obscures that distinction.

A separate offline probe ran
`ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da`
with `--buffer_size 1` and `--buffer_size 5000`, `--fork 1`, the GRCh38 cache
and the same options as the
corpus provenance. For the isolated record and the 11 affected GRCh38 records,
respectively, the two runs had identical multisets of 80 and 687 tuples keyed
by input, transcript, VEP `variant_allele`, and sorted terms. This rules out
buffer size as the cause in those probes; it does not decide how the two
endpoint alleles should map to DuckVEP's event API.

This audit compares literal-SNP location and allele spelling only. It retains
no-transcript rows but does not establish parity for them. It does not compare
regulatory or other non-transcript features, load a reference FASTA, or establish
parity for unsupported and ambiguous multiallelic inputs. Full VEP parity and
the historical issue baseline require further adjudication. Run the diagnostic
using the instructions in [README.md](README.md).
