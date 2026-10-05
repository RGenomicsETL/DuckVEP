# Species structural evidence against VEP 116

The structural preparation builders and executable Ensembl VEP 116 differentials cover the closed species matrix
([#5](https://github.com/RGenomicsETL/DuckVEP/issues/5)): **mouse GRCm39** and two non-vertebrates, **fly BDGP6.54** and
**Arabidopsis TAIR10**. The oracle is the digest-pinned image
`ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11` with the pinned indexed caches and
FASTA files of [mouse](mouse_grcm39_sources.md), [fly](fly_bdgp654_sources.md) and [plant](plant_tair10_sources.md). The comparator is
[`species_sv_vep116_differential.R`](../../test/duckvep/conformance/species_sv_vep116_differential.R).

```sh
DUCKVEP_EXTENSION_FILE=/path/to/immutable/duckvep.duckdb_extension \
Rscript test/duckvep/conformance/species_sv_vep116_differential.R \
  SPECIES ASSEMBLY CACHE_VERSION MODEL_DB FASTA CACHE_DIR OUTPUT_DIR [SEED] [PER_STRATUM]
# mus_musculus GRCm39 116 ... 11639 10   drosophila_melanogaster BDGP6.54 116 ... 11654 10
# arabidopsis_thaliana TAIR10 63 ... 11663 10
```

## Corpus

Everything is derived from the pinned model with a seeded hash, at most ten events per stratum (`PER_STRATUM`):

- **DEL, DUP, DUP:TANDEM (TDUP), INV**, exact-span symbolic (`END=`, `SVTYPE=`), over 13 geometry states sampled per type and
  transcript strand: exon start/end crossings, exact and interior exons, intron interiors, spans across an intron and two exons,
  transcript left/right crossings, exact transcript spans, left/right flanks (1-5 kb), and CDS start/end crossings.
- **INS**, symbolic `<INS>` with `END`=`POS`, at exon boundaries, exon interiors, introns, and CDS start/end (one VEP predicate family,
  insertion overlaps).
- **BND**, reciprocal pairs (`MATEID`, `EVENT`), both physical records sent to VEP and DuckVEP separately, four bracket orientations,
  same- and cross-chromosome mates, 12 endpoint landmarks (transcript start/mid/end, flanks, CDS ends, exon starts/ends/mids), sampled per
  landmark, orientation and local-transcript strand. VEP runs with `--buffer_size 1`.
- **Structural HGVS**: 1-40 bp exact-span DEL/DUP/TDUP/INV (400 events per species) at exon edges, exon and intron interiors, and CDS
  start, with reference sequences from the pinned FASTA.

`corpus_manifest.csv` in each directory gives the events per (type, state, strand); `events.csv` lists every event.
Comparison is the set of normalized consequence terms per (event, transcript). DuckVEP annotates the same typed events
(`structural_type`, `copy_change`, `mate_seq_region`/`mate_position` for BND) with `upstream_distance = downstream_distance = 5000` to match
VEP's `--distance 5000`.

## Results

Exact means equal term sets for the (event, transcript) pair; `oracle_only`, `duckvep_only` and `terms_differ` are counterexamples. Directories
under [sv_species/](sv_species/) hold every stratum table (`sv_type_state.csv`, `sv_type_source_strand.csv`, `biotype.csv`, `so_term.csv`),
`disagreements.csv`, the VCFs, and hashes; the multi-megabyte oracle JSON files are identified by SHA-256 only.

**Admitted: every DEL, DUP, TDUP, INV and INS pair in all three species (all 0 counterexamples), all structural-HGVS strings on the
literal-equivalent edits, and all geometry-builder and BND-identity checks.** BND has 37 nonexact pairs in the pooled runs, all in three
records that are exact when run alone; they are retained below and not admitted as pooled agreement.

## Tables

### Mus musculus GRCm39 (Ensembl 116, seed 11639)

Model SHA-256 `ae39ffc9e647d0a096a13b737599d480fae2b94d3a202938a8eeaa44208c62be`; 1159 symbolic SV records, 576 BND records (288 reciprocal pairs), 400 small exact-span events for HGVS. Corpus manifest: [mouse_grcm39_seed11639/corpus_manifest.csv](sv_species/mouse_grcm39_seed11639/corpus_manifest.csv); records: [sv.vcf](sv_species/mouse_grcm39_seed11639/sv.vcf), [bnd.vcf](sv_species/mouse_grcm39_seed11639/bnd.vcf), [hgvs-literal.vcf](sv_species/mouse_grcm39_seed11639/hgvs-literal.vcf).

Oracle output SHA-256: SV `2f83c07f2c5173d81d39fb4509192645eb0e416086ea893a7506bfd3b9967da1`, BND `3d97369344c7a66980562121fa3346793fd102614141a3e3b6fff37a0da29a02`.

Consequence-term sets per (event, transcript), pooled oracle run:

| SV type | exact | oracle_only | duckvep_only | terms_differ | pairs | events |
| --- | --- | --- | --- | --- | --- | --- |
| BND | 46,050 | 0 | 28 | 0 | 46,078 | 576 |
| DEL | 10,780 | 0 | 0 | 0 | 10,780 | 259 |
| DUP | 10,323 | 0 | 0 | 0 | 10,323 | 260 |
| INS | 4,744 | 0 | 0 | 0 | 4,744 | 120 |
| INV | 10,101 | 0 | 0 | 0 | 10,101 | 260 |
| TDUP | 10,600 | 0 | 0 | 0 | 10,600 | 260 |
| **All** | 92,598 | 0 | 28 | 0 | 92,626 | 1,735 |

By SV type and transcript strand:

| SV type | transcript strand | exact | oracle_only | duckvep_only | terms_differ | pairs | events |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BND | -1 | 19,808 | 0 | 28 | 0 | 19,836 | 534 |
| BND | 1 | 26,242 | 0 | 0 | 0 | 26,242 | 544 |
| DEL | -1 | 5,309 | 0 | 0 | 0 | 5,309 | 202 |
| DEL | 1 | 5,471 | 0 | 0 | 0 | 5,471 | 195 |
| DUP | -1 | 4,802 | 0 | 0 | 0 | 4,802 | 192 |
| DUP | 1 | 5,521 | 0 | 0 | 0 | 5,521 | 192 |
| INS | -1 | 2,180 | 0 | 0 | 0 | 2,180 | 83 |
| INS | 1 | 2,564 | 0 | 0 | 0 | 2,564 | 83 |
| INV | -1 | 5,358 | 0 | 0 | 0 | 5,358 | 186 |
| INV | 1 | 4,743 | 0 | 0 | 0 | 4,743 | 204 |
| TDUP | -1 | 4,907 | 0 | 0 | 0 | 4,907 | 196 |
| TDUP | 1 | 5,693 | 0 | 0 | 0 | 5,693 | 209 |

Every record with a pooled disagreement, re-run alone with `--buffer_size 1` (see the BND note below):

| kind | id | pooled_pairs | isolated_pairs | isolated_exact | isolated_nonexact |
| --- | --- | --- | --- | --- | --- |
| bnd | bnd186b | 28 | 29 | 29 | 0 |

Structural HGVS on the literal-equivalent edits (builder `hgvs_g` vs VEP `--hgvsg --shift_hgvs 0`; DuckVEP `hgvsc` vs VEP default `--hgvs`, per event over all transcripts):

| SV type | transcript-source strand | genomic events | genomic equal | transcript events | transcript equal |
| --- | --- | --- | --- | --- | --- |
| DEL | -1 | 50 | 50 | 50 | 50 |
| DEL | 1 | 50 | 50 | 50 | 50 |
| DUP | -1 | 50 | 50 | 50 | 50 |
| DUP | 1 | 50 | 50 | 50 | 50 |
| INV | -1 | 49 | 49 | 49 | 49 |
| INV | 1 | 47 | 47 | 47 | 47 |
| TDUP | -1 | 50 | 50 | 50 | 50 |
| TDUP | 1 | 50 | 50 | 50 | 50 |

Builder statuses: unsupported (inversion_not_reducible) 4; supported 396. VEP's default 3' shift moves 73 of the 396 supported genomic strings, retained in [hgvs_genomic.csv](sv_species/mouse_grcm39_seed11639/hgvs_genomic.csv) (`shift1_differs`).

`duckvep_prepare_sv_geometry_sql` nominal start/end equal VEP's `start`/`end` for all 1159 symbolic records (0 disagreements). `duckvep_prepare_breakend_pairs_sql` on the BND records: reciprocal 576.

### Drosophila melanogaster BDGP6.54 (Ensembl 116, seed 11654)

Model SHA-256 `d84f53924eac8fab69fc8601d38ff650ff5d8bd3d6901e5510742806361114de`; 1158 symbolic SV records, 576 BND records (288 reciprocal pairs), 400 small exact-span events for HGVS. Corpus manifest: [fly_bdgp654_seed11654/corpus_manifest.csv](sv_species/fly_bdgp654_seed11654/corpus_manifest.csv); records: [sv.vcf](sv_species/fly_bdgp654_seed11654/sv.vcf), [bnd.vcf](sv_species/fly_bdgp654_seed11654/bnd.vcf), [hgvs-literal.vcf](sv_species/fly_bdgp654_seed11654/hgvs-literal.vcf).

Oracle output SHA-256: SV `0c5f143b59012c18e5e33a171c347f5405267cd90f3607f53da55d3dd568d948`, BND `798c0d78ec834c55cd08e90bff2c17f5b5c45bd15ba57c3cf0b8eb72e2102954`.

Consequence-term sets per (event, transcript), pooled oracle run:

| SV type | exact | oracle_only | duckvep_only | terms_differ | pairs | events |
| --- | --- | --- | --- | --- | --- | --- |
| BND | 12,916 | 0 | 0 | 0 | 12,916 | 576 |
| DEL | 2,890 | 0 | 0 | 0 | 2,890 | 259 |
| DUP | 2,893 | 0 | 0 | 0 | 2,893 | 260 |
| INS | 1,434 | 0 | 0 | 0 | 1,434 | 120 |
| INV | 2,678 | 0 | 0 | 0 | 2,678 | 259 |
| TDUP | 2,845 | 0 | 0 | 0 | 2,845 | 260 |
| **All** | 25,656 | 0 | 0 | 0 | 25,656 | 1,734 |

By SV type and transcript strand:

| SV type | transcript strand | exact | oracle_only | duckvep_only | terms_differ | pairs | events |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BND | -1 | 7,138 | 0 | 0 | 0 | 7,138 | 564 |
| BND | 1 | 5,778 | 0 | 0 | 0 | 5,778 | 556 |
| DEL | -1 | 1,576 | 0 | 0 | 0 | 1,576 | 234 |
| DEL | 1 | 1,314 | 0 | 0 | 0 | 1,314 | 224 |
| DUP | -1 | 1,562 | 0 | 0 | 0 | 1,562 | 233 |
| DUP | 1 | 1,331 | 0 | 0 | 0 | 1,331 | 220 |
| INS | -1 | 734 | 0 | 0 | 0 | 734 | 104 |
| INS | 1 | 700 | 0 | 0 | 0 | 700 | 100 |
| INV | -1 | 1,401 | 0 | 0 | 0 | 1,401 | 226 |
| INV | 1 | 1,277 | 0 | 0 | 0 | 1,277 | 224 |
| TDUP | -1 | 1,613 | 0 | 0 | 0 | 1,613 | 228 |
| TDUP | 1 | 1,232 | 0 | 0 | 0 | 1,232 | 226 |

No pooled disagreement; no isolated re-run was needed.

Structural HGVS on the literal-equivalent edits (builder `hgvs_g` vs VEP `--hgvsg --shift_hgvs 0`; DuckVEP `hgvsc` vs VEP default `--hgvs`, per event over all transcripts):

| SV type | transcript-source strand | genomic events | genomic equal | transcript events | transcript equal |
| --- | --- | --- | --- | --- | --- |
| DEL | -1 | 50 | 50 | 50 | 50 |
| DEL | 1 | 50 | 50 | 50 | 50 |
| DUP | -1 | 50 | 50 | 50 | 50 |
| DUP | 1 | 49 | 49 | 49 | 49 |
| INV | -1 | 48 | 48 | 48 | 48 |
| INV | 1 | 49 | 49 | 49 | 49 |
| TDUP | -1 | 50 | 50 | 50 | 50 |
| TDUP | 1 | 50 | 50 | 50 | 50 |

Builder statuses: unsupported (duplication_adjacent_repeat) 1; unsupported (inversion_identity) 1; unsupported (inversion_not_reducible) 2; supported 396. VEP's default 3' shift moves 81 of the 396 supported genomic strings, retained in [hgvs_genomic.csv](sv_species/fly_bdgp654_seed11654/hgvs_genomic.csv) (`shift1_differs`).

`duckvep_prepare_sv_geometry_sql` nominal start/end equal VEP's `start`/`end` for all 1158 symbolic records (0 disagreements). `duckvep_prepare_breakend_pairs_sql` on the BND records: reciprocal 576.

### Arabidopsis thaliana TAIR10 (Ensembl Plants 63 / VEP 116, seed 11663)

Model SHA-256 `521c53d0331a1954eadfacf111314ae74c6eeadff92af29c48dad468c87485a6`; 1157 symbolic SV records, 576 BND records (288 reciprocal pairs), 400 small exact-span events for HGVS. Corpus manifest: [plant_tair10_seed11663/corpus_manifest.csv](sv_species/plant_tair10_seed11663/corpus_manifest.csv); records: [sv.vcf](sv_species/plant_tair10_seed11663/sv.vcf), [bnd.vcf](sv_species/plant_tair10_seed11663/bnd.vcf), [hgvs-literal.vcf](sv_species/plant_tair10_seed11663/hgvs-literal.vcf).

Oracle output SHA-256: SV `ec489e6eb1f42a1b43dfc19f077403874513f9a147ab825436f96f079b8685a9`, BND `cc097ab3529b50e5f0d26d2a9975299f481761835f9622d9df119a6818c1cfcb`.

Consequence-term sets per (event, transcript), pooled oracle run:

| SV type | exact | oracle_only | duckvep_only | terms_differ | pairs | events |
| --- | --- | --- | --- | --- | --- | --- |
| BND | 9,526 | 0 | 9 | 0 | 9,535 | 576 |
| DEL | 2,239 | 0 | 0 | 0 | 2,239 | 258 |
| DUP | 2,369 | 0 | 0 | 0 | 2,369 | 260 |
| INS | 990 | 0 | 0 | 0 | 990 | 120 |
| INV | 2,225 | 0 | 0 | 0 | 2,225 | 260 |
| TDUP | 2,330 | 0 | 0 | 0 | 2,330 | 259 |
| **All** | 19,679 | 0 | 9 | 0 | 19,688 | 1,733 |

By SV type and transcript strand:

| SV type | transcript strand | exact | oracle_only | duckvep_only | terms_differ | pairs | events |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BND | -1 | 4,826 | 0 | 6 | 0 | 4,832 | 572 |
| BND | 1 | 4,700 | 0 | 3 | 0 | 4,703 | 572 |
| DEL | -1 | 1,151 | 0 | 0 | 0 | 1,151 | 237 |
| DEL | 1 | 1,088 | 0 | 0 | 0 | 1,088 | 236 |
| DUP | -1 | 1,261 | 0 | 0 | 0 | 1,261 | 240 |
| DUP | 1 | 1,108 | 0 | 0 | 0 | 1,108 | 225 |
| INS | -1 | 458 | 0 | 0 | 0 | 458 | 105 |
| INS | 1 | 532 | 0 | 0 | 0 | 532 | 111 |
| INV | -1 | 1,066 | 0 | 0 | 0 | 1,066 | 233 |
| INV | 1 | 1,159 | 0 | 0 | 0 | 1,159 | 225 |
| TDUP | -1 | 1,162 | 0 | 0 | 0 | 1,162 | 233 |
| TDUP | 1 | 1,168 | 0 | 0 | 0 | 1,168 | 242 |

Every record with a pooled disagreement, re-run alone with `--buffer_size 1` (see the BND note below):

| kind | id | pooled_pairs | isolated_pairs | isolated_exact | isolated_nonexact |
| --- | --- | --- | --- | --- | --- |
| bnd | bnd130b | 5 | 11 | 11 | 0 |
| bnd | bnd230b | 4 | 12 | 12 | 0 |

Structural HGVS on the literal-equivalent edits (builder `hgvs_g` vs VEP `--hgvsg --shift_hgvs 0`; DuckVEP `hgvsc` vs VEP default `--hgvs`, per event over all transcripts):

| SV type | transcript-source strand | genomic events | genomic equal | transcript events | transcript equal |
| --- | --- | --- | --- | --- | --- |
| DEL | -1 | 50 | 50 | 50 | 50 |
| DEL | 1 | 50 | 50 | 50 | 50 |
| DUP | -1 | 50 | 50 | 50 | 50 |
| DUP | 1 | 49 | 49 | 49 | 49 |
| INV | -1 | 48 | 48 | 48 | 48 |
| INV | 1 | 49 | 49 | 49 | 49 |
| TDUP | -1 | 50 | 50 | 50 | 50 |
| TDUP | 1 | 50 | 50 | 50 | 50 |

Builder statuses: unsupported (duplication_adjacent_repeat) 1; unsupported (inversion_not_reducible) 3; supported 396. VEP's default 3' shift moves 78 of the 396 supported genomic strings, retained in [hgvs_genomic.csv](sv_species/plant_tair10_seed11663/hgvs_genomic.csv) (`shift1_differs`).

`duckvep_prepare_sv_geometry_sql` nominal start/end equal VEP's `start`/`end` for all 1157 symbolic records (0 disagreements). `duckvep_prepare_breakend_pairs_sql` on the BND records: reciprocal 576.


## The pooled-BND counterexamples

In the mouse and Arabidopsis pooled BND runs, three mate records (`bnd186b` in mouse; `bnd130b` and `bnd230b` in Arabidopsis) had transcripts that
DuckVEP reports but VEP did not: 28, 5 and 4 `duckvep_only` pairs. All three are the generated mate (`b`) record of a pair whose `a` record is in the same file; most other `b` records agree. Re-run alone (`isolated_rerun.csv`), VEP reports exactly the DuckVEP transcript sets
(29/29, 11/11 and 12/12 pairs exact). The mouse case: in the pooled run VEP returned one chr1 transcript for `bnd186b` (chr2:113,877,763) and 29 for its partner `bnd186a`, 28 of them
on chr2 at `bnd186b`'s own locus; run alone, `bnd186b` returns those chr2 transcripts with `3_prime_UTR_variant&feature_truncation`, as DuckVEP
does. This is the record-dependence that [ERRATA.md](../../ERRATA.md) records for VEP's chromosome-blind BND input-buffer interval tree. The
observation here is that `--buffer_size 1` on a multi-record file did not remove it when a record's mate coordinate is another record's own
position; the mechanism beyond that is not established by this differential. DuckVEP annotates each physical record independently. The claim made
is therefore *record-isolated* BND agreement: the 3 disagreeing records are exact against isolated VEP runs, and the rest are exact in the pooled runs.

## Structural HGVS scope

VEP 116 emits no HGVS for symbolic alleles or breakends, so the domain is the same as in
[`structural_hgvs_vep116_differential.R`](../../test/duckvep/conformance/structural_hgvs_vep116_differential.R): the builder's
literal-equivalent edit is annotated by VEP and its genomic string compared with the builder's `hgvs_g` under `--shift_hgvs 0` (the builder
declares `normalization = 'none'`), and the transcript strings are compared against DuckVEP `hgvs` on the same edit. VEP's default 3' shifting
changes 73, 81 and 78 genomic strings (mouse, fly, Arabidopsis): these are retained counterexamples to any normalized-HGVS claim, not evidence
against the unshifted claim. A few events per species are refused by the builder (`inversion_not_reducible`, plus one `duplication_adjacent_repeat`
in Arabidopsis); the random corpus is too small in repeat-adjacent duplications to replace the targeted human control, so refusal behavior is
not re-validated per species.

## Models, pins and provenance

- Extension: `duckvep.duckdb_extension` SHA-256 `bbc690cc3fac46f1651a8c6854030fe8bf9e938572eb6288a153dca6970cedf3`, built from
  `7a741c0` with an INT64 overflow guard in `duckvep_prepare_sv_geometry_sql` and hash-join mate lookup in
  `duckvep_prepare_breakend_pairs_sql`. These fixes do not affect the results in these tables. Every run loads an immutable copy.
- Mouse: the existing pinned model (`model_sha256 ae39ffc9...`, receipt in `duckvep_model_receipts.csv`).
- Fly and Arabidopsis models were built from the pinned sources with `scripts/build_species_model.R`. The staged `core.duckdb` bytes differ
  from the recorded snapshot digest (DuckDB files are not byte-reproducible), and `model_sha256` covers `source_manifest_sha256`, so the raw
  hashes differ (`d84f5392...` for fly, `521c53d0...` for Arabidopsis). Recomputing `duckvep_model_receipt_sql` over the rebuilt tables with the
  *pinned* manifest digest reproduces the recorded hashes `e6deea1a...` (fly) and `5b3fac7c...` (Arabidopsis) exactly, with matching
  counts, so the row content is identical. The FASTA, FASTA index, and indexed-cache archive hashes match the source pins. Staging the three inputs and caches used about 11 GB.
