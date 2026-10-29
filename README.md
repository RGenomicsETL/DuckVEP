
<!-- README.md is generated from README.Rmd using duckknit. Please edit that file. -->

# DuckVEP

[![Extension CI](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/MainDistributionPipeline.yml/badge.svg)](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/MainDistributionPipeline.yml)
[![Site](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/pages.yml/badge.svg)](https://rgenomicsetl.github.io/DuckVEP/)

**Variant effect prediction inside DuckDB.** Load an Ensembl transcript model once; from then on annotation is a relation you filter, join and aggregate like any table, next to your VCFs, ClinVar or population frequencies. Rows carry Sequence Ontology consequences and impact, cDNA/CDS/protein positions, amino acids, NMD predictions and HGVS. Phased calls give whole-haplotype consequences. Everything is checked against executable Ensembl VEP.

[Install](#install) · [Try it](#try-it) · [Function reference](docs/functions.md) · [Conformance](benchmarks/duckvep_conformance.md) · [Errata](ERRATA.md) · [Design](design/duckvep.md) · [R package](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/)

## Install

DuckVEP is a C extension on DuckDB’s stable C API v1. CI builds and tests it against DuckDB 1.5.6, and the same build loads in DuckDB 2.0.0. Build it from source and load it unsigned:

```sh
git clone --recurse-submodules https://github.com/RGenomicsETL/DuckVEP
cd DuckVEP
make configure release test_release
duckdb -unsigned -c "LOAD 'build/release/duckvep.duckdb_extension'"
```

`LOAD` registers functions and nothing else: it downloads nothing and creates no tables. Transcript models come from Ensembl data you supply, or from a snapshot you saved earlier. Prebuilt Ensembl GRCh38 and GRCh37 models are published in the [duckvep-models](https://huggingface.co/datasets/RgenomicsETL/duckvep-models) dataset on Hugging Face; download one, or attach it in place with `httpfs` as `hf://datasets/RgenomicsETL/duckvep-models/...`. The extension links its own htslib (for indexed FASTA) and cgranges. From R, [Rduckvep](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/) builds the same sources offline and adds connection, model and haplotype helpers.

A second host built on DuckDB’s C API v2 is available as a preview ([docs/v2-host.md](docs/v2-host.md)); making it the default is [\#48](https://github.com/RGenomicsETL/DuckVEP/issues/48).

## Try it

A model is loaded from three relations: regions, transcripts (with CDS and flanking sequence) and exons. In practice they come from an Ensembl core database through `query(duckvep_ensembl_regions_sql(...))` and `query(duckvep_ensembl_transcripts_sql(...))`; here a [one-transcript fixture](test/data/duckvep/readme.sql) stands in:

```sql
SELECT loaded FROM duckvep_model_load('demo',
  'SELECT * FROM readme_regions ORDER BY seq_region',
  'SELECT * FROM readme_transcripts ORDER BY seq_region, transcript_start',
  'SELECT * FROM readme_exons ORDER BY transcript_index, exon_start');
#> ┌─────────┐
#> │ loaded  │
#> │ boolean │
#> ├─────────┤
#> │ true    │
#> └─────────┘
```

The model stays resident in the session. Annotate seven variants, from the 5′ UTR to the intron, and join the results back to the input:

```sql
CREATE TABLE demo_events AS
SELECT row_number() OVER (ORDER BY position, alternate)::UBIGINT AS event_index,
       1::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate,
       NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM (VALUES (110, 'C', 'T'), (124, 'T', 'C'), (125, 'A', 'G'), (130, 'C', 'CA'),
             (134, 'C', 'A'), (151, 'G', 'A'), (175, 'A', 'G')) v(position, reference, alternate);
SELECT e.position, e.reference || '>' || e.alternate AS change,
       a.consequence, a.impact, a.transcript_hgvs, a.protein_hgvs, a.nmd_prediction
FROM query(duckvep_annotate_sql('demo_events', 'demo', struct_pack(rich := true, hgvs := true))) a
JOIN demo_events e USING (event_index)
ORDER BY e.position;
#>
#> ┌──────────┬─────────┬──────────────────────┬──────────┬─────────────────┬──────────────┬────────────────┐
#> │ position │ change  │     consequence      │  impact  │ transcript_hgvs │ protein_hgvs │ nmd_prediction │
#> │  uint64  │ varchar │       varchar        │ varchar  │     varchar     │   varchar    │    varchar     │
#> ├──────────┼─────────┼──────────────────────┼──────────┼─────────────────┼──────────────┼────────────────┤
#> │      110 │ C>T     │ 5_prime_UTR_variant  │ MODIFIER │ NULL            │ NULL         │ NULL           │
#> │      124 │ T>C     │ missense_variant     │ MODERATE │ c.5T>C          │ p.Val2Ala    │ NULL           │
#> │      125 │ A>G     │ synonymous_variant   │ LOW      │ c.6A>G          │ p.Val2=      │ NULL           │
#> │      130 │ C>CA    │ frameshift_variant   │ HIGH     │ NULL            │ NULL         │ escaping       │
#> │      134 │ C>A     │ stop_gained          │ HIGH     │ c.15C>A         │ p.Tyr5Ter    │ escaping       │
#> │      151 │ G>A     │ splice_donor_variant │ HIGH     │ NULL            │ NULL         │ unresolved     │
#> │      175 │ A>G     │ intron_variant       │ MODIFIER │ NULL            │ NULL         │ NULL           │
#> └──────────┴─────────┴──────────────────────┴──────────┴─────────────────┴──────────────┴────────────────┘
```

The frameshift and the stop gain escape nonsense-mediated decay because they fall in the first CDS positions (`nmd_escape_early_cds`).

The empty HGVS fields are deliberate: this model was loaded without a reference FASTA, and those notations need flanking sequence to shift and render. Every row says why:

```sql
SELECT e.position, a.transcript_hgvs_status, a.transcript_hgvs_reason
FROM query(duckvep_annotate_sql('demo_events', 'demo', struct_pack(hgvs := true))) a
JOIN demo_events e USING (event_index)
WHERE a.transcript_hgvs IS NULL
ORDER BY e.position;
#> ┌──────────┬────────────────────────┬────────────────────────┐
#> │ position │ transcript_hgvs_status │ transcript_hgvs_reason │
#> │  uint64  │        varchar         │        varchar         │
#> ├──────────┼────────────────────────┼────────────────────────┤
#> │      110 │ unresolved             │ missing_reference      │
#> │      130 │ unresolved             │ missing_reference      │
#> │      151 │ unresolved             │ missing_reference      │
#> │      175 │ unresolved             │ missing_reference      │
#> └──────────┴────────────────────────┴────────────────────────┘
```

Phased genotypes give whole-haplotype consequences. The changes at 124 and 125 sit in the same codon (Val2): alone, one is missense and the other synonymous; on the same haplotype they make one codon, GTA to GCG, classified once:

```sql
CREATE TABLE demo_calls AS
SELECT event_index, seq_region, position, reference, alternate,
       1::UINTEGER alt_index, 0::UINTEGER transcript_index, 0::UBIGINT sample_index,
       [1, 0]::INTEGER[] alleles, [false, true]::BOOLEAN[] phase_before, NULL::BIGINT phase_set
FROM demo_events WHERE position IN (124, 125);
SELECT carrier_count, length(contributors) AS contributors, prediction_status,
       haplotype_consequences, haplotype_impact
FROM duckvep_haplotypes('SELECT * FROM demo_calls', 'demo');
#>
#> ┌───────────────┬──────────────┬───────────────────┬────────────────────────┬──────────────────┐
#> │ carrier_count │ contributors │ prediction_status │ haplotype_consequences │ haplotype_impact │
#> │    uint32     │    int64     │      varchar      │       varchar[]        │     varchar      │
#> ├───────────────┼──────────────┼───────────────────┼────────────────────────┼──────────────────┤
#> │             1 │            2 │ predicted         │ [missense_variant]     │ MODERATE         │
#> └───────────────┴──────────────┴───────────────────┴────────────────────────┴──────────────────┘
```

For VCF or BCF input, [`duckvep_coding_calls(model, path)`](docs/functions.md#duckvep_coding_calls) builds that calls relation directly: it maps `CHROM` to the model’s `seq_region_name`, finds the coding transcripts each record touches and decodes GT and PS only for those.

## How close to VEP

Across 9 recorded corpora, executable Ensembl VEP and DuckVEP agree on **1,486,561 of 1,486,561 transcript and feature pairs**; 0 are unresolved and 0 differ. The corpora cover GRCh37 and GRCh38, *Plasmodium falciparum* genetic codes, small variants, exact structural events, paired breakends, and regulatory and motif features. HGVS against VEP `--hgvs`: 357,806 of 357,806 HGVSc and 357,806 of 357,806 HGVSp pairs agree, 56,998 of them from ClinVar chromosome 21. The newest of these runs is dated 2026-09-07; changes merged after it are not covered until the corpora are rerun.

The [conformance report](benchmarks/duckvep_conformance.md) has the methods and term-level results; the [consequence](test/duckvep/conformance/data/conformance_history.csv) and [HGVS](test/duckvep/conformance/data/hgvs_history.csv) ledgers hold every run with its source revision.

<details>
<summary>
Per-corpus counts
</summary>

| Corpus                                 | Assembly        | Oracle  |   Pairs |   Exact | Unresolved | Different | Run date   |
|:---------------------------------------|:----------------|:--------|--------:|--------:|-----------:|----------:|:-----------|
| dbSNP 157 windows                      | GRCh38          | VEP 116 |  73,620 |  73,620 |          0 |         0 | 2026-07-14 |
| GIAB HG002 small variants              | GRCh38          | VEP 116 |  54,905 |  54,905 |          0 |         0 | 2026-07-14 |
| ClinVar coding                         | GRCh38          | VEP 116 | 287,836 | 287,836 |          0 |         0 | 2026-07-16 |
| ClinVar, all chromosomes               | GRCh38          | VEP 116 | 316,397 | 316,397 |          0 |         0 | 2026-07-16 |
| GRCh37 cache corpus                    | GRCh37          | VEP 116 | 486,464 | 486,464 |          0 |         0 | 2026-07-15 |
| *P. falciparum* (Ensembl Genomes 63)   | GCA_000002765v3 | VEP 116 |  40,732 |  40,732 |          0 |         0 | 2026-07-15 |
| Paired breakends, multi-chromosome     | GRCh38          | VEP 116 |  91,428 |  91,428 |          0 |         0 | 2026-07-16 |
| GIAB + regulatory and motif features   | GRCh38          | VEP 116 |  14,955 |  14,955 |          0 |         0 | 2026-07-16 |
| Exact structural variants + regulation | GRCh38          | VEP 116 | 120,224 | 120,224 |          0 |         0 | 2026-07-16 |

</details>

## Whole haplotypes

VEP annotates each record on its own. DuckVEP also rebuilds each haplotype’s transcript and classifies the whole product, following how VEP reads the same data wherever VEP defines an answer:

- **Phase.** `strict` uses GT and PS as written. `vep_compat` reads genotype slots in order, as VEP’s haplotype mode does.
- **Overlapping records** apply at their own coordinates. Records that agree compose once; two that claim one base or insertion site differently are an `edit_conflict`.
- **The product.** The edited CDS is translated, read through the 3′ UTR when the stop is lost, and classified against the reference protein. Curated translations (selenocysteine, initial Met, readthrough) and RNA editing programmes are part of the reference, as in Ensembl. NMD follows the 50-base exon-junction rule.
- **Splicing is the model’s**, as in VEP. A record that touches no exon leaves the product alone.
- **Per-record terms.** Some paths have no whole product to rebuild: an edit across a splice site or a CDS end, or an in-frame stop in the reference CDS. Those rows carry each record’s VEP terms, with reason `per_record_terms`.

On GIAB HG002 with the Ensembl 116 model and strict phase, 157,597 paths get a whole-product prediction and 212 carry per-record terms. The other 328 are `incomplete_input`: several unphased heterozygous sites on one transcript, which the calls leave undefined. [`duckvep_haplotype_arrangements`](docs/functions.md#duckvep_haplotype_arrangements) enumerates their alternatives, and `vep_compat` reads them in slot order.

[`bcftools csq`](https://github.com/samtools/bcftools) is the executable check where its model and phase rules are comparable. Single-record paths are also compared with DuckVEP’s per-variant annotation, which is exact against VEP; the [contract](design/duckvep_haplotype_contract.md) lists the rules, the evidence and the named differences.

## Speed

Measured on one host (Intel i5-13500, cores pinned) with the Ensembl 116 GRCh38 model; the receipts below record the inputs, revisions and every run.

| Workload                                                        | DuckVEP                                       | Comparison                                                                  |
|-----------------------------------------------------------------|-----------------------------------------------|-----------------------------------------------------------------------------|
| Annotate GIAB HG002 (4,095,611 alleles), consequence codes only | 1,090,418 alleles/s on 1 core; 4,031,113 on 4 |                                                                             |
| Same, with rich columns and HGVS                                | 174,772 alleles/s on 1 core; 626,528 on 4     |                                                                             |
| HG002 end to end on 1 core: read VCF, annotate, write output    | 35.0 s, 2.7 GiB peak                          | [vep-rs](https://github.com/natera-open-source/vep-rs) 41.9 s, 2.8 GiB peak |
| HG002 whole haplotypes on 1 core, model load included           | 8.7 s from the CLI, 9.1 s from R              | `bcftools csq` 17.1 s                                                       |
| Load the model                                                  | 0.55 s from a snapshot (1.3 GiB)              | 2.33 s from relations (1.48 s on 8 threads)                                 |

On the transcripts both tools annotate, DuckVEP and vep-rs agree on 34,146,531 of 34,148,222 consequence tuples; executable VEP sides with DuckVEP on 1,689 of the 1,691 disagreements. Snapshots are memory-mapped, so processes restoring the same file share its pages. Native memory is charged against an enforced budget: exceeding it is an error, never a truncated result.

Receipts: [annotation throughput](benchmarks/duckvep_throughput.md), [vep-rs comparison](benchmarks/benchmark_duckvep_vep_rs.md), [haplotypes against csq](benchmarks/data/haplotype_scale/hg002_2417ac9/), [model load](benchmarks/data/duckvep_model_load.csv), [multi-job scale runs](docs/scale-runner.md).

## Where it differs from VEP

Differences are recorded with a witness and a verdict (a VEP convention DuckVEP follows, a gap in DuckVEP, or a likely upstream erratum) in the [errata](ERRATA.md). The main ones:

- **Executable VEP is the oracle, not Ensembl’s published release files.** In release 116, `X/Y:276322 G>A` is published as `intergenic_variant`, while executable VEP emits three `5_prime_UTR_variant` rows per chromosome; DuckVEP follows the executable.
- **Breakends are evaluated one event at a time.** VEP 116’s buffered breakend path can drop valid transcript pairs in multi-chromosome batches, so the oracle runs with `--buffer_size 1`.
- **Structural events** use their nominal span, with no `--max_sv_size` cutoff; `CIPOS`/`CIEND` stay on the row as metadata.
- **Circular regions** use lifted intervals, so origin-crossing transcripts and features keep their overlaps and HGVS. VEP cannot check these; they are tested by rotation and against a linearized model.
- **Unknown is a value.** When a result can’t be computed (no reference sequence, a reference mismatch, an unsupported symbolic allele) the row says so with a status and reason.

Untested species and releases, and genomic or structural HGVS beyond the [exact-span profile](docs/structural-identity-hgvs.md), are outside the claim.

## From R

`Rduckvep` runs the same builders on your DuckDB connection: builder functions return SQL, and `rduckvep_annotate()` runs it and returns a data frame. Queries can read anything on that connection, TEMP tables included.

```r
con <- Rduckvep::rduckvep_connect()
# After loading a model named "demo" and creating "demo_events":
results <- Rduckvep::rduckvep_annotate(con, "demo_events", "demo", hgvs = TRUE)
DBI::dbDisconnect(con, shutdown = TRUE)
```

## With the rest of the stack

DuckVEP needs no other extension, but it is meant to share a query with them: [DuckHTS](https://github.com/RGenomicsETL/duckhts) reads the VCF and [DuckClinVarbitration](https://github.com/RGenomicsETL/DuckClinVarbitration) supplies arbitrated ClinVar decisions. Not evaluated here, since it needs a full Ensembl model and both extensions:

```sql
-- The model's regions under the VCF's own contig names: DuckHTS reads the header, and
-- duckhts_contig_key maps each name (chr1, chrM) to Ensembl's (1, MT). Breakend mates
-- written in ALT (A[chr2:200[) resolve through the same rows.
CREATE VIEW cohort_regions AS
SELECT seq_region, seq_region_name FROM grch38_regions
UNION
SELECT r.seq_region, h.id FROM read_hts_header('cohort.vcf.gz') h
JOIN grch38_regions r ON r.seq_region_name = duckhts_contig_key(h.id)
WHERE h.record_type = 'contig';

-- One row per ALT allele, read as VEP reads it (DuckHTS reads the VCF)
CREATE VIEW cohort AS SELECT * FROM read_bcf('cohort.vcf.gz');
CREATE TABLE alleles AS
SELECT * FROM query(duckvep_vcf_events_sql('cohort', 'cohort_regions'));
CREATE VIEW annotatable AS SELECT * FROM alleles WHERE event_status = 'ok' ORDER BY event_index;

-- High-impact consequences next to their arbitrated ClinVar decision
SELECT r.seq_region_name AS contig, e.source_pos, e.source_ref, e.source_alt,
       a.consequence, a.protein_hgvs, d.policy_classification, d.gold_stars
FROM query(duckvep_annotate_sql('annotatable', 'human_116_grch38', struct_pack(rich := true, hgvs := true))) a
JOIN annotatable e USING (event_index)
JOIN grch38_regions r USING (seq_region)
LEFT JOIN clinvar_vcf v                              -- DuckClinVarbitration
  ON v.assembly = 'GRCh38' AND v.contig = r.seq_region_name AND v.position = e.source_pos
 AND v.reference = e.source_ref AND v.alternate = e.source_alt
LEFT JOIN clinvar_policy_allele_decisions d USING (allele_id)
WHERE a.impact = 'HIGH';
```

`event_reason` says why a record was left out, for example a gVCF `A,<NON_REF>` record that VEP skips, or a contig the model does not have. The structural preparation builders (breakend identity, confidence bounds, structural HGVS, ExpansionHunter repeats) read the same `read_bcf` rows, with their typed `INFO_*` columns, once each row has an `event_index`.

## Documentation

- **Reference:** [SQL functions](docs/functions.md) with signatures, options, columns and checked examples; [R package](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/); [C API hosts](docs/v2-host.md)
- **Design:** [implementation](design/duckvep.md); [haplotype contract](design/duckvep_haplotype_contract.md); [structural identity and HGVS](docs/structural-identity-hgvs.md); [corpus workflow](design/duckvep_corpus_workflow.md)
- **Evidence:** [conformance](benchmarks/duckvep_conformance.md) and its [harness](test/duckvep/conformance/README.md); [errata](ERRATA.md); [haplotypes](benchmarks/duckvep_haplotypes.md)
- **Performance:** [throughput](benchmarks/duckvep_throughput.md), [vep-rs comparison](benchmarks/benchmark_duckvep_vep_rs.md), [scale runner](docs/scale-runner.md)

The [project site](https://rgenomicsetl.github.io/DuckVEP/) publishes these pages.

## Origins

> Variant annotation take time, too much time, variant_myth try to solve this problem using modern and parallel tools.
>
> — [variant_myth](https://github.com/natir/variant_myth)

DuckVEP started from that sentence and asks the same question inside the database; thank you to variant_myth’s author for asking it. Correctness rests on two other projects whose executables make a runnable oracle possible: [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) and, for whole haplotypes, [`bcftools csq`](https://github.com/samtools/bcftools). DuckVEP was developed inside [DuckHTS](https://github.com/RGenomicsETL/duckhts) and extracted with its Git history.

## License

GPL-2.0-or-later; see [LICENSE](LICENSE). Vendored htslib and cgranges keep their own licences. Ensembl data and VEP are © EMBL-EBI under their own terms.
