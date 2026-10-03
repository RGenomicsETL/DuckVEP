
<!-- README.md is generated from README.Rmd using duckknit. Please edit that file. -->

# DuckVEP

[![Extension CI](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/MainDistributionPipeline.yml/badge.svg)](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/MainDistributionPipeline.yml)
[![Site](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/pages.yml/badge.svg)](https://rgenomicsetl.github.io/DuckVEP/)

**Variant effect prediction as a SQL relation: load an Ensembl release into DuckDB once, then ask your variants questions interactively.**

> Variant annotation take time, too much time, variant_myth try to solve this problem using modern and parallel tools.
>
> — [variant_myth](https://github.com/natir/variant_myth)

That sentence is where this project began. variant_myth answers it with a fast, parallel annotator: it reads GFF3, FASTA and VCF and writes Parquet, one row per variant and transcript. DuckVEP takes the idea one step further. If annotation ends up as a table anyway, it can be a table you query directly.

Load a transcript model into DuckDB once. Annotation then becomes a relational operator, `FROM query(duckvep_annotate_sql(...))`, that you filter, join and aggregate in the same statement as everything else. When the question changes, you change the query and run it again in seconds; you don’t rerun a pipeline. Variant effect prediction stops being a batch job you wait for, and becomes something you explore.

What comes back are typed rows, not a text file to parse:
- Sequence Ontology consequences and impact;
- cDNA, CDS and protein positions and amino acids;
- NMD prediction;
- HGVS.

Those rows join with anything DuckDB reads: a VCF through [DuckHTS](https://github.com/RGenomicsETL/duckhts), Parquet, DuckLake, ClinVar, population frequencies. There is no plugin system and no new file format.

Speed is only useful if the answers are right. So DuckVEP’s reference is executable Ensembl VEP itself, pinned to one release (today VEP 116 with its matching core and variation databases). Every claim below is measured against it. Every result that differs stays in the count, and every deliberate difference is written down. Moving to a new VEP release means re-running that ledger, not re-arguing the rules.

## Try it

A model is loaded from ordinary relations: regions, transcripts with their CDS and flanking sequence, and exons. In production these come from an Ensembl core database through `query(duckvep_ensembl_regions_sql(...))` and `query(duckvep_ensembl_transcripts_sql(...))`. Here a [one-transcript fixture](test/data/duckvep/readme.sql) stands in:

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

The model now stays in the session. Variants are just rows. Here are seven across the transcript, from the 5′ UTR to the intron, annotated and joined back to their input in one query:

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

Some HGVS fields are empty, and that isn’t a parsing failure. This model was loaded without a reference FASTA, and those notations need flanking sequence to shift and render. DuckVEP says so rather than guess, and asking why is just another query:

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

Phased genotypes become whole-haplotype consequences, again as a relation. The changes at 124 and 125 fall in the same codon (Val2). On their own, one is missense and the other synonymous. Carried on the same haplotype, they are read together as one codon, GTA to GCG, and classified once:

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

To build that calls relation from a whole-genome VCF, `duckvep_coding_transcripts(model, seq_region, position, reference, alternate)` looks each record up in the model’s interval index. It returns the transcripts whose coding sequence the record touches, using the same normalization as the annotation builder, so an indel’s anchor base alone never counts. Records outside coding sequence return an empty list, so an `unnest()` over it reduces a genome to its few coding records before any per-record work. In R it is `Rduckvep::rduckvep_coding_transcripts()`.

## How close to VEP

Transcript and feature consequence pairs compared with executable VEP, newest tested run per corpus:

| Corpus                                 | Assembly        | Oracle  |   Pairs |   Exact | Unresolved | Different |
|:---------------------------------------|:----------------|:--------|--------:|--------:|-----------:|----------:|
| dbSNP 157 windows                      | GRCh38          | VEP 116 |  73,620 |  73,620 |          0 |         0 |
| GIAB HG002 small variants              | GRCh38          | VEP 116 |  54,905 |  54,905 |          0 |         0 |
| ClinVar coding                         | GRCh38          | VEP 116 | 287,836 | 287,836 |          0 |         0 |
| ClinVar, all chromosomes               | GRCh38          | VEP 116 | 316,397 | 316,397 |          0 |         0 |
| GRCh37 cache corpus                    | GRCh37          | VEP 116 | 486,464 | 486,464 |          0 |         0 |
| *P. falciparum* (Ensembl Genomes 63)   | GCA_000002765v3 | VEP 116 |  40,732 |  40,732 |          0 |         0 |
| Paired breakends, multi-chromosome     | GRCh38          | VEP 116 |  91,428 |  91,428 |          0 |         0 |
| GIAB + regulatory and motif features   | GRCh38          | VEP 116 |  14,955 |  14,955 |          0 |         0 |
| Exact structural variants + regulation | GRCh38          | VEP 116 | 120,224 | 120,224 |          0 |         0 |

That is **1,486,561 of 1,486,561** pairs identical to VEP. They cover:
- human GRCh38 and GRCh37, and *Plasmodium falciparum* with its own genetic codes;
- small variants, exact structural events and paired breakends;
- regulatory and motif features.

The release gate doesn’t bargain: one discordant, missing, extra or unresolved pair fails the run.

HGVS is held to the same standard: 357,806 HGVSc and 357,806 HGVSp transcript pairs match VEP `--hgvs` string for string (56,998 of them from ClinVar chromosome 21), with 0 discordant.

The [conformance report](benchmarks/duckvep_conformance.md) has:
- the methods;
- the per-term and per-impact tables;
- the generated state-exploration campaigns;
- the statistical caveats.

Ledger revisions from before the extraction from DuckHTS resolve through the [commit map](design/duckhts-commit-map.txt).

## Where it differs, and why

Compatibility is claimed only where it has been measured, and differences are kept, not hidden. The [compatibility and errata record](ERRATA.md) separates three verdicts, each with its witness:
- a VEP convention DuckVEP must follow;
- a gap in DuckVEP;
- a likely upstream erratum.

The main cases:

- **Published release annotations are not the oracle.** Ensembl’s release VCFs sometimes disagree with executable VEP. In release 116, `X/Y:276322 G>A` is published as `intergenic_variant`, while the executable emits three `5_prime_UTR_variant` rows per chromosome. DuckVEP follows the executable, and the [PAR witnesses](test/duckvep/conformance/README.md) pin both.
- **Breakends are evaluated one event at a time.** VEP 116’s buffered breakend path uses a chromosome-blind interval tree and can drop valid transcript pairs in multi-chromosome batches. The oracle runs with `--buffer_size 1` to isolate that; DuckVEP has no such batching effect ([paired-breakend differential](benchmarks/duckvep_conformance.md#paired-breakend-differential)).
- **Imprecise structural variants use the nominal span**, as VEP’s registered predicates do. `CIPOS` and `CIEND` stay on the row as metadata rather than being dropped.
- **No silent size limits.** VEP skips structural events above `--max_sv_size` (5 kb by default). DuckVEP annotates every exact span, and the oracle runs with a 10 Mb limit so the comparison covers them.
- **Circular regions run on lifted intervals.**
  - A circular sequence region with an object that crosses its origin runs on a lifted linear copy of the model. Such an object can be a transcript, exon, regulatory or motif feature. On the lifted copy, flanks, splice sites, CDS positions, HGVS shifts and regulatory overlaps all cross the origin consistently.
  - This is separate from the mitochondrial codon table, which only chooses a translation rule. Human MT has no wrapped object, so it keeps VEP’s linear behaviour byte for byte.
  - VEP models origin-crossing transcripts as reversed intervals, so it can’t be the oracle for them. Those are checked instead by rotation equivariance and against a linear model ([design](design/duckvep.md), [survey and differential](benchmarks/data/circular_source_survey.md)).
- **Unknown is a value.** When a result can’t be computed, the row carries `duckvep_status` and a reason instead of a guess. Examples are a missing reference sequence, a reference mismatch, or an unsupported symbolic allele.
- **Whole-haplotype consequences are checked against `bcftools csq`, not VEP.**
  - VEP 116 does not define them. The coding-only contract is [`duckvep-coding-v1`](design/duckvep_haplotype_contract.md), with independent base-R goldens and a declared list of every divergence from csq.
  - DuckVEP is 1.56× faster than csq end to end on HG002 (the 2× goal is tracked in [\#34](https://github.com/RGenomicsETL/DuckVEP/issues/34)).
- **Outside the claim:**
  - species and releases that haven’t been tested;
  - compound (cis) HGVS presentation, where the errata record the disagreements that remain;
  - genomic (`g.`) and structural HGVS.

  See the [design contract](design/duckvep.md) and the [roadmap](https://github.com/RGenomicsETL/DuckVEP/issues).

## Speed

Interactive only works if the engine keeps up.
- **The public SQL builder** annotates **1,034,768 alleles per second on one core** in compact output. That is every model-addressable allele of GIAB HG002 (4,095,611), with all 1,383,580 regulatory and motif features resident, straight through `FROM query(duckvep_annotate_sql(...))` ([throughput report](benchmarks/duckvep_throughput.md)).
- **End to end against FastVEP**, VCF in and table out, including DuckVEP’s own coordinate sort: DuckVEP finished **2.55×** faster on one core (65 s against 164 s) and **2.12×** faster on four. That is for a compact output contract on the full GIAB HG002 GRCh38 benchmark (4,096,123 ALT alleles). Complete-field comparisons, CSQ with HGVS, and every caveat are in the [FastVEP benchmark](benchmarks/benchmark_duckvep_fastvep.md).
- **Memory** is bounded, not hoped for. Each job’s native memory is charged against an enforced budget, and exceeding it is an explicit capacity error, never a truncated result. Three concurrent 5-million-allele gnomAD jobs are certified on one 20-thread host ([scale runner](docs/scale-runner.md)).

## From R

`Rduckvep` runs the same native builders on your DuckDB connection. Builder functions return SQL, and `rduckvep_annotate()` executes an annotation and returns a data frame. The relations you pass stay visible in that connection, including TEMP tables and uncommitted rows.

```r
con <- Rduckvep::rduckvep_connect()
# After loading a model named "demo" and creating "demo_events":
sql <- Rduckvep::rduckvep_annotate_sql(con, "demo_events", "demo", hgvs = TRUE)
results <- Rduckvep::rduckvep_annotate(con, "demo_events", "demo", hgvs = TRUE)
DBI::dbDisconnect(con, shutdown = TRUE)
```

## With the rest of the stack

DuckVEP needs no other extension, but it is meant to sit in a query next to them. With [DuckHTS](https://github.com/RGenomicsETL/duckhts) reading the VCF and [DuckClinVarbitration](https://github.com/RGenomicsETL/DuckClinVarbitration) supplying arbitrated ClinVar decisions, annotation and clinical evidence meet in one query. This one isn’t evaluated here, because it needs a full Ensembl model and both extensions:

```sql
-- One row per ALT allele, with the model's region index (DuckHTS reads the VCF)
CREATE TABLE alleles AS
SELECT row_number() OVER (ORDER BY r.seq_region, v.POS, alt_allele)::UBIGINT AS event_index,
       r.seq_region, v.POS::UBIGINT AS position, v.REF AS reference, alt_allele AS alternate,
       NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM read_bcf('cohort.vcf.gz') v CROSS JOIN unnest(v.ALT) AS t(alt_allele)
JOIN grch38_regions r   -- Ensembl names contigs 1..22, X, Y, MT; UCSC-style VCFs say chr1, chrM
  ON r.seq_region_name = CASE WHEN v.CHROM IN ('chrM', 'M') THEN 'MT'
                              ELSE regexp_replace(v.CHROM, '^chr', '') END;

-- High-impact consequences next to their arbitrated ClinVar decision
SELECT r.seq_region_name AS contig, e.position, e.reference, e.alternate,
       a.consequence, a.protein_hgvs, d.policy_classification, d.gold_stars
FROM query(duckvep_annotate_sql('alleles', 'human_116_grch38', struct_pack(rich := true, hgvs := true))) a
JOIN alleles e USING (event_index)
JOIN grch38_regions r USING (seq_region)
LEFT JOIN clinvar_vcf v                              -- DuckClinVarbitration
  ON v.assembly = 'GRCh38' AND v.contig = r.seq_region_name AND v.position = e.position
 AND v.reference = e.reference AND v.alternate = e.alternate
LEFT JOIN clinvar_policy_allele_decisions d USING (allele_id)
WHERE a.impact = 'HIGH';
```

## Install and build

```sh
git clone --recurse-submodules https://github.com/RGenomicsETL/DuckVEP
cd DuckVEP
make configure release test_release
duckdb -unsigned -c "LOAD 'build/release/duckvep.duckdb_extension'"
```

DuckHTS 1.5.2 and earlier still ship the DuckVEP functions they once bundled, under the same SQL names. With such a DuckHTS loaded in the same database, models load into one extension and annotation looks in the other, and queries fail with `unknown model name`. Use a DuckHTS release newer than 1.5.2 next to DuckVEP.

The extension uses the stable DuckDB C API (tested on DuckDB 1.5 and the 2.0 pre-release) and links its own htslib (for indexed reference FASTA) and cgranges. From R, [Rduckvep](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/) builds the same sources offline and adds connection, model and haplotype helpers.

## Documentation

- [Function reference](docs/functions.md): every public SQL function, its signature, options, columns and a checked example
- [Site](https://rgenomicsetl.github.io/DuckVEP/) and [R package reference](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/)
- [Compatibility and errata](ERRATA.md): every known difference from VEP, classified, with witnesses
- [Design and implementation contract](design/duckvep.md): model build, the sweep, NMD, phased edits, structural events, HGVS, supplementary annotation
- [Haplotype contract `duckvep-coding-v1`](design/duckvep_haplotype_contract.md) and its [scale qualification](benchmarks/data/haplotype_scale/README.md)
- Reports: [conformance](benchmarks/duckvep_conformance.md), [throughput](benchmarks/duckvep_throughput.md), [FastVEP comparison](benchmarks/benchmark_duckvep_fastvep.md), [haplotypes](benchmarks/duckvep_haplotypes.md), [scale runner](docs/scale-runner.md)
- [Corpus workflow](design/duckvep_corpus_workflow.md) and [conformance harness](test/duckvep/conformance/README.md)

DuckVEP grew up inside [DuckHTS](https://github.com/RGenomicsETL/duckhts) and was extracted with its history.

## License

GPL-2.0-or-later; see [LICENSE](LICENSE). Vendored htslib and cgranges keep their own licences. Ensembl data and VEP are © EMBL-EBI under their own terms.
