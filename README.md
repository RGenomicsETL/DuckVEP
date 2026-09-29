
<!-- README.md is generated from README.Rmd using duckknit. Please edit that file. -->

# DuckVEP

[![Extension CI](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/MainDistributionPipeline.yml/badge.svg)](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/MainDistributionPipeline.yml)
[![Site](https://github.com/RGenomicsETL/DuckVEP/actions/workflows/pages.yml/badge.svg)](https://rgenomicsetl.github.io/DuckVEP/)

**A variant consequence engine for DuckDB, with Ensembl VEP semantics and a ledger of every place it differs.**

Variant effect prediction usually means a Perl pipeline that reads a VCF and writes a text file. DuckVEP makes it a relational operator. An Ensembl release is compiled once into an immutable, shared transcript model; sorted alleles then stream through a native C sweep that emits typed rows: Sequence Ontology consequences, impact, cDNA/CDS/protein positions, amino acids, NMD prediction and HGVS. Those rows join with anything DuckDB can read (a VCF through [DuckHTS](https://github.com/RGenomicsETL/duckhts), Parquet, DuckLake, ClinVar, population frequencies) without a plugin system or a new file format.

The semantics are not a reinterpretation of the SO rules. They are defined by executable Ensembl VEP: every claim is measured against a pinned release of the real tool (today VEP 116, with its matching core and variation databases), every result that differs is kept in the denominator, and the places where DuckVEP deliberately behaves differently are written down. Moving to a new VEP release means re-running the ledger, not re-arguing the rules.

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

That is **1,486,561 of 1,486,561** pairs identical to VEP, across human GRCh38 and GRCh37 and *Plasmodium falciparum* (with its own genetic codes), small variants, exact structural events, paired breakends, and regulatory and motif features. The release gate is deterministic: a single discordant, missing, extra or unresolved pair fails the run. HGVS is held to the same standard: 357,806 HGVSc and 357,806 HGVSp transcript pairs match VEP `--hgvs` string for string (56,998 of them from ClinVar chromosome 21), with 0 discordant.

Methods, per-term and per-impact tables, generated state-exploration campaigns and statistical caveats are in the [conformance report](benchmarks/duckvep_conformance.md). Ledger revisions predate the extraction from DuckHTS and are resolved through the [commit map](design/duckhts-commit-map.txt).

## Where it differs, and why

Compatibility is claimed only where it is measured, and differences are recorded rather than hidden. The [compatibility and errata record](ERRATA.md) keeps three verdicts apart (an observed VEP convention DuckVEP must follow, a DuckVEP gap, and a potential upstream erratum), each with its witness; the main cases:

- **Published release annotations are not the oracle.** Ensembl’s release VCFs disagree with executable VEP at some sites (in release 116, `X/Y:276322 G>A` is published as `intergenic_variant`, while the executable emits three `5_prime_UTR_variant` rows per chromosome). DuckVEP follows the executable; the [PAR witnesses](test/duckvep/conformance/README.md) pin both.
- **Breakends are evaluated one event at a time.** VEP 116’s buffered breakend path uses a chromosome-blind interval tree and can drop valid transcript pairs in multi-chromosome batches. The oracle runs with `--buffer_size 1` to isolate that; DuckVEP has no such batching effect ([paired-breakend differential](benchmarks/duckvep_conformance.md#paired-breakend-differential)).
- **Imprecise structural variants use the nominal span**, exactly as VEP’s registered predicates do, while `CIPOS`/`CIEND` stay on the row as metadata instead of being dropped.
- **No silent size limits.** VEP skips structural events above `--max_sv_size` (5 kb by default); DuckVEP annotates every exact span, and the oracle is run with a 10 Mb limit so the comparison covers them.
- **Circular regions execute on lifted intervals; MT is not "circular support".** A circular sequence region that carries an origin-crossing object (a transcript, exon, regulatory or motif feature with `start > end`) runs on a lifted linear copy of the model, so flanks, splice sites, CDS positions, HGVS 3' shifts and regulatory overlaps cross the origin consistently. Circular-coordinate support is independent of the mitochondrial codon table, which only selects a translation rule. Public crossing transcripts exist (19 in Ensembl Genomes 63, mostly bacterial) but VEP models them as reversed-bound intervals, so it is not an oracle for them: on three genomes with a VEP cache DuckVEP equals VEP on all other transcripts except flank rows that exist only through the origin, and the crossing object is proved by rotation equivariance and agreement with a linear model ([design](design/duckvep.md), [survey and differential](benchmarks/data/circular_source_survey.md)). A circular region without a wrapped object, such as human MT, keeps VEP's linear behavior byte for byte.
- **Unknown is a value.** When a result cannot be computed (missing reference sequence, reference mismatch, an unsupported symbolic allele) the row carries `duckvep_status` and a reason instead of a guess.
- **Outside the claim:** untested species and releases, phased multi-record haplotypes against an executable oracle, compound (cis) HGVS presentation, where the errata record the retained disagreements, and genomic (`g.`) and structural HGVS. See the [design contract](design/duckvep.md) and the [roadmap](https://github.com/RGenomicsETL/DuckVEP/issues).

## Speed

On the full GIAB HG002 GRCh38 benchmark (4,096,123 ALT alleles, VCF in, uncompressed table out, including DuckVEP’s explicit coordinate sort), DuckVEP finished **2.55×** faster than FastVEP on one core (65 s against 164 s) and **2.12×** on four, on a compact output contract. Complete-field comparisons, including CSQ with HGVS, and every caveat are in the [FastVEP benchmark](benchmarks/benchmark_duckvep_fastvep.md). The annotation kernel alone sustains 448,698 variants per second on one core with rich output ([throughput report](benchmarks/duckvep_throughput.md)).

## Try it

A model is loaded from ordinary relations: regions, transcripts with their CDS and flanking sequence, and exons. In production these come from an Ensembl core database through `query(duckvep_ensembl_regions_sql(...))` and `query(duckvep_ensembl_transcripts_sql(...))`; here a [one-transcript fixture](test/data/duckvep/readme.sql) stands in:

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

Alleles are rows too. Seven variants across the transcript, from the 5′ UTR to the intron:

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

The frameshift and the stop-gain escape nonsense-mediated decay because they fall in the first CDS positions (`nmd_escape_early_cds`). The empty HGVS fields are not failures to parse: this model was loaded without a reference FASTA, and those notations need flanking sequence to shift and render, so DuckVEP says so instead of guessing:

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

Phased genotypes become whole-haplotype consequences, again as a relation:

```sql
CREATE TABLE demo_calls AS
SELECT event_index, seq_region, position, reference, alternate,
       1::UINTEGER alt_index, 0::UINTEGER transcript_index, 0::UBIGINT sample_index,
       [1, 0]::INTEGER[] alleles, [false, true]::BOOLEAN[] phase_before, NULL::BIGINT phase_set
FROM demo_events WHERE position IN (124, 125);
SELECT transcript_index, carrier_count, sequence_status, length(contributors) AS contributors
FROM duckvep_haplotypes('SELECT * FROM demo_calls', 'demo');
#>
#> ┌──────────────────┬───────────────┬─────────────────┬──────────────┐
#> │ transcript_index │ carrier_count │ sequence_status │ contributors │
#> │      uint32      │    uint32     │     varchar     │    int64     │
#> ├──────────────────┼───────────────┼─────────────────┼──────────────┤
#> │                0 │             1 │ ok              │            2 │
#> └──────────────────┴───────────────┴─────────────────┴──────────────┘
```

## R interface

`Rduckvep` runs the same native builders on the caller’s DuckDB connection. Builder functions return SQL; `rduckvep_annotate()` executes annotation and returns a data frame. The supplied relation names remain visible in that connection, including TEMP tables and uncommitted rows.

```r
con <- Rduckvep::rduckvep_connect()
# After loading a model named "demo" and creating "demo_events":
sql <- Rduckvep::rduckvep_annotate_sql(con, "demo_events", "demo", hgvs = TRUE)
results <- Rduckvep::rduckvep_annotate(con, "demo_events", "demo", hgvs = TRUE)
DBI::dbDisconnect(con, shutdown = TRUE)
```

## With the rest of the stack

DuckVEP needs no other extension, but it is built to sit in a query next to them. With [DuckHTS](https://github.com/RGenomicsETL/duckhts) reading the VCF and [DuckClinVarbitration](https://github.com/RGenomicsETL/DuckClinVarbitration) supplying arbitrated ClinVar decisions, annotation and clinical evidence meet in one query (not evaluated here: it needs a full Ensembl model and both extensions):

```sql
-- One row per ALT allele, with the model's region index (DuckHTS reads the VCF)
CREATE TABLE alleles AS
SELECT row_number() OVER (ORDER BY r.seq_region, v.POS, alt_allele)::UBIGINT AS event_index,
       r.seq_region, v.POS::UBIGINT AS position, v.REF AS reference, alt_allele AS alternate,
       NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM read_bcf('cohort.vcf.gz') v CROSS JOIN unnest(v.ALT) AS t(alt_allele)
JOIN grch38_regions r ON r.seq_region_name = v.CHROM;

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

The extension uses the stable DuckDB C API (tested on DuckDB 1.5 and the 2.0 pre-release) and links its own htslib (for indexed reference FASTA) and cgranges. From R, [Rduckvep](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/) builds the same sources offline and provides connection, model and haplotype helpers.

## Documentation

- [Site](https://rgenomicsetl.github.io/DuckVEP/) and [R package reference](https://rgenomicsetl.github.io/DuckVEP/Rduckvep/)
- [Compatibility and errata](ERRATA.md): every known difference from VEP, classified, with witnesses
- [Design and implementation contract](design/duckvep.md): model build, the sweep, NMD, phased edits, structural events, HGVS, supplementary annotation
- Reports: [conformance](benchmarks/duckvep_conformance.md), [throughput](benchmarks/duckvep_throughput.md), [FastVEP comparison](benchmarks/benchmark_duckvep_fastvep.md), [haplotypes](benchmarks/duckvep_haplotypes.md)
- [Corpus workflow](design/duckvep_corpus_workflow.md) and [conformance harness](test/duckvep/conformance/README.md)

DuckVEP was developed inside [DuckHTS](https://github.com/RGenomicsETL/duckhts) and extracted with its history.

## License

GPL-2.0-or-later; see [LICENSE](LICENSE). Vendored htslib and cgranges keep their own licences. Ensembl data and VEP are © EMBL-EBI under their own terms.
