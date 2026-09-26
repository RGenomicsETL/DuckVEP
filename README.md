
# DuckVEP

DuckVEP implements Ensembl VEP 116 consequence semantics as a resident, relation-first DuckDB kernel. Models are loaded from ordered relations; annotation, transcript HGVS, and phased haplotypes remain queryable relations. The extension uses htslib for indexed reference access, not a DuckHTS runtime dependency.

## Evidence

| Evidence                                 | Pairs / variants |   Exact | Different | Unresolved | Variants/s | Receipt                                                                      |
|:-----------------------------------------|-----------------:|--------:|----------:|-----------:|-----------:|:-----------------------------------------------------------------------------|
| final-dbsnp                              |           73,620 |  73,620 |         0 |          0 |          — | [conformance history](test/duckvep/conformance/data/conformance_history.csv) |
| final-grch37                             |          486,464 | 486,464 |         0 |          0 |          — | [conformance history](test/duckvep/conformance/data/conformance_history.csv) |
| plasmodium-falciparum-63                 |           40,732 |  40,732 |         0 |          0 |          — | [conformance history](test/duckvep/conformance/data/conformance_history.csv) |
| GRCh38 HGVSc suffix agreement (present)  |           44,871 |  44,871 |         0 |          0 |          — | [executable-VEP comparison](benchmarks/data/duckvep_fastvep/conformance.csv) |
| GRCh38 HGVSp suffix agreement (present)  |           20,782 |  20,782 |         0 |          0 |          — | [executable-VEP comparison](benchmarks/data/duckvep_fastvep/conformance.csv) |
| GRCh38 sorted rich-output, single thread |          100,957 |       — |         — |          — |    448,698 | [timing receipt](benchmarks/data/duckvep_throughput.csv)                     |

Each corpus result applies to the named receipt and model, not an arbitrary deployment sample. See the [conformance report](benchmarks/duckvep_conformance.md) and [throughput report](benchmarks/duckvep_throughput.md) for methods and limitations.

## Query a model

The [small model setup](test/data/duckvep/readme.sql) loads a transcript and exons and stages a variant from a [committed Parquet fixture](test/data/duckvep/minimal_bcsq.parquet). The README runs against a locally built extension in a DuckDB CLI session.

``` sql
SELECT event_index, consequence_mask, transcript_hgvs, transcript_hgvs_status
FROM duckvep_annotate('readme_events', 'readme', hgvs := true);
#> ┌─────────────┬──────────────────┬─────────────────┬────────────────────────┐
#> │ event_index │ consequence_mask │ transcript_hgvs │ transcript_hgvs_status │
#> │   uint64    │      uint64      │     varchar     │        varchar         │
#> ├─────────────┼──────────────────┼─────────────────┼────────────────────────┤
#> │           1 │             8192 │ c.5T>C          │ supported              │
#> └─────────────┴──────────────────┴─────────────────┴────────────────────────┘
```

Calls are relations too:

``` sql
CREATE TABLE readme_calls AS SELECT event_index, seq_region, position,
  reference, alternate, 1::UINTEGER alt_index, 0::UINTEGER transcript_index,
  0::UBIGINT sample_index, [1,0]::INTEGER[] alleles,
  [false,true]::BOOLEAN[] phase_before, NULL::BIGINT phase_set
FROM readme_events;
SELECT transcript_index, carrier_count, sequence_status,
  length(contributors) AS contributors
FROM duckvep_haplotypes('SELECT * FROM readme_calls', 'readme');
#>
#> WARNING:
#> Deprecated lambda arrow (->) detected. Please transition to the new lambda syntax, i.e.., lambda x, i: x + i, before DuckDB's next release.
#> Use SET lambda_syntax='ENABLE_SINGLE_ARROW' to revert to the deprecated behavior.
#> For more information, see https://duckdb.org/docs/stable/sql/functions/lambda.html.
#>
#> ┌──────────────────┬───────────────┬─────────────────┬──────────────┐
#> │ transcript_index │ carrier_count │ sequence_status │ contributors │
#> │      uint32      │    uint32     │     varchar     │    int64     │
#> ├──────────────────┼───────────────┼─────────────────┼──────────────┤
#> │                0 │             1 │ ok              │            1 │
#> └──────────────────┴───────────────┴─────────────────┴──────────────┘
```

The [design](design/duckvep.md) describes model and consequence contracts. [Rduckvep](r/Rduckvep) builds the same sources offline and provides a connection and haplotype front end; DuckHTS is optional for VCF reading.
