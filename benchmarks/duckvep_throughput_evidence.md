# Throughput: methods and measurement scope

[The visual report](duckvep_throughput.md) separates resident annotation,
dense-region threading, whole-genome composition and interval joins. Each chart
uses a complete measured configuration, with source revisions retained in the CSV.
The chart source pins exact receipt revisions; the render fails if paired inputs
or output checks disagree.

## Resident annotation

The timer covers typed result production, list expansion and checksum aggregation.
Model loading, input staging and warm-up are outside the measured passes. Compact,
rich, HGVS and fused rich+HGVS are different output contracts; a count-only query is
not interchangeable with a materialized file write.

The public-relation chart uses the full literal HG002 relation, 5,000-base transcript
halo, 1,383,580 regulatory/motif features, five measured passes and 100,000 warm-up
alleles. Revision `660a4ed9` (2026-09-30) supplies all four output modes at both
one and four threads. Source corpus, staged input, physical model, allele count and
output-row count agree across those rows. Each output mode retains its own checksum
and full-row fingerprint contract.

The full-core chart uses revision `03803bd8` (2026-07-21): a complete same-revision
pair of ClinVar and HG002 corpora, each with compact, rich and cumulative-HGVS
measurements on one core.
Its Ensembl 116 GRCh38 model contains 644,427 transcripts, 380,818 regulatory
features and 1,002,762 motif features. Source, extension, physical/logical model,
reference, region-ordinal and staged-corpus digests bind each receipt. Public-row
XOR/sum fingerprints have schema version 2.

Whiskers are the observed minimum and maximum pass times, not confidence intervals.
Input alleles per second and output rows per second answer different questions:
a corpus with greater transcript fan-out produces more work per input allele.
Comparisons across source revisions require identical input and output contracts.

## Dense-region parallelism

The panel contains 517,097 ClinVar alleles in 318 annotation-dense tiles. It tests
transcript halos of 0, 5,000, 10,000 and 50,000 bases. One ordered partition runs on
CPU 2; four run on CPUs 2,4,6,8. Five passes follow 100,000 warm-up alleles.

At each distance, both configurations retain the same allele count, output-row
count, model and staged-input digests, and XOR/sum fingerprints. These checks
establish row-multiset equality; parallel emission order is not the contract.
The measurements cover four partitions on this panel, not arbitrary thread counts.

Peak RSS is GNU time's process high-water mark. It includes loading, DuckDB and
allocator/page effects; it is not an attribution of C allocations to workers.
The separate [model-resource receipt](data/duckvep_model_resources.csv) describes a
transcript-only model and must not be pooled with regulation-enabled process RSS.

## Whole-genome composition

The DeepVariant HG002 40x PCR-free GRCh38 WGS campaign scans and canonicalizes
literal alleles on chromosomes 1–22, X, Y and MT. It requests rich consequences,
HGVSc/HGVSp and core regulation, then joins dated ClinVar, ClinvArbitration,
AlphaMissense, gnomAD v2.1.1 gene constraint and Ensembl regulatory intervals.
The complete relation is written as ZSTD Parquet.

These are warm-page-cache, unpinned integration measurements. Four writers,
a 4 GB DuckDB buffer-manager limit and per-worker output files bound the materialized
query. Native immutable-model allocations are outside that buffer-manager ceiling.
The statement timer excludes setup when measuring a query; GNU time records process RSS.

Annotation-only and all-provider queries overlap in work and are not additive stages.
The [output fingerprints](data/duckvep_human_annotation_fingerprint.csv) agree in row
count and order-independent XOR/sum over every projected column for the two retained
composition configurations. The CSV records both configurations and their resources.

## Interval joins

The equivalence campaign checks every matched allele's stable IDs and SO terms.
All three methods return 745,252 overlap pairs for 414,813 alleles, with zero mismatched
allele rows.

- A chromosome equality join produces same-chromosome candidates before residual
  range filtering; its physical plan is `HASH_JOIN`.
- Two packed RegionKey inequalities encode chromosome and half-open overlap and
  select `IE_JOIN` without a string equality.
- A cgranges bulk index supports arbitrary literal contig names and has its own
  measured time/memory trade-off.

`EXPLAIN` identifies the plan. RSS is available for IEJoin and cgranges, not for the
chromosome-hash receipt. An absent measurement is not zero.

## Receipts

- [Resident annotation history](data/duckvep_throughput.csv): all recorded surfaces,
  source revisions, pass ranges, hardware, input/model identities and output checks.
- [Dense threading](data/duckvep_dense_threading.csv): halos, partitions, exact row
  counts, process RSS and cross-thread fingerprints.
- [Model resources](data/duckvep_model_resources.csv): load/drop RSS and planned scratch.
- [Whole-genome statements](data/duckvep_human_annotation.csv) and
  [full-output fingerprints](data/duckvep_human_annotation_fingerprint.csv).
- [Interval equivalence](data/duckvep_human_interval_validation.csv).

## Reproduction

Build the pinned extension and provide the model and staged corpus identified by
the receipt. A resident compact measurement has this form:

```sh
taskset -c 2 Rscript benchmarks/duckvep_throughput.R \
  --database /path/to/model.duckdb \
  --variants-database /path/to/staged-corpus.duckdb \
  --variants-table bench_variants --corpus-source /path/to/source.vcf.gz \
  --workload-name ensembl116_grch38_full_corpus_regulation \
  --regulatory --variants 4095611 --passes 5 --warmup 100000 \
  --threads 1 --input-partitions 1 --transcript-distance 5000 \
  --output compact --fingerprint /path/to/full-public-row-fingerprint.csv
```

For HGVS, supply `--reference-fasta` and select `--output hgvs`. The runner records
source, binary, model, reference, corpus and public-row identities. Rendering the
report reads the committed receipts and performs no annotation benchmark:

```sh
Rscript benchmarks/scripts/render_benchmarks.R throughput
```
