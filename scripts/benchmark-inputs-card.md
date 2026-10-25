---
license: other
license_name: giab-nist-and-gnomad-cc0
pretty_name: DuckVEP benchmark inputs
tags:
- genomics
- variant-effect-prediction
- benchmark
- duckdb
---

# DuckVEP benchmark inputs

The variant inputs of the [DuckVEP](https://github.com/RGenomicsETL/DuckVEP) benchmarks, as the VCFs the runners read and as Parquet with one file per chromosome. Each Parquet relation has one row per ALT allele and was checked, row for row, against `bcftools query` of the VCF it comes from. `SHA256SUMS` lists every file.

| Directory | Records | ALT alleles | Used by |
|---|---:|---:|---|
| `hg002/` | 4,023,088 | 4,070,522 | vep-rs comparison, haplotypes against `bcftools csq`, scale qualification |
| `qual5m/` | 5,000,000 | 5,047,434 | 5M-record scale qualification |
| `controls/dense`, `controls/low_sharing` | 4,000 and 19,065 | 4,000 and 19,065 | scale qualification stress controls |
| `throughput/` | 4,048,342 | 4,095,611 | annotation throughput |

## Columns

The Parquet files under `hg002/alleles/`, `qual5m/alleles/` and `controls/` have:

| Column | Type | Meaning |
|---|---|---|
| `record_index` | UBIGINT | 0-based ordinal of the record in the VCF beside it |
| `chrom` | VARCHAR | Ensembl contig name (`1`..`22`, `X`, `Y`); also the partition directory |
| `position` | UBIGINT | VCF `POS` |
| `reference`, `alternate` | VARCHAR | `REF` and one `ALT` allele |
| `alt_index` | UTINYINT | 1-based index of that allele in `ALT` |
| `filter` | VARCHAR | VCF `FILTER` |
| `gt`, `ps` | VARCHAR | The sample's `GT` and `PS` as written; `PS` may be a number or a label such as `PATMAT` or `HOMVAR` |
| `origin` | VARCHAR | `qual5m` only: `giab_hg002` or `gnomad_v4.1_exomes` |

`throughput/alleles/` has `record_index`, `chrom`, `position`, `reference`, `alternate` and `alt_index`, where `record_index` refers to GIAB's `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`.

Read a relation in DuckDB with `httpfs`:

```sql
SELECT chrom, count(*) AS alleles
FROM read_parquet('hf://datasets/REPO_ID/hg002/alleles/*/*.parquet', hive_partitioning = false)
GROUP BY chrom ORDER BY chrom;
```

To annotate with DuckVEP, join `chrom` to a model's `model_regions.seq_region_name` to get `seq_region`, and number the rows as `event_index`.

## Provenance

These files are derived from data of the National Institute of Standards and Technology (NIST) Genome in a Bottle Consortium and of gnomAD. They were changed as follows.

- **`hg002/hg002.ens.vcf.gz`** is GIAB's [`HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/) (SHA-256 `d027d7f0…2bf6`). 25,254 chr6 MHC records whose sample column has more fields than `FORMAT` declares were removed, and `chr1`..`chr22` were renamed `1`..`22`. [`benchmarks/prepare_hg002.sh`](https://github.com/RGenomicsETL/DuckVEP/blob/main/benchmarks/prepare_hg002.sh) repeats these steps.
- **`qual5m/qual5m.vcf.gz`** is that file plus 976,912 gnomAD v4.1 exome sites (`FILTER` PASS, absent from HG002, chosen by a fixed hash) with synthetic phased genotypes. **`controls/`** holds gnomAD v4.1 exome sites with synthetic genotypes: 4,000 on one haplotype of the TTN coding sequence (`dense`), and 19,065 sites placed one per MANE transcript (`low_sharing`). [`benchmarks/haplotype_scale/prepare_inputs.R`](https://github.com/RGenomicsETL/DuckVEP/blob/main/benchmarks/haplotype_scale/prepare_inputs.R) records the rules.
- **`throughput/`** holds the literal ALT alleles of GIAB's [`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/) (SHA-256 `adb4d4a5…175c`), with `chr` removed from contig names.

## Terms and citation

GIAB data are NIST data: not subject to copyright in the United States, and free to modify and redistribute when the changes are stated, as above, and NIST is acknowledged as the source ([NIST data license](https://www.nist.gov/open/license)). For the v4.2.1 benchmark, cite [Wagner et al.](https://doi.org/10.1101/2020.07.24.212712) and, for the sequencing data, [Zook et al. 2016](https://doi.org/10.1038/sdata.2016.25).

This dataset includes data from the [gnomAD](https://gnomad.broadinstitute.org) v4.1 release. gnomAD's primary data are available under the [Creative Commons Zero Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/); please acknowledge the gnomAD project. Only site coordinates and alleles are used here, with synthetic genotypes.

HG002 is a GIAB reference sample consented for public and commercial redistribution.
