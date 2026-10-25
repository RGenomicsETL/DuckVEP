# Haplotype scale measurements

The [haplotype contract](../../../design/duckvep_haplotype_contract.md) asks for complete output within fixed memory caps, and for whole-genome prediction at least twice as fast as pinned `bcftools csq` on the same input, model and core. Two runners measure this on GIAB HG002 v4.2.1 (4,023,088 records, 4,070,522 ALT alleles) with the Ensembl 116 GRCh38 model. Every timed process is fresh and pinned to one core.

## DuckVEP against csq

[`run_hg002_csq.sh`](../../haplotype_scale/run_hg002_csq.sh) alternates csq (`-p a -Ou -o /dev/null`), the DuckDB CLI and R for three rounds on core 6. DuckVEP reads the VCF with its native reader (`duckvep_coding_calls`), stages the calls and writes the complete `duckvep_haplotypes` output to Parquet; model load is inside the timer. Each `hg002_<revision>/` directory holds `process.tsv` (every process), `identity.txt` (input, model, tool and extension hashes, and the output checksum) and `status_counts.tsv`.

| `hg002_2417ac9` | wall s (runs) | median | csq median / median |
|---|---|---|---|
| `bcftools csq` | 17.05 / 17.51 / 16.87 | 17.05 | |
| DuckVEP, DuckDB CLI | 8.80 / 8.74 / 8.71 | 8.74 | 1.95 |
| DuckVEP, from R | 9.41 / 9.10 / 9.12 | 9.12 | 1.87 |

This run is just under the two-times requirement. Inside DuckVEP's 8.7 s, model load takes 2.3 s, reading and staging 3.1 s, and prediction with Parquet output 3.2 s.

## Qualification

[`run_qualification.sh`](../../haplotype_scale/run_qualification.sh) runs every process inside a systemd scope capped at 16 GiB with no swap, a 4 GiB native budget, a DuckDB `memory_limit` of 8 GB, a 32 GiB temp directory and one DuckDB thread. It starts a run only when the one-minute load is at most 4, and alternates csq with DuckVEP. Its R worker has two pipelines:

- **Mode B** reads the VCF text with DuckDB's `read_csv`, discovers coding transcripts with `duckvep_coding_transcripts` and parses genotypes in SQL.
- **Mode A** starts from the calls relation already staged in Parquet, so it times prediction alone.

A cold run is a fresh process; a warm pass is a second execution in the same process. Both modes write the same output, checked by a full-output checksum. [`summarize.R`](../../haplotype_scale/summarize.R) turns the raw receipts into `runs.csv` and `summary.csv`, and [`tables.R`](../../haplotype_scale/tables.R) prints the tables below.

### HG002

| invocation | wall s (runs) | median wall s | median user s | peak RSS KiB | memory.peak MiB | load at start |
|---|---|---|---|---|---|---|
| `-p a -Ob -o file.bcf` | 28.53 / 28.42 / 29.10 | 28.53 | 27.64 | 811,312 | 860 | 1.1 / 1.2 / 1.2 |
| `-p a -Ou -o /dev/null` | 17.01 / 17.09 / 16.93 | 17.01 | 16.35 | 811,312 | 792 | 0.9 / 1.3 / 1.4 |

| mode | cold process wall s (runs) | cold wall median | load | stage | predict + Parquet | warm pass median | warm stage | warm predict |
|---|---|---|---|---|---|---|---|---|
| B | 11.51 / 11.74 / 11.47 | 11.51 | 2.29 | 5.67 | 3.10 | 8.64 | 5.50 | 3.12 |
| A | 6.09 / 6.16 / 6.06 | 6.09 | 2.31 | 0.03 | 3.25 | 3.19 | 0.03 | 3.16 |

Times are seconds; load, stage and predict are the worker's own timers. The output has 158,137 rows: 157,597 whole-product predictions, 212 per-record predictions and 328 `incomplete_input` paths. Of the 4,023,088 records, 25,969 touch coding sequence and become 247,374 calls; the rest reach no call.

Split by the committed csq accounting (mode A, cold):

| partition | records with calls | calls | output rows | predicted | predict + Parquet s (median) | load s |
|---|---|---|---|---|---|---|
| inside | 23865 | 239300 | 153654 | 153654 | 3.18 | 2.30 |
| outside | 2104 | 8074 | 6462 | 6292 | 0.13 | 2.30 |

### 5M records on the MANE model

The input is HG002 plus 976,912 gnomAD v4.1 exome sites chosen by a fixed hash: 5,000,000 records and 5,047,434 ALT alleles, one sample, all phased. The model is the MANE Select subset (19,297 transcripts). [`prepare_inputs.R`](../../haplotype_scale/prepare_inputs.R) builds the inputs and records the rule.

| mode | cold wall s (runs) | median | load | stage | predict | warm pass | output checksum (audited runs) |
|---|---|---|---|---|---|---|---|
| B | 15.96 / 16.06 / 15.80 | 15.96 | 0.93 | 5.80 | 8.74 | 14.39 | 336410809266679938646157 |
| A | 10.28 / 10.34 / 10.18 | 10.28 | 0.92 | 0.04 | 8.85 | 8.81 | 336410809266679938646157 |

294,000 records reach a call. The 36,443 output rows are 35,682 whole-product predictions, 331 per-record predictions, 279 `edit_conflict` paths and 151 `incomplete_input` paths. At most 916 events share one transcript.

### Controls

`dense` places 4,000 variants on one haplotype of the longest MANE coding sequence (TTN); `low_sharing` places one variant on each MANE transcript.

| control | median wall s | predict s | calls | output rows | carriers | translated bases | max events per transcript | transcripts | output bytes |
|---|---|---|---|---|---|---|---|---|---|
| dense_B | 3.18 | 1.79 | 4000 | 2 | 2 | 215914 | 4000 | 1 | 553338 |
| low_sharing_B | 1.98 | 0.49 | 19116 | 19099 | 19099 | 33382784 | 2 | 19087 | 27884526 |
| dense_A | 3.16 | 1.80 | 4000 | 2 | 2 | 215914 | 4000 | 1 | 553338 |
| low_sharing_A | 1.93 | 0.49 | 19116 | 19099 | 19099 | 33382784 | 2 | 19087 | 27884526 |

### Caps observed

Maximum over the three cold processes of each job. No job spilled to disk.

| job | memory.peak GiB (16 cap) | max RSS GiB | native high-water MiB, load (total) | native high-water MiB per owner, pass | native pass total MiB (4096 cap) | DuckDB spill peak bytes |
|---|---|---|---|---|---|---|
| hg002_B | 4.14 | 4.18 | 1,595 | model 1,317, index 14, workspace 346, scratch 0, emit 0, control 0 | 1,677 | 0 |
| hg002_A | 3.90 | 3.93 | 1,595 | model 1,317, index 14, workspace 346, scratch 0, emit 0, control 0 | 1,677 | 0 |
| qual5m_B | 3.14 | 3.11 | 98 | model 75, index 0, workspace 344, scratch 0, emit 0, control 0 | 418 | 0 |
| qual5m_A | 2.56 | 2.56 | 98 | model 75, index 0, workspace 344, scratch 0, emit 0, control 0 | 418 | 0 |
| dense_B | 2.34 | 2.38 | 98 | model 75, index 0, workspace 1,110, scratch 0, emit 0, control 0 | 1,185 | 0 |
| low_sharing_B | 1.99 | 2.02 | 98 | model 75, index 0, workspace 344, scratch 0, emit 0, control 0 | 418 | 0 |

### Failures are explicit

[`failures.R`](../../haplotype_scale/failures.R) checks that exhausting a budget or capacity raises an error and leaves the connection usable (`failures.tsv`):

- A 64 MiB native budget refuses the model load, publishes no model, and the load succeeds after the budget is restored.
- A budget of the resident model plus 64 MiB refuses `duckvep_haplotypes`; no table or output file is published, and native bytes return to the resident value.
- `max_alignment_cells := 100000` refuses the HG002 calls with the cell count it needed; nothing is published.
- The builder's default capacities complete the HG002 job (158,137 rows).

## Files

`runs.csv` has one row per process and pass, `summary.csv` the medians and ratios, `inputs.tsv` and `input_counts.tsv` the input hashes and counts, and `failures.tsv` the failure controls. Regenerate with:

```sh
bash benchmarks/prepare_hg002.sh /root/duckvep/data/hg002-csq
Rscript benchmarks/haplotype_scale/prepare_inputs.R WORK
WORK=WORK OUT=RESULTS bash benchmarks/haplotype_scale/run_qualification.sh
Rscript benchmarks/haplotype_scale/summarize.R RESULTS
Rscript benchmarks/haplotype_scale/tables.R
OUT=benchmarks/data/haplotype_scale/hg002_REVISION bash benchmarks/haplotype_scale/run_hg002_csq.sh
```

Large inputs and outputs stay outside Git.
