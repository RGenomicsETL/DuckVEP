# Paired single-core HG002 run

Measured 2026-10-05 with the biological implementation at `bbe2ec4` and the
integrated benchmark runner. `run_receipt.tsv` records executable and runner
hashes; `preparation.receipt.tsv` and `preparation_outputs.sha256` bind the inputs
and generated snapshot/mappings. Paths are measurement locators.

Three cold-process pairs, always vep-rs followed by DuckVEP, were pinned to
CPU 19 on the i5-13500: one E-core, maximum clock 3.5 GHz. Outputs were written
to disk, and other workloads were active on the host. These conditions differ
from the historical memory-backed measurements; the numbers are descriptive,
not a general speed qualification.

| Tool | Median wall (s) | Median CPU (s) | Median RSS (KiB) |
|---|---:|---:|---:|
| DuckVEP | 58.09 | 57.98 | 2,778,116 |
| vep-rs | 73.26 | 70.83 | 2,896,592 |

`stage_timings.tsv` contains one separate five-query probe after the paired runs:
compact count 5.788 s, rich count 12.991 s, rich-column touch 18.684 s,
rich joined Parquet 50.436 s, compact Parquet 28.054 s. These are cumulative
query variants, not additive components; extension setup and model restore are
outside their timers.

The common-transcript comparison contains 34,148,222 tuples per tool:
34,146,531 match, with 1,691 differences on each side. Executable VEP 116 agrees
with DuckVEP on 1,689, with vep-rs on one, and with neither on one.
`adjudication.txt` retains the complete counts, including missing transcript rows
in the expanded VEP replay.

The vep-rs cache is `homo_sapiens/116_GRCh38`, containing `transcripts/`.
The runner checks that layout and requires nonempty annotation output before
accepting timings. Its synthetic regression test also exercises one-sided
mismatches and absent VEP rows.

Reproduce with `benchmark_duckvep_vep_rs.sh --threads 1 --runs 3 --cpu-affinity 19`
and the input paths/environment described in its header. Large outputs remain
outside Git in `/root/duckvep/work/luna-six/artifacts/integration/veprs-{work,out}`.
