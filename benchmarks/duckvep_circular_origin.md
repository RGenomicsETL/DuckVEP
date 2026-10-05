# Origin-focused circular workload

`benchmarks/duckvep_circular_origin.py` measures the lifted-interval execution of circular regions
on a synthetic 100 kb circle with 400 transcripts (70% coding, one to five exons) and 120 regulation features packed
within 15 kb of the origin, and 355,446 events, 70% of them within 4 kb of the origin (60% SNVs, deletions, insertions and
replacements up to 12 bases). Every run is a fresh child process that loads an immutable, hashed copy of the extension
(`extension_sha256` in the CSV), so `peak_rss_mb` (kernel `VmHWM`) belongs to that run alone; `loaded_rss_mb` is the
same process after loading the model and events and before annotating. Each configuration is the median of three
passes of `duckvep_annotate_sql` with 1 kb flanks, aggregated in DuckDB so the 12.2 M result rows are not
materialized.

Three frames annotate the same events:

- **lifted**: the origin lies inside the object cloud, so most objects wrap and the region runs on the lifted model;
- **lifted_rot**: the same world rotated by 777 bases, a different wrap pattern;
- **linear**: rotated by half the circle so that nothing wraps; the model runs on the ordinary linear kernel. This is
  the control for the cost of lifting.

Events within 300 bases of any frame's seam are excluded from all three, so the linear control can hold them.

**Output equality.** The checksum is a sum of a hash over (event, transcript or feature stable identifier, consequence
mask, region mask, impact, cDNA/CDS/protein position, HGVSc, HGVSp, HGVS shift, NMD), not over frame-specific
ordinals. It is identical across 1, 4 and 8 threads and across all three frames: `112207876093400761556032700` without HGVS and `112183060793327480582475232` with
HGVS (12,163,925 rows each). Rotation equivariance therefore holds at this scale on the real annotation path, and
the lifted result equals the linear-kernel result on the same world.

| Frame | HGVS | Threads | Events/s | vs linear, 1 thread | Peak RSS (MB) |
| --- | --- | ---: | ---: | ---: | ---: |
| lifted | no | 1 | 139,006 | 0.83 | 326.5 |
| lifted | no | 4 | 505,652 |  | 391.1 |
| lifted | no | 8 | 875,773 |  | 457.4 |
| lifted | yes | 1 | 56,500 | 0.80 | 349.7 |
| lifted | yes | 4 | 209,939 |  | 488.7 |
| lifted | yes | 8 | 352,332 |  | 643.2 |
| lifted_rot | no | 1 | 138,527 | 0.83 | 331.5 |
| lifted_rot | no | 4 | 519,079 |  | 398.0 |
| lifted_rot | no | 8 | 863,064 |  | 478.6 |
| lifted_rot | yes | 1 | 56,886 | 0.81 | 350.3 |
| lifted_rot | yes | 4 | 209,177 |  | 478.6 |
| lifted_rot | yes | 8 | 358,026 |  | 651.2 |
| linear | no | 1 | 166,818 | 1.00 | 331.9 |
| linear | no | 4 | 603,951 |  | 399.2 |
| linear | no | 8 | 1,027,782 |  | 476.3 |
| linear | yes | 1 | 70,402 | 1.00 | 357.0 |
| linear | yes | 4 | 266,996 |  | 481.7 |
| linear | yes | 8 | 455,212 |  | 675.5 |

Host: 13th Gen Intel(R) Core(TM) i5-13500, DuckDB Python 1.5.2, revision `a88399998df62955f3e72ed57b71a9b99d3ff7c9`. Lifting costs about a fifth of one-thread throughput (0.83 without HGVS, 0.80
with): each object is admitted at three images and every lifted row is then resolved to one row per event and object.
The profile was not decomposed further. Peak memory is within a few percent of the linear control,
because result aggregation, not the model, dominates the resident set at this size. Scaling from 1 to 8 threads is
6.3x without HGVS and 6.2x with HGVS on the lifted frame.

Reproduce with `python3 benchmarks/duckvep_circular_origin.py --extension build/release/extension/duckvep/duckvep.duckdb_extension
--out benchmarks/data/duckvep_circular_origin.csv`.
