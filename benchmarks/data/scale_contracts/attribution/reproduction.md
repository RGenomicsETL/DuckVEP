# Frozen scale-contract reproduction

Source revision: `2f4d8526d5e29103731073859c06b792dbd1f89d` (already on origin/main; the extension sources and bundled R extension sources in this checkout match that commit). Extension binary: immutable `/tmp/scale3-restored/duckvep.duckdb_extension`, SHA-256 `5e77dc2c175901afc635121b0958729e8b6a59ae880e82a4ffe5e77fb191b61e`. Runtime: DuckDB R 1.5.5, Linux. Command: `DUCKVEP_SCALE_EXTENSION=/tmp/scale3-restored/duckvep.duckdb_extension Rscript benchmarks/scale_contracts.R`.

The frozen `inputs.tsv`, `eligibility.tsv`, `schemas.tsv`, and `receipts.tsv` checks passed: compact, compact_regulation, complete17, and 41,942 retained cells, model SHA-256 `0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b`. The public annotation paths use native SQL builders through `FROM query(duckvep_annotate_sql(...))` and the projected SQL builder. This verifies the restored origin/main-equivalent extension, **not** the rejected keyset probe. No new publishable receipt or remote operation was produced.

## Scatter loader reproduction

The unordered scatter loader (`096b945`, binary SHA-256 `b831bd915f138d2b3fb0d267663b8ddb02c3ab72c4916b56d6ed1394def0fc11`, run from an immutable copy) passed the same `benchmarks/scale_contracts.R` check unchanged: compact, compact_regulation, complete17, 41,942 retained cells, model SHA-256 `0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b`. No new publishable receipt was produced.

## Scatter loader rebased onto circular-topology main

After rebasing `scale-slice3` onto origin/main (which carries the #19 circular contract), the merged loader (immutable copy, SHA-256 `4907665be5961850d2e70c1f18c035a3a20bab38b70c9b61762d5ce87b42f813`) again passed `benchmarks/scale_contracts.R` unchanged: compact, compact_regulation, complete17, 41,942 retained cells, model SHA-256 `0301905915b038b6b1f5db04d7dda1041ac5e30a7732ec700dd82b1616402c1b`. No new publishable receipt was produced.
