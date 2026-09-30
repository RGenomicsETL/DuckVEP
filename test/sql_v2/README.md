# v2 host tests

- `equality_cases.sql`: queries run on both hosts; results are compared with
  `equality_golden.json`, the recorded output of the v1 host (`run_v2_tests.py --record`).
  A case named `... error` must fail with a `duckvep_*:` message; every other case must succeed.
- `v2_native.sql`: v2-only assertions (registration kinds, no catalog objects, NULL
  invariants of struct results, table-function rescans).
- Run with `DUCKVEP_V2_DUCKDB=<DuckDB CLI built at the pinned revision> make test_v2`.
  The PyPI `duckdb>=2.0.0.dev0` wheels are different snapshots of a preview ABI that is not
  frozen and cannot load this extension; see `docs/v2-host.md`.
