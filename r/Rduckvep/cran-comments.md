# Rduckvep 0.1.0 — CRAN submission comments

## Test environments

* local Ubuntu 24.04.3 LTS, R 4.6.0 (x86_64-pc-linux-gnu)
  `R CMD check --as-cran --no-manual` (CRAN incoming feasibility OK)
* GitHub Actions "Check Rduckvep": ubuntu-latest, macos-latest, windows-latest,
  windows-11-arm, wasm
* win-builder (devel and release) — TODO before submission

## R CMD check results

```
Status: OK
0 errors | 0 warnings | 0 notes
```

The local run used `--no-manual` and disabled remote incoming checks. A real
submission will additionally carry the standard `New submission` note.

## Notes for the reviewer

* **New submission.** Rduckvep is the R front end of the DuckVEP DuckDB
  extension; no reverse dependencies.
* **Installed size.** 10.4 Mb, of which 9.4 Mb is `inst/duckvep_extension`,
  the C sources of the bundled extension. `configure` compiles them from the
  vendored tree during installation; it performs no network access and
  installs nothing. `SystemRequirements: cmake, GNU make, zlib`.
* **Bundled third-party code.** HTSlib sources are vendored under
  `inst/duckvep_extension/third_party`; copyright and licences are recorded in
  `inst/COPYRIGHT`.
* **Suggests.** `Rduckhts` is on CRAN; `tinytest` runs the test suite.
* **Runs offline.** Examples and tests use a small bundled model and require no
  network, reference FASTA or Ensembl data.
* **Licence.** GPL (>= 2), matching the bundled DuckVEP extension.
