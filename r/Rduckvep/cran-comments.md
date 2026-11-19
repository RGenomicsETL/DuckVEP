# Rduckvep 2.0.0.0.1.0 — CRAN submission comments

## Test environments

* local Ubuntu 24.04.3 LTS, R 4.6.2 (x86_64-pc-linux-gnu), R `duckdb` 2.0.0:
  `R CMD check --as-cran`, manual and remote incoming checks included
* GitHub Actions "Check Rduckvep": ubuntu-latest, macos-latest, windows-latest,
  windows-11-arm, wasm
* win-builder (devel and release) — TODO before submission

## R CMD check results

```
Status: 1 NOTE
0 errors | 0 warnings | 1 note
```

```
* checking CRAN incoming feasibility ... NOTE
New submission
```

The version `2.0.0.0.1.0` is the bundled extension version followed by the
package revision.

## Notes for the reviewer

* **New submission.** Rduckvep is the R front end of the DuckVEP DuckDB
  extension; no reverse dependencies.
* **Installed size.** 10.6 Mb, of which 9.6 Mb is `inst/duckvep_extension`,
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
