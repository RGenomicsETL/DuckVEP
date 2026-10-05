# Rduckvep

R front end for [DuckVEP](https://rgenomicsetl.github.io/DuckVEP), a variant consequence engine for DuckDB with Ensembl VEP semantics.

The package bundles DuckVEP extension sources and builds them offline for the installed `duckdb` R package. It provides:

- `rduckvep_connect()` and `rduckvep_load()` for a DuckDB connection with the extension loaded;
- `rduckvep_annotate()` for annotation results as a data frame, and `rduckvep_annotate_sql()` for the SQL builder;
- `rduckvep_haplotypes()`, `rduckvep_coding_calls()` and `rduckvep_coding_transcripts()` for phased coding calls;
- wrappers for Ensembl model builders, transcript projection, model receipts, loss-of-function calls, breakend evidence and structural HGVS.

Model loading uses `duckvep_model_load()` through DBI on the caller's connection. Its relations can include temporary tables and uncommitted rows. See the [project page](https://rgenomicsetl.github.io/DuckVEP/) for compatibility evidence, differences from VEP and worked examples.

```r
con <- Rduckvep::rduckvep_connect()
# After loading a model and creating an events relation:
Rduckvep::rduckvep_annotate(con, "events", "model", hgvs = TRUE)
DBI::dbGetQuery(con, "SELECT function_name FROM duckdb_functions() WHERE function_name LIKE 'duckvep%'")
```

Install from a source clone with `Rscript r/Rduckvep/bootstrap.R .` followed by `R CMD INSTALL r/Rduckvep`. [DuckHTS](https://github.com/RGenomicsETL/duckhts) (`Rduckhts`) is optional for VCF/BCF input.
