# Rduckvep

R front end for [DuckVEP](https://rgenomicsetl.github.io/DuckVEP/), a variant consequence engine for DuckDB with Ensembl VEP semantics.

The package bundles the DuckVEP extension sources and builds them for the installed `duckdb` R package, offline. It provides:

- `rduckvep_connect()` and `rduckvep_load()`: a DuckDB connection with the extension loaded;
- `rduckvep_haplotypes()`: whole-haplotype consequences from phased calls;
- `rduckvep_annotate()`: annotated events as a data frame;
- `rduckvep_annotate_sql()`, `rduckvep_transcript_projection_sql()`, `rduckvep_ensembl_regions_sql()`, `rduckvep_ensembl_transcripts_sql()`, `rduckvep_ensembl_regulation_features_sql()` and `rduckvep_model_receipt_sql()`: SQL from the extension's native builders.

Model loading remains SQL through DBI: `duckvep_model_load()` compiles model relations on that connection. See the [project page](https://rgenomicsetl.github.io/DuckVEP/) for the evidence, the documented differences from VEP and worked examples.

```r
con <- Rduckvep::rduckvep_connect()
# After loading a model and creating an events relation:
# Rduckvep::rduckvep_annotate(con, "events", "model", hgvs = TRUE)
DBI::dbGetQuery(con, "SELECT function_name FROM duckdb_functions() WHERE function_name LIKE 'duckvep%'")
```

Install from source in a clone: `Rscript r/Rduckvep/bootstrap.R .` then `R CMD INSTALL r/Rduckvep`. [DuckHTS](https://github.com/RGenomicsETL/duckhts) (Rduckhts) is optional, for reading VCF/BCF.
