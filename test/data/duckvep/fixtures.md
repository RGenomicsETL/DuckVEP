# Standalone DuckVEP SQL fixtures

The Parquet files are snapshots produced with DuckDB's `COPY ... TO ... (FORMAT PARQUET)` after loading a local DuckHTS extension build. The source-tree revision of that binary was not independently verified; the committed fixture bytes are the test inputs. The fixture sources are committed here; the extension under test reads Parquet with DuckDB's built-in reader.

| Parquet file | Source | DuckHTS reader |
|:--|:--|:--|
| `ensembl_release_consequences.parquet` | `ensembl_release_consequences.vcf` | `read_bcf` |
| `minimal_bcsq.parquet` | `minimal_bcsq.vcf` | `read_bcf` |
| `minimal_gff.parquet` | `minimal.gff3` | `read_gff(..., attributes_map := true)` |
| `geno_phase_partial.parquet` | `../geno_phase_partial.vcf` | `read_geno` |
| `geno_vcf44.parquet` | `../geno_vcf44.vcf` | `read_geno` |
| `ce_1000_5999.parquet` | `../ce.fa`, `CHROMOSOME_I:1000-5999` | `fasta_nuc(..., bin_width := 5000, include_seq := TRUE)`; `string_agg(seq, '' ORDER BY start)` |

The VCF, GFF and FASTA inputs are not extension dependencies. `readme.sql` defines an in-memory model from a Parquet variant. The optional DuckHTS integration script reads its VCF at runtime. Set `DUCKHTS_EXTENSION` to a core-only DuckHTS build that does not register the same `duckvep_*` functions: two independent model registries behind identical SQL function names cannot interoperate.
