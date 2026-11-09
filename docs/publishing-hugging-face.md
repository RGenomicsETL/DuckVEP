# Publish datasets to Hugging Face

This procedure stages data locally. You review the files and run the upload command yourself. The scripts do not upload files in `--dry-run` mode.

## Prepare

Use Linux with GNU command-line utilities. Use a checkout that contains the DuckVEP extension and publication scripts. Build the extension with `make release`, or select a compatible local artifact. Install the DuckDB command-line client, `sha256sum`, and the Hugging Face `hf` CLI. Benchmark inputs also require `bcftools` and `zcat`.

Make sure that the extension loads in your DuckDB client. Set `DUCKDB` and `DUCKVEP_EXTENSION` if the defaults do not point to the correct files:

```bash
cd /path/to/DuckVEP
export DUCKDB=/path/to/duckdb
export DUCKVEP_EXTENSION="$PWD/build/release/duckvep.duckdb_extension"
"$DUCKDB" -unsigned :memory: -c "LOAD '$DUCKVEP_EXTENSION'; SELECT duckvep_allele_geometry(100, 'A', 'ATG').interbase;"
```

The destination is a Hugging Face **dataset** repository. Confirm that you can redistribute every input. Review the repository visibility before you upload. A public repository stays public when you upload to it.

## Stage model files

Each source database must contain `model_regions`, `model_transcripts`, and one `model_receipt` row. The script also exports `duckvep_regulation_features` and `model_rna_edits` when those tables are present. Set `STAGE` to a new, empty directory. Do not put source files in the stage directory.

```bash
export REPO=RgenomicsETL/duckvep-models
export STAGE="$(mktemp -d)"

bash scripts/publish-models.sh --dry-run "$REPO" \
  /path/to/homo_sapiens_116_GRCh37_final.duckdb \
  /path/to/homo_sapiens_116_GRCh38_final.duckdb
```

You can give one model file or several model files. The script compares native snapshots loaded from each source and export. It writes a dataset card and `SHA256SUMS` in `STAGE`. It needs temporary disk space for the exported database and two snapshots.

The comparison covers the load queries in the generated card. The script includes circular regions, cDNA, and RNA authority fields when the source supplies them. It loads RNA edits from `model_rna_edits`. It does not derive edits from other source tables. It does not supply a reference FASTA. Check that this load profile is appropriate for the intended model before publication.

By default, the generated load query does not assert complete transcript coverage. Add `--complete-coverage` only when the model contains the complete transcript set required for intergenic annotation. The option is a publisher assertion, not a completeness test.

To include GRCh37 MANE files, add `--mane-grch37 /path/to/mane-grch37`. That directory must contain `mane_grch37_mapping.parquet`, `grch37_transcript_authorities.parquet`, `mane_grch37_receipt.csv`, and `SHA256SUMS`. The mapping model hash must match the exported model receipt.

## Review the stage

Check the script output. A failed run is not ready for upload. Start a new stage if you repeat a failed run.

```bash
find "$STAGE" -type f -print | sort
(cd "$STAGE" && sha256sum -c SHA256SUMS)
less "$STAGE/README.md"
(cd "$STAGE" && sha256sum README.md SHA256SUMS)
```

Check the model assemblies, release numbers, dataset card, license, and provenance. Record the `README.md` and `SHA256SUMS` hashes with your release notes. The manifest covers staged data files. It does not cover the card or the manifest itself.

## Upload explicitly

Authenticate with a token that has write access to the destination. Keep the token out of scripts and shell history.

```bash
hf auth login
hf auth whoami
hf repos create "$REPO" --repo-type dataset --private --exist-ok
hf upload "$REPO" "$STAGE" . --repo-type dataset \
  --commit-message "Publish DuckVEP dataset" --exclude '*.wal'
```

The upload adds or replaces files. It does not remove other files from the current repository or from earlier revisions. The checksum check does not identify extra remote files. Review those files and repository visibility separately. `--exist-ok` does not change a repository's visibility.

For an independent download check, select a directory with enough disk space:

```bash
export VERIFY_DIR="$(mktemp -d)"
hf download "$REPO" --repo-type dataset --local-dir "$VERIFY_DIR"
(cd "$VERIFY_DIR" && sha256sum -c SHA256SUMS)
sha256sum "$STAGE/README.md" "$VERIFY_DIR/README.md"
```

Check that the two card hashes agree. Record the Hugging Face commit or revision with the release evidence.

## Stage benchmark inputs

The benchmark publisher expects the specific GIAB and gnomAD-derived inputs described in `scripts/benchmark-inputs-card.md`. Check their identities and provenance before publication. Set the input paths and use a different empty stage:

```bash
export REPO=RgenomicsETL/duckvep-benchmark-inputs
export STAGE="$(mktemp -d)"
export HG002_VCF=/path/to/hg002.ens.vcf.gz
export QUAL5M_VCF=/path/to/qual5m.vcf.gz
export CONTROLS_DIR=/path/to/controls
export CORPUS_DB=/path/to/hg002_corpus.duckdb
export MODEL_DB=/path/to/model.duckdb
export CORPUS_VCF=/path/to/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

bash scripts/publish-benchmark-inputs.sh --dry-run "$REPO"
(cd "$STAGE" && sha256sum -c SHA256SUMS)
less "$STAGE/README.md"
```

The controls directory must contain `dense.vcf.gz` and `low_sharing.vcf.gz`. Each VCF must have a `.tbi` index. The script compares Parquet rows with an independent VCF parse and checks the throughput corpus. Exact comparisons can require substantial time, memory, and temporary disk space.

To run the small local publication tests, install `bgzip` and `tabix`, then run `make test-publication`. These tests use temporary fixtures and do not upload files.

Review the staged card and source provenance before you run the upload commands above. The scripts do not check that an arbitrary replacement input has the GIAB or gnomAD provenance in the card.
