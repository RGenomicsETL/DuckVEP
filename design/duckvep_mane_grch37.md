# MANE v1.5 mapped to GRCh37

Status: data release
[`mane-grch37-v1.5-e116`](https://github.com/RGenomicsETL/DuckVEP/releases/tag/mane-grch37-v1.5-e116),
published 2026-09-28 (MANE v1.5 against the Ensembl 116 GRCh37 model, #6). This is a data
release, not a software release; the extension and R package ship separately.

## What this relation is

Current MANE is defined only on GRCh38; a GRCh37 model has no native MANE data (see
[the GRCh37 Ensembl site](https://grch37.ensembl.org/index.html) and
[Ensembl's MANE description](https://www.ensembl.org/info/genome/genebuild/mane.html)).
This release is an external association that audits MANE v1.5 RefSeq transcripts against
DuckVEP's Ensembl 116 GRCh37 core model. Every row produced by it must be labelled
**"MANE mapped to GRCh37,"** never native MANE.

Two cold Parquet relations ship in the release, built by
[`scripts/build_mane_grch37.R`](../scripts/build_mane_grch37.R):

- `mane_grch37_mapping.parquet` — one row per MANE v1.5 record (19,437 rows). It carries the
  exact RefSeq and Ensembl accessions, the GRCh37.p13 target geometry, the per-gate
  validation results, a reason-coded `mapping_status`, and the `transcript_index` of the
  model transcript it applies to.
- `grch37_transcript_authorities.parquet` — native GRCh37 facts for every filtered model
  transcript: Ensembl canonical, GENCODE-19 Basic, and `no_retained_canonical`.
  `gencode_primary` is always `false`, because GENCODE Primary isn't defined for GRCh37.

Neither relation is loaded into the DuckVEP consequence model. They never change native
consequence masks, ordinals, HGVS, or VEP-parity results; a caller joins them onto already
expanded consequence rows. See
[the design contract's GRCh37 section](duckvep.md#grch37-transcript-selection-and-external-mane-mapping)
for the full construction and validation policy.

## Status tiers

| `mapping_status`         |  rows | what is claimed |
|:--------------------------|------:|:-----------------|
| `exact_model_match`       |   740 | The full transcript is identical: exon chain, CDS and phase, spliced sequence, and translation. |
| `cds_exact_utr_differs`   | 15,241 | **Coding region only.** CDS segments and phases are identical, the translation equals the RefSeq protein, and there are no reference differences within the CDS. Coding consequences and c./p. HGVS are equivalent. UTR `c.-N`/`c.*N` and `n.` positions, and UTR or non-coding exon consequences, are **not** claimed. |
| all other statuses        | 3,456 | No association; each row keeps its reason code (`ambiguous_target_locus`, `cds_phase_mismatch`, `geometry_mismatch`, `refseq_only_no_gencode19_match`, `sequence_mismatch`, `target_reference_unavailable`, `target_transcript_absent`, `translation_mismatch`). |

The 3,456 unmapped rows break down as: `refseq_only_no_gencode19_match` 2,020,
`geometry_mismatch` 1,281, `target_transcript_absent` 70, `cds_phase_mismatch` 23,
`target_reference_unavailable` 38, `ambiguous_target_locus` 8, `sequence_mismatch` 12, and
`translation_mismatch` 4, per the committed
[receipt ledger](../benchmarks/data/mane_grch37_receipts.csv).

## The model this release applies to

The mapping was validated against, and applies only to, the Ensembl 116 GRCh37 model with

```text
model_sha256 = 21e113d9148132491bc935f3d1b0ec7d50663f450b62e346cb1f0f447de0b290
```

That is the deterministic hash `query(duckvep_model_receipt_sql(...))` computes over the model's
declared source, transcripts, exons, and reference sequence. A caller joins on
`model_sha256` together with `transcript_index`; a model built from a different core dump or
reference FASTA has a different `model_sha256`, was not validated against this release, and
the join below simply returns no rows for it rather than a silent mismatch.

## Caller recipe

[`scripts/mane_grch37_caller.sql`](../scripts/mane_grch37_caller.sql) is the reference join. It runs *after* consequence
expansion, aggregates each tier into a list per `(model_sha256, transcript_index)` instead of
duplicating consequence rows one-to-many, and joins the native GRCh37 transcript-authority
facts alongside it. Replace `OUTPUT_DIR` with the directory holding the downloaded release
Parquet files.

```sql
-- Run after expanding consequences against the native GRCh37 model.
-- Replace OUTPUT_DIR with the receipt-matched external Parquet directory.
-- consequence_rows has model_sha256, transcript_index and caller event columns.
WITH mapped AS (
  SELECT model_sha256, transcript_index,
         list(struct_pack(status := mane_status, refseq_nuc := refseq_nuc,
                          refseq_prot := refseq_prot, mapping_label := mapping_label,
                          source_digest := mane_sha256, target_assembly := target_assembly)
              ORDER BY mane_status, refseq_nuc)
           FILTER (WHERE mapping_status = 'exact_model_match') AS mane_mapped_to_grch37,
         list(struct_pack(status := mane_status, refseq_nuc := refseq_nuc,
                          refseq_prot := refseq_prot, mapping_label := mapping_label,
                          source_digest := mane_sha256, target_assembly := target_assembly)
              ORDER BY mane_status, refseq_nuc)
           FILTER (WHERE mapping_status = 'cds_exact_utr_differs') AS mane_coding_region_only
  FROM read_parquet('OUTPUT_DIR/mane_grch37_mapping.parquet')
  WHERE mapping_status IN ('exact_model_match', 'cds_exact_utr_differs')
  GROUP BY model_sha256, transcript_index
)
SELECT c.*, n.canonical AS ensembl_grch37_canonical,
       n.gencode_basic AS gencode19_basic,
       n.no_retained_canonical, n.model_sha256 AS native_model_sha256,
       mapped.mane_mapped_to_grch37, mapped.mane_coding_region_only
FROM consequence_rows AS c
JOIN read_parquet('OUTPUT_DIR/grch37_transcript_authorities.parquet') AS n
  ON n.model_sha256 = c.model_sha256 AND n.transcript_index = c.transcript_index
LEFT JOIN mapped
  ON mapped.model_sha256 = c.model_sha256 AND mapped.transcript_index = c.transcript_index;
```

`mane_coding_region_only` is usable only for coding consequences or coding c./p. HGVS; do
not read UTR or non-coding positions from it. If several MANE rows reference one transcript
they arrive as a list rather than duplicated consequence rows. Canonical, GENCODE-19 Basic,
and mapped MANE are distinct source-attributed facts — selecting one representative
transcript remains an explicit caller policy, not VEP `--pick`.

## Download and verify

| Asset | SHA-256 |
|:------|:--------|
| [`mane_grch37_mapping.parquet`](https://github.com/RGenomicsETL/DuckVEP/releases/download/mane-grch37-v1.5-e116/mane_grch37_mapping.parquet) | `dc6c369c3b61cc4d232bb7aa84e8dce736de544b3f4c3c5cac49b9cf1aab3861` |
| [`grch37_transcript_authorities.parquet`](https://github.com/RGenomicsETL/DuckVEP/releases/download/mane-grch37-v1.5-e116/grch37_transcript_authorities.parquet) | `e893d63b9a878d9c60c3b704184ff9edf43868273746c75bf418999cfa34687b` |
| [`mane_grch37_receipt.csv`](https://github.com/RGenomicsETL/DuckVEP/releases/download/mane-grch37-v1.5-e116/mane_grch37_receipt.csv) | `3e444e9007e824692a533e4733f4eb96d47d0fc4c0c8b51af714ca0a27e1131f` |
| [`SHA256SUMS`](https://github.com/RGenomicsETL/DuckVEP/releases/download/mane-grch37-v1.5-e116/SHA256SUMS) | (lists the three digests above) |

None of these Parquet files are committed to the repository; they are external release
artifacts, referenced here by their release download links and checksums only. After
downloading `SHA256SUMS` alongside the other three files into the same directory, verify
them with:

```sh
sha256sum -c SHA256SUMS
```

## Rebuilding

The release is reproducible from pinned, checksum-verified sources, not redistributed
upstream data. From a checkout of this repository:

```sh
sh scripts/stage_mane_grch37.sh STAGING_DIR
Rscript scripts/build_mane_grch37.R STAGING_DIR MODEL.duckdb OUTPUT_DIR
```

`scripts/stage_mane_grch37.sh` downloads the MANE v1.5 summary, the NCBI GRCh37.p13
assembly report/GFF/RNA/protein FASTAs, and the Ensembl GRCh37 primary-assembly FASTA,
checking each against a pinned SHA-256 before use; it requires `curl` and `samtools`.
`scripts/build_mane_grch37.R` then builds the association against `MODEL.duckdb` (an
Ensembl 116 GRCh37 core model with the `model_sha256` above) and writes the two Parquet
relations plus `mane_grch37_receipt.csv` into `OUTPUT_DIR`; it requires the R packages
`DBI`, `duckdb`, `data.table`, `Biostrings`, `Rsamtools`, `GenomicRanges`, and `digest`. The
script uses plain R FASTA readers and DuckDB SQL — it loads neither DuckHTS nor DuckVEP into
its DuckDB connections, and all network access stays in the staging step.

`make test_mane_grch37` runs the offline policy fixtures without any staged data.
`make test_mane_grch37 MANE_GRCH37_OUTPUT=OUTPUT_DIR` additionally checks a rebuilt
`OUTPUT_DIR` against the last row of the checked-in
[receipt ledger](../benchmarks/data/mane_grch37_receipts.csv),
including the relation checksum `4a1f3a03e9af72e35dad742636c1c9e4c8a0e0e3b03043d921aa34ec1777dace`
recorded there for this release.
