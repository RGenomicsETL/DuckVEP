# Function reference

DuckVEP registers 25 public SQL functions: 21 scalar functions and 4 table functions. This page lists all of them, grouped by purpose. Internal helpers whose names start with `_duckvep_` or `__duckvep_` are implementation details and are not documented. `scripts/check-function-docs.py` (`make check-function-docs`) fails when this page misses or adds a public function relative to `duckdb_functions()`, and it runs every example below.

Every `sql` example on this page runs against the fixture model in `test/data/duckvep/readme.sql`, in the order it appears on the page: later examples use tables created by earlier ones. Blocks marked `sql no-run` are illustrative fragments or need external data.

## Conventions

**Builders return SQL text.** Most of the functions ending in `_sql` are scalar functions that return a SQL string, not rows. The string is executed in the caller's connection with `query()`, so the relations you name stay visible, including temporary tables and uncommitted rows:

```sql no-run
FROM query(duckvep_annotate_sql('demo_events', 'demo', {hgvs: true}))
```

**Relation names are strings.** The builders quote them: a name may be qualified with one schema (`schema.table`), and quotes inside a name are escaped. An empty name, a NULL name or an invalid qualified name is an error.

**Options go in a trailing STRUCT.** Optional settings are the fields of one final STRUCT argument (`{hgvs: true}` or `struct_pack(hgvs := true)`), never named function parameters: the stable DuckDB C API used by this extension has no named scalar arguments, so `hgvs := true` written as a bare argument does not bind. An unknown field, or a field with the wrong type, is an error that names the option. A field set to NULL keeps the default. A builder that takes no options rejects any field.

**Table functions take named parameters.** `duckvep_model_load` and `duckvep_haplotypes` are table functions and use ordinary `name := value` parameters.

**Unknown is a value.** When a result cannot be computed (a missing reference sequence, a reference mismatch, an unsupported allele) the row carries a status and a reason instead of a guess.

**Positions** are one-based unless a column name ends in `0`, which marks a zero-based coordinate.

## Index

| Section | Function | Kind | Purpose |
| --- | --- | --- | --- |
| [Model: build and load](#model-build-and-load) | [`duckvep_ensembl_regions_sql`](#duckvep_ensembl_regions_sql) | scalar | SQL that turns an Ensembl core schema and a tiled FASTA into the regions relation. |
| | [`duckvep_ensembl_transcripts_sql`](#duckvep_ensembl_transcripts_sql) | scalar | SQL that builds the transcripts relation, with CDS and flank sequence and exons. |
| | [`duckvep_ensembl_regulation_features_sql`](#duckvep_ensembl_regulation_features_sql) | scalar | SQL that builds regulatory and motif features from an Ensembl funcgen schema. |
| | [`duckvep_model_receipt_sql`](#duckvep_model_receipt_sql) | scalar | SQL that summarizes a prepared model with counts and content hashes. |
| | [`duckvep_model_load`](#duckvep_model_load) | table | Compile relations into a named, immutable, resident model. |
| | [`duckvep_model_drop`](#duckvep_model_drop) | scalar | Unload a model and return its memory. |
| [Annotation: SQL builders](#annotation-sql-builders) | [`duckvep_annotate_sql`](#duckvep_annotate_sql) | scalar | SQL that annotates a relation of alleles with consequences, impact and HGVS. |
| | [`duckvep_annotate_projected_sql`](#duckvep_annotate_projected_sql) | scalar | SQL that annotates ordered small variants and returns projected-edit facts. |
| | [`duckvep_transcript_projection_sql`](#duckvep_transcript_projection_sql) | scalar | SQL that presents transcript positions, codons and peptides for annotated events. |
| [Structural and repeat preparation](#structural-and-repeat-preparation) | [`duckvep_prepare_sv_geometry_sql`](#duckvep_prepare_sv_geometry_sql) | scalar | SQL that normalizes VCF structural geometry and insertion provenance. |
| | [`duckvep_prepare_breakend_pairs_sql`](#duckvep_prepare_breakend_pairs_sql) | scalar | SQL that validates BND mate and event identity. |
| | [`duckvep_prepare_breakend_fusion_sql`](#duckvep_prepare_breakend_fusion_sql) | scalar | SQL that joins BND identity with endpoint genes. |
| | [`duckvep_prepare_structural_hgvs_sql`](#duckvep_prepare_structural_hgvs_sql) | scalar | SQL that describes exact-span DEL, DUP and INV alleles as genomic HGVS. |
| | [`duckvep_prepare_expansionhunter_sql`](#duckvep_prepare_expansionhunter_sql) | scalar | SQL that prepares ExpansionHunter v5 repeat alleles. |
| | [`duckvep_repeat_alleles`](#duckvep_repeat_alleles) | scalar | Expand ordered repeat units and counts into reference and alternate sequences. |
| [Haplotypes](#haplotypes) | [`duckvep_phase_call`](#duckvep_phase_call) | scalar | Assign genotype slots to haplotype lanes and phase sets. |
| | [`duckvep_coding_transcripts`](#duckvep_coding_transcripts) | scalar | List the transcripts whose coding sequence a VCF record overlaps. |
| | [`duckvep_haplotypes`](#duckvep_haplotypes) | table | Replay phased calls into whole-haplotype CDS, protein and consequence rows. |
| [Geometry and helpers](#geometry-and-helpers) | [`duckvep_allele_geometry`](#duckvep_allele_geometry) | scalar | Normalized coordinates of one small allele. |
| | [`duckvep_breakend_geometry`](#duckvep_breakend_geometry) | scalar | Parse a breakend ALT into its mate coordinate and replacement sequence. |
| [Resource control](#resource-control) | [`duckvep_native_budget`](#duckvep_native_budget) | table | Report native memory use per owner. |
| | [`duckvep_native_budget_set`](#duckvep_native_budget_set) | scalar | Set the process-wide native memory ceiling. |
| | [`duckvep_native_budget_reset_high_water`](#duckvep_native_budget_reset_high_water) | scalar | Restart the high-water marks. |
| | [`duckvep_worker_limits_set`](#duckvep_worker_limits_set) | scalar | Set the annotation worker count and per-worker byte limits. |
| [Vocabulary](#vocabulary) | [`duckvep_so_terms`](#duckvep_so_terms) | table | The Sequence Ontology terms, bit positions, impacts and ranks. |

---

## Model: build and load

An Ensembl release is compiled once into an immutable, resident transcript model that annotation functions then read by name. A model is loaded from ordinary relations: regions, transcripts (with their CDS and flanking sequence) and exons. In production the relations come from an Ensembl core database through the three `duckvep_ensembl_*_sql` builders; the examples here use the one-transcript fixture from `test/data/duckvep/readme.sql`.

<a id="duckvep_ensembl_regions_sql"></a>

### duckvep_ensembl_regions_sql

Builds the SQL that produces the region relation for an assembly from an Ensembl core schema and a table of reference sequence chunks.

Signatures:

```text
duckvep_ensembl_regions_sql(core_schema VARCHAR, reference_chunks_table VARCHAR, assembly VARCHAR) -> VARCHAR
duckvep_ensembl_regions_sql(core_schema VARCHAR, reference_chunks_table VARCHAR, assembly VARCHAR, options STRUCT) -> VARCHAR
```

Parameters:

| Parameter | Description |
| --- | --- |
| `core_schema` | Schema holding the Ensembl core tables. |
| `reference_chunks_table` | Relation with columns `chrom`, `start`, `end` and `seq`, for example from `fasta_nuc(..., include_seq := true)`. |
| `assembly` | Assembly name, for example `'GRCh38'`. |

Options: `species_id` (INTEGER, default 1) selects the species in a multi-species core schema.

Returns: a SQL string. The relation it produces has the columns `seq_region` UINTEGER (dense region ordinal), `seq_region_name` VARCHAR, `sequence_length` UBIGINT, `circular` BOOLEAN, `source_seq_region_id` BIGINT, `source_coord_system_id` BIGINT, `coord_system_name` VARCHAR and `coord_system_rank` BIGINT. The Ensembl compiler requires each region to match a reference region of the same name and length.

```sql
CREATE TABLE prepared_regions AS
SELECT * FROM query(duckvep_ensembl_regions_sql(
  'duckvep_core', 'duckvep_reference_chunks', 'GRCh38'));
```

<a id="duckvep_ensembl_transcripts_sql"></a>

### duckvep_ensembl_transcripts_sql

Builds the SQL that produces the transcript relation, including CDS, pre- and post-CDS transcript flank sequence and the exon list, from an Ensembl core schema and reference chunks.

Signatures:

```text
duckvep_ensembl_transcripts_sql(core_schema VARCHAR, reference_chunks_table VARCHAR, assembly VARCHAR) -> VARCHAR
duckvep_ensembl_transcripts_sql(core_schema VARCHAR, reference_chunks_table VARCHAR, assembly VARCHAR, options STRUCT) -> VARCHAR
```

Parameters are those of [`duckvep_ensembl_regions_sql`](#duckvep_ensembl_regions_sql). Options: `species_id` (INTEGER, default 1).

Returns: a SQL string. The relation it produces carries the model transcript columns (`transcript_index`, `seq_region`, `transcript_start`, `transcript_end`, `strand`, `gene_index`, `transcript_flags`, `cds_start`, `cds_end`, `cds_sequence`, `codon_table`, `pre_cds_sequence`, `post_cds_sequence`) followed by source identity and metadata columns (`seq_region_name`, `sequence_length`, `source_seq_region_id`, `source_transcript_id`, `transcript_stable_id`, `transcript_version`, `transcript_biotype`, `mane_select_refseq`, `mane_plus_clinical_refseq`, `source_gene_id`, `gene_stable_id`, `gene_version`, `gene_biotype`, `source_translation_id`, `translation_stable_id`, `translation_version`, `sequence_withheld_reason`), the list columns `exons`, `mature_mirna_regions` and `peptide_edits`, and `circular` and `origin_crossing` flags. Stable identifiers stay DuckDB columns; the resident model holds only numeric ordinals.

```sql
CREATE TABLE prepared_transcripts AS
SELECT * FROM query(duckvep_ensembl_transcripts_sql(
  'duckvep_core', 'duckvep_reference_chunks', 'GRCh38'));
```

<a id="duckvep_ensembl_regulation_features_sql"></a>

### duckvep_ensembl_regulation_features_sql

Builds the SQL that produces regulatory-feature and motif-feature rows from an Ensembl funcgen schema, aligned to a prepared regions relation. Epigenetically modified regions (EMAR) are removed, as VEP 116 does before it builds regulatory overlap objects.

Signatures:

```text
duckvep_ensembl_regulation_features_sql(funcgen_schema VARCHAR, regions_table VARCHAR) -> VARCHAR
duckvep_ensembl_regulation_features_sql(funcgen_schema VARCHAR, regions_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `funcgen_schema` is the Ensembl funcgen schema and `regions_table` the prepared regions relation. This builder accepts no options; a non-empty STRUCT is an error.

Returns: a SQL string. The relation it produces has the hot columns `regulation_feature_index` UINTEGER, `seq_region` UINTEGER, `feature_start` UINTEGER, `feature_end` UINTEGER, `strand` TINYINT and `feature_kind` UTINYINT (1 for a RegulatoryFeature, 2 for a MotifFeature), and metadata columns `feature_class`, `source_feature_table`, `source_feature_id`, `stable_id`, `source_feature_type_id`, `source_binding_matrix_id`, `source_regulatory_build_id`, `feature_name`, `feature_so_accession`, `feature_so_term`, `score`, `circular` and `origin_crossing`. Only the hot columns enter the resident model, through `interval_feature_query` of [`duckvep_model_load`](#duckvep_model_load).

```sql
CREATE TABLE prepared_regulation AS
SELECT * FROM query(duckvep_ensembl_regulation_features_sql(
  'duckvep_funcgen', 'prepared_regions'));
```

<a id="duckvep_model_receipt_sql"></a>

### duckvep_model_receipt_sql

Builds the SQL that summarizes a prepared model as one receipt row: source identity, counts and content hashes. A reused model can be checked by reproducing its stored receipt.

Signatures:

```text
duckvep_model_receipt_sql(regions_table, transcripts_table, source_name, source_version, assembly,
                          source_manifest_sha256, reference_sha256, transcript_filter) -> VARCHAR
duckvep_model_receipt_sql(regions_table, transcripts_table, source_name, source_version, assembly,
                          source_manifest_sha256, reference_sha256, transcript_filter, options STRUCT) -> VARCHAR
```

All eight required parameters are VARCHAR. `regions_table` and `transcripts_table` name the prepared relations. `source_name`, `source_version` and `assembly` identify the source (for example `'Ensembl'`, `'116'`, `'GRCh38'`). `source_manifest_sha256` and `reference_sha256` are digests supplied by the caller, and `transcript_filter` describes the transcript selection.

Options: `regulation_features_table` (VARCHAR) names a prepared regulation relation to include in the counts.

Returns: a SQL string. The relation it produces has one row with `source_name`, `source_version`, `assembly`, `source_manifest_sha256`, `reference_sha256`, `transcript_filter`, `model_sha256`, `region_count`, `reference_base_count`, `transcript_count`, `gene_count`, `coding_transcript_count`, `sequence_backed_transcript_count`, `sequence_withheld_transcript_count`, `exon_membership_count`, `cds_base_count`, `mature_mirna_transcript_count`, `mature_mirna_segment_count`, `peptide_edit_count`, `transcript_flank_base_count`, `regulation_feature_count`, `regulatory_region_count`, `motif_feature_count`, `circular_region_count`, `circular_regions` and `topology_sha256`.

```sql
SELECT region_count, transcript_count, gene_count, model_sha256
FROM query(duckvep_model_receipt_sql(
  'prepared_regions', 'prepared_transcripts',
  'Ensembl', '116', 'GRCh38',
  repeat('a', 64), repeat('b', 64), 'fixture transcripts',
  {regulation_features_table: 'prepared_regulation'}));
```

<a id="duckvep_model_load"></a>

### duckvep_model_load

Reads relations through a private connection, validates and narrows every value, builds independent transcript and regulation/motif interval indexes, and only then publishes the named immutable model. A failed load publishes nothing. Several models can coexist in one database instance. Ordinals are meaningful only within their model.

Signature:

```text
duckvep_model_load(name VARCHAR, regions_query VARCHAR, transcripts_query VARCHAR, exons_query VARCHAR
                   [, transcript_coverage_complete := BOOLEAN]
                   [, interval_feature_query := VARCHAR]
                   [, peptide_edit_query := VARCHAR]
                   [, mature_mirna_query := VARCHAR]
                   [, reference_fasta := VARCHAR]) -> TABLE(loaded BOOLEAN)
```

Parameters (positional):

| Parameter | Description |
| --- | --- |
| `name` | Non-empty model name. Loading an existing name is an error. |
| `regions_query` | A SELECT returning `seq_region`, and optionally `sequence_length`, `seq_region_name` and `circular`. |
| `transcripts_query` | A SELECT returning either the 11 CDS-only columns or the 13-column form ending in `pre_cds_sequence` and `post_cds_sequence`. Only the complete form can resolve length-changing edits that cross the CDS start or end; otherwise those return `missing_transcript_flank`. |
| `exons_query` | A SELECT returning `transcript_index`, `exon_start`, `exon_end`, `exon_cdna_start`, `exon_cdna_end`, `phase` and `end_phase`. |

Named parameters:

| Name | Type | Description |
| --- | --- | --- |
| `transcript_coverage_complete` | BOOLEAN | Default false. Only a model loaded with true reports `intergenic_variant` where no loaded transcript lies; a partial model returns an unresolved result there. NULL is an error. |
| `interval_feature_query` | VARCHAR | Regulatory and motif features: feature ordinal, region ordinal, inclusive start, inclusive end and kind (1 RegulatoryFeature, 2 MotifFeature), ordered by region, start and ordinal. |
| `peptide_edit_query` | VARCHAR | Curated peptide edits: transcript ordinal, one-based protein position and an uppercase replacement amino acid, unique and ordered by transcript and position. |
| `mature_mirna_query` | VARCHAR | Mature miRNA ranges: transcript ordinal, inclusive genomic start and end, ordered by transcript and start. |
| `reference_fasta` | VARCHAR | Path to an indexed reference FASTA. Needed to shift and render HGVS notation. Requires `seq_region`, `sequence_length` and `seq_region_name` in the region query. |

A NULL value for an optional query or the FASTA path is the same as omitting it. Without a reference FASTA, HGVS fields that need flanking sequence come back with a `missing_reference` reason (see [`duckvep_annotate_sql`](#duckvep_annotate_sql)).

Returns: one row, `loaded BOOLEAN`. A load that exceeds the native budget publishes nothing and leaves every loaded model usable.

```sql
SELECT loaded FROM duckvep_model_load('demo',
  'SELECT * FROM readme_regions ORDER BY seq_region',
  'SELECT * FROM readme_transcripts ORDER BY seq_region, transcript_start',
  'SELECT * FROM readme_exons ORDER BY transcript_index, exon_start');
```

<a id="duckvep_model_drop"></a>

### duckvep_model_drop

Unloads a named model and returns its model and index bytes to the native budget.

Signature:

```text
duckvep_model_drop(name VARCHAR) -> BOOLEAN
```

Returns true when the model was dropped. Returns false when no model has that name or when the model is pinned by a running query. An empty name is an error.

```sql
SELECT loaded FROM duckvep_model_load('scratch',
  'SELECT * FROM readme_regions ORDER BY seq_region',
  'SELECT * FROM readme_transcripts ORDER BY seq_region, transcript_start',
  'SELECT * FROM readme_exons ORDER BY transcript_index, exon_start');
SELECT duckvep_model_drop('scratch') AS dropped,
       duckvep_model_drop('scratch') AS dropped_again;
```

---

## Annotation: SQL builders

Annotation runs as SQL emitted by a builder and executed with `query()`. The event relation stays in the caller's transaction. The three builders below share one events relation.

**Events relation.** One row per ALT allele, already in global coordinate order (the builder never hides an `ORDER BY`; input that is out of order is rejected by the native sorted-stream checks) with the columns:

| Column | Description |
| --- | --- |
| `event_index` | Unsigned identity, used to join results back to the source rows. |
| `seq_region` | Model-local region ordinal. |
| `position` | One-based position. |
| `reference`, `alternate` | Literal REF and ALT alleles. |
| `end_position`, `structural_type`, `copy_change` | Nullable single-locus structural span, type and copy direction. |
| `mate_seq_region`, `mate_position` | Nullable mate coordinates for a paired breakend. |

The relation is classified as a small variant, an exact structural event or a paired breakend, and dispatched to the matching native lane. Wide provenance (genotypes, raw ALT, confidence intervals) stays in the caller's table and is joined back by `event_index`.

The examples use seven variants across the fixture transcript, from the 5' UTR to the intron:

```sql
CREATE TABLE demo_events AS
SELECT row_number() OVER (ORDER BY position, alternate)::UBIGINT AS event_index,
       1::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate,
       NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM (VALUES (110, 'C', 'T'), (124, 'T', 'C'), (125, 'A', 'G'), (130, 'C', 'CA'),
             (134, 'C', 'A'), (151, 'G', 'A'), (175, 'A', 'G')) v(position, reference, alternate);
```

<a id="duckvep_annotate_sql"></a>

### duckvep_annotate_sql

Builds the SQL that annotates every allele of an events relation against a loaded model, one row per allele and overlapping transcript or feature. The result is typed columns, not text: Sequence Ontology consequences, impact, cDNA/CDS/protein positions, amino acids, nonsense-mediated decay (NMD) prediction and, on request, HGVS.

Signatures:

```text
duckvep_annotate_sql(events_table VARCHAR, model_name VARCHAR) -> VARCHAR
duckvep_annotate_sql(events_table VARCHAR, model_name VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `events_table` names the events relation and `model_name` a loaded model. A NULL argument is an error.

Options:

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `hgvs` | BOOLEAN | false | Also compute transcript and protein HGVS. |
| `rich` | BOOLEAN | false | Fill the decoded text columns (`consequence`, `impact`, `region`, amino acids, `nmd_prediction`, `overlap_object`, `duckvep_status`, `duckvep_reason`). The compact form leaves them NULL and reports numeric codes only. |
| `upstream_distance` | INTEGER | 5000 | Bases upstream of a transcript in which it is still a candidate. 0 disables the direction. |
| `downstream_distance` | INTEGER | 5000 | The same, downstream. |

The distances extend candidate admission beyond transcript endpoints. They do not cap an allele span or clip an event that overlaps a transcript.

Returns: the same fixed columns for all three event families.

| Columns | Description |
| --- | --- |
| `event_index` UBIGINT, `duckvep_event_kind` VARCHAR | Event identity and event family (`small_variant` for literal small alleles). |
| `transcript_index` UINTEGER, `gene_index` UINTEGER | Overlapped transcript and gene ordinals. |
| `consequence_mask` UBIGINT, `region_mask` UINTEGER, `impact_code` UTINYINT, `status_code` UTINYINT, `reason_code` UTINYINT | Numeric results. Decode `consequence_mask` with [`duckvep_so_terms`](#duckvep_so_terms). |
| `cdna_position`, `cds_position`, `protein_position` UINTEGER | Projected positions. |
| `reference_amino_acid_code`, `alternate_amino_acid_code` UTINYINT | Numeric amino acids. |
| `nmd_prediction_code` UTINYINT, `nmd_escape_reasons` UTINYINT | Numeric NMD result and escape reasons. |
| `regulation_feature_index` UINTEGER, `overlap_object_code` UTINYINT | Overlapped regulation or motif feature, when the model has one. |
| `transcript_hgvs`, `protein_hgvs` VARCHAR, `hgvs_shift` UINTEGER | HGVS text and shift, when `hgvs` is true. |
| `transcript_hgvs_status`, `transcript_hgvs_reason`, `protein_hgvs_status`, `protein_hgvs_reason` VARCHAR | Status and reason of each HGVS string. |
| `consequence`, `impact`, `region`, `reference_amino_acid`, `alternate_amino_acid`, `nmd_prediction` VARCHAR | Decoded text (with `rich`). |
| `nmd_escape_intronless`, `nmd_escape_early_cds`, `nmd_escape_last_exon`, `nmd_escape_penultimate_exon_end` BOOLEAN | Each NMD escape reason as a boolean. |
| `overlap_object` VARCHAR | Decoded overlap object (with `rich`). |
| `duckvep_status`, `duckvep_reason` VARCHAR | Whether the result was computed and why not (with `rich`). |

Join results to the source table by `event_index`. An HGVS string that could not be rendered is NULL and its status and reason say why, for example `unresolved` with `missing_reference` when the model was loaded without a reference FASTA.

```sql
SELECT e.position, e.reference || '>' || e.alternate AS change,
       a.consequence, a.impact, a.transcript_hgvs, a.protein_hgvs, a.nmd_prediction
FROM query(duckvep_annotate_sql('demo_events', 'demo', {rich: true, hgvs: true})) a
JOIN demo_events e USING (event_index)
ORDER BY e.position;
```

The compact form reports masks; `duckvep_so_terms()` decodes them:

```sql
SELECT a.event_index, list(s.consequence ORDER BY s.severity_rank) AS consequences
FROM query(duckvep_annotate_sql('demo_events', 'demo')) a
JOIN duckvep_so_terms() s ON (a.consequence_mask & s.consequence_mask) <> 0
GROUP BY a.event_index
ORDER BY a.event_index;
```

Without a reference FASTA the HGVS fields say why they are empty:

```sql
SELECT e.position, a.transcript_hgvs_status, a.transcript_hgvs_reason
FROM query(duckvep_annotate_sql('demo_events', 'demo', {hgvs: true})) a
JOIN demo_events e USING (event_index)
WHERE a.transcript_hgvs IS NULL
ORDER BY e.position;
```

An illustrative production query, with a VCF read by DuckHTS, needs external data and is not run:

```sql no-run
CREATE TABLE alleles AS
SELECT row_number() OVER (ORDER BY r.seq_region, v.POS, alt_allele)::UBIGINT AS event_index,
       r.seq_region, v.POS::UBIGINT AS position, v.REF AS reference, alt_allele AS alternate,
       NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM read_bcf('cohort.vcf.gz') v CROSS JOIN unnest(v.ALT) AS t(alt_allele)
JOIN grch38_regions r ON r.seq_region_name = v.CHROM;

SELECT * FROM query(duckvep_annotate_sql('alleles', 'human_116_grch38', {rich: true, hgvs: true}));
```

<a id="duckvep_annotate_projected_sql"></a>

### duckvep_annotate_projected_sql

Builds the SQL for a dedicated presentation contract: ordered, literal small-variant events are annotated once and returned together with their source columns, the normalized allele geometry and the projected-edit facts, in one relation. It keeps the projected presentation separate from the compact and rich relations of [`duckvep_annotate_sql`](#duckvep_annotate_sql) without a second projection pass.

Signatures:

```text
duckvep_annotate_projected_sql(events_table VARCHAR, model_name VARCHAR) -> VARCHAR
duckvep_annotate_projected_sql(events_table VARCHAR, model_name VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: as for `duckvep_annotate_sql`. The model name must be non-empty. Options: `upstream_distance` and `downstream_distance` (INTEGER, default 5000). Any other field, including `hgvs`, is an error.

Returns: the events columns, `geometry` (the STRUCT of [`duckvep_allele_geometry`](#duckvep_allele_geometry)), `transcript_index`, `gene_index`, the decoded text columns `consequence`, `impact`, `region`, `status`, `reason`, `reference_amino_acid`, `alternate_amino_acid`, `nmd_prediction`, `overlap_object`, the four `nmd_escape_*` booleans, `cdna_position`, `cds_position`, `protein_position`, `regulation_feature_index`, the numeric codes of `duckvep_annotate_sql`, and the projected-edit columns `output_allele`, `interbase`, `cdna_start`, `cdna_end`, `cds_start`, `cds_end`, `protein_start`, `protein_end`, `exon_first`, `exon_last`, `exon_total`, `intron_first`, `intron_last`, `intron_total`, `transcript_distance`, `cds_start_nf`, `cds_end_nf`, `reference_amino_acids`, `alternate_amino_acids`, `reference_codons` and `alternate_codons`.

```sql
SELECT event_index, consequence, impact, status, output_allele, cds_start
FROM query(duckvep_annotate_projected_sql('demo_events', 'demo', {upstream_distance: 0}))
ORDER BY event_index;
```

<a id="duckvep_transcript_projection_sql"></a>

### duckvep_transcript_projection_sql

Builds the SQL that presents, for events already annotated, the transcript-level coordinates, amino acids and codons as separate typed columns, from the events, their annotation relation and a transcript relation. This is the independent reference presentation for the fields that `duckvep_annotate_sql` reports as numeric codes.

Signatures:

```text
duckvep_transcript_projection_sql(events_table VARCHAR, annotations_table VARCHAR, transcripts_table VARCHAR) -> VARCHAR
duckvep_transcript_projection_sql(events_table VARCHAR, annotations_table VARCHAR, transcripts_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `events_table` is the events relation, `annotations_table` a materialized `duckvep_annotate_sql` result (it needs `event_index`, `transcript_index` and `consequence_mask`), and `transcripts_table` a transcript relation with the model transcript columns plus the list columns `exons` and `peptide_edits`, as [`duckvep_ensembl_transcripts_sql`](#duckvep_ensembl_transcripts_sql) produces. Names must be non-empty. The builder accepts no options.

Returns: `event_index`, `transcript_index`, `output_allele`, `interbase`, `cdna_start`, `cdna_end`, `cds_start`, `cds_end`, `protein_start`, `protein_end`, `exon_first`, `exon_last`, `exon_total`, `intron_first`, `intron_last`, `intron_total`, `transcript_distance`, `cds_start_nf`, `cds_end_nf`, `reference_amino_acids`, `alternate_amino_acids`, `reference_codons` and `alternate_codons`. Codon and peptide fields are strings; an insertion or deletion can give empty or partial values, and events outside the coding sequence give NULL.

```sql
CREATE TABLE readme_projection_transcripts AS
SELECT t.*,
       (SELECT list(struct_pack(exon_start := e.exon_start, exon_end := e.exon_end,
                                exon_cdna_start := e.exon_cdna_start, exon_cdna_end := e.exon_cdna_end,
                                phase := e.phase, end_phase := e.end_phase) ORDER BY e.exon_start)
        FROM readme_exons e WHERE e.transcript_index = t.transcript_index) AS exons,
       []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] AS peptide_edits
FROM readme_transcripts t;

CREATE TABLE demo_annotations AS
SELECT * FROM query(duckvep_annotate_sql('demo_events', 'demo'));

SELECT event_index, reference_amino_acids, alternate_amino_acids, reference_codons, alternate_codons
FROM query(duckvep_transcript_projection_sql(
  'demo_events', 'demo_annotations', 'readme_projection_transcripts'))
ORDER BY event_index;
```

---

## Structural and repeat preparation

These builders normalize raw VCF records into typed preparation rows, one row per input record, with a status and a stable reason. They return SQL text and take a trailing options STRUCT like every builder. Only `duckvep_prepare_structural_hgvs_sql` has an option; the others take no options. `docs/structural-identity-hgvs.md` explains the identity and HGVS rules and the VEP 116 evidence behind them.

<a id="duckvep_prepare_sv_geometry_sql"></a>

### duckvep_prepare_sv_geometry_sql

Builds the SQL that prepares structural geometry and insertion provenance from raw records. Symbolic alleles keep their breakpoint confidence bounds separately from the nominal coordinates. Literal anchored insertions carry their inserted bases; a symbolic `INFO/SEQ` is source provenance, not a literal insertion.

Signatures:

```text
duckvep_prepare_sv_geometry_sql(input_table VARCHAR) -> VARCHAR
duckvep_prepare_sv_geometry_sql(input_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `input_table` names a relation with `event_index`, `pos`, `ref`, `alt` and `info`.

Returns: one row per input event with `event_index`, `status`, `mode`, `nominal_start`, `nominal_end`, `outer_start`, `inner_start`, `inner_end`, `outer_end`, `inserted_sequence` and `source_sequence`. `status` is `ok`, or `unsupported_geometry` for missing, inconsistent or extreme coordinates; such a row does not fail the batch. `mode` is `structural`, `literal_insertion` or `unsupported`.

```sql
CREATE TABLE sv_records AS SELECT * FROM (VALUES
  (0::BIGINT, 13546123::BIGINT, 'G', '<INS>', 'SVTYPE=INS;END=13546123;CIPOS=-2,0;CIEND=0,2;SEQ=ATG'),
  (1, 100, 'A', '<DEL>', 'END=220;CIPOS=-5,5'),
  (2, 500, 'G', 'GATG', '.')
) t(event_index, pos, ref, alt, info);

SELECT * FROM query(duckvep_prepare_sv_geometry_sql('sv_records')) ORDER BY event_index;
```

<a id="duckvep_prepare_breakend_pairs_sql"></a>

### duckvep_prepare_breakend_pairs_sql

Builds the SQL that validates breakend (BND) mate and event identity across raw records. It returns one row per physical VCF record and never merges a record with its mate. ALT syntax alone never proves a fusion, a phase or inserted-only sequence.

Signatures:

```text
duckvep_prepare_breakend_pairs_sql(events_table VARCHAR) -> VARCHAR
duckvep_prepare_breakend_pairs_sql(events_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `events_table` names a relation with `event_index`, `chrom`, `pos`, `id`, `ref`, `alt` and `info`.

Returns: `event_index`, `id`, `chrom`, `pos`, `record_kind`, `status`, `reason`, `mate_id`, `event_id`, `mate_event_index`, `declared_mate_chrom`, `declared_mate_position`, `id_reciprocal`, `coordinate_reciprocal`, `orientation_reciprocal`, `insert_agree`, `event_agree`, `event_record_count`, `inserted_length`, `inserted_sequence`, `pair_key`, `phase_status` and `fusion_status`. `record_kind` is `paired_breakend`, `single_breakend`, `malformed_alt` or `not_breakend`. `status` is `reciprocal`, `unproven`, `conflict`, `invalid` or `not_applicable`, and `reason` is a stable identifier. `fusion_status` is always `not_asserted` and `phase_status` always `not_evaluated`. The evidence columns stay separate so partial agreement is visible.

```sql
CREATE TABLE bnd_records AS SELECT * FROM (VALUES
  (0::INTEGER, '1', 100::BIGINT, 'a1', 'A', 'ATT[2:200[', 'MATEID=a2;EVENT=e1'),
  (1, '2', 200, 'a2', 'C', ']1:100]AAC', 'MATEID=a1;EVENT=e1'),
  (2, '3', 10, 's1', 'G', '.G', 'SVTYPE=BND')
) t(event_index, chrom, pos, id, ref, alt, info);

SELECT event_index, record_kind, status, reason, mate_event_index, pair_key
FROM query(duckvep_prepare_breakend_pairs_sql('bnd_records'))
ORDER BY event_index;
```

<a id="duckvep_prepare_breakend_fusion_sql"></a>

### duckvep_prepare_breakend_fusion_sql

Builds the SQL that joins breakend identity with caller-supplied endpoint genes and reports partner-gene evidence per physical record. It never asserts a fusion, a reading frame or a phase.

Signatures:

```text
duckvep_prepare_breakend_fusion_sql(pairs_table VARCHAR, genes_table VARCHAR) -> VARCHAR
duckvep_prepare_breakend_fusion_sql(pairs_table VARCHAR, genes_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `pairs_table` holds a materialized result of `duckvep_prepare_breakend_pairs_sql`, and `genes_table` has `event_index` and `gene_id`, one row per gene overlapped by that physical endpoint.

Returns: `event_index`, `mate_event_index`, `pair_key`, `status`, `reason`, `endpoint_genes` VARCHAR[], `mate_endpoint_genes` VARCHAR[], `fusion_asserted` and `phase_status`. `status` is `identity_unproven` (with the identity reason), `endpoint_without_gene`, `shared_gene_endpoints`, `candidate_orientation_conflict` or `candidate_partner_genes`. `fusion_asserted` is always false.

```sql
CREATE TABLE bnd_pairs AS SELECT * FROM query(duckvep_prepare_breakend_pairs_sql('bnd_records'));
CREATE TABLE bnd_genes AS SELECT * FROM (VALUES (0, 'G1'), (1, 'G2')) t(event_index, gene_id);

SELECT event_index, mate_event_index, status, endpoint_genes, mate_endpoint_genes, fusion_asserted
FROM query(duckvep_prepare_breakend_fusion_sql('bnd_pairs', 'bnd_genes'))
ORDER BY event_index;
```

<a id="duckvep_prepare_structural_hgvs_sql"></a>

### duckvep_prepare_structural_hgvs_sql

Builds the SQL that describes exact-span symbolic DEL, DUP, DUP:TANDEM and INV alleles as unshifted genomic HGVS, with an equivalent literal edit that the small-variant path can turn into transcript HGVS. Every other allele is `unavailable` or `unsupported` with a stable reason. VEP 116 itself emits no HGVS for symbolic structural alleles or breakends.

Signatures:

```text
duckvep_prepare_structural_hgvs_sql(events_table VARCHAR, reference_table VARCHAR) -> VARCHAR
duckvep_prepare_structural_hgvs_sql(events_table VARCHAR, reference_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `events_table` has `event_index`, `chrom`, `pos`, `ref`, `alt` and `info`, and `reference_table` has `event_index` and `reference_sequence`, the reference from `pos` through `END`.

Options: `max_span` (INTEGER, default 5000, range 1 to 60000) bounds the described span; a longer span is `unsupported` with reason `span_capacity`.

Returns: `event_index`, `hgvs_status` (`supported`, `unsupported` or `unavailable`), `hgvs_reason`, `edit`, `hgvs_g`, `normalization` (`none`: no 3' shifting), `transcript_hgvs_route` (`literal_equivalent` when supported), `literal_position`, `literal_reference` and `literal_alternate`.

```sql
CREATE TABLE shgvs_events AS SELECT * FROM (VALUES
  (0::INTEGER, '21', 7467462::DOUBLE, 'T', '<DEL>', 'END=7467465'),
  (1, '21', 7467462, 'T', 'N[21:7467900[', 'SVTYPE=BND')
) t(event_index, chrom, pos, ref, alt, info);
CREATE TABLE shgvs_reference AS SELECT 0::INTEGER AS event_index, 'TCAG' AS reference_sequence;

SELECT event_index, hgvs_status, hgvs_reason, hgvs_g, literal_reference, literal_alternate
FROM query(duckvep_prepare_structural_hgvs_sql('shgvs_events', 'shgvs_reference', {max_span: 100}))
ORDER BY event_index;
```

<a id="duckvep_prepare_expansionhunter_sql"></a>

### duckvep_prepare_expansionhunter_sql

Builds the SQL that prepares ExpansionHunter v5 repeat alleles. Each input event keeps its source fields, a preparation reason and an exactness decision, so a summary is never promoted to an exact sequence.

Signatures:

```text
duckvep_prepare_expansionhunter_sql(input_table VARCHAR, reference_table VARCHAR) -> VARCHAR
duckvep_prepare_expansionhunter_sql(input_table VARCHAR, reference_table VARCHAR, options STRUCT) -> VARCHAR
```

Parameters: `input_table` has `event_index`, `info`, `format`, `sample`, `ref`, `alt` and the one-based `alt_index`, and `reference_table` has `event_index` and the literal `reference_sequence` from POS+1 through END.

Returns: `event_index`, `status`, `reason`, `source` (a STRUCT echoing the input fields and the reference sequence), `reference_components` and `alternate_components` (each a STRUCT of `unit` and `count`) and `sequence_exact`. `status` values include `ok`, `summary_only`, `incomplete` and `invalid`. Duplicate reference rows never duplicate output rows: the event is `invalid` with reason `ambiguous_reference_sequence`.

```sql
CREATE TABLE eh_records AS SELECT * FROM (VALUES
  (0::INTEGER, 'END=2008;REF=1;RL=3;RU=CAG', 'GT:SO:REPCN:REPCI', '1/2:SPANNING/SPANNING:2/10:2-2/10-10', 'C', '<STR2>,<STR10>', 1::DOUBLE),
  (1, 'END=2008;REF=1;RL=3;RU=CAG', 'GT:SO:REPCN:REPCI', '1/2:SPANNING/SPANNING:2/10:2-2/10-10', 'C', '<STR2>,<STR10>', 2)
) t(event_index, info, format, "sample", ref, alt, alt_index);
CREATE TABLE eh_reference AS SELECT event_index, 'CAG' AS reference_sequence FROM eh_records;

SELECT event_index, status, reason, reference_components, alternate_components, sequence_exact
FROM query(duckvep_prepare_expansionhunter_sql('eh_records', 'eh_reference'))
ORDER BY event_index;
```

<a id="duckvep_repeat_alleles"></a>

### duckvep_repeat_alleles

Expands reference and alternate repeat descriptions, each an ordered list of `(unit, count)`, into literal sequences as one event fact. It is a scalar function evaluated per row. A required `sequence_exact` assertion separates exact descriptions from summaries. The result is not a copy-number inference.

Signatures:

```text
duckvep_repeat_alleles(reference_components, alternate_components, sequence_exact BOOLEAN) -> STRUCT
duckvep_repeat_alleles(reference_components, alternate_components, sequence_exact BOOLEAN, options STRUCT) -> STRUCT
```

Parameters:

| Parameter | Description |
| --- | --- |
| `reference_components`, `alternate_components` | A list of STRUCT(`unit` VARCHAR, `count`), any numeric count type. Units are non-empty IUPAC DNA and their case is preserved. An empty list is an empty allele. A NULL list, or a NULL element or field, means the sequence is unavailable. |
| `sequence_exact` | Must not be NULL. False makes the result `summary_only`. |

Options: `max_allele_bases` (integer from 0 to 2147483647, default 5000). Each complete allele must fit, or the call is an error naming the required length; the limit belongs to this call.

Returns: `STRUCT(reference VARCHAR, alternate VARCHAR, reference_length UBIGINT, alternate_length UBIGINT, length_change BIGINT, length_direction VARCHAR, status VARCHAR)`. `status` is `ok`, `summary_only`, `incomplete_input` or `nonintegral_count`. Unless it is `ok`, every other field is NULL. `length_direction` is `GAIN`, `LOSS` or `NEUTRAL`, by base length.

```sql
SELECT duckvep_repeat_alleles(
  [{unit: 'CAG', count: 10}], [{unit: 'CAG', count: 14}], true) AS exact,
  duckvep_repeat_alleles(
  [{unit: 'CAG', count: 10}], [{unit: 'CAG', count: 14}], false) AS summary;
```

---

## Haplotypes

Phased genotypes become whole-haplotype consequences: one row per occupied shared path, with the edited CDS and protein, the carriers and every contributing event. The v1 contract covers coding sequence, and the supported domain of its predictions is documented in `design/duckvep.md`.

<a id="duckvep_phase_call"></a>

### duckvep_phase_call

Assigns each slot of one genotype to a haplotype lane and phase scope, as explicit typed genotype preparation.

Signatures:

```text
duckvep_phase_call(alleles, phase_before) -> STRUCT[]
duckvep_phase_call(alleles, phase_before, options STRUCT) -> STRUCT[]
```

Parameters: `alleles` is a list of integer allele indices (NULL for a missing call) and `phase_before` a list of booleans of the same length: a slot is true when it is phased to the slot before it. Either argument may be NULL: a NULL `alleles` gives a NULL result, and under the strict policy a NULL `phase_before` leaves every slot unphased. Lists of different lengths, an empty genotype or a ploidy above 65,535 are errors.

Options: `phase_set` (INTEGER, the record's phase-set label; without it phased slots carry no label) and `phase_policy` (VARCHAR, `'strict'` by default, or `'vep116_compat'` for VEP 116's called-slot order).

Returns: one element per genotype slot, `STRUCT(input_slot USMALLINT, allele_index INTEGER, haplotype_lane USMALLINT, ploidy USMALLINT, phase_set BIGINT, phase_scope VARCHAR, status VARCHAR)[]`. `status` is `called`, `unphased` or `missing`; `phase_scope` is `phase_set`, `all_phase_sets`, `allele_slot` or `unresolved`. A slot with no lane has a NULL `haplotype_lane`.

```sql
SELECT a.input_slot, a.allele_index, a.haplotype_lane, a.phase_set, a.phase_scope, a.status
FROM (SELECT unnest(duckvep_phase_call([0, 1], [false, true], {phase_set: 7})) AS a)
ORDER BY a.input_slot;
```

<a id="duckvep_coding_transcripts"></a>

### duckvep_coding_transcripts

The discovery step for haplotypes: one lookup in the resident transcript interval index returns the transcripts whose coding sequence a VCF event touches, so a whole-genome VCF reduces to its coding records before any per-record work. The event is normalized exactly as the annotation builder normalizes it: a shared anchor base creates no overlap, an insertion is an interbase point, and an intron of at most 13 bases inside the CDS counts as coding. The pairs are exactly those to which `duckvep_annotate_sql` assigns the CDS region bit.

Signature:

```text
duckvep_coding_transcripts(model VARCHAR, seq_region, position, reference VARCHAR, alternate VARCHAR) -> UINTEGER[]
```

Parameters: `model` is a loaded model name. `seq_region` and `position` are integers of any width and signedness (one-based position). `reference` and `alternate` are one VCF record and one ALT allele with its shared anchor base.

Returns: the ascending model transcript ordinals, as `UINTEGER[]`. The list is empty for a record in an intron, UTR, flank or intergenic sequence, for alleles that are not literal bases and for events without a difference. Models with wrapped circular objects are refused, as for `duckvep_haplotypes`.

```sql
SELECT event_index, unnest(duckvep_coding_transcripts('demo', seq_region, position, reference, alternate)) AS transcript_index
FROM demo_events
ORDER BY event_index;
```

<a id="duckvep_haplotypes"></a>

### duckvep_haplotypes

Replays phased calls through the loaded model and returns one row per occupied shared path. DuckDB derives complete phase-set domains and sorts the calls; the function consumes the flat calls relation as a query. Incomplete calls and projection or edit failures keep their provenance without inventing sequence. Candidate selection is explicit in the input: `duckvep_coding_transcripts` supplies the `transcript_index` column.

Signature:

```text
duckvep_haplotypes(calls_query VARCHAR, model_name VARCHAR
                   [, phase_policy := VARCHAR] [, input_mode := VARCHAR] [, hgvs := BOOLEAN]
                   [, capacity limits...]) -> TABLE
```

Parameters (positional): `calls_query` is one non-empty SELECT and `model_name` a loaded model without wrapped circular objects.

Input in the default `alt_events` mode is one row per `event_index`, `transcript_index` and `sample_index`, with the columns `seq_region`, `position`, `reference`, `alternate`, `alt_index`, `alleles`, `phase_before` and a nullable `phase_set`. `event_index` identifies an individual ALT event. In `source_records` mode the columns are `event_index` (the whole source record), `seq_region`, `position`, `reference`, `alternates` (the complete source ALT list), `transcript_index`, `sample_index` and `gt` (the original VCF text).

Named parameters:

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `phase_policy` | VARCHAR | `'strict'` | `'strict'` interprets GT and PS strictly; `'vep116_compat'` follows VEP 116 called-slot order. |
| `input_mode` | VARCHAR | `'alt_events'` | `'alt_events'` for decoded per-ALT calls, or `'source_records'` for raw GT text, which requires `phase_policy := 'vep116_compat'`. |
| `hgvs` | BOOLEAN | false | Request bounded protein HGVS for supported completed paths (`hgvsp`, `hgvsp_status`). |

Capacity limits are positive integers. A path that would exceed one is an error, never a truncated result. Defaults: `max_active_events` 16384, `max_active_transcripts` 4096, `max_active_carriers` 65536, `max_active_prefixes` 262144, `max_active_projections` 262144, `max_allele_bytes` 8388608, `max_leaf_events` 4096, `max_leaf_edits` 65536, `max_sequence_bases` 1048576, `max_ploidy` 64 (at most 65535), `max_phase_sets` 1024, `max_alignment_cells` 16777216, `max_leaf_differences` 65536, `max_hgvs_operations` 65536, `max_hgvs_bytes` 1048576, `max_hgvs_reference_bytes` 262144 and `workspace_limit` 268435456 bytes.

Returns:

| Column | Type | Description |
| --- | --- | --- |
| `transcript_index` | UINTEGER | Model transcript ordinal. |
| `cds`, `protein` | VARCHAR | The edited coding sequence and its translation. |
| `sequence_flags`, `evidence_flags` | UINTEGER, UTINYINT | Sequence and evidence bit flags. |
| `projection_status`, `sequence_status` | VARCHAR | Whether the edits projected and whether the sequence was rebuilt (`ok`, or an explicit reason). |
| `edit_count`, `carrier_count` | UBIGINT, UINTEGER | Edits in the path and carriers of it. |
| `carriers` | STRUCT[] | `sample_index`, `phase_set`, `haplotype_lane`, `ploidy`. |
| `contributors`, `contributor_provenance` | STRUCT[] | Every source event of the path with its raw alleles, evidence, projection status and, in the provenance list, `alt_index`, `role` and edit count. |
| `coding_blocks`, `normalized_edits` | STRUCT[] | Physical edits grouped by shared codon or displaced frame, and the differing CDS edit islands with their source event. |
| `cds_differences`, `protein_differences` | STRUCT[] | Aligned differing runs (`ref_start0`, `alt_start0`, `reference`, `alternate`, `alignment_start0`). |
| `stop_in_displaced_frame` | BOOLEAN | Whether the first stop codon overlaps a frame-displaced span. |
| `hgvsp`, `hgvsp_status` | VARCHAR | Protein HGVS and its status (`not_requested` when `hgvs` is false). |
| `prediction_policy`, `prediction_status`, `prediction_reason` | VARCHAR | The versioned contract (`duckvep-coding-v1`) and whether a path is in its supported domain. |
| `carrier_predictions` | STRUCT[] | The per-carrier keyed result, with impact, consequences and NMD. |
| `haplotype_consequences`, `haplotype_impact` | VARCHAR[], VARCHAR | The whole-protein Sequence Ontology set and IMPACT of the edited sequence, NULL unless `prediction_status` is `predicted`. |
| `nmd_rule`, `nmd_prediction`, `nmd_stop_position`, `nmd_junction_position`, `nmd_contributors` | VARCHAR, VARCHAR, UBIGINT, UBIGINT, UBIGINT[] | The whole-haplotype NMD prediction under rule `ejc50-v1` and its evidence. |
| `nominal_length_diff` | BIGINT | Signed sum of projected replacement ALT-minus-REF lengths. |

`nmd_prediction` is an EJC-distance heuristic on the whole haplotype, not a union of single-allele results and not the VEP NMD plugin.

```sql
CREATE TABLE demo_calls AS
SELECT event_index, seq_region, position, reference, alternate,
       1::UINTEGER alt_index, 0::UINTEGER transcript_index, 0::UBIGINT sample_index,
       [1, 0]::INTEGER[] alleles, [false, true]::BOOLEAN[] phase_before, NULL::BIGINT phase_set
FROM demo_events WHERE position IN (124, 125);

SELECT transcript_index, carrier_count, sequence_status, length(contributors) AS contributors,
       protein_differences
FROM duckvep_haplotypes('SELECT * FROM demo_calls', 'demo', hgvs := true);
```

---

## Geometry and helpers

<a id="duckvep_allele_geometry"></a>

### duckvep_allele_geometry

Returns the normalized geometry of one small allele: the uploaded, VEP-feature and minimized-edit coordinates, kept apart. All intervals are zero-based and half-open; an insertion is an explicit interbase site rather than a fake base.

Signature:

```text
duckvep_allele_geometry(position UBIGINT, reference VARCHAR, alternate VARCHAR) -> STRUCT
```

Parameters: a one-based position that fits UINTEGER, and distinct, non-empty REF and ALT alleles of at most 65,535 bases in `A`, `C`, `G`, `T` and `N` (either case). Any NULL argument gives a NULL result. Other input, or identical alleles, is an error.

Returns: `STRUCT(kind_code UTINYINT, interbase BOOLEAN, anchor_side_code UTINYINT, raw_start0 UBIGINT, raw_end0 UBIGINT, feature_start0 UBIGINT, feature_end0 UBIGINT, edit_start0 UBIGINT, edit_end0 UBIGINT, insertion_boundary0 UBIGINT, reference_difference_offset USMALLINT, reference_difference_length USMALLINT, alternate_difference_offset USMALLINT, alternate_difference_length USMALLINT)`. `kind_code` and `anchor_side_code` are numeric codes. `insertion_boundary0` is NULL unless `interbase`. The two difference offsets and lengths locate the differing region inside REF and ALT.

```sql
SELECT duckvep_allele_geometry(100, 'AAC', 'ATC') AS substitution,
       duckvep_allele_geometry(100, 'A', 'ATG') AS insertion;
```

<a id="duckvep_breakend_geometry"></a>

### duckvep_breakend_geometry

Parses a VCF breakend ALT into the mate coordinate, orientation and replacement sequence.

Signature:

```text
duckvep_breakend_geometry(alt VARCHAR) -> STRUCT
```

Parameter: the ALT text. A NULL or non-breakend ALT (a plain allele, `<DEL>`) gives a NULL STRUCT whose children are NULL, so `TRY` and projections are safe. A malformed breakend (unmatched brackets, an empty replacement, a mate position above UBIGINT) is an error; wrap it in `TRY` to turn that into NULL.

Returns: `STRUCT(mate_chrom VARCHAR, mate_position UBIGINT, local_join_after BOOLEAN, mate_extends_right BOOLEAN, replacement_sequence VARCHAR)`. A single breakend (a leading or trailing dot) has no mate, so `mate_chrom`, `mate_position` and `mate_extends_right` are NULL.

```sql
SELECT alt, TRY(duckvep_breakend_geometry(alt)) AS geometry
FROM (VALUES ('A[chr2:321682['), ('.A'), ('<DEL>'), ('N[1:5')) t(alt);
```

---

## Resource control

Every native owner (model arrays, interval indexes, reference readers, workspaces, per-worker result and text arenas, SQL builders) allocates through one process-wide atomic budget, 4 GiB by default. A refusal is an explicit `capacity error: ... budget exceeded` naming the requested bytes, the bytes in use and the limit; nothing is truncated. A model load that exceeds the budget publishes nothing. Annotation admits at most 6 concurrent workers by default, each holding a 128 MiB native scratch lease and a 256 MiB emitted-output allowance, and an idle worker keeps at most 64 MiB. htslib's own buffers are not routed through the budget; a fixed reservation stands in for each open FASTA index.

<a id="duckvep_native_budget"></a>

### duckvep_native_budget

Reports native memory use for every owner and in total.

Signature:

```text
duckvep_native_budget() -> TABLE
```

Returns: one row per owner (`model`, `index`, `reference`, `workspace`, `scratch`, `emit`, `control`) and a `total` row, with `owner` VARCHAR and `current_bytes`, `high_water_bytes`, `limit_bytes`, `charges` and `refusals` (all UBIGINT). `charges` and `refusals` are process-wide counters reported on the `total` row and zero elsewhere.

```sql
SELECT owner, current_bytes, high_water_bytes, limit_bytes
FROM duckvep_native_budget()
ORDER BY owner;
```

<a id="duckvep_native_budget_set"></a>

### duckvep_native_budget_set

Sets the process-wide ceiling in bytes.

Signature:

```text
duckvep_native_budget_set(bytes BIGINT) -> BIGINT
```

Returns the new ceiling. A non-positive or NULL value is an error, and so is a value below the bytes already charged.

```sql
SELECT duckvep_native_budget_set(4294967296) AS limit_bytes;
```

<a id="duckvep_native_budget_reset_high_water"></a>

### duckvep_native_budget_reset_high_water

Restarts the high-water marks reported by `duckvep_native_budget()`.

Signature:

```text
duckvep_native_budget_reset_high_water() -> BOOLEAN
```

Always returns true.

```sql
SELECT duckvep_native_budget_reset_high_water() AS reset;
```

<a id="duckvep_worker_limits_set"></a>

### duckvep_worker_limits_set

Sets the annotation worker limits.

Signature:

```text
duckvep_worker_limits_set(workers BIGINT, scratch_bytes BIGINT, emit_bytes BIGINT, idle_bytes BIGINT) -> BOOLEAN
```

Parameters: `workers` is the maximum number of concurrent annotation workers (1 to 1024). `scratch_bytes` is the native scratch each worker holds, `emit_bytes` its emitted-output allowance and `idle_bytes` how much an idle worker keeps; each is non-negative. An allele or vector over its lease is a capacity error. Anything out of range, or NULL, is an error. Returns true.

```sql
SELECT duckvep_worker_limits_set(6, 134217728, 268435456, 67108864) AS applied;
```

---

## Vocabulary

<a id="duckvep_so_terms"></a>

### duckvep_so_terms

Lists the Sequence Ontology consequence terms DuckVEP can report. `consequence_mask` in the annotation results is a bit set over these terms.

Signature:

```text
duckvep_so_terms() -> TABLE
```

Returns one row per term with `bit_index` UTINYINT, `consequence_mask` UBIGINT (the single bit for the term), `consequence` VARCHAR, `impact_code` UTINYINT, `impact` VARCHAR (`MODIFIER`, `LOW`, `MODERATE` or `HIGH`), `severity_rank` UTINYINT and `evaluator_tier` UTINYINT.

```sql
SELECT bit_index, consequence, impact, severity_rank
FROM duckvep_so_terms()
ORDER BY severity_rank
LIMIT 5;
```
