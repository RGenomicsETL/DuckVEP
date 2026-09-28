# DuckVEP

## Breaking change: native SQL builders

DuckVEP registers functions without modifying the database catalog during `LOAD`. Preparation, annotation and projection run SQL emitted by native scalar builders in the caller's connection. Existing macro invocations must be migrated; `LOAD` does not install compatibility macros.

| Before | After |
| --- | --- |
| `FROM duckvep_annotate('t', 'model', hgvs := true)` | `FROM query(duckvep_annotate_sql('t', 'model', {hgvs: true}))` |
| `FROM duckvep_ensembl_regions('core', 'ref', 'GRCh38')` | `FROM query(duckvep_ensembl_regions_sql('core', 'ref', 'GRCh38'))` |
| `FROM duckvep_ensembl_transcripts('core', 'ref', 'GRCh38')` | `FROM query(duckvep_ensembl_transcripts_sql('core', 'ref', 'GRCh38'))` |
| `FROM duckvep_ensembl_regulation_features('funcgen', 'regions')` | `FROM query(duckvep_ensembl_regulation_features_sql('funcgen', 'regions'))` |
| `FROM duckvep_transcript_projection('events', 'annotations', 'transcripts')` | `FROM query(duckvep_transcript_projection_sql('events', 'annotations', 'transcripts'))` |
| `FROM duckvep_model_receipt('regions', 'transcripts', ...)` | `FROM query(duckvep_model_receipt_sql('regions', 'transcripts', ...))` |

`duckvep_annotate_sql()` requires the model name. Optional builder settings are fields of a STRUCT rather than named function parameters. `Rduckvep` exposes matching `rduckvep_*_sql()` wrappers and `rduckvep_annotate()` for a data-frame result.
