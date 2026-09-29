# Offline species excerpts

`scripts/extract_species_fixture.R PF_MODEL TETRA_MODEL EXTENSION_COPY PF_FASTA TETRA_FASTA` extracts source-derived records here and mirrors the Parquet and indexed FASTA excerpts to `r/Rduckvep/inst/extdata/species/` for installed-package tests. It requires staged models from the [Protists 63 Tetrahymena sources](../../../../benchmarks/data/tetrahymena_jcvi_tta1_sources.md) and [P. falciparum differential](../../../../benchmarks/data/plasmodium_seed11663_differential.md), plus the committed Ensembl 116 GRCh37 core excerpt in `test/data/duckvep/ensembl_core/grch37/`. The GRCh37 native builder loads the provided immutable extension copy. The extraction runs offline given those staged models.

| Product | Source transcript(s) | Regions | Code | Witness |
| --- | --- | ---: | ---: | --- |
| *P. falciparum* GCA000002765v3 | PF3D7_MIT01400.1 | 1 | 4 | mitochondrial TGG → TGA |
| *T. thermophila* JCVI-TTA1-2.2 | EAR80522; EAR80553 | 2 | 6 | TAA → CAA; TAG → CAG |
| Human GRCh37 | ENST00000400678 | 1 | 1 | model coexistence |

The extractor copies the model's source-derived region and transcript columns, including CDS, flanks, and exons, and uses `samtools faidx` to take complete contigs from the pinned FASTAs. The small indexed Tetrahymena FASTA is pinned by name in an additional model during the coexistence tests. It reindexes only the model-local `seq_region` and `transcript_index` ordinals to start at zero and retains `source_seq_region_id` and `transcript_stable_id` for traceability. These excerpts are test inputs, not whole-genome differential receipts. SQL and installed-package R tests load three models at once and verify nonstandard codons. The two mirrored fixture sets together occupy less than 200 KB.
