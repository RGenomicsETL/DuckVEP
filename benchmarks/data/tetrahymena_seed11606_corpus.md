# Tetrahymena thermophila seed-11606 corpus

`Rscript scripts/stage_species_corpus.R MODEL_DB FASTA OUTPUT.vcf 11606 tetrahymena_thermophila 8` samples 26,176 distinct FASTA-checked variants (SNV, MNV, insertions and deletions) across eight biotypes and adds 12 table-6 witnesses. The VCF SHA-256 is `5a80f6f5f1a7064365c7d053a9a661075e3e040fe86c45b6c4a72b79b8382938`, and the adjacent provenance TSV SHA-256 is `ae387d821cc2ac05ed375018dfd4627940263bb43351b5a5c7e890d67f3b7eff`.

The witnesses mutate source TAA or TAG codons to CAA or CAG. Eleven witnesses carry TAA and one carries TAG in the original coding sequence, including `EAR80522` (GG663223:143 T>C, TAA) and `EAR80553` (GG663148:295 T>C, TAG). In the table-6 model these are synonymous glutamine substitutions; under table 1 the same changes lose stops. The [offline fixtures](../../test/data/duckvep/fixtures.md) and coexistence tests exercise both readings without downloading the cache.
