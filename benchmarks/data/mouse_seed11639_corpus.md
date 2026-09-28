# Mus musculus GRCm39 model-derived corpus (seed 11639)

```sh
Rscript scripts/stage_species_corpus.R \
  /tmp/duckvep-mouse/mouse.duckdb /tmp/duckvep-mouse/mouse.fa \
  /tmp/duckvep-mouse/mouse-seed11639.vcf 11639 mus_musculus 8
```

The generator samples up to eight transcripts per contig, biotype, and strand from the pinned model. Each sampled transcript contributes reference-checked SNV, MNV, insertion, and deletion records; single-exon coding transcripts also supply codon witnesses. The output has 22,204 variants: 5,547 of each variant shape and 16 codon witnesses (12 table-1 TGG→TGA; four mitochondrial table-2 TGG→TGA). It covers all 43 modeled biotypes and both strands across 38 contigs, including MT. Thirty-three biotypes occur on both strands; the other ten have transcripts on only one strand in the source model.

| Output | SHA-256 |
| --- | --- |
| `mouse-seed11639.vcf` | `bb0c04e37245ac4dbc253f590e0d87c366fbdb75253528a60606883b3444072f` |
| `mouse-seed11639.vcf.provenance.tsv` | `ef314b1f0ed62f4facbabab9eeb870d31c850504743e712a57aa104957dd1bee` |

The same generalized generator with species `plasmodium_falciparum` and its seed 11663 reproduces the pinned P. falciparum VCF and provenance byte-for-byte.
