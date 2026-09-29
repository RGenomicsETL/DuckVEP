# Arabidopsis thaliana TAIR10 model-derived corpus (seed 11663)

```sh
Rscript scripts/stage_species_corpus.R \
  /tmp/duckvep-plant/plant.duckdb /tmp/duckvep-plant/plant.fa \
  /tmp/duckvep-plant/plant-seed11663.vcf 11663 arabidopsis_thaliana 8
```

The generator samples up to eight transcripts per contig/biotype/strand group, checks reference alleles against the indexed FASTA and produces 2,667 unique variants: 657 each of SNV, MNV, insertion and deletion, plus 39 single-exon codon witnesses. The corpus covers all nine model biotypes, both strands and all seven contigs, including `Mt` (148 variants) and `Pt` (190 variants). Every modeled region has a FASTA sequence; seed and sampling limit are fixed above.

Nine table-1 nuclear witnesses mutate in-frame TGG→TGA (`stop_gained`). Table 11 has the same ordinary amino-acid assignments as table 1 but a different set of permitted initiation codons: GTG can initiate under table 11 but not under table 1. Of 30 table-11 plastid witnesses, four mutate an *observed GTG start* to GCG; 26 change an ATG start to GTG. All 30 yield `start_lost` in both engines. The [start witness table](plant_seed11663_differential/plastid_start_witnesses.csv) retains the original and alternate codons, transcripts, alleles and both consequences; four GTG→GCG rows have exact `start_lost` matches. This tests initiation at the table-11-specific GTG codon rather than claiming that the amino-acid translation of table 11 differs from table 1.

| Output | SHA-256 |
| --- | --- |
| `plant-seed11663.vcf` | `ccc8614b3e218fc690a7d12d6c9c4838a0f219bc06a27d03d4babb4b5cc11bca` |
| `plant-seed11663.vcf.provenance.tsv` | `a1d86e8af0b1e4bec723045e84faf44363a7c5a38ddc6a9fd9c1ba347ff8d85a` |
