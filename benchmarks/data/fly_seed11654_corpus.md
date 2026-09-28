# Drosophila melanogaster BDGP6.54 model-derived corpus (seed 11654)

```sh
Rscript scripts/stage_species_corpus.R \
  /tmp/duckvep-fly/fly.duckdb /tmp/duckvep-fly/fly.fa \
  /tmp/duckvep-fly/fly-seed11654.vcf 11654 drosophila_melanogaster 8
```

The generator samples up to eight transcripts per contig/biotype/strand group, checks reference alleles against the indexed FASTA, and produces SNVs, MNVs, insertions, and deletions. The VCF contains 3,527 unique variants: 878 of each allele shape and 15 single-exon codon witnesses. It covers all ten model biotypes, both strands, and 29 contigs including `mitochondrion_genome`.

Twelve nuclear table-1 witnesses mutate an in-frame TGG to TGA (`stop_gained`). Three mitochondrial table-5 witnesses mutate TGG to TGA on the appropriate strand (`synonymous_variant`, Trp→Trp instead of the table-1 stop). The mitochondrial witness IDs and designated transcripts are `DM002542`/`FBtr0100857` (positive strand), `DM002607`/`FBtr0433501` and `DM002640`/`FBtr0433499` (negative strand); the [witness pair table](fly_seed11654_differential/codon_witness_pairs.csv) records oracle and DuckVEP terms.

| Output | SHA-256 |
| --- | --- |
| `fly-seed11654.vcf` | `8dcd7cb98341d8be202f022d9ac1865a3e21d905fc91ad97523aae583bce2534` |
| `fly-seed11654.vcf.provenance.tsv` | `3ad98b068f39e1af3efe09f5d965cdddedc2db3ebddcf05c24cf9ec45731228f` |
