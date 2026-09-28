# Release-63 P. falciparum transcript-pair differential

Run `scripts/run_plasmodium_vep116_docker.sh /tmp/duckvep-species/cache /tmp/duckvep-species/plasmodium.fa /tmp/duckvep-species/species-seed11663.vcf /tmp/duckvep-species/oracle-seed11663.json`, then `Rscript scripts/compare_plasmodium_corpus.R /tmp/duckvep-species/plasmodium.duckdb /tmp/duckvep-species/species-seed11663.vcf /tmp/duckvep-species/species-seed11663.vcf.provenance.tsv /tmp/duckvep-species/oracle-seed11663.json /tmp/duckvep-species/differential-seed11663`.

Oracle: `ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11`, VEP 116.0, indexed Protists release-63 cache `plasmodium_falciparum/63_GCA000002765v3` (archive SHA-256 `2ed9cdafff5e96a4c1430fdd3993b38504ef9ed35e314d213fd009ab6fea36ab`). Oracle JSON: 5,366 records, SHA-256 `69e0bb4a7a8eca7da9f871fe7561a4acaeea80d265ffc6ebc54f8ff560634c3a`. DuckVEP: rebuilt `d9c7056823755f83d5b189b9e59fdd1129ebe3e342095453e0d2208f3081d81b` model, release extension. The comparison unions SO terms independently per `(VCF ID, transcript stable ID)` and full-outer-joins the pair sets; it checks the 43 designated witness transcripts separately. Distance: 5,000 bp. No regulatory or HGVS comparison.

| Dimension | Exact pairs | Oracle-only | DuckVEP-only | Different SO sets |
| --- | ---: | ---: | ---: | ---: |
| All | 32,131 | 0 | 0 | 0 |
| Code table 1 | 16,299 | 0 | 0 | 0 |
| Code table 4 | 513 | 0 | 0 | 0 |
| Code table 11 | 3,122 | 0 | 0 | 0 |
| Code table not assigned | 12,197 | 0 | 0 | 0 |
| Nuclear | 19,422 | 0 | 0 | 0 |
| Mitochondrial | 6,382 | 0 | 0 | 0 |
| Apicoplast | 6,327 | 0 | 0 | 0 |
| ncRNA | 531 | 0 | 0 | 0 |
| protein_coding | 19,934 | 0 | 0 | 0 |
| pseudogene | 1,200 | 0 | 0 | 0 |
| rRNA | 6,519 | 0 | 0 | 0 |
| snoRNA | 524 | 0 | 0 | 0 |
| snRNA | 40 | 0 | 0 | 0 |
| tRNA | 3,383 | 0 | 0 | 0 |
| Plus strand | 17,620 | 0 | 0 | 0 |
| Minus strand | 14,511 | 0 | 0 | 0 |
| SNV | 8,647 | 0 | 0 | 0 |
| MNV | 7,829 | 0 | 0 | 0 |
| Insertion | 7,829 | 0 | 0 | 0 |
| Deletion | 7,826 | 0 | 0 | 0 |

All 21 observed SO terms agree pairwise, including `stop_gained` (112), `synonymous_variant` (127), `start_lost` (41), `missense_variant` (957), and `frameshift_variant` (1,163). Each count is a number of pairs carrying the term; pairs can have multiple terms. The per-term and disagreement strata are exported to `so_term.csv`, `codon_table.csv`, `contig_class.csv`, `biotype.csv`, `strand.csv`, `allele_shape.csv` and `disagreements.csv` in the untracked output directory.

Designated TGG→TGA witnesses have `stop_gained` in 11 table-1 transcripts and `synonymous_variant` in all 3 table-4 transcripts: table 1 changes W→stop, while table 4 retains W. The 29 table-11 ATG→GTG initiator witnesses are classified `start_lost` by both annotators; this tests their agreement on start-site consequences, not alternative-start retention. No claim of table-11-specific amino-acid behavior follows from that classification.
