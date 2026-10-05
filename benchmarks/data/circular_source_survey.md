# Ensembl 116 circular-region transcript survey

Source: public Ensembl release-116 MySQL dumps at
`https://ftp.ensembl.org/pub/release-116/mysql/`. The checked databases are
`drosophila_melanogaster_core_63_116_11` and
`saccharomyces_cerevisiae_core_63_116_4`. For each database, the survey read
`attrib_type.txt.gz`, `seq_region_attrib.txt.gz`, `seq_region.txt.gz`,
`transcript.txt.gz`, `exon.txt.gz`, and `exon_transcript.txt.gz`. The downloaded
dumps are not model inputs or checked-in fixtures.

| Database | Circular region (`seq_region_id`, length) | Transcripts on region | Inverted transcript bounds | Exon-rank origin crossings |
| --- | --- | ---: | ---: | ---: |
| D. melanogaster | `mitochondrion_genome` (1070, 19524) | 38 | 0 | 0 |
| S. cerevisiae | `Mito` (34, 85779) | 55 | 0 | 0 |

Both databases declare `attrib_type.code = 'circular_seq'` with
`attrib_type_id = 316` and exactly one `seq_region_attrib.value = '1'` for
the region above. The survey joined the attribute to its region, counted
transcripts whose `seq_region_start > seq_region_end`, then joined exons by
`exon_transcript.rank` and counted adjacent rank pairs whose genomic starts
reverse the expected direction for their transcript strand. No wrapped
transcript was found in these two release/species combinations. The survey does
not cover all public Ensembl releases or species; the wider survey below
reports the broader database set and its crossing transcripts.

SHA-256 of the compressed `transcript.txt.gz` downloads, respectively:
`8aa0bf8182e8e6db9d3f6c8be809d39f6e72d3d34570b9ca5fe3657e650f0d59`
and `57a9b22eb4183d893d303ebec3da006c54d5c3df207933e762ec6570fc321b2e`.
SHA-256 of `seq_region_attrib.txt.gz`, respectively:
`9a352ddd9753e9ade957f4bec9e1e7f69033cb` and
`0cf0ee6a200b65eea7031e8c591cabc67cc4c777d8b53cd594af7990e143fe16`.

## Wider Ensembl survey

`scripts/survey_circular_transcripts.py` reads `attrib_type`, `seq_region_attrib` and, only for
databases that declare a `circular_seq = 1` region, `transcript` from the public MySQL dumps,
and counts transcripts on those regions with `seq_region_start > seq_region_end`. The full
tables are `circular_survey_release116.txt` and `circular_survey_genomes63.txt` in this directory.

| Source | Core databases | With a circular region | Inverted transcripts |
| --- | ---: | ---: | ---: |
| Ensembl release 116 (`ftp.ensembl.org/pub/release-116/mysql/`) | 359 | 2 (D. melanogaster, S. cerevisiae) | 0 |
| Ensembl Genomes 63 plants | 267 | 12 | 1 (`daucus_carota_core_63_116_1`) |
| Ensembl Genomes 63 fungi | 81 | 19 | 0 |
| Ensembl Genomes 63 protists | 56 | 4 | 0 |
| Ensembl Genomes 63 metazoa | 386 | 112 | 0 |
| Ensembl Genomes 63 bacteria (129 collections) | 129 | 129 | 18 (`bacteria_0_collection_core_63_116_1`) |

Every Ensembl Genomes bacteria collection declares circular chromosomes and plasmids, but only
`bacteria_0_collection` holds transcripts with inverted bounds. The 19 inverted transcripts are real, public, origin-crossing
objects: 18 in ten bacterial and archaeal genomes of `bacteria_0_collection` (for example
`AAC68473` of *Chlamydia trachomatis* D/UW-3/CX, `AAR38856` of *Nanoarchaeum equitans* Kin4-M, three genes of
*Pyrococcus horikoshii* OT3, five plasmid and chromosome genes of *Synechocystis* sp. PCC 6803) and one in
the *Daucus carota* plastid, `KZM81246`, a trans-spliced rps12 whose exons are ranked 1 to 3 across the
inverted-repeat boundary so that the circular builder classifies it as origin-crossing. Every one has
`seq_region_start > seq_region_end` and passes the exon-rank validation described above.

## Executable VEP differential

A VEP 116 image (`scripts/run_species_vep116_docker.sh`, pinned by digest) and an Ensembl Genomes 63 VEP cache
exist for three of these genomes, so `scripts/circular_vep_differential.py` annotates the same events with both
tools, offline: carrot plastid `Pt` (155,848 bp; `KZM81246`, minus strand), *C. trachomatis* chromosome
(1,042,519 bp; `AAC68473`, plus strand) and *N. equitans* chromosome (490,885 bp; `AAR38856`, minus strand). Events
are every base within 60 bases of the origin and 120 bases of each exon edge of the crossing transcript, sparse
tiling elsewhere, and small deletions and insertions (6,268 to 8,100 events per genome, 5 kb flanks). The
receipts, with source, reference, cache and extension hashes, are in `circular_vep_differential/*.json`.
The HGVS run leaves out events that touch the crossing transcript's exons, because VEP's transcript mapper aborts
on that transcript (`Bio::EnsEMBL::Mapper::map_insert`) for an insertion in its translation.

| Genome | Transcripts | Other-transcript pairs | Identical SO terms | DuckVEP only | VEP only | Different | HGVSc identical | HGVSp identical |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| *D. carota* Pt | 152 | 83,892 | 81,365 | 2,527 | 0 | 0 | 4,714 / 4,738 | 3,192 / 3,197 |
| *C. trachomatis* | 937 | 59,323 | 56,381 | 2,942 | 0 | 0 | 4,862 / 4,862 | 4,792 / 4,793 |
| *N. equitans* | 553 | 80,060 | 78,435 | 1,624 | 0 | 1 | 5,918 / 5,918 | 5,884 / 5,888 |

What the comparison establishes, and what it does not:

- **Unwrapped transcripts agree with VEP wherever topology is not involved.** Every DuckVEP-only pair (2,527,
  2,942 and 1,624) is an `upstream_gene_variant` or `downstream_gene_variant` that exists only through the
  origin: its linear gap exceeds 5 kb and its circular gap does not. VEP, which has no circular semantics,
  cannot report them. There is no VEP-only pair and no different term set caused by the lift. HGVSc differs in
  24 carrot pairs, all within 1,100 bases of the origin, where DuckVEP's 3' shift window continues around the
  origin and VEP's is clipped at the sequence end.
- **Residual disagreements are not caused by lifting.** The one different SO pair (*N. equitans* insertion at the
  first base of `AAR39122`, VEP `inframe_insertion&stop_retained_variant`, DuckVEP
  `coding_sequence_variant&inframe_insertion`) and the six interior HGVSp differences (one to four per genome)
  also occur with the extension at revision `5b6d195`, which has no circular lifting. Each receipt's
  `lifted_vs_prelift_linear` block compares the lifted result to that extension on the unwrapped transcripts:
  no interior difference in any genome, 24 near-origin differences in carrot (the HGVS windows above) and none
  elsewhere, plus the DuckVEP-only flank rows.
- **The crossing transcript itself is not oracle-proved.** VEP does not implement circular topology for a
  transcript stored with `start > end`. For `KZM81246`, `AAC68473` and `AAR38856` it reports no row for most events
  inside the transcript, or `intergenic_variant` as a transcript consequence, and it labels part of the carrot
  transcript's first exon `3_prime_UTR_variant`; DuckVEP reports the exon, intron, splice and coding
  consequences of the physical circular transcript (the last column of each receipt lists the pairs). This
  is a difference in what is modelled, not a matched or refuted result, so the differential neither confirms nor
  refutes DuckVEP on the origin-crossing object. That behavior is proved by rotation equivariance, by
  agreement with an ordinary linear model in which each object is placed mid-sequence, and by the properties in
  `test/duckvep/property/duckvep_prop_circular.c`, not by VEP.
