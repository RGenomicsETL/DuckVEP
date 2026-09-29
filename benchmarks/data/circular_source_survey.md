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
transcript was found in these two release/species combinations. This is not
a survey of all public Ensembl releases or species; a public crossing
transcript suitable for an executable VEP differential remains unidentified.

SHA-256 of the compressed `transcript.txt.gz` downloads, respectively:
`8aa0bf8182e8e6db9d3f6c8be809d39f6e72d3d34570b9ca5fe3657e650f0d59`
and `57a9b22eb4183d893d303ebec3da006c54d5c3df207933e762ec6570fc321b2e`.
SHA-256 of `seq_region_attrib.txt.gz`, respectively:
`9a352ddd9753e9ade957f4bec9e1e7f69033cb` and
`0cf0ee6a200b65eea7031e8c591cabc67cc4c777d8b53cd594af7990e143fe16`.
