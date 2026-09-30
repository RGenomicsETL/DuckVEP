CREATE SCHEMA duckvep_core;
CREATE TABLE duckvep_reference_chunks(chrom VARCHAR, "start" BIGINT, "end" BIGINT, seq VARCHAR);
INSERT INTO duckvep_reference_chunks VALUES
  ('1', 0, 20, 'ATGAAACCCCGGGTAACCCC'),
  ('1', 20, 40, 'CCTTTAGGGCATAAAAAAAA');
CREATE TABLE duckvep_core.coord_system(
  coord_system_id BIGINT, species_id BIGINT, name VARCHAR, version VARCHAR, rank BIGINT
);
INSERT INTO duckvep_core.coord_system VALUES
  (1, 1, 'chromosome', 'GRCh38', 1),
  (2, 2, 'chromosome', 'GRCh38', 1);
CREATE TABLE duckvep_core.seq_region(
  seq_region_id BIGINT, name VARCHAR, coord_system_id BIGINT, length BIGINT
);
INSERT INTO duckvep_core.seq_region VALUES
  (10, '1', 1, 40),
  (11, '1', 2, 40);
CREATE TABLE duckvep_core.seq_region_attrib(
  seq_region_id BIGINT, attrib_type_id BIGINT, value BIGINT
);
CREATE TABLE duckvep_core.gene(
  gene_id BIGINT, biotype VARCHAR, is_current BIGINT, stable_id VARCHAR, version BIGINT
);
INSERT INTO duckvep_core.gene VALUES
  (20, 'protein_coding', 1, 'ENSG_TEST_1', 1),
  (21, 'protein_coding', 1, 'ENSG_TEST_2', 1);
CREATE TABLE duckvep_core.transcript(
  transcript_id BIGINT, gene_id BIGINT, seq_region_id BIGINT,
  seq_region_start BIGINT, seq_region_end BIGINT, seq_region_strand BIGINT,
  biotype VARCHAR, is_current BIGINT, stable_id VARCHAR, version BIGINT
);
INSERT INTO duckvep_core.transcript VALUES
  (30, 20, 10, 1, 18, 1, 'protein_coding', 1, 'ENST_TEST_FORWARD', 1),
  (31, 21, 10, 21, 32, -1, 'nonsense_mediated_decay', 1, 'ENST_TEST_REVERSE', 1),
  (32, 21, 10, 21, 32, -1, 'protein_coding', 1, 'ENST_TEST_EDITED', 1);
CREATE TABLE duckvep_core.exon(
  exon_id BIGINT, seq_region_id BIGINT, seq_region_start BIGINT,
  seq_region_end BIGINT, seq_region_strand BIGINT, phase BIGINT,
  end_phase BIGINT, is_current BIGINT, stable_id VARCHAR
);
INSERT INTO duckvep_core.exon VALUES
  (40, 10, 1, 6, 1, -1, 0, 1, 'ENSE_TEST_1'),
  (41, 10, 11, 18, 1, 0, 0, 1, 'ENSE_TEST_2'),
  (42, 10, 21, 32, -1, 0, 0, 1, 'ENSE_TEST_3');
CREATE TABLE duckvep_core.exon_transcript(exon_id BIGINT, transcript_id BIGINT, rank BIGINT);
INSERT INTO duckvep_core.exon_transcript VALUES
  (40, 30, 1), (41, 30, 2), (42, 31, 1), (42, 32, 1);
CREATE TABLE duckvep_core.translation(
  translation_id BIGINT, transcript_id BIGINT, seq_start BIGINT,
  start_exon_id BIGINT, seq_end BIGINT, end_exon_id BIGINT,
  stable_id VARCHAR, version BIGINT
);
INSERT INTO duckvep_core.translation VALUES
  (50, 30, 1, 40, 6, 41, 'ENSP_TEST_FORWARD', 1),
  (51, 31, 1, 42, 9, 42, 'ENSP_TEST_REVERSE', 1),
  (52, 32, 1, 42, 9, 42, 'ENSP_TEST_EDITED', 1);
CREATE TABLE duckvep_core.attrib_type(attrib_type_id BIGINT, code VARCHAR);
INSERT INTO duckvep_core.attrib_type VALUES
  (60, 'MANE_Select'), (61, 'gencode_basic'), (62, 'ccds_transcript'),
  (63, 'gencode_primary'), (64, 'initial_met'), (65, '_rna_edit'),
  (66, 'codon_table'), (68, 'amino_acid_sub'), (69, 'circular_seq');
INSERT INTO duckvep_core.seq_region_attrib VALUES (10, 66, 2);
CREATE TABLE duckvep_core.transcript_attrib(
  transcript_id BIGINT, attrib_type_id BIGINT, value VARCHAR
);
INSERT INTO duckvep_core.transcript_attrib VALUES
  (30, 60, 'NM_TEST'), (30, 61, '1'), (30, 62, 'CCDS_TEST'), (31, 63, '1'),
  (32, 65, '4 5 A');
CREATE TABLE duckvep_core.translation_attrib(
  translation_id BIGINT, attrib_type_id BIGINT, value VARCHAR
);
INSERT INTO duckvep_core.translation_attrib VALUES
  (50, 68, '2 2 U'),
  (51, 65, '1 2'),
  (52, 64, '1 1 M');
CREATE SCHEMA duckvep_funcgen;
CREATE TABLE duckvep_funcgen.feature_type(
  feature_type_id BIGINT, name VARCHAR, so_accession VARCHAR, so_term VARCHAR
);
INSERT INTO duckvep_funcgen.feature_type VALUES
  (70, 'Promoter', 'SO:0000167', 'promoter'),
  (71, 'EMAR', 'SO:0001720', 'epigenetically_modified_region');
CREATE TABLE duckvep_funcgen.regulatory_feature(
  regulatory_feature_id BIGINT, feature_type_id BIGINT, seq_region_id BIGINT,
  seq_region_strand BIGINT, seq_region_start BIGINT, seq_region_end BIGINT,
  stable_id VARCHAR, regulatory_build_id BIGINT
);
INSERT INTO duckvep_funcgen.regulatory_feature VALUES
  (80, 70, 10, 0, 2, 6, 'ENSR_TEST_1', 1),
  (81, 71, 10, 0, 20, 22, 'ENSR_TEST_EMAR', 1);
CREATE TABLE duckvep_funcgen.motif_feature(
  motif_feature_id BIGINT, binding_matrix_id BIGINT, seq_region_id BIGINT,
  seq_region_start BIGINT, seq_region_end BIGINT, seq_region_strand BIGINT,
  score DOUBLE, stable_id VARCHAR
);
INSERT INTO duckvep_funcgen.motif_feature VALUES
  (90, 100, 10, 11, 14, 1, 8.5, 'ENSM_TEST_1');
