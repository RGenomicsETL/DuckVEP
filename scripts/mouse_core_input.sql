SELECT CASE WHEN count(*) = 1 THEN 1
  ELSE error('mouse core miRNA quarantine does not match the pinned snapshot') END
FROM core_snapshot.source_core.transcript_attrib
WHERE transcript_id = 144081 AND attrib_type_id = 15 AND value = '35-56';
CREATE SCHEMA core_input;
CREATE VIEW core_input.attrib_type AS FROM core_snapshot.source_core.attrib_type;
CREATE VIEW core_input.coord_system AS FROM core_snapshot.source_core.coord_system;
CREATE VIEW core_input.seq_region AS FROM core_snapshot.source_core.seq_region;
CREATE VIEW core_input.seq_region_attrib AS FROM core_snapshot.source_core.seq_region_attrib;
CREATE VIEW core_input.gene AS FROM core_snapshot.source_core.gene;
CREATE VIEW core_input.transcript AS FROM core_snapshot.source_core.transcript;
CREATE VIEW core_input.translation AS FROM core_snapshot.source_core.translation;
CREATE VIEW core_input.translation_attrib AS FROM core_snapshot.source_core.translation_attrib;
CREATE VIEW core_input.exon AS FROM core_snapshot.source_core.exon;
CREATE VIEW core_input.exon_transcript AS FROM core_snapshot.source_core.exon_transcript;
CREATE VIEW core_input.transcript_attrib AS
  SELECT * FROM core_snapshot.source_core.transcript_attrib
  WHERE NOT (transcript_id = 144081 AND attrib_type_id = 15 AND value = '35-56');
