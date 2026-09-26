CREATE TABLE readme_regions AS SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER)) t(seq_region);
CREATE TABLE readme_transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region,
  100::UBIGINT transcript_start, 250::UBIGINT transcript_end, 1::TINYINT strand,
  0::UINTEGER gene_index, 3::UBIGINT transcript_flags,
  120::UBIGINT cds_start, 240::UBIGINT cds_end,
  'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence,
  1::UTINYINT codon_table,
  'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence,
  'ACGTACGTAC'::BLOB post_cds_sequence;
CREATE TABLE readme_exons AS SELECT * FROM (VALUES
  (0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT),
  (0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)
) t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase);
CREATE TABLE readme_events AS SELECT row_number() OVER ()::UBIGINT AS event_index,
  1::UINTEGER AS seq_region, POS::UBIGINT AS position,
  REF AS reference, ALT[1] AS alternate,
  NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
  NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
  NULL::UBIGINT AS mate_position
FROM read_parquet('test/data/duckvep/minimal_bcsq.parquet')
WHERE POS = 124;
SELECT loaded FROM duckvep_model_load('readme',
  'SELECT * FROM readme_regions ORDER BY seq_region',
  'SELECT * FROM readme_transcripts ORDER BY seq_region, transcript_start',
  'SELECT * FROM readme_exons ORDER BY transcript_index, exon_start');
