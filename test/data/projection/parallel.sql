CREATE TABLE projection_events AS
SELECT i::UBIGINT AS event_index, 1::UINTEGER AS seq_region,
       13::UBIGINT AS position, 'A'::VARCHAR AS reference, 'C'::VARCHAR AS alternate
FROM range(34892) r(i);
CREATE TABLE projection_transcripts AS
SELECT i::UINTEGER AS transcript_index, 1::UINTEGER AS seq_region,
       10::UBIGINT AS transcript_start, 30::UBIGINT AS transcript_end,
       1::TINYINT AS strand, 3::UBIGINT AS transcript_flags,
       10::UBIGINT AS cds_start, 30::UBIGINT AS cds_end,
       'ATGAAACCCGGGTAA'::BLOB AS cds_sequence, 1::UTINYINT AS codon_table,
       ''::BLOB AS post_cds_sequence,
       [{exon_start:10::UBIGINT,exon_end:18::UBIGINT,exon_cdna_start:1::UBIGINT,
         exon_cdna_end:9::UBIGINT,phase:0::TINYINT,end_phase:0::TINYINT},
        {exon_start:25::UBIGINT,exon_end:30::UBIGINT,exon_cdna_start:10::UBIGINT,
         exon_cdna_end:15::UBIGINT,phase:0::TINYINT,end_phase:0::TINYINT}] AS exons,
       []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] AS peptide_edits
FROM range(173111) r(i);
CREATE TABLE projection_annotations AS
SELECT (i % 34892)::UBIGINT AS event_index, (i % 173111)::UINTEGER AS transcript_index,
       1::UBIGINT AS consequence_mask
FROM range(692444) r(i);
