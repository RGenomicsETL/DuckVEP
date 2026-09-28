-- Run after expanding consequences against the native GRCh37 model.
-- Replace OUTPUT_DIR with the receipt-matched external Parquet directory.
-- consequence_rows has model_sha256, transcript_index and caller event columns.
WITH mapped AS (
  SELECT model_sha256, transcript_index,
         list(struct_pack(status := mane_status, refseq_nuc := refseq_nuc,
                          refseq_prot := refseq_prot, source_digest := mane_sha256,
                          target_assembly := target_assembly)
              ORDER BY mane_status, refseq_nuc) AS mane_mapped_to_grch37
  FROM read_parquet('OUTPUT_DIR/mane_grch37_mapping.parquet')
  WHERE mapping_status = 'exact_model_match'
  GROUP BY model_sha256, transcript_index
)
SELECT c.*, n.canonical AS ensembl_grch37_canonical,
       n.gencode_basic AS gencode19_basic,
       n.no_retained_canonical, n.model_sha256 AS native_model_sha256,
       mapped.mane_mapped_to_grch37
FROM consequence_rows AS c
JOIN read_parquet('OUTPUT_DIR/grch37_transcript_authorities.parquet') AS n
  ON n.model_sha256 = c.model_sha256 AND n.transcript_index = c.transcript_index
LEFT JOIN mapped
  ON mapped.model_sha256 = c.model_sha256 AND mapped.transcript_index = c.transcript_index;
