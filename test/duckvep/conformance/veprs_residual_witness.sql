-- Run from the repository root after loading DuckVEP and setting the SQL variable
-- comparison_dir to a directory produced by benchmark_duckvep_vep_rs.sh.
-- The snapshot and mappings are read-only; all witness tables are temporary.

SET threads=1;
SELECT duckvep_model_restore('residual_witness',
  getvariable('comparison_dir') || '/model.dvsnap');

CREATE TEMP TABLE residual_events AS
SELECT row_number() OVER (ORDER BY r.seq_region, v.pos, v.id)::UBIGINT AS event_index,
       v.id AS vcf_id,
       r.seq_region::UINTEGER AS seq_region,
       v.pos::UBIGINT AS "position",
       v.ref AS reference,
       v.alt AS alternate
FROM read_csv('test/duckvep/conformance/veprs_residual_witness.vcf',
       delim='\t', header=false, skip=2, quote='', escape='', auto_detect=false,
       columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR',
                'alt':'VARCHAR','qual':'VARCHAR','filter':'VARCHAR','info':'VARCHAR'}) v
JOIN read_parquet(getvariable('comparison_dir') || '/regions.parquet') r
  ON r.seq_region_name = v.chrom;

CREATE TEMP TABLE residual_alleles AS
SELECT event_index, seq_region, "position", reference, alternate,
       NULL::UBIGINT AS end_position,
       NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change,
       NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM residual_events;

SELECT e.vcf_id, t.transcript_stable_id, a.consequence, a.region,
       a.protein_position, a.reference_amino_acid, a.alternate_amino_acid,
       a.duckvep_status, a.duckvep_reason
FROM query(duckvep_annotate_sql(
       'residual_alleles', 'residual_witness', struct_pack(rich := true))) a
JOIN residual_events e USING (event_index)
JOIN read_parquet(getvariable('comparison_dir') || '/transcript_ids.parquet') t
  USING (transcript_index)
WHERE (e.vcf_id, t.transcript_stable_id) IN (
  ('e284879', 'ENST00000696609'),
  ('e391762', 'ENST00001011173'))
ORDER BY e.vcf_id, t.transcript_stable_id;
