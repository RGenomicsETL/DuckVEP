-- Run from the repository root on a disposable copy of the GRCh37 model.
LOAD 'build/release/duckvep.duckdb_extension';
CREATE TEMP TABLE rebuilt_regions AS
SELECT * FROM query(duckvep_ensembl_regions_sql('ensembl_core', 'reference_chunks', 'GRCh37'));
CREATE TEMP TABLE rebuilt_transcripts AS
SELECT * FROM query(duckvep_ensembl_transcripts_sql('ensembl_core', 'reference_chunks', 'GRCh37'));
CREATE TEMP TABLE rebuilt_receipt AS SELECT model_sha256, region_count, transcript_count
FROM query(duckvep_model_receipt_sql(
  'rebuilt_regions', 'rebuilt_transcripts',
  'Ensembl GRCh37', '116', 'GRCh37',
  'c379d7bf8a50991ea2e5133a15e6870f545b2a1db55e191e7d5624f47067ee97',
  '0a43b56dec40debae976d6e70cac68ea6ed874f9fb7c8c814363702ff1d47865',
  'VEP 116 core: is_current=1, non-empty stable ID, biotype!=artifact, no readthrough_tra; GRCh37 primary-assembly FASTA regions'));
SELECT model_sha256, region_count, transcript_count FROM rebuilt_receipt;
SELECT CASE WHEN model_sha256 = '21e113d9148132491bc935f3d1b0ec7d50663f450b62e346cb1f0f447de0b290'
    THEN 'GRCh37 model hash matches' ELSE error('GRCh37 model hash mismatch') END
FROM rebuilt_receipt;
