WITH profiles AS (
  SELECT unnest([
    {id:'include',blocks:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],rna:'ATGAAACCCAAATTTCCCTAAATGAAATAAAAA',rescue:false},
    {id:'skip',blocks:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],rna:'ATGAAACCCTTTCCCTAAATGAAATAAAAA',rescue:false},
    {id:'retain',blocks:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:10::UBIGINT,end1:12::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],rna:'ATGAAACCCGGGAAATTTCCCTAAATGAAATAAAAA',rescue:false},
    {id:'rescue',blocks:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],rna:'ATGAAACCCAAATTTCCCTAAATGAAATAAAAA',rescue:true}
  ]) profile
), strands AS (SELECT unnest([1::TINYINT,-1::TINYINT]) strand),
scenario AS (
  SELECT profile.id, strand, duckvep_product_scenario({
    reference:{start1:4294967297::UBIGINT,
      dna:CASE WHEN strand=1 THEN 'ATGAAACCCGGGAAAGGGTTTCCCTAAGGGATGAAATAAGGGAAA'
               ELSE _duckvep_revcomp('ATGAAACCCGGGAAAGGGTTTCCCTAAGGGATGAAATAAGGGAAA') END,
      strand:strand},
    exons:CASE WHEN strand=1 THEN
      list_transform([{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],
        lambda e: {start1:e.start1+4294967296::UBIGINT,end1:e.end1+4294967296::UBIGINT})
      ELSE list_transform(list_reverse([{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}]),
        lambda e: {start1:46::UBIGINT-e.end1+4294967296::UBIGINT,end1:46::UBIGINT-e.start1+4294967296::UBIGINT}) END,
    blocks:CASE WHEN strand=1 THEN
      list_transform(profile.blocks,lambda e: {start1:e.start1+4294967296::UBIGINT,end1:e.end1+4294967296::UBIGINT})
      ELSE list_transform(profile.blocks,lambda e: {start1:46::UBIGINT-e.end1+4294967296::UBIGINT,end1:46::UBIGINT-e.start1+4294967296::UBIGINT}) END,
    sources:CASE WHEN strand=1 THEN
      [{source_id:22::UBIGINT,start1:29::UBIGINT+4294967296::UBIGINT,reference:'G',alternate:'C'},
       {source_id:11::UBIGINT,start1:(CASE WHEN profile.rescue THEN 2 ELSE 13 END)::UBIGINT+4294967296::UBIGINT,reference:CASE WHEN profile.rescue THEN 'T' ELSE 'A' END,alternate:'C'}]
      ELSE [{source_id:22::UBIGINT,start1:17::UBIGINT+4294967296::UBIGINT,reference:'C',alternate:'G'},
            {source_id:11::UBIGINT,start1:(CASE WHEN profile.rescue THEN 44 ELSE 33 END)::UBIGINT+4294967296::UBIGINT,reference:CASE WHEN profile.rescue THEN 'A' ELSE 'T' END,alternate:'G'}] END,
    programme:{id:profile.id,evidence_id:'declared',
      assumption_id:CASE WHEN profile.rescue THEN 'alternative-start' ELSE 'splice' END,
      expected_reference_rna:profile.rna,annotated_start1:1::UINTEGER,
      alternative_start1:(CASE WHEN profile.rescue THEN 22 ELSE 0 END)::UINTEGER,
      alternative_trigger_id:CASE WHEN profile.rescue THEN 'annotated-start-lost' ELSE '' END,
      genetic_code:1::UTINYINT}
  }) p FROM profiles CROSS JOIN strands
)
SELECT id, strand, typeof(p) AS result_type, to_json(p) AS full_typed_fields
FROM scenario ORDER BY id,strand;
