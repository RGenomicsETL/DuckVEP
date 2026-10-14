WITH fixture AS (
 SELECT {
  reference:{start1:1::UBIGINT,dna:'ATGAAACCCGGGAAAGGGTTTCCCTAAGGGATGAAATAAGGGAAA',strand:1::TINYINT},
  exons:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],
  blocks:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],
  sources:[{source_id:0::UBIGINT,start1:13::UBIGINT,reference:'A',alternate:'C'},
           {source_id:22::UBIGINT,start1:29::UBIGINT,reference:'G',alternate:'C'}],
  programme:{id:'include',evidence_id:'declared',assumption_id:'splice',
    expected_reference_rna:'ATGAAACCCAAATTTCCCTAAATGAAATAAAAA',
    annotated_start1:1::UINTEGER,alternative_start1:0::UINTEGER,
    alternative_trigger_id:'',genetic_code:1::UTINYINT}
 } s
), rows AS (
 SELECT i AS row_id, CASE WHEN i=2 THEN struct_update(s,sources := s.sources[:1]) ELSE s END s
 FROM fixture, range(1,4) AS indices(i)
), products AS (
 SELECT row_id,duckvep_product_scenario_v1(s) p FROM rows
)
SELECT row_id,list_count(p.sources),list_count(p.source_effects),p.sources[1].source_id,
       p.source_effects[1].source_id,p.sources[1].reference,p.sources[1].alternate,
       p.source_effects[2].source_id FROM products ORDER BY row_id;
