WITH subject AS (
  SELECT repeat('ATG',400)||'TAA' AS dna,
         list_transform(range(0,257),lambda i: {
           source_id:(i+1)::UBIGINT,start1:(4+i*3)::UBIGINT,
           reference:'A',alternate:'C'}) AS sources
)
SELECT p.verified_sources, list_count(p.sources) AS original_sources,
       list_count(p.source_effects) AS effects,
       p.source_effects[257].source_id AS last_id,
       p.source_effects[257].changed_rna_bases AS last_rna_bases
FROM (SELECT duckvep_product_scenario_v1({
  reference:{start1:1::UBIGINT,dna:subject.dna,strand:1::TINYINT},
  exons:[{start1:1::UBIGINT,end1:1203::UBIGINT}],
  blocks:[{start1:1::UBIGINT,end1:1203::UBIGINT}],
  sources:subject.sources,
  programme:{id:'257-disjoint',evidence_id:'declared',assumption_id:'fixed-length',
    expected_reference_rna:subject.dna,annotated_start1:1::UINTEGER,
    alternative_start1:0::UINTEGER,alternative_trigger_id:'',genetic_code:1::UTINYINT}
}) p FROM subject);
