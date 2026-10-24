SELECT p.source_effects[1].genomic_start1 AS genomic_start1,
       p.source_effects[1].genomic_end1 AS genomic_end1,
       p.reference_start1 AS reference_start1
FROM (SELECT duckvep_product_scenario({
 reference:{start1:18446744073709551607::UBIGINT,dna:'ATGAAATAA',strand:1::TINYINT},
 exons:[{start1:18446744073709551607::UBIGINT,end1:18446744073709551615::UBIGINT}],
 blocks:[{start1:18446744073709551607::UBIGINT,end1:18446744073709551615::UBIGINT}],
 sources:[{source_id:11::UBIGINT,start1:18446744073709551610::UBIGINT,reference:'A',alternate:'C'}],
 programme:{id:'edge',evidence_id:'declared',assumption_id:'edge',
   expected_reference_rna:'ATGAAATAA',annotated_start1:1::UINTEGER,
   alternative_start1:0::UINTEGER,alternative_trigger_id:'',genetic_code:1::UTINYINT}
}) p);
