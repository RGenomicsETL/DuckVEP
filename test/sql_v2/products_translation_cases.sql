WITH cases AS (
  SELECT unnest([
    {id:'post_stop',dna:'ATGAAATAATTT',pos:10,ref:'T',alt:'C',init:0,refp:'MK*',altp:'MK*'},
    {id:'stop_lost',dna:'ATGAAATAATTTTGA',pos:7,ref:'T',alt:'C',init:0,refp:'MK*',altp:'MKQF*'},
    {id:'no_stop',dna:'ATGAAATAATTT',pos:7,ref:'T',alt:'C',init:0,refp:'MK*',altp:'MKQF'},
    {id:'annotated_ctg',dna:'CTGAAATAA',pos:4,ref:'A',alt:'G',init:0,refp:'MK*',altp:'ME*'},
    {id:'alternative_ctg',dna:'ATGAAACTGAAATAA',pos:2,ref:'T',alt:'C',init:7,refp:'MKLK*',altp:'MK*'}
  ]) c
), results AS (
  SELECT c, duckvep_product_scenario_v1({
    reference:{start1:1::UBIGINT,dna:c.dna,strand:1::TINYINT},
    exons:[{start1:1::UBIGINT,end1:length(c.dna)::UBIGINT}],
    blocks:[{start1:1::UBIGINT,end1:length(c.dna)::UBIGINT}],
    sources:[{source_id:0::UBIGINT,start1:c.pos::UBIGINT,reference:c.ref,alternate:c.alt}],
    programme:{id:c.id,evidence_id:'declared',assumption_id:'joint-assignment',
      expected_reference_rna:c.dna,annotated_start1:1::UINTEGER,
      alternative_start1:c.init::UINTEGER,
      alternative_trigger_id:CASE WHEN c.init>0 THEN 'annotated-start-lost' ELSE '' END,
      genetic_code:1::UTINYINT}
  }) p FROM cases
)
SELECT c.id,
  CASE WHEN p.reference_peptide=c.refp AND p.alternate_peptide=c.altp
        AND p.conditional_nmd='unknown' AND p.allele_evidence='hypothetical_assignment'
       THEN 'exact' ELSE error('conditional product translation mismatch: ' || c.id) END AS verified,
  to_json(p) AS full_typed_fields
FROM results ORDER BY c.id;
