CREATE TEMP VIEW product_fixture AS SELECT {
 reference:{start1:1::UBIGINT,dna:'ATGAAACCCGGGAAAGGGTTTCCCTAAGGGATGAAATAAGGGAAA',strand:1::TINYINT},
 exons:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],
 blocks:[{start1:1::UBIGINT,end1:9::UBIGINT},{start1:13::UBIGINT,end1:15::UBIGINT},{start1:19::UBIGINT,end1:27::UBIGINT},{start1:31::UBIGINT,end1:39::UBIGINT},{start1:43::UBIGINT,end1:45::UBIGINT}],
 sources:[{source_id:11::UBIGINT,start1:13::UBIGINT,reference:'A',alternate:'C'},
          {source_id:22::UBIGINT,start1:29::UBIGINT,reference:'G',alternate:'C'}],
 programme:{id:'include',evidence_id:'declared',assumption_id:'splice',
            expected_reference_rna:'ATGAAACCCAAATTTCCCTAAATGAAATAAAAA',
            annotated_start1:1::UINTEGER,alternative_start1:0::UINTEGER,
            alternative_trigger_id:'',genetic_code:1::UTINYINT}
} AS s;
SELECT (duckvep_product_scenario_v1(s)).verified_sources AS valid FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := []::STRUCT(source_id UBIGINT,start1 UBIGINT,reference VARCHAR,alternate VARCHAR)[]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := [NULL::STRUCT(source_id UBIGINT,start1 UBIGINT,reference VARCHAR,alternate VARCHAR)]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  exons := []::STRUCT(start1 UBIGINT,end1 UBIGINT)[])) FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  blocks := [NULL::STRUCT(start1 UBIGINT,end1 UBIGINT)])) FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := [struct_update(s.sources[1],reference := NULL::VARCHAR),s.sources[2]]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  exons := [struct_update(s.exons[1],end1 := 0::UBIGINT),s.exons[2],s.exons[3],s.exons[4],s.exons[5]]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := NULL::STRUCT(source_id UBIGINT,start1 UBIGINT,reference VARCHAR,alternate VARCHAR)[]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := [s.sources[1],struct_update(s.sources[2],source_id := 11::UBIGINT)]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := [s.sources[1],struct_update(s.sources[2],start1 := 13::UBIGINT)]))
  FROM product_fixture;
SELECT duckvep_product_scenario_v1(struct_update(s,
  sources := [struct_update(s.sources[1],alternate := 'CC'),s.sources[2]]))
  FROM product_fixture;
SELECT (duckvep_product_scenario_v1(s)).verified_sources AS valid_after_errors FROM product_fixture LIMIT 1;
SELECT (duckvep_product_scenario_v1(struct_update(s, sources := s.sources[:1]))).verified_sources AS one_source,
       (duckvep_product_scenario_v1(s)).verified_sources AS two_sources FROM product_fixture;
