#include "duckvep_property.h"
#include "core/duckvep_core_products.h"
#include "duckvep_budget.h"
#include <string.h>

static const duckvep_product_interval_t exons[]={{1,9},{13,15},{19,27},{31,39},{43,45}};
static const duckvep_product_interval_t skip[]={{1,9},{19,27},{31,39},{43,45}};
static const duckvep_product_interval_t retain[]={{1,9},{10,12},{13,15},{19,27},{31,39},{43,45}};
static const duckvep_product_source_t sources[]={
    {11,13,{"A",1},{"C",1}}, {22,29,{"G",1},{"C",1}}};
static duckvep_core_product_fields_t scenario(void){
    return (duckvep_core_product_fields_t){
        .reference_start1=1,.dna={"ATGAAACCCGGGAAAGGGTTTCCCTAAGGGATGAAATAAGGGAAA",45},
        .strand=1,.exons=exons,.exon_count=5,.blocks=exons,.block_count=5,
        .sources=sources,.source_count=2,
        .id={"include",7},.evidence_id={"declared",8},.assumption_id={"splice",6},
        .expected_rna={"ATGAAACCCAAATTTCCCTAAATGAAATAAAAA",33},
        .annotated_start1=1,.genetic_code=1
    };
}
TEST product_core_owner_fault_reuse_and_geometry(void){
    duckvep_core_product_fields_t f=scenario(), original=f;
    duckvep_core_product_result_t out={0};
    duckvep_budget_stats_t before,after;
    duckvep_budget_stats(&before);
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    ASSERT_EQ(2u,out.result.verified_sources);
    ASSERT_EQ(11u,out.effects[0].source_id);
    ASSERT_EQ(1u,out.effects[0].changed_rna_bases);
    ASSERT_EQ(22u,out.effects[1].source_id);
    ASSERT_EQ(0u,out.effects[1].changed_rna_bases);
    ASSERT(memcmp(&f,&original,sizeof f)==0);
    duckvep_core_product_dispose(&out);
    duckvep_budget_stats(&after);
    ASSERT_EQ(before.total_current,after.total_current);
    duckvep_fault_arm(UINT64_MAX);
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    uint64_t attempts=duckvep_fault_allocations();
    duckvep_core_product_dispose(&out);
    duckvep_fault_arm(0);
    for(uint64_t n=1;n<=attempts;++n){
        duckvep_fault_arm(n);
        ASSERT_EQ(DUCKVEP_PROFILE_CAPACITY,duckvep_core_product_scenario(&f,&out));
        ASSERT_EQ(1u,duckvep_fault_fired());
        ASSERT_EQ(NULL,out.dna);
        duckvep_budget_stats(&after);
        ASSERT_EQ(before.total_current,after.total_current);
        duckvep_fault_arm(0);
    }
    uint64_t original_limit=duckvep_budget_limit();
    duckvep_budget_reset_high_water();
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    duckvep_budget_stats(&after);
    uint64_t exact=after.total_high_water;
    duckvep_core_product_dispose(&out);
    ASSERT(exact>before.total_current);
    ASSERT(duckvep_budget_set_limit(exact));
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    duckvep_core_product_dispose(&out);
    ASSERT(duckvep_budget_set_limit(exact-1u));
    ASSERT_EQ(DUCKVEP_PROFILE_CAPACITY,duckvep_core_product_scenario(&f,&out));
    ASSERT_EQ(NULL,out.dna);
    ASSERT(duckvep_budget_set_limit(original_limit));
    duckvep_product_source_t swapped[]={sources[1],sources[0]};
    f.sources=swapped;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    ASSERT_EQ(22u,out.effects[0].source_id);
    ASSERT_EQ(11u,out.effects[1].source_id);
    duckvep_core_product_dispose(&out);
    f=scenario();f.blocks=skip;f.block_count=4;
    f.expected_rna=(duckvep_product_text_t){"ATGAAACCCTTTCCCTAAATGAAATAAAAA",30};
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    ASSERT_EQ(0u,out.effects[0].changed_rna_bases);
    duckvep_core_product_dispose(&out);
    f=scenario();f.blocks=retain;f.block_count=6;
    f.expected_rna=(duckvep_product_text_t){"ATGAAACCCGGGAAATTTCCCTAAATGAAATAAAAA",36};
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    duckvep_core_product_dispose(&out);
    f=scenario();
    duckvep_product_source_t zero_sources[]={sources[0],sources[1]};
    zero_sources[0].source_id=0;
    f.sources=zero_sources;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_core_product_scenario(&f,&out));
    ASSERT_EQ(0u,out.effects[0].source_id);
    ASSERT_EQ(22u,out.effects[1].source_id);
    duckvep_core_product_dispose(&out);
    zero_sources[1].source_id=0;
    ASSERT_EQ(DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,duckvep_core_product_scenario(&f,&out));
    f=scenario();
    duckvep_product_source_t duplicated[]={sources[0],sources[0]};
    f.sources=duplicated;
    ASSERT_EQ(DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,duckvep_core_product_scenario(&f,&out));
    ASSERT_EQ(NULL,out.dna);
    duckvep_budget_stats(&after);
    ASSERT_EQ(before.total_current,after.total_current);
    PASS();
}
