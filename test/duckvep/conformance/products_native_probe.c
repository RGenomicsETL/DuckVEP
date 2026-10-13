#include "duckvep_product_profile.h"
#include <stdio.h>
#include <string.h>

/* A small stdin-free bridge: every source and RNA assumption is an argv value.
 * The product implementation, not this bridge, owns mutation and translation. */
int main(int argc, char **argv) {
    if (argc != 6) return 2;
    const uint8_t *dna=(const uint8_t *)argv[1];
    const uint8_t *expected=(const uint8_t *)argv[3];
    const uint8_t *ref=(const uint8_t *)argv[4];
    const uint8_t *alt=(const uint8_t *)argv[5];
    duckvep_profile_interval_t exons[]={{1,9},{13,15},{19,27}};
    duckvep_profile_rna_block_v1_t blocks[3]={{{1,9},1},{{19,27},10}};
    if (!strcmp(argv[2],"include")) {
        blocks[1]=(duckvep_profile_rna_block_v1_t){{13,15},10};
        blocks[2]=(duckvep_profile_rna_block_v1_t){{19,27},13};
    } else if (strcmp(argv[2],"skip")) return 2;
    duckvep_profile_v1_input_t geometry={
        .reference_start1=1,.reference_dna=dna,.dna_length=strlen(argv[1]),
        .transcript={1,27},.exons=exons,.exon_count=3,.strand=1,.model_complete=1};
    duckvep_profile_programme_v1_t p={
        .identity=argv[2],.evidence_id="user-mapped",.assumption_id=argv[2],
        .blocks=blocks,.block_count=!strcmp(argv[2],"include")?3u:2u,
        .reference_rna=expected,.rna_length=strlen(argv[3]),
        .cds_start1=1,.cds_stop_end1=(uint32_t)strlen(argv[3]),.annotated_start1=1};
    duckvep_haplotype_edit_t edit={7,(uint32_t)strlen(argv[4]),ref,
                                  (uint32_t)strlen(argv[5]),alt,1};
    uint64_t source_id=11;
    duckvep_profile_execution_t query={&geometry,&edit,&source_id,1,&p,1};
    uint8_t dna_out[32],ref_rna[32],alt_rna[32],ref_peptide[32],alt_peptide[32];
    duckvep_profile_source_effect_t effect;
    duckvep_profile_workspace_t w={dna_out,ref_rna,alt_rna,ref_peptide,alt_peptide,
                                   32,32,32,&effect,1};
    duckvep_profile_execution_result_t result;
    duckvep_profile_status_t status=duckvep_product_execute_v1(&query,&w,&result);
    if (status != DUCKVEP_PROFILE_OK) {
        printf("status=%d\n",status);
        return 0;
    }
    printf("status=0 product=%s ref_rna=%.*s alt_rna=%.*s ref_peptide=%.*s alt_peptide=%.*s ref_stop=%zu alt_stop=%zu junction=%u source=%llu/rna%u\n",
        result.identity,(int)result.rna_length,ref_rna,(int)result.rna_length,alt_rna,
        (int)result.reference_peptide_length,ref_peptide,
        (int)result.alternate_peptide_length,alt_peptide,
        result.reference_stop1,result.alternate_stop1,result.last_junction1,
        (unsigned long long)effect.source_id,effect.rna_first1);
    return 0;
}
