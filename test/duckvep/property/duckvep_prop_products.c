#include "duckvep_property.h"
#include "duckvep_product_profile.h"

typedef struct {
    uint8_t dna[161], rna[161], alt_dna[161], ref_rna[161], alt_rna[161];
    uint8_t ref_peptide[64], alt_peptide[64];
    duckvep_profile_interval_t exons[3];
    duckvep_profile_rna_block_v1_t blocks[3];
    duckvep_profile_v1_input_t geometry;
    duckvep_profile_programme_v1_t programme;
    duckvep_haplotype_edit_t edits[2];
    uint64_t ids[2];
    duckvep_profile_source_effect_t effects[2];
    duckvep_profile_execution_t request;
    duckvep_profile_workspace_t workspace;
    duckvep_profile_execution_result_t result;
    duckvep_profile_v1_fact_t facts[8];
    size_t count;
} product_scene;

static void product_setup(product_scene *s, int strand, int include_exon) {
    memset(s, 0, sizeof *s);
    memset(s->dna, 'A', sizeof s->dna);
    s->exons[0] = (duckvep_profile_interval_t){1, 90};
    s->exons[1] = (duckvep_profile_interval_t){101, 121};
    s->exons[2] = (duckvep_profile_interval_t){132, 161};
    memset(s->dna + 100, 'C', 21);
    memcpy(s->dna, "ATGATG", 6);
    memcpy(s->dna + 131, "TAA", 3);
    if (strand == -1) {
        memset(s->dna, 'T', 90);
        memset(s->dna + 131, 'T', 30);
        memcpy(s->dna + 158, "CAT", 3);
        memcpy(s->dna + 87, "TTA", 3); /* reverse complement TAA */
    }
    s->geometry = (duckvep_profile_v1_input_t){
        .reference_start1=1, .reference_dna=s->dna, .dna_length=161,
        .transcript={1,161}, .exons=s->exons, .exon_count=3,
        .strand=(int8_t)strand, .model_complete=1};
    if (strand == 1) {
        s->blocks[0]=(duckvep_profile_rna_block_v1_t){{1,90},1};
        s->blocks[1]=(duckvep_profile_rna_block_v1_t){{101,121},91};
        s->blocks[2]=(duckvep_profile_rna_block_v1_t){{132,161},112};
        if (!include_exon) s->blocks[1]=(duckvep_profile_rna_block_v1_t){{132,161},91};
    } else {
        s->blocks[0]=(duckvep_profile_rna_block_v1_t){{132,161},1};
        s->blocks[1]=(duckvep_profile_rna_block_v1_t){{101,121},31};
        s->blocks[2]=(duckvep_profile_rna_block_v1_t){{1,90},52};
        if (!include_exon) s->blocks[1]=(duckvep_profile_rna_block_v1_t){{1,90},31};
    }
    size_t count = include_exon ? 3u : 2u;
    size_t offset = 0;
    for (size_t i=0;i<count;i++) {
        duckvep_profile_interval_t b=s->blocks[i].genomic_source;
        for (uint32_t p=b.start1;p<=b.end1;p++) {
            uint8_t v=s->dna[strand==1 ? p-1 : b.end1-(p-b.start1)-1];
            if (strand==-1) v=v=='A'?'T':v=='T'?'A':v=='C'?'G':'C';
            s->rna[offset++]=v;
        }
    }
    s->programme=(duckvep_profile_programme_v1_t){
        .identity=include_exon?"included":"skipped", .evidence_id="declared-source",
        .assumption_id=include_exon?"include-exon":"skip-exon",
        .blocks=s->blocks, .block_count=count, .reference_rna=s->rna,
        .rna_length=offset, .cds_start1=1, .cds_stop_end1=(uint32_t)offset,
        .annotated_start1=1};
    s->workspace=(duckvep_profile_workspace_t){
        s->alt_dna,s->ref_rna,s->alt_rna,s->ref_peptide,s->alt_peptide,
        sizeof s->alt_dna,sizeof s->alt_rna,sizeof s->ref_peptide,
        s->effects,2};
    s->request=(duckvep_profile_execution_t){
        &s->geometry,s->edits,s->ids,0,&s->programme,DUCKVEP_CODON_TABLE_STANDARD};
}

TEST product_two_named_mappings_and_junction_geometry(void) {
    for (int strand=-1;strand<=1;strand+=2) {
        product_scene s;
        product_setup(&s,strand,0);
        uint32_t pos=strand==1?22:150;
        static const uint8_t stop_ref[]="AAA", stop_alt[]="TAA";
        static const uint8_t minus_ref[]="TTT", minus_alt[]="TTA";
        /* Reverse genomic TTT at 150..152 becomes an RNA stop with TTA. */
        s.edits[0]=(duckvep_haplotype_edit_t){pos,3,
            strand==1?stop_ref:minus_ref,3,strand==1?stop_alt:minus_alt,1};
        s.ids[0]=37; s.request.source_count=1;
        ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_position_v1(
            &s.request,&s.workspace,&s.result,s.facts,8,&s.count));
        ASSERT_EQ(1u,s.count);
        ASSERT_EQ(strand==1?22u:152u,s.facts[0].pos1);
        ASSERT_EQ(strand==1?90u:30u,s.result.last_junction1);
        uint8_t skipped[161]; memcpy(skipped,s.alt_rna,s.result.rna_length);
        uint8_t skip_peptide[64]; memcpy(skip_peptide,s.alt_peptide,s.result.alternate_peptide_length);
        size_t skip_peptide_len=s.result.alternate_peptide_length;
        size_t skip_len=s.result.rna_length;
        product_setup(&s,strand,1);
        s.edits[0]=(duckvep_haplotype_edit_t){pos,3,
            strand==1?stop_ref:minus_ref,3,strand==1?stop_alt:minus_alt,1};
        s.ids[0]=37; s.request.source_count=1;
        ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
        ASSERT_EQ(skip_len+21,s.result.rna_length);
        ASSERT(memcmp(skipped,s.alt_rna,skip_len)!=0);
        ASSERT_EQ(skip_peptide_len,s.result.alternate_peptide_length);
        ASSERT_EQ(0,memcmp(skip_peptide,s.alt_peptide,skip_peptide_len));
        ASSERT_EQ(strand==1?8u:4u,s.result.alternate_peptide_length);
        ASSERT_EQ(0,strcmp((const char *)s.alt_peptide,strand==1?"MMKKKKK*":"MKK*"));
        ASSERT_EQ(s.result.reference_stop1,s.result.reference_peptide_length);
        ASSERT_EQ(0u,s.ref_peptide[s.result.reference_peptide_length]);
        ASSERT_EQ(1u,s.result.verified_sources);
        ASSERT_EQ(37u,s.effects[0].source_id);
        ASSERT_EQ(1u,s.effects[0].changed_rna_bases);
        if (strand==1) {
            ASSERT_EQ(8u,s.result.alternate_stop1);
            ASSERT_EQ(111u,s.result.last_junction1);
        }
    }
    PASS();
}

TEST product_stop_and_junction_geometry(void) {
    product_scene s;
    product_setup(&s,1,0);
    memset(s.dna,'A',90);
    memcpy(s.dna+1,"ATG",3);
    memcpy(s.rna,s.dna,90);
    s.programme.cds_start1=2;
    s.programme.annotated_start1=2;
    static const uint8_t ref[]="AAA", alt[]="TAA";
    s.edits[0]=(duckvep_haplotype_edit_t){38,3,ref,3,alt,1};
    s.ids[0]=17;
    s.request.source_count=1;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(90u,s.result.last_junction1);
    ASSERT_EQ(13u,s.result.alternate_stop1);
    s.edits[0].cds_start=35;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(12u,s.result.alternate_stop1);
    memcpy(s.dna,"ATG",3);
    memcpy(s.rna,s.dna,90);
    s.programme.cds_start1=1;
    s.programme.annotated_start1=1;
    s.edits[0].cds_start=37;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(13u,s.result.alternate_stop1);
    PASS();
}

TEST product_exon_specific_source_effects(void) {
    product_scene s;
    product_setup(&s,1,0);
    static const uint8_t ref[]="AAA", alt[]="TAA";
    static const uint8_t middle_ref[]="C", middle_alt[]="T";
    s.edits[0]=(duckvep_haplotype_edit_t){101,1,middle_ref,1,middle_alt,1};
    s.edits[1]=(duckvep_haplotype_edit_t){22,3,ref,3,alt,1};
    s.ids[0]=71; s.ids[1]=37;
    s.request.source_count=2;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(71u,s.effects[0].source_id);
    ASSERT_EQ(1u,s.effects[0].changed_dna);
    ASSERT_EQ(0u,s.effects[0].changed_rna_bases);
    ASSERT_EQ(0u,s.effects[0].rna_first1);
    ASSERT_EQ(22u,s.effects[1].rna_first1);
    s.programme.identity="included";
    s.programme.assumption_id="include-exon";
    s.blocks[1]=(duckvep_profile_rna_block_v1_t){{101,121},91};
    s.blocks[2]=(duckvep_profile_rna_block_v1_t){{132,161},112};
    s.programme.block_count=3;
    s.programme.rna_length=141;
    s.programme.cds_stop_end1=141;
    s.workspace.rna_capacity=161;
    for (size_t i=0;i<21;i++) s.rna[90+i]='C';
    for (size_t i=0;i<30;i++) s.rna[111+i]=s.dna[131+i];
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(91u,s.effects[0].rna_first1);
    ASSERT_EQ(1u,s.effects[0].changed_rna_bases);
    ASSERT_EQ(22u,s.effects[1].rna_first1);
    PASS();
}

TEST product_start_loss_uses_only_declared_alternative(void) {
    product_scene s; product_setup(&s,1,0);
    static const uint8_t ref[]="T", alt[]="C";
    s.edits[0]=(duckvep_haplotype_edit_t){2,1,ref,1,alt,1};
    s.ids[0]=1; s.request.source_count=1;
    s.programme.alternative_start1=4;
    s.programme.trigger_id="annotated-start-lost";
    s.programme.counterfactual_id="reference-annotated-start-retained";
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(1u,s.result.annotated_start_lost);
    ASSERT_EQ(1u,s.result.alternative_start_used);
    ASSERT_EQ(1u,s.result.reference_start1);
    ASSERT_EQ(4u,s.result.alternate_start1);
    ASSERT_EQ('M',s.ref_peptide[0]);
    ASSERT_EQ('M',s.alt_peptide[0]);
    s.programme.trigger_id="scan-kozak";
    duckvep_profile_execution_result_t before=s.result;
    ASSERT_EQ(DUCKVEP_PROFILE_UNSUPPORTED_PRODUCT,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(&before,&s.result,sizeof before));
    s.programme.trigger_id="annotated-start-lost";
    s.programme.alternative_start1=7;
    ASSERT_EQ(DUCKVEP_PROFILE_UNSUPPORTED_PRODUCT,duckvep_product_execute_v1(
        &s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(&before,&s.result,sizeof before));
    s.programme.alternative_start1=4;
    static const uint8_t other_ref[]="A", other_alt[]="C";
    s.edits[0]=(duckvep_haplotype_edit_t){22,1,other_ref,1,other_alt,1};
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(0u,s.result.alternative_start_used);
    PASS();
}

TEST product_sources_capacity_atomic_and_reuse(void) {
    product_scene s; product_setup(&s,1,0);
    static const uint8_t ref[]="AAA", alt[]="TAA", wrong[]="CCC";
    s.edits[0]=(duckvep_haplotype_edit_t){22,3,ref,3,alt,1};
    s.ids[0]=99; s.request.source_count=1;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    duckvep_profile_execution_result_t before=s.result;
    duckvep_profile_source_effect_t before_effect=s.effects[0];
    s.edits[0].ref=wrong;
    ASSERT_EQ(DUCKVEP_PROFILE_REFERENCE_MISMATCH,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(&before,&s.result,sizeof before));
    ASSERT_EQ(0,memcmp(&before_effect,&s.effects[0],sizeof before_effect));
    s.edits[0].ref=ref;
    s.workspace.rna_capacity=119;
    ASSERT_EQ(DUCKVEP_PROFILE_SHORT_BUFFER,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    s.workspace.rna_capacity=161;
    s.workspace.source_effect_capacity=0;
    ASSERT_EQ(DUCKVEP_PROFILE_SHORT_BUFFER,duckvep_product_execute_v1(
        &s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(&before_effect,&s.effects[0],sizeof before_effect));
    s.workspace.source_effect_capacity=2;
    s.workspace.peptide_capacity=1;
    ASSERT_EQ(DUCKVEP_PROFILE_SHORT_BUFFER,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    s.workspace.peptide_capacity=64;
    s.workspace.rna_capacity=120;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    s.workspace.rna_capacity=161;
    s.request.source_count=2;
    s.edits[1]=s.edits[0]; s.ids[1]=100;
    ASSERT_EQ(DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    s.request.source_count=1;
    s.rna[21]='G';
    ASSERT_EQ(DUCKVEP_PROFILE_REFERENCE_MISMATCH,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    s.rna[21]='A';
    s.dna[70]='N';
    ASSERT_EQ(DUCKVEP_PROFILE_SEQUENCE_UNAVAILABLE,duckvep_product_execute_v1(
        &s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(&before,&s.result,sizeof before));
    s.dna[70]='A';
    s.blocks[1].rna_start1=92;
    ASSERT_EQ(DUCKVEP_PROFILE_PARTIAL_MODEL,duckvep_product_execute_v1(
        &s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(&before,&s.result,sizeof before));
    s.blocks[1].rna_start1=91;
    s.geometry.reference_witness=(const uint8_t *)"not-an-authority";
    s.geometry.final_dna=(const uint8_t *)"not-an-authority";
    uint8_t immutable[161]; memcpy(immutable,s.dna,sizeof immutable);
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_v1(&s.request,&s.workspace,&s.result));
    ASSERT_EQ(0,memcmp(immutable,s.dna,sizeof immutable));
    ASSERT_EQ(0,memcmp(&before,&s.result,sizeof before));
    s.facts[0].pos1=UINT32_MAX;
    ASSERT_EQ(DUCKVEP_PROFILE_SHORT_BUFFER,duckvep_product_execute_position_v1(
        &s.request,&s.workspace,&s.result,s.facts,0,&s.count));
    ASSERT_EQ(0u,s.count);
    ASSERT_EQ(UINT32_MAX,s.facts[0].pos1);
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_position_v1(
        &s.request,&s.workspace,&s.result,s.facts,1,&s.count));
    ASSERT_EQ(1u,s.count);
    static const uint8_t other_ref[]="A", other_alt[]="C";
    s.edits[1]=(duckvep_haplotype_edit_t){4,1,other_ref,1,other_alt,1};
    s.ids[1]=101; s.request.source_count=2;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_execute_position_v1(
        &s.request,&s.workspace,&s.result,s.facts,2,&s.count));
    ASSERT_EQ(2u,s.count);
    ASSERT_EQ(2u,s.result.verified_sources);
    ASSERT_EQ(101u,s.effects[1].source_id);
    ASSERT_EQ(4u,s.effects[1].rna_first1);
    ASSERT_EQ(4u,s.facts[0].pos1);
    ASSERT_EQ(22u,s.facts[1].pos1);
    s.edits[1].cds_start=24;
    ASSERT_EQ(DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,duckvep_product_execute_v1(
        &s.request,&s.workspace,&s.result));
    s.edits[1].cds_start=4;
    s.ids[1]=99;
    ASSERT_EQ(DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,duckvep_product_execute_v1(
        &s.request,&s.workspace,&s.result));
    PASS();
}

TEST product_position_fact_contract(void) {
    const uint8_t reference[]="ACGTACGT", changed[]="ATGTATGT";
    const duckvep_profile_interval_t exons[]={{1,2},{5,8}};
    const duckvep_profile_feature_t features[]={{{2,3},7}};
    duckvep_profile_v1_input_t in={
        .reference_start1=1,.reference_dna=reference,.dna_length=8,
        .reference_witness=reference,.witness_length=8,.final_dna=reference,.final_length=8,
        .transcript={1,8},.strand=1,.exons=exons,.exon_count=2,.model_complete=1,
        .features=features,.feature_count=1,.joint_status=DUCKVEP_PROFILE_OK};
    duckvep_profile_v1_fact_t facts[4]={{999,0,0}};
    size_t n=99;
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_position_v1(&in,facts,4,&n));
    ASSERT_EQ(0u,n);
    in.final_dna=changed;
    ASSERT_EQ(DUCKVEP_PROFILE_SHORT_BUFFER,duckvep_product_position_v1(&in,facts,2,&n));
    ASSERT_EQ(0u,n);
    ASSERT_EQ(999u,facts[0].pos1);
    ASSERT_EQ(DUCKVEP_PROFILE_OK,duckvep_product_position_v1(&in,facts,4,&n));
    ASSERT_EQ(3u,n);
    ASSERT_EQ(2u,facts[0].pos1);
    ASSERT_EQ(0u,facts[0].flags);
    ASSERT_EQ(7u,facts[1].feature_id);
    ASSERT_EQ(DUCKVEP_PROFILE_REGULATORY_INTERVAL,facts[1].flags);
    ASSERT_EQ(6u,facts[2].pos1);
    ASSERT_EQ(0u,facts[2].flags);
    in.reference_witness=changed;
    ASSERT_EQ(DUCKVEP_PROFILE_REFERENCE_MISMATCH,duckvep_product_position_v1(&in,facts,4,&n));
    in.reference_witness=reference;in.model_complete=0;
    ASSERT_EQ(DUCKVEP_PROFILE_PARTIAL_MODEL,duckvep_product_position_v1(&in,facts,4,&n));
    in.model_complete=1;in.joint_status=DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION;
    ASSERT_EQ(DUCKVEP_PROFILE_AMBIGUOUS_COMPOSITION,duckvep_product_position_v1(&in,facts,4,&n));
    PASS();
}
