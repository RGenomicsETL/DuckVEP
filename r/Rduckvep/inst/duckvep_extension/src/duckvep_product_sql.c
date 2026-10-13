#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "core/duckvep_core_products.h"
#include "kernel/src/duckvep_budget.h"
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

static bool valid(duckdb_vector v, idx_t row) {
    uint64_t *bits = duckdb_vector_get_validity(v);
    return !bits || duckdb_validity_row_is_valid(bits, row);
}
static duckvep_product_text_t text_at(duckdb_vector v, idx_t row) {
    duckdb_string_t *str = &((duckdb_string_t *)duckdb_vector_get_data(v))[row];
    return (duckvep_product_text_t){duckdb_string_t_data(str), duckdb_string_t_length(*str)};
}
static bool list_span(duckdb_vector v, idx_t row, idx_t *start, size_t *count) {
    if (!valid(v, row)) return false;
    duckdb_list_entry e = ((duckdb_list_entry *)duckdb_vector_get_data(v))[row];
    idx_t size = duckdb_list_vector_get_size(v);
    if (e.offset > size || e.length > size - e.offset || e.length > SIZE_MAX) return false;
    *start = e.offset; *count = (size_t)e.length;
    return true;
}
static bool intervals(duckdb_vector list, idx_t row, duckvep_product_interval_t *dst,
                      idx_t *start, size_t *count) {
    if (!list_span(list, row, start, count) || !*count || !dst) return false;
    duckdb_vector child = duckdb_list_vector_get_child(list);
    duckdb_vector a = duckdb_struct_vector_get_child(child, 0);
    duckdb_vector b = duckdb_struct_vector_get_child(child, 1);
    for (size_t i = 0; i < *count; ++i) {
        idx_t at = *start + i;
        if (!valid(child, at) || !valid(a, at) || !valid(b, at)) return false;
        dst[i] = (duckvep_product_interval_t){((uint64_t *)duckdb_vector_get_data(a))[at],
                                               ((uint64_t *)duckdb_vector_get_data(b))[at]};
    }
    return true;
}
static void product_scalar(duckdb_function_info info, duckdb_data_chunk chunk, duckdb_vector output) {
    duckdb_vector input = duckdb_data_chunk_get_vector(chunk, 0);
    duckdb_vector reference = duckdb_struct_vector_get_child(input, 0);
    duckdb_vector exon_list = duckdb_struct_vector_get_child(input, 1);
    duckdb_vector block_list = duckdb_struct_vector_get_child(input, 2);
    duckdb_vector source_list = duckdb_struct_vector_get_child(input, 3);
    duckdb_vector programme = duckdb_struct_vector_get_child(input, 4);
    duckdb_vector ref[3], prog[8], results[PRODUCT_OUTPUT_COUNT];
    for (size_t i=0;i<3;++i) ref[i]=duckdb_struct_vector_get_child(reference,i);
    for (size_t i=0;i<8;++i) prog[i]=duckdb_struct_vector_get_child(programme,i);
    for (size_t i=0;i<PRODUCT_OUTPUT_COUNT;++i) results[i]=duckdb_struct_vector_get_child(output,i);
    idx_t written=0;
    for (idx_t row=0;row<duckdb_data_chunk_get_size(chunk);++row) {
        duckvep_core_product_fields_t f={0};
        duckvep_product_interval_t *exons=NULL,*blocks=NULL;
        duckvep_product_source_t *sources=NULL;
        duckvep_core_product_result_t owned={0};
        const char *failure="invalid scenario shape or source";
        idx_t first=0; size_t count=0;
        if (!valid(input,row)||!valid(reference,row)||!valid(programme,row)) goto fail;
        for(size_t i=0;i<3;++i) if(!valid(ref[i],row)) goto fail;
        for(size_t i=0;i<8;++i) if(!valid(prog[i],row)) goto fail;
        f.reference_start1=((uint64_t *)duckdb_vector_get_data(ref[0]))[row];
        f.dna=text_at(ref[1],row);
        f.strand=((int8_t *)duckdb_vector_get_data(ref[2]))[row];
        if(!list_span(exon_list,row,&first,&count)||!count||count>SIZE_MAX/sizeof *exons) goto fail;
        exons=duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,count*sizeof *exons);
        if(!exons) {failure="allocation limit";goto fail;}
        if(!intervals(exon_list,row,exons,&first,&count)) goto fail;
        f.exons=exons;f.exon_count=count;
        if(!list_span(block_list,row,&first,&count)||!count||count>SIZE_MAX/sizeof *blocks) goto fail;
        blocks=duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,count*sizeof *blocks);
        if(!blocks) {failure="allocation limit";goto fail;}
        if(!intervals(block_list,row,blocks,&first,&count)) goto fail;
        f.blocks=blocks;f.block_count=count;
        if(!list_span(source_list,row,&first,&count)||!count||count>SIZE_MAX/sizeof *sources) goto fail;
        sources=duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,count*sizeof *sources);
        if(!sources) {failure="allocation limit";goto fail;}
        duckdb_vector sv=duckdb_list_vector_get_child(source_list);
        duckdb_vector sc[4];for(size_t i=0;i<4;++i) sc[i]=duckdb_struct_vector_get_child(sv,i);
        for(size_t i=0;i<count;++i){
            idx_t at=first+i;
            if(!valid(sv,at)) goto fail;
            for(size_t j=0;j<4;++j) if(!valid(sc[j],at)) goto fail;
            sources[i]=(duckvep_product_source_t){
                ((uint64_t *)duckdb_vector_get_data(sc[0]))[at],
                ((uint64_t *)duckdb_vector_get_data(sc[1]))[at],
                text_at(sc[2],at),text_at(sc[3],at)};
        }
        f.sources=sources;f.source_count=count;
        f.id=text_at(prog[0],row);f.evidence_id=text_at(prog[1],row);
        f.assumption_id=text_at(prog[2],row);f.expected_rna=text_at(prog[3],row);
        f.annotated_start1=((uint32_t *)duckdb_vector_get_data(prog[4]))[row];
        f.alternative_start1=((uint32_t *)duckdb_vector_get_data(prog[5]))[row];
        f.trigger_id=text_at(prog[6],row);f.genetic_code=((uint8_t *)duckdb_vector_get_data(prog[7]))[row];
        duckvep_profile_status_t status=duckvep_core_product_scenario(&f,&owned);
        if(status!=DUCKVEP_PROFILE_OK){
            const char *messages[]={"ok","invalid scenario shape or source","genomic or RNA reference mismatch",
                "sequence unavailable","ambiguous composition","partial model or programme",
                "short workspace","allocation limit","unsupported product or initiation"};
            failure=messages[status];goto fail;
        }
        for(size_t which=0;which<2;++which){
            duckdb_vector list=results[which];
            if(count>UINT64_MAX-written || duckdb_list_vector_reserve(list,written+count)!=DuckDBSuccess ||
               duckdb_list_vector_set_size(list,written+count)!=DuckDBSuccess){failure="output list allocation failed";goto fail;}
            ((duckdb_list_entry *)duckdb_vector_get_data(list))[row]=(duckdb_list_entry){written,count};
            duckdb_vector item=duckdb_list_vector_get_child(list);
            duckdb_vector part[7];size_t fields=which?7:4;
            for(size_t j=0;j<fields;++j) part[j]=duckdb_struct_vector_get_child(item,j);
            for(size_t i=0;i<count;++i){
                idx_t at=written+i;
                if(!which){
                    ((uint64_t *)duckdb_vector_get_data(part[0]))[at]=sources[i].source_id;
                    ((uint64_t *)duckdb_vector_get_data(part[1]))[at]=sources[i].start1;
                    duckdb_vector_assign_string_element_len(part[2],at,sources[i].reference.ptr,sources[i].reference.length);
                    duckdb_vector_assign_string_element_len(part[3],at,sources[i].alternate.ptr,sources[i].alternate.length);
                }else{
                    duckvep_core_product_effect_t *e=&owned.effects[i];
                    ((uint64_t *)duckdb_vector_get_data(part[0]))[at]=e->source_id;
                    ((bool *)duckdb_vector_get_data(part[1]))[at]=e->changed_dna!=0;
                    ((uint32_t *)duckdb_vector_get_data(part[2]))[at]=e->changed_rna_bases;
                    ((uint64_t *)duckdb_vector_get_data(part[3]))[at]=e->genomic_start1;
                    ((uint64_t *)duckdb_vector_get_data(part[4]))[at]=e->genomic_end1;
                    ((uint32_t *)duckdb_vector_get_data(part[5]))[at]=e->rna_first1;
                    ((uint32_t *)duckdb_vector_get_data(part[6]))[at]=e->rna_last1;
                }
            }
        }
        written+=count;
        for(size_t i=2;i<PRODUCT_OUTPUT_COUNT;++i){
            const char *text;size_t length;uint64_t number;
            (void)duckvep_core_product_output(&f,&owned,i,&text,&length,&number);
            if(text) duckdb_vector_assign_string_element_len(results[i],row,text,length);
            else if(!strcmp(duckvep_product_output_types[i],"UBIGINT"))
                ((uint64_t *)duckdb_vector_get_data(results[i]))[row]=number;
            else if(!strcmp(duckvep_product_output_types[i],"TINYINT"))
                ((int8_t *)duckdb_vector_get_data(results[i]))[row]=(int8_t)number;
            else if(!strcmp(duckvep_product_output_types[i],"UTINYINT"))
                ((uint8_t *)duckdb_vector_get_data(results[i]))[row]=(uint8_t)number;
            else ((uint32_t *)duckdb_vector_get_data(results[i]))[row]=(uint32_t)number;
        }
        duckvep_core_product_dispose(&owned);
        duckvep_budget_free(exons);duckvep_budget_free(blocks);duckvep_budget_free(sources);
        continue;
fail:
        duckvep_core_product_dispose(&owned);
        duckvep_budget_free(exons);duckvep_budget_free(blocks);duckvep_budget_free(sources);
        duckdb_scalar_function_set_error(info,failure);
        return;
    }
}
static duckdb_logical_type record(const duckdb_type *types,const char *const *names,size_t n){
    duckdb_logical_type fields[26];
    for(size_t i=0;i<n;++i) fields[i]=duckdb_create_logical_type(types[i]);
    duckdb_logical_type result=duckdb_create_struct_type(fields,(const char **)names,n);
    for(size_t i=0;i<n;++i) duckdb_destroy_logical_type(&fields[i]);
    return result;
}
static duckdb_logical_type interval_type(void){
    const duckdb_type types[]={DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_UBIGINT};
    const char *names[]={"start1","end1"};return record(types,names,2);
}
static duckdb_logical_type source_type(void){
    const duckdb_type types[]={DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_VARCHAR,DUCKDB_TYPE_VARCHAR};
    const char *names[]={"source_id","start1","reference","alternate"};return record(types,names,4);
}
static duckdb_logical_type effect_type(void){
    const duckdb_type types[]={DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_BOOLEAN,DUCKDB_TYPE_UINTEGER,
        DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_UINTEGER,DUCKDB_TYPE_UINTEGER};
    const char *names[]={"source_id","changed_dna","changed_rna_bases","genomic_start1",
        "genomic_end1","rna_first1","rna_last1"};return record(types,names,7);
}
bool duckvep_register_product_scenario(duckdb_connection connection){
    duckdb_logical_type ref,interval,source,exons,blocks,sources,programme,arg,result,effect,effects;
    const duckdb_type rt[]={DUCKDB_TYPE_UBIGINT,DUCKDB_TYPE_VARCHAR,DUCKDB_TYPE_TINYINT};
    const char *rn[]={"start1","dna","strand"};ref=record(rt,rn,3);
    interval=interval_type();exons=duckdb_create_list_type(interval);blocks=duckdb_create_list_type(interval);
    source=source_type();sources=duckdb_create_list_type(source);
    const duckdb_type pt[]={DUCKDB_TYPE_VARCHAR,DUCKDB_TYPE_VARCHAR,DUCKDB_TYPE_VARCHAR,
        DUCKDB_TYPE_VARCHAR,DUCKDB_TYPE_UINTEGER,DUCKDB_TYPE_UINTEGER,DUCKDB_TYPE_VARCHAR,DUCKDB_TYPE_UTINYINT};
    const char *pn[]={"id","evidence_id","assumption_id","expected_reference_rna",
        "annotated_start1","alternative_start1","alternative_trigger_id","genetic_code"};
    programme=record(pt,pn,8);
    duckdb_logical_type args[]={ref,exons,blocks,sources,programme};
    const char *an[]={"reference","exons","blocks","sources","programme"};
    arg=duckdb_create_struct_type(args,an,5);
    effect=effect_type();effects=duckdb_create_list_type(effect);
    duckdb_logical_type outputs[PRODUCT_OUTPUT_COUNT]={sources,effects};
    for(size_t i=2;i<PRODUCT_OUTPUT_COUNT;++i){
        duckdb_type type=!strcmp(duckvep_product_output_types[i],"VARCHAR")?DUCKDB_TYPE_VARCHAR:
            !strcmp(duckvep_product_output_types[i],"UBIGINT")?DUCKDB_TYPE_UBIGINT:
            !strcmp(duckvep_product_output_types[i],"TINYINT")?DUCKDB_TYPE_TINYINT:
            !strcmp(duckvep_product_output_types[i],"UTINYINT")?DUCKDB_TYPE_UTINYINT:DUCKDB_TYPE_UINTEGER;
        outputs[i]=duckdb_create_logical_type(type);
    }
    result=duckdb_create_struct_type(outputs,(const char **)duckvep_product_output_names,PRODUCT_OUTPUT_COUNT);
    duckdb_scalar_function fun=duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fun,"duckvep_product_scenario_v1");
    duckdb_scalar_function_add_parameter(fun,arg);
    duckdb_scalar_function_set_return_type(fun,result);
    duckdb_scalar_function_set_special_handling(fun);
    duckdb_scalar_function_set_function(fun,product_scalar);
    duckdb_state state=duckdb_register_scalar_function(connection,fun);
    duckdb_destroy_scalar_function(&fun);duckdb_destroy_logical_type(&result);
    for(size_t i=2;i<PRODUCT_OUTPUT_COUNT;++i) duckdb_destroy_logical_type(&outputs[i]);
    duckdb_destroy_logical_type(&effects);duckdb_destroy_logical_type(&effect);
    duckdb_destroy_logical_type(&arg);duckdb_destroy_logical_type(&programme);
    duckdb_destroy_logical_type(&sources);duckdb_destroy_logical_type(&source);
    duckdb_destroy_logical_type(&blocks);duckdb_destroy_logical_type(&exons);
    duckdb_destroy_logical_type(&interval);duckdb_destroy_logical_type(&ref);
    return state==DuckDBSuccess;
}
