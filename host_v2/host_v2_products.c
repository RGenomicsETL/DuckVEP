#include "host_v2_columns.h"
#include "duckvep_host.h"
#include "core/duckvep_core_products.h"
#include <string.h>

static bool refresh(column *c, duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail=NULL;
    DUCKDB_V2_ERROR status=duckdb_v2_vector_get_view(c->vector,&c->view,&detail);
    if(status!=DUCKDB_V2_ERROR_NONE){copy_duckdb_error(*error,status,detail);
        (void)duckdb_v2_error_info_destroy(&detail);return false;}
    (void)duckdb_v2_error_info_destroy(&detail);return true;
}
static bool open_fields(column *parent,column *fields,size_t count,
                        duckdb_v2_error_info_handle *error){
    for(size_t i=0;i<count;++i) if(!column_child(parent,i,&fields[i],error)) return false;
    if(!refresh(parent,error)) return false;
    for(size_t i=0;i<count;++i) if(!refresh(&fields[i],error)) return false;
    return true;
}
static bool span(const column *list,const column *elements,idx_t row,size_t *at,size_t *count){
    if(!column_valid(list,row))return false;
    duckdb_v2_list_entry entry=list_at(list,row);
    if(entry.offset>elements->view.count || entry.length>elements->view.count-entry.offset ||
       entry.length>SIZE_MAX)return false;
    *at=(size_t)entry.offset;*count=(size_t)entry.length;return true;
}
static bool cell(const column *c,idx_t at,duckvep_cell_t *value){
    fill_cell(c,at,value);return value->valid;
}
static duckvep_product_text_t as_text(duckvep_cell_t c){
    return (duckvep_product_text_t){c.text,c.text_length};
}
static void product_exec(duckdb_v2_scalar_function_exec_info_handle info,
                         duckdb_v2_context_handle context,duckdb_v2_error_info_handle *error){
    (void)context;
    duckdb_v2_error_info_handle detail=NULL;
    duckdb_v2_vector_handle arg=NULL,out=NULL;
    idx_t rows=0;
    column root,top[5],reference[3],programme[8],exon_item,block_item,source_item;
    column exons[2],blocks[2],sources[4];
    v2_call call={0};call.error=error;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info,&rows,&detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info,0,&arg,&detail));
    if(!column_open(arg,true,&root,error) || !open_fields(&root,top,5,error))goto cleanup;
    if(!open_fields(&top[0],reference,3,error)||!open_fields(&top[4],programme,8,error))goto cleanup;
    if(!column_child(&top[1],0,&exon_item,error)||!open_fields(&exon_item,exons,2,error) ||
       !column_child(&top[2],0,&block_item,error)||!open_fields(&block_item,blocks,2,error) ||
       !column_child(&top[3],0,&source_item,error)||!open_fields(&source_item,sources,4,error))goto cleanup;
    if(!refresh(&root,error))goto cleanup;
    for(size_t i=0;i<5;++i)if(!refresh(&top[i],error))goto cleanup;
    for(size_t i=0;i<3;++i)if(!refresh(&reference[i],error))goto cleanup;
    for(size_t i=0;i<8;++i)if(!refresh(&programme[i],error))goto cleanup;
    column *items[]={&exon_item,&block_item,&source_item};
    for(size_t i=0;i<3;++i)if(!refresh(items[i],error))goto cleanup;
    for(size_t i=0;i<2;++i)if(!refresh(&exons[i],error)||!refresh(&blocks[i],error))goto cleanup;
    for(size_t i=0;i<4;++i)if(!refresh(&sources[i],error))goto cleanup;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info,&out,&detail));
    if(!v2_open_writable(&call,&call.output,out,rows,false))goto cleanup;
    v2_vec *result[PRODUCT_OUTPUT_COUNT];
    for(size_t i=0;i<PRODUCT_OUTPUT_COUNT;++i)
        result[i]=v2_child_of(&call.output,i,rows,true);
    if(call.failed)goto cleanup;
    for(idx_t row=0;row<rows;++row){
        duckvep_core_product_fields_t f={0};
        duckvep_core_product_result_t owned={0};
        duckvep_product_interval_t *xi=NULL,*bi=NULL;
        duckvep_product_source_t *si=NULL;
        duckvep_cell_t r[3],p[8],c[4];
        const char *message="invalid scenario shape or source";
        size_t at=0,n=0;
        if(!column_valid(&root,row)||!column_valid(&top[0],row)||!column_valid(&top[4],row))goto row_fail;
        for(size_t i=0;i<3;++i)if(!cell(&reference[i],row,&r[i]))goto row_fail;
        for(size_t i=0;i<8;++i)if(!cell(&programme[i],row,&p[i]))goto row_fail;
        f.reference_start1=r[0].u;f.dna=as_text(r[1]);f.strand=(int8_t)r[2].i;
        if(!span(&top[1],&exon_item,row,&at,&n)||!n||n>SIZE_MAX/sizeof *xi)goto row_fail;
        xi=duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,n*sizeof *xi);
        if(!xi){message="allocation limit";goto row_fail;}
        for(size_t i=0;i<n;++i){
            if(!column_valid(&exon_item,at+i)||!cell(&exons[0],at+i,&c[0])||
               !cell(&exons[1],at+i,&c[1]))goto row_fail;
            xi[i]=(duckvep_product_interval_t){c[0].u,c[1].u};
        }
        f.exons=xi;f.exon_count=n;
        if(!span(&top[2],&block_item,row,&at,&n)||!n||n>SIZE_MAX/sizeof *bi)goto row_fail;
        bi=duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,n*sizeof *bi);
        if(!bi){message="allocation limit";goto row_fail;}
        for(size_t i=0;i<n;++i){
            if(!column_valid(&block_item,at+i)||!cell(&blocks[0],at+i,&c[0])||
               !cell(&blocks[1],at+i,&c[1]))goto row_fail;
            bi[i]=(duckvep_product_interval_t){c[0].u,c[1].u};
        }
        f.blocks=bi;f.block_count=n;
        if(!span(&top[3],&source_item,row,&at,&n)||!n||n>SIZE_MAX/sizeof *si)goto row_fail;
        si=duckvep_budget_malloc(DUCKVEP_OWNER_WORKSPACE,n*sizeof *si);
        if(!si){message="allocation limit";goto row_fail;}
        for(size_t i=0;i<n;++i){
            if(!column_valid(&source_item,at+i))goto row_fail;
            for(size_t j=0;j<4;++j)if(!cell(&sources[j],at+i,&c[j]))goto row_fail;
            si[i]=(duckvep_product_source_t){c[0].u,c[1].u,as_text(c[2]),as_text(c[3])};
        }
        f.sources=si;f.source_count=n;
        f.id=as_text(p[0]);f.evidence_id=as_text(p[1]);f.assumption_id=as_text(p[2]);
        f.expected_rna=as_text(p[3]);f.annotated_start1=(uint32_t)p[4].u;
        f.alternative_start1=(uint32_t)p[5].u;f.trigger_id=as_text(p[6]);
        f.genetic_code=(uint8_t)p[7].u;
        duckvep_profile_status_t status=duckvep_core_product_scenario(&f,&owned);
        if(status!=DUCKVEP_PROFILE_OK){
            const char *messages[]={"ok","invalid scenario shape or source","genomic or RNA reference mismatch",
                "sequence unavailable","ambiguous composition","partial model or programme",
                "short workspace","allocation limit","unsupported product or initiation"};
            message=messages[status];goto row_fail;
        }
        for(size_t k=0;k<2;++k){
            duckdb_v2_list_entry entry;
            if(!duckvep_h_list_extend(result[k],n,&entry)){
                message="output list size overflow";goto row_fail;
            }
            ((duckdb_v2_list_entry *)result[k]->data)[row]=entry;
            v2_vec *item=duckvep_h_list_child(result[k]);
            size_t field_count=k?7:4;
            v2_vec *parts[7];
            for(size_t j=0;j<field_count;++j)parts[j]=duckvep_h_struct_child(item,j);
            if(call.failed){message="output list allocation failed";goto row_fail;}
            for(size_t i=0;i<n;++i){
                size_t pos=entry.offset+i;
                if(k){
                    duckvep_core_product_effect_t *e=&owned.effects[i];
                    ((uint64_t *)parts[0]->data)[pos]=e->source_id;
                    ((uint8_t *)parts[1]->data)[pos]=e->changed_dna!=0;
                    ((uint32_t *)parts[2]->data)[pos]=e->changed_rna_bases;
                    ((uint64_t *)parts[3]->data)[pos]=e->genomic_start1;
                    ((uint64_t *)parts[4]->data)[pos]=e->genomic_end1;
                    ((uint32_t *)parts[5]->data)[pos]=e->rna_first1;
                    ((uint32_t *)parts[6]->data)[pos]=e->rna_last1;
                }else{
                    ((uint64_t *)parts[0]->data)[pos]=si[i].source_id;
                    ((uint64_t *)parts[1]->data)[pos]=si[i].start1;
                    duckvep_h_assign_string(parts[2],pos,si[i].reference.ptr,si[i].reference.length);
                    duckvep_h_assign_string(parts[3],pos,si[i].alternate.ptr,si[i].alternate.length);
                }
            }
        }
        for(size_t i=2;i<PRODUCT_OUTPUT_COUNT;++i){
            const char *text;size_t length;uint64_t value;
            (void)duckvep_core_product_output(&f,&owned,i,&text,&length,&value);
            if(text)duckvep_h_assign_string(result[i],row,text,length);
            else if(!strcmp(duckvep_product_output_types[i],"UBIGINT"))
                ((uint64_t *)result[i]->data)[row]=value;
            else if(!strcmp(duckvep_product_output_types[i],"TINYINT"))
                ((int8_t *)result[i]->data)[row]=(int8_t)value;
            else if(!strcmp(duckvep_product_output_types[i],"UTINYINT"))
                ((uint8_t *)result[i]->data)[row]=(uint8_t)value;
            else ((uint32_t *)result[i]->data)[row]=(uint32_t)value;
        }
        duckvep_core_product_dispose(&owned);
        duckvep_budget_free(xi);duckvep_budget_free(bi);duckvep_budget_free(si);
        if(call.failed)goto cleanup;
        continue;
row_fail:
        duckvep_core_product_dispose(&owned);
        duckvep_budget_free(xi);duckvep_budget_free(bi);duckvep_budget_free(si);
        if(!call.failed)report(*error,message);
        goto cleanup;
    }
    v2_finish(&call);
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

bool host_v2_register_products(duckdb_v2_extension_handle extension,duckdb_v2_context_handle context,
                               duckdb_v2_error_info_handle *error){
    const char *const args[]={duckvep_product_input_type};
    const char *const parameters[]={"scenario"};
    return register_scalar(extension,context,"duckvep_product_scenario",args,parameters,
                           1,duckvep_product_output_type,product_exec,error);
}
