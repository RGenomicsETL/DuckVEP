/* Model reloads between streamed annotation chunks must not mix generations. */
#include "duckdb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int run(duckdb_connection con, const char *sql) {
    duckdb_result result;
    int ok = duckdb_query(con, sql, &result) == DuckDBSuccess;
    if (!ok) fprintf(stderr, "%s\n", duckdb_result_error(&result));
    duckdb_destroy_result(&result);
    return ok;
}

static int load_model(duckdb_connection con, int v2) {
    if (!v2)
        return run(con, "SELECT loaded FROM duckvep_model_load('generation',"
                        "'SELECT * FROM regions','SELECT * FROM transcripts','SELECT * FROM exons')");
    duckdb_result script;
    if (duckdb_query(con, "SELECT unnest(duckvep_model_load_sql('generation',"
                         "'SELECT * FROM regions','SELECT * FROM transcripts','SELECT * FROM exons'))",
                     &script) != DuckDBSuccess) {
        fprintf(stderr, "%s\n", duckdb_result_error(&script));
        duckdb_destroy_result(&script);
        return 0;
    }
    int ok = 1;
    for (idx_t row = 0; ok && row < duckdb_row_count(&script); ++row) {
        char *statement = duckdb_value_varchar(&script, 0, row);
        ok = statement != NULL && run(con, statement);
        duckdb_free(statement);
    }
    duckdb_destroy_result(&script);
    return ok;
}

static int stream_case(duckdb_connection con, duckdb_connection control, int v2, int projected, int reload) {
    duckdb_prepared_statement prepared = NULL;
    duckdb_result result, dropped;
    const char *sql = projected ?
        "SELECT * FROM query(duckvep_annotate_projected_sql('events','generation'))" :
        "SELECT * FROM query(duckvep_annotate_sql('events','generation'))";
    if (!load_model(control, v2)) return 0;
    if (duckdb_prepare(con, sql, &prepared) != DuckDBSuccess) {
        fprintf(stderr, "Annotation prepare failed: %s\n", duckdb_prepare_error(prepared));
        return 0;
    }
    if (duckdb_execute_prepared_streaming(prepared, &result) != DuckDBSuccess ||
        !duckdb_result_is_streaming(result)) {
        fprintf(stderr, "Annotation did not produce a streaming result: %s\n", duckdb_result_error(&result));
        return 0;
    }
    duckdb_data_chunk chunk = duckdb_fetch_chunk(result);
    if (!chunk || !duckdb_data_chunk_get_size(chunk)) return 0;
    idx_t rows = duckdb_data_chunk_get_size(chunk);
    duckdb_destroy_data_chunk(&chunk);
    if (duckdb_query(control, "SELECT duckvep_model_drop('generation')", &dropped) != DuckDBSuccess ||
        !duckdb_value_boolean(&dropped, 0, 0)) {
        fprintf(stderr, "Model was not droppable between chunks\n");
        return 0;
    }
    duckdb_destroy_result(&dropped);
    if (reload && !load_model(control, v2)) return 0;
    while ((chunk = duckdb_fetch_chunk(result)) != NULL) {
        rows += duckdb_data_chunk_get_size(chunk);
        duckdb_destroy_data_chunk(&chunk);
    }
    const char *error = duckdb_result_error(&result);
    int ok = rows < 1000000 && error &&
        (strstr(error, "dropped or reloaded") || (!reload && strstr(error, "not loaded")));
    printf("projected=%d reload=%d consumed=%llu error=%s\n", projected, reload,
           (unsigned long long)rows, error ? error : "NONE");
    duckdb_destroy_result(&result);
    duckdb_destroy_prepare(&prepared);
    if (reload) ok = run(control, "SELECT duckvep_model_drop('generation')") && ok;
    return ok;
}

int main(int argc, char **argv) {
    if (argc != 2 || strchr(argv[1], '\'')) return 2;
    duckdb_config config = NULL;
    duckdb_database db = NULL;
    duckdb_connection con = NULL, control = NULL;
    char *error = NULL;
    if (duckdb_create_config(&config) != DuckDBSuccess ||
        duckdb_set_config(config, "allow_unsigned_extensions", "true") != DuckDBSuccess ||
        duckdb_set_config(config, "threads", "1") != DuckDBSuccess ||
        duckdb_open_ext(NULL, &db, config, &error) != DuckDBSuccess ||
        duckdb_connect(db, &con) != DuckDBSuccess || duckdb_connect(db, &control) != DuckDBSuccess)
        return 3;
    char sql[4096];
    if (snprintf(sql, sizeof(sql), "LOAD '%s'", argv[1]) >= (int)sizeof(sql) || !run(con, sql)) return 4;
    if (!run(con,
        "CREATE TABLE regions AS SELECT 1::UINTEGER seq_region;"
        "CREATE TABLE transcripts AS SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region,"
        "100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,"
        "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,111::UBIGINT cds_end,"
        "'ATGAAAGCCTAA'::BLOB cds_sequence,1::UTINYINT codon_table,"
        "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence;"
        "CREATE TABLE exons AS SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,"
        "111::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,"
        "0::TINYINT phase,0::TINYINT end_phase;"
        "CREATE TABLE events AS SELECT i::UBIGINT event_index,1::UINTEGER seq_region,"
        "104::UBIGINT position,'A' reference,'C' alternate,NULL::UBIGINT end_position,"
        "NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,"
        "NULL::UBIGINT mate_position FROM range(1000000) t(i)")) return 5;
    duckdb_result functions;
    if (duckdb_query(con, "SELECT count(*) FROM duckdb_functions() WHERE function_name='duckvep_model_publish'",
                     &functions) != DuckDBSuccess) return 1;
    int v2 = duckdb_value_int64(&functions, 0, 0) > 0;
    duckdb_destroy_result(&functions);
    for (int projected = 0; projected <= 1; ++projected)
        for (int reload = 1; reload >= 0; --reload)
            if (!stream_case(con, control, v2, projected, reload)) return 1;
    duckdb_disconnect(&control);
    duckdb_disconnect(&con);
    duckdb_close(&db);
    duckdb_destroy_config(&config);
    duckdb_free(error);
    return 0;
}
