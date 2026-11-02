/* Linux subprocess probe: annotation must recover after host allocation failure. */
#define _POSIX_C_SOURCE 200809L
#include "duckdb.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/resource.h>
#include <unistd.h>

static int run(duckdb_connection c, const char *sql) {
    duckdb_result r;
    int ok = duckdb_query(c, sql, &r) == DuckDBSuccess;
    if (!ok) fprintf(stderr, "query error: %s\n", duckdb_result_error(&r));
    duckdb_destroy_result(&r);
    return ok;
}

int main(int argc, char **argv) {
    if (argc < 3) return 2;
    setvbuf(stdout, NULL, _IONBF, 0);
    duckdb_config cfg = NULL;
    duckdb_database db = NULL;
    duckdb_connection con = NULL;
    duckdb_create_config(&cfg);
    duckdb_set_config(cfg, "allow_unsigned_extensions", "true");
    duckdb_set_config(cfg, "threads", "1");
    char *err = NULL;
    if (duckdb_open_ext(NULL, &db, cfg, &err) == DuckDBError || duckdb_connect(db, &con) == DuckDBError) return 3;
    char sql[4096];
    snprintf(sql, sizeof(sql), "LOAD '%s'", argv[1]);
    if (!run(con, sql)) return 4;
    int n = argc > 3 ? atoi(argv[3]) : 10000;
    if (!run(con, "SELECT duckvep_worker_limits_set(1,134217728,134217728,0); CREATE TABLE regions AS SELECT 1::UINTEGER seq_region;")) return 5;
    snprintf(sql, sizeof(sql), "CREATE TABLE transcripts AS SELECT i::UINTEGER transcript_index, 1::UINTEGER seq_region, 100::UBIGINT transcript_start, 111::UBIGINT transcript_end, 1::TINYINT strand, i::UINTEGER gene_index, 3::UBIGINT transcript_flags, 100::UBIGINT cds_start, 111::UBIGINT cds_end, 'ATGAAAGCCTAA'::BLOB cds_sequence, 1::UTINYINT codon_table, ''::BLOB pre_cds_sequence, ''::BLOB post_cds_sequence FROM range(%d) r(i)", n);
    if (!run(con, sql) || !run(con, "CREATE TABLE exons AS SELECT transcript_index, 100::UBIGINT exon_start, 111::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, 12::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase FROM transcripts; SELECT loaded FROM duckvep_model_load('m','SELECT * FROM regions','SELECT * FROM transcripts','SELECT * FROM exons'); CREATE TABLE events AS SELECT 1::UBIGINT event_index, 1::UINTEGER seq_region, 104::UBIGINT position, 'A' reference, 'C' alternate, NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change, NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position;")) return 6;
    duckdb_prepared_statement prepared = NULL;
    const char *annotation = "SELECT count(*) FROM query(duckvep_annotate_sql('events','m',{'upstream_distance':0,'downstream_distance':0}))";
    if (duckdb_prepare(con, annotation, &prepared) == DuckDBError) {
        fprintf(stderr,"prepare error %s\n",duckdb_prepare_error(prepared));return 7;
    }
    struct rlimit original, limited;
    getrlimit(RLIMIT_AS, &original);
    FILE *f = fopen("/proc/self/statm", "r");
    unsigned long pages = 0;
    if (!f || fscanf(f, "%lu", &pages) != 1) return 8;
    fclose(f);
    long slack = atol(argv[2]);
    limited = original;
    limited.rlim_cur = pages * (unsigned long)sysconf(_SC_PAGESIZE) + (rlim_t)slack * 1024;
    printf("pid=%d n=%d slack_kb=%ld vm_kb=%lu\n",getpid(),n,slack,pages * (unsigned long)sysconf(_SC_PAGESIZE)/1024);
    if (slack >= 0 && setrlimit(RLIMIT_AS, &limited) != 0) return 9;
    duckdb_result result;
    duckdb_state state = duckdb_execute_prepared(prepared, &result);
    if (setrlimit(RLIMIT_AS, &original) != 0) return 10;
    printf("address-space limit restored\n");
    printf("annotation=%s %s\n",state == DuckDBSuccess ? "ok" : "error", state == DuckDBSuccess ? "" : duckdb_result_error(&result));
    if (state == DuckDBSuccess) printf("rows=%lld\n", (long long)duckdb_value_int64(&result,0,0));
    duckdb_destroy_result(&result);
    if (argc > 4) {
        printf("recovery-start\n");
        state = duckdb_execute_prepared(prepared, &result);
        printf("recovery=%s\n", state == DuckDBSuccess ? "ok" : duckdb_result_error(&result));
        int recovered = state == DuckDBSuccess && duckdb_value_int64(&result, 0, 0) == n;
        duckdb_destroy_result(&result);
        if (!recovered) return 11;
    }
    duckdb_destroy_prepare(&prepared);
    printf("drop-start\n");
    state = duckdb_query(con, "SELECT duckvep_model_drop('m')", &result);
    int dropped = state == DuckDBSuccess && duckdb_value_boolean(&result, 0, 0);
    printf("drop=%s\n", dropped ? "true" : "false");
    duckdb_destroy_result(&result);
    duckdb_disconnect(&con);
    duckdb_close(&db);
    duckdb_destroy_config(&cfg);
    return dropped ? 0 : 12;
}
