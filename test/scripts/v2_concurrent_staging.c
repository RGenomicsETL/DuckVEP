/* Stress two exact stage keys from two connections to one database. */
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "duckdb.h"

static const char *regions_a =
    "COPY (SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER)) t(seq_region) ORDER BY seq_region) "
    "TO 'x' (FORMAT duckvep_stage, MODEL 'stage-race', LOAD_KEY 'A', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)";
static const char *regions_b =
    "COPY (SELECT * FROM (VALUES (0::UINTEGER), (1::UINTEGER), (2::UINTEGER)) t(seq_region) ORDER BY seq_region) "
    "TO 'x' (FORMAT duckvep_stage, MODEL 'stage-race', LOAD_KEY 'B', RELATION 'regions', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)";
static const char *transcripts =
    "COPY (SELECT 0::UINTEGER transcript_index, 1::UINTEGER seq_region, "
    "100::UBIGINT transcript_start, 250::UBIGINT transcript_end, 1::TINYINT strand, "
    "0::UINTEGER gene_index, 3::UBIGINT transcript_flags, 120::UBIGINT cds_start, 240::UBIGINT cds_end, "
    "'ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA'::BLOB cds_sequence, "
    "1::UTINYINT codon_table, 'TACGTACGTACGTACGTACG'::BLOB pre_cds_sequence, 'ACGTACGTAC'::BLOB post_cds_sequence) "
    "TO 'x' (FORMAT duckvep_stage, MODEL 'stage-race', LOAD_KEY '%s', RELATION 'transcripts', USE_TMP_FILE FALSE)";
static const char *exons =
    "COPY (SELECT * FROM (VALUES "
    "(0::UINTEGER, 100::UBIGINT, 150::UBIGINT, 1::UBIGINT, 51::UBIGINT, 0::TINYINT, 0::TINYINT), "
    "(0::UINTEGER, 200::UBIGINT, 250::UBIGINT, 52::UBIGINT, 102::UBIGINT, 0::TINYINT, 0::TINYINT)) "
    "t(transcript_index, exon_start, exon_end, exon_cdna_start, exon_cdna_end, phase, end_phase) "
    "ORDER BY transcript_index, exon_start) TO 'x' (FORMAT duckvep_stage, MODEL 'stage-race', "
    "LOAD_KEY '%s', RELATION 'exons', USE_TMP_FILE FALSE, PRESERVE_ORDER TRUE)";

static int query(duckdb_connection con, const char *sql, const char *expected_error) {
    duckdb_result result;
    duckdb_state status = duckdb_query(con, sql, &result);
    const char *error = duckdb_result_error(&result);
    int ok = expected_error ? status == DuckDBError && error && strstr(error, expected_error)
                            : status == DuckDBSuccess;
    if (!ok) {
        fprintf(stderr, "unexpected result for %s: %s\n", sql, error ? error : "success");
    }
    duckdb_destroy_result(&result);
    return ok;
}

static int query_keyed(duckdb_connection con, const char *sql_format, const char *load_key) {
    char sql[4096];
    int length = snprintf(sql, sizeof sql, sql_format, load_key);
    if (length < 0 || (size_t)length >= sizeof sql) {
        fputs("keyed staging SQL is too long\n", stderr);
        return 0;
    }
    return query(con, sql, NULL);
}

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t condition;
    int ready;
    int start;
} gate;

typedef struct {
    duckdb_connection connection;
    gate *start;
    const char *regions;
    int iterations;
    int ok;
} worker;

static void *stage_regions(void *arg) {
    worker *w = arg;
    pthread_mutex_lock(&w->start->mutex);
    w->start->ready++;
    pthread_cond_broadcast(&w->start->condition);
    while (!w->start->start) {
        pthread_cond_wait(&w->start->condition, &w->start->mutex);
    }
    pthread_mutex_unlock(&w->start->mutex);
    w->ok = 1;
    for (int i = 0; i < w->iterations && w->ok; i++) {
        w->ok = query(w->connection, w->regions, NULL);
    }
    return NULL;
}

static int stress(duckdb_connection connections[2], const char *first, const char *second, int iterations) {
    gate start = {PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, 0, 0};
    worker workers[2] = {{connections[0], &start, first, iterations, 0},
                         {connections[1], &start, second, iterations, 0}};
    pthread_t threads[2];
    int ok = 0;
    if (pthread_create(&threads[0], NULL, stage_regions, &workers[0]) != 0) goto cleanup;
    if (pthread_create(&threads[1], NULL, stage_regions, &workers[1]) != 0) {
        pthread_mutex_lock(&start.mutex);
        start.start = 1;
        pthread_cond_broadcast(&start.condition);
        pthread_mutex_unlock(&start.mutex);
        pthread_join(threads[0], NULL);
        goto cleanup;
    }
    pthread_mutex_lock(&start.mutex);
    while (start.ready < 2) pthread_cond_wait(&start.condition, &start.mutex);
    start.start = 1;
    pthread_cond_broadcast(&start.condition);
    pthread_mutex_unlock(&start.mutex);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    ok = workers[0].ok && workers[1].ok;
cleanup:
    pthread_cond_destroy(&start.condition);
    pthread_mutex_destroy(&start.mutex);
    return ok;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s /path/to/duckvep.duckdb_extension\n", argv[0]);
        return 2;
    }
    duckdb_database db = NULL;
    duckdb_connection main_connection = NULL, worker_connections[2] = {NULL, NULL};
    duckdb_config config = NULL;
    char *error = NULL;
    int ok = 0;
    if (duckdb_create_config(&config) == DuckDBError ||
        duckdb_set_config(config, "allow_unsigned_extensions", "true") == DuckDBError ||
        duckdb_open_ext(NULL, &db, config, &error) == DuckDBError ||
        duckdb_connect(db, &main_connection) == DuckDBError ||
        duckdb_connect(db, &worker_connections[0]) == DuckDBError ||
        duckdb_connect(db, &worker_connections[1]) == DuckDBError) {
        fprintf(stderr, "could not open DuckDB: %s\n", error ? error : "unknown error");
        goto cleanup;
    }
    /* Escape the extension path as a SQL string. */
    size_t path_size = strlen(argv[1]);
    char *load = malloc(2 * path_size + 16);
    if (!load) goto cleanup;
    char *p = load;
    memcpy(p, "LOAD '", 6); p += 6;
    for (const char *s = argv[1]; *s; s++) {
        *p++ = *s;
        if (*s == '\'') *p++ = '\'';
    }
    *p++ = '\''; *p = '\0';
    int loaded = query(main_connection, load, NULL);
    free(load);
    if (!loaded) goto cleanup;
    if (!query(main_connection, regions_b, NULL) ||
        !query_keyed(main_connection, transcripts, "B") || !query_keyed(main_connection, exons, "B") ||
        !query(main_connection, "SELECT duckvep_model_publish('stage-race', {'load_key': 'B'})", NULL) ||
        !query(main_connection, "SET VARIABLE expected_fingerprint = "
               "(SELECT CAST(_duckvep_model_fingerprint('stage-race') AS VARCHAR))", NULL) ||
        !query(main_connection, "SELECT duckvep_model_drop('stage-race')", NULL)) goto cleanup;

    if (!stress(worker_connections, regions_a, regions_b, 100)) goto cleanup;

    /* A failed A publish consumes only A; B remains complete and publishes its distinct rows. */
    ok = query_keyed(main_connection, transcripts, "A") && query_keyed(main_connection, exons, "A") &&
         query(main_connection, "SELECT duckvep_model_publish('stage-race', {'load_key': 'A', "
               "'reference_fasta': '/dev/null/duckvep-missing.fa'})", "reference_fasta requires") &&
         query_keyed(main_connection, transcripts, "B") && query_keyed(main_connection, exons, "B") &&
         query(main_connection, "SELECT duckvep_model_publish('stage-race', {'load_key': 'B'})", NULL) &&
         query(main_connection, "SELECT CASE WHEN CAST(_duckvep_model_fingerprint('stage-race') AS VARCHAR) = "
               "getvariable('expected_fingerprint') THEN true ELSE error('concurrent model rows differ') END", NULL) &&
         query(main_connection, "SELECT duckvep_model_drop('stage-race')", NULL);
    if (!ok) goto cleanup;
    puts("concurrent exact-key staging: B survives failed A publish with the serial fingerprint");
    const char *unkeyed = "COPY (SELECT 0::UINTEGER seq_region) TO 'x' "
        "(FORMAT duckvep_stage, MODEL 'atomic-race', RELATION 'regions', USE_TMP_FILE FALSE)";
    ok = stress(worker_connections, unkeyed, unkeyed, 1000) &&
        query(main_connection, "SELECT duckvep_model_publish('atomic-race')", "transcripts") &&
        query(main_connection, "SELECT duckvep_model_publish('atomic-race')", "regions");
    if (ok) puts("concurrent unkeyed replacement: one stage consumed by the first failed publish");
cleanup:
    if (error) duckdb_free(error);
    if (worker_connections[0]) duckdb_disconnect(&worker_connections[0]);
    if (worker_connections[1]) duckdb_disconnect(&worker_connections[1]);
    if (main_connection) duckdb_disconnect(&main_connection);
    if (db) duckdb_close(&db);
    if (config) duckdb_destroy_config(&config);
    return ok ? 0 : 1;
}
