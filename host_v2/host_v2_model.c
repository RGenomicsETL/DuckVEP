/* The model sink of the v2 host: duckvep_model_load_sql (the statements to run),
 * the COPY format duckvep_stage (spillable staging of a caller's typed SELECT,
 * run in the caller's own transaction), duckvep_model_publish (build and install
 * from the staged rows), duckvep_model_drop and an internal fingerprint.
 *
 * v1 loads through private connections and stays as it is. v2 cannot run the
 * caller's queries on another connection, so the caller's COPY statements feed
 * DuckDB-managed column collections (they spill under the memory limit) and
 * publish scans them through the same src/core loaders the v1 host uses. A COPY
 * that fails or is cancelled destroys its rows; publish consumes the staging,
 * whether it succeeds or not. */
#include "host_v2_columns.h"
#include "host_v2_stage.h"

#include "core/duckvep_core_model.h"
#include "core/duckvep_core_snapshot.h"
#include "core/duckvep_core_model_script.h"

#include <pthread.h>


/* ---------------------------------------------------------------------------
 * State shared by the model functions of one database
 * ------------------------------------------------------------------------- */

void host_v2_stage_destroy(stage *s) {
    if (!s || __atomic_sub_fetch(&s->refs, 1, __ATOMIC_ACQ_REL) != 0) {
        return;
    }
    if (s->collection) {
        (void)duckdb_v2_column_data_collection_destroy(&s->collection);
    }
    for (size_t i = 0; i < s->columns; ++i) {
        free(s->names[i]);
        (void)duckdb_v2_logical_type_destroy(&s->types[i]);
    }
    free(s->model);
    free(s->relation);
    free(s->hap_model);
    free(s);
}

void host_v2_state_retain(model_state *state) {
    pthread_mutex_lock(&state->lock);
    state->references++;
    pthread_mutex_unlock(&state->lock);
}

void host_v2_state_release(void *pointer) {
    model_state *state = pointer;
    bool last;
    if (!state) {
        return;
    }
    pthread_mutex_lock(&state->lock);
    last = --state->references == 0;
    pthread_mutex_unlock(&state->lock);
    if (!last) {
        return;
    }
    while (state->stages) {
        stage *next = state->stages->next;
        host_v2_stage_destroy(state->stages);
        state->stages = next;
    }
    duckvep_registry_release(state->registry);
    pthread_mutex_destroy(&state->lock);
    free(state);
}

/* Removes and returns the stage of (model, relation), or NULL. */
stage *host_v2_stage_take(model_state *state, const char *model, const char *relation) {
    stage **link;
    pthread_mutex_lock(&state->lock);
    for (link = &state->stages; *link; link = &(*link)->next) {
        if (strcmp((*link)->model, model) == 0 && strcmp((*link)->relation, relation) == 0) {
            stage *found = *link;
            *link = found->next;
            found->next = NULL;
            pthread_mutex_unlock(&state->lock);
            return found;
        }
    }
    pthread_mutex_unlock(&state->lock);
    return NULL;
}

static char *copy_text(const char *text, size_t length);

stage *host_v2_stage_acquire(model_state *state, const char *key, const char *relation) {
    stage *found = NULL;
    pthread_mutex_lock(&state->lock);
    for (stage *s = state->stages; s; s = s->next) {
        if (strcmp(s->model, key) == 0 && strcmp(s->relation, relation) == 0) {
            __atomic_add_fetch(&s->refs, 1, __ATOMIC_ACQ_REL);
            found = s;
            break;
        }
    }
    pthread_mutex_unlock(&state->lock);
    return found;
}

bool host_v2_stage_peek_job(model_state *state, const char *job, const char *relation, char **model,
                            duckvep_hap_config_t *config) {
    bool found = false;
    pthread_mutex_lock(&state->lock);
    for (stage *s = state->stages; s; s = s->next) {
        if (strcmp(s->model, job) == 0 && strcmp(s->relation, relation) == 0 && s->hap_model) {
            *model = copy_text(s->hap_model, strlen(s->hap_model));
            *config = s->hap;
            found = *model != NULL;
            break;
        }
    }
    pthread_mutex_unlock(&state->lock);
    return found;
}

static void state_discard_model(model_state *state, const char *model) {
    stage *s;
    for (size_t i = 0; i < 6; ++i) {
        while ((s = host_v2_stage_take(state, model, duckvep_model_relations[i])) != NULL) {
            host_v2_stage_destroy(s);
        }
    }
}

static bool model_exists(model_state *state, const char *name) {
    bool found;
    pthread_mutex_lock(&state->registry->mutex);
    found = duckvep_registry_find_locked(state->registry, name) != NULL;
    pthread_mutex_unlock(&state->registry->mutex);
    return found;
}

static duckvep_ctype_t ctype_of(DUCKDB_V2_LOGICAL_TYPE_ID id) {
    switch (id) {
    case DUCKDB_V2_LOGICAL_TYPE_ID_BOOLEAN: return DUCKVEP_CT_BOOLEAN;
    case DUCKDB_V2_LOGICAL_TYPE_ID_TINYINT: return DUCKVEP_CT_TINYINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UTINYINT: return DUCKVEP_CT_UTINYINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UINTEGER: return DUCKVEP_CT_UINTEGER;
    case DUCKDB_V2_LOGICAL_TYPE_ID_UBIGINT: return DUCKVEP_CT_UBIGINT;
    case DUCKDB_V2_LOGICAL_TYPE_ID_VARCHAR: return DUCKVEP_CT_VARCHAR;
    case DUCKDB_V2_LOGICAL_TYPE_ID_BLOB: return DUCKVEP_CT_BLOB;
    default: return DUCKVEP_CT_OTHER;
    }
}

static char *copy_text(const char *text, size_t length) {
    char *copy = malloc(length + 1);
    if (copy) {
        memcpy(copy, text, length);
        copy[length] = '\0';
    }
    return copy;
}

/* ---------------------------------------------------------------------------
 * COPY ... TO ... (FORMAT duckvep_stage, MODEL '<name>', RELATION '<relation>')
 * ------------------------------------------------------------------------- */

typedef struct {
    model_state *state;
    stage *layout; /* the schema, with no collection */
} stage_bind;

typedef struct {
    stage_bind *bind;
    duckdb_v2_column_data_collection_handle collection;
} stage_init;

typedef struct {
    duckdb_v2_column_data_collection_handle collection;
} stage_batch;

static void stage_bind_destroy(void *pointer) {
    stage_bind *bind = pointer;
    if (bind) {
        host_v2_stage_destroy(bind->layout);
        free(bind);
    }
}

static void stage_init_destroy(void *pointer) {
    stage_init *init = pointer;
    if (init) {
        /* A COPY that did not finalize (failed or cancelled) drops its rows here. */
        if (init->collection) {
            (void)duckdb_v2_column_data_collection_destroy(&init->collection);
        }
        free(init);
    }
}

static void stage_batch_destroy(void *pointer) {
    stage_batch *batch = pointer;
    if (batch) {
        if (batch->collection) {
            (void)duckdb_v2_column_data_collection_destroy(&batch->collection);
        }
        free(batch);
    }
}

static bool identifier_is(const duckdb_v2_identifier_t *name, const char *text) {
    size_t length = strlen(text);
    if (name->len != length) {
        return false;
    }
    for (size_t i = 0; i < length; ++i) {
        char c = name->ptr[i];
        if (c >= 'A' && c <= 'Z') {
            c = (char)(c + ('a' - 'A'));
        }
        if (c != text[i]) {
            return false;
        }
    }
    return true;
}

static void stage_bind_exec(duckdb_v2_copy_to_bind_info_handle info, duckdb_v2_context_handle context,
                            duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    stage_bind *bind = NULL;
    duckdb_v2_value_handle value = NULL;
    idx_t options = 0, columns = 0;
    char *model = NULL, *relation = NULL;
    hap_bind *hap = NULL;
    bool owned = false, job = false;
    DUCKDB_CALL(duckdb_v2_copy_to_bind_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_copy_to_bind_get_option_count(info, &options, &detail));
    hap = host_v2_hap_bind_create();
    if (!hap) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        goto cleanup;
    }
    for (idx_t i = 0; i < options; ++i) {
        duckdb_v2_identifier_t name;
        duckdb_v2_str text;
        bool is_model, is_relation, consumed = false;
        DUCKDB_CALL(duckdb_v2_copy_to_bind_get_option_name(info, i, &name, &detail));
        is_model = identifier_is(&name, "model");
        is_relation = identifier_is(&name, "relation");
        DUCKDB_CALL(duckdb_v2_copy_to_bind_get_option_value(info, i, &value, &detail));
        if (!is_model && !is_relation) {
            /* The options of a haplotype job (JOB, STAGE, PHASE_POLICY, ...); anything else is DuckDB's own. */
            if (!host_v2_hap_bind_option(hap, &name, value, &consumed, error)) {
                goto cleanup;
            }
            DUCKDB_CALL(duckdb_v2_value_destroy(&value));
            continue;
        }
        DUCKDB_CALL(duckdb_v2_value_get_varchar(value, &text, &detail));
        if (is_model) {
            free(model);
            model = copy_text(text.ptr, text.len);
        } else {
            free(relation);
            relation = copy_text(text.ptr, text.len);
        }
        DUCKDB_CALL(duckdb_v2_value_destroy(&value));
    }
    job = host_v2_hap_bind_is_job(hap);
    if (!job) {
        if (!model || !*model || !relation) {
            INPUT_ERROR("duckvep_stage: MODEL and RELATION options are required");
        }
        {
            bool known = false;
            for (size_t i = 0; i < 6; ++i) {
                known = known || strcmp(relation, duckvep_model_relations[i]) == 0;
            }
            if (!known) {
                INPUT_ERROR("duckvep_stage: RELATION must be regions, transcripts, exons, mature_mirna, "
                            "peptide_edits or interval_features");
            }
        }
        if (model_exists(state, model)) {
            INPUT_ERROR("duckvep_stage: model name already exists");
        }
    }
    DUCKDB_CALL(duckdb_v2_copy_to_bind_get_column_count(info, &columns, &detail));
    if (columns == 0 || columns > MAX_STAGE_COLUMNS) {
        INPUT_ERROR("duckvep_stage: the query must return between 1 and 20 columns");
    }
    bind = calloc(1, sizeof(*bind));
    if (bind) {
        bind->layout = calloc(1, sizeof(*bind->layout));
        if (bind->layout) {
            bind->layout->refs = 1;
        }
    }
    if (!bind || !bind->layout) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        goto cleanup;
    }
    bind->state = state;
    if (!job) {
        bind->layout->model = model;
        bind->layout->relation = relation;
        model = relation = NULL;
    }
    for (idx_t i = 0; i < columns; ++i) {
        duckdb_v2_identifier_t name;
        duckdb_v2_logical_type_handle type = NULL;
        DUCKDB_V2_LOGICAL_TYPE_ID id;
        DUCKDB_CALL(duckdb_v2_copy_to_bind_get_column_name(info, i, &name, &detail));
        DUCKDB_CALL(duckdb_v2_copy_to_bind_get_column_type(info, i, &type, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_get_id(type, &id, &detail));
        bind->layout->names[i] = copy_text(name.ptr, name.len);
        bind->layout->types[i] = type;
        bind->layout->ctypes[i] = ctype_of(id);
        bind->layout->columns = i + 1;
        if (!bind->layout->names[i]) {
            set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
            goto cleanup;
        }
    }
    if (job && !host_v2_hap_bind_finish(hap, state, context, bind->layout, model, error)) {
        goto cleanup;
    }
    {
        duckdb_v2_opaque data = {bind, stage_bind_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_copy_to_bind_set_bind_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        stage_bind_destroy(bind);
    }
    host_v2_hap_bind_destroy(hap);
    free(model);
    free(relation);
    (void)duckdb_v2_value_destroy(&value);
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void stage_init_exec(duckdb_v2_copy_to_init_info_handle info, duckdb_v2_context_handle context,
                            duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    stage_bind *bind = NULL;
    stage_init *init = calloc(1, sizeof(*init));
    bool owned = false;
    if (!init) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_copy_to_init_get_bind_data(info, (void **)&bind, &detail));
    init->bind = bind;
    {
        duckdb_v2_opaque data = {init, stage_init_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_copy_to_init_set_init_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        stage_init_destroy(init);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void stage_batch_exec(duckdb_v2_copy_to_batch_info_handle info, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    stage_batch *batch = calloc(1, sizeof(*batch));
    bool owned = false;
    if (!batch) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        return;
    }
    DUCKDB_CALL(duckdb_v2_copy_to_batch_take_input(info, &batch->collection, &detail));
    {
        duckdb_v2_opaque data = {batch, stage_batch_destroy, NULL};
        DUCKDB_CALL(duckdb_v2_copy_to_batch_set_batch_data(info, &data, &detail));
        owned = true;
    }
cleanup:
    if (!owned) {
        stage_batch_destroy(batch);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* Ordered: each batch is merged into the COPY's own collection, so the row
 * order the loaders require is the query's order. */
static void stage_flush_exec(duckdb_v2_copy_to_flush_info_handle info, duckdb_v2_context_handle context,
                             duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    stage_init *init = NULL;
    stage_batch *batch = NULL;
    DUCKDB_CALL(duckdb_v2_copy_to_flush_get_init_data(info, (void **)&init, &detail));
    DUCKDB_CALL(duckdb_v2_copy_to_flush_get_batch_data(info, (void **)&batch, &detail));
    if (!init->collection) {
        stage *layout = init->bind->layout;
        DUCKDB_CALL(duckdb_v2_column_data_collection_create_with_context(
            context, layout->types, layout->columns, &init->collection, &detail));
    }
    DUCKDB_CALL(duckdb_v2_column_data_collection_combine(init->collection, &batch->collection, &detail));
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* Only a COPY that ran to its end publishes its rows to the staging area. */
static void stage_finalize_exec(duckdb_v2_copy_to_finalize_info_handle info,
                                duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    stage_init *init = NULL;
    stage *entry = NULL, *previous;
    stage_bind *bind;
    DUCKDB_CALL(duckdb_v2_copy_to_finalize_get_init_data(info, (void **)&init, &detail));
    bind = init->bind;
    entry = calloc(1, sizeof(*entry));
    if (entry) {
        entry->refs = 1;
    }
    if (!entry) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        goto cleanup;
    }
    entry->model = copy_text(bind->layout->model, strlen(bind->layout->model));
    entry->relation = copy_text(bind->layout->relation, strlen(bind->layout->relation));
    if (bind->layout->hap_model) {
        entry->hap_model = copy_text(bind->layout->hap_model, strlen(bind->layout->hap_model));
        entry->hap = bind->layout->hap;
    }
    if (!entry->model || !entry->relation || (bind->layout->hap_model && !entry->hap_model)) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
        goto cleanup;
    }
    for (size_t i = 0; i < bind->layout->columns; ++i) {
        entry->names[i] = copy_text(bind->layout->names[i], strlen(bind->layout->names[i]));
        entry->ctypes[i] = bind->layout->ctypes[i];
        entry->columns = i + 1;
        if (!entry->names[i]) {
            set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_stage: out of memory");
            goto cleanup;
        }
        DUCKDB_CALL(duckdb_v2_logical_type_copy(bind->layout->types[i], &entry->types[i], &detail));
    }
    if (init->collection) {
        DUCKDB_CALL(duckdb_v2_column_data_collection_row_count(init->collection, &entry->rows, &detail));
    }
    entry->collection = init->collection;
    init->collection = NULL;
    /* Staging the same relation again replaces the earlier rows. */
    previous = host_v2_stage_take(bind->state, entry->model, entry->relation);
    host_v2_stage_destroy(previous);
    pthread_mutex_lock(&bind->state->lock);
    entry->next = bind->state->stages;
    bind->state->stages = entry;
    pthread_mutex_unlock(&bind->state->lock);
    entry = NULL;
cleanup:
    host_v2_stage_destroy(entry);
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * Row sources over staged rows
 * ------------------------------------------------------------------------- */

typedef struct {
    stage *staged;
    duckdb_v2_context_handle context;
    duckdb_v2_error_info_handle *error;
    duckdb_v2_column_data_collection_shared_scan_state_handle shared;
    duckdb_v2_column_data_collection_worker_scan_state_handle worker;
    duckdb_v2_data_chunk_handle chunk;
    duckvep_batch_t batch;
    duckvep_col_t columns[MAX_STAGE_COLUMNS];
} staged_source;

static int ss_count(void *context, size_t *rows, char *error, size_t error_size) {
    staged_source *source = context;
    (void)error;
    (void)error_size;
    *rows = (size_t)source->staged->rows;
    return 1;
}

static void ss_close(void *context) {
    staged_source *source = context;
    (void)duckdb_v2_column_data_collection_worker_scan_state_destroy(&source->worker);
    (void)duckdb_v2_column_data_collection_shared_scan_state_destroy(&source->shared);
    (void)duckdb_v2_data_chunk_destroy(&source->chunk);
}

static int ss_open(void *context, char *error, size_t error_size) {
    staged_source *source = context;
    duckdb_v2_error_info_handle detail = NULL;
    stage *s = source->staged;
    int ok = 0;
    if (!s->collection) {
        return 1;
    }
    if (duckdb_v2_column_data_collection_shared_scan_state_create(s->collection, &source->shared, &detail) ||
        duckdb_v2_column_data_collection_worker_scan_state_create(s->collection, &source->worker, &detail) ||
        duckdb_v2_data_chunk_create_with_context(source->context, s->types, s->columns, &source->chunk, &detail)) {
        duckdb_v2_str text = {"could not scan the staged rows", 30};
        if (detail) {
            (void)duckdb_v2_error_info_get_text(detail, &text);
        }
        (void)snprintf(error, error_size, "%.*s", (int)text.len, text.ptr);
        ss_close(source);
        goto done;
    }
    ok = 1;
done:
    (void)duckdb_v2_error_info_destroy(&detail);
    return ok;
}

static size_t ss_columns(void *context) {
    return ((staged_source *)context)->staged->columns;
}

static const char *ss_name(void *context, size_t column) {
    return ((staged_source *)context)->staged->names[column];
}

static duckvep_ctype_t ss_type(void *context, size_t column) {
    return ((staged_source *)context)->staged->ctypes[column];
}

static void ss_release(duckvep_batch_t *batch) {
    (void)batch; /* the chunk is reused by the next scan */
}

static duckvep_batch_t *ss_next(void *context) {
    staged_source *source = context;
    duckdb_v2_error_info_handle detail = NULL;
    stage *s = source->staged;
    bool produced = false;
    idx_t rows = 0;
    duckvep_batch_t *result = NULL;
    if (!s->collection || !source->chunk) {
        return NULL;
    }
    if (duckdb_v2_column_data_collection_scan(s->collection, source->shared, source->worker, source->chunk,
                                              &produced, &detail) ||
        !produced || duckdb_v2_data_chunk_get_size(source->chunk, &rows, &detail)) {
        goto done;
    }
    for (size_t column = 0; column < s->columns; ++column) {
        duckdb_v2_vector_handle vector = NULL;
        duckdb_v2_vector_view view;
        if (duckdb_v2_data_chunk_get_vector(source->chunk, column, &vector, &detail) ||
            duckdb_v2_vector_flatten(vector, &detail) ||
            duckdb_v2_vector_get_view(vector, &view, &detail)) {
            goto done;
        }
        source->columns[column].data = (void *)view.data;
        source->columns[column].validity = view.validity;
    }
    source->batch.rows = rows;
    source->batch.columns = source->columns;
    source->batch.context = source;
    source->batch.release = ss_release;
    result = &source->batch;
done:
    (void)duckdb_v2_error_info_destroy(&detail);
    return result;
}

static void staged_source_init(staged_source *source, duckvep_source_t *out, stage *staged,
                               duckdb_v2_context_handle context) {
    memset(source, 0, sizeof(*source));
    source->staged = staged;
    source->context = context;
    out->context = source;
    out->count = ss_count;
    out->open = ss_open;
    out->column_count = ss_columns;
    out->column_name = ss_name;
    out->column_type = ss_type;
    out->next = ss_next;
    out->close = ss_close;
}

/* ---------------------------------------------------------------------------
 * duckvep_model_publish(name [, {reference_fasta, transcript_coverage_complete}])
 * ------------------------------------------------------------------------- */

static void publish_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                         duckdb_v2_error_info_handle *error) {
    static const char *const keys[] = {"reference_fasta", "transcript_coverage_complete"};
    static const duckvep_core_option_kind_t kinds[] = {DUCKVEP_CORE_OPTION_TEXT, DUCKVEP_CORE_OPTION_BOOLEAN};
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    duckdb_v2_vector_handle arguments[2] = {0};
    column name_column;
    options option;
    duckdb_v2_vector_handle output = NULL;
    void *output_data = NULL;
    uint64_t *output_validity = NULL;
    idx_t rows = 0;
    uint32_t argc = 0;
    stage *staged[6] = {0};
    char *name = NULL, *reference = NULL;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    for (uint32_t i = 0; i < argc; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
    }
    if (!column_open(arguments[0], false, &name_column, error) ||
        !options_open(argc == 2 ? arguments[1] : NULL, keys, kinds, 2, &option, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &output_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &output_validity, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_cell_t cell;
        int coverage = 0;
        staged_source sources[6];
        duckvep_source_t source_views[6];
        duckvep_model_sources_t model_sources;
        char message[DUCKVEP_SQL_ERROR_SIZE + 256];
        bool ok;
        name = column_valid(&name_column, row) ? NULL : NULL;
        if (column_valid(&name_column, row)) {
            duckdb_v2_str text = string_at(&name_column.view, row);
            name = duckvep_core_string_copy(text.ptr, text.len);
        }
        if (!name || !*name) {
            INPUT_ERROR("duckvep_model_publish: name must be a non-empty string");
        }
        if (!options_check(&option, row, error)) {
            goto cleanup;
        }
        if (option_cell(&option, 0, row, &cell) && cell.valid) {
            reference = duckvep_core_string_copy(cell.text, cell.text_length);
            if (!reference || !*reference) {
                INPUT_ERROR("duckvep_model_publish: reference_fasta must be a non-empty string");
            }
        }
        if (option_cell(&option, 1, row, &cell)) {
            if (!cell.valid) {
                INPUT_ERROR("duckvep_model_publish: transcript_coverage_complete cannot be NULL");
            }
            coverage = cell.boolean;
        }
        /* The staging is consumed by this call, whatever its outcome. */
        for (size_t i = 0; i < 6; ++i) {
            staged[i] = host_v2_stage_take(state, name, duckvep_model_relations[i]);
        }
        for (size_t i = 0; i < 3; ++i) {
            if (!staged[i]) {
                (void)snprintf(message, sizeof message,
                               "duckvep_model_publish: relation '%s' is not staged for this model "
                               "(run the statements of duckvep_model_load_sql first)",
                               duckvep_model_relations[i]);
                INPUT_ERROR(message);
            }
        }
        memset(&model_sources, 0, sizeof(model_sources));
        for (size_t i = 0; i < 6; ++i) {
            if (staged[i]) {
                staged_source_init(&sources[i], &source_views[i], staged[i], context);
            }
        }
        model_sources.regions = &source_views[0];
        model_sources.transcripts = &source_views[1];
        model_sources.exons = &source_views[2];
        model_sources.mature_mirna = staged[3] ? &source_views[3] : NULL;
        model_sources.peptide_edits = staged[4] ? &source_views[4] : NULL;
        model_sources.interval_features = staged[5] ? &source_views[5] : NULL;
        model_sources.reference_fasta = reference;
        model_sources.transcript_coverage_complete = coverage;
        ok = duckvep_core_model_install(state->registry, name, &model_sources, "duckvep_model_publish",
                                        message, sizeof message) != 0;
        for (size_t i = 0; i < 6; ++i) {
            host_v2_stage_destroy(staged[i]);
            staged[i] = NULL;
        }
        if (!ok) {
            INPUT_ERROR(message);
        }
        ((bool *)output_data)[row] = true;
        mark_valid(output_validity, row);
        duckvep_budget_free(name);
        duckvep_budget_free(reference);
        name = reference = NULL;
    }
cleanup:
    for (size_t i = 0; i < 6; ++i) {
        host_v2_stage_destroy(staged[i]);
    }
    duckvep_budget_free(name);
    duckvep_budget_free(reference);
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_model_drop(name): removes the model and any staged rows
 * ------------------------------------------------------------------------- */

static void drop_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                      duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    duckdb_v2_vector_view views[1];
    duckdb_v2_vector_handle output = NULL;
    void *output_data = NULL;
    uint64_t *output_validity = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 1, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &output_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &output_validity, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_model_entry_t *entry = NULL, *previous = NULL;
        char *name = NULL;
        if (row_is_valid(&views[0], row)) {
            duckdb_v2_str text = string_at(&views[0], row);
            name = duckvep_core_string_copy(text.ptr, text.len);
        }
        if (!name || !*name) {
            duckvep_budget_free(name);
            INPUT_ERROR("duckvep_model_drop: name must be a non-empty string");
        }
        state_discard_model(state, name);
        pthread_mutex_lock(&state->registry->mutex);
        for (entry = state->registry->models; entry; previous = entry, entry = entry->next) {
            if (strcmp(entry->name, name) == 0) {
                break;
            }
        }
        if (entry && entry->pins == 0) {
            if (previous) {
                previous->next = entry->next;
            } else {
                state->registry->models = entry->next;
            }
            entry->next = NULL;
        } else {
            entry = NULL;
        }
        pthread_mutex_unlock(&state->registry->mutex);
        duckvep_budget_free(name);
        ((bool *)output_data)[row] = entry != NULL;
        mark_valid(output_validity, row);
        duckvep_model_entry_destroy(entry);
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_model_save(name, path) and duckvep_model_restore(name, path): model
 * snapshots (src/core/duckvep_core_snapshot.c), the same code as on v1
 * ------------------------------------------------------------------------- */

static void snapshot_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_error_info_handle *error,
                          bool restore) {
    const char *label = restore ? "duckvep_model_restore" : "duckvep_model_save";
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    duckdb_v2_vector_view views[2];
    duckdb_v2_vector_handle output = NULL;
    void *output_data = NULL;
    uint64_t *output_validity = NULL;
    char message[DUCKVEP_SQL_ERROR_SIZE + 256];
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 2, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &output_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &output_validity, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        char *name = NULL, *path = NULL;
        int done = 0;
        if (row_is_valid(&views[0], row) && row_is_valid(&views[1], row)) {
            duckdb_v2_str name_text = string_at(&views[0], row), path_text = string_at(&views[1], row);
            name = duckvep_core_string_copy(name_text.ptr, name_text.len);
            path = duckvep_core_string_copy(path_text.ptr, path_text.len);
        }
        message[0] = '\0';
        if (!name || !*name || !path || !*path) {
            (void)snprintf(message, sizeof message, "%s: name and path must be non-empty strings", label);
        } else if (restore) {
            done = duckvep_core_model_snapshot_install(state->registry, name, path, message, sizeof message);
        } else {
            duckvep_model_entry_t *entry = duckvep_registry_pin(state->registry, name);
            if (!entry) {
                (void)snprintf(message, sizeof message, "%s: unknown model name", label);
            } else {
                done = duckvep_core_model_snapshot_save(&entry->model, path, message, sizeof message);
                duckvep_registry_unpin(state->registry, entry);
            }
        }
        duckvep_budget_free(name);
        duckvep_budget_free(path);
        if (!done) {
            report(*error, message);
            goto cleanup;
        }
        ((bool *)output_data)[row] = true;
        mark_valid(output_validity, row);
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

static void save_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                      duckdb_v2_error_info_handle *error) {
    (void)context;
    snapshot_exec(info, error, false);
}

static void restore_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                         duckdb_v2_error_info_handle *error) {
    (void)context;
    snapshot_exec(info, error, true);
}

/* ---------------------------------------------------------------------------
 * _duckvep_model_fingerprint(name): FNV-1a of the loaded arrays (a diagnostic
 * for comparing a v1 and a v2 load of the same inputs), NULL for no such model
 * ------------------------------------------------------------------------- */

static void fingerprint_exec(duckdb_v2_scalar_function_exec_info_handle info,
                             duckdb_v2_context_handle context, duckdb_v2_error_info_handle *error) {
    (void)context;
    duckdb_v2_error_info_handle detail = NULL;
    model_state *state = NULL;
    duckdb_v2_vector_view views[1];
    duckdb_v2_vector_handle output = NULL;
    uint64_t *data = NULL;
    idx_t rows = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_user_data(info, (void **)&state, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    if (!load_views(info, 1, views, error)) {
        goto cleanup;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, (void **)&data, &detail));
    for (idx_t row = 0; row < rows; ++row) {
        duckvep_model_entry_t *entry = NULL;
        char *name = NULL;
        if (row_is_valid(&views[0], row)) {
            duckdb_v2_str text = string_at(&views[0], row);
            name = duckvep_core_string_copy(text.ptr, text.len);
        }
        if (name) {
            entry = duckvep_registry_pin(state->registry, name);
        }
        duckvep_budget_free(name);
        if (!entry) {
            DUCKDB_CALL(duckdb_v2_vector_set_null(output, row, &detail));
            continue;
        }
        data[row] = duckvep_core_model_fingerprint(&entry->model);
        duckvep_registry_unpin(state->registry, entry);
    }
cleanup:
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * duckvep_model_load_sql(name, regions_query, transcripts_query, exons_query [, options])
 *   -> VARCHAR[]: the statements to run, in order
 * options: mature_mirna_query, peptide_edit_query, interval_feature_query,
 * reference_fasta (VARCHAR) and transcript_coverage_complete (BOOLEAN)
 * ------------------------------------------------------------------------- */

typedef struct {
    duckvep_sql_text statements[DUCKVEP_MODEL_SCRIPT_MAX];
    size_t count;
} script;

static char *argument_copy(const column *c, idx_t row) {
    duckdb_v2_str text;
    if (!column_valid(c, row)) {
        return NULL;
    }
    text = string_at(&c->view, row);
    return duckvep_core_string_copy(text.ptr, text.len);
}

static void load_sql_exec(duckdb_v2_scalar_function_exec_info_handle info, duckdb_v2_context_handle context,
                          duckdb_v2_error_info_handle *error) {
    (void)context;
    static const char *const keys[] = {"mature_mirna_query", "peptide_edit_query", "interval_feature_query",
                                       "reference_fasta", "transcript_coverage_complete"};
    static const duckvep_core_option_kind_t kinds[] = {
        DUCKVEP_CORE_OPTION_TEXT, DUCKVEP_CORE_OPTION_TEXT, DUCKVEP_CORE_OPTION_TEXT,
        DUCKVEP_CORE_OPTION_TEXT, DUCKVEP_CORE_OPTION_BOOLEAN};
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_vector_handle arguments[5] = {0};
    column strings[4];
    options option;
    duckdb_v2_vector_handle output = NULL, child = NULL;
    void *list_data = NULL, *child_data = NULL;
    uint64_t *list_validity = NULL;
    duckdb_v2_arena_handle arena = NULL;
    idx_t rows = 0, at = 0;
    uint32_t argc = 0;
    script *scripts = NULL;
    size_t total = 0;
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_row_count(info, &rows, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg_count(info, &argc, &detail));
    for (uint32_t i = 0; i < argc; ++i) {
        DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_arg(info, i, &arguments[i], &detail));
    }
    for (uint32_t i = 0; i < 4; ++i) {
        if (!column_open(arguments[i], false, &strings[i], error)) {
            goto cleanup;
        }
    }
    if (!options_open(argc == 5 ? arguments[4] : NULL, keys, kinds, 5, &option, error)) {
        goto cleanup;
    }
    scripts = calloc(rows ? rows : 1, sizeof(*scripts));
    if (!scripts) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep_model_load_sql: out of memory");
        goto cleanup;
    }
    for (idx_t row = 0; row < rows; ++row) {
        char *values[4] = {0};
        char *optional[3] = {0};
        char *reference = NULL;
        const char *queries[6] = {0};
        int coverage = -1;
        duckvep_cell_t cell;
        const char *failure = NULL;
        bool built = false;
        if (!options_check(&option, row, error)) {
            goto cleanup;
        }
        for (uint32_t i = 0; i < 4; ++i) {
            values[i] = argument_copy(&strings[i], row);
            if (!values[i] || !*values[i]) {
                failure = "duckvep_model_load_sql: arguments must be non-empty strings";
            }
        }
        for (size_t o = 0; !failure && o < 3; ++o) {
            if (option_cell(&option, o, row, &cell) && cell.valid) {
                optional[o] = duckvep_core_string_copy(cell.text, cell.text_length);
                if (!optional[o] || !*optional[o]) {
                    failure = "duckvep_model_load_sql: optional queries must be non-empty strings";
                }
            }
        }
        if (!failure && option_cell(&option, 3, row, &cell) && cell.valid) {
            reference = duckvep_core_string_copy(cell.text, cell.text_length);
            if (!reference || !*reference) {
                failure = "duckvep_model_load_sql: reference_fasta must be a non-empty string";
            }
        }
        if (!failure && option_cell(&option, 4, row, &cell)) {
            if (!cell.valid) {
                failure = "duckvep_model_load_sql: transcript_coverage_complete cannot be NULL";
            } else {
                coverage = cell.boolean ? 1 : 0;
            }
        }
        if (!failure) {
            queries[0] = values[1];
            queries[1] = values[2];
            queries[2] = values[3];
            queries[3] = optional[0];
            queries[4] = optional[1];
            queries[5] = optional[2];
            built = duckvep_core_model_script(values[0], queries, reference, coverage, scripts[row].statements,
                                              &scripts[row].count);
            if (!built) {
                failure = "duckvep_model_load_sql: out of memory";
            }
        }
        for (size_t i = 0; i < 4; ++i) {
            duckvep_budget_free(values[i]);
        }
        for (size_t i = 0; i < 3; ++i) {
            duckvep_budget_free(optional[i]);
        }
        duckvep_budget_free(reference);
        if (failure) {
            report(*error, failure);
            goto cleanup;
        }
        total += scripts[row].count;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_exec_get_result(info, &output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flatten(output, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(output, rows, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(output, &list_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(output, &list_validity, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_child(output, 0, &child, &detail));
    DUCKDB_CALL(duckdb_v2_vector_set_size(child, total, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_data_mutable(child, &child_data, &detail));
    DUCKDB_CALL(duckdb_v2_vector_get_arena(child, &arena, &detail));
    {
        uint64_t *child_validity = NULL;
        DUCKDB_CALL(duckdb_v2_vector_flat_get_validity_mutable(child, &child_validity, &detail));
        for (size_t i = 0; i < total; ++i) {
            mark_valid(child_validity, i);
        }
    }
    for (idx_t row = 0; row < rows; ++row) {
        ((duckdb_v2_list_entry *)list_data)[row] = (duckdb_v2_list_entry){at, scripts[row].count};
        mark_valid(list_validity, row);
        for (size_t i = 0; i < scripts[row].count; ++i, ++at) {
            DUCKDB_CALL(write_string(arena, &((duckdb_v2_bytes *)child_data)[at], scripts[row].statements[i].data,
                                     scripts[row].statements[i].length, &detail));
        }
    }
cleanup:
    if (scripts) {
        for (idx_t row = 0; row < rows; ++row) {
            for (size_t i = 0; i < scripts[row].count; ++i) {
                duckvep_sql_free(&scripts[row].statements[i]);
            }
        }
        free(scripts);
    }
    (void)duckdb_v2_error_info_destroy(&detail);
}

/* ---------------------------------------------------------------------------
 * Registration
 * ------------------------------------------------------------------------- */

static bool register_model_scalar(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                                  model_state *state, const char *name, const char *const *parameter_types,
                                  const char *const *parameter_names, idx_t parameter_count,
                                  const char *return_type,
                                  duckdb_v2_scalar_function_exec_callback_fn callback,
                                  duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_scalar_function_handle function = NULL;
    duckdb_v2_function_signature_handle signature = NULL;
    duckdb_v2_logical_type_handle type = NULL;
    duckdb_v2_str text = string_view(name);
    bool success = false;
    bool retained = false;
    DUCKDB_CALL(duckdb_v2_scalar_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_name(function, &text, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_get_signature(function, &signature, &detail));
    for (idx_t i = 0; i < parameter_count; ++i) {
        DUCKDB_CALL(make_type(context, parameter_types[i], &type, &detail));
        DUCKDB_CALL(duckdb_v2_function_signature_add_parameter(signature, string_view(parameter_names[i]), type,
                                                               NULL, &detail));
        DUCKDB_CALL(duckdb_v2_logical_type_destroy(&type));
    }
    DUCKDB_CALL(make_type(context, return_type, &type, &detail));
    DUCKDB_CALL(duckdb_v2_function_signature_set_return_type(signature, type, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(
        function, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY, DUCKDB_V2_FUNCTION_PROPERTY_STABILITY_VOLATILE,
        &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_set_property(
        function, DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING, DUCKDB_V2_FUNCTION_PROPERTY_NULL_HANDLING_SPECIAL,
        &detail));
    if (state) {
        duckdb_v2_opaque data = {state, host_v2_state_release, NULL};
        host_v2_state_retain(state);
        retained = true;
        DUCKDB_CALL(duckdb_v2_scalar_function_set_user_data(function, &data, &detail));
        retained = false;
    }
    DUCKDB_CALL(duckdb_v2_scalar_function_set_exec_callback(function, callback, &detail));
    DUCKDB_CALL(duckdb_v2_scalar_function_register(function, &detail));
    success = true;
cleanup:
    if (retained) {
        host_v2_state_release(state);
    }
    (void)duckdb_v2_logical_type_destroy(&type);
    (void)duckdb_v2_scalar_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

static bool register_stage_format(duckdb_v2_extension_handle extension, model_state *state,
                                  duckdb_v2_error_info_handle *error) {
    duckdb_v2_error_info_handle detail = NULL;
    duckdb_v2_copy_function_handle function = NULL;
    duckdb_v2_identifier_t name = {"duckvep_stage", 13};
    bool success = false;
    bool retained = false;
    DUCKDB_CALL(duckdb_v2_copy_function_create_with_extension(extension, &function, &detail));
    DUCKDB_CALL(duckdb_v2_copy_function_set_name(function, &name, &detail));
    {
        duckdb_v2_opaque data = {state, host_v2_state_release, NULL};
        host_v2_state_retain(state);
        retained = true;
        DUCKDB_CALL(duckdb_v2_copy_function_set_user_data(function, &data, &detail));
        retained = false;
    }
    DUCKDB_CALL(duckdb_v2_copy_to_set_bind_callback(function, stage_bind_exec, &detail));
    DUCKDB_CALL(duckdb_v2_copy_to_set_init_callback(function, stage_init_exec, &detail));
    DUCKDB_CALL(duckdb_v2_copy_to_set_batch_callback(function, stage_batch_exec, &detail));
    DUCKDB_CALL(duckdb_v2_copy_to_set_flush_callback(function, stage_flush_exec, &detail));
    DUCKDB_CALL(duckdb_v2_copy_to_set_finalize_callback(function, stage_finalize_exec, &detail));
    DUCKDB_CALL(duckdb_v2_copy_function_register(function, &detail));
    success = true;
cleanup:
    if (retained) {
        host_v2_state_release(state);
    }
    (void)duckdb_v2_copy_function_destroy(&function);
    (void)duckdb_v2_error_info_destroy(&detail);
    return success;
}

bool host_v2_register_model(duckdb_v2_extension_handle extension, duckdb_v2_context_handle context,
                            duckdb_v2_error_info_handle *error) {
    static const char *const publish_types[] = {"VARCHAR", "ANY"};
    static const char *const publish_names[] = {"name", "options"};
    static const char *const one_type[] = {"VARCHAR"};
    static const char *const one_name[] = {"name"};
    static const char *const snapshot_types[] = {"VARCHAR", "VARCHAR"};
    static const char *const snapshot_names[] = {"name", "path"};
    static const char *const load_types[] = {"VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "ANY"};
    static const char *const load_names[] = {"name", "regions_query", "transcripts_query", "exons_query",
                                             "options"};
    model_state *state = calloc(1, sizeof(*state));
    bool ok = false;
    if (!state) {
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep: out of memory");
        return false;
    }
    pthread_mutex_init(&state->lock, NULL);
    state->references = 1;
    state->registry = duckvep_core_registry_create();
    if (!state->registry) {
        pthread_mutex_destroy(&state->lock);
        free(state);
        set_error(*error, DUCKDB_V2_ERROR_RESOURCE_OUT_OF_MEMORY, "duckvep: out of memory");
        return false;
    }
    ok = register_stage_format(extension, state, error) &&
         register_model_scalar(extension, context, state, "duckvep_model_publish", publish_types, publish_names,
                               1, "BOOLEAN", publish_exec, error) &&
         register_model_scalar(extension, context, state, "duckvep_model_publish", publish_types, publish_names,
                               2, "BOOLEAN", publish_exec, error) &&
         register_model_scalar(extension, context, state, "duckvep_model_drop", one_type, one_name, 1, "BOOLEAN",
                               drop_exec, error) &&
         register_model_scalar(extension, context, state, "duckvep_model_save", snapshot_types, snapshot_names, 2,
                               "BOOLEAN", save_exec, error) &&
         register_model_scalar(extension, context, state, "duckvep_model_restore", snapshot_types,
                               snapshot_names, 2, "BOOLEAN", restore_exec, error) &&
         register_model_scalar(extension, context, state, "_duckvep_model_fingerprint", one_type, one_name, 1,
                               "UBIGINT", fingerprint_exec, error) &&
         register_model_scalar(extension, context, NULL, "duckvep_model_load_sql", load_types, load_names, 4,
                               "VARCHAR[]", load_sql_exec, error) &&
         register_model_scalar(extension, context, NULL, "duckvep_model_load_sql", load_types, load_names, 5,
                               "VARCHAR[]", load_sql_exec, error);
    ok = ok && host_v2_register_annotate(extension, context, state, error);
    ok = ok && host_v2_register_haplotypes(extension, context, state, error);
    host_v2_state_release(state); /* the registrations hold their own references */
    return ok;
}

duckvep_registry_t *host_v2_registry_of(void *user_data) {
    return ((model_state *)user_data)->registry;
}
