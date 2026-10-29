sql_quote_string <- function(con, value) {
  as.character(DBI::dbQuoteString(con, value))
}

build_param_str <- function(params) {
  paste0(",", paste(paste0(names(params), ":=", unlist(params)), collapse = ","))
}

.duckvep_check_table_target <- function(con, table_name, overwrite) {
  if (!is.null(table_name) &&
      (!is.character(table_name) || length(table_name) != 1L ||
       is.na(table_name) || !nzchar(table_name))) {
    stop("table_name must be one nonempty string", call. = FALSE)
  }
  if (!is.null(table_name) && !overwrite && DBI::dbExistsTable(con, table_name)) {
    stop("Table already exists: ", table_name, call. = FALSE)
  }
}

# The C API v2 host registers caller-side staging builders in place of the v1 table functions.
.duckvep_is_v2 <- function(con) {
  nrow(DBI::dbGetQuery(con, "SELECT 1 FROM duckdb_functions()
    WHERE function_name = 'duckvep_model_publish' LIMIT 1")) > 0L
}

# A v2 staging job name unique within this process; staged jobs are shared by every
# connection to the database, so a fixed name would let concurrent calls replace each other.
.duckvep_job_counter <- new.env(parent = emptyenv())
.duckvep_job_name <- function(prefix) {
  count <- (if (is.null(.duckvep_job_counter$n)) 0 else .duckvep_job_counter$n) + 1
  .duckvep_job_counter$n <- count
  paste0(prefix, "_", Sys.getpid(), "_", format(count, scientific = FALSE), "_",
         format(as.numeric(Sys.time()) * 1e6, scientific = FALSE))
}

# CREATE OR REPLACE swaps the table in only once the query succeeds, so a failed query
# leaves the existing table in place.
.duckvep_create_table <- function(con, table_name, query, overwrite) {
  name <- as.character(DBI::dbQuoteIdentifier(con, table_name))
  create <- if (overwrite) "CREATE OR REPLACE TABLE" else "CREATE TABLE"
  DBI::dbExecute(con, paste(create, name, "AS", query))
}
