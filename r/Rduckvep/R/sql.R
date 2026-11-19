sql_quote_string <- function(con, value) {
  as.character(DBI::dbQuoteString(con, value))
}

build_param_str <- function(params) {
  paste0(",", paste(paste0(names(params), ":=", unlist(params)), collapse = ","))
}

# Admits the string arguments of an exported function. Each is passed by name,
# and the error names the first one that is not one nonempty string.
.duckvep_check_strings <- function(...) {
  values <- list(...)
  for (name in names(values)) {
    value <- values[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(name, " must be one nonempty string", call. = FALSE)
    }
  }
}

.duckvep_check_table_target <- function(con, table_name, overwrite) {
  if (!isTRUE(overwrite) && !isFALSE(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  if (is.null(table_name)) return(invisible())
  .duckvep_check_strings(table_name = table_name)
  if (!overwrite && DBI::dbExistsTable(con, table_name)) {
    stop("Table already exists: ", table_name, call. = FALSE)
  }
}

# The rows of `query`, or a table of them when `table_name` is given.
.duckvep_query_result <- function(con, query, table_name, overwrite) {
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  name <- as.character(DBI::dbQuoteIdentifier(con, table_name))
  create <- if (overwrite) "CREATE OR REPLACE TABLE" else "CREATE TABLE"
  DBI::dbExecute(con, paste(create, name, "AS", query))
  invisible(TRUE)
}

.duckvep_is_v2 <- function(con) {
  nrow(DBI::dbGetQuery(con, "SELECT 1 FROM duckdb_functions()
    WHERE function_name = 'duckvep_model_publish' LIMIT 1")) > 0L
}
