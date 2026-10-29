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

# CREATE OR REPLACE swaps the table in only once the query succeeds, so a failed query
# leaves the existing table in place.
.duckvep_create_table <- function(con, table_name, query, overwrite) {
  name <- as.character(DBI::dbQuoteIdentifier(con, table_name))
  create <- if (overwrite) "CREATE OR REPLACE TABLE" else "CREATE TABLE"
  DBI::dbExecute(con, paste(create, name, "AS", query))
}
