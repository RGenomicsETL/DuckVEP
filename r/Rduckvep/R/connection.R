#' Load the bundled DuckVEP extension
#' @param con An open DuckDB connection with unsigned extension loading enabled.
#' @param extension_path Optional extension artifact for this DuckDB version.
#' @return Invisibly, the connection.
#' @export
rduckvep_load <- function(con, extension_path = NULL) {
  if (is.null(extension_path)) {
    extension_path <- system.file(
      "duckvep_extension", "artifacts",
      paste0("v", as.character(utils::packageVersion("duckdb"))),
      "duckvep.duckdb_extension", package = "Rduckvep"
    )
  }
  if (!length(extension_path) || !nzchar(extension_path) || !file.exists(extension_path)) {
    stop("DuckVEP artifact for installed DuckDB is not available", call. = FALSE)
  }
  DBI::dbExecute(con, paste("LOAD", as.character(DBI::dbQuoteString(con, extension_path))))
  invisible(con)
}

#' Open DuckDB with DuckVEP loaded
#' @param dbdir Database path, or `:memory:`.
#' @param extension_path Optional locally built extension artifact.
#' @return A DBI connection. Disconnect with `shutdown = TRUE`.
#' @export
rduckvep_connect <- function(dbdir = ":memory:", extension_path = NULL) {
  config <- list(allow_unsigned_extensions = "true",
                 autoinstall_known_extensions = "false",
                 autoload_known_extensions = "false")
  args <- list(dbdir = dbdir, config = config)
  supported <- names(formals(duckdb::duckdb))
  if ("allow_extensions" %in% supported) args$allow_extensions <- TRUE
  if ("shared_home" %in% supported) args$shared_home <- FALSE
  con <- DBI::dbConnect(do.call(duckdb::duckdb, args))
  tryCatch(rduckvep_load(con, extension_path), error = function(error) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop(error)
  })
  con
}
