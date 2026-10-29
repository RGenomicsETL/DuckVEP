#' Load the bundled DuckVEP extension
#' @param con An open DuckDB connection with unsigned extension loading enabled.
#' @param extension_path Optional extension artifact for this DuckDB version.
#' @return Invisibly, the connection.
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbGetQuery dbQuoteIdentifier dbQuoteString dbExistsTable
#' @export
rduckvep_load <- function(con, extension_path = NULL) {
  if (is.null(extension_path)) {
    extension_path <- .duckvep_bundled_artifact(as.character(utils::packageVersion("duckdb")))
  }
  if (!length(extension_path) || !nzchar(extension_path) || !file.exists(extension_path)) {
    stop("DuckVEP artifact for installed DuckDB is not available", call. = FALSE)
  }
  DBI::dbExecute(con, paste("LOAD", as.character(DBI::dbQuoteString(con, extension_path))))
  invisible(con)
}

# The artifact built for the installed duckdb, else the newest one bundled. Artifacts declare the
# stable C API v1.2.0, so one built against an older duckdb still loads after duckdb is upgraded;
# DuckDB itself refuses a footer it cannot load.
.duckvep_bundled_artifact <- function(duckdb_version) {
  artifacts <- system.file("duckvep_extension", "artifacts", package = "Rduckvep")
  exact <- file.path(artifacts, paste0("v", duckdb_version), "duckvep.duckdb_extension")
  if (nzchar(artifacts) && file.exists(exact)) return(exact)
  found <- Sys.glob(file.path(artifacts, "v*", "duckvep.duckdb_extension"))
  if (!nzchar(artifacts) || !length(found)) return("")
  versions <- package_version(sub("^v", "", basename(dirname(found))), strict = FALSE)
  if (all(is.na(versions))) return(found[[1L]])
  found[order(versions, decreasing = TRUE, na.last = TRUE)][[1L]]
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
