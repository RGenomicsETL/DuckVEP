duckvep_stable_artifacts <- function(root, platform) {
  if (!nzchar(root) || !dir.exists(root)) return(character())
  paths <- file.path(list.dirs(root, recursive = FALSE), "duckvep.duckdb_extension")
  paths <- sort(paths[file.exists(paths)])
  stable <- vapply(paths, function(path) {
    tryCatch({
      con <- file(path, "rb")
      on.exit(close(con))
      seek(con, where = -512, origin = "end")
      footer <- readBin(con, "raw", n = 512L)
      if (length(footer) != 512L) return(FALSE)
      field <- function(i) {
        bytes <- footer[seq.int((i - 1L) * 32L + 1L, i * 32L)]
        rawToChar(bytes[bytes != as.raw(0L)])
      }
      identical(field(8L), "4") && identical(field(7L), platform) &&
        identical(field(4L), "C_STRUCT") && grepl("^v1\\.[0-9]+\\.[0-9]+$", field(6L))
    }, error = function(e) FALSE)
  }, logical(1))
  paths[stable]
}

#' Load the bundled DuckVEP extension
#' @param con An open DuckDB connection with unsigned extension loading enabled.
#' @param extension_path Optional locally built extension artifact.
#' @details If an exact-version artifact is unavailable, the loader tries bundled
#'   stable v1 artifacts for the same DuckDB platform. DuckDB checks their API
#'   compatibility. Version-locked artifacts are not used as fallbacks.
#' @return Invisibly, the connection.
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbGetQuery dbQuoteIdentifier dbQuoteString dbExistsTable
#' @export
rduckvep_load <- function(con, extension_path = NULL) {
  if (!is.null(extension_path)) {
    if (length(extension_path) == 0L || !nzchar(extension_path) || !file.exists(extension_path)) {
      stop("DuckVEP artifact for installed DuckDB is not available", call. = FALSE)
    }
    DBI::dbExecute(con, paste("LOAD", as.character(DBI::dbQuoteString(con, extension_path))))
    return(invisible(con))
  }
  root <- system.file("duckvep_extension", "artifacts", package = "Rduckvep")
  exact <- file.path(root, paste0("v", as.character(utils::packageVersion("duckdb"))),
                     "duckvep.duckdb_extension")
  if (nzchar(root) && file.exists(exact)) {
    DBI::dbExecute(con, paste("LOAD", as.character(DBI::dbQuoteString(con, exact))))
    return(invisible(con))
  }
  platform <- DBI::dbGetQuery(con, "PRAGMA platform")$platform[[1L]]
  candidates <- duckvep_stable_artifacts(root, platform)
  if (length(candidates) == 0L) {
    stop("DuckVEP artifact for installed DuckDB is not available", call. = FALSE)
  }
  for (candidate in candidates) {
    result <- tryCatch(
      DBI::dbExecute(con, paste("LOAD", as.character(DBI::dbQuoteString(con, candidate)))),
      error = identity
    )
    if (!inherits(result, "error")) return(invisible(con))
  }
  stop("No compatible stable DuckVEP artifact could be loaded: ", conditionMessage(result), call. = FALSE)
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
