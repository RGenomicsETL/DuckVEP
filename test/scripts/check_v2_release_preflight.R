#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))]
if (length(script_arg) != 1L) stop("run this file with Rscript", call. = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."), mustWork = TRUE)

usage <- function() {
  cat(paste0(
    "Usage: Rscript test/scripts/check_v2_release_preflight.R [options]\n\n",
    "Checks a concrete v2 release qualification. It performs no network access unless\n",
    "--live-status is supplied. ADMITTED requires live upstream status, an actual\n",
    "extension and pinned CLI, an explicitly selected R library, and an HG002 run\n",
    "directory produced by scripts/check_v2_hg002.sh.\n\n",
    "  --extension PATH       v2 extension artifact to inspect and load\n",
    "  --v2-cli PATH          DuckDB CLI used to inspect and load the artifact\n",
    "  --r-library PATH       library path containing the R duckdb package to query\n",
    "  --hg002-dir PATH       completed check_v2_hg002.sh output directory\n",
    "  --live-status          explicitly query DuckDB GitHub and CRAN metadata\n",
    "  --offline-status PATH  captured status JSON; useful for diagnostics only\n",
    "  --selftest             run offline negative controls\n"
  ))
}

parse_args <- function(values) {
  result <- list(extension = NULL, v2_cli = NULL, r_library = NULL, hg002_dir = NULL,
                 live_status = FALSE, offline_status = NULL, selftest = FALSE)
  while (length(values) > 0L) {
    key <- values[[1L]]
    if (key %in% c("-h", "--help")) {
      usage()
      quit(status = 0L)
    }
    if (key %in% c("--live-status", "--selftest")) {
      result[[gsub("-", "_", sub("^--", "", key))]] <- TRUE
      values <- values[-1L]
      next
    }
    field <- gsub("-", "_", sub("^--", "", key))
    if (!field %in% names(result) || length(values) < 2L) {
      stop(sprintf("unknown option or missing value: %s", key), call. = FALSE)
    }
    result[[field]] <- values[[2L]]
    values <- values[-c(1L, 2L)]
  }
  result
}

run <- function(command, arguments = character(), env = character()) {
  output <- suppressWarnings(system2(command, shQuote(arguments), stdout = TRUE, stderr = TRUE, env = env))
  list(ok = identical(attr(output, "status"), NULL) || attr(output, "status") == 0L,
       output = output)
}

source_matches <- function(actual, pinned) {
  nchar(actual) >= 10L && startsWith(pinned, actual)
}

read_receipt <- function(path) {
  table <- tryCatch(read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                               stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(table) || ncol(table) != 2L || anyDuplicated(table[[1L]])) return(NULL)
  setNames(table[[2L]], table[[1L]])
}

check_budget <- function(path) {
  table <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  required <- c("current_bytes", "high_water_bytes", "limit_bytes")
  if (is.null(table) || nrow(table) == 0L || !all(required %in% names(table)) ||
      !all(vapply(table[required], is.numeric, logical(1L)))) return(FALSE)
  bytes <- as.matrix(table[required])
  all(is.finite(bytes)) && all(bytes >= 0) && all(bytes == floor(bytes)) &&
    all(table$current_bytes <= table$high_water_bytes) &&
    all(table$high_water_bytes <= table$limit_bytes)
}

check_exact_rows <- function(path) {
  table <- tryCatch(read.csv(path, colClasses = "character"), error = function(e) NULL)
  if (is.null(table) || nrow(table) != 1L ||
      !identical(names(table), c("left_rows", "right_rows", "row_differences"))) return(FALSE)
  counts <- unlist(table, use.names = FALSE)
  all(!is.na(counts)) && all(grepl("^(0|[1-9][0-9]*)$", counts)) &&
    identical(table$left_rows[[1L]], table$right_rows[[1L]]) &&
    identical(table$row_differences[[1L]], "0")
}

check_fingerprints <- function(directory) {
  paths <- file.path(directory, c("v1-model-fingerprint.csv", "v2-model-fingerprint.csv"))
  tables <- lapply(paths, function(path) tryCatch(read.csv(path, colClasses = "character"), error = function(e) NULL))
  !any(vapply(tables, is.null, logical(1L))) &&
    all(vapply(tables, function(x) identical(names(x), "fingerprint") && nrow(x) == 1L &&
      !is.na(x$fingerprint[[1L]]) && nzchar(x$fingerprint[[1L]]), logical(1L))) &&
    identical(tables[[1L]]$fingerprint[[1L]], tables[[2L]]$fingerprint[[1L]])
}

check_input_hashes <- function(path, required) {
  lines <- readLines(path, warn = FALSE)
  fields <- regexec("^([0-9a-f]{64})  (.+)$", lines)
  matches <- regmatches(lines, fields)
  if (!length(lines) || any(lengths(matches) != 3L)) return(FALSE)
  expected <- vapply(matches, `[[`, "", 2L)
  inputs <- vapply(matches, `[[`, "", 3L)
  if (any(!file.exists(c(inputs, required))) || anyDuplicated(inputs) ||
      !all(normalizePath(required) %in% normalizePath(inputs))) return(FALSE)
  observed <- run("sha256sum", c("--", inputs))
  if (!observed$ok || length(observed$output) != length(inputs)) return(FALSE)
  observed_hashes <- sub("  .*$", "", observed$output)
  identical(unname(expected), unname(observed_hashes))
}

hg002_checks <- function(directory, pin, extension, cli_source) {
  checks <- c("HG002 output directory supplied" = !is.null(directory))
  if (is.null(directory) || !dir.exists(directory)) {
    checks["HG002 output directory exists"] <- FALSE
    return(checks)
  }
  directory <- normalizePath(directory)
  receipt_path <- file.path(directory, "receipt.tsv")
  receipt <- if (file.exists(receipt_path)) read_receipt(receipt_path) else NULL
  checks["HG002 receipt is a valid key-value table"] <- !is.null(receipt)
  if (is.null(receipt)) return(checks)
  required <- c("status", "v2_pinned_revision", "v2_duckdb_source", "v1_duckdb", "v2_duckdb",
                "v1_extension", "v2_extension", "vcf", "model", "v1_parquet", "v2_parquet",
                "v1_v2_exact_multiset")
  checks["HG002 receipt has qualification fields"] <- all(required %in% names(receipt))
  if (!all(required %in% names(receipt))) return(checks)
  checks["HG002 receipt status is PASS"] <- identical(receipt[["status"]], "PASS")
  checks["HG002 receipt uses the repository SDK revision"] <- identical(receipt[["v2_pinned_revision"]], pin)
  checks["HG002 runtime source matches the supplied CLI"] <- source_matches(receipt[["v2_duckdb_source"]], pin) &&
    identical(receipt[["v2_duckdb_source"]], cli_source)
  checks["HG002 receipt is bound to the supplied extension"] <- !is.null(extension) && file.exists(extension) &&
    identical(normalizePath(receipt[["v2_extension"]], mustWork = FALSE), normalizePath(extension))
  checks["HG002 records exact v1-v2 parity"] <- identical(receipt[["v1_v2_exact_multiset"]], "PASS") &&
    check_exact_rows(file.path(directory, "v1-v2-rows.csv"))
  checks["HG002 records equal model fingerprints"] <- check_fingerprints(directory)
  budget_files <- file.path(directory, c("v1-budget-after-load.csv", "v1-budget-final.csv",
                                          "v2-budget-after-load.csv", "v2-budget-after-capture.csv",
                                          "v2-budget-final.csv"))
  checks["HG002 budget snapshots stay within their limits"] <- all(file.exists(budget_files)) &&
    all(vapply(budget_files, check_budget, logical(1L)))
  concrete <- c(receipt[["vcf"]], receipt[["model"]], receipt[["v1_parquet"]], receipt[["v2_parquet"]])
  checks["HG002 comparison, model, and output paths exist"] <- all(file.exists(concrete)) &&
    all(file.info(concrete)$size > 0L)
  checks["HG002 input hashes match concrete inputs"] <- file.exists(file.path(directory, "inputs.sha256")) &&
    check_input_hashes(file.path(directory, "inputs.sha256"),
      unname(receipt[c("v1_duckdb", "v2_duckdb", "v1_extension", "v2_extension", "vcf", "model")]))
  checks
}

r_runtime <- function(library_path, pin, extension) {
  if (is.null(library_path)) return(list(ok = FALSE, detail = "--r-library was not supplied"))
  if (!dir.exists(library_path)) return(list(ok = FALSE, detail = "R library does not exist"))
  library_path <- normalizePath(library_path)
  old_paths <- .libPaths()
  on.exit(.libPaths(old_paths), add = TRUE)
  .libPaths(c(library_path, old_paths))
  tryCatch({
    if (!requireNamespace("DBI", quietly = TRUE) ||
        !requireNamespace("duckdb", lib.loc = library_path, quietly = TRUE))
      stop("duckdb must be installed in --r-library and DBI must be available")
    expected_path <- normalizePath(find.package("duckdb", lib.loc = library_path))
    if (!identical(normalizePath(getNamespaceInfo("duckdb", "path")), expected_path))
      stop("loaded duckdb namespace is not from --r-library")
    version <- as.character(utils::packageVersion("duckdb", lib.loc = library_path))
    if (!identical(version, "2.0.0")) stop("R duckdb package version is ", version, "; expected 2.0.0")
    driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
      config = list(allow_unsigned_extensions = "true"))
    connection <- DBI::dbConnect(driver)
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    runtime <- DBI::dbGetQuery(connection, "SELECT version() AS version, source_id FROM pragma_version()")
    ready <- nrow(runtime) == 1L && identical(sub("^v", "", runtime$version[[1L]]), "2.0.0") &&
      source_matches(runtime$source_id[[1L]], pin)
    if (!ready) stop("R duckdb runtime does not match the released SDK pin")
    if (is.null(extension) || !file.exists(extension)) stop("--extension was not supplied or does not exist")
    DBI::dbExecute(connection, paste("LOAD", DBI::dbQuoteString(connection, normalizePath(extension))))
    loaded <- DBI::dbGetQuery(connection,
      "SELECT loaded FROM duckdb_extensions() WHERE extension_name = 'duckvep'")
    list(ok = nrow(loaded) == 1L && isTRUE(loaded$loaded[[1L]]),
      detail = sprintf("package=%s runtime=%s source=%s; extension load inspected", version,
        runtime$version[[1L]], runtime$source_id[[1L]]))
  }, error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
}

live_status <- function() {
  release <- run("gh", c("api", "repos/duckdb/duckdb/releases/tags/v2.0.0"))
  if (!release$ok) return(list(ok = FALSE, detail = paste(release$output, collapse = "\n")))
  release_json <- tryCatch(jsonlite::fromJSON(paste(release$output, collapse = "\n")), error = function(e) NULL)
  cran <- tryCatch(utils::available.packages(repos = c(CRAN = "https://cloud.r-project.org")), error = function(e) NULL)
  ok <- !is.null(release_json) && identical(release_json$tag_name, "v2.0.0") && !is.null(cran) &&
    "duckdb" %in% rownames(cran) && identical(unname(cran["duckdb", "Version"]), "2.0.0")
  list(ok = ok, detail = if (is.null(cran) || !"duckdb" %in% rownames(cran)) "DuckDB release found; CRAN duckdb unavailable" else
    sprintf("GitHub tag=%s CRAN duckdb=%s", release_json$tag_name, cran["duckdb", "Version"]))
}

write_bad_footer <- function(path) {
  fields <- c("4", "linux", "v2.0.0", "duckvep", "C_STRUCT_UNSTABLE", "", "", "")
  bytes <- raw(512L)
  for (i in seq_along(fields)) {
    value <- charToRaw(fields[[length(fields) - i + 1L]])
    start <- (i - 1L) * 32L + 1L
    bytes[start + seq_along(value) - 1L] <- value
  }
  writeBin(c(charToRaw("not-an-extension"), bytes), path)
}

selftest <- function() {
  pin <- jsonlite::fromJSON(file.path(root, "duckvep-package.json"))$v2_host
  sdk <- run("python3", file.path(root, "scripts", "fetch-v2-sdk.py"))
  stopifnot(sdk$ok, source_matches(pin$duckdb_sdk_revision, pin$duckdb_sdk_revision),
            !source_matches("deadbeef00", pin$duckdb_sdk_revision))
  directory <- tempfile("release-preflight-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  footer <- file.path(directory, "wrong-footer.duckdb_extension")
  write_bad_footer(footer)
  old <- Sys.getenv("DUCKVEP_V2_EXTENSION", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("DUCKVEP_V2_EXTENSION") else Sys.setenv(DUCKVEP_V2_EXTENSION = old), add = TRUE)
  Sys.setenv(DUCKVEP_V2_EXTENSION = footer)
  footer_check <- run("python3", c(file.path(root, "test", "scripts", "check_v2_host.py"), "footer"))
  stopifnot(!footer_check$ok, any(grepl("C_STRUCT_UNSTABLE", footer_check$output, fixed = TRUE)))
  missing <- hg002_checks(tempfile("missing-hg002-"), pin$duckdb_sdk_revision, footer, pin$duckdb_sdk_revision)
  stopifnot(!all(missing))
  rows <- file.path(directory, "comparison.csv")
  for (values in c("2,2,1", "2,3,0", "-1,-1,0", "9007199254740992,9007199254740993,0")) {
    writeLines(c("left_rows,right_rows,row_differences", values), rows)
    stopifnot(!check_exact_rows(rows))
  }
  writeLines(c("left_rows,right_rows,row_differences", "2,2,0"), rows)
  stopifnot(check_exact_rows(rows))
  fingerprint_paths <- file.path(directory, c("v1-model-fingerprint.csv", "v2-model-fingerprint.csv"))
  for (path in fingerprint_paths) writeLines(c("fingerprint", "8650913055114199712"), path)
  stopifnot(check_fingerprints(directory))
  writeLines(c("fingerprint", "8650913055114199713"), fingerprint_paths[[2L]])
  stopifnot(!check_fingerprints(directory))
  budget <- file.path(directory, "budget.csv")
  writeLines(c("current_bytes,high_water_bytes,limit_bytes", "1,2,3"), budget)
  stopifnot(check_budget(budget))
  writeLines(c("current_bytes,high_water_bytes,limit_bytes", "-1,2,3"), budget)
  stopifnot(!check_budget(budget))
  hashes <- file.path(directory, "inputs.sha256")
  hash <- run("sha256sum", c("--", footer))
  stopifnot(hash$ok)
  writeLines(hash$output, hashes)
  stopifnot(check_input_hashes(hashes, footer), !check_input_hashes(hashes, c(footer, rows)))
  cat("RELEASE-PREFLIGHT selftest: source/ABI, missing qualification, exact counts, 64-bit fingerprints, budgets and hash coverage passed\n")
}

options <- parse_args(args)
if (options$selftest) {
  selftest()
  quit(status = 0L)
}

pin <- jsonlite::fromJSON(file.path(root, "duckvep-package.json"))$v2_host
checks <- list()
add_check <- function(name, passed, detail) {
  checks[[name]] <<- list(passed = isTRUE(passed), detail = detail)
}

sdk <- run("python3", file.path(root, "scripts", "fetch-v2-sdk.py"))
add_check("repository SDK manifest and headers", sdk$ok, paste(sdk$output, collapse = "\n"))
add_check("repository v2 host is a release", identical(pin$status, "release"), sprintf("manifest status=%s", pin$status))

if (!is.null(options$offline_status)) {
  status <- tryCatch(jsonlite::fromJSON(options$offline_status), error = function(e) NULL)
  add_check("offline status is diagnostic only", FALSE,
            if (is.null(status)) "offline status JSON is invalid" else "offline captured status cannot admit a release")
}
if (options$live_status) {
  status <- live_status()
  add_check("explicit live DuckDB and CRAN status", status$ok, status$detail)
} else {
  add_check("explicit live DuckDB and CRAN status", FALSE, "not requested; use --live-status")
}

cli_source <- NULL
if (!is.null(options$v2_cli) && file.exists(options$v2_cli) && file.access(options$v2_cli, 1L) == 0L) {
  version <- run(options$v2_cli, c("-no-init", "-batch", "-list", "-noheader", "-c",
                                   "SELECT version() || '|' || source_id FROM pragma_version()"))
  fields <- if (version$ok && length(version$output) == 1L) strsplit(version$output[[1L]], "|", fixed = TRUE)[[1L]] else character()
  cli_source <- if (length(fields) == 2L) fields[[2L]] else NULL
  add_check("supplied CLI is DuckDB 2.0.0 at the SDK revision", length(fields) == 2L &&
              identical(sub("^v", "", fields[[1L]]), "2.0.0") && source_matches(fields[[2L]], pin$duckdb_sdk_revision),
            paste(version$output, collapse = "\n"))
} else {
  add_check("supplied CLI is DuckDB 2.0.0 at the SDK revision", FALSE, "--v2-cli is missing, unreadable, or not executable")
}

if (!is.null(options$extension) && file.exists(options$extension)) {
  old <- Sys.getenv("DUCKVEP_V2_EXTENSION", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("DUCKVEP_V2_EXTENSION") else Sys.setenv(DUCKVEP_V2_EXTENSION = old), add = TRUE)
  Sys.setenv(DUCKVEP_V2_EXTENSION = normalizePath(options$extension))
  footer <- run("python3", c(file.path(root, "test", "scripts", "check_v2_host.py"), "footer"))
  add_check("actual extension footer has stable v2 ABI", footer$ok, paste(footer$output, collapse = "\n"))
  if (!is.null(options$v2_cli) && file.exists(options$v2_cli)) {
    escaped <- gsub("'", "''", normalizePath(options$extension), fixed = TRUE)
    load <- run(options$v2_cli, c("-unsigned", "-no-init", "-batch", "-bail", "-c", sprintf("LOAD '%s';", escaped)))
    add_check("actual extension loads in the supplied CLI", load$ok, paste(load$output, collapse = "\n"))
  } else {
    add_check("actual extension loads in the supplied CLI", FALSE, "--v2-cli is not available")
  }
} else {
  add_check("actual extension footer has stable v2 ABI", FALSE, "--extension is missing or unreadable")
  add_check("actual extension loads in the supplied CLI", FALSE, "no extension artifact supplied")
}

r_check <- r_runtime(options$r_library, pin$duckdb_sdk_revision, options$extension)
add_check("selected R duckdb is 2.0.0 at the SDK revision and loads the extension", r_check$ok, r_check$detail)

if (is.null(cli_source)) cli_source <- ""
hg002 <- hg002_checks(options$hg002_dir, pin$duckdb_sdk_revision, options$extension, cli_source)
for (name in names(hg002)) add_check(name, hg002[[name]], "")

for (name in names(checks)) {
  check <- checks[[name]]
  cat(sprintf("%s: %s%s\n", if (check$passed) "PASS" else "BLOCKED", name,
              if (nzchar(check$detail)) paste0(" — ", check$detail) else ""))
}
all_pass <- all(vapply(checks, `[[`, logical(1L), "passed"))
if (all_pass) {
  cat("ADMITTED: release inputs and qualification artifacts were inspected\n")
  quit(status = 0L)
}
ready <- checks[["repository v2 host is a release"]]$passed &&
  checks[["explicit live DuckDB and CRAN status"]]$passed
outcome <- "BLOCKED: release evidence is incomplete or inconsistent"
if (ready) outcome <- "READY_TO_QUALIFY: upstream releases are visible; artifact and HG002 qualification is incomplete"
cat(outcome, "\n", sep = "")
quit(status = 1L)
