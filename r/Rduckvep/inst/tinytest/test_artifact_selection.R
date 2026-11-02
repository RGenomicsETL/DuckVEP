library(tinytest)

root <- tempfile("duckvep-artifacts-")
dir.create(root)
args <- list(config = list(allow_unsigned_extensions = "true",
                           autoinstall_known_extensions = "false",
                           autoload_known_extensions = "false"))
if ("shared_home" %in% names(formals(duckdb::duckdb))) args$shared_home <- FALSE
con <- DBI::dbConnect(do.call(duckdb::duckdb, args))
platform <- DBI::dbGetQuery(con, "PRAGMA platform")$platform[[1L]]
original <- system.file(
  "duckvep_extension", "artifacts", paste0("v", as.character(utils::packageVersion("duckdb"))),
  "duckvep.duckdb_extension", package = "Rduckvep"
)
dir.create(file.path(root, "v0.0.1"))
stable <- file.path(root, "v0.0.1", "duckvep.duckdb_extension")
expect_true(file.copy(original, stable))

write_metadata <- getFromNamespace("duckvep_append_extension_metadata", "Rduckvep")
make_artifact <- function(version, abi, api, target_platform) {
  path <- file.path(root, version, "duckvep.duckdb_extension")
  dir.create(dirname(path))
  library_file <- tempfile("dummy-library-")
  writeBin(as.raw(1:4), library_file)
  write_metadata(library_file, path, duckdb_platform = target_platform,
                 duckdb_version = api, extension_version = version, abi_type = abi)
  unlink(library_file)
  path
}
unstable <- make_artifact("v2.0.0", "C_STRUCT_UNSTABLE", "v2.0.0", platform)
other_platform <- make_artifact("v3.0.0", "C_STRUCT", "v1.2.0", "different_platform")
v2_stable <- make_artifact("v4.0.0", "C_STRUCT", "v2.0.0", platform)
dir.create(file.path(root, "v5.0.0"))
writeBin(as.raw(1:4), file.path(root, "v5.0.0", "duckvep.duckdb_extension"))

candidates <- Rduckvep:::duckvep_stable_artifacts(root, platform)
expect_identical(candidates, stable)
expect_false(any(c(unstable, other_platform, v2_stable) %in% candidates))
expect_identical(Rduckvep:::duckvep_stable_artifacts(root, "not_this_platform"), character())
expect_identical(Rduckvep:::duckvep_stable_artifacts(file.path(root, "missing"), platform), character())
expect_identical(rduckvep_load(con, candidates[[1L]]), con)
expect_true(DBI::dbGetQuery(con, "SELECT duckvep_allele_geometry(100, 'A', 'ATG').interbase AS is_interbase")$is_interbase)
DBI::dbDisconnect(con, shutdown = TRUE)
unlink(root, recursive = TRUE)
