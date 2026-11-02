targets_file <- local({
  frame_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile) || length(frame$ofile) != 1L) return("")
    as.character(frame$ofile)
  }, character(1))
  frame_scripts <- vapply(sys.frames(), function(frame) {
    if (!exists("script", envir = frame, inherits = FALSE)) return("")
    script <- get("script", envir = frame, inherits = FALSE)
    if (!is.character(script) || length(script) != 1L) return("")
    script
  }, character(1))
  argument_files <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  )
  candidates <- c(frame_files, frame_scripts, argument_files)
  candidates <- candidates[nzchar(candidates) & basename(candidates) == "_targets.R"]
  if (length(candidates) == 0L) {
    stop("could not locate _targets.R; set DUCKVEP_ROOT to a DuckVEP checkout")
  }
  normalizePath(candidates[[length(candidates)]], mustWork = TRUE)
})
source(file.path(dirname(targets_file), "R", "root.R"))
duckvep_root <- duckvep_targets_root(targets_file)

library(targets)

source(file.path(duckvep_root, "scripts", "duckvep_evidence.R"))
source(file.path(duckvep_root, "pipelines", "duckvep", "R", "pipeline.R"))

campaign_manifest <- Sys.getenv(
  "DUCKVEP_TARGETS_CAMPAIGNS",
  unset = file.path(duckvep_root, "pipelines", "duckvep", "campaigns.tsv")
)
campaign_manifest <- duckvep_targets_path(campaign_manifest, duckvep_root)
output_root <- Sys.getenv(
  "DUCKVEP_TARGETS_OUTPUT",
  unset = file.path(duckvep_root, ".tmp", "duckvep_targets", "results")
)
output_root <- duckvep_targets_path(output_root, duckvep_root)
targets_store <- Sys.getenv(
  "DUCKVEP_TARGETS_STORE",
  unset = file.path(duckvep_root, ".tmp", "duckvep_targets", "store")
)
targets_store <- duckvep_targets_path(targets_store, duckvep_root)
extension_receipt <- file.path(
  duckvep_root,
  ".tmp",
  "duckvep_targets",
  "extension.tsv"
)
default_micromamba <- Sys.getenv("MICROMAMBA", unset = "micromamba")
default_vep_prefix <- Sys.getenv("VEP_PREFIX")

tar_option_set(
  packages = "blit",
  memory = "transient",
  garbage_collection = TRUE,
  error = "continue",
  workspace_on_error = TRUE
)

list(
  tar_target(
    duckvep_source_revision,
    duckvep_targets_source_revision(duckvep_root),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    duckvep_campaign_manifest,
    normalizePath(campaign_manifest, mustWork = TRUE),
    format = "file"
  ),
  tar_target(
    duckvep_runtime_defaults,
    list(
      micromamba = default_micromamba,
      vep_prefix = default_vep_prefix
    )
  ),
  tar_target(
    duckvep_campaigns,
    duckvep_targets_read_campaigns(
      duckvep_campaign_manifest,
      duckvep_root,
      duckvep_runtime_defaults$micromamba,
      duckvep_runtime_defaults$vep_prefix,
      output_root,
      c(
        file.path(
          duckvep_root,
          "build",
          "release",
          "duckvep.duckdb_extension"
        ),
        extension_receipt
      ),
      targets_store
    ),
    iteration = "list"
  ),
  tar_target(
    duckvep_extension_bundle,
    duckvep_targets_build_extension(
      duckvep_root,
      duckvep_source_revision,
      extension_receipt
    ),
    format = "file"
  ),
  tar_target(
    duckvep_campaign,
    duckvep_campaigns,
    pattern = map(duckvep_campaigns),
    iteration = "list"
  ),
  tar_target(
    duckvep_campaign_inputs,
    duckvep_targets_campaign_inputs(duckvep_campaign),
    pattern = map(duckvep_campaign),
    format = "file"
  ),
  tar_target(
    duckvep_cache_state,
    duckvep_targets_cache_state(duckvep_campaign),
    pattern = map(duckvep_campaign),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    duckvep_campaign_evidence,
    duckvep_targets_run_campaign(
      duckvep_campaign,
      duckvep_campaign_inputs,
      duckvep_cache_state,
      duckvep_extension_bundle,
      duckvep_root
    ),
    pattern = map(
      duckvep_campaign,
      duckvep_campaign_inputs,
      duckvep_cache_state
    ),
    format = "file"
  )
)
