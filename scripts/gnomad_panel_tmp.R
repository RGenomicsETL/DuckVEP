# Shared by the panel builders. DuckDB spill goes to $DUCKVEP_SCALE_TMP (default
# /root/duckvep/data/scale-tmp), never under the staged gnomAD tree, and is removed when the
# process ends. The spill quota is capped so that at least DUCKVEP_MIN_FREE_GIB (50) GiB stay
# free on the filesystem; the build refuses to start if less than 2 GiB of quota would remain.
panel_spill_setup <- function(execute, q) {
  scale_tmp <- Sys.getenv("DUCKVEP_SCALE_TMP", "/root/duckvep/data/scale-tmp")
  dir.create(scale_tmp, recursive = TRUE, showWarnings = FALSE)
  spill <- tempfile("panel-spill-", tmpdir = scale_tmp)
  dir.create(spill)
  reg.finalizer(globalenv(), function(e) unlink(spill, recursive = TRUE), onexit = TRUE)
  parse_gib <- function(x) {
    n <- as.numeric(sub("[^0-9.].*$", "", x))
    unit <- toupper(sub("^[0-9.]+", "", x))
    n * c(GIB = 1, GB = 1, MIB = 1 / 1024, MB = 1 / 1024, TIB = 1024, TB = 1024)[[unit]]
  }
  requested <- parse_gib(Sys.getenv("DUCKVEP_PANEL_TEMP", "12GiB"))
  free <- as.numeric(strsplit(system2("df", c("-PB1", scale_tmp), stdout = TRUE)[2L], "[[:space:]]+")[[1L]][4L]) / 2^30
  allowed <- min(requested, floor(free - as.numeric(Sys.getenv("DUCKVEP_MIN_FREE_GIB", "50"))))
  if (allowed < 2) stop(sprintf("only %.1f GiB free; the panel build needs the disk to stay above %s GiB free",
    free, Sys.getenv("DUCKVEP_MIN_FREE_GIB", "50")))
  execute(paste0("SET temp_directory=", q(spill)))
  execute(paste0("SET max_temp_directory_size=", q(paste0(allowed, "GiB"))))
  message(sprintf("spill %s, quota %g GiB (%.1f GiB free)", spill, allowed, free))
  invisible(spill)
}
