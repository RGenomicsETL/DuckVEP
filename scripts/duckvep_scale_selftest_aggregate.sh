#!/usr/bin/env bash
# Aggregator self-test: synthetic run directories, no DuckVEP run. A clean receipt set must exit 0;
# each contradiction or gap in the evidence must exit 3 with its reason in summary.md.
# Usage: scripts/duckvep_scale_selftest_aggregate.sh
set -euo pipefail
SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/duckvep-aggregate-selftest.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
export AGG="$SCRIPT_DIR/duckvep_scale_aggregate.R" WORK
rc=0
Rscript --vanilla - <<'RSCRIPT' || rc=$?
agg <- Sys.getenv("AGG"); work <- Sys.getenv("WORK")
w <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
make <- function(name, mut = function(x) x) {
  d <- file.path(work, name); dir.create(d)
  run <- c(run_id = name, host = "h", kernel = "k", cores = "8", mem_total_bytes = "68719476736", jobs = "2",
    panel = "1M", panel_file = "p", panel_sha256 = "aa", model_sha256 = "bb", extension_sha256 = "cc",
    panel_sha256_end = "aa", model_sha256_end = "bb", extension_sha256_end = "cc", git_revision = "g", git_dirty = "0",
    cgroup_mode = "systemd", enforced = "yes", regulation = "0", threads_per_job = "6", modes = "compact,complete17",
    limit_rows = "all", memory_max_bytes = "17179869184", duckdb_memory_limit = "8GB", max_temp_bytes = "34359738368",
    native_budget_mib = "4096", workers = "6", scratch_mib = "128", emit_mib = "256", wall_s = "1", loadavg_at_start = "0 0 0")
  res <- list(); cg <- list()
  for (j in 1:2) {
    res[[j]] <- data.frame(mode = c("compact", "complete17"), outcome = "ok", reason = "", panel_rows = "1000",
      mode_s = c("1.5", "2.5"), alleles_per_s = "500", out_rows = c("13", "14"), out_bytes = "100",
      hash_sum = c("111", "222"), hash_xor = c("333", "444"), rss_peak_kib = "1", native_hw_total = "1048576",
      stringsAsFactors = FALSE)
    cg[[j]] <- data.frame(exit_code = "0", wall_s = "5", memory_peak = "1000000", memory_max = "17179869184",
      memory_swap_max = "0", memory_swap_peak = "0", oom = "0", oom_kill = "0", event_max = "0", event_high = "0",
      spill_peak_bytes = "0", cgroup = "x", stringsAsFactors = FALSE)
  }
  s <- mut(list(run = run, res = res, cg = cg))
  w(data.frame(key = names(s$run), value = unname(s$run)), file.path(d, "run.tsv"))
  for (j in 1:2) {
    jd <- file.path(d, paste0("job-", j)); dir.create(jd)
    w(s$res[[j]], file.path(jd, "result.tsv")); w(s$cg[[j]], file.path(jd, "cgroup.tsv"))
    w(data.frame(load_s = "1"), file.path(jd, "load.tsv"))
  }
  d
}
fail <- 0L
expect <- function(name, code, pattern, mut = function(x) x) {
  d <- make(name, mut)
  out <- suppressWarnings(system2("Rscript", c(shQuote(agg), shQuote(d)), stdout = TRUE, stderr = TRUE))
  got <- attr(out, "status"); if (is.null(got)) got <- 0L
  sm <- paste(readLines(file.path(d, "summary.md")), collapse = "\n")
  ok <- identical(as.integer(got), code) && (is.null(pattern) || grepl(pattern, sm, fixed = TRUE))
  cat(if (ok) "PASS  " else "FAIL  ", name, " (exit ", got, if (!is.null(pattern)) paste0(", expects \"", pattern, "\"") else "", ")\n", sep = "")
  if (!ok) fail <<- 1L
}
expect("clean receipts certify", 0L, "Certified: **yes**")
expect("checksums disagree", 3L, "checksums disagree in complete17", function(s) { s$res[[2]]$hash_xor[2] <- "999"; s })
expect("missing native_hw_total", 3L, "native_hw_total missing for job 2", function(s) { s$res[[2]]$native_hw_total <- NA; s })
expect("counter missing in some rows only", 3L, "memory.peak missing for job 1", function(s) { s$cg[[1]]$memory_peak <- NA; s })
expect("swap.max not 0", 3L, "swap not disabled", function(s) { s$cg[[1]]$memory_swap_max <- "max"; s })
expect("oom_kill above 0", 3L, "oom_kill=1", function(s) { s$cg[[2]]$oom_kill <- "1"; s })
expect("model hash changed", 3L, "model changed during the run", function(s) { s$run[["model_sha256_end"]] <- "zz"; s })
expect("end hashes missing", 3L, "panel hash missing at start or end", function(s) { s$run <- s$run[names(s$run) != "panel_sha256_end"]; s })
expect("ceilings not enforced", 3L, "ceilings not enforced", function(s) { s$run[["enforced"]] <- "no: ceilings not enforced"; s })
quit(save = "no", status = fail)
RSCRIPT
if [[ "$rc" == 0 ]]; then echo "aggregate self-test: all checks passed"; else echo "aggregate self-test: FAILED"; fi
exit "$rc"
