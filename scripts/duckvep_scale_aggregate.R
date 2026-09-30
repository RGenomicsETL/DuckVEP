#!/usr/bin/env Rscript
# Merge the per-job result, cgroup and load files of one scale run into receipt.csv and
# summary.md. Classifies every (job, mode) as ok, capacity_error or failed; an OOM kill
# (memory.events oom_kill > 0) is always failed.
#
#   Rscript scripts/duckvep_scale_aggregate.R RUN_DIR
#
# Exit status: 0 when every job and mode is ok and met every ceiling, 3 otherwise.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: Rscript scripts/duckvep_scale_aggregate.R RUN_DIR")
dir <- normalizePath(args[[1L]], mustWork = TRUE)
read_tsv <- function(path) if (file.exists(path)) read.delim(path, colClasses = "character",
  na.strings = "", check.names = FALSE, quote = "") else NULL
run <- read_tsv(file.path(dir, "run.tsv"))
info <- setNames(run$value, run$key)
modes <- strsplit(info[["modes"]], ",", fixed = TRUE)[[1L]]
jobs <- seq_len(as.integer(info[["jobs"]]))
num <- function(x) suppressWarnings(as.numeric(x))
first <- function(x, name, default = NA_character_) {
  if (is.null(x) || !name %in% names(x) || nrow(x) == 0L) default else x[[name]][1L]
}
getinfo <- function(k, default = NA_character_) if (k %in% names(info)) info[[k]] else default
enforced <- identical(info[["enforced"]], "yes")
max_temp <- num(info[["max_temp_bytes"]])
budget_bytes <- num(info[["native_budget_mib"]]) * 1048576
rows <- list()
for (j in jobs) {
  jd <- file.path(dir, paste0("job-", j))
  res <- read_tsv(file.path(jd, "result.tsv"))
  cg <- read_tsv(file.path(jd, "cgroup.tsv"))
  ld <- read_tsv(file.path(jd, "load.tsv"))
  for (m in modes) {
    r <- if (is.null(res)) NULL else res[res$mode == m, , drop = FALSE]
    if (!is.null(r) && nrow(r) == 0L) r <- NULL
    outcome <- first(r, "outcome", "failed")
    reason <- first(r, "reason", "")
    exit_code <- first(cg, "exit_code")
    if (is.null(r)) reason <- paste0("no result reported (job exit code ", exit_code,
      if (is.null(cg)) ", wrapper lost" else "", "; see job.log)")
    oom_kill <- num(first(cg, "oom_kill", NA))
    peak <- num(first(cg, "memory_peak", NA)); mmax <- num(first(cg, "memory_max", NA))
    spill <- num(first(cg, "spill_peak_bytes", NA))
    hw_total <- num(first(r, "native_hw_total", NA))
    checks <- character()
    if (!is.na(oom_kill) && oom_kill > 0) { outcome <- "failed"; checks <- c(checks, paste0("cgroup oom_kill=", oom_kill)) }
    if (!is.na(peak) && !is.na(mmax) && peak > mmax) checks <- c(checks, "memory.peak above memory.max")
    if (!is.na(spill) && spill > max_temp) checks <- c(checks, "spill above quota")
    if (!is.na(hw_total) && hw_total > budget_bytes) checks <- c(checks, "native budget exceeded")
    # Mandatory evidence: a row that is ok but lacks any of these cannot certify.
    if (outcome == "ok") {
      missing <- c(if (is.na(hw_total)) "native_hw_total", if (is.na(spill)) "spill peak",
        if (is.na(num(first(r, "mode_s", NA)))) "mode seconds",
        if (anyNA(c(first(r, "out_rows", NA), first(r, "hash_sum", NA), first(r, "hash_xor", NA)))) "checksum")
      if (enforced) {
        swapmax <- first(cg, "memory_swap_max", NA)
        missing <- c(missing, if (is.na(peak)) "memory.peak", if (is.na(mmax)) "memory.max",
          if (is.na(swapmax)) "memory.swap.max", if (is.na(oom_kill)) "oom_kill")
        if (!is.na(swapmax) && !identical(trimws(swapmax), "0")) checks <- c(checks, "swap not disabled")
      }
      if (length(missing)) checks <- c(checks, paste0(missing, " missing for job ", j))
    }
    if (enforced && (is.na(peak) || is.na(mmax))) checks <- unique(c(checks, "cgroup counters unavailable"))
    if (length(checks) && outcome == "ok") outcome <- "failed"
    if (length(checks)) reason <- paste(c(if (!is.na(reason) && nzchar(reason)) reason, checks), collapse = "; ")
    row <- data.frame(run_id = info[["run_id"]], job = j, mode = m, outcome = outcome, reason = reason,
      ceilings = if (enforced) "enforced" else "ceilings not enforced",
      ceilings_met = outcome == "ok" && enforced,
      panel = info[["panel"]], regulation = info[["regulation"]], regulation_features = first(r, "regulation_features"), panel_rows = first(r, "panel_rows"), threads = info[["threads_per_job"]],
      load_s = first(ld, "load_s"), mode_s = first(r, "mode_s"), job_wall_s = first(cg, "wall_s"),
      alleles_per_s = first(r, "alleles_per_s"), out_rows = first(r, "out_rows"),
      out_bytes = first(r, "out_bytes"), hash_sum = first(r, "hash_sum"), hash_xor = first(r, "hash_xor"),
      memory_peak_bytes = first(cg, "memory_peak"), memory_max_bytes = first(cg, "memory_max"),
      memory_swap_max = first(cg, "memory_swap_max"), memory_swap_peak_bytes = first(cg, "memory_swap_peak"), oom_kill = first(cg, "oom_kill"),
      spill_peak_bytes = first(cg, "spill_peak_bytes"), rss_peak_kib = first(r, "rss_peak_kib"),
      stringsAsFactors = FALSE)
    for (o in c("model", "index", "reference", "workspace", "scratch", "emit", "control", "total")) {
      row[[paste0("native_load_hw_", o)]] <- first(ld, paste0("native_load_hw_", o))
      row[[paste0("native_hw_", o)]] <- first(r, paste0("native_hw_", o))
    }
    rows[[length(rows) + 1L]] <- row
  }
}
receipt <- do.call(rbind, rows)
write.csv(receipt, file.path(dir, "receipt.csv"), row.names = FALSE, na = "")

# Per-mode aggregates over the concurrent jobs.
ok <- receipt$outcome == "ok"
gib <- function(x) sprintf("%.2f", num(x) / 1073741824)
mib <- function(x) sprintf("%.1f", num(x) / 1048576)
pct <- function(x, p) if (length(x) && !all(is.na(x))) unname(quantile(x, p, na.rm = TRUE, type = 1)) else NA_real_
agg <- do.call(rbind, lapply(modes, function(m) {
  s <- receipt[receipt$mode == m, ]
  good <- s$outcome == "ok"
  t <- num(s$mode_s[good]); w <- num(s$job_wall_s)
  agree <- length(unique(paste(s$out_rows[good], s$hash_sum[good], s$hash_xor[good]))) <= 1L && any(good)
  data.frame(mode = m, jobs_ok = sum(good), jobs = nrow(s),
    aggregate_alleles_per_s = if (any(good)) round(sum(num(s$panel_rows[good])) / max(t)) else NA_real_,
    mode_s_p50 = pct(t, 0.5), mode_s_max = if (any(good)) max(t) else NA_real_,
    job_wall_s_p50 = pct(w, 0.5), job_wall_s_max = if (any(!is.na(w))) max(w, na.rm = TRUE) else NA_real_,
    checksums_agree = agree, stringsAsFactors = FALSE)
}))
# Certification: every job and mode ok, checksums agree within each mode (2 or more ok jobs), the
# ceilings were enforced, and the frozen artifacts did not change. Anything else is exit 3 with a reason.
cert_reasons <- character()
if (!enforced) cert_reasons <- c(cert_reasons, "ceilings not enforced")
if (any(receipt$outcome != "ok")) cert_reasons <- c(cert_reasons, "not every job and mode is ok")
for (m in modes) {
  s <- receipt[receipt$mode == m & receipt$outcome == "ok", ]
  if (nrow(s) >= 2L && length(unique(paste(s$out_rows, s$hash_sum, s$hash_xor))) > 1L)
    cert_reasons <- c(cert_reasons, paste0("checksums disagree in ", m))
}
for (a in c("model", "panel", "extension")) {
  start <- getinfo(paste0(a, "_sha256")); end <- getinfo(paste0(a, "_sha256_end"))
  if (is.na(start) || is.na(end)) cert_reasons <- c(cert_reasons, paste0(a, " hash missing at start or end"))
  else if (!identical(start, end)) cert_reasons <- c(cert_reasons, paste0(a, " changed during the run"))
}
all_met <- length(cert_reasons) == 0L
counts <- table(factor(receipt$outcome, c("ok", "capacity_error", "failed")))
hashes <- function(a) paste0(getinfo(paste0(a, "_sha256"), "missing"), " / ", getinfo(paste0(a, "_sha256_end"), "missing"))
tbl <- function(df) {
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1L, function(r) paste0("| ", paste(ifelse(is.na(r), "", r), collapse = " | "), " |")))
}
lines <- c(paste0("# DuckVEP scale run ", info[["run_id"]]), "",
  if (!enforced) c("**ceilings not enforced** (`--cgroup none`): the process memory and swap ceilings were not applied.", ""),
  paste0("Certified: **", if (all_met) "yes" else "no", "**. Outcomes: ok ", counts[["ok"]],
    ", capacity_error ", counts[["capacity_error"]], ", failed ", counts[["failed"]], "."), "",
  if (length(cert_reasons)) c("Not certified:", paste0("- ", cert_reasons), ""),
  "## Run", "",
  tbl(data.frame(item = c("host", "cores", "RAM (GiB)", "kernel", "jobs x threads", "panel", "panel rows",
    "regulation (resident features)", "panel sha256 (start / end)", "model sha256 (start / end)", "extension sha256 (start / end)", "git revision (dirty files)", "cgroup mode", "wall (s)",
    "load average at start"),
    value = c(info[["host"]], info[["cores"]], gib(info[["mem_total_bytes"]]), info[["kernel"]],
      paste0(info[["jobs"]], " x ", info[["threads_per_job"]]), paste0(info[["panel"]], " (", info[["limit_rows"]], ")"),
      first(receipt[ok, , drop = FALSE], "panel_rows", "n/a"),
      if (identical(info[["regulation"]], "1")) paste0("on (", first(receipt, "regulation_features"), ")") else "off", 
      hashes("panel"), hashes("model"), hashes("extension"), paste0(info[["git_revision"]], " (", info[["git_dirty"]], ")"),
      info[["cgroup_mode"]], info[["wall_s"]], info[["loadavg_at_start"]]), stringsAsFactors = FALSE)), "",
  "## Ceilings", "",
  paste0("Process memory.max ", gib(info[["memory_max_bytes"]]), " GiB with memory.swap.max 0; DuckDB memory_limit ",
    info[["duckdb_memory_limit"]], ", max_temp_directory_size ", gib(info[["max_temp_bytes"]]),
    " GiB; native budget ", info[["native_budget_mib"]], " MiB, ", info[["workers"]], " workers x ",
    info[["scratch_mib"]], " MiB scratch + ", info[["emit_mib"]], " MiB emit."), "",
  "## Aggregate", "", tbl(agg), "",
  "Aggregate alleles/s is the sum of the jobs' panel rows over the slowest job's annotation-and-write time; it excludes model load.", "",
  "## Per job", "",
  tbl(data.frame(job = receipt$job, mode = receipt$mode, outcome = receipt$outcome,
    load_s = receipt$load_s, mode_s = receipt$mode_s, alleles_per_s = receipt$alleles_per_s,
    out_rows = receipt$out_rows, out_MiB = mib(receipt$out_bytes),
    memory_peak_GiB = gib(receipt$memory_peak_bytes), oom_kill = receipt$oom_kill,
    spill_peak_MiB = mib(receipt$spill_peak_bytes), native_hw_total_MiB = mib(receipt$native_hw_total),
    native_hw_model_MiB = mib(receipt$native_hw_model), stringsAsFactors = FALSE)), "",
  "## Checksums", "",
  tbl(data.frame(job = receipt$job, mode = receipt$mode, out_rows = receipt$out_rows,
    hash_sum = receipt$hash_sum, hash_xor = receipt$hash_xor, stringsAsFactors = FALSE)), "",
  "The checksum is the row count with the sum and XOR of DuckDB `hash()` over every output row, independent of row order.")
bad <- receipt[receipt$outcome != "ok", ]
if (nrow(bad)) lines <- c(lines, "## Non-ok jobs", "",
  paste0("- job ", bad$job, " ", bad$mode, ": **", bad$outcome, "**: ", bad$reason))
writeLines(lines, file.path(dir, "summary.md"))
quit(save = "no", status = if (all_met) 0L else 3L)
