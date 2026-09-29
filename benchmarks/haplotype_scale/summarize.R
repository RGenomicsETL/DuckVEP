#!/usr/bin/env Rscript
# Turns the raw receipts of run_qualification.sh into the committed CSVs and prints the tables and the verdict.
#
#   Rscript benchmarks/haplotype_scale/summarize.R RAW_DIR [OUTPUT_DIR]
#
# RAW_DIR holds process.tsv (one row per fresh process from capped_run.sh) and <label>.worker.tsv (one row per pass).
# OUTPUT_DIR (default benchmarks/data/haplotype_scale) receives runs.csv (one row per process and pass) and summary.csv
# (medians and ratios). The gate is judged on mode B, cold process, model load included, wall clock of the whole process
# (R and DuckDB start-up included) against the median wall clock of the pinned csq on the same input and core.
args <- commandArgs(trailingOnly = TRUE)
raw <- args[[1L]]
here <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
out <- if (length(args) >= 2L) args[[2L]] else file.path(dirname(here), "..", "data", "haplotype_scale")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
process <- read.delim(file.path(raw, "process.tsv"), stringsAsFactors = FALSE)
workers <- do.call(rbind, lapply(list.files(raw, "\\.worker\\.tsv$", full.names = TRUE), function(f) {
  x <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE); x
}))
extension_sha <- if (file.exists(file.path(raw, "ext", "sha256"))) readLines(file.path(raw, "ext", "sha256"), warn = FALSE)[1L] else NA
group_of <- function(label) sub("_[0-9]+$", "", label)
runs <- process
runs$group <- group_of(runs$label)
runs$tool <- ifelse(grepl("^csq", runs$group), "bcftools csq", "duckvep")
if (!is.null(workers) && nrow(workers)) {
  cold <- workers[workers$pass == "cold", ]
  warm <- workers[workers$pass != "cold", ]
  names(cold)[names(cold) == "label"] <- "worker_label"
  merged <- merge(runs, cold, by.x = "label", by.y = "worker_label", all.x = TRUE)
  merged$pass <- ifelse(is.na(merged$pass), NA, "cold")
  warm_rows <- NULL
  if (nrow(warm)) {
    warm$group <- group_of(warm$label)
    warm_rows <- merge(warm, runs[, c("label", "load1_at_start", "core", "sibling_busy_pct")], by = "label")
  }
  runs <- merged
  if (!is.null(warm_rows)) {
    fill <- setdiff(names(runs), names(warm_rows))
    for (column in fill) warm_rows[[column]] <- NA
    runs <- rbind(runs, warm_rows[, names(runs)])
  }
}
# Native budget per owner during the pass (MiB), from the BUDGET line the worker prints: current and high-water.
budget_of <- function(label) {
  file <- file.path(raw, paste0(label, ".out"))
  if (!file.exists(file)) return(NULL)
  line <- grep("^BUDGET\\t", readLines(file, warn = FALSE), value = TRUE)
  if (!length(line)) return(NULL)
  pairs <- strsplit(strsplit(sub("^BUDGET\\t", "", line[length(line)]), ";", fixed = TRUE)[[1L]], "=", fixed = TRUE)
  setNames(as.numeric(vapply(pairs, `[`, "", 2L)), vapply(pairs, `[`, "", 1L))
}
budgets <- lapply(runs$label, budget_of)
columns <- unique(unlist(lapply(budgets, names)))
for (column in columns) runs[[paste0(column, "_mib")]] <- vapply(budgets, function(b) if (is.null(b) || is.na(b[column])) NA_real_ else b[[column]], 0)
runs$extension_sha256 <- extension_sha
runs$host_note <- "i5-13500 (P-core 6, SMT sibling 7), 1 DuckDB thread, memory_limit 8GB, native budget 4 GiB, MemoryMax 16G"
write.csv(runs, file.path(out, "runs.csv"), row.names = FALSE)
med <- function(x) if (length(x[!is.na(x)])) median(x, na.rm = TRUE) else NA_real_
group_rows <- function(g, pass = NULL) {
  x <- runs[runs$group == g, ]
  if (!is.null(pass) && "pass" %in% names(x)) x <- x[!is.na(x$pass) & x$pass == pass, ]
  x
}
gate_csq <- med(group_rows("csq_Ou")$wall_s)
diagnostic_csq <- med(group_rows("csq_Ob")$wall_s)
summary <- data.frame(measure = character(), value = numeric(), note = character(), stringsAsFactors = FALSE)
add <- function(measure, value, note = "") summary[nrow(summary) + 1L, ] <<- list(measure, value, note)
add("csq_median_wall_s_Ou_devnull", gate_csq, "bcftools csq -p a -Ou -o /dev/null, the gate baseline (matches the recorded 22.2 s log: same warnings, 811296 KiB RSS)")
add("csq_median_wall_s_Ob_file", diagnostic_csq, "same with -Ob to a file: BCF compression adds about 12 s")
for (m in c("B", "A")) {
  g <- paste0("hg002_", m)
  cold <- group_rows(g, "cold"); warm <- group_rows(paste0(g, "W"), "warm"); coldw <- group_rows(paste0(g, "W"), "cold")
  add(paste0(g, "_cold_process_wall_s_median"), med(cold$wall_s), "cold-only fresh processes, whole process incl. R/DuckDB start-up")
  add(paste0(g, "_cold_internal_total_s_median"), med(cold$total_s), "model load + stage + predict")
  add(paste0(g, "_cold_load_s_median"), med(cold$load_s), "duckvep_model_load")
  add(paste0(g, "_cold_stage_s_median"), med(cold$stage_s), if (m == "B") "decode, discovery, calls" else "read the ordered calls")
  add(paste0(g, "_cold_predict_s_median"), med(cold$predict_s), "duckvep_haplotypes + Parquet")
  add(paste0(g, "_warm_pass_s_median"), med(warm$pass_s), "second pass in the same process, model resident")
  add(paste0(g, "_warm_stage_s_median"), med(warm$stage_s), "")
  add(paste0(g, "_warm_predict_s_median"), med(warm$predict_s), "")
  add(paste0(g, "W_cold_pass_internal_total_s_median"), med(coldw$total_s), "cold pass of the runs that also did a warm pass")
}
b_wall <- med(group_rows("hg002_B", "cold")$wall_s)
add("gate_ratio_B_cold_wall_over_csq_Ou", b_wall / gate_csq, "the signed gate is <= 0.5")
add("ratio_B_cold_wall_over_csq_Ob", b_wall / diagnostic_csq, "against csq writing BCF")
add("ratio_B_cold_internal_over_csq_Ou", med(group_rows("hg002_B", "cold")$total_s) / gate_csq, "without R/DuckDB start-up")
add("ratio_B_warm_over_csq_Ou", med(group_rows("hg002_BW", "warm")$pass_s) / gate_csq, "warm execution only")
bz <- group_rows("hg002_Bbgzip", "cold")
if (nrow(bz)) {
  add("hg002_Bbgzip_cold_process_wall_s_median", med(bz$wall_s), "diagnostic: VCF inflated by bgzip -dc (libdeflate) through a FIFO on the same core")
  add("hg002_Bbgzip_cold_stage_s_median", med(bz$stage_s), "")
  add("ratio_Bbgzip_cold_wall_over_csq_Ou", med(bz$wall_s) / gate_csq, "diagnostic only; the gate uses DuckDB's own decoder")
}
add("ratio_A_cold_wall_over_csq_Ou", med(group_rows("hg002_A", "cold")$wall_s) / gate_csq, "mode A")
add("ratio_A_warm_over_csq_Ou", med(group_rows("hg002_AW", "warm")$pass_s) / gate_csq, "mode A warm")
add("gate_met", as.numeric(b_wall / gate_csq <= 0.5), "1 when mode B cold wall <= 0.5 x csq median")
add("gate_target_wall_s", 0.5 * gate_csq, "")
for (g in c("hg002_inside_A", "hg002_outside_A")) {
  x <- group_rows(g, "cold")
  add(paste0(g, "_predict_s_median"), med(x$predict_s), "duckvep_haplotypes + Parquet on that partition of the calls")
}
for (g in c("qual5m_B", "qual5m_A", "dense_B", "dense_A", "low_sharing_B", "low_sharing_A")) {
  cold <- group_rows(g, "cold"); warm <- group_rows(paste0(sub("_([AB])$", "_\\1W", g)), "warm")
  add(paste0(g, "_cold_process_wall_s_median"), med(cold$wall_s), "")
  add(paste0(g, "_cold_internal_total_s_median"), med(cold$total_s), "")
  add(paste0(g, "_cold_load_s_median"), med(cold$load_s), "")
  add(paste0(g, "_cold_stage_s_median"), med(cold$stage_s), "")
  add(paste0(g, "_cold_predict_s_median"), med(cold$predict_s), "")
  if (nrow(warm)) add(paste0(g, "_warm_pass_s_median"), med(warm$pass_s), "")
  add(paste0(g, "_cgroup_peak_bytes_max"), suppressWarnings(max(as.numeric(cold$cgroup_peak_bytes), na.rm = TRUE)), "memory.peak")
  add(paste0(g, "_max_rss_kib_max"), suppressWarnings(max(as.numeric(cold$max_rss_kib), na.rm = TRUE)), "GNU time")
  add(paste0(g, "_spill_peak_bytes_max"), suppressWarnings(max(as.numeric(cold$spill_peak_bytes), na.rm = TRUE)), "sampled temp directory size")
}
for (g in sort(unique(sub("_[0-9]+$", "", runs$label[grepl("^ab_", runs$label)])))) {
  x <- group_rows(g, "cold")
  add(paste0(g, "_load_s_median"), med(x$load_s), "model load (mode A, cold)")
  add(paste0(g, "_predict_s_median"), med(x$predict_s), "duckvep_haplotypes + Parquet")
}
for (g in sort(unique(sub("_[0-9]+$", "", runs$label[grepl("^abx_", runs$label)])))) {
  x <- group_rows(g, "cold")
  add(paste0(g, "_load_s_median"), med(x$load_s), "abx: A/B of a pipeline choice, cold")
  add(paste0(g, "_stage_s_median"), med(x$stage_s), "")
  add(paste0(g, "_process_wall_s_median"), med(x$wall_s), "")
}
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 4)
cat(sprintf("\nVERDICT: mode B cold wall median %.2f s vs csq -Ou median %.2f s: ratio %.3f (gate <= 0.500): %s\n",
  b_wall, gate_csq, b_wall / gate_csq, if (b_wall / gate_csq <= 0.5) "MET" else "NOT MET"))
