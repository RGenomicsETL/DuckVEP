#!/usr/bin/env Rscript
# Markdown tables of the committed receipts (benchmarks/data/haplotype_scale), the numbers quoted in its README.
#
#   Rscript benchmarks/haplotype_scale/tables.R [DATA_DIR] [RAW_DIR]
args <- commandArgs(trailingOnly = TRUE)
here <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
data <- if (length(args) >= 1L) args[[1L]] else file.path(dirname(here), "..", "data", "haplotype_scale")
raw <- if (length(args) >= 2L) args[[2L]] else NA
runs <- read.csv(file.path(data, "runs.csv"), stringsAsFactors = FALSE)
summary <- read.csv(file.path(data, "summary.csv"), stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
num <- function(x) suppressWarnings(as.numeric(x))
med <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
fmt <- function(x, digits = 2) ifelse(is.na(x), "-", formatC(x, format = "f", digits = digits, big.mark = ","))
md <- function(df) {
  cat("| ", paste(names(df), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(df)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) cat("| ", paste(as.character(unlist(df[i, ])), collapse = " | "), " |\n", sep = "")
  cat("\n")
}
audit_of <- function(text) {
  if (is.na(text) || !nzchar(text)) return(list())
  parts <- strsplit(strsplit(text, ";", fixed = TRUE)[[1L]], "=", fixed = TRUE)
  setNames(lapply(parts, function(p) paste(p[-1L], collapse = "=")), vapply(parts, `[`, "", 1L))
}

cat("### csq baseline (one core, 3 fresh runs each)\n\n")
csq <- runs[runs$tool == "bcftools csq", ]
md(do.call(rbind, lapply(split(csq, csq$group), function(x) data.frame(
  invocation = c(csq_Ou = "`-p a -Ou -o /dev/null`", csq_Ob = "`-p a -Ob -o file.bcf`")[[x$group[1L]]],
  `wall s (runs)` = paste(fmt(num(x$wall_s)), collapse = " / "), `median wall s` = fmt(med(num(x$wall_s))),
  `median user s` = fmt(med(num(x$user_s))), `peak RSS KiB` = fmt(max(num(x$max_rss_kib)), 0), `memory.peak MiB` = fmt(max(num(x$cgroup_peak_bytes)) / 1048576, 0),
  `load at start` = paste(fmt(num(x$load1_at_start), 1), collapse = " / "), check.names = FALSE))))

cat("### HG002 on the full Ensembl 116 model (single core)\n\n")
rows <- list()
for (m in c("B", "A")) {
  cold <- runs[runs$group == paste0("hg002_", m) & !is.na(runs$pass) & runs$pass == "cold", ]
  warm <- runs[runs$group == paste0("hg002_", m, "W") & !is.na(runs$pass) & runs$pass == "warm", ]
  rows[[length(rows) + 1L]] <- data.frame(mode = m,
    `cold process wall s (runs)` = paste(fmt(num(cold$wall_s)), collapse = " / "), `cold wall median` = fmt(med(num(cold$wall_s))),
    `load` = fmt(med(num(cold$load_s))), `stage` = fmt(med(num(cold$stage_s))), `predict + Parquet` = fmt(med(num(cold$predict_s))),
    `warm pass median` = fmt(med(num(warm$pass_s))), `warm stage` = fmt(med(num(warm$stage_s))), `warm predict` = fmt(med(num(warm$predict_s))), check.names = FALSE)
}
md(do.call(rbind, rows))
for (name in c("gate_ratio_B_cold_wall_over_csq_Ou", "ratio_B_cold_wall_over_csq_Ob", "ratio_B_cold_internal_over_csq_Ou", "ratio_B_warm_over_csq_Ou",
  "ratio_A_cold_wall_over_csq_Ou", "ratio_A_warm_over_csq_Ou", "gate_target_wall_s"))
  cat("- ", name, ": ", fmt(summary$value[summary$measure == name], 3), "\n", sep = "")
for (g in c("hg002_BW", "hg002_AW")) {
  x <- runs[runs$group == g & !is.na(runs$pass) & runs$pass == "warm", ]
  a <- audit_of(x$audit[nrow(x)])
  cat("- ", g, " audit: ", paste(sprintf("%s=%s", names(a), unlist(a)), collapse = "; "), "\n", sep = "")
}
cat("\n")

cat("### Inside and outside the csq domain (mode A, cold, calls partitioned by the committed accounting)\n\n")
part <- lapply(c("inside", "outside"), function(p) {
  x <- runs[runs$group == paste0("hg002_", p, "_A") & !is.na(runs$pass) & runs$pass == "cold", ]
  a <- audit_of(x$audit[nrow(x)])
  data.frame(partition = p, `records with calls` = a$records_in_scope %||% "-", calls = a$calls %||% "-", `output rows` = a$rows %||% "-", predicted = a$predicted %||% "-",
    `predict + Parquet s (median)` = fmt(med(num(x$predict_s))), `load s` = fmt(med(num(x$load_s))), check.names = FALSE)
})
md(do.call(rbind, part))

cat("### Caps observed (maximum over the three cold processes of each job)\n\n")
owners <- c("model", "index", "workspace", "scratch", "emit", "control", "total")
cap_rows <- lapply(c("hg002_B", "hg002_A", "qual5m_B", "qual5m_A", "dense_B", "low_sharing_B"), function(g) {
  x <- runs[runs$group == g & !is.na(runs$pass) & runs$pass == "cold", ]
  if (!nrow(x)) return(NULL)
  hw <- vapply(owners, function(o) max(num(x[[paste0("pass_hw_", o, "_mib")]]), na.rm = TRUE), 0)
  lw <- vapply(owners, function(o) max(num(x[[paste0("load_hw_", o)]]) / 1048576, na.rm = TRUE), 0)
  data.frame(job = g, `memory.peak GiB (16 cap)` = fmt(max(num(x$cgroup_peak_bytes)) / 2^30), `max RSS GiB` = fmt(max(num(x$max_rss_kib)) / 2^20),
    `native high-water MiB, load (total)` = fmt(lw[["total"]], 0),
    `native high-water MiB per owner, pass` = paste(sprintf("%s %s", owners[-length(owners)], fmt(hw[-length(hw)], 0)), collapse = ", "),
    `native pass total MiB (4096 cap)` = fmt(hw[["total"]], 0), `DuckDB spill peak bytes` = fmt(max(num(x$spill_peak_bytes)), 0), check.names = FALSE)
})
md(do.call(rbind, cap_rows))

cat("### 5M MANE qualification (3 fresh processes per mode)\n\n")
counts <- if (!is.na(raw) && file.exists(file.path(raw, "input_counts.tsv"))) read.delim(file.path(raw, "input_counts.tsv")) else NULL
if (!is.null(counts)) md(counts)
q5 <- lapply(c("qual5m_B", "qual5m_A"), function(g) {
  x <- runs[runs$group == g & !is.na(runs$pass) & runs$pass == "cold", ]
  w <- runs[runs$group == sub("_([AB])$", "_\\1W", g) & !is.na(runs$pass) & runs$pass == "warm", ]
  wa <- runs[runs$group == sub("_([AB])$", "_\\1W", g) & !is.na(runs$pass) & runs$pass == "warm", ]
  a <- audit_of(wa$audit[nrow(wa)])
  data.frame(mode = sub("qual5m_", "", g), `cold wall s (runs)` = paste(fmt(num(x$wall_s)), collapse = " / "), `median` = fmt(med(num(x$wall_s))),
    load = fmt(med(num(x$load_s))), stage = fmt(med(num(x$stage_s))), predict = fmt(med(num(x$predict_s))), `warm pass` = fmt(med(num(w$pass_s))),
    `output checksum (audited runs)` = paste(unique(vapply(wa$audit, function(t) audit_of(t)$checksum %||% "-", "")), collapse = ","), check.names = FALSE)
})
md(do.call(rbind, q5))
for (g in c("qual5m_BW", "qual5m_AW")) {
  x <- runs[runs$group == g & !is.na(runs$pass) & runs$pass == "warm", ]
  a <- audit_of(x$audit[nrow(x)])
  cat("- ", g, " audit: ", paste(sprintf("%s=%s", names(a), unlist(a)), collapse = "; "), "\n", sep = "")
}
cat("\n### Controls (mode B, 3 fresh processes)\n\n")
ctl <- lapply(c("dense_B", "low_sharing_B", "dense_A", "low_sharing_A"), function(g) {
  x <- runs[runs$group == g & !is.na(runs$pass) & runs$pass == "cold", ]
  a <- audit_of(x$audit[nrow(x)])
  data.frame(control = g, `median wall s` = fmt(med(num(x$wall_s))), `predict s` = fmt(med(num(x$predict_s))), calls = a$calls %||% "-", `output rows` = a$rows %||% "-",
    carriers = a$carriers %||% "-", `translated bases` = a$translated_bases %||% "-", `max events per transcript` = a$max_events_per_transcript %||% "-", transcripts = a$transcripts %||% "-",
    `output bytes` = a$output_bytes %||% "-", check.names = FALSE)
})
md(do.call(rbind, ctl))

cat("### Optimizations (mode A cold on the HG002 calls, alternating variants, median of 3)\n\n")
ab <- runs[grepl("^ab_", runs$label) & !is.na(runs$pass) & runs$pass == "cold", ]
if (nrow(ab)) md(do.call(rbind, lapply(split(ab, ab$group), function(x) data.frame(variant = sub("^ab_", "", x$group[1L]),
  `load s` = fmt(med(num(x$load_s))), `predict + Parquet s` = fmt(med(num(x$predict_s))), `process wall s` = fmt(med(num(x$wall_s))),
  `checksum` = paste(unique(vapply(x$audit, function(t) audit_of(t)$checksum %||% "-", "")), collapse = ","), check.names = FALSE))))

cat("### Pipeline choices (cold, alternating, median of 3)\n\n")
abx <- runs[grepl("^abx_", runs$label) & !is.na(runs$pass) & runs$pass == "cold", ]
if (nrow(abx)) md(do.call(rbind, lapply(split(abx, abx$group), function(x) data.frame(choice = sub("^abx_", "", x$group[1L]),
  `process wall s` = fmt(med(num(x$wall_s))), `load s` = fmt(med(num(x$load_s))), `stage s` = fmt(med(num(x$stage_s))),
  `max RSS GiB` = fmt(max(num(x$max_rss_kib)) / 2^20), check.names = FALSE))))
