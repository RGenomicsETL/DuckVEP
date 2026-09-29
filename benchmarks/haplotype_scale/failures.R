#!/usr/bin/env Rscript
# Explicit capacity failures on the haplotype path (#2 slice 7): each must raise an explicit capacity error, publish
# nothing, return the native budget to its baseline and leave the connection and the loaded model usable.
#
#   Rscript benchmarks/haplotype_scale/failures.R EXTENSION CALLS.parquet OUT_DIR RECEIPT.tsv
#
# CALLS.parquet is the ordered HG002 calls relation written by worker.R --mode stage (full Ensembl 116 model).
# Controls, in one process so that reuse after a failure is real:
#   model_load_budget     native budget 64 MiB before duckvep_model_load
#   haplotype_budget      budget = resident model + 64 MiB before duckvep_haplotypes (workspace allocation)
#   alignment_capacity    default max_alignment_cells (16,777,216) on the calls of the longest transcripts
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: Rscript failures.R EXTENSION CALLS.parquet OUT_DIR RECEIPT.tsv")
extension <- normalizePath(args[[1L]], mustWork = TRUE); calls <- normalizePath(args[[2L]], mustWork = TRUE)
out_dir <- args[[3L]]; receipt <- args[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model_path <- "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb"
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true", threads = "1")))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
q <- function(x) as.character(dbQuoteString(con, x))
sql <- function(s) invisible(dbExecute(con, s))
get <- function(s) dbGetQuery(con, s)
lines <- character()
note <- function(control, key, value) {
  value <- gsub("[\t\n\r]+", " ", as.character(value))
  lines[[length(lines) + 1L]] <<- paste(control, key, value, sep = "\t")
  cat(control, key, value, "\n", sep = "\t")
}
attempt <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
# R releases a failed statement's bind state at its next garbage collection, so collect before reading the budget.
total_bytes <- function() { invisible(gc()); invisible(gc()); get("SELECT current_bytes AS b FROM duckvep_native_budget() WHERE owner = 'total'")$b }
# COPY ... TO creates its target when it starts. The runner protocol writes <out>.partial and renames it after success, so
# a failed statement publishes nothing at the final name; the leftover partial file is reported and removed.
copy_to <- function(query, target) {
  partial <- paste0(target, ".partial")
  message <- attempt(sql(sprintf("COPY (%s) TO %s (FORMAT parquet)", query, q(partial))))
  if (is.na(message)) file.rename(partial, target)
  list(message = message, partial_bytes = if (file.exists(partial)) file.size(partial) else 0, partial = partial)
}
sql(paste("LOAD", q(extension)))
sql("SET memory_limit = '8GB'")
sql(paste0("ATTACH ", q(model_path), " AS m (READ_ONLY)"))
sql("CREATE TABLE regions AS SELECT seq_region::BIGINT AS seq_region, seq_region_name, sequence_length FROM m.model_regions")
sql(paste0("CREATE TABLE calls AS SELECT * FROM read_parquet(", q(calls), ")"))
load_sql <- "SELECT loaded FROM duckvep_model_load('hap', 'SELECT seq_region::UINTEGER AS seq_region, sequence_length FROM regions ORDER BY seq_region',
 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM m.model_transcripts',
 'SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase FROM m.model_transcripts, unnest(exons) u(e)',
 mature_mirna_query := 'SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM m.model_transcripts, unnest(mature_mirna_regions) u(x)',
 peptide_edit_query := 'SELECT transcript_index, x.protein_position, x.alternate_amino_acid FROM m.model_transcripts, unnest(peptide_edits) u(x)',
 transcript_coverage_complete := TRUE)"
predict_sql <- function(cells = 268435456, workspace = 1073741824, source = "SELECT * FROM calls")
  sprintf("SELECT * FROM duckvep_haplotypes(%s, 'hap', max_alignment_cells := %.0f, workspace_limit := %.0f)", q(source), cells, workspace)
small_query <- "SELECT count(*) AS n FROM duckvep_haplotypes('SELECT * FROM calls WHERE event_index < 0', 'hap')"

# 1. Model load under a 64 MiB budget.
sql("SELECT duckvep_native_budget_set(67108864)")
message1 <- attempt(get(load_sql))
note("model_load_budget", "error", message1)
note("model_load_budget", "explicit_capacity_error", grepl("capacity error", message1, fixed = TRUE))
note("model_load_budget", "bytes_charged_after", total_bytes())
note("model_load_budget", "model_published", is.na(attempt(get(small_query))) )
sql("SELECT duckvep_native_budget_set(4294967296)")
note("model_load_budget", "reload_after_failure", isTRUE(get(load_sql)$loaded))
invisible(get(small_query))
resident <- total_bytes()
note("model_load_budget", "resident_bytes_after_reload", resident)

# 2. Haplotype builder under resident + 64 MiB: the workspace allocation is refused.
target <- file.path(out_dir, "must_not_exist.parquet")
if (file.exists(target)) file.remove(target)
sql(sprintf("SELECT duckvep_native_budget_set(%.0f)", resident + 67108864))
r2 <- copy_to(predict_sql(), target); message2 <- r2$message
note("haplotype_budget", "partial_file_bytes_left_by_copy", r2$partial_bytes); if (file.exists(r2$partial)) invisible(file.remove(r2$partial))
note("haplotype_budget", "error", message2)
note("haplotype_budget", "explicit_capacity_error", grepl("capacity error|budget exceeded|workspace", message2))
note("haplotype_budget", "output_published", file.exists(target))
message2b <- attempt(sql(paste("CREATE TABLE must_not_exist AS", predict_sql())))
note("haplotype_budget", "table_published", nrow(get("SELECT 1 FROM duckdb_tables() WHERE table_name = 'must_not_exist'")) > 0)
note("haplotype_budget", "bytes_charged_after", total_bytes())
note("haplotype_budget", "bytes_equal_resident", total_bytes() == resident)
sql("SELECT duckvep_native_budget_set(4294967296)")
note("haplotype_budget", "connection_reusable", identical(get(small_query)$n, 0))

# 3. A small alignment capacity refuses, explicitly, the long coding sequences of the calls; the builder's own default (16,777,216
# cells) no longer refuses them since capacity follows the band an alignment needs, so it is recorded as a success.
target3 <- file.path(out_dir, "must_not_exist_alignment.parquet")
if (file.exists(target3)) file.remove(target3)
r3 <- copy_to(predict_sql(cells = 100000, workspace = 268435456), target3); message3 <- r3$message
note("alignment_capacity", "cells", 100000)
note("alignment_capacity", "error", message3)
note("alignment_capacity", "explicit_error", !is.na(message3))
note("alignment_capacity", "names_limit", grepl("max_alignment_cells", message3, fixed = TRUE))
note("alignment_capacity", "partial_file_bytes_left_by_copy", r3$partial_bytes); if (file.exists(r3$partial)) invisible(file.remove(r3$partial))
note("alignment_capacity", "output_published", file.exists(target3))
note("alignment_capacity", "bytes_equal_resident", total_bytes() == resident)
note("alignment_capacity", "connection_reusable", identical(get(small_query)$n, 0))
target4 <- file.path(out_dir, "builder_defaults.parquet")
if (file.exists(target4)) file.remove(target4)
r4 <- copy_to(sprintf("SELECT * FROM duckvep_haplotypes(%s, 'hap')", q("SELECT * FROM calls")), target4)
note("builder_defaults", "succeeds_with_max_alignment_cells_16777216", is.na(r4$message) && file.exists(target4))
note("builder_defaults", "output_rows", if (file.exists(target4)) get(sprintf("SELECT count(*) AS n FROM read_parquet(%s)", q(target4)))$n else NA)
if (file.exists(target4)) file.remove(target4)
writeLines(c("control\tkey\tvalue", lines), receipt)
