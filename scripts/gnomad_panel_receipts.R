#!/usr/bin/env Rscript
# Write the committed receipts for the panels under $DUCKVEP_GNOMAD_ROOT/panels and refresh the
# stable links the scale runner resolves (--panel 1M, 5M, 25M, 100M, exome2M, structural).
#
#   Rscript scripts/gnomad_panel_receipts.R GENOME_PANEL_DIR EXOME_PANEL_DIR [STRUCTURAL_PARQUET]
#
# Panels are data and are never committed; benchmarks/data/scale_contracts/panels/ holds only
#   panel-receipts.tsv   one row per panel: rows, file SHA-256, content checksum, status, sources
#   panel-shards.tsv     per source shard: set, chromosome, parts, bytes, alleles, manifest SHA-256
# A panel that cannot be built yet (too few distinct alleles until every genome shard lands) is
# recorded with its status and no checksum, and its link is removed.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("usage: Rscript scripts/gnomad_panel_receipts.R GENOME_DIR EXOME_DIR [STRUCTURAL_PARQUET]")
root <- Sys.getenv("DUCKVEP_GNOMAD_ROOT", "/root/duckvep/data/gnomad-v4.1")
panels_root <- file.path(root, "panels")
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), ".."))
out_dir <- file.path(repo, "benchmarks/data/scale_contracts/panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
genome_dir <- normalizePath(args[[1L]]); exome_dir <- normalizePath(args[[2L]])
structural <- if (length(args) >= 3L) normalizePath(args[[3L]]) else file.path(root, "structural-controls.parquet")
sha256 <- function(path) substr(system2("sha256sum", path, stdout = TRUE), 1L, 64L)
con <- dbConnect(duckdb(shared_home = FALSE))
q <- function(x) as.character(dbQuoteString(con, x))
content <- function(path) {
  r <- dbGetQuery(con, paste0("SELECT count(*)::VARCHAR AS rows, sum(hash(t)::HUGEINT)::VARCHAR AS hash_sum, ",
    "bit_xor(hash(t))::VARCHAR AS hash_xor FROM read_parquet(", q(path), ") t"))
  r
}
sources_of <- function(dir) read.delim(file.path(dir, "sources.tsv"), colClasses = "character")
shard_of <- function(files) sub("^.*/((genomes|exomes)-chr[0-9XY]+)/part-.*$", "\\1", files)
description <- function(dir) {
  s <- sources_of(dir)
  shards <- unique(shard_of(s$file))
  list(sets = paste(sort(unique(sub("-chr.*$", "", shards))), collapse = ","),
    chromosomes = paste(sub("^.*-chr", "", shards)[order(suppressWarnings(as.integer(sub("^.*-chr", "", shards))),
      sub("^.*-chr", "", shards))], collapse = ","),
    parts = nrow(s), sources_sha256 = sha256(file.path(dir, "sources.tsv")))
}
git_revision <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) "unknown")
duckdb_version <- dbGetQuery(con, "SELECT version() AS v")$v
rows <- list()
add <- function(panel, kind, quota, dir, status, path, d, extra = list()) {
  c <- if (status == "complete") content(path) else list(rows = "0", hash_sum = NA, hash_xor = NA)
  rows[[length(rows) + 1L]] <<- data.frame(panel = panel, kind = kind, quota = quota, status = status,
    rows = c$rows, file_sha256 = if (status == "complete") sha256(path) else NA_character_,
    hash_sum = c$hash_sum, hash_xor = c$hash_xor, sets = d$sets, chromosomes = d$chromosomes,
    source_parts = d$parts, sources_sha256 = d$sources_sha256, artifact = basename(dir),
    builder_revision = git_revision, duckdb = duckdb_version, stringsAsFactors = FALSE)
}
link <- function(name, target) {
  path <- file.path(panels_root, name)
  unlink(path)
  if (!is.null(target)) invisible(file.symlink(target, path))
}

g <- description(genome_dir)
gp <- read.delim(file.path(genome_dir, "panels.tsv"), colClasses = "character")
counts <- read.delim(file.path(genome_dir, "counts.tsv"), colClasses = "character")
label <- c("1000000" = "1M", "5000000" = "5M", "25000000" = "25M", "100000000" = "100M")
for (i in seq_len(nrow(gp))) {
  path <- file.path(genome_dir, paste0("panel-", gp$quota[i], ".parquet"))
  add(paste0("genomes-", label[[gp$quota[i]]]), "genome distinct alleles", gp$quota[i], genome_dir, gp$status[i], path, g)
  link(paste0("genomes-", label[[gp$quota[i]]], ".parquet"), if (gp$status[i] == "complete") path else NULL)
}
e <- description(exome_dir)
ep <- read.delim(file.path(exome_dir, "panel.tsv"), colClasses = "character")
epath <- file.path(exome_dir, "exome-2m.parquet")
add("exomes-2M", "exome coding/splice enriched", "2000000", exome_dir, ep$status[1L], epath, e)
link("exomes-2M.parquet", if (ep$status[1L] == "complete") epath else NULL)
if (file.exists(structural)) {
  sd <- list(sets = "synthetic", chromosomes = "1-22,X,Y", parts = 1L, sources_sha256 = NA_character_)
  add("structural-controls", "synthetic structural controls", "12000", dirname(structural), "complete", structural, sd)
  link("structural-controls.parquet", structural)
}
receipts <- do.call(rbind, rows)
avail_path <- file.path(exome_dir, "availability.tsv")
bounds_path <- file.path(exome_dir, "bounds.tsv")
receipts$detail <- ""
if (file.exists(avail_path)) {
  a <- read.delim(avail_path, colClasses = "character")
  b <- read.delim(bounds_path, colClasses = "character")
  receipts$detail[receipts$panel == "exomes-2M"] <- paste0(
    paste0(a$bin, ": required ", a$required, ", available ", a$available, collapse = "; "),
    "; distinct literal alleles by class (snv, indel, mnv, all): ", paste(b$source_upper_bound, collapse = ","),
    "; gnomAD sites contain no MNVs, so the 100k CDS MNV quota moved to cds_indel and MNV coverage comes from the ClinVar and indel/MNV controls")
}
is_genome <- startsWith(receipts$panel, "genomes-")
receipts$source_alleles <- ifelse(is_genome, counts$source_alts, "")
receipts$distinct_alleles <- ifelse(is_genome, counts$distinct_alts, "")
receipts[] <- lapply(receipts, function(x) { x <- as.character(x); x[is.na(x) | x == ""] <- "-"; x })
write.table(receipts, file.path(out_dir, "panel-receipts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# One row per source shard that fed a panel.
shard_rows <- do.call(rbind, lapply(c(genome_dir, exome_dir), function(dir) {
  s <- sources_of(dir)
  s$shard <- shard_of(s$file)
  do.call(rbind, lapply(split(s, s$shard), function(x) {
    shard_dir <- file.path(root, x$shard[1L])
    cnt <- read.delim(file.path(shard_dir, "counts.tsv"), colClasses = "character")
    manifest <- tempfile()
    write.table(x[order(x$file), c("file", "bytes", "sha256")], manifest, sep = "\t", row.names = FALSE, quote = FALSE)
    data.frame(panel_artifact = basename(dir), shard = x$shard[1L], parts = nrow(x),
      part_bytes = sum(as.numeric(x$bytes)), alleles = cnt$alts, literal = cnt$literal,
      manifest_sha256 = sha256(manifest), stringsAsFactors = FALSE)
  }))
}))
write.table(shard_rows, file.path(out_dir, "panel-shards.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
cat(file.path(out_dir, "panel-receipts.tsv"), "\n")
