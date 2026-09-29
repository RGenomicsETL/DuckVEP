#!/usr/bin/env Rscript
# Real BND corpora: mate/event identity and endpoint-gene evidence.
#
#  * FusionCatcher v1.20 test call set (GRCh38, Ensembl 98 gene ids, 17 published
#    fusions as 44 physical BND records) pinned in
#    test/duckvep/conformance/data/breakend_fusion/PROVENANCE.tsv;
#  * GRIDSS DO52605T example excerpt (GRCh37), a spec-conformant caller output
#    with inserted sequence, single breakends and one contradictory mate pair.
# Identity is checked by the native builders; endpoint genes come from
# duckvep_annotate over each physical record (records are never merged). The
# expected gene at each FusionCatcher endpoint is the caller's own GENS5/GENS3
# INFO field. Loads an immutable copy of the extension; writes no receipt.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("r/Rduckvep/R/builders.R")
source("r/Rduckvep/R/structural_identity.R")
data_dir <- "test/duckvep/conformance/data/breakend_fusion"
model <- Sys.getenv("DUCKVEP_GRCH38_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
results_file <- file.path(data_dir, "fusioncatcher_v1.20_endpoint_genes.tsv")

provenance <- read.delim(file.path(data_dir, "PROVENANCE.tsv"), stringsAsFactors = FALSE)
pin <- function(fixture, field) provenance$value[provenance$fixture == fixture & provenance$field == field]
sha256 <- function(path) tolower(strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " ")[[1L]][1L])
fc_path <- file.path(data_dir, "fusioncatcher_v1.20_test.vcf")
gr_path <- file.path(data_dir, "gridss_DO52605T_excerpt.vcf")
stopifnot(identical(sha256(fc_path), pin("fusioncatcher_v1.20_test.vcf", "fixture_sha256")),
          identical(sha256(gr_path), pin("gridss_DO52605T_excerpt.vcf", "fixture_sha256")))

read_vcf <- function(path) {
  lines <- grep("^[^#]", readLines(path), value = TRUE)
  x <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
  data.frame(event_index = seq_len(nrow(x)) - 1L, chrom = x[, 1L], pos = as.numeric(x[, 2L]),
             id = x[, 3L], ref = x[, 4L], alt = x[, 5L], info = x[, 8L], stringsAsFactors = FALSE)
}
directory <- tempfile("bnd-fusion-")
dir.create(directory)
binary <- file.path(directory, "duckvep.duckdb_extension")
stopifnot(file.copy(Sys.getenv("DUCKVEP_EXTENSION_FILE", "build/release/duckvep.duckdb_extension"), binary))
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
on.exit({dbDisconnect(con, shutdown = TRUE); unlink(directory, recursive = TRUE)}, add = TRUE)
q <- function(value) as.character(dbQuoteString(con, value))
invisible(dbExecute(con, paste("LOAD", q(binary))))
info_field <- function(info, key)
  vapply(strsplit(info, ";", fixed = TRUE), function(tokens) {
    hit <- grep(paste0("^", key, "="), tokens, value = TRUE)
    if (length(hit)) sub(paste0("^", key, "="), "", hit[1L]) else NA_character_
  }, "")

## ---- GRIDSS excerpt: spec-conformant pairs, insertions, a contradictory pair ----
gridss <- read_vcf(gr_path)
dbWriteTable(con, "gridss", gridss, temporary = TRUE)
g <- rduckvep_prepare_breakend_pairs(con, "gridss")
stopifnot(nrow(g) == nrow(gridss), all(g$event_index == gridss$event_index))
gr_reason <- table(g$reason)
stopifnot(identical(as.integer(gr_reason[c("reciprocal", "mate_coordinate_conflict", "single_breakend")]),
                    c(6L, 2L, 3L)), sum(gr_reason) == 11L)
insertion_pair <- g[g$id %in% c("gridss28fb_8844o", "gridss28fb_8844h"), ]
stopifnot(all(insertion_pair$status == "reciprocal"), all(insertion_pair$inserted_length == 28L),
          all(!is.na(insertion_pair$inserted_sequence)),
          all(g$fusion_status == "not_asserted"), all(g$phase_status == "not_evaluated"))
contradictory <- g[g$reason == "mate_coordinate_conflict", ]
stopifnot(all(contradictory$id_reciprocal), !any(contradictory$coordinate_reciprocal),
          all(is.na(contradictory$pair_key)))
# Both physical records of a pair stay separate rows with the same pair key.
stopifnot(all(table(g$pair_key[!is.na(g$pair_key)]) == 2L))

## ---- FusionCatcher: real known fusions -----------------------------------------
fc <- read_vcf(fc_path)
stopifnot(nrow(fc) == 44L)
dbWriteTable(con, "fc", fc, temporary = TRUE)
p <- rduckvep_prepare_breakend_pairs(con, "fc")
stopifnot(nrow(p) == 44L, all(p$event_index == fc$event_index),
          all(p$record_kind == "paired_breakend"),
          all(p$id_reciprocal), all(p$coordinate_reciprocal), all(p$event_agree),
          all(p$fusion_status == "not_asserted"))
# The caller writes 5' and 3' records with the same bracket family. By the VCF
# BND mate rule that is a contradiction for every pair except one.
identity_reasons <- table(p$reason)
stopifnot(identical(as.integer(identity_reasons[c("mate_orientation_conflict", "reciprocal")]),
                    c(42L, 2L)))
reciprocal_pair <- p[p$status == "reciprocal", ]
stopifnot(identical(sort(reciprocal_pair$id), c("BRD4--NUTM1__3", "BRD4--NUTM1__5")),
          length(unique(reciprocal_pair$pair_key)) == 1L)
dbWriteTable(con, "fc_pairs", p, temporary = TRUE)

invisible(dbExecute(con, paste("ATTACH", q(model), "AS m (READ_ONLY)")))
queries <- c(
  "SELECT seq_region, sequence_length FROM m.bench_regions ORDER BY seq_region",
  paste("SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
    "codon_table, pre_cds_sequence, post_cds_sequence FROM m.bench_transcripts",
    "ORDER BY transcript_index"),
  paste("SELECT transcript_index, exon_start, exon_end, exon_cdna_start,",
    "exon_cdna_end, phase, end_phase FROM m.bench_exons",
    "ORDER BY transcript_index, exon_cdna_start"))
load <- paste0("SELECT loaded FROM duckvep_model_load('breakend_fusion',",
  paste(vapply(queries, q, ""), collapse = ","), ")")
stopifnot(isTRUE(dbGetQuery(con, load)$loaded))
invisible(dbExecute(con, paste(
  "CREATE TEMP TABLE fc_events AS SELECT f.event_index::UBIGINT event_index,",
  "r.seq_region::UINTEGER seq_region, f.pos::UBIGINT AS \"position\", NULL::VARCHAR reference,",
  "f.alt alternate, NULL::UBIGINT end_position, NULL::VARCHAR structural_type,",
  "NULL::VARCHAR copy_change, r2.seq_region::UINTEGER mate_seq_region,",
  "try_cast(regexp_extract(f.alt, ':([0-9]+)', 1) AS UBIGINT) mate_position",
  "FROM fc f JOIN m.bench_regions r ON r.chrom = f.chrom",
  "JOIN m.bench_regions r2 ON r2.chrom = regexp_extract(f.alt, '[][]([^:]+):', 1)",
  "WHERE f.pos <= r.sequence_length AND",
  "try_cast(regexp_extract(f.alt, ':([0-9]+)', 1) AS UBIGINT) <= r2.sequence_length")))
# The caller wrote chr14:190287409 (a chr4-scale coordinate; chr14 is 107,043,718 bp)
# for DUX4--IGH@ event 17. DuckVEP refuses such endpoints, so both physical records
# are excluded from annotation and reported as unmatched, not repaired.
skipped <- setdiff(fc$id, dbGetQuery(con, "SELECT f.id FROM fc f JOIN fc_events e USING (event_index)")$id)
stopifnot(identical(sort(skipped), c("DUX4--IGH@__3__17", "DUX4--IGH@__5__17")))
invisible(dbExecute(con, paste(
  "CREATE TEMP TABLE fc_genes AS SELECT DISTINCT a.event_index, t.gene_stable_id AS gene_id",
  "FROM query(duckvep_annotate_sql('fc_events', 'breakend_fusion', struct_pack(rich := true,",
  "upstream_distance := 5000, downstream_distance := 5000))) a",
  "JOIN m.model_transcripts t USING (transcript_index)",
  "WHERE a.region IS NOT NULL AND a.region NOT IN ('upstream', 'downstream')")))
f <- rduckvep_prepare_breakend_fusion(con, "fc_pairs", "fc_genes")
stopifnot(nrow(f) == 44L, all(f$event_index == fc$event_index), all(!f$fusion_asserted),
          all(f$phase_status == "unproven"))
genes <- dbGetQuery(con, "SELECT event_index, list(gene_id ORDER BY gene_id) AS genes FROM fc_genes GROUP BY event_index")
endpoint_genes <- setNames(genes$genes, genes$event_index)
fc$end <- info_field(fc$info, "END")
fc$expected_gene <- ifelse(fc$end == "5", info_field(fc$info, "GENS5"), info_field(fc$info, "GENS3"))
fc$fusion <- info_field(fc$info, "FUSION")
found <- vapply(seq_len(nrow(fc)), function(i)
  fc$expected_gene[i] %in% endpoint_genes[[as.character(fc$event_index[i])]], TRUE)
composite <- fc$expected_gene == "ENSG09000000014" | fc$expected_gene == "ENSG09000001017"
result <- data.frame(id = fc$id, fusion = fc$fusion, end = fc$end, expected_gene = fc$expected_gene,
  endpoint_genes = vapply(seq_len(nrow(fc)), function(i)
    paste(endpoint_genes[[as.character(fc$event_index[i])]], collapse = ","), ""),
  expected_gene_at_endpoint = found, identity_reason = p$reason, fusion_status = f$status,
  stringsAsFactors = FALSE)
write.table(result, results_file, sep = "\t", quote = FALSE, row.names = FALSE)
# IGH@ is a FusionCatcher locus label, not an Ensembl gene: its ids cannot be found
# in the Ensembl model and are retained as unmatched, not scored.
stopifnot(all(!found[composite]), all(!found[fc$id %in% skipped]))
scored <- !composite & !(fc$id %in% skipped)
cat(sprintf("FusionCatcher: %d records, %d scored endpoints, %d expected genes found; %d IGH@ composite ids and %d out-of-contig records unmatched\n",
            nrow(fc), sum(scored), sum(found[scored]), sum(composite), length(skipped)))
# Four scored records do not overlap their declared gene: three CRLF2 records whose
# chrX coordinates the caller wrote away from CRLF2 (105.8 Mb; 1.2287 Mb, 16 kb past
# its end), and the MALT1 5' record 1.5 kb upstream of MALT1. Retained, not repaired.
stopifnot(identical(sort(fc$id[scored & !found]),
  c("IGH@--CRLF2__3__5", "IGH@--CRLF2__3__6", "IGH@--CRLF2__3__7", "MALT1--IGH@__5")),
  sum(found[scored]) == 33L, sum(scored) == 37L)
status_counts <- table(f$status)
cat("fusion builder statuses:", paste(names(status_counts), status_counts, collapse = ", "), "\n")
# Only the one orientation-consistent pair can reach candidate_partner_genes.
stopifnot(all(f$status[f$status == "candidate_partner_genes"] == "candidate_partner_genes"),
          !any(f$status[p$reason == "mate_orientation_conflict"] == "candidate_partner_genes"),
          all(f$status[p$reason == "mate_orientation_conflict"] %in%
                c("candidate_orientation_conflict", "endpoint_without_gene", "shared_gene_endpoints")))

## ---- Failure controls on the real records ---------------------------------------
control <- function(rows) {
  dbWriteTable(con, "bnd_control", rows, temporary = TRUE, overwrite = TRUE)
  rduckvep_prepare_breakend_pairs(con, "bnd_control")
}
base <- fc[fc$fusion == "BRD4--NUTM1", c("event_index", "chrom", "pos", "id", "ref", "alt", "info")]
base$event_index <- 0:1
stopifnot(all(control(base)$reason == "reciprocal"))
missing_mate <- control(base[1L, ])
stopifnot(identical(missing_mate$reason, "mate_not_found"))
duplicated_mate <- control(rbind(base, transform(base[2L, ], event_index = 2L)))
stopifnot(all(duplicated_mate$reason[1:3] == c("duplicate_id", "duplicate_id", "duplicate_id")) ||
          all(duplicated_mate$status == "conflict"))
shifted_mate <- base
shifted_mate$pos[2L] <- shifted_mate$pos[2L] + 1
stopifnot(identical(control(shifted_mate)$reason[1L], "mate_coordinate_conflict"))
wrong_event <- base
wrong_event$info[2L] <- sub("EVENT=[^;]*", "EVENT=OTHER", wrong_event$info[2L])
stopifnot(identical(control(wrong_event)$reason[1L], "event_conflict"))
no_mateid <- base
no_mateid$info[1L] <- gsub("MATEID=[^;]*;", "", no_mateid$info[1L])
stopifnot(identical(control(no_mateid)$reason, c("missing_mateid", "mate_missing_mateid")))
cat("BND controls: 6 mutations of the BRD4--NUTM1 pair each fail with the expected stable reason\n")
