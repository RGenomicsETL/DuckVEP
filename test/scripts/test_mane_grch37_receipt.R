#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DBI); library(data.table); library(digest) })
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("usage: Rscript test/scripts/test_mane_grch37_receipt.R OUTPUT_DIR RECEIPT.csv")
out <- args[1]
receipt <- fread(args[2])
stopifnot(nrow(receipt) >= 1L)
receipt <- tail(receipt, 1L)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = ":memory:")
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
relation <- as.data.table(dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
  dbQuoteString(con, file.path(out, "mane_grch37_mapping.parquet")), ") ORDER BY mane_row")))
native <- as.data.table(dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
  dbQuoteString(con, file.path(out, "grch37_transcript_authorities.parquet")), ")")))
canonical <- tempfile("mane-receipt-")
on.exit(unlink(canonical), add = TRUE)
fwrite(relation, canonical, sep = "\t", na = "\\N", quote = TRUE)
counts <- relation[, .(n = .N), by = mapping_status][order(mapping_status)]
status_counts <- paste(paste(counts$mapping_status, counts$n, sep = ":"), collapse = ";")
accepted <- relation[mapping_status == "exact_model_match"]
coding <- relation[mapping_status == "cds_exact_utr_differs"]
joinable <- relation[!is.na(transcript_index)]
# Each admitted association retains its native transcript authorities.
matched_native <- merge(joinable[, .(transcript_index, canonical, gencode_basic)],
  native[, .(transcript_index, canonical, gencode_basic)], by = "transcript_index",
  suffixes = c("_mapped", "_native"))
stopifnot(nrow(relation) == receipt$row_count, nrow(relation) == 19437L,
  !anyDuplicated(relation$mane_row), !anyNA(relation$mapping_status),
  relation[grepl("^NM_", refseq_nuc) & exact_refseq, .N] == 19306L,
  relation[grepl("^NR_", refseq_nuc) & exact_refseq, .N] == 61L,
  sum(relation$exact_refseq) == receipt$exact_refseq_count,
  identical(status_counts, receipt$status_counts),
  digest(canonical, algo = "sha256", file = TRUE) == receipt$relation_sha256,
  all(relation$model_sha256 == receipt$model_sha256),
  all(relation$source_manifest_sha256 == receipt$source_manifest_sha256),
  all(relation$reference_sha256 == receipt$reference_sha256),
  all(native$model_sha256 == receipt$model_sha256),
  nrow(native) == 195379L,
  uniqueN(native$source_gene_id) == receipt$native_gene_count,
  sum(native$canonical) == receipt$retained_canonical_gene_count,
  sum(native$gencode_basic) == receipt$gencode_basic_transcript_count,
  uniqueN(native[no_retained_canonical == TRUE]$source_gene_id) == receipt$no_retained_canonical_gene_count,
  !any(native$native_mane), !any(native$gencode_primary),
  all(is.na(relation$transcript_index) ==
      !relation$mapping_status %in% c("exact_model_match", "cds_exact_utr_differs")),
  all(accepted$mapping_label == "MANE mapped to GRCh37"),
  all(accepted$exact_refseq & accepted$exon_chain_match & accepted$utr_exon_chain_match &
      accepted$cds_phase_match & accepted$spliced_sequence_match &
      accepted$translated_sequence_match),
  nrow(accepted) == 740L,
  nrow(coding) > 0L,
  all(coding$mapping_label == "MANE mapped to GRCh37 (coding region only)"),
  !anyNA(coding$refseq_rna_sha256), !anyNA(coding$refseq_protein_sha256),
  all(grepl("^(NM|NR)_[0-9]+\\.[0-9]+$", coding$refseq_nuc)),
  all(coding$exact_refseq & !coding$utr_exon_chain_match & coding$cds_phase_match &
      coding$translated_sequence_match & coding$cds_reference_difference_bases == 0L),
  !anyNA(coding[, .(exact_refseq, exon_chain_match, utr_exon_chain_match,
                    cds_phase_match, spliced_sequence_match, translated_sequence_match,
                    cds_reference_difference_bases)]),
  nrow(matched_native) == nrow(joinable),
  all(matched_native$canonical_mapped == matched_native$canonical_native),
  all(matched_native$gencode_basic_mapped == matched_native$gencode_basic_native),
  all(c("+", "-") %in% accepted[mane_status == "MANE Select"]$target_strand),
  all(c("+", "-") %in% accepted[mane_status == "MANE Plus Clinical"]$target_strand))
cat("MANE GRCh37 full-release receipt and Parquet: passed\n")
print(counts)
