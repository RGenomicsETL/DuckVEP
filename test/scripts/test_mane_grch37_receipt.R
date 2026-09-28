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

# Named witnesses: specific rows of the real release, so a regression in a
# known case is caught by name, not only by status counts. Chosen
# deterministically as the lowest refseq_nuc (plain lexicographic order,
# identical between DuckDB's default collation and data.table's radix sort)
# within each category below.
witness <- function(refseq, status, enst, enst_version, index, strand) {
  row <- relation[refseq_nuc == refseq]
  stopifnot(nrow(row) == 1L,
    identical(row$mapping_status, status),
    identical(row$candidate_enst, enst),
    identical(row$candidate_enst_version, enst_version),
    identical(row$transcript_index, index),
    identical(row$target_strand, strand))
}

# MANE Select, exact_model_match, one witness per strand.
witness("NM_000394.4", "exact_model_match", "ENST00000291554", 2L, 13812L, "+")
witness("NM_000079.4", "exact_model_match", "ENST00000348749", 5L, 33707L, "-")
# MANE Select, cds_exact_utr_differs, one witness per strand.
witness("NM_000015.3", "cds_exact_utr_differs", "ENST00000286479", 3L, 152504L, "+")
witness("NM_000014.6", "cds_exact_utr_differs", "ENST00000318602", 7L, 122870L, "-")
# MANE Plus Clinical, accepted in either tier (lowest refseq_nuc across both
# exact_model_match and cds_exact_utr_differs Plus Clinical rows).
witness("NM_000248.4", "cds_exact_utr_differs", "ENST00000394351", 3L, 108374L, "+")
stopifnot(relation[refseq_nuc == "NM_000248.4"]$mane_status == "MANE Plus Clinical")
# One witness per rejection status (every mapping_status other than the two
# accepted tiers).
witness("NM_000451.4", "ambiguous_target_locus", NA_character_, NA_integer_, NA_integer_, "+")
witness("NM_001004457.2", "cds_phase_mismatch", "ENST00000373688", 2L, NA_integer_, "+")
witness("NM_000036.3", "geometry_mismatch", "ENST00000520113", 2L, NA_integer_, "-")
witness("NM_000026.4", "refseq_only_no_gencode19_match", NA_character_, NA_integer_, NA_integer_, "+")
witness("NR_002728.4", "sequence_mismatch", "ENST00000597346", 1L, NA_integer_, "-")
witness("NM_001005513.1", "target_reference_unavailable", NA_character_, NA_integer_, NA_integer_, "+")
witness("NM_000131.5", "target_transcript_absent", "ENST00000375581", 3L, NA_integer_, NA_character_)
witness("NM_000792.7", "translation_mismatch", "ENST00000361921", 3L, NA_integer_, "+")
stopifnot(setequal(relation[!mapping_status %in% c("exact_model_match", "cds_exact_utr_differs")]$mapping_status,
  c("ambiguous_target_locus", "cds_phase_mismatch", "geometry_mismatch",
    "refseq_only_no_gencode19_match", "sequence_mismatch", "target_reference_unavailable",
    "target_transcript_absent", "translation_mismatch")))

# MANE v1.5 has no mitochondrial gene records at all (no source row with
# GRCh38_chr chrMT or a symbol prefixed MT-), so this release has no row
# placed on the RefSeq mitochondrial contig NC_012920 either. There is no
# witness row to name for this category; this assertion documents the
# absence directly from the release output instead.
stopifnot(sum(grepl("^NC_012920", relation$target_sequence_accession)) == 0L)

cat("MANE GRCh37 full-release receipt and Parquet: passed\n")
print(counts)
