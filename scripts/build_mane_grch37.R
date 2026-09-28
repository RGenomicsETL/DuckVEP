#!/usr/bin/env Rscript
# Build the cold MANE v1.5 -> GRCh37.p13 association; never load it into the consequence model.
suppressPackageStartupMessages({
  library(DBI)
  library(data.table)
  library(Biostrings)
  library(Rsamtools)
  library(GenomicRanges)
  library(digest)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script_file), "mane_grch37_policy.R"))
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("usage: Rscript scripts/build_mane_grch37.R STAGING_DIR MODEL.duckdb OUTPUT_DIR")
stage <- normalizePath(args[1]); model_path <- normalizePath(args[2]); out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- c(mane = "MANE.GRCh38.v1.5.summary.txt.gz",
           report = "GCF_000001405.25_GRCh37.p13_assembly_report.txt",
           gff = "GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
           rna = "GCF_000001405.25_GRCh37.p13_rna.fna.gz",
           protein = "GCF_000001405.25_GRCh37.p13_protein.faa.gz",
           reference = "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
           reference_fasta = "Homo_sapiens.GRCh37.dna.primary_assembly.fa")
paths <- setNames(file.path(stage, files), names(files))
if (any(!file.exists(paths))) stop("missing staging inputs: ", paste(files[!file.exists(paths)], collapse = ", "))
sha <- setNames(vapply(paths, digest, character(1), algo = "sha256", file = TRUE), names(files))
expected <- c(mane = "d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c",
  report = "a6cf8300aa2cef9188590bad2d9d54a5909f6d4d2da3f22aa5c5ba2fda1adab3",
  gff = "5fcadac26be5d82a1f1c52e33cc5047247f719a60b50182c6ba298bda77cc80f",
  rna = "0f48eeebe6e5631ff8b346e41db8a9e0d7a3506e1acc9df3fcdcdac8139beb2b",
  protein = "7697d816c1a1639ed304e0fc0fff6867d55df91a5e50f6b6f7287055f389f514",
  reference = "0a43b56dec40debae976d6e70cac68ea6ed874f9fb7c8c814363702ff1d47865",
  reference_fasta = "3a3872e7bdd1532fdbfbc3afd70c04ff8a12a109319ae050c91d962f32e16a4f")
if (!identical(sha, expected)) stop("staged source digest mismatch: ", paste(names(sha)[sha != expected], collapse = ", "))
# The core is copied from Ensembl's public MySQL homo_sapiens_core_116_37
# (ensembldb.ensembl.org:3337, anonymous), as in
# test/scripts/prepare_duckvep_ensembl_fixture.sql. The model's
# source_manifest_sha256 is carried forward: the complete July-definition
# model hash was independently reproduced on this rebuilt model, while the
# current hash excludes the retired post_cds_bases field. The manifest cannot
# be independently re-hashed from a model after the source dumps are removed.
con <- dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = model_path, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
receipt <- as.data.table(dbGetQuery(con, "SELECT * FROM model_receipt"))
stopifnot(nrow(receipt) == 1L, receipt$assembly == "GRCh37", receipt$reference_sha256 == sha["reference"])
model_hash <- receipt$model_sha256
if (!identical(model_hash, "21e113d9148132491bc935f3d1b0ec7d50663f450b62e346cb1f0f447de0b290")) stop("unexpected model hash: ", model_hash)

mane <- fread(cmd = paste("gzip -dc", shQuote(paths["mane"])), sep = "\t", header = TRUE, check.names = FALSE)
setnames(mane, "#NCBI_GeneID", "NCBI_GeneID")
stopifnot(nrow(mane) == 19437L, !anyDuplicated(mane$RefSeq_nuc),
  all(grepl("^(NM|NR)_[0-9]+\\.[0-9]+$", mane$RefSeq_nuc)))
mane[, `:=`(mane_row = .I, enst_root = sub("\\.[0-9]+$", "", Ensembl_nuc))]

report <- fread(paths["report"], skip = "# Sequence-Name", sep = "\t")
setnames(report, sub("^# ", "", names(report)))
report <- report[!is.na(`RefSeq-Accn`) & `RefSeq-Accn` != "na",
                 .(target_sequence_accession = `RefSeq-Accn`,
                   model_contig = `Sequence-Name`, assembly_unit = `Assembly-Unit`,
                   sequence_length = as.integer(`Sequence-Length`))]
stopifnot(!anyDuplicated(report$target_sequence_accession))

# Extract versioned release accessions with their RNA, exon and CDS rows.
# A stable Ensembl root cannot create a missing RefSeq target.
wanted <- unique(mane$RefSeq_nuc)
con_gff <- gzfile(paths["gff"], "rt")
on.exit(close(con_gff), add = TRUE)
chunks <- list(); n <- 0L
repeat {
  lines <- readLines(con_gff, n = 50000L)
  if (length(lines) == 0L) break
  lines <- lines[!startsWith(lines, "#")]
  lines <- lines[grepl("\t(exon|CDS)\t|\tID=rna-", lines)]
  if (length(lines) == 0L) next
  fields <- tstrsplit(lines, "\t", fixed = TRUE, keep = c(1L, 3L, 4L, 5L, 7L, 8L, 9L))
  transcript_row <- grepl("(?:^|;)ID=rna-", fields[[7]], perl = TRUE)
  parent <- sub(".*(?:^|;)Parent=rna-([^;,]+).*", "\\1", fields[[7]], perl = TRUE)
  # NCBI appends -2 to the placement ID when one RNA occurs on two loci.
  # Parent ownership is retained while the accession remains version-exact.
  parent <- sub("-[0-9]+$", "", parent)
  ids <- ifelse(transcript_row,
    sub(".*(?:^|;)transcript_id=([^;]+).*", "\\1", fields[[7]], perl = TRUE), parent)
  keep <- ids %chin% wanted
  if (!any(keep)) next
  n <- n + 1L
  chunks[[n]] <- data.table(accession = ids[keep], target_sequence_accession = fields[[1]][keep],
                            feature = fields[[2]][keep], start = as.integer(fields[[3]][keep]),
                            end = as.integer(fields[[4]][keep]), strand = fields[[5]][keep],
                            phase = fields[[6]][keep], attributes = fields[[7]][keep])
}
gff <- rbindlist(chunks)
rm(chunks)
gff <- merge(gff, report, by = "target_sequence_accession", all.x = TRUE, sort = FALSE)
model_contigs <- dbGetQuery(con, "SELECT seq_region_name FROM model_regions")[[1]]
# A candidate interval must lie on a report-declared region backed by the model.
gff[, target_ok := !is.na(model_contig) & model_contig %chin% model_contigs &
      start >= 1L & end <= sequence_length & start <= end]
setkey(gff, accession)

# Biostrings indexes the accession up to the first space in the FASTA header.
read_accessions <- function(path, accessions, kind = "DNA") {
  fa <- if (kind == "AA") readAAStringSet(path) else readDNAStringSet(path)
  names(fa) <- sub(" .*", "", names(fa))
  fa <- fa[names(fa) %in% accessions]
  setNames(as.character(fa), names(fa))
}
rna <- read_accessions(paths["rna"], wanted)
protein <- read_accessions(paths["protein"], unique(mane$RefSeq_prot), kind = "AA")

roots <- unique(mane$enst_root)
quoted <- paste(DBI::dbQuoteString(con, roots), collapse = ",")
model <- as.data.table(dbGetQuery(con, paste0(
  "SELECT t.transcript_index, t.transcript_stable_id, t.transcript_version, t.seq_region_name, t.strand," ,
  " t.transcript_flags, t.cds_start, t.cds_end, t.codon_table, t.exons,",
  " t.source_transcript_id = g.canonical_transcript_id AS canonical,",
  " CAST(t.pre_cds_sequence AS VARCHAR) pre, CAST(t.cds_sequence AS VARCHAR) cds,",
  " CAST(t.post_cds_sequence AS VARCHAR) post FROM model_transcripts t",
  " JOIN ensembl_core.gene g ON g.gene_id = t.source_gene_id",
  " WHERE t.transcript_stable_id IN (", quoted, ")")))
setkey(model, transcript_stable_id)

native <- as.data.table(dbGetQuery(con, paste(
  "SELECT t.transcript_index, t.transcript_stable_id, t.transcript_version,",
  "t.gene_stable_id, t.source_gene_id,",
  "t.source_transcript_id = g.canonical_transcript_id AS canonical,",
  "(t.transcript_flags & 4096) <> 0 AS gencode_basic,",
  "false AS gencode_primary, false AS native_mane,",
  "count(*) FILTER (WHERE t.source_transcript_id = g.canonical_transcript_id)",
  "OVER (PARTITION BY t.source_gene_id) = 0 AS no_retained_canonical,",
  "r.model_sha256 FROM model_transcripts t",
  "JOIN ensembl_core.gene g ON g.gene_id = t.source_gene_id",
  "CROSS JOIN model_receipt r ORDER BY t.transcript_index")))
stopifnot(dbGetQuery(con, "SELECT count(*) n FROM ensembl_core.gene")$n == 64102L,
  dbGetQuery(con, "SELECT count(*) FILTER (WHERE (transcript_flags & 11264) <> 0) n FROM model_transcripts")$n == 0L,
  nrow(native) == 195379L, uniqueN(native$source_gene_id) == 57343L,
  sum(native$canonical) == 57298L, sum(native$gencode_basic) == 99260L,
  uniqueN(native[no_retained_canonical == TRUE]$source_gene_id) == 45L,
  all(!native$gencode_primary), all(!native$native_mane))

# Fetch target reference exon sequences from the indexed primary-assembly FASTA.
# The model and NCBI annotation are checked against the SAME target assembly.
fa <- FaFile(paths["reference_fasta"])
open(fa)
on.exit(close(fa), add = TRUE)
positions <- function(x) paste(x$start, x$end, sep = "-")
geometry <- function(x, strand) {
  ordering <- order(if (strand == "+") x$start else -x$end)
  x <- x[ordering]
  paste(positions(x), collapse = ";")
}
# GFF phase is bases skipped at the 5' end of each CDS segment. Ensembl exon
# phase describes the codon position after the previous coding segment.
cds_geometry <- function(x, strand) {
  if (!nrow(x)) return("")
  ordering <- order(if (strand == "+") x$start else -x$end)
  x <- x[ordering]
  paste(paste(positions(x), x$phase, sep = ":"), collapse = ";")
}
# Compact per-accession GFF records prevent a source row from being merged with
# another isoform's exons or with a version-colliding annotation.
refseq_roots <- setNames(mane$enst_root, mane$RefSeq_nuc)
targets <- lapply(split(gff, by = "accession", keep.by = TRUE), function(x) {
  location_count <- nrow(x[!feature %chin% c("exon", "CDS")])
  root_candidates <- mane_candidates(model, refseq_roots[[x$accession[1]]])
  candidate_contig <- if (nrow(root_candidates) == 1L) root_candidates$seq_region_name else NULL
  x <- mane_supported_loci(x, candidate_contig)
  tx <- x[!feature %chin% c("exon", "CDS")]
  ex <- x[feature == "exon"]
  cds <- x[feature == "CDS"]
  # One accession must define one strand, locus and complete exon chain.
  if (nrow(tx) != 1L || nrow(ex) == 0L || any(!x$target_ok) ||
      uniqueN(x$target_sequence_accession) != 1L || uniqueN(x$strand) != 1L) {
    locus <- if (nrow(tx) > 0L) tx[1] else x[1]
    reason <- if (uniqueN(x$target_sequence_accession) > 1L || nrow(tx) > 1L)
      "ambiguous_target_locus" else "target_reference_unavailable"
    return(list(valid = FALSE, accession = locus$target_sequence_accession,
                start = locus$start, end = locus$end, strand = locus$strand,
                target_location_count = location_count, reason = reason))
  }
  strand <- tx$strand[1]
  exon_order <- order(if (strand == "+") ex$start else -ex$end)
  cds_order <- order(if (strand == "+") cds$start else -cds$end)
  ex <- ex[exon_order]
  cds <- cds[cds_order]
  ranges <- GRanges(rep(tx$model_contig[1], nrow(ex)), IRanges(ex$start, ex$end))
  bases <- as.character(scanFa(fa, param = ranges))
  if (strand == "-") bases <- as.character(reverseComplement(DNAStringSet(bases)))
  cdna <- paste0(bases, collapse = "")
  list(valid = TRUE, contig = tx$model_contig[1], accession = tx$target_sequence_accession[1],
       start = tx$start[1], end = tx$end[1], strand = strand,
       target_location_count = location_count,
       exons = geometry(ex, strand), cds = cds_geometry(cds, strand),
       cdna = cdna, cds_rows = cds, exon_rows = ex)
})
rm(gff)

classify <- function(row) {
  refseq <- row$RefSeq_nuc
  target <- targets[[refseq]]
  candidate <- mane_candidates(model, row$enst_root)
  result <- list(mane_row = row$mane_row, mane_release = "1.5", mane_status = row$MANE_status,
    refseq_nuc = refseq, refseq_prot = row$RefSeq_prot,
    ensembl_nuc = row$Ensembl_nuc, ensembl_prot = row$Ensembl_prot,
    source_assembly = "GRCh38", target_assembly = "GRCh37.p13",
    mane_source_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz",
    target_annotation_source_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
    target_sequence_accession = NA_character_, target_start = NA_integer_, target_end = NA_integer_,
    target_strand = NA_character_, target_exons = NA_character_, target_cds = NA_character_,
    target_location_count = 0L,
    target_annotation_release = "GCF_000001405.25", candidate_enst = NA_character_,
    candidate_enst_version = NA_integer_, transcript_index = NA_integer_,
    canonical = NA, gencode_basic = NA, gencode_primary = FALSE,
    exact_refseq = FALSE, exon_chain_match = NA, utr_exon_chain_match = NA,
    cds_phase_match = NA, spliced_sequence_match = NA, translated_sequence_match = NA,
    reference_difference = NA, reference_difference_bases = NA_integer_,
    cds_reference_difference_bases = NA_integer_,
    target_reference_sha256 = NA_character_, refseq_rna_sha256 = NA_character_,
    refseq_protein_sha256 = NA_character_, mapping_status = mane_target_status(target),
    mapping_label = NA_character_,
    model_sha256 = model_hash, source_manifest_sha256 = receipt$source_manifest_sha256,
    reference_sha256 = receipt$reference_sha256, mane_sha256 = sha["mane"],
    report_sha256 = sha["report"], gff_sha256 = sha["gff"],
    rna_sha256 = sha["rna"], protein_sha256 = sha["protein"])
  if (nrow(candidate) == 1L) {
    result$candidate_enst <- candidate$transcript_stable_id
    result$candidate_enst_version <- as.integer(candidate$transcript_version)
    result$canonical <- candidate$canonical
    result$gencode_basic <- bitwAnd(as.integer(candidate$transcript_flags), 4096L) != 0L
  }
  if (is.null(target)) return(result)
  result$exact_refseq <- TRUE
  result$target_location_count <- target$target_location_count
  ref_rna <- unname(rna[refseq])
  ref_protein <- unname(protein[row$RefSeq_prot])
  if (!is.na(ref_rna))
    result$refseq_rna_sha256 <- digest(ref_rna, algo = "sha256", serialize = FALSE)
  if (!is.na(ref_protein))
    result$refseq_protein_sha256 <- digest(ref_protein, algo = "sha256", serialize = FALSE)
  result$target_sequence_accession <- target$accession
  result$target_start <- target$start; result$target_end <- target$end
  result$target_strand <- target$strand
  if (!isTRUE(target$valid)) {
    result$mapping_status <- target$reason
    return(result)
  }
  result$target_exons <- target$exons; result$target_cds <- target$cds
  result$target_reference_sha256 <- digest(target$cdna, algo = "sha256", serialize = FALSE)
  result$mapping_status <- "refseq_only_no_gencode19_match"
  if (!nrow(candidate)) return(result)
  gate <- mane_pair_gate(target, candidate,
    if (is.na(ref_rna)) NULL else ref_rna,
    if (is.na(ref_protein)) NULL else ref_protein)
  for (name in names(gate$evidence)) result[[name]] <- gate$evidence[[name]]
  result$mapping_status <- gate$status
  if (!gate$status %in% c("exact_model_match", "cds_exact_utr_differs")) return(result)
  result$mapping_label <- if (gate$status == "exact_model_match")
    "MANE mapped to GRCh37" else "MANE mapped to GRCh37 (coding region only)"
  selected <- gate$selected
  result$candidate_enst <- selected$transcript_stable_id
  result$candidate_enst_version <- as.integer(selected$transcript_version)
  result$transcript_index <- as.integer(selected$transcript_index)
  result$canonical <- selected$canonical
  result$gencode_basic <- bitwAnd(as.integer(selected$transcript_flags), 4096L) != 0L
  result
}
rows <- lapply(seq_len(nrow(mane)), function(i) classify(mane[i]))
mapping <- rbindlist(rows)
stopifnot(nrow(mapping) == nrow(mane), !anyDuplicated(mapping$mane_row),
          all(mapping$gencode_primary == FALSE),
          all(is.na(mapping$transcript_index) ==
        !mapping$mapping_status %in% c("exact_model_match", "cds_exact_utr_differs")))
setorder(mapping, mane_row)
# Stable checksum of ordered UTF-8 TSV rows, including column names and explicit NA.
checksum_file <- tempfile("mane-mapping-")
fwrite(mapping, checksum_file, sep = "\t", na = "\\N", quote = TRUE)
checksum <- digest(checksum_file, algo = "sha256", file = TRUE)
unlink(checksum_file)
counts <- mapping[, .(n = .N), by = mapping_status][order(mapping_status)]
receipt_out <- data.table(run_date = "2026-09-28", mane_release = "1.5",
  mane_source_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz",
  ncbi_source_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/",
  ensembl_core_source = "ensembldb.ensembl.org:3337/homo_sapiens_core_116_37",
  reference_source_url = "https://ftp.ensembl.org/pub/grch37/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
  model_sha256 = model_hash, source_manifest_sha256 = receipt$source_manifest_sha256,
  reference_sha256 = sha["reference"], reference_fasta_sha256 = sha["reference_fasta"],
  mane_sha256 = sha["mane"], report_sha256 = sha["report"],
  gff_sha256 = sha["gff"], rna_sha256 = sha["rna"], protein_sha256 = sha["protein"],
  row_count = nrow(mapping), native_gene_count = uniqueN(native$source_gene_id),
  retained_canonical_gene_count = sum(native$canonical),
  no_retained_canonical_gene_count = uniqueN(native[no_retained_canonical == TRUE]$source_gene_id),
  gencode_basic_transcript_count = sum(native$gencode_basic),
  exact_refseq_count = sum(mapping$exact_refseq),
  status_counts = paste(paste(counts$mapping_status, counts$n, sep = ":"), collapse = ";"),
  relation_sha256 = checksum)
# Export with a separate DuckDB process/connection from all native readers.
parquet <- file.path(out, "mane_grch37_mapping.parquet")
export <- dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = ":memory:")
dbWriteTable(export, "mapping", as.data.frame(mapping))
dbExecute(export, paste0("COPY mapping TO ", dbQuoteString(export, parquet), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
dbWriteTable(export, "native", as.data.frame(native))
dbExecute(export, paste0("COPY native TO ", dbQuoteString(export, file.path(out, "grch37_transcript_authorities.parquet")),
                         " (FORMAT PARQUET, COMPRESSION ZSTD)"))
dbDisconnect(export, shutdown = TRUE)
fwrite(receipt_out, file.path(out, "mane_grch37_receipt.csv"), quote = TRUE)
print(counts)
print(receipt_out)
