#!/usr/bin/env Rscript
# Structural HGVS domain against executable VEP 116 (digest-pinned Docker).
#
# VEP 116 emits no HGVS for symbolic structural alleles or breakends (checked
# below). The supported domain is therefore validated on the equivalent literal
# edit that the native builder derives from an exact-span symbolic allele:
#   * genomic strings: builder hgvs_g vs VEP --hgvsg --shift_hgvs 0 on the
#     literal edit (the builder declares normalization = 'none');
#   * transcript strings: DuckVEP hgvs on the literal edit vs VEP --hgvs (its
#     default 3' shifting) on the same edit;
#   * failure controls: --shift_hgvs 1 genomic counterexamples are retained.
# Loads an immutable copy of the built extension; writes no receipt.
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(jsonlite) })
source("r/Rduckvep/R/builders.R")
source("r/Rduckvep/R/structural_hgvs.R")
fasta <- "/root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
cache <- "/root/.cache/duckhts/vep/cache-grch38-chr21"
model <- "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb"
counterexamples <- "test/duckvep/conformance/data/structural_hgvs_vep116_counterexamples.tsv"
refused_file <- "test/duckvep/conformance/data/structural_hgvs_vep116_refused.tsv"
n_events <- 600L # random events; targeted duplications are appended
directory <- tempfile("structural-hgvs-")
dir.create(directory)
binary <- file.path(directory, "duckvep.duckdb_extension")
stopifnot(file.copy("build/release/duckvep.duckdb_extension", binary))
con <- dbConnect(duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
on.exit({dbDisconnect(con, shutdown = TRUE); unlink(directory, recursive = TRUE)}, add = TRUE)
q <- function(value) as.character(dbQuoteString(con, value))
dbExecute(con, paste("LOAD", q(binary)))
dbExecute(con, paste("ATTACH", q(model), "AS m (READ_ONLY)"))
region <- dbGetQuery(con, "SELECT seq_region FROM m.bench_regions WHERE chrom='21'")$seq_region

# Deterministic corpus: symbolic DEL/DUP/DUP:TANDEM/INV starts inside chr21 exons.
exons <- dbGetQuery(con, paste(
  "SELECT e.exon_start, e.exon_end FROM m.bench_exons e JOIN m.bench_transcripts t",
  "USING (transcript_index) WHERE t.seq_region =", region,
  "AND e.exon_end - e.exon_start >= 40 ORDER BY t.transcript_index, e.exon_start"))
set.seed(4116L)
fetch <- function(from, to) {
  out <- system2("samtools", c("faidx", "-n", "100000", shQuote(fasta),
    paste0("21:", from, "-", to)), stdout = TRUE)
  stopifnot(length(out) == 2L * length(from))
  out[seq(2L, length(out), by = 2L)]
}
draw <- function(n) {
  pick <- sample(nrow(exons), n, replace = TRUE)
  data.frame(pos = vapply(seq_len(n), function(i)
    as.numeric(exons$exon_start[pick[i]] + sample(0:20, 1L)), 0),
    length = sample(c(1:12, 13:40), n, replace = TRUE))
}
random <- draw(n_events)
random$alt <- rep(c("<DEL>", "<DUP>", "<DUP:TANDEM>", "<INV>"), length.out = n_events)
# Targeted duplications whose first duplicated base equals the next reference
# base but whose unit is not repeated, plus fully repeated units (excluded).
pool <- draw(3000L)
pool$alt <- "<DUP>"
pool_seq <- fetch(pool$pos, pool$pos + 2 * pool$length)
span <- substr(pool_seq, 2L, pool$length + 1L)
flank <- substr(pool_seq, pool$length + 2L, 2L * pool$length + 1L)
partial <- which(substr(span, 1L, 1L) == substr(flank, 1L, 1L) & span != flank & pool$length >= 2L)
full <- which(span == flank & pool$length >= 2L)
stopifnot(length(partial) >= 60L)
targeted <- pool[c(head(partial, 60L), head(full, 10L)), ]
corpus <- rbind(random, targeted)
pos <- corpus$pos
lengths <- corpus$length
kinds <- corpus$alt
end <- pos + lengths
n_events <- nrow(corpus)
refseq <- fetch(pos, end + lengths)
stopifnot(all(nchar(refseq) == 2L * lengths + 1L))
# DEL and INV need no flank; every second one is supplied without it.
short <- kinds %in% c("<DEL>", "<INV>") & seq_len(n_events) %% 2L == 0L
refseq[short] <- substr(refseq[short], 1L, lengths[short] + 1L)
events <- data.frame(event_index = seq_len(n_events) - 1L, chrom = "21", pos = pos,
  ref = substr(refseq, 1L, 1L), alt = kinds, info = paste0("END=", end),
  stringsAsFactors = FALSE)
refs <- data.frame(event_index = events$event_index, reference_sequence = refseq)
dbWriteTable(con, "sh_events", events, temporary = TRUE)
dbWriteTable(con, "sh_refs", refs, temporary = TRUE)
built <- rduckvep_prepare_structural_hgvs(con, "sh_events", "sh_refs")
stopifnot(nrow(built) == n_events, all(built$event_index == events$event_index))
ok <- built$hgvs_status == "supported"
inversion_reasons <- table(built$hgvs_reason[!ok])
cat("supported:", sum(ok), "of", n_events, "; unsupported reasons:",
    paste(names(inversion_reasons), inversion_reasons, collapse = ", "), "\n")
stopifnot(all(built$hgvs_reason[!ok] %in% c("inversion_not_reducible", "inversion_identity",
                                            "duplication_adjacent_repeat")),
          sum(ok) > 0.6 * n_events, sum(built$hgvs_reason == "duplication_adjacent_repeat", na.rm = TRUE) > 0L)

runner <- function(vcf, out, ...) {
  rc <- system2("scripts/run_species_vep116_docker.sh", c("homo_sapiens", "GRCh38", "116",
    cache, fasta, vcf, out, ...))
  stopifnot(identical(rc, 0L))
  lapply(readLines(out), fromJSON, simplifyVector = FALSE)
}
write_vcf <- function(path, chrom, position, ref, alt, id, info = ".") {
  rows <- data.frame(position = position, ref = ref, alt = alt, id = id, info = info,
                     stringsAsFactors = FALSE)
  rows <- rows[order(rows$position, rows$id), ]
  position <- rows$position; ref <- rows$ref; alt <- rows$alt; id <- rows$id; info <- rows$info
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=21,length=46709983>",
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(chrom, format(position, scientific = FALSE, trim = TRUE), id, ref, alt, ".", "PASS", info,
          sep = "\t")), path)
}
supported <- built[ok, ]
literal_vcf <- file.path(directory, "literal.vcf")
write_vcf(literal_vcf, "21", supported$literal_position, supported$literal_reference,
          supported$literal_alternate, paste0("e", supported$event_index))
by_id <- function(records) {
  names(records) <- vapply(records, function(x) x$id, "")
  records
}
transcript_field <- function(record, field) {
  out <- vapply(record$transcript_consequences, function(t)
    if (is.null(t[[field]])) NA_character_ else t[[field]], "")
  names(out) <- vapply(record$transcript_consequences, function(t) t$transcript_id, "")
  out
}
vep0 <- by_id(runner(literal_vcf, file.path(directory, "shift0.json"), "--hgvs", "--hgvsg",
                     "--shift_hgvs", "0"))
vep1 <- by_id(runner(literal_vcf, file.path(directory, "shift1.json"), "--hgvs", "--hgvsg"))
stopifnot(length(vep0) == nrow(supported), length(vep1) == nrow(supported))

# 1. Genomic strings, unshifted: every event that VEP reports on a transcript.
compared <- 0L
mismatches <- list()
for (i in seq_len(nrow(supported))) {
  g <- unique(na.omit(transcript_field(vep0[[paste0("e", supported$event_index[i])]], "hgvsg")))
  if (!length(g)) next
  if (!(length(g) == 1L && identical(g, supported$hgvs_g[i]))) {
    mismatches[[length(mismatches) + 1L]] <- data.frame(event_index = supported$event_index[i],
      allele = events$alt[supported$event_index[i] + 1L], builder = supported$hgvs_g[i],
      vep = paste(g, collapse = ","), literal = paste(supported$literal_reference[i],
        supported$literal_alternate[i], sep = ">"))
    next
  }
  compared <- compared + 1L
}
if (length(mismatches)) { print(head(do.call(rbind, mismatches), 20L)) }
stopifnot(!length(mismatches))
stopifnot(compared > 0.9 * nrow(supported))

# 2. Failure control: VEP's default 3' shift moves some descriptions, so the
# builder's declared normalization = 'none' is not a claim about shifted output.
shifted <- do.call(rbind, lapply(seq_len(nrow(supported)), function(i) {
  g <- unique(na.omit(transcript_field(vep1[[paste0("e", supported$event_index[i])]], "hgvsg")))
  if (length(g) == 1L && !identical(g, supported$hgvs_g[i]))
    data.frame(event_index = supported$event_index[i], allele = events$alt[supported$event_index[i] + 1L],
               builder_hgvs_g = supported$hgvs_g[i], vep116_shifted_hgvsg = g,
               stringsAsFactors = FALSE)
}))
stopifnot(!is.null(shifted), nrow(shifted) > 0L)
shifted <- shifted[order(shifted$event_index), ]
write.table(head(shifted, 25L), counterexamples, sep = "\t", quote = FALSE, row.names = FALSE)

# 2b. Exclusion controls: what VEP does on the events the builder refuses.
refused <- built[!ok, ]
complement <- function(x) chartr("ACGT", "TGCA", x)
reverse <- function(x) vapply(strsplit(x, ""), function(y) paste(rev(y), collapse = ""), "")
idx <- refused$event_index + 1L
refused_span <- substr(refseq[idx], 2L, lengths[idx] + 1L)
is_inversion <- kinds[idx] == "<INV>"
refused_vcf <- file.path(directory, "refused.vcf")
write_vcf(refused_vcf, "21", ifelse(is_inversion, pos[idx] + 1, end[idx]),
  ifelse(is_inversion, refused_span, substr(refseq[idx], lengths[idx] + 1L, lengths[idx] + 1L)),
  ifelse(is_inversion, reverse(complement(refused_span)),
    paste0(substr(refseq[idx], lengths[idx] + 1L, lengths[idx] + 1L), refused_span)),
  paste0("r", refused$event_index))
refused_vep <- by_id(runner(refused_vcf, file.path(directory, "refused.json"), "--hgvsg",
                            "--shift_hgvs", "0"))
refused_g <- vapply(seq_len(nrow(refused)), function(i) {
  g <- unique(na.omit(transcript_field(refused_vep[[paste0("r", refused$event_index[i])]], "hgvsg")))
  if (length(g) == 1L) g else NA_character_
}, "")
naive <- ifelse(is_inversion, NA_character_, paste0("21:g.", pos[idx] + 1,
  ifelse(lengths[idx] > 1L, paste0("_", end[idx]), ""), "dup"))
adjacent <- refused$hgvs_reason == "duplication_adjacent_repeat"
# VEP places an adjacent repeated duplication at the later copy, never at the
# nominal span the builder would have written.
stopifnot(sum(adjacent) > 0L, all(is.na(refused_g[adjacent]) | refused_g[adjacent] != naive[adjacent]),
          sum(!is.na(refused_g[adjacent]) & refused_g[adjacent] != naive[adjacent]) > 0L)
# VEP does not describe a palindromic or single-base-changing inversion as inv.
inverted_refused <- refused$hgvs_reason %in% c("inversion_not_reducible", "inversion_identity")
stopifnot(sum(inverted_refused) > 0L,
          !any(grepl("inv$", refused_g[inverted_refused]), na.rm = TRUE))

write.table(data.frame(event_index = refused$event_index, allele = kinds[idx],
  chrom_pos = paste0("21:", pos[idx]), end = end[idx], builder_reason = refused$hgvs_reason,
  literal_reference = ifelse(is_inversion, refused_span, substr(refseq[idx], lengths[idx] + 1L, lengths[idx] + 1L)),
  literal_alternate = ifelse(is_inversion, reverse(complement(refused_span)),
    paste0(substr(refseq[idx], lengths[idx] + 1L, lengths[idx] + 1L), refused_span)),
  vep116_shift0_hgvsg = refused_g, stringsAsFactors = FALSE),
  refused_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# 3. Transcript strings: DuckVEP hgvs on the literal edit vs VEP default hgvs.
queries <- c(
  "SELECT seq_region, sequence_length, chrom AS seq_region_name FROM m.bench_regions ORDER BY seq_region",
  paste("SELECT transcript_index, seq_region, transcript_start, transcript_end,",
    "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence,",
    "codon_table, pre_cds_sequence, post_cds_sequence FROM m.bench_transcripts",
    "ORDER BY transcript_index"),
  paste("SELECT transcript_index, exon_start, exon_end, exon_cdna_start,",
    "exon_cdna_end, phase, end_phase FROM m.bench_exons",
    "ORDER BY transcript_index, exon_cdna_start"))
load <- paste0("SELECT loaded FROM duckvep_model_load('structural_hgvs',",
  paste(vapply(queries, q, ""), collapse = ","), ", reference_fasta := ", q(fasta), ")")
stopifnot(isTRUE(dbGetQuery(con, load)$loaded))
dbWriteTable(con, "sh_literal", data.frame(
  event_index = supported$event_index, seq_region = region,
  position = supported$literal_position, reference = supported$literal_reference,
  alternate = supported$literal_alternate), temporary = TRUE)
invisible(dbExecute(con, paste(
  "CREATE TEMP TABLE sh_annotate_events AS SELECT event_index::UBIGINT event_index,",
  "seq_region::UINTEGER seq_region, \"position\"::UBIGINT AS \"position\", reference, alternate,",
  "NULL::UBIGINT end_position, NULL::VARCHAR structural_type, NULL::VARCHAR copy_change,",
  "NULL::UINTEGER mate_seq_region, NULL::UBIGINT mate_position FROM sh_literal")))
duck <- dbGetQuery(con, paste(
  "SELECT a.event_index, t.transcript_stable_id AS transcript_id,",
  "t.transcript_stable_id || '.' || t.transcript_version || ':' || a.transcript_hgvs AS hgvsc,",
  "a.transcript_hgvs_status AS status",
  "FROM query(duckvep_annotate_sql('sh_annotate_events', 'structural_hgvs',",
  "struct_pack(hgvs := true, upstream_distance := 5000, downstream_distance := 5000))) a",
  "JOIN m.model_transcripts t USING (transcript_index) WHERE a.transcript_hgvs IS NOT NULL",
  "ORDER BY event_index, transcript_id"))
duck_pairs <- split(paste(duck$transcript_id, duck$hgvsc), duck$event_index)
transcript_pairs <- 0L
for (i in seq_len(nrow(supported))) {
  id <- paste0("e", supported$event_index[i])
  vh <- transcript_field(vep1[[id]], "hgvsc")
  vh <- vh[!is.na(vh)]
  expected <- sort(paste(names(vh), vh))
  actual <- sort(duck_pairs[[as.character(supported$event_index[i])]])
  stopifnot(identical(as.character(expected), as.character(actual)))
  transcript_pairs <- transcript_pairs + length(expected)
}
stopifnot(transcript_pairs > 0L)

# 4. Symbolic alleles and a BND pair: VEP 116 emits no HGVS of any kind.
symbolic_vcf <- file.path(directory, "symbolic.vcf")
sym <- events[events$event_index < 40L, ]
write_vcf(symbolic_vcf, "21", sym$pos, "N", sym$alt, paste0("s", sym$event_index),
          sym$info)
symbolic <- runner(symbolic_vcf, file.path(directory, "symbolic.json"), "--hgvs", "--hgvsg")
bnd_vcf <- file.path(directory, "bnd.vcf")
write_vcf(bnd_vcf, "21", c(7467462, 7467900), "N", c("N[21:7467900[", "]21:7467462]N"),
          c("bnd1", "bnd2"), "SVTYPE=BND")
bnd <- runner(bnd_vcf, file.path(directory, "bnd.json"), "--hgvs", "--hgvsg")
emitted <- vapply(c(symbolic, bnd), function(r) sum(vapply(
  c(r$transcript_consequences, r$intergenic_consequences), function(t)
    any(grepl("^hgvs", names(t))), TRUE)), 0L)
stopifnot(length(symbolic) == 40L, length(bnd) == 2L, all(emitted == 0L),
          sum(vapply(c(symbolic, bnd), function(r) length(r$transcript_consequences), 0L)) > 0L)

# 5. Builder failure controls, none of which VEP could validate.
controls <- data.frame(event_index = 0:11, chrom = "21", pos = 7467462,
  ref = "T", alt = c("<DEL>", "<DEL>", "<DEL>", "<INS>", "<CNV>", "<STR7>", "<DUP>",
                     "N[21:7467900[", "<DEL>", "<INV>", "<DEL>", "TCAG"),
  info = c("END=7467465;CIPOS=-1,1", "SVTYPE=DEL", "END=7467465", "END=7467462",
           "END=7467465", "END=7467465", "END=7467465;END=7467466", "SVTYPE=BND",
           "END=7467465", "END=7467465", "END=7467465", "."),
  stringsAsFactors = FALSE)
control_refs <- data.frame(event_index = c(2L, 9L, 10L, 6L),
  reference_sequence = c("TCAN", "TCAG", "ACAG", "TCAG"))
dbWriteTable(con, "sh_controls", controls, temporary = TRUE)
dbWriteTable(con, "sh_control_refs", control_refs, temporary = TRUE)
c1 <- rduckvep_prepare_structural_hgvs(con, "sh_controls", "sh_control_refs")
expected_controls <- c("unsupported:imprecise", "unavailable:missing_end",
  "unavailable:reference_alphabet", "unsupported:symbolic_insertion", "unsupported:copy_number",
  "unsupported:repeat_expansion", "unavailable:duplicate_end", "unsupported:breakend",
  "unavailable:missing_reference_sequence", "unsupported:inversion_not_reducible",
  "unavailable:reference_anchor_mismatch", "unsupported:literal_allele")
stopifnot(identical(paste0(c1$hgvs_status, ":", c1$hgvs_reason), expected_controls),
          all(is.na(c1$hgvs_g)))
cat(sprintf(paste0("VEP 116: %d exact-span DEL/DUP/INV events, %d genomic strings equal (shift 0),",
  " %d transcript strings equal, %d shift-1 counterexamples retained, %d refused events checked,",
  " symbolic and BND emit no HGVS\n"),
  nrow(supported), compared, transcript_pairs, nrow(shifted), nrow(refused)))
