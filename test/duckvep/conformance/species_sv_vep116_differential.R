#!/usr/bin/env Rscript
# Species structural evidence against executable VEP 116 (digest-pinned Docker).
#
#   species_sv_vep116_differential.R SPECIES ASSEMBLY CACHE_VERSION MODEL_DB FASTA \
#       CACHE_DIR OUTPUT_DIR [SEED] [PER_STRATUM]
#
# The model is one of the pinned species models built by
# scripts/build_species_model.R (tables model_regions and model_transcripts).
# Everything is derived from that model with a seeded hash, so the corpus is
# reproducible:
#   * exact-span symbolic DEL, DUP, DUP:TANDEM and INV over 16 geometry states
#     (exon-edge crossings, exact/interior exons, introns, multi-exon spans,
#     transcript-edge crossings, CDS-edge crossings, flanks), sampled per
#     type x state x transcript strand;
#   * symbolic <INS> at exon/intron/CDS edges;
#   * BND pairs: both physical records of a reciprocal pair, four bracket
#     orientations, same- and cross-chromosome mates, endpoints at transcript
#     landmarks; VEP runs with --buffer_size 1 (its interval tree is
#     chromosome-blind across records);
#   * a small-span DEL/DUP/DUP:TANDEM/INV corpus for the structural HGVS domain.
# Per (event, transcript) consequence-term sets from VEP and DuckVEP are
# compared; nonexact pairs are written to disagreements.csv and the script
# exits nonzero. Structural HGVS is checked on the literal-equivalent edit as in
# structural_hgvs_vep116_differential.R. The extension is loaded from an
# immutable copy of DUCKVEP_EXTENSION_FILE (default build/release/...).
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(jsonlite) })
source("r/Rduckvep/R/builders.R")
source("r/Rduckvep/R/structural_geometry.R")
source("r/Rduckvep/R/structural_identity.R")
source("r/Rduckvep/R/structural_hgvs.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7L || length(args) > 9L) {
  stop("usage: species_sv_vep116_differential.R SPECIES ASSEMBLY CACHE_VERSION MODEL_DB FASTA CACHE_DIR OUTPUT_DIR [SEED] [PER_STRATUM]", call. = FALSE)
}
species <- args[[1L]]; assembly <- args[[2L]]; cache_version <- args[[3L]]
model <- normalizePath(args[[4L]], mustWork = TRUE)
fasta <- normalizePath(args[[5L]], mustWork = TRUE)
cache <- normalizePath(args[[6L]], mustWork = TRUE)
out <- args[[7L]]
seed <- if (length(args) >= 8L) as.integer(args[[8L]]) else 4116L
per <- if (length(args) >= 9L) as.integer(args[[9L]]) else 10L
stopifnot(!is.na(seed), !is.na(per), per > 0L, file.exists(paste0(fasta, ".fai")))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
extension <- Sys.getenv("DUCKVEP_EXTENSION_FILE", "build/release/duckvep.duckdb_extension")
work <- tempfile("species-sv-"); dir.create(work)
binary <- file.path(work, "duckvep.duckdb_extension")
stopifnot(file.copy(extension, binary))
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
on.exit({dbDisconnect(con, shutdown = TRUE); unlink(work, recursive = TRUE)}, add = TRUE)
q <- function(x) as.character(dbQuoteString(con, x))
run <- function(sql) invisible(dbExecute(con, sql))
run(paste("LOAD", q(binary)))
run(paste("ATTACH", q(model), "AS model (READ_ONLY)"))
run("SET threads=4")
sha256 <- function(path) tolower(strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " ")[[1L]][1L])

## ---- model -----------------------------------------------------------------------
run("CREATE TABLE mregions AS SELECT seq_region, sequence_length, seq_region_name FROM model.model_regions")
run("CREATE TEMP TABLE mnames AS SELECT transcript_index, transcript_stable_id AS tx,
     transcript_version AS version, strand, transcript_biotype AS biotype FROM model.model_transcripts")
fai <- read.delim(paste0(fasta, ".fai"), header = FALSE, stringsAsFactors = FALSE)
dbWriteTable(con, "fa_order", data.frame(chrom = fai$V1, ord = seq_len(nrow(fai))), temporary = TRUE)
stopifnot(dbGetQuery(con, "SELECT count(*) n FROM mregions r LEFT JOIN fa_order f
  ON f.chrom = r.seq_region_name WHERE f.chrom IS NULL")$n == 0L)
# All transcripts must be on the pinned model; transcript_coverage_complete follows the species differential.
load_model <- function() {
  run("CREATE TABLE l_transcripts AS SELECT transcript_index, seq_region, transcript_start, transcript_end,
       strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence,
       post_cds_sequence FROM model.model_transcripts")
  run("CREATE TABLE l_exons AS SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start,
       e.exon_cdna_end, e.phase, e.end_phase FROM model.model_transcripts, unnest(exons) AS x(e)")
  run("CREATE TABLE l_mirna AS SELECT transcript_index, r.mature_mirna_start, r.mature_mirna_end
       FROM model.model_transcripts, unnest(mature_mirna_regions) AS x(r)")
  sql <- paste0("SELECT loaded FROM duckvep_model_load(", q("sv_species"), ",",
    q("SELECT seq_region, sequence_length, seq_region_name FROM mregions ORDER BY seq_region"), ",",
    q(paste("SELECT * FROM l_transcripts ORDER BY seq_region, transcript_start, transcript_index")), ",",
    q("SELECT * FROM l_exons ORDER BY transcript_index, exon_cdna_start"), ",",
    "mature_mirna_query := ", q("SELECT * FROM l_mirna ORDER BY transcript_index, mature_mirna_start"),
    ", transcript_coverage_complete := true, reference_fasta := ", q(fasta), ")")
  stopifnot(isTRUE(dbGetQuery(con, sql)$loaded))
}
load_model()
model_sha <- dbGetQuery(con, "SELECT model_sha256 FROM model.model_receipt")$model_sha256

## ---- seeded corpus ---------------------------------------------------------------
run(sprintf("CREATE TEMP TABLE g_tx AS SELECT transcript_index, seq_region, seq_region_name AS chrom,
  transcript_start::BIGINT AS ts, transcript_end::BIGINT AS te, strand::INTEGER AS strand,
  cds_start::BIGINT AS cs, cds_end::BIGINT AS ce, sequence_length::BIGINT AS slen
  FROM model.model_transcripts WHERE transcript_start > 6000
  AND transcript_end + 6000 < sequence_length AND transcript_end - transcript_start < 400000"))
run("CREATE TEMP TABLE g_ex AS SELECT * FROM (
  SELECT t.transcript_index, e.exon_start::BIGINT AS es, e.exon_end::BIGINT AS ee,
    lead(e.exon_start) OVER w AS next_es
  FROM model.model_transcripts t, unnest(exons) AS x(e)
  WINDOW w AS (PARTITION BY t.transcript_index ORDER BY e.exon_start)) x
  WHERE transcript_index IN (SELECT transcript_index FROM g_tx)")
h <- function(k, n) sprintf("CAST(hash(x.transcript_index, x.es, %d, %d) %% %d AS BIGINT)", seed, k, n)
# Span states: (state, start, end, filter). Offsets are seeded hashes of the exon.
exon_states <- list(
  c("exon_start_cross", sprintf("x.es - 1 - %s", h(1, 40)), sprintf("x.es + %s", h(2, 40)), "true"),
  c("exon_end_cross", sprintf("x.ee - %s", h(1, 40)), sprintf("x.ee + 1 + %s", h(2, 40)), "true"),
  c("exon_exact", "x.es", "x.ee", "true"),
  c("exon_interior", sprintf("x.es + 2 + %s", h(1, 4)), sprintf("x.ee - 2 - %s", h(2, 4)),
    "x.ee - x.es >= 14"),
  c("intron_interior", sprintf("x.ee + 2 + %s", h(1, 4)), sprintf("x.next_es - 2 - %s", h(2, 4)),
    "x.next_es - x.ee >= 20"),
  c("multi_exon", sprintf("x.ee - LEAST(%s, x.ee - x.es)", h(1, 20)),
    sprintf("x.next_es + LEAST(%s, 5)", h(2, 20)), "x.next_es IS NOT NULL"))
tx_states <- list(
  c("transcript_left_cross", sprintf("t.ts - 1 - CAST(hash(t.transcript_index, %d, 3) %% 2000 AS BIGINT)", seed),
    sprintf("t.ts + CAST(hash(t.transcript_index, %d, 4) %% 500 AS BIGINT)", seed), "true"),
  c("transcript_right_cross", sprintf("t.te - CAST(hash(t.transcript_index, %d, 4) %% 500 AS BIGINT)", seed),
    sprintf("t.te + 1 + CAST(hash(t.transcript_index, %d, 3) %% 2000 AS BIGINT)", seed), "true"),
  c("transcript_exact", "t.ts", "t.te", "true"),
  c("flank_left", sprintf("t.ts - 3000 - CAST(hash(t.transcript_index, %d, 5) %% 1000 AS BIGINT)", seed),
    sprintf("t.ts - 3000 - CAST(hash(t.transcript_index, %d, 5) %% 1000 AS BIGINT) + CAST(hash(t.transcript_index, %d, 6) %% 400 AS BIGINT)", seed, seed), "true"),
  c("flank_right", sprintf("t.te + 1000 + CAST(hash(t.transcript_index, %d, 5) %% 1000 AS BIGINT)", seed),
    sprintf("t.te + 1000 + CAST(hash(t.transcript_index, %d, 5) %% 1000 AS BIGINT) + CAST(hash(t.transcript_index, %d, 6) %% 400 AS BIGINT)", seed, seed), "true"),
  c("cds_start_cross", sprintf("t.cs - CAST(hash(t.transcript_index, %d, 7) %% 30 AS BIGINT)", seed),
    sprintf("t.cs + CAST(hash(t.transcript_index, %d, 8) %% 30 AS BIGINT)", seed), "t.cs IS NOT NULL"),
  c("cds_end_cross", sprintf("t.ce - CAST(hash(t.transcript_index, %d, 7) %% 30 AS BIGINT)", seed),
    sprintf("t.ce + CAST(hash(t.transcript_index, %d, 8) %% 30 AS BIGINT)", seed), "t.ce IS NOT NULL"))
span_sql <- c(
  vapply(exon_states, function(s) sprintf("SELECT x.transcript_index, '%s' AS state, %s AS s, %s AS e
    FROM g_ex x WHERE %s", s[1L], s[2L], s[3L], s[4L]), ""),
  vapply(tx_states, function(s) sprintf("SELECT t.transcript_index, '%s' AS state, %s AS s, %s AS e
    FROM g_tx t WHERE %s", s[1L], s[2L], s[3L], s[4L]), ""))
ins_states <- list(
  c("ins_before_exon", "x.es - 1", "true"), c("ins_exon_end", "x.ee", "true"),
  c("ins_exon_interior", sprintf("x.es + %s %% GREATEST(x.ee - x.es + 1, 1)", h(1, 100000)), "true"),
  c("ins_intron", sprintf("x.ee + 1 + %s %% GREATEST(x.next_es - x.ee - 1, 1)", h(1, 100000)),
    "x.next_es - x.ee > 2"))
ins_sql <- c(vapply(ins_states, function(s) sprintf("SELECT x.transcript_index, '%s' AS state, %s AS s, %s AS e
    FROM g_ex x WHERE %s", s[1L], s[2L], s[2L], s[3L]), ""),
  "SELECT t.transcript_index, 'ins_cds_start' AS state, t.cs AS s, t.cs AS e FROM g_tx t WHERE t.cs IS NOT NULL",
  "SELECT t.transcript_index, 'ins_cds_end' AS state, t.ce AS s, t.ce AS e FROM g_tx t WHERE t.ce IS NOT NULL")
run(paste0("CREATE TEMP TABLE g_geom AS ", paste(span_sql, collapse = " UNION ALL ")))
run(paste0("CREATE TEMP TABLE g_ins AS ", paste(ins_sql, collapse = " UNION ALL ")))
types <- "(VALUES ('DEL','<DEL>','LOSS'), ('DUP','<DUP>','GAIN'), ('TDUP','<DUP:TANDEM>','GAIN'), ('INV','<INV>','NEUTRAL')) o(stype, alt, copy_change)"
run(sprintf("CREATE TEMP TABLE g_pick AS SELECT * FROM (
  SELECT g.transcript_index, o.stype, o.alt, o.copy_change, g.state, t.strand, t.chrom, t.seq_region,
    g.s, g.e, row_number() OVER (PARTITION BY o.stype, g.state, t.strand
      ORDER BY hash(g.transcript_index, g.s, g.e, o.stype, %d)) AS rk
  FROM g_geom g JOIN g_tx t USING (transcript_index) CROSS JOIN %s
  WHERE g.s >= 2 AND g.e >= g.s AND g.e <= t.slen AND g.e - g.s + 1 <= 500000) WHERE rk <= %d", seed, types, per))
run(sprintf("CREATE TEMP TABLE g_pick_ins AS SELECT * FROM (
  SELECT g.transcript_index, 'INS' AS stype, '<INS>' AS alt, 'UNKNOWN' AS copy_change, g.state, t.strand, t.chrom,
    t.seq_region, g.s, g.e, row_number() OVER (PARTITION BY g.state, t.strand
      ORDER BY hash(g.transcript_index, g.s, %d)) AS rk
  FROM g_ins g JOIN g_tx t USING (transcript_index)
  WHERE g.s >= 2 AND g.s < t.slen) WHERE rk <= %d", seed, per))
# One row per distinct physical event; earliest stratum wins, so strata stay disjoint.
run("CREATE TEMP TABLE g_events0 AS SELECT * FROM (
  SELECT *, row_number() OVER (PARTITION BY chrom, s, e, stype ORDER BY state, transcript_index) AS dup
  FROM (SELECT * FROM g_pick UNION ALL BY NAME SELECT * FROM g_pick_ins)) WHERE dup = 1")
run("CREATE TEMP TABLE g_events AS SELECT (row_number() OVER (ORDER BY f.ord, e.s, e.e, e.stype) - 1)::UBIGINT
  AS event_index, e.*, f.ord AS chrom_order FROM g_events0 e JOIN fa_order f ON f.chrom = e.chrom")

## ---- BND pairs ---------------------------------------------------------------------
pt_states <- list(c("transcript_start", "t.ts"), c("transcript_mid", "(t.ts + t.te) // 2"),
  c("transcript_end", "t.te"), c("upstream_flank", "t.ts - 1"), c("downstream_flank", "t.te + 1"),
  c("cds_start", "t.cs"), c("cds_end", "t.ce"))
ex_pt <- list(c("exon_start", "x.es"), c("exon_end", "x.ee"), c("before_exon_start", "x.es - 1"),
  c("after_exon_end", "x.ee + 1"), c("exon_mid", "(x.es + x.ee) // 2"))
pt_sql <- c(vapply(pt_states, function(s) sprintf("SELECT t.transcript_index, '%s' AS state, %s AS pos
  FROM g_tx t WHERE %s IS NOT NULL", s[1L], s[2L], s[2L]), ""),
  vapply(ex_pt, function(s) sprintf("SELECT x.transcript_index, '%s' AS state, %s AS pos FROM g_ex x", s[1L], s[2L]), ""))
run(paste0("CREATE TEMP TABLE g_pt0 AS ", paste(pt_sql, collapse = " UNION ALL ")))
run("CREATE TEMP TABLE g_pt AS SELECT (row_number() OVER (ORDER BY p.transcript_index, p.state, p.pos) - 1) AS rn,
  p.transcript_index, p.state, p.pos, t.chrom, t.seq_region, t.strand, t.slen FROM g_pt0 p
  JOIN g_tx t USING (transcript_index) WHERE p.pos >= 2 AND p.pos < t.slen")
bper <- max(2L, per %/% 3L)
run(sprintf("CREATE TEMP TABLE b_pick AS SELECT * FROM (
  SELECT l.rn AS lrn, o.orientation, l.state AS local_state, l.strand, row_number() OVER (
    PARTITION BY l.state, o.orientation, l.strand ORDER BY hash(l.rn, o.orientation, %d)) AS rk
  FROM g_pt l CROSS JOIN (VALUES ('N[M['), ('N]M]'), (']M]N'), ('[M[N')) o(orientation)) WHERE rk <= %d", seed, bper))
run(sprintf("CREATE TEMP TABLE b_pairs AS SELECT (row_number() OVER (ORDER BY hash(b.lrn, b.orientation, %d)) - 1) AS pair_index,
  b.orientation, b.local_state, l.chrom AS lchrom, l.seq_region AS lreg, l.pos AS lpos, l.strand AS lstrand,
  m.state AS mate_state, m.chrom AS mchrom, m.seq_region AS mreg, m.pos AS mpos, m.strand AS mstrand
  FROM b_pick b JOIN g_pt l ON l.rn = b.lrn
  JOIN g_pt m ON m.rn = CAST(hash(b.lrn, %d, 9) %% (SELECT count(*) FROM g_pt) AS BIGINT)
  WHERE NOT (l.chrom = m.chrom AND l.pos = m.pos)", seed, seed))
# The mate record: N[M[ <-> ]L]N, ]M]N <-> N[L[, N]M] <-> N]L], [M[N <-> [L[N.
run("CREATE TEMP TABLE b_records AS
  SELECT pair_index, 'a' AS side, format('bnd{}a', pair_index) AS id, format('bnd{}b', pair_index) AS mate_id,
    lchrom AS chrom, lreg AS seq_region, lpos AS pos, mchrom AS mate_chrom, mreg AS mate_seq_region, mpos AS mate_pos,
    orientation AS shape, local_state AS state, lstrand AS strand FROM b_pairs
  UNION ALL
  SELECT pair_index, 'b', format('bnd{}b', pair_index), format('bnd{}a', pair_index), mchrom, mreg, mpos, lchrom, lreg, lpos,
    CASE orientation WHEN 'N[M[' THEN ']M]N' WHEN ']M]N' THEN 'N[M[' ELSE orientation END, mate_state, mstrand FROM b_pairs")
run("ALTER TABLE b_records ADD COLUMN alt VARCHAR");
run("UPDATE b_records SET alt = CASE shape
  WHEN 'N[M[' THEN 'N[' || mate_chrom || ':' || mate_pos || '['
  WHEN 'N]M]' THEN 'N]' || mate_chrom || ':' || mate_pos || ']'
  WHEN ']M]N' THEN ']' || mate_chrom || ':' || mate_pos || ']N'
  ELSE '[' || mate_chrom || ':' || mate_pos || '[N' END")
run("CREATE TEMP TABLE b_events AS SELECT (row_number() OVER (ORDER BY f.ord, r.pos, r.id) - 1)::UBIGINT AS event_index, r.*
  FROM b_records r JOIN fa_order f ON f.chrom = r.chrom")

## ---- VCFs ----------------------------------------------------------------------------
vcf_header <- c("##fileformat=VCFv4.2",
  "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
  "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">",
  "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"Mate id\">",
  "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"Event id\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
write_query_vcf <- function(path, sql) {
  writeLines(c(vcf_header, dbGetQuery(con, sql)$line), path)
}
sv_vcf <- file.path(out, "sv.vcf")
write_query_vcf(sv_vcf, "SELECT format('{}\t{}\t{}\tN\t{}\t.\tPASS\tEND={};SVTYPE={}',
  chrom, CASE WHEN stype = 'INS' THEN s ELSE s - 1 END, format('sv{}', event_index), alt,
  CASE WHEN stype = 'INS' THEN s ELSE e END,
  CASE stype WHEN 'TDUP' THEN 'DUP' ELSE stype END) AS line FROM g_events ORDER BY event_index")
bnd_vcf <- file.path(out, "bnd.vcf")
write_query_vcf(bnd_vcf, "SELECT format('{}\t{}\t{}\tN\t{}\t.\tPASS\tSVTYPE=BND;MATEID={};EVENT=bnd{}', chrom, pos, id, alt,
  mate_id, pair_index) AS line FROM b_events ORDER BY event_index")
n_sv <- dbGetQuery(con, "SELECT count(*) n FROM g_events")$n
n_bnd <- dbGetQuery(con, "SELECT count(*) n FROM b_events")$n
cat("corpus: ", n_sv, " symbolic SV records, ", n_bnd, " BND records (", n_bnd / 2, " pairs)\n", sep = "")
stopifnot(n_sv > 0L, n_bnd > 0L)

## ---- VEP -----------------------------------------------------------------------------
vep <- function(vcf, json, ..., buffer = NULL) {
  if (file.exists(json)) file.remove(json)
  if (!is.null(buffer)) Sys.setenv(VEP_BUFFER_SIZE = buffer) else Sys.unsetenv("VEP_BUFFER_SIZE")
  rc <- system2("scripts/run_species_vep116_docker.sh", c(species, assembly, cache_version, cache, fasta, vcf, json, ...))
  stopifnot(identical(rc, 0L))
  lapply(readLines(json), fromJSON, simplifyVector = FALSE)
}
sv_json <- file.path(out, "sv-oracle.json")
bnd_json <- file.path(out, "bnd-oracle.json")
vep_sv <- vep(sv_vcf, sv_json, "--max_sv_size", "10000000")
vep_bnd <- vep(bnd_vcf, bnd_json, "--max_sv_size", "10000000", buffer = 1L)
stopifnot(length(vep_sv) == n_sv, length(vep_bnd) == n_bnd)

## ---- DuckVEP annotation ----------------------------------------------------------------
run("CREATE TEMP TABLE a_sv AS SELECT event_index, seq_region::UINTEGER AS seq_region,
  s::UBIGINT AS \"position\", 'N' AS reference, alt AS alternate, e::UBIGINT AS end_position,
  stype AS structural_type, copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
  FROM g_events")
run("CREATE TEMP TABLE a_bnd AS SELECT event_index, seq_region::UINTEGER AS seq_region, pos::UBIGINT AS \"position\",
  NULL::VARCHAR AS reference, alt AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
  NULL::VARCHAR AS copy_change, mate_seq_region::UINTEGER AS mate_seq_region, mate_pos::UBIGINT AS mate_position
  FROM b_events")
annotate <- function(table) dbGetQuery(con, sprintf(
  "SELECT a.event_index, n.tx, a.consequence, a.duckvep_status FROM query(duckvep_annotate_sql(%s, 'sv_species',
   struct_pack(rich := true, upstream_distance := 5000, downstream_distance := 5000))) a
   JOIN mnames n USING (transcript_index)", q(table)))
duck_sv <- annotate("a_sv"); duck_bnd <- annotate("a_bnd")
oracle_rows <- function(records, index_of) {
  rows <- lapply(records, function(r) lapply(r$transcript_consequences, function(t)
    data.frame(event_index = index_of(r$id), tx = t$transcript_id,
      terms = paste(sort(unique(unlist(t$consequence_terms))), collapse = "&"))))
  do.call(rbind, unlist(rows, recursive = FALSE))
}
o_sv <- oracle_rows(vep_sv, function(id) as.numeric(sub("^sv", "", id)))
o_bnd <- oracle_rows(vep_bnd, function(id) {
  idx <- match(id, dbGetQuery(con, "SELECT id FROM b_events ORDER BY event_index")$id)
  idx - 1
})
duck_sv$kind <- "sv"; duck_bnd$kind <- "bnd"; o_sv$kind <- "sv"; o_bnd$kind <- "bnd"
dbWriteTable(con, "duck_rows", rbind(duck_sv, duck_bnd), temporary = TRUE)
dbWriteTable(con, "oracle_rows", rbind(o_sv, o_bnd), temporary = TRUE)
# Canonical term sets per (event, transcript).
run("CREATE TEMP TABLE oracle_pairs AS SELECT kind, event_index::BIGINT AS event_index, tx,
    list_aggregate(list_sort(list_distinct(list(term))), 'string_agg', '&') AS consequence
  FROM (SELECT kind, event_index, tx, unnest(string_split(terms, '&')) AS term FROM oracle_rows) GROUP BY ALL")
run("CREATE TEMP TABLE duck_pairs AS SELECT kind, event_index::BIGINT AS event_index, tx,
    list_aggregate(list_sort(list_distinct(list(term))), 'string_agg', '&') AS consequence,
    string_agg(DISTINCT duckvep_status, '&') AS status
  FROM (SELECT kind, event_index, tx, duckvep_status, unnest(string_split(consequence, '&')) AS term FROM duck_rows) GROUP BY ALL")
run("CREATE TEMP TABLE differential AS SELECT coalesce(o.kind, d.kind) AS kind, coalesce(o.event_index, d.event_index) AS event_index,
  coalesce(o.tx, d.tx) AS tx, o.consequence AS oracle_terms, d.consequence AS duckvep_terms, d.status,
  CASE WHEN o.tx IS NULL THEN 'duckvep_only' WHEN d.tx IS NULL THEN 'oracle_only'
       WHEN o.consequence != d.consequence THEN 'terms_differ' ELSE 'exact' END AS verdict
  FROM oracle_pairs o FULL OUTER JOIN duck_pairs d USING (kind, event_index, tx)")
run("CREATE TEMP TABLE strata AS
  SELECT 'sv' AS kind, event_index::BIGINT AS event_index, stype AS sv_type, state, strand AS source_strand, format('sv{}', event_index) AS id FROM g_events
  UNION ALL SELECT 'bnd', event_index::BIGINT, 'BND', state || '/' || shape, strand, id FROM b_events")
run("CREATE TEMP TABLE annotated AS SELECT d.*, s.sv_type, s.state, s.source_strand, s.id, n.strand AS pair_strand, n.biotype
  FROM differential d JOIN strata s USING (kind, event_index) LEFT JOIN mnames n ON n.tx = d.tx")
stratify <- function(key, file) {
  res <- dbGetQuery(con, sprintf("SELECT %s, count(*) FILTER (verdict = 'exact') AS exact,
    count(*) FILTER (verdict = 'oracle_only') AS oracle_only, count(*) FILTER (verdict = 'duckvep_only') AS duckvep_only,
    count(*) FILTER (verdict = 'terms_differ') AS terms_differ, count(*) AS pairs,
    count(DISTINCT event_index || kind) AS events FROM annotated GROUP BY ALL ORDER BY ALL", key))
  write.csv(res, file.path(out, file), row.names = FALSE); res
}
by_type <- stratify("sv_type", "sv_type.csv")
by_type_strand <- stratify("sv_type, pair_strand", "sv_type_pair_strand.csv")
by_state <- stratify("sv_type, state", "sv_type_state.csv")
by_source_strand <- stratify("sv_type, source_strand", "sv_type_source_strand.csv")
by_biotype <- stratify("biotype", "biotype.csv")
overall <- dbGetQuery(con, "SELECT verdict, count(*) n FROM differential GROUP BY verdict ORDER BY verdict")
write.csv(overall, file.path(out, "verdicts.csv"), row.names = FALSE)
write.csv(dbGetQuery(con, "SELECT status, sv_type, count(*) n FROM (SELECT d.status, s.sv_type FROM duck_pairs d
  JOIN strata s USING (kind, event_index)) GROUP BY ALL ORDER BY ALL"), file.path(out, "duckvep_status.csv"), row.names = FALSE)
write.csv(dbGetQuery(con, "SELECT sv_type, term, verdict, count(*) n FROM (SELECT sv_type, verdict,
  unnest(string_split(coalesce(oracle_terms, duckvep_terms), '&')) AS term FROM annotated) GROUP BY ALL ORDER BY ALL"),
  file.path(out, "so_term.csv"), row.names = FALSE)
disagreements <- dbGetQuery(con, "SELECT * FROM annotated WHERE verdict <> 'exact' ORDER BY id, tx")
write.csv(disagreements, file.path(out, "disagreements.csv"), row.names = FALSE, na = "")
write.csv(dbGetQuery(con, "SELECT stype AS sv_type, state, strand, count(*) events FROM g_events GROUP BY ALL
  UNION ALL SELECT 'BND', state || '/' || shape, strand, count(*) FROM b_events GROUP BY ALL ORDER BY 1, 2, 3"),
  file.path(out, "corpus_manifest.csv"), row.names = FALSE)
write.csv(dbGetQuery(con, "SELECT 'sv' AS kind, id, chrom, s AS event_start, e AS event_end, stype, state, strand,
    transcript_index FROM (SELECT format('sv{}', event_index) AS id, * FROM g_events)
  UNION ALL SELECT 'bnd', id, chrom, pos, mate_pos, 'BND', state || '/' || shape, strand, NULL FROM b_events ORDER BY 1, 2"),
  file.path(out, "events.csv"), row.names = FALSE)

## ---- isolated re-runs of every disagreeing record ------------------------------------------
# VEP's BND interval tree is chromosome-blind and outlives the record it was built for, so a
# mate record placed at another record's mate coordinate can change the executable oracle's
# answer. Each record that disagrees in the pooled runs is therefore re-run alone (buffer 1)
# and both verdicts are retained; the pooled counterexamples are never removed.
bad <- unique(disagreements[c("kind", "id")])
isolated <- data.frame()
if (nrow(bad) > 0L) {
  stopifnot(nrow(bad) <= 60L)
  source_lines <- list(sv = readLines(sv_vcf), bnd = readLines(bnd_vcf))
  for (i in seq_len(nrow(bad))) {
    id <- bad$id[i]; kind <- bad$kind[i]
    one <- file.path(work, paste0("isolated-", id, ".vcf"))
    writeLines(c(grep("^#", source_lines[[kind]], value = TRUE),
      grep(paste0("\t", id, "\t"), source_lines[[kind]], value = TRUE)), one)
    rec <- vep(one, file.path(work, paste0("isolated-", id, ".json")), "--max_sv_size", "10000000", buffer = 1L)
    o <- do.call(rbind, lapply(rec[[1L]]$transcript_consequences, function(t) data.frame(
      tx = t$transcript_id, terms = paste(sort(unique(unlist(t$consequence_terms))), collapse = "&"))))
    ev <- dbGetQuery(con, paste("SELECT * FROM annotated WHERE id =", q(id)))
    d <- setNames(ev$duckvep_terms, ev$tx)
    d <- d[!is.na(d)]
    all_tx <- union(names(d), o$tx)
    same <- vapply(all_tx, function(tx) tx %in% names(d) && tx %in% o$tx &&
      identical(sort(strsplit(d[[tx]], "&")[[1L]]), sort(strsplit(o$terms[o$tx == tx][1L], "&")[[1L]])), TRUE)
    isolated <- rbind(isolated, data.frame(kind = kind, id = id, pooled_pairs = sum(ev$verdict != "exact"),
      isolated_pairs = length(all_tx), isolated_exact = sum(same), isolated_nonexact = sum(!same)))
  }
}
write.csv(isolated, file.path(out, "isolated_rerun.csv"), row.names = FALSE)

## ---- geometry and identity builders on the same records ----------------------------------
sv_in <- dbGetQuery(con, "SELECT event_index, CASE WHEN stype = 'INS' THEN s ELSE s - 1 END::DOUBLE AS pos, 'N' AS ref, alt,
  'END=' || CASE WHEN stype = 'INS' THEN s ELSE e END || ';SVTYPE=' || CASE stype WHEN 'TDUP' THEN 'DUP' ELSE stype END AS info
  FROM g_events ORDER BY event_index")
dbWriteTable(con, "sv_geometry_in", sv_in, temporary = TRUE)
geometry <- rduckvep_prepare_sv_geometry(con, "sv_geometry_in")
vep_span <- data.frame(event_index = as.numeric(sub("^sv", "", vapply(vep_sv, `[[`, "", "id"))),
  vep_start = vapply(vep_sv, function(r) as.numeric(r$start), 0), vep_end = vapply(vep_sv, function(r) as.numeric(r$end), 0))
g <- merge(geometry, vep_span, by = "event_index")
geometry_mismatch <- g[g$status != "ok" | g$nominal_start != g$vep_start | g$nominal_end != g$vep_end, ]
write.csv(geometry_mismatch, file.path(out, "geometry_disagreements.csv"), row.names = FALSE)
bnd_in <- dbGetQuery(con, "SELECT event_index, chrom, pos::DOUBLE AS pos, id, 'N' AS ref, alt,
  'SVTYPE=BND;MATEID=' || mate_id || ';EVENT=bnd' || pair_index AS info FROM b_events ORDER BY event_index")
dbWriteTable(con, "bnd_identity_in", bnd_in, temporary = TRUE)
identity <- rduckvep_prepare_breakend_pairs(con, "bnd_identity_in")
identity_reasons <- table(identity$reason)
write.csv(as.data.frame(identity_reasons), file.path(out, "bnd_identity_reasons.csv"), row.names = FALSE)

## ---- structural HGVS on literal-equivalent edits ------------------------------------------
# Small exact-span events around exon/CDS/intron landmarks (1-40 bp).
hg_states <- list(
  c("exon_start_cross", "x.es - 1 - CAST(hash(x.transcript_index, x.es, S, 1) % 5 AS BIGINT)"),
  c("exon_end_cross", "x.ee - CAST(hash(x.transcript_index, x.es, S, 1) % 5 AS BIGINT)"),
  c("exon_interior", "x.es + 3 + CAST(hash(x.transcript_index, x.es, S, 1) % GREATEST(x.ee - x.es - 50, 1) AS BIGINT)"),
  c("intron_interior", "x.ee + 4 + CAST(hash(x.transcript_index, x.es, S, 1) % GREATEST(x.next_es - x.ee - 50, 1) AS BIGINT)"))
hg_sql <- vapply(hg_states, function(s) sprintf("SELECT x.transcript_index, '%s' AS state, %s AS s,
  1 + CAST(hash(x.transcript_index, x.es, %d, 2) %% 40 AS BIGINT) AS len FROM g_ex x WHERE %s",
  s[1L], sub("\\bS\\b", as.character(seed), s[2L]), seed,
  if (s[1L] == "intron_interior") "x.next_es - x.ee > 60" else if (s[1L] == "exon_interior") "x.ee - x.es > 60" else "true"), "")
hg_sql <- c(hg_sql, sprintf("SELECT t.transcript_index, 'cds_start_cross' AS state,
  t.cs - CAST(hash(t.transcript_index, %d, 1) %% 5 AS BIGINT) AS s,
  1 + CAST(hash(t.transcript_index, %d, 2) %% 40 AS BIGINT) AS len FROM g_tx t WHERE t.cs IS NOT NULL", seed, seed))
run(paste0("CREATE TEMP TABLE h_geom AS ", paste(hg_sql, collapse = " UNION ALL ")))
run(sprintf("CREATE TEMP TABLE h_events0 AS SELECT * FROM (
  SELECT o.stype, o.alt, h.transcript_index, h.state, t.strand, t.chrom, t.seq_region, h.s AS pos, h.len,
    row_number() OVER (PARTITION BY o.stype, h.state, t.strand ORDER BY hash(h.transcript_index, h.s, h.len, o.stype, %d)) AS rk
  FROM h_geom h JOIN g_tx t USING (transcript_index) CROSS JOIN %s
  WHERE h.s >= 2 AND h.s + 2 * h.len + 1 < t.slen) WHERE rk <= %d", seed, types, as.integer(per)))
run("CREATE TEMP TABLE h_events AS SELECT (row_number() OVER (ORDER BY f.ord, e.pos, e.stype) - 1)::INTEGER AS event_index, e.*
  FROM h_events0 e JOIN fa_order f ON f.chrom = e.chrom")
hev <- dbGetQuery(con, "SELECT * FROM h_events ORDER BY event_index")
regions <- hev$chrom
fetch <- function(chrom, from, to) {
  spec <- paste0(chrom, ":", format(from, scientific = FALSE, trim = TRUE), "-", format(to, scientific = FALSE, trim = TRUE))
  listing <- file.path(work, "regions.txt"); writeLines(spec, listing)
  outp <- system2("samtools", c("faidx", "-n", "100000", "-r", listing, shQuote(fasta)), stdout = TRUE)
  stopifnot(length(outp) == 2L * length(spec))
  outp[seq(2L, length(outp), by = 2L)]
}
end <- hev$pos + hev$len
refseq <- toupper(fetch(hev$chrom, hev$pos, end + hev$len))
stopifnot(all(nchar(refseq) == 2L * hev$len + 1L))
n_ambig <- sum(grepl("[^ACGT]", refseq))
events <- data.frame(event_index = hev$event_index, chrom = hev$chrom, pos = hev$pos, ref = substr(refseq, 1L, 1L),
  alt = hev$alt, info = paste0("END=", end), stringsAsFactors = FALSE)
refs <- data.frame(event_index = hev$event_index, reference_sequence = refseq)
dbWriteTable(con, "sh_events", events, temporary = TRUE)
dbWriteTable(con, "sh_refs", refs, temporary = TRUE)
built <- rduckvep_prepare_structural_hgvs(con, "sh_events", "sh_refs")
stopifnot(nrow(built) == nrow(events), all(built$event_index == events$event_index))
ok <- built$hgvs_status == "supported"
hg_reason <- as.data.frame(table(status = built$hgvs_status, reason = built$hgvs_reason, useNA = "ifany"))
write.csv(hg_reason[hg_reason$Freq > 0, ], file.path(out, "hgvs_builder_statuses.csv"), row.names = FALSE)
supported <- built[ok, ]
lookup <- hev[match(supported$event_index, hev$event_index), ]
write_lit_vcf <- function(path, chrom, position, ref, alt, id) {
  ord <- order(match(chrom, fai$V1), position, id)
  writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(chrom[ord], format(position[ord], scientific = FALSE, trim = TRUE), id[ord], ref[ord], alt[ord], ".", "PASS", ".", sep = "\t")), path)
}
literal_vcf <- file.path(out, "hgvs-literal.vcf")
write_lit_vcf(literal_vcf, lookup$chrom, supported$literal_position, supported$literal_reference,
  supported$literal_alternate, paste0("e", supported$event_index))
vep0 <- vep(literal_vcf, file.path(out, "hgvs-shift0-oracle.json"), "--hgvs", "--hgvsg", "--shift_hgvs", "0")
vep1 <- vep(literal_vcf, file.path(out, "hgvs-shift1-oracle.json"), "--hgvs", "--hgvsg")
by_id <- function(records) { names(records) <- vapply(records, `[[`, "", "id"); records }
vep0 <- by_id(vep0); vep1 <- by_id(vep1)
field <- function(record, name) {
  v <- vapply(record$transcript_consequences, function(t) if (is.null(t[[name]])) NA_character_ else t[[name]], "")
  names(v) <- vapply(record$transcript_consequences, `[[`, "", "transcript_id"); v
}
genomic <- do.call(rbind, lapply(seq_len(nrow(supported)), function(i) {
  id <- paste0("e", supported$event_index[i])
  g0 <- unique(na.omit(field(vep0[[id]], "hgvsg"))); g1 <- unique(na.omit(field(vep1[[id]], "hgvsg")))
  data.frame(event_index = supported$event_index[i], sv_type = lookup$stype[i], state = lookup$state[i],
    strand = lookup$strand[i], builder = supported$hgvs_g[i], vep_shift0 = if (length(g0)) paste(g0, collapse = ",") else NA_character_,
    vep_shift1 = if (length(g1)) paste(g1, collapse = ",") else NA_character_, stringsAsFactors = FALSE)
}))
genomic$verdict <- ifelse(is.na(genomic$vep_shift0), "vep_no_hgvs",
  ifelse(genomic$builder == genomic$vep_shift0, "exact", "differs"))
genomic$shift1_differs <- !is.na(genomic$vep_shift1) & genomic$vep_shift1 != genomic$builder
write.csv(genomic, file.path(out, "hgvs_genomic.csv"), row.names = FALSE, na = "")
dbWriteTable(con, "sh_literal", data.frame(event_index = supported$event_index, chrom = lookup$chrom,
  position = supported$literal_position, reference = supported$literal_reference, alternate = supported$literal_alternate),
  temporary = TRUE)
run("CREATE TEMP TABLE sh_annotate AS SELECT l.event_index::UBIGINT AS event_index, r.seq_region::UINTEGER AS seq_region,
  l.\"position\"::UBIGINT AS \"position\", l.reference, l.alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
  NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
  FROM sh_literal l JOIN mregions r ON r.seq_region_name = l.chrom")
duck_h <- dbGetQuery(con, "SELECT a.event_index, n.tx, n.tx || coalesce('.' || n.version, '') || ':' || a.transcript_hgvs AS hgvsc
  FROM query(duckvep_annotate_sql('sh_annotate', 'sv_species', struct_pack(hgvs := true, upstream_distance := 5000,
  downstream_distance := 5000))) a JOIN mnames n USING (transcript_index) WHERE a.transcript_hgvs IS NOT NULL")
duck_lookup <- split(duck_h$hgvsc, duck_h$event_index)
transcript_rows <- lapply(seq_len(nrow(supported)), function(i) {
  id <- paste0("e", supported$event_index[i])
  vh <- field(vep1[[id]], "hgvsc"); vh <- vh[!is.na(vh)]
  actual <- duck_lookup[[as.character(supported$event_index[i])]]
  # DuckVEP strings carry tx.version:...; VEP hgvsc is already versioned.
  expected_s <- sort(unname(vh)); actual_s <- if (is.null(actual)) character() else sort(actual)
  data.frame(event_index = supported$event_index[i], sv_type = lookup$stype[i], state = lookup$state[i],
    strand = lookup$strand[i], vep_n = length(expected_s), duckvep_n = length(actual_s),
    exact = identical(expected_s, actual_s), vep_only = paste(setdiff(expected_s, actual_s), collapse = "|"),
    duckvep_only = paste(setdiff(actual_s, expected_s), collapse = "|"), stringsAsFactors = FALSE)
})
transcript <- do.call(rbind, transcript_rows)
write.csv(transcript, file.path(out, "hgvs_transcript.csv"), row.names = FALSE, na = "")
write.csv(transcript[!transcript$exact, ], file.path(out, "hgvs_transcript_disagreements.csv"), row.names = FALSE)
genomic$ok <- genomic$verdict == "exact"
transcript$ok <- transcript$exact
write.csv(aggregate(cbind(events = 1L, exact = ok) ~ sv_type + strand, data = transform(genomic, events = 1L), FUN = sum),
  file.path(out, "hgvs_genomic_by_type_strand.csv"), row.names = FALSE)
write.csv(aggregate(cbind(events = 1L, exact = ok) ~ sv_type + strand, data = transform(transcript, events = 1L), FUN = sum),
  file.path(out, "hgvs_transcript_by_type_strand.csv"), row.names = FALSE)

## ---- receipt ----------------------------------------------------------------------------
manifest <- data.frame(item = c("species", "assembly", "cache_version", "seed", "per_stratum", "model_sha256", "extension_sha256",
  "sv_vcf_sha256", "bnd_vcf_sha256", "sv_oracle_json_sha256", "bnd_oracle_json_sha256", "hgvs_literal_vcf_sha256",
  "hgvs_shift0_oracle_sha256", "hgvs_shift1_oracle_sha256", "sv_records", "bnd_records", "hgvs_events", "hgvs_ambiguous_reference_events"),
  value = c(species, assembly, cache_version, seed, per, model_sha, sha256(binary), sha256(sv_vcf), sha256(bnd_vcf),
    sha256(sv_json), sha256(bnd_json), sha256(literal_vcf), sha256(file.path(out, "hgvs-shift0-oracle.json")),
    sha256(file.path(out, "hgvs-shift1-oracle.json")), n_sv, n_bnd, nrow(hev), n_ambig))
write.csv(manifest, file.path(out, "manifest.csv"), row.names = FALSE)

## ---- summary and verdict ------------------------------------------------------------------
cat("\nconsequence-term pairs by verdict\n"); print(overall)
cat("\nby SV type\n"); print(by_type)
cat("\nDuckVEP status:\n"); print(dbGetQuery(con, "SELECT status, count(*) n FROM duck_pairs GROUP BY 1"))
cat("\ngeometry builder vs VEP nominal start/end mismatches:", nrow(geometry_mismatch), "of", nrow(g), "\n")
cat("BND identity reasons:", paste(names(identity_reasons), identity_reasons, collapse = ", "), "\n")
cat("HGVS builder statuses:\n"); print(table(built$hgvs_status, built$hgvs_reason, useNA = "ifany"))
cat("HGVS genomic (shift 0):", paste(names(table(genomic$verdict)), table(genomic$verdict), collapse = ", "),
    "; shift-1 counterexamples:", sum(genomic$shift1_differs), "\n")
cat("HGVS transcript strings exact events:", sum(transcript$exact), "of", nrow(transcript), "\n")
if (nrow(isolated)) { cat("\nisolated re-runs of pooled disagreements:\n"); print(isolated) }
exact_all <- (nrow(disagreements) == 0L || all(isolated$isolated_nonexact == 0L)) && nrow(geometry_mismatch) == 0L &&
  all(genomic$verdict != "differs") && all(transcript$exact)
quit(status = if (exact_all) 0L else 1L)
