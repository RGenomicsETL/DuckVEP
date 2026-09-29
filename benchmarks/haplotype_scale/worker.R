#!/usr/bin/env Rscript
# One fresh process of the #2 slice 7 scale qualification: model load, then mode A or B through the public
# builder, materialized to Parquet. Appends one TSV receipt row per pass (cold, warm) to --receipt.
#
#   Rscript benchmarks/haplotype_scale/worker.R --extension EXT --model full|mane --mode A|B|stage \
#     --calls CALLS.parquet | --vcf INPUT.vcf.gz --out OUT.parquet --receipt R.tsv --label LABEL [options]
#
# Modes
#   B      the unsorted VCF text: decode (DuckDB read_csv of the bgzipped VCF), transcript discovery
#          (duckvep_coding_transcripts), call construction and staging, then duckvep_haplotypes into Parquet.
#   A      the preordered calls relation (--calls, written by --mode stage) staged from Parquet, then
#          duckvep_haplotypes into Parquet. The builder's own ordering sort is included.
#   stage  mode B's staging only, writing the ordered calls to --calls (untimed input preparation for mode A).
# Options: --warm (a second pass in the same process after the cold pass), --native-budget BYTES (default 4 GiB),
# --memory-limit (default 8GB), --spill DIR, --max-spill (default 32GiB), --calls-filter SQL (restrict calls in
# mode A, for ad-hoc partitions), --passes N (one cold pass and N - 1 warm ones), --order-load.
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NA_character_) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) return(default)
  value <- args[i + 1L]; args <<- args[-c(i, i + 1L)]; value
}
flag <- function(name) { i <- match(paste0("--", name), args); if (is.na(i)) FALSE else { args <<- args[-i]; TRUE } }
extension <- normalizePath(opt("extension"), mustWork = TRUE)
model_kind <- opt("model", "full")
mode <- opt("mode", "B")
vcf <- opt("vcf"); calls_path <- opt("calls"); out <- opt("out"); receipt <- opt("receipt"); label <- opt("label", mode)
warm <- flag("warm")
passes <- as.integer(opt("passes", if (warm) "2" else "1"))   # 1 cold pass, then passes - 1 warm passes
budget_bytes <- opt("native-budget", "4294967296")
memory_limit <- opt("memory-limit", "8GB"); max_spill <- opt("max-spill", "32GiB")
spill <- opt("spill", file.path(tempdir(), "duckvep-spill")); calls_filter <- opt("calls-filter", "")
audit <- !flag("no-audit")
model_path <- opt("model-file", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
discovery <- opt("discovery", "scalar")   # scalar: duckvep_coding_transcripts; annotate: the annotation builder (region_mask & 16), the pre-scalar method, for the A/B
decoder <- opt("decoder", "duckdb")   # duckdb: read_csv inflates the .vcf.gz itself; bgzip: `bgzip -dc` (htslib, libdeflate) feeds it through a FIFO
order_load <- flag("order-load")   # add ORDER BY to the model queries of the full model (the stored model is already in load order)
alignment_cells <- opt("max-alignment-cells", "268435456"); workspace <- opt("workspace-limit", "1073741824")
stopifnot(mode %in% c("A", "B", "stage"), model_kind %in% c("full", "mane"), decoder %in% c("duckdb", "bgzip"), discovery %in% c("scalar", "annotate"))

con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true", threads = "1")))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
q <- function(x) as.character(dbQuoteString(con, x))
sql <- function(s) invisible(dbExecute(con, s))
get <- function(s) dbGetQuery(con, s)
elapsed <- function(code) unname(system.time(force(code))[["elapsed"]])
status <- function(key) {
  line <- grep(paste0("^", key, ":"), readLines("/proc/self/status", warn = FALSE), value = TRUE)
  as.numeric(sub(paste0("^", key, ":[[:space:]]*([0-9]+).*"), "\\1", line))
}
sql(paste("LOAD", q(extension)))
dir.create(spill, recursive = TRUE, showWarnings = FALSE)
sql(paste("SET memory_limit =", q(memory_limit)))
sql(paste("SET temp_directory =", q(spill)))
sql(paste("SET max_temp_directory_size =", q(max_spill)))
sql(paste0("ATTACH ", q(model_path), " AS m (READ_ONLY)"))
sql(sprintf("SELECT duckvep_native_budget_set(%s)", budget_bytes))

model_load <- function() {
  if (model_kind == "mane") {
    # MANE Select subset: the model's own MANE annotation, re-indexed densely in (seq_region, start, index) order.
    # Both the transcript and the dependent exon/miRNA/peptide queries go through the same map.
    sql("CREATE TABLE mane_map AS SELECT transcript_index AS source_index, (row_number() OVER (ORDER BY seq_region,
      transcript_start, transcript_index) - 1)::UINTEGER AS transcript_index FROM m.model_transcripts
      WHERE mane_select_refseq IS NOT NULL")
    transcripts <- "SELECT k.transcript_index, t.seq_region, t.transcript_start, t.transcript_end, t.strand, t.gene_index,
      t.transcript_flags, t.cds_start, t.cds_end, t.cds_sequence, t.codon_table, t.pre_cds_sequence, t.post_cds_sequence
      FROM m.model_transcripts t JOIN mane_map k ON k.source_index = t.transcript_index ORDER BY k.transcript_index"
    exons <- "SELECT k.transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase
      FROM m.model_transcripts t JOIN mane_map k ON k.source_index = t.transcript_index, unnest(t.exons) u(e)
      ORDER BY k.transcript_index, e.exon_cdna_start"
    mirna <- "SELECT k.transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM m.model_transcripts t
      JOIN mane_map k ON k.source_index = t.transcript_index, unnest(t.mature_mirna_regions) u(x) ORDER BY 1, 2"
    edits <- "SELECT k.transcript_index, x.protein_position, x.alternate_amino_acid FROM m.model_transcripts t
      JOIN mane_map k ON k.source_index = t.transcript_index, unnest(t.peptide_edits) u(x) ORDER BY 1, 2"
  } else {
    # The compiled model is stored in load order (transcript_index dense, sorted by region and start, exons by cDNA
    # position), and the loader rejects rows that arrive out of order, so ORDER BY (a sort of 1.1 GB of sequence) is
    # optional here; --order-load restores it.
    by <- function(x) if (order_load) paste("ORDER BY", x) else ""
    transcripts <- paste("SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index,
      transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence
      FROM m.model_transcripts", by("seq_region, transcript_start, transcript_index"))
    exons <- paste("SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start, e.exon_cdna_end, e.phase, e.end_phase
      FROM m.model_transcripts, unnest(exons) u(e)", by("transcript_index, e.exon_cdna_start"))
    mirna <- paste("SELECT transcript_index, x.mature_mirna_start, x.mature_mirna_end FROM m.model_transcripts,
      unnest(mature_mirna_regions) u(x)", by("1, 2"))
    edits <- paste("SELECT transcript_index, x.protein_position, x.alternate_amino_acid FROM m.model_transcripts,
      unnest(peptide_edits) u(x)", by("1, 2"))
  }
  sql("CREATE TABLE regions AS SELECT seq_region::BIGINT AS seq_region, seq_region_name, sequence_length FROM m.model_regions")
  get(paste0("SELECT loaded FROM duckvep_model_load('hap', 'SELECT seq_region::UINTEGER AS seq_region, sequence_length
    FROM regions ORDER BY seq_region', ", q(transcripts), ", ", q(exons), ",
    mature_mirna_query := ", q(mirna), ", peptide_edit_query := ", q(edits), ", transcript_coverage_complete := TRUE)"))
}

header_lines <- function(path) {
  n <- as.integer(system2("sh", c("-c", shQuote(paste0("zcat ", shQuote(path), " | head -n 5000 | grep -c '^#'"))), stdout = TRUE))
  stopifnot(is.finite(n), n > 0L)
  n
}
vcf_columns <- "columns={'chrom':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','qual':'VARCHAR',
  'filter':'VARCHAR','info':'VARCHAR','fmt':'VARCHAR','sample':'VARCHAR'}"
vcf_source <- function(path) {
  skip <- header_lines(path)
  if (decoder == "duckdb")
    return(sprintf("read_csv(%s, delim='\\t', header=false, skip=%d, auto_detect=false, quote='', escape='', strict_mode=false, %s, compression='gzip')",
      q(path), skip, vcf_columns))
  # bgzip -dc inflates with libdeflate several times faster than the miniz inside DuckDB. It runs on the same core (the
  # affinity is inherited), so this only changes the inflate implementation, not the cores.
  fifo <- tempfile("vcf-", fileext = ".fifo", tmpdir = spill)
  stopifnot(system2("mkfifo", shQuote(fifo)) == 0L)
  system(sprintf("bgzip -dc %s > %s", shQuote(path), shQuote(fifo)), wait = FALSE)
  sprintf("read_csv(%s, delim='\\t', header=false, skip=%d, auto_detect=false, quote='', escape='', strict_mode=false, %s, compression='none')",
    q(fifo), skip, vcf_columns)
}

stage_b <- function() {
  sql("DROP TABLE IF EXISTS cand")
  if (discovery == "scalar") {
    # One pass over the VCF text: each ALT allele is an event, and only events that overlap the CDS of a model transcript survive the scalar.
    sql(paste0("CREATE TABLE cand AS SELECT * FROM (SELECT v.record_index, a.i, r.seq_region, v.pos, v.ref, a.alt AS alt, v.sample,
      duckvep_coding_transcripts('hap', r.seq_region, v.pos, v.ref, a.alt) AS tx
      FROM (SELECT row_number() OVER () AS record_index, chrom, pos, ref, alt, sample FROM ", vcf_source(vcf), ") v
      JOIN regions r ON r.seq_region_name = v.chrom, unnest(string_split(v.alt, ',')) WITH ORDINALITY a(alt, i) WHERE v.alt <> '.') WHERE len(tx) > 0"))
  } else {
    # The pre-scalar method: every ALT allele becomes an event, the annotation builder lists the transcripts each event
    # overlaps (twelve per intronic event), and those whose CDS bit is set are kept.
    sql("DROP TABLE IF EXISTS ev_all"); sql("DROP VIEW IF EXISTS ordered_events_all")
    sql(paste0("CREATE TABLE ev_all AS SELECT ((v.record_index << 6) | (a.i - 1))::UBIGINT AS event_index, v.record_index, r.seq_region::UINTEGER AS seq_region,
      v.pos::UBIGINT AS position, v.ref AS reference, a.alt AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
      NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position, v.sample
      FROM (SELECT row_number() OVER () AS record_index, chrom, pos, ref, alt, sample FROM ", vcf_source(vcf), ") v
      JOIN regions r ON r.seq_region_name = v.chrom, unnest(string_split(v.alt, ',')) WITH ORDINALITY a(alt, i) WHERE v.alt <> '.'"))
    sql("CREATE VIEW ordered_events_all AS SELECT event_index, seq_region, position, reference, alternate, end_position, structural_type, copy_change, mate_seq_region, mate_position
      FROM ev_all ORDER BY seq_region, position, event_index")
    sql("DROP TABLE IF EXISTS calls")
    sql("CREATE TABLE calls AS SELECT e.event_index::BIGINT AS event_index, e.seq_region::INTEGER AS seq_region, e.position::BIGINT AS position, e.reference,
      e.alternate, (e.event_index & 63)::INTEGER + 1 AS alt_index, h.transcript_index::INTEGER AS transcript_index, 0::INTEGER AS sample_index,
      list_transform(string_split_regex(split_part(e.sample, ':', 1), '[/|]'), lambda x: try_cast(x AS INTEGER)) AS alleles,
      [false] || list_transform(regexp_extract_all(split_part(e.sample, ':', 1), '[/|]'), lambda s: s = '|') AS phase_before,
      try_cast(list_last(string_split(e.sample, ':')) AS BIGINT) AS phase_set
      FROM ev_all e JOIN (SELECT DISTINCT event_index, transcript_index FROM query(duckvep_annotate_sql('ordered_events_all', 'hap')) WHERE region_mask & 16 <> 0) h USING (event_index)
      ORDER BY seq_region, position, event_index, transcript_index")
    sql("DROP TABLE ev_all")
    return(invisible())
  }
  sql("DROP TABLE IF EXISTS calls")
  # One call per (ALT allele, transcript); GT '1|0', '0/1', '1|2' become alleles and per-lane phase flags. PS is the
  # phase set only when it is an integer; labels such as PATMAT and HOMVAR mean no phase set (whole-domain phase).
  sql("CREATE TABLE calls AS SELECT ((c.record_index << 6) | (c.i - 1))::BIGINT AS event_index, c.seq_region::INTEGER AS seq_region,
    c.pos::BIGINT AS position, c.ref AS reference, c.alt AS alternate, c.i::INTEGER AS alt_index, t.transcript_index::INTEGER AS transcript_index,
    0::INTEGER AS sample_index, list_transform(string_split_regex(split_part(c.sample, ':', 1), '[/|]'), lambda x: try_cast(x AS INTEGER)) AS alleles,
    [false] || list_transform(regexp_extract_all(split_part(c.sample, ':', 1), '[/|]'), lambda s: s = '|') AS phase_before,
    try_cast(list_last(string_split(c.sample, ':')) AS BIGINT) AS phase_set
    FROM cand c, unnest(c.tx) t(transcript_index)
    WHERE CASE WHEN c.i > 64 THEN error('more than 64 ALT alleles in one record') ELSE true END
    ORDER BY seq_region, position, event_index, transcript_index")
}

stage_a <- function() {
  sql("DROP TABLE IF EXISTS calls")
  where <- if (nzchar(calls_filter)) paste("WHERE", calls_filter) else ""
  sql(paste0("CREATE TABLE calls AS SELECT * FROM read_parquet(", q(calls_path), ") ", where))
}

predict <- function() {
  if (file.exists(out)) file.remove(out)
  partial <- paste0(out, ".partial")   # COPY creates its target when it starts: publish by rename after success only
  on.exit(if (file.exists(partial)) file.remove(partial), add = TRUE)
  sql(sprintf("COPY (SELECT * FROM duckvep_haplotypes('SELECT * FROM calls', 'hap', max_alignment_cells := %s, workspace_limit := %s))
    TO %s (FORMAT parquet)", alignment_cells, workspace, q(partial)))
  stopifnot(file.rename(partial, out))
}

budget_snapshot <- function(tag) {
  b <- get("SELECT owner, current_bytes, high_water_bytes FROM duckvep_native_budget()")
  setNames(c(b$current_bytes, b$high_water_bytes), c(paste0(tag, "_cur_", b$owner), paste0(tag, "_hw_", b$owner)))
}
spill_bytes <- function() {
  x <- tryCatch(get("SELECT coalesce(sum(size), 0)::BIGINT AS b FROM duckdb_temporary_files()")$b, error = function(e) NA)
  x
}

t_load <- elapsed(model_load())
load_snapshot <- budget_snapshot("load")
rss_after_load <- status("VmHWM")
row_header <- c("label", "mode", "model", "pass", "load_s", "stage_s", "predict_s", "pass_s", "total_s", "rss_hwm_kib", "duckdb_spill_bytes")
if (!is.na(receipt) && !file.exists(receipt)) {
  cat(paste(c(row_header, names(load_snapshot), "audit"), collapse = "\t"), "\n", file = receipt, sep = "")
}
metrics <- list()
pass_index <- 0L
for (pass in c("cold", rep("warm", passes - 1L))) {
  pass_index <- pass_index + 1L
  sql("SELECT duckvep_native_budget_reset_high_water()")
  if (mode == "stage") {
    t_stage <- elapsed(stage_b()); t_predict <- NA
    elapsed(sql(sprintf("COPY (SELECT * FROM calls) TO %s (FORMAT parquet)", q(calls_path))))
  } else {
    t_stage <- elapsed(if (mode == "B") stage_b() else stage_a())
    t_predict <- elapsed(predict())
  }
  snap <- budget_snapshot("pass")
  pass_s <- sum(c(t_stage, t_predict), na.rm = TRUE)
  total <- if (pass == "cold") t_load + pass_s else pass_s
  audit_text <- ""
  if (audit && mode != "stage" && pass_index == passes) {
    a <- get(sprintf("SELECT count(*) AS rows, coalesce(sum(len(carriers)), 0) AS carriers,
      coalesce(sum(length(cds)), 0) AS translated_bases, sum(hash(t)::HUGEINT)::VARCHAR AS checksum,
      count(*) FILTER (WHERE prediction_status = 'predicted') AS predicted FROM read_parquet(%s) t", q(out)))
    c <- get("SELECT count(*) AS calls, count(DISTINCT (event_index, transcript_index)) AS projections,
      count(DISTINCT event_index) AS events, count(DISTINCT (seq_region, position, reference)) AS sources,
      count(DISTINCT transcript_index) AS transcripts FROM calls")
    peaks <- get("SELECT max(n) AS max_events_per_transcript FROM (SELECT count(DISTINCT event_index) AS n FROM calls GROUP BY transcript_index)")
    statuses <- get(sprintf("SELECT string_agg(prediction_status || '/' || prediction_reason || ':' || n, '|' ORDER BY n DESC) AS status_counts
      FROM (SELECT prediction_status, prediction_reason, count(*) AS n FROM read_parquet(%s) GROUP BY ALL)", q(out)))
    records <- get("SELECT count(DISTINCT event_index >> 6) AS records_in_scope FROM calls")
    audit_text <- paste(sprintf("%s=%s", c(names(a), names(c), names(peaks), names(records), names(statuses), "output_bytes"),
      c(unlist(lapply(a, as.character)), unlist(lapply(c, as.character)), unlist(lapply(peaks, as.character)),
        unlist(lapply(records, as.character)), unlist(lapply(statuses, as.character)), file.size(out))), collapse = ";")
  }
  row <- c(label, mode, model_kind, pass, sprintf("%.3f", if (pass == "cold") t_load else 0), sprintf("%.3f", t_stage),
    sprintf("%.3f", t_predict), sprintf("%.3f", pass_s), sprintf("%.3f", total), sprintf("%.0f", status("VmHWM")),
    sprintf("%s", spill_bytes()))
  line <- paste(c(row, sprintf("%.0f", load_snapshot), audit_text), collapse = "\t")
  if (!is.na(receipt)) cat(line, "\n", file = receipt, sep = "", append = TRUE)
  cat(paste0("PASS\t", pass, "\tload=", sprintf("%.3f", t_load), "\tstage=", sprintf("%.3f", t_stage), "\tpredict=",
    sprintf("%.3f", t_predict), "\n"))
  cat("BUDGET\t", paste(sprintf("%s=%.1f", names(snap), snap / 1048576), collapse = ";"), "\n", sep = "")
  if (nzchar(audit_text)) cat("AUDIT\t", audit_text, "\n", sep = "")
}
