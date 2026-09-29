library(tinytest)
library(DBI)

local({
  con <- rduckvep_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fixture <- system.file("extdata", "species", package = "Rduckvep")
  expect_true(nzchar(fixture))
  q <- function(x) as.character(dbQuoteString(con, x))
  products <- c("plasmodium", "tetrahymena", "human_grch37")
  names <- c("pf", "tetra", "human")
  for (i in seq_along(products)) {
    for (relation in c("regions", "transcripts")) {
      dbExecute(con, paste0("CREATE TABLE ", names[[i]], "_", relation,
                            " AS FROM read_parquet(",
                            q(file.path(fixture, paste0(products[[i]], "_", relation, ".parquet"))),
                            ")"))
    }
  }
  tables <- vapply(names, function(name) dbGetQuery(con,
    paste0("SELECT min(codon_table) n FROM ", name, "_transcripts"))$n, integer(1))
  expect_identical(unname(tables), c(4L, 6L, 1L))
  expect_equal(as.integer(dbGetQuery(con,
    "SELECT list_sort(list(transcript_index)) AS n FROM tetra_transcripts")$n[[1L]]),
    c(0L, 1L))

  load_model <- function(name, transcript = paste0(name, "_transcripts"),
                         reference_fasta = NULL) {
    regions <- if (name %in% c("tetra-standard", "tetra-pinned"))
      "tetra_regions" else paste0(name, "_regions")
    transcripts <- paste0("SELECT transcript_index, seq_region, transcript_start, transcript_end, ",
      "strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, ",
      "codon_table, pre_cds_sequence, post_cds_sequence FROM ", transcript,
      " ORDER BY seq_region, transcript_start, transcript_index")
    exons <- paste0("SELECT transcript_index, exon.exon_start, exon.exon_end, ",
      "exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM ",
      transcript, ", unnest(exons) u(exon) ORDER BY transcript_index, exon.exon_cdna_start")
    dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load(", q(name), ",",
      q(paste0("SELECT seq_region, sequence_length, seq_region_name FROM ", regions,
               " ORDER BY seq_region")),
      ",", q(transcripts), ",", q(exons),
      if (is.null(reference_fasta)) "" else paste0(",reference_fasta := ", q(reference_fasta)),
      ",transcript_coverage_complete := true)"))$loaded
  }
  for (name in names) expect_true(load_model(name))
  dbExecute(con, "CREATE TABLE tetra_standard AS SELECT * REPLACE (1::UTINYINT AS codon_table) FROM tetra_transcripts")
  expect_true(load_model("tetra-standard", "tetra_standard"))
  expect_true(load_model("tetra-pinned", "tetra_transcripts",
                         file.path(fixture, "tetrahymena.fa")))

  dbExecute(con, paste(
    "CREATE TABLE witness_events AS SELECT e.*, NULL::UBIGINT AS end_position,",
    "NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,",
    "NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM (VALUES",
    "(1::UBIGINT, 0::UINTEGER, 1241::UBIGINT, 'C', 'T'),",
    "(2::UBIGINT, (SELECT seq_region FROM tetra_transcripts WHERE transcript_stable_id='EAR80522'), 143::UBIGINT, 'T', 'C'),",
    "(3::UBIGINT, (SELECT seq_region FROM tetra_transcripts WHERE transcript_stable_id='EAR80553'), 295::UBIGINT, 'T', 'C')",
    ") e(event_index, seq_region, position, reference, alternate)"))
  dbExecute(con, "CREATE TABLE pf_events AS SELECT * FROM witness_events WHERE event_index=1")
  dbExecute(con, "CREATE TABLE tetra_events AS SELECT * FROM witness_events WHERE event_index IN (2,3)")
  dbExecute(con, paste(
    "CREATE TABLE human_events AS SELECT 4::UBIGINT AS event_index,",
    "0::UINTEGER AS seq_region, 1::UBIGINT AS position, 'A' AS reference,",
    "'C' AS alternate, NULL::UBIGINT AS end_position,",
    "NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,",
    "NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position"))
  annotate <- function(name) {
    events <- if (name == "pf") "pf_events" else if (name == "human")
      "human_events" else "tetra_events"
    dbGetQuery(con, paste0("SELECT a.event_index, s.consequence FROM query(duckvep_annotate_sql(",
      q(events), ",", q(name), ")) a JOIN duckvep_so_terms() s ON ",
      "(a.consequence_mask & s.consequence_mask) != 0 ORDER BY a.event_index"))
  }
  expect_identical(annotate("pf")$consequence[[1L]], "synonymous_variant")
  expect_identical(annotate("tetra")$consequence, rep("synonymous_variant", 2L))
  expect_identical(annotate("tetra-standard")$consequence, rep("stop_lost", 2L))
  expect_identical(annotate("tetra-pinned"), annotate("tetra"))
  human_before <- annotate("human")
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('pf') dropped")$dropped)
  expect_equal(nrow(annotate("tetra")), 2L)
  expect_identical(annotate("human"), human_before)

  receipt <- function(name, assembly) dbGetQuery(con, paste0(
    "SELECT model_sha256, transcript_count, assembly FROM query(duckvep_model_receipt_sql(",
    q(paste0(name, "_regions")), ",", q(paste0(name, "_transcripts")),
    ",'Ensembl','116',", q(assembly), ",repeat('a',64),repeat('b',64),",
    "'source-derived excerpt'))"))
  hashes <- vapply(seq_along(names), function(i) {
    record <- receipt(names[[i]], c("GCA000002765v3", "JCVI-TTA1-2.2", "GRCh37")[[i]])
    expect_equal(record$transcript_count, c(1L, 2L, 1L)[[i]])
    record$model_sha256
  }, character(1))
  expect_equal(length(unique(hashes)), 3L)
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('tetra-standard') dropped")$dropped)
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('tetra-pinned') dropped")$dropped)
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('tetra') dropped")$dropped)
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('human') dropped")$dropped)
})
