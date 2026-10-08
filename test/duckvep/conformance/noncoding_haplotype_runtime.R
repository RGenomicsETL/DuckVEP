run_full_noncoding_producer <- function(connection, sql_quote) {
  execute <- function(sql) DBI::dbExecute(connection, sql)
  execute("CREATE SCHEMA nc_core")
  execute("CREATE TABLE nc_reference(chrom VARCHAR, start BIGINT, \"end\" BIGINT, seq VARCHAR)")
  execute("INSERT INTO nc_reference VALUES ('1', 0, 40, 'NNNNNNNNNNAAAAACCCCCNNNNNNNNNNGGGGGTTTTT')")
  execute("CREATE TABLE nc_core.coord_system(coord_system_id BIGINT, species_id BIGINT, name VARCHAR, version VARCHAR, rank BIGINT)")
  execute("INSERT INTO nc_core.coord_system VALUES (1, 1, 'chromosome', 'NC', 1)")
  execute("CREATE TABLE nc_core.seq_region(seq_region_id BIGINT, name VARCHAR, coord_system_id BIGINT, length BIGINT)")
  execute("INSERT INTO nc_core.seq_region VALUES (1, '1', 1, 40)")
  execute("CREATE TABLE nc_core.seq_region_attrib(seq_region_id BIGINT, attrib_type_id BIGINT, value BIGINT)")
  execute("CREATE TABLE nc_core.gene(gene_id BIGINT, biotype VARCHAR, is_current BIGINT, stable_id VARCHAR, version BIGINT)")
  execute("INSERT INTO nc_core.gene VALUES (1, 'lncRNA', 1, 'ENSG_NC', 1)")
  execute(paste0("CREATE TABLE nc_core.transcript(transcript_id BIGINT, gene_id BIGINT, seq_region_id BIGINT, ",
    "seq_region_start BIGINT, seq_region_end BIGINT, seq_region_strand BIGINT, biotype VARCHAR, ",
    "is_current BIGINT, stable_id VARCHAR, version BIGINT)"))
  execute("CREATE TABLE nc_core.exon(exon_id BIGINT, seq_region_id BIGINT, seq_region_start BIGINT, seq_region_end BIGINT, seq_region_strand BIGINT, phase BIGINT, end_phase BIGINT, is_current BIGINT, stable_id VARCHAR)")
  execute("CREATE TABLE nc_core.exon_transcript(exon_id BIGINT, transcript_id BIGINT, rank BIGINT)")
  execute("CREATE TABLE nc_core.translation(translation_id BIGINT, transcript_id BIGINT, seq_start BIGINT, start_exon_id BIGINT, seq_end BIGINT, end_exon_id BIGINT, stable_id VARCHAR, version BIGINT)")
  execute("CREATE TABLE nc_core.attrib_type(attrib_type_id BIGINT, code VARCHAR)")
  execute("CREATE TABLE nc_core.transcript_attrib(transcript_id BIGINT, attrib_type_id BIGINT, value VARCHAR)")
  execute("CREATE TABLE nc_core.translation_attrib(translation_id BIGINT, attrib_type_id BIGINT, value VARCHAR)")
  make_calls <- function(rows) {
    execute("DROP TABLE IF EXISTS nc_calls")
    execute(paste0("CREATE TABLE nc_calls AS SELECT * FROM (VALUES ", paste(rows, collapse = ","),
      ") t(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,",
      "sample_index,phase_set,alleles,phase_before)"))
  }
  for (strand in c(1L, -1L)) {
    suffix <- if (strand == 1L) "plus" else "minus"
    model <- paste0("full_noncoding_", suffix)
    execute("DELETE FROM nc_core.transcript")
    execute("DELETE FROM nc_core.exon")
    execute("DELETE FROM nc_core.exon_transcript")
    execute(sprintf("INSERT INTO nc_core.transcript VALUES (1, 1, 1, 11, 40, %d, 'lncRNA', 1, 'ENST_NC_%s', 1)", strand, suffix))
    if (strand == 1L) {
      execute("INSERT INTO nc_core.exon VALUES (11, 1, 11, 20, 1, -1, 0, 1, 'ENSE_NC_1'), (12, 1, 31, 40, 1, 0, 0, 1, 'ENSE_NC_2')")
      execute("INSERT INTO nc_core.exon_transcript VALUES (11, 1, 1), (12, 1, 2)")
      valid <- c("(1,0,16,'C','A',1,0,0,12,[1,1],[false,true])", "(2,0,36,'T','C',1,0,0,12,[1,1],[false,true])")
      splice <- "(1,0,18,'C','A',1,0,0,12,[1,1],[false,true])"
      external_edge <- "(1,0,13,'A','C',1,0,0,12,[1,1],[false,true])"
      conflict <- c("(1,0,16,'C','A',1,0,0,12,[1,1],[false,true])", "(2,0,16,'C','G',1,0,0,12,[1,1],[false,true])")
    } else {
      execute("INSERT INTO nc_core.exon VALUES (11, 1, 31, 40, -1, -1, 0, 1, 'ENSE_NC_1'), (12, 1, 11, 20, -1, 0, 0, 1, 'ENSE_NC_2')")
      execute("INSERT INTO nc_core.exon_transcript VALUES (11, 1, 1), (12, 1, 2)")
      valid <- c("(1,0,35,'G','T',1,0,0,12,[1,1],[false,true])", "(2,0,15,'A','G',1,0,0,12,[1,1],[false,true])")
      splice <- "(1,0,33,'G','T',1,0,0,12,[1,1],[false,true])"
      external_edge <- "(1,0,38,'T','G',1,0,0,12,[1,1],[false,true])"
      conflict <- c("(1,0,35,'G','T',1,0,0,12,[1,1],[false,true])", "(2,0,35,'G','C',1,0,0,12,[1,1],[false,true])")
    }
    region_sql <- DBI::dbGetQuery(connection, "SELECT duckvep_ensembl_regions_sql('nc_core', 'nc_reference', 'NC') AS query_text")$query_text[[1L]]
    execute(paste0("CREATE OR REPLACE TABLE nc_regions AS FROM query(", sql_quote(region_sql), ")"))
    transcript_sql <- DBI::dbGetQuery(connection, "SELECT duckvep_ensembl_transcripts_sql('nc_core', 'nc_reference', 'NC') AS query_text")$query_text[[1L]]
    execute(paste0("CREATE OR REPLACE TABLE nc_transcripts AS FROM query(", sql_quote(transcript_sql), ")"))
    prepared <- DBI::dbGetQuery(connection, "SELECT hex(cdna_sequence) cdna, cds_sequence FROM nc_transcripts")
    stopifnot(nrow(prepared) == 1L, prepared$cdna == "4141414141434343434347474747475454545454", is.na(prepared$cds_sequence[[1L]]))
    transcript_query <- "SELECT transcript_index,seq_region,transcript_start,transcript_end,strand,gene_index,transcript_flags,cds_start,cds_end,cds_sequence,codon_table,pre_cds_sequence,post_cds_sequence,cdna_sequence FROM nc_transcripts"
    exon_query <- paste0("SELECT transcript_index,e.exon_start,e.exon_end,e.exon_cdna_start,e.exon_cdna_end,e.phase,e.end_phase ",
      "FROM nc_transcripts,unnest(exons) x(e)")
    load <- function(name, query = transcript_query) DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(",
      sql_quote(name), ",", sql_quote("SELECT seq_region,sequence_length FROM nc_regions"), ",", sql_quote(query), ",", sql_quote(exon_query), ")"))
    stopifnot(load(model)$loaded[[1L]] == 1L)
    relation <- paste0("duckvep_haplotypes(", sql_quote("SELECT * FROM nc_calls ORDER BY seq_region,position,event_index"), ",", sql_quote(model), ")")
    result <- function() DBI::dbGetQuery(connection, paste0("SELECT cds,protein,prediction_status,array_to_string(haplotype_consequences,',') consequences,edit_count,",
      "array_to_string(list_sort(list_transform(contributors,lambda c: c.event_index)),',') source_ids FROM ", relation))
    make_calls(valid)
    observed <- result()
    stopifnot(nrow(observed) == 1L, observed$prediction_status == "predicted",
      observed$consequences == "non_coding_transcript_exon_variant", is.na(observed$cds[[1L]]),
      is.na(observed$protein[[1L]]), observed$edit_count == 0L, observed$source_ids == "1,2")
    capacity <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ", sub("\\)$", ",max_sequence_bases:=19)", relation)))
    stopifnot(nrow(capacity) == 1L, capacity$prediction_status == "unsupported_context")
    make_calls(conflict)
    stopifnot(result()$prediction_status == "unsupported_context")
    make_calls(splice)
    stopifnot(result()$prediction_status == "unsupported_context")
    make_calls(external_edge)
    stopifnot(result()$prediction_status == "predicted")
    make_calls(valid)
    stopifnot(result()$prediction_status == "predicted")
    uncertain <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM duckvep_haplotypes(",
      sql_quote("SELECT * REPLACE ([0,1]::INTEGER[] AS alleles, [false,false]::BOOLEAN[] AS phase_before) FROM nc_calls"), ",", sql_quote(model), ")"))
    stopifnot(nrow(uncertain) > 0L, all(uncertain$prediction_status != "predicted"))
    old_model <- paste0(model, "_old")
    old_query <- sub(",cdna_sequence FROM nc_transcripts$", " FROM nc_transcripts", transcript_query)
    stopifnot(load(old_model, old_query)$loaded[[1L]] == 1L)
    old_result <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ", sub(sql_quote(model), sql_quote(old_model), relation, fixed = TRUE)))
    stopifnot(nrow(old_result) == 1L, old_result$prediction_status == "unsupported_context")
    bad_query <- sub("cdna_sequence FROM nc_transcripts$", "substr(CAST(cdna_sequence AS VARCHAR),1,19)::BLOB cdna_sequence FROM nc_transcripts", transcript_query)
    bad <- try(load(paste0(model, "_bad"), bad_query), silent = TRUE)
    stopifnot(inherits(bad, "try-error"))
    for (name in c(model, old_model)) DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", sql_quote(name), ")"))
  }
}

main <- function(root) {
  fixtures <- list(
    utr5 = list(term = "5_prime_UTR_variant", positions = c(2L, 4L), alternate = "C", other = "utr3"),
    utr3 = list(term = "3_prime_UTR_variant", positions = c(17L, 19L), alternate = "T", other = "utr5"))
  axes <- list(
    list(name = "plus", strand = 1L, cds_start = 106L, spliced_start = 112L,
      exons = "(100,103,1,4),(110,125,5,20)"),
    list(name = "minus", strand = -1L, cds_start = 105L, spliced_start = 105L,
      exons = "(122,125,1,4),(100,115,5,20)"))
  complement <- c(A = "T", C = "G", G = "C", T = "A")
  driver <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE)
  connection <- DBI::dbConnect(driver)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  sql_quote <- function(x) as.character(DBI::dbQuoteString(connection, x))
  extension <- Sys.getenv("DUCKVEP_EXTENSION",
    file.path(root, "build/release/duckvep.duckdb_extension"))
  DBI::dbExecute(connection, paste("LOAD", sql_quote(normalizePath(extension, mustWork = TRUE))))
  checked <- 0L
  for (axis in axes) for (feature in names(fixtures)) {
    strand <- axis$strand
    fixture <- fixtures[[feature]]
    model <- paste0("noncoding_", feature, "_", axis$name)
    cds_start <- axis$cds_start
    cds_end <- cds_start + 8L
    transcript <- sprintf(paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
      "100::UBIGINT transcript_start,119::UBIGINT transcript_end,%d::TINYINT strand,",
      "0::UINTEGER gene_index,3::UBIGINT transcript_flags,%d::UBIGINT cds_start,",
      "%d::UBIGINT cds_end,'ATGAAATAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
      "'AAAAAA'::BLOB pre_cds_sequence,'CCCCC'::BLOB post_cds_sequence"), strand, cds_start, cds_end)
    exons <- paste0("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,",
      "119::UBIGINT exon_end,1::UBIGINT exon_cdna_start,20::UBIGINT exon_cdna_end,",
      "0::TINYINT phase,0::TINYINT end_phase")
    loaded <- DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(",
      sql_quote(model), ",", sql_quote("SELECT 0::UINTEGER seq_region"), ",",
      sql_quote(transcript), ",", sql_quote(exons), ")"))
    stopifnot(loaded$loaded[[1L]] == 1L)
    positions <- fixture$positions
    reference <- "AAAAAAATGAAATAACCCCC"
    alternate <- strsplit(reference, "", fixed = TRUE)[[1L]]
    refs <- alternate[positions]
    alts <- rep(fixture$alternate, 2L)
    alternate[positions] <- alts
    changed <- which(strsplit(reference, "", fixed = TRUE)[[1L]] != alternate)
    stopifnot(identical(changed, positions), all(changed < 7L | changed > 15L))
    genomic <- if (strand == 1L) 99L + positions else 120L - positions
    if (strand == -1L) { refs <- unname(complement[refs]); alts <- unname(complement[alts]) }
    rows <- vapply(seq_along(positions), function(i) sprintf(paste0("(%d::UBIGINT,0::UINTEGER,",
      "%d::UBIGINT,%s::VARCHAR,%s::VARCHAR,1::UINTEGER,0::UINTEGER,0::UINTEGER,12::BIGINT,",
      "[1,1]::INTEGER[],[false,true]::BOOLEAN[])"), i, genomic[[i]],
      sql_quote(refs[[i]]), sql_quote(alts[[i]])), character(1L))
    create_calls <- paste0("CREATE OR REPLACE TABLE noncoding_calls AS SELECT * FROM (VALUES ",
      paste(rows[order(genomic)], collapse = ","), ") t(event_index,seq_region,position,reference,",
      "alternate,alt_index,transcript_index,sample_index,phase_set,alleles,phase_before)")
    DBI::dbExecute(connection, create_calls)
    relation <- paste0("duckvep_haplotypes(",
      sql_quote("SELECT * FROM noncoding_calls ORDER BY seq_region,position,event_index"), ",",
      sql_quote(model), ")")
    result <- DBI::dbGetQuery(connection, paste0("SELECT cds,protein,prediction_status,",
      "array_to_string(haplotype_consequences,',') consequences,nmd_prediction,edit_count,",
      "array_to_string(list_sort(list_transform(contributors,lambda c: c.event_index)),',') source_ids,",
      "array_to_string(list_distinct(list_transform(contributors,lambda c: c.projection_status)),',') projections ",
      "FROM ", relation))
    stopifnot(nrow(result) == 1L, result$prediction_status == "predicted",
      result$consequences == fixture$term,
      result$cds == "ATGAAATAA", result$protein == "MK*", result$edit_count == 0,
      result$nmd_prediction == "not_applicable", result$source_ids == "1,2",
      result$projections == "outside_cds")
    capacity <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ",
      sub("\\)$", ",max_sequence_bases:=19)", relation)))
    stopifnot(nrow(capacity) == 1L, all(capacity$prediction_status == "unsupported_context"))
    ambiguous <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM duckvep_haplotypes(",
      sql_quote(paste0("SELECT event_index,seq_region,position,reference,alternate,",
        "alt_index,transcript_index,sample_index,phase_set,[0,1]::INTEGER[] alleles,",
        "[false,false]::BOOLEAN[] phase_before FROM noncoding_calls")), ",", sql_quote(model), ")"))
    stopifnot(nrow(ambiguous) > 0L, all(ambiguous$prediction_status != "predicted"))
    DBI::dbExecute(connection, "UPDATE noncoding_calls SET alternate=alternate || 'A' WHERE event_index=2")
    length_change <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ", relation))
    stopifnot(nrow(length_change) == 1L, length_change$prediction_status == "unsupported_context")
    DBI::dbExecute(connection, create_calls)
    unchanged <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,",
      "array_to_string(haplotype_consequences,',') consequences FROM duckvep_haplotypes(",
      sql_quote("SELECT * REPLACE ([0,0]::INTEGER[] AS alleles) FROM noncoding_calls ORDER BY position"),
      ",", sql_quote(model), ")"))
    stopifnot(nrow(unchanged) == 0L)
    other <- fixtures[[fixture$other]]
    other_position <- other$positions[[1L]]
    other_ref <- substr(reference, other_position, other_position)
    other_alt <- other$alternate
    if (strand == -1L) { other_ref <- complement[[other_ref]]; other_alt <- complement[[other_alt]] }
    DBI::dbExecute(connection, sprintf(paste0("UPDATE noncoding_calls SET position=%d,",
      "reference=%s,alternate=%s WHERE event_index=2"),
      if (strand == 1L) 99L + other_position else 120L - other_position,
      sql_quote(other_ref), sql_quote(other_alt)))
    mixed <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,",
      "array_to_string(haplotype_consequences,',') consequences FROM ", relation))
    stopifnot(nrow(mixed) == 1L, mixed$prediction_status == "unsupported_context", is.na(mixed$consequences))
    DBI::dbExecute(connection, create_calls)
    replayed <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,",
      "array_to_string(haplotype_consequences,',') consequences FROM ", relation))
    stopifnot(nrow(replayed) == 1L, replayed$prediction_status == "predicted",
      replayed$consequences == fixture$term)
    if (feature == "utr5") {
      spliced_model <- paste0(model, "_spliced")
      spliced_start <- axis$spliced_start
      spliced_transcript <- sub("119::UBIGINT transcript_end", "125::UBIGINT transcript_end",
        transcript, fixed = TRUE)
      spliced_transcript <- sub(sprintf("%d::UBIGINT cds_start,%d::UBIGINT cds_end", cds_start, cds_end),
        sprintf("%d::UBIGINT cds_start,%d::UBIGINT cds_end", spliced_start, spliced_start + 8L),
        spliced_transcript, fixed = TRUE)
      spliced_exons <- paste0("SELECT 0::UINTEGER transcript_index,s::UBIGINT exon_start,",
        "e::UBIGINT exon_end,c::UBIGINT exon_cdna_start,d::UBIGINT exon_cdna_end,",
        "0::TINYINT phase,0::TINYINT end_phase FROM (VALUES ", axis$exons, ") t(s,e,c,d)")
      spliced_loaded <- DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(",
        sql_quote(spliced_model), ",", sql_quote("SELECT 0::UINTEGER seq_region"), ",",
        sql_quote(spliced_transcript), ",", sql_quote(spliced_exons), ")"))
      stopifnot(spliced_loaded$loaded == 1L)
      if (strand == -1L) DBI::dbExecute(connection,
        "UPDATE noncoding_calls SET position=position+6")
      spliced <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ",
        sub(sql_quote(model), sql_quote(spliced_model), relation, fixed = TRUE)))
      stopifnot(nrow(spliced) == 1L, spliced$prediction_status == "unsupported_context")
      DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", sql_quote(spliced_model), ")"))
      DBI::dbExecute(connection, create_calls)
    }
    missing_model <- paste0(model, "_no_flanks")
    missing_transcript <- sub(",1::UTINYINT codon_table,.*$", ",1::UTINYINT codon_table", transcript)
    missing_loaded <- DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(",
      sql_quote(missing_model), ",", sql_quote("SELECT 0::UINTEGER seq_region"), ",",
      sql_quote(missing_transcript), ",", sql_quote(exons), ")"))
    stopifnot(missing_loaded$loaded == 1L)
    missing <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ",
      sub(sql_quote(model), sql_quote(missing_model), relation, fixed = TRUE)))
    stopifnot(nrow(missing) == 1L, missing$prediction_status == "unsupported_context")
    for (name in c(model, missing_model))
      DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", sql_quote(name), ")"))
    checked <- checked + 1L
  }
  run_full_noncoding_producer(connection, sql_quote)
  cat("Native final-allele UTR:", checked,
    "strand/feature cases; raw operands, contributors, reference calls, capacity, phase, mixed-feature, length-change, splice and missing-metadata guards verified\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: noncoding_haplotype_runtime.R [repository]")
main(normalizePath(if (length(args) == 1L) args[[1L]] else ".", mustWork = TRUE))
