#!/usr/bin/env Rscript

main <- function() {
  extension <- Sys.getenv("DUCKVEP_EXTENSION", "build/release/duckvep.duckdb_extension")
  extension <- normalizePath(extension, mustWork = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))

  cases <- list(
    list(
      name = "sec", cds = "ATGTGAGCCTAA", codon_table = 1L,
      edit_code = "_selenocysteine", edit_position = 2L, edit_amino_acid = "U",
      reference_protein = "MUA*",
      variants = list(
        list(position = 8L, reference = "C", alternate = "T", prediction = "MUV*", protein = "M*"),
        list(position = 5L, reference = "G", alternate = "C", prediction = "MSA*", protein = "MSA*")
      )
    ),
    list(
      name = "readthrough", cds = "ATGTAAGCCTAA", codon_table = 1L,
      edit_code = "_stop_codon_rt", edit_position = 2L, edit_amino_acid = "W",
      reference_protein = "MWA*",
      variants = list(
        list(position = 8L, reference = "C", alternate = "T", prediction = "MWV*", protein = "M*"),
        list(position = 5L, reference = "A", alternate = "C", prediction = "MSA*", protein = "MSA*")
      )
    ),
    list(
      name = "initial_met", cds = "GTG GCC GCT TAA", codon_table = 1L,
      edit_code = "initial_met", edit_position = 1L, edit_amino_acid = "M",
      reference_protein = "MAA*",
      variants = list(
        list(position = 5L, reference = "C", alternate = "T", prediction = "MVA*", protein = "VVA*"),
        list(position = 2L, reference = "T", alternate = "C", prediction = "AAA*", protein = "AAA*")
      )
    ),
    list(
      name = "mitochondrial_initial_met", cds = "ATA GCC GCT TAA", codon_table = 2L,
      edit_code = "amino_acid_sub", edit_position = 1L, edit_amino_acid = "M",
      reference_protein = "MAA*",
      variants = list(
        list(position = 5L, reference = "C", alternate = "T", prediction = "MVA*", protein = "MVA*"),
        list(position = 1L, reference = "A", alternate = "G", prediction = "VAA*", protein = "VAA*")
      )
    ),
    list(
      name = "early_alternate_stop", cds = "ATGCAATGATAA", codon_table = 1L,
      edit_code = "_selenocysteine", edit_position = 3L, edit_amino_acid = "U",
      reference_protein = "MQU*",
      variants = list(
        list(position = 4L, reference = "C", alternate = "T", prediction = "M*", protein = "M*")
      )
    )
  )

  model_queries <- function(cds, codon_table, transcript_flags = 0L) {
    cds <- gsub(" ", "", cds)
    cds_length <- nchar(cds)
    c(
      paste0(
        "SELECT 0::UINTEGER seq_region, ", cds_length,
        "::UBIGINT sequence_length, 'chr1'::VARCHAR seq_region_name"
      ),
      paste0(
        "SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region, ",
        "1::UBIGINT transcript_start, ", cds_length,
        "::UBIGINT transcript_end, 1::TINYINT strand, 0::UINTEGER gene_index, ",
        transcript_flags, "::UBIGINT transcript_flags, 1::UBIGINT cds_start, ",
        cds_length, "::UBIGINT cds_end, ", DBI::dbQuoteString(con, gsub(" ", "", cds)),
        "::BLOB cds_sequence, ", codon_table, "::UTINYINT codon_table"
      ),
      paste0(
        "SELECT 0::UINTEGER transcript_index, 1::UBIGINT exon_start, ",
        cds_length, "::UBIGINT exon_end, 1::UBIGINT exon_cdna_start, ",
        cds_length, "::UBIGINT exon_cdna_end, 0::TINYINT phase, 0::TINYINT end_phase"
      )
    )
  }

  calls_query <- function(variant) {
    paste0(
      "SELECT 1::UBIGINT record_index, 1::UBIGINT event_index, 0::UINTEGER seq_region, ",
      variant$position, "::UBIGINT \"position\", ",
      DBI::dbQuoteString(con, variant$reference), "::VARCHAR reference, ",
      DBI::dbQuoteString(con, variant$alternate), "::VARCHAR alternate, ",
      "1::UINTEGER alt_index, 0::UINTEGER transcript_index, 0::UINTEGER sample_index, ",
      "[0,1]::INTEGER[] alleles, [TRUE]::BOOLEAN[] phase_before, 1::BIGINT phase_set"
    )
  }

  haplotype_query <- function(calls, model) {
    paste0(
      "SELECT prediction_status,prediction_reason,prediction_reference_protein,",
      "prediction_protein,protein,nmd_prediction,",
      "haplotype_consequences IS NULL AS no_unconditional_consequence ",
      "FROM duckvep_haplotypes(", DBI::dbQuoteString(con, calls), ",",
      DBI::dbQuoteString(con, model), ",phase_policy:='strict')"
    )
  }

  sec_model <- "curated_sec"
  sec_model_queries <- NULL
  for (case in cases) {
    model <- paste0("curated_", case$name)
    queries <- model_queries(case$cds, case$codon_table)
    if (case$name == "sec") sec_model_queries <- queries
    typed_peptide_query <- paste0(
      "SELECT 0::UINTEGER transcript_index, ", case$edit_position,
      "::UINTEGER protein_position, ", DBI::dbQuoteString(con, case$edit_amino_acid),
      "::VARCHAR alternate_amino_acid, ", DBI::dbQuoteString(con, case$edit_code),
      "::VARCHAR edit_code"
    )
    loaded <- DBI::dbGetQuery(con, paste0(
      "SELECT loaded FROM duckvep_model_load(", DBI::dbQuoteString(con, model), ",",
      paste(DBI::dbQuoteString(con, queries), collapse = ","),
      ",peptide_edit_query:=", DBI::dbQuoteString(con, typed_peptide_query), ")"
    ))
    stopifnot(identical(loaded$loaded, TRUE))

    for (variant in case$variants) {
      observed <- DBI::dbGetQuery(con, haplotype_query(calls_query(variant), model))
      stopifnot(
        nrow(observed) == 1L,
        identical(observed$prediction_status, "conditional_assumed_recoding"),
        identical(observed$prediction_reason, "assumed_recoding_programme"),
        identical(observed$prediction_reference_protein, case$reference_protein),
        identical(observed$prediction_protein, variant$prediction),
        identical(observed$protein, variant$protein),
        identical(observed$nmd_prediction, "unknown"),
        isTRUE(observed$no_unconditional_consequence)
      )
    }
  }

  sec_calls <- calls_query(cases[[1L]]$variants[[1L]])
  snapshot <- tempfile(fileext = ".duckvep")
  on.exit(unlink(snapshot), add = TRUE)
  fingerprint <- DBI::dbGetQuery(con,
    paste0("SELECT _duckvep_model_fingerprint(", DBI::dbQuoteString(con, sec_model), ")"))[[1L]]
  stopifnot(identical(DBI::dbGetQuery(con, paste0(
    "SELECT duckvep_model_save(", DBI::dbQuoteString(con, sec_model), ",",
    DBI::dbQuoteString(con, snapshot), ")"
  ))[[1L]], TRUE))
  stopifnot(identical(DBI::dbGetQuery(con, paste0(
    "SELECT duckvep_model_restore('snapshot_curated_sec',",
    DBI::dbQuoteString(con, snapshot), ")"
  ))[[1L]], TRUE))
  stopifnot(identical(fingerprint, DBI::dbGetQuery(con,
    "SELECT _duckvep_model_fingerprint('snapshot_curated_sec')")[[1L]]))
  restored <- DBI::dbGetQuery(con, haplotype_query(sec_calls, "snapshot_curated_sec"))
  stopifnot(
    nrow(restored) == 1L,
    identical(restored$prediction_status, "conditional_assumed_recoding"),
    identical(restored$prediction_reference_protein, "MUA*"),
    identical(restored$prediction_protein, "MUV*"),
    identical(restored$protein, "M*"),
    identical(restored$nmd_prediction, "unknown"),
    isTRUE(restored$no_unconditional_consequence)
  )

  legacy_query <- paste(
    "SELECT 0::UINTEGER transcript_index, 2::UINTEGER protein_position,",
    "'U'::VARCHAR alternate_amino_acid"
  )
  legacy_model <- "legacy_curated_fixture"
  legacy_loaded <- DBI::dbGetQuery(con, paste0(
    "SELECT loaded FROM duckvep_model_load(", DBI::dbQuoteString(con, legacy_model), ",",
    paste(DBI::dbQuoteString(con, sec_model_queries), collapse = ","),
    ",peptide_edit_query:=", DBI::dbQuoteString(con, legacy_query), ")"
  ))
  stopifnot(identical(legacy_loaded$loaded, TRUE))
  legacy <- DBI::dbGetQuery(con, haplotype_query(sec_calls, legacy_model))
  stopifnot(
    nrow(legacy) == 1L,
    identical(legacy$prediction_status, "unsupported_context"),
    identical(legacy$prediction_reason, "untyped_curated_metadata"),
    is.na(legacy$prediction_reference_protein),
    is.na(legacy$prediction_protein),
    identical(legacy$protein, "M*")
  )

  for (typed in c(FALSE, TRUE)) {
    rna_model <- paste0("rna_edit_fixture_", typed)
    rna_edit_query <- paste0("SELECT *, '_selenocysteine'::VARCHAR edit_code FROM (", legacy_query, ")")
    rna_loaded <- DBI::dbGetQuery(con, paste0(
      "SELECT loaded FROM duckvep_model_load(", DBI::dbQuoteString(con, rna_model), ",",
      paste(DBI::dbQuoteString(con, model_queries("ATGTGAGCCTAA", 1L, 323L)), collapse = ","),
      if (typed) paste0(",peptide_edit_query:=", DBI::dbQuoteString(con, rna_edit_query)) else "", ")"
    ))
    stopifnot(identical(rna_loaded$loaded, TRUE))
    rna_edit <- DBI::dbGetQuery(con, haplotype_query(sec_calls, rna_model))
    stopifnot(
      nrow(rna_edit) == 1L,
      identical(rna_edit$prediction_status, "unsupported_context"),
      identical(rna_edit$prediction_reason,
        if (typed) "unsupported_curated_edit" else "untyped_curated_metadata"),
      is.na(rna_edit$prediction_reference_protein),
      is.na(rna_edit$prediction_protein)
    )
  }

  invalid_calls <- calls_query(list(position = 8L, reference = "A", alternate = "T"))
  invalid <- DBI::dbGetQuery(con, haplotype_query(invalid_calls, sec_model))
  stopifnot(
    nrow(invalid) == 1L,
    identical(invalid$prediction_status, "unsupported_context"),
    identical(invalid$prediction_reason, "reference_mismatch"),
    is.na(invalid$prediction_reference_protein),
    is.na(invalid$prediction_protein)
  )
  replay <- DBI::dbGetQuery(con, haplotype_query(sec_calls, sec_model))
  stopifnot(
    nrow(replay) == 1L,
    identical(replay$prediction_status, "conditional_assumed_recoding"),
    identical(replay$prediction_reference_protein, "MUA*"),
    identical(replay$prediction_protein, "MUV*")
  )

  producer_setup <- c(
    "CREATE SCHEMA curated_core",
    "CREATE TABLE curated_reference_chunks(chrom VARCHAR, \"start\" BIGINT, \"end\" BIGINT, seq VARCHAR)",
    "INSERT INTO curated_reference_chunks VALUES ('chr1', 0, 12, 'ATGTGAGCCTAA')",
    "CREATE TABLE curated_core.coord_system(coord_system_id BIGINT, species_id BIGINT, name VARCHAR, version VARCHAR, rank BIGINT)",
    "INSERT INTO curated_core.coord_system VALUES (1, 1, 'chromosome', 'GRCh38', 1)",
    "CREATE TABLE curated_core.seq_region(seq_region_id BIGINT, name VARCHAR, coord_system_id BIGINT, length BIGINT)",
    "INSERT INTO curated_core.seq_region VALUES (10, 'chr1', 1, 12)",
    "CREATE TABLE curated_core.seq_region_attrib(seq_region_id BIGINT, attrib_type_id BIGINT, value BIGINT)",
    "CREATE TABLE curated_core.gene(gene_id BIGINT, biotype VARCHAR, is_current BIGINT, stable_id VARCHAR, version BIGINT)",
    "INSERT INTO curated_core.gene VALUES (20, 'protein_coding', 1, 'ENSG_CURATED', 1)",
    paste(
      "CREATE TABLE curated_core.transcript(transcript_id BIGINT, gene_id BIGINT, seq_region_id BIGINT,",
      "seq_region_start BIGINT, seq_region_end BIGINT, seq_region_strand BIGINT, biotype VARCHAR,",
      "is_current BIGINT, stable_id VARCHAR, version BIGINT)"
    ),
    "INSERT INTO curated_core.transcript VALUES (30, 20, 10, 1, 12, 1, 'protein_coding', 1, 'ENST_CURATED', 1)",
    paste(
      "CREATE TABLE curated_core.exon(exon_id BIGINT, seq_region_id BIGINT, seq_region_start BIGINT,",
      "seq_region_end BIGINT, seq_region_strand BIGINT, phase BIGINT, end_phase BIGINT,",
      "is_current BIGINT, stable_id VARCHAR)"
    ),
    "INSERT INTO curated_core.exon VALUES (40, 10, 1, 12, 1, -1, 0, 1, 'ENSE_CURATED')",
    "CREATE TABLE curated_core.exon_transcript(exon_id BIGINT, transcript_id BIGINT, rank BIGINT)",
    "INSERT INTO curated_core.exon_transcript VALUES (40, 30, 1)",
    paste(
      "CREATE TABLE curated_core.translation(translation_id BIGINT, transcript_id BIGINT, seq_start BIGINT,",
      "start_exon_id BIGINT, seq_end BIGINT, end_exon_id BIGINT, stable_id VARCHAR, version BIGINT)"
    ),
    "INSERT INTO curated_core.translation VALUES (50, 30, 1, 40, 12, 40, 'ENSP_CURATED', 1)",
    "CREATE TABLE curated_core.attrib_type(attrib_type_id BIGINT, code VARCHAR)",
    "INSERT INTO curated_core.attrib_type VALUES (60, '_selenocysteine')",
    "CREATE TABLE curated_core.transcript_attrib(transcript_id BIGINT, attrib_type_id BIGINT, value VARCHAR)",
    "CREATE TABLE curated_core.translation_attrib(translation_id BIGINT, attrib_type_id BIGINT, value VARCHAR)",
    "INSERT INTO curated_core.translation_attrib VALUES (50, 60, '2 2 U')"
  )
  for (statement in producer_setup) DBI::dbExecute(con, statement)

  regions_sql <- DBI::dbGetQuery(con,
    "SELECT duckvep_ensembl_regions_sql('curated_core','curated_reference_chunks','GRCh38')")[[1L]]
  transcripts_sql <- DBI::dbGetQuery(con,
    "SELECT duckvep_ensembl_transcripts_sql('curated_core','curated_reference_chunks','GRCh38')")[[1L]]
  DBI::dbExecute(con, paste0(
    "CREATE TABLE producer_curated_transcripts AS SELECT * FROM query(",
    DBI::dbQuoteString(con, transcripts_sql), ")"
  ))
  producer_edit <- DBI::dbGetQuery(con, paste(
    "SELECT edit.protein_position, edit.alternate_amino_acid, edit.edit_code",
    "FROM producer_curated_transcripts, LATERAL unnest(peptide_edits) AS u(edit)"
  ))
  stopifnot(
    nrow(producer_edit) == 1L,
    producer_edit$protein_position == 2,
    identical(producer_edit$alternate_amino_acid, "U"),
    identical(producer_edit$edit_code, "_selenocysteine")
  )

  producer_edit_query <- paste(
    "SELECT transcript_index, edit.protein_position, edit.alternate_amino_acid, edit.edit_code",
    "FROM producer_curated_transcripts, LATERAL unnest(peptide_edits) AS u(edit)",
    "ORDER BY transcript_index, edit.protein_position"
  )
  producer_queries <- c(
    paste0(
      "SELECT seq_region, sequence_length FROM query(",
      DBI::dbQuoteString(con, regions_sql), ") ORDER BY seq_region"
    ),
    paste(
      "SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index,",
      "transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence",
      "FROM producer_curated_transcripts ORDER BY seq_region, transcript_start, transcript_index"
    ),
    paste(
      "SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start,",
      "exon.exon_cdna_end, exon.phase, exon.end_phase FROM producer_curated_transcripts,",
      "LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start"
    )
  )
  producer_model <- "producer_curated_fixture"
  producer_loaded <- DBI::dbGetQuery(con, paste0(
    "SELECT loaded FROM duckvep_model_load(", DBI::dbQuoteString(con, producer_model), ",",
    paste(DBI::dbQuoteString(con, producer_queries), collapse = ","),
    ",peptide_edit_query:=", DBI::dbQuoteString(con, producer_edit_query),
    ",transcript_coverage_complete:=true)"
  ))
  stopifnot(identical(producer_loaded$loaded, TRUE))
  producer_prediction <- DBI::dbGetQuery(con,
    haplotype_query(sec_calls, producer_model))
  stopifnot(
    nrow(producer_prediction) == 1L,
    identical(producer_prediction$prediction_status, "conditional_assumed_recoding"),
    identical(producer_prediction$prediction_reference_protein, "MUA*"),
    identical(producer_prediction$prediction_protein, "MUV*"),
    identical(producer_prediction$protein, "M*")
  )
}

main()
