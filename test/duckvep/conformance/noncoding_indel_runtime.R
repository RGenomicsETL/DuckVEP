run_noncoding_indels <- function(connection, quote) {
  execute <- function(sql) DBI::dbExecute(connection, sql)
  rc <- reverse_complement
  reference <- "AAAAAAATGAAATAACCCCC"
  for (lower in c(FALSE, TRUE)) for (strand in c(1L, -1L)) for (feature in c("nc", "utr5", "utr3")) {
    model <- paste0("rna_indel_", feature, "_", if (strand == 1L) "plus" else "minus", if (lower) "_lower")
    coding <- feature != "nc"
    start <- if (strand == 1L) 106L else 105L
    stored <- if (lower) tolower(reference) else reference
    transcripts <- sprintf(paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
      "100::UBIGINT transcript_start,119::UBIGINT transcript_end,%d::TINYINT strand,",
      "0::UINTEGER gene_index,3::UBIGINT transcript_flags,%s::UBIGINT cds_start,%s::UBIGINT cds_end,",
      "%s::BLOB cds_sequence,%s::UTINYINT codon_table,%s::BLOB pre_cds_sequence,",
      "%s::BLOB post_cds_sequence,%s::BLOB cdna_sequence"), strand,
      if (coding) start else "NULL", if (coding) start + 8L else "NULL",
      if (coding) quote(substr(stored, 7L, 15L)) else "NULL", if (coding) "1" else "NULL",
      if (coding) quote(substr(stored, 1L, 6L)) else "NULL",
      if (coding) quote(substr(stored, 16L, 20L)) else "NULL", quote(stored))
    exons <- paste0("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,119::UBIGINT exon_end,",
      "1::UBIGINT exon_cdna_start,20::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase")
    loaded <- DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(", quote(model), ",",
      quote("SELECT 0::UINTEGER seq_region"), ",", quote(transcripts), ",", quote(exons), ")"))
    stopifnot(loaded$loaded == 1L)
    fingerprint_sql <- paste0("SELECT _duckvep_model_fingerprint(", quote(model), ") AS fingerprint")
    fingerprint <- DBI::dbGetQuery(connection, fingerprint_sql)
    cases <- rna_edit_cases(if (feature == "utr3") "C" else "A")
    cases$endpoint <- list(list(position = if (feature == "utr3") 20L else 1L,
      reference = if (feature == "utr3") "C" else "A", alternate = "GT"))
    for (name in names(cases)) {
      edits <- cases[[name]]
      if (feature == "utr3" && name != "endpoint")
        edits <- lapply(edits, function(e) {
          e$position <- if (e$position == 5L) 19L else e$position + 15L
          e
        })
      expected <- compose_rna(reference, edits)
      if (coding) {
        shift <- if (feature == "utr5") expected$delta else 0L
        stopifnot(substr(expected$sequence, 7L + shift, 15L + shift) == "ATGAAATAA")
      }
      rows <- vapply(seq_along(edits), function(i) {
        e <- edits[[i]]
        position <- if (strand == 1L) 99L + e$position else 121L - e$position - nchar(e$reference)
        ref <- if (strand == 1L) e$reference else rc(e$reference)
        alt <- if (strand == 1L) e$alternate else rc(e$alternate)
        sprintf("(%d,0,%d,%s,%s,1,0,0,12,[1,1],[false,true])", i, position, quote(ref), quote(alt))
      }, character(1L))
      create <- paste0("CREATE OR REPLACE TABLE rna_calls AS SELECT * FROM (VALUES ", paste(rows, collapse = ","),
        ") t(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,sample_index,phase_set,alleles,phase_before)")
      execute(create)
      relation <- paste0("duckvep_haplotypes(", quote("SELECT * FROM rna_calls ORDER BY position,event_index"), ",", quote(model), ")")
      observed <- DBI::dbGetQuery(connection, paste0("SELECT cds,protein,prediction_status,nominal_length_diff,edit_count,nmd_prediction,carrier_count,",
        "list_bool_and(list_transform(carriers,lambda c: c.sample_index=0 AND c.phase_set IS NULL AND c.ploidy=2)) carrier_keys,",
        "array_to_string(list_sort(list_transform(carriers,lambda c: c.haplotype_lane)),',') lanes,",
        "array_to_string(haplotype_consequences,',') consequences,",
        "array_to_string(list_sort(list_transform(contributors,lambda c: c.event_index)),',') ids,",
        "array_to_string(list_sort(list_transform(contributors,lambda c: c.reference || '>' || c.alternate)),',') operands FROM ", relation))
      term <- switch(feature, nc = "non_coding_transcript_exon_variant", utr5 = "5_prime_UTR_variant", utr3 = "3_prime_UTR_variant")
      stopifnot(nrow(observed) == 1L, observed$prediction_status == "predicted", observed$edit_count == 0,
        if (coding) observed$nominal_length_diff == 0 else is.na(observed$nominal_length_diff),
        observed$carrier_count == 2L, observed$carrier_keys, observed$lanes == "1,2",
        observed$ids == paste(seq_along(edits), collapse = ","))
      if (expected$sequence == reference) stopifnot(is.na(observed$consequences) || observed$consequences == "")
      else stopifnot(observed$consequences == term)
      expected_operands <- vapply(edits, function(e) if (strand == 1L) paste0(e$reference, ">", e$alternate)
        else paste0(rc(e$reference), ">", rc(e$alternate)), character(1L))
      stopifnot(observed$operands == paste(sort(expected_operands), collapse = ","))
      if (coding) stopifnot(observed$cds == "ATGAAATAA", observed$protein == "MK*")
      else stopifnot(is.na(observed$cds[[1L]]), is.na(observed$protein[[1L]]), observed$nmd_prediction == "not_applicable")
      short <- max(nchar(reference), nchar(expected$sequence)) - 1L
      refused <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,prediction_reason FROM ", sub("\\)$", paste0(",max_sequence_bases:=", short, ")"), relation)))
      stopifnot(refused$prediction_status == "unsupported_context", refused$prediction_reason == "sequence_capacity")
      execute("UPDATE rna_calls SET reference='N' || substr(reference,2) WHERE event_index=1")
      refused <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ", relation))
      stopifnot(refused$prediction_status != "predicted")
      execute(create)
      stopifnot(DBI::dbGetQuery(connection, paste0("SELECT prediction_status FROM ", relation))$prediction_status == "predicted",
        DBI::dbGetQuery(connection, "SELECT bool_and(phase_set=12 AND alleles=[1,1] AND phase_before=[false,true]) unchanged FROM rna_calls")$unchanged,
        identical(DBI::dbGetQuery(connection, fingerprint_sql), fingerprint))
    }
    if (!coding) {
      # The raw one-slot GT admits a conditional empty-ALT second slot.
      calls <- paste0("SELECT 1::UBIGINT event_index,0::UINTEGER seq_region,100::UBIGINT AS position,",
        quote(if (strand == 1L) reference else rc(reference)), "::VARCHAR reference,['A']::VARCHAR[] alternates,",
        "0::UINTEGER transcript_index,0::UINTEGER sample_index,'1'::VARCHAR gt")
      observed <- DBI::dbGetQuery(connection, paste0("SELECT cds,protein,prediction_status,prediction_reason,nmd_prediction,",
        "array_to_string(haplotype_consequences,',') consequence,contributors[1].alt_index IS NULL undefined,",
        "contributors[1].alternate='' deleted FROM duckvep_haplotypes(", quote(calls), ",", quote(model),
        ",input_mode:='source_records',phase_policy:='vep_compat') WHERE carriers[1].haplotype_lane=2"))
      stopifnot(nrow(observed) == 1L, observed$prediction_status == "unsupported_context",
        observed$prediction_reason == "non_strict_phase_policy", is.na(observed$consequence),
        observed$undefined, observed$deleted, is.na(observed$cds[[1L]]), is.na(observed$protein[[1L]]), observed$nmd_prediction == "unknown",
        identical(DBI::dbGetQuery(connection, fingerprint_sql), fingerprint))
    }
    DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", quote(model), ")"))
  }
  # A stop at cDNA 15 lies 65 bases before the original final junction (80).
  # Deleting 20 bases exclusively in the 3' UTR moves J to 60: EJC50 escapes.
  transcripts <- paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,100::UBIGINT transcript_start,",
    "249::UBIGINT transcript_end,1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
    "106::UBIGINT cds_start,114::UBIGINT cds_end,'ATGAAATAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "'AAAAAA'::BLOB pre_cds_sequence,repeat('C',115)::BLOB post_cds_sequence")
  exons <- paste0("SELECT 0::UINTEGER transcript_index,s::UBIGINT exon_start,e::UBIGINT exon_end,",
    "a::UBIGINT exon_cdna_start,z::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase ",
    "FROM (VALUES (100,179,1,80),(200,249,81,130)) t(s,e,a,z)")
  DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load('utr_nmd',", quote("SELECT 0::UINTEGER seq_region"),
    ",", quote(transcripts), ",", quote(exons), ")"))
  for (deleted in c(FALSE, TRUE)) {
    calls <- paste0("SELECT 1::UBIGINT event_index,0::UINTEGER seq_region,120::UBIGINT AS position,",
      quote(if (deleted) strrep("C",21L) else "C"), "::VARCHAR reference,",
      quote(if (deleted) "C" else "CAAAAAAAAAA"), "::VARCHAR alternate,1::UINTEGER alt_index,",
      "0::UINTEGER transcript_index,0::UINTEGER sample_index,12::BIGINT phase_set,[1,1]::INTEGER[] alleles,[false,true]::BOOLEAN[] phase_before")
    observed <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,nmd_prediction,nmd_stop_position,nmd_junction_position,cds,protein FROM duckvep_haplotypes(", quote(calls), ",'utr_nmd')"))
    stopifnot(observed$prediction_status == "predicted", observed$cds == "ATGAAATAA", observed$protein == "MK*",
      observed$nmd_stop_position == 15, observed$nmd_junction_position == if (deleted) 60 else 90,
      observed$nmd_prediction == if (deleted) "escape" else "trigger")
  }
  unknown_transcripts <- sub("3::UBIGINT transcript_flags", "259::UBIGINT transcript_flags", transcripts, fixed = TRUE)
  DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load('utr_nmd_unknown',", quote("SELECT 0::UINTEGER seq_region"),
    ",", quote(unknown_transcripts), ",", quote(exons), ")"))
  unknown <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,prediction_reason,nmd_prediction FROM duckvep_haplotypes(", quote(calls), ",'utr_nmd_unknown')"))
  stopifnot(unknown$prediction_status == "predicted", unknown$nmd_prediction == "unknown", unknown$prediction_reason == "nmd_geometry_unavailable")
  positional <- DBI::dbGetQuery(connection, paste0("SELECT cds,protein,array_to_string(haplotype_consequences,',') terms,",
    "prediction_reference_protein,prediction_protein FROM duckvep_haplotypes(", quote(calls), ",'utr_nmd_unknown')"))
  stopifnot(positional$cds == "ATGAAATAA", positional$protein == "MK*",
    positional$terms == "3_prime_UTR_variant", is.na(positional$prediction_reference_protein),
    is.na(positional$prediction_protein))
  raw_cdna <- paste0("'AAAAAAATGAAATAA", strrep("C", 115L), "'::BLOB cdna_sequence")
  withheld <- paste0("SELECT * REPLACE (NULL::UBIGINT AS cds_start,NULL::UBIGINT AS cds_end,",
    "NULL::BLOB AS cds_sequence,NULL::UTINYINT AS codon_table,",
    "NULL::BLOB AS pre_cds_sequence,NULL::BLOB AS post_cds_sequence),", raw_cdna,
    " FROM (", unknown_transcripts, ")")
  supplied <- paste0("SELECT *,", raw_cdna, ",7::UBIGINT rna_cds_cdna_start,",
    "15::UBIGINT rna_cds_cdna_end,1::UTINYINT rna_codon_table FROM (", unknown_transcripts, ")")
  programme <- paste0("SELECT 0::UINTEGER transcript_index,10::UINTEGER cdna_start,10::UINTEGER cdna_end,",
    "'A'::VARCHAR expected_substrate,'G'::VARCHAR alternate_sequence,'substitution'::VARCHAR operation_kind,",
    "'Ensembl116:_rna_edit'::VARCHAR programme_name,'transcript_attrib:SeqEdit'::VARCHAR evidence_name,",
    "0::UTINYINT boundary_overrides")
  for (kind in c("withheld", "supplied")) {
    rna_calls <- paste0("SELECT * REPLACE (",
      if (kind == "withheld") "109::UBIGINT AS position,'A'::VARCHAR AS reference,'G'::VARCHAR AS alternate" else
        "'C'::VARCHAR AS reference,'T'::VARCHAR AS alternate", ") FROM (", calls, ")")
    model <- paste0("utr_rna_", kind)
    transcript_query <- if (kind == "withheld") withheld else supplied
    DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(", quote(model), ",",
      quote("SELECT 0::UINTEGER seq_region"), ",", quote(transcript_query), ",", quote(exons),
      if (kind == "supplied") paste0(",rna_edit_query := ", quote(programme)), ")"))
    refused <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,prediction_reason,nmd_prediction,",
      "haplotype_consequences IS NULL AS no_terms,prediction_reference_protein,prediction_protein,cds,protein ",
      "FROM duckvep_haplotypes(", quote(rna_calls), ",", quote(model), ")"))
    stopifnot(refused$prediction_status == "unsupported_context",
      refused$prediction_reason == "unsupported_curated_edit", refused$nmd_prediction == "unknown",
      refused$no_terms, is.na(refused$prediction_reference_protein), is.na(refused$prediction_protein))
    if (kind == "withheld") stopifnot(is.na(refused$cds), is.na(refused$protein)) else
      stopifnot(refused$cds == "ATGAAATAA", refused$protein == "MK*")
    DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", quote(model), ")"))
  }
  DBI::dbGetQuery(connection, "SELECT duckvep_model_drop('utr_nmd_unknown')")
  DBI::dbGetQuery(connection, "SELECT duckvep_model_drop('utr_nmd')")
  for (lower in c(FALSE, TRUE)) for (strand in c(1L, -1L)) {
    model <- paste0("utr_nmd_restore_", strand == -1L, "_", lower)
    starts <- if (strand == 1L) c(100L, 200L) else c(200L, 100L)
    ends <- starts + c(64L, 29L)
    rna <- paste0("AAAAAAATGAAATAA", strrep("C", 80L))
    stored <- if (lower) tolower(rna) else rna
    transcripts <- sprintf(paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,100::UBIGINT transcript_start,",
      "%d::UBIGINT transcript_end,%d::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
      "%d::UBIGINT cds_start,%d::UBIGINT cds_end,%s::BLOB cds_sequence,1::UTINYINT codon_table,",
      "%s::BLOB pre_cds_sequence,%s::BLOB post_cds_sequence"), max(ends), strand,
      if (strand == 1L) 106L else 250L, if (strand == 1L) 114L else 258L,
      quote(substr(stored, 7L, 15L)), quote(substr(stored, 1L, 6L)), quote(substr(stored, 16L, 95L)))
    exons <- sprintf(paste0("SELECT 0::UINTEGER transcript_index,s::UBIGINT exon_start,e::UBIGINT exon_end,",
      "a::UBIGINT exon_cdna_start,z::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase ",
      "FROM (VALUES (%d,%d,1,65),(%d,%d,66,95)) t(s,e,a,z)"), starts[1L], ends[1L], starts[2L], ends[2L])
    stopifnot(DBI::dbGetQuery(connection, paste0("SELECT loaded FROM duckvep_model_load(", quote(model), ",",
      quote("SELECT 0::UINTEGER seq_region"), ",", quote(transcripts), ",", quote(exons), ")"))$loaded == 1L)
    for (kind in c("same_exon", "same_exon_terminal", "transfer_forward", "transfer_reverse")) {
      restored_geometry <- kind %in% c("same_exon", "same_exon_terminal")
      positions <- if (kind == "same_exon_terminal") c(75L, 95L) else
        c(25L, if (kind == "same_exon") 35L else 75L)
      refs <- if (kind %in% c("transfer_reverse", "same_exon_terminal")) c("CC", "C") else c("C", "CC")
      alts <- rev(refs)
      edits <- lapply(seq_along(positions), function(i) list(position = positions[i], reference = refs[i], alternate = alts[i]))
      stopifnot(compose_rna(rna, edits)$sequence == rna)
      # Exon lengths are 65 and 30; the stop occupies cDNA 13..15.
      junction <- 65L + sum((nchar(alts) - nchar(refs))[positions <= 65L])
      rows <- vapply(seq_along(positions), function(i) {
        exon <- if (positions[i] <= 65L) 1L else 2L
        offset <- positions[i] - if (exon == 1L) 1L else 66L
        position <- if (strand == 1L) starts[exon] + offset else ends[exon] - offset - nchar(refs[i]) + 1L
        sprintf("(%d,0,%d,%s,%s,1,0,0,12,[1,1],[false,true])", i, position,
          quote(if (strand == 1L) refs[i] else rc(refs[i])), quote(if (strand == 1L) alts[i] else rc(alts[i])))
      }, character(1L))
      calls <- paste0("SELECT * FROM (VALUES ", paste(rows, collapse = ","),
        ") t(event_index,seq_region,position,reference,alternate,alt_index,transcript_index,sample_index,phase_set,alleles,phase_before) ORDER BY position,event_index")
      observed <- DBI::dbGetQuery(connection, paste0("SELECT prediction_status,cds,protein,nmd_prediction,nmd_stop_position,nmd_junction_position,",
        "len(haplotype_consequences) terms,len(contributors) sources,",
        "array_to_string(list_sort(list_transform(contributors,lambda c: c.event_index)),',') ids,",
        "list_bool_and(list_transform(contributors,lambda c: c.projection_status='outside_cds')) outside_cds ",
        "FROM duckvep_haplotypes(", quote(calls), ",", quote(model), ")"))
      stopifnot(observed$prediction_status == "predicted", observed$terms == 0L,
        observed$cds == "ATGAAATAA", observed$protein == "MK*", observed$sources == 2L,
        observed$ids == "1,2", observed$outside_cds)
      if (restored_geometry) stopifnot(observed$nmd_prediction == "not_applicable",
        is.na(observed$nmd_stop_position), is.na(observed$nmd_junction_position)) else
        stopifnot(observed$nmd_stop_position == 15L, observed$nmd_junction_position == junction,
          observed$nmd_prediction == if (junction - 15L > 50L) "trigger" else "escape")
    }
    DBI::dbGetQuery(connection, paste0("SELECT duckvep_model_drop(", quote(model), ")"))
  }
  for (strand in c(1L, -1L)) for (last_length in 406:408) {
    starts <- if (strand == 1L) c(100L, 200L) else c(700L, 100L)
    ends <- starts + c(80L, last_length) - 1L
    cds <- if (strand == 1L) c(206L, 214L) else c(ends[2L] - 14L, ends[2L] - 6L)
    tx <- paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,100::UBIGINT transcript_start,",
      max(ends), "::UBIGINT transcript_end,", strand, "::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
      cds[1L], "::UBIGINT cds_start,", cds[2L], "::UBIGINT cds_end,",
      "'ATGAAATAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
      quote(strrep("A", 86L)), "::BLOB pre_cds_sequence,", quote(strrep("C", last_length - 15L)), "::BLOB post_cds_sequence")
    exons <- sprintf("SELECT 0::UINTEGER transcript_index,s::UBIGINT exon_start,e::UBIGINT exon_end,a::UBIGINT exon_cdna_start,z::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase FROM (VALUES (%d,%d,1,80),(%d,%d,81,%d)) v(s,e,a,z)",
      starts[1L], ends[1L], starts[2L], ends[2L], 80L + last_length)
    model <- paste("terminal_nmd", strand, last_length, sep = "_")
    load <- sprintf("SELECT loaded FROM duckvep_model_load(%s,'SELECT 0::UINTEGER seq_region',%s,%s)", quote(model), quote(tx), quote(exons))
    stopifnot(isTRUE(DBI::dbGetQuery(connection, load)$loaded))
    position <- if (strand == 1L) ends[2L] else starts[2L]
    ref <- if (strand == 1L) "C" else "G"
    calls <- sprintf("SELECT 1::UBIGINT event_index,0::UINTEGER seq_region,%d::BIGINT AS position,%s AS reference,'CG' AS alternate,1 AS alt_index,0 AS transcript_index,0 AS sample_index,[0,1] AS alleles,[false,true] AS phase_before,7::BIGINT AS phase_set", position, quote(ref))
    result <- DBI::dbGetQuery(connection, sprintf("SELECT * FROM duckvep_haplotypes(%s,%s)", quote(calls), quote(model)))
    stopifnot(nrow(result) == 1L, result$prediction_status == "predicted", result$nmd_prediction == "escape",
      result$nmd_stop_position == 95, result$nmd_junction_position == 80,
      grepl("long_exon", result$nmd_exceptions, fixed = TRUE) == (last_length + 1L > 407L),
      result$cds == "ATGAAATAA", result$protein == "MK*")
    DBI::dbGetQuery(connection, sprintf("SELECT duckvep_model_drop(%s)", quote(model)))
  }
  cat("RNA indel SQL matrix: case/strand/feature replay, ablation boundary, restored RNA and exon-aware EJC50/long-exon exceptions verified\n")
}
