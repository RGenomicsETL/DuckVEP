library(tinytest)
library(DBI)

# Coding-v1 eligibility and provenance (slice 2 of #2), the same-codon classifier (slice 3) and the frame/stop-gain
# classifier (slice 4) and the start/stop classifier (slice 5). prediction_status says whether a path is inside the supported domain and whether a
# classifier decided it (`predicted`); since slice 5 no eligible path is left pending (`eligible_classifier_pending`).
# Transcript 0 is + strand with exons 100-108 (ATGGCTGCT) and 120-128 (GAAGGTTAA).
# Transcripts 1-6 are single-exon variants that each break one v1 CDS requirement.
local({
  con <- rduckvep_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste(
    "SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,100::UBIGINT transcript_start,",
    "CASE WHEN i=0 THEN 128 ELSE 117 END::UBIGINT transcript_end,1::TINYINT strand,",
    "i::UINTEGER gene_index,CASE i WHEN 4 THEN 259 WHEN 5 THEN 19 ELSE 3 END::UBIGINT transcript_flags,",
    "100::UBIGINT cds_start,CASE WHEN i=0 THEN 128 ELSE 117 END::UBIGINT cds_end,",
    "CASE i WHEN 1 THEN 'CTGGCTGCTGAAGGTTAA' WHEN 2 THEN 'ATGGCTTAAGAAGGTTAA'",
    "WHEN 3 THEN 'ATGGCTGCTGAAGGTGGT' ELSE 'ATGGCTGCTGAAGGTTAA' END::BLOB cds_sequence,",
    "CASE WHEN i=6 THEN 2 ELSE 1 END::UTINYINT codon_table,''::BLOB pre_cds_sequence,",
    "''::BLOB post_cds_sequence FROM range(8) t(i)")
  exons <- paste(
    "SELECT * FROM (SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,",
    "(CASE WHEN i=0 THEN 108 ELSE 117 END)::UBIGINT exon_end,1::UBIGINT exon_cdna_start,",
    "(CASE WHEN i=0 THEN 9 ELSE 18 END)::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM range(8) t(i) UNION ALL SELECT 0::UINTEGER,120::UBIGINT,128::UBIGINT,10::UBIGINT,",
    "18::UBIGINT,0::TINYINT,0::TINYINT) ORDER BY transcript_index,exon_start")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('elig',",
    "'SELECT i::UINTEGER seq_region FROM range(8) t(i)',", dbQuoteString(con, tx), ",",
    dbQuoteString(con, exons), ")"))$loaded)
  dbExecute(con, paste("CREATE TABLE elig_events(event_index BIGINT, seq_region INT, position BIGINT,",
    "reference VARCHAR, alternate VARCHAR, alt_index INT, transcript_index INT)"))
  dbExecute(con, paste("INSERT INTO elig_events VALUES",
    "(1,0,103,'G','A',1,0),(2,0,104,'C','T',1,0),(3,0,123,'G','C',1,0),(5,0,102,'GGCT','G',1,0),",
    "(6,0,103,'G','T',1,0),(10,0,103,'GCT','TAG',1,0),(11,0,112,'A','G',1,0),",
    "(16,0,103,'GC','AT',1,0),(20,0,106,'G','T',2,0),(50,0,108,'T','C',1,0)"))
  dbExecute(con, "INSERT INTO elig_events SELECT 100+i,i,103,'G','A',1,i FROM range(1,7) t(i)")
  dbExecute(con, paste("CREATE TABLE elig_gt(event_index BIGINT, sample_index INT, alleles INTEGER[],",
    "phase_before BOOLEAN[], phase_set BIGINT)"))
  dbExecute(con, paste("INSERT INTO elig_gt VALUES",
    "(1,1,[1,0],[false,true],NULL),(2,1,[1,0],[false,true],NULL),",
    "(1,2,[1,0],[false,true],NULL),(2,2,[0,1],[false,true],NULL),",
    "(1,3,[1,NULL],[false,true],NULL),(1,4,[0,1],[false,false],NULL),",
    "(1,5,[1,0],[false,true],10),(3,5,[1,0],[false,true],20),",
    "(1,7,[1,0],[false,true],NULL),(5,7,[1,0],[false,true],NULL),",
    "(1,8,[1,0],[false,true],NULL),(6,8,[1,0],[false,true],NULL),",
    "(11,12,[1,0],[false,true],NULL),(3,12,[1,0],[false,true],NULL),",
    "(10,16,[1,0],[false,true],NULL),(3,16,[1,0],[false,true],NULL),",
    "(1,17,[1],[false],NULL),(16,30,[1,0],[false,true],NULL),",
    "(20,40,[2,0],[false,true],NULL),(50,50,[1,0],[false,true],NULL)"))
  dbExecute(con, "INSERT INTO elig_gt SELECT 100+i,20+i,[1,0],[false,true],NULL FROM range(1,7) t(i)")
  dbExecute(con, paste("CREATE TABLE elig_calls AS SELECT e.event_index,e.seq_region,e.position,",
    "e.reference,e.alternate,e.alt_index,e.transcript_index,g.sample_index,g.alleles,",
    "g.phase_before,g.phase_set FROM elig_events e JOIN elig_gt g USING(event_index)"))

  result <- rduckvep_haplotypes(con, "SELECT * FROM elig_calls", "elig")
  # The public relation keeps raw and conditional prediction operands separate.
  expect_identical(ncol(result), 36L)
  expect_identical(names(result)[16:19], c("hgvsc", "hgvsc_status", "hgvsp", "hgvsp_status"))
  expect_identical(names(result)[20:31], c("prediction_policy", "prediction_status",
    "prediction_reason", "contributor_provenance", "normalized_edits", "carrier_predictions",
    "haplotype_consequences", "haplotype_impact", "nmd_rule", "nmd_prediction", "nmd_stop_position",
    "nmd_junction_position"))
  expect_identical(names(result)[32:36], c("nmd_contributors", "nmd_exceptions",
    "prediction_reference_protein", "prediction_protein", "nominal_length_diff"))
  expect_true(all(result$prediction_policy == "duckvep-coding"))

  keyed <- do.call(rbind, lapply(result$carrier_predictions, function(x) x))
  keyed <- keyed[order(keyed$sample_index, keyed$phase_set, keyed$haplotype_lane), ]
  rownames(keyed) <- NULL
  keyed$phase_set <- as.numeric(keyed$phase_set)
  status <- c(eligible = "eligible_classifier_pending", predicted = "predicted", incomplete = "incomplete_input",
    conflict = "edit_conflict", overlap = "unsupported_overlap", context = "unsupported_context")
  expected <- data.frame(
    sample_index = c(1, 2, 2, 3, 3, 4, 5, 5, 7, 8, 12, 16, 17, 21:26, 30, 40, 50),
    phase_set = c(rep(NA, 6), 10, 20, rep(NA, 14)),
    haplotype_lane = c(1, 1, 2, 1, 2, 2, 1, 1, rep(1, 14)),
    prediction_status = unname(status[c("predicted", "predicted", "predicted", "predicted", "incomplete",
      "predicted", "incomplete", "incomplete", "overlap", "conflict", "context",
      "predicted", "predicted", rep("context", 4), "predicted", "predicted", "predicted", "predicted", "predicted")]),
    prediction_reason = c("supported_domain", "supported_domain", "supported_domain",
      "supported_domain", "missing_call", "supported_domain",
      "unresolved_cross_ps_phase", "unresolved_cross_ps_phase", "overlapping_edits",
      "contradictory_edits", "outside_cds", "supported_domain", "supported_domain",
      "noncanonical_start", "internal_stop", "noncanonical_stop", "untyped_curated_metadata",
      "supported_domain", "supported_domain", "supported_domain", "supported_domain",
      "supported_domain"),
    stringsAsFactors = FALSE)
  # Sample 40 carries ALT ordinal 2 on lane 1 only.
  expected$phase_set <- as.numeric(expected$phase_set)
  expect_equal(keyed[, names(expected)], expected, check.attributes = FALSE)

  # Row summaries never claim eligibility or a classification for an ineligible carrier.
  for (i in seq_len(nrow(result))) {
    per_carrier <- result$carrier_predictions[[i]]
    expect_identical(nrow(per_carrier), nrow(result$carriers[[i]]))
    expect_identical(per_carrier$sample_index, result$carriers[[i]]$sample_index)
    ok <- c("eligible_classifier_pending", "predicted")
    expect_identical(result$prediction_status[i] %in% ok, all(per_carrier$prediction_status %in% ok))
  }

  # Same-codon and frame/stop-gain classifiers: the combined haplotype is classified, and the row-level SO
  # set/IMPACT are NULL unless the row is predicted. A cis SNV pair and an MNV in one codon give the same
  # missense; a trans pair is two paths; a wobble SNV is synonymous (LOW); a created stop before the
  # terminator is stop_gained (HIGH) and the SNV after it stays a contributor (slice 4 transition: this path
  # was eligible_classifier_pending at slice 3), so nothing in this fixture is pending any more.
  so <- vapply(result$haplotype_consequences, function(x) if (is.null(x)) NA_character_ else
    paste(x, collapse = ","), "")
  expect_identical(is.na(so), result$prediction_status != "predicted")
  expect_identical(is.na(result$haplotype_impact), result$prediction_status != "predicted")
  expect_true(all(so[result$cds == "ATGATTGCTGAAGGTTAA" & !is.na(result$cds)] == "missense_variant"))
  expect_true(all(result$haplotype_impact[!is.na(so) & so == "missense_variant"] == "MODERATE"))
  syn <- which(so == "synonymous_variant")
  expect_equal(length(syn), 1L)
  expect_identical(result$cds[syn], "ATGGCTGCCGAAGGTTAA")
  expect_identical(result$haplotype_impact[syn], "LOW")
  stopped <- which(so %in% "stop_gained")
  expect_equal(length(stopped), 1L)
  expect_identical(result$haplotype_impact[stopped], "HIGH")
  expect_equal(nrow(result$contributor_provenance[[stopped]]), 2L)
  expect_equal(sum(result$prediction_status == "eligible_classifier_pending"), 0L)

  # NMD (slice 6, rule ejc50). Transcript 0's exons are cDNA 1-9 and 10-18, so J = 9. The created stop TAG is
  # codon 2 (S = 6): J - S = 3, escape, with the stop's own edit attributed. unknown is a row that is not predicted or
  # a start loss; every decided row without a premature stop is not_applicable and has no positions or attribution.
  expect_true(all(result$nmd_rule == "ejc50"))
  expect_identical(result$nmd_prediction[stopped], "escape")
  expect_equal(c(result$nmd_stop_position[stopped], result$nmd_junction_position[stopped]), c(6, 9))
  expect_equal(length(result$nmd_contributors[[stopped]]), 1L)
  prov <- result$contributor_provenance[[stopped]]
  expect_identical(prov$role[prov$event_index %in% result$nmd_contributors[[stopped]]], "applied")
  expect_identical(result$nmd_prediction[result$prediction_status != "predicted"], rep("unknown", sum(result$prediction_status != "predicted")))
  expect_identical(result$nmd_prediction[result$prediction_status == "predicted" & seq_len(nrow(result)) != stopped &
    !grepl("stop_lost|start_lost", so)], rep("not_applicable", sum(result$prediction_status == "predicted" &
    seq_len(nrow(result)) != stopped & !grepl("stop_lost|start_lost", so))))
  expect_true(all(is.na(result$nmd_stop_position[-stopped])) && all(is.na(result$nmd_junction_position[-stopped])))
  expect_true(all(vapply(result$nmd_contributors[-stopped], is.null, NA)))

  # Per-carrier consequences (slice 4): every carrier whose own keyed status is predicted carries the shared
  # edited sequence's set and IMPACT, so an ineligible carrier of a row does not hide it; other carriers are NULL.
  for (i in seq_len(nrow(result))) {
    per_carrier <- result$carrier_predictions[[i]]
    expect_true(all(c("haplotype_impact", "haplotype_consequences") %in% names(per_carrier)))
    ok <- per_carrier$prediction_status == "predicted"
    expect_true(all(c("nmd_prediction", "nmd_stop_position", "nmd_junction_position") %in% names(per_carrier)))
    expect_true(all(per_carrier$nmd_prediction[!ok] == "unknown") && all(is.na(per_carrier$nmd_stop_position[!ok])))
    if (all(ok)) expect_true(all(per_carrier$nmd_prediction == result$nmd_prediction[i]))
    expect_identical(is.na(per_carrier$haplotype_impact), !ok) # every predicted path here has a non-empty set
    expect_identical(vapply(per_carrier$haplotype_consequences, is.null, NA), !ok)
    if (all(ok)) {
      expect_true(all(vapply(per_carrier$haplotype_consequences, identical, NA, result$haplotype_consequences[[i]])))
      expect_true(all(is.na(per_carrier$haplotype_impact) | per_carrier$haplotype_impact == result$haplotype_impact[i]))
    }
  }
  mixed <- which(vapply(result$carrier_predictions, function(x) any(x$prediction_status == "predicted") &&
    any(x$prediction_status != "predicted"), NA))
  expect_true(length(mixed) > 0L)
  expect_true(all(is.na(result$haplotype_impact[mixed])))
  expect_true(all(vapply(mixed, function(i) any(!is.na(result$carrier_predictions[[i]]$haplotype_impact)), NA)))

  # Every contributor keeps its operands, evidence, ALT ordinal, role and normalized edits.
  for (i in seq_len(nrow(result))) {
    old <- result$contributors[[i]]
    prov <- result$contributor_provenance[[i]]
    expect_identical(sort(prov$event_index), sort(old$event_index))
    expect_identical(prov$reference[order(prov$event_index)], old$reference[order(old$event_index)])
    expect_identical(prov$alternate[order(prov$event_index)], old$alternate[order(old$event_index)])
    expect_true(all(result$normalized_edits[[i]]$event_index %in% prov$event_index))
    expect_equal(nrow(result$normalized_edits[[i]]), result$edit_count[i])
  }
  roles <- function(sample, lane = 1L) {
    i <- which(vapply(result$carriers, function(x) any(x$sample_index == sample &
      x$haplotype_lane == lane) && nrow(x) == 1L, NA))
    expect_equal(length(i), 1L)
    prov <- result$contributor_provenance[[i]]
    paste(prov$event_index[order(prov$event_index)], prov$alt_index[order(prov$event_index)],
      prov$role[order(prov$event_index)], sep = ":", collapse = " ")
  }
  expect_identical(roles(1), "1:1:applied 2:1:applied")
  expect_identical(roles(7), "1:1:unapplied 5:1:unapplied")
  expect_identical(roles(8), "1:1:unapplied 6:1:unapplied")
  expect_identical(roles(12), "3:1:unapplied 11:1:omitted")
  expect_identical(roles(16), "3:1:post_stop 10:1:applied")
  expect_identical(roles(40), "20:2:applied")
  # A conflicting leaf still lists its edits, ascending and without a block.
  conflict <- result$normalized_edits[[which(result$prediction_reason == "contradictory_edits")]]
  expect_identical(conflict$edit_index, c(0, 1))
  expect_true(all(is.na(conflict$block_index)))

  # Reference lanes are not leaves in strict decoded input; the same sequence from an MNV
  # and from a cis SNV pair stays two leaves that keep their own contributors.
  lanes <- do.call(rbind, result$carriers)
  expect_equal(sum(lanes$sample_index == 1 & lanes$haplotype_lane == 2), 0L)
  same <- which(result$cds == "ATGATTGCTGAAGGTTAA")
  expect_equal(sort(vapply(same, function(i) nrow(result$contributor_provenance[[i]]), 0L)), c(1L, 2L))

  # vep_compat is a compatibility interpretation, outside the strict v1 domain.
  compat <- rduckvep_haplotypes(con, "SELECT * FROM elig_calls WHERE sample_index IN (1,2,3,4)",
    "elig", phase_policy = "vep_compat")
  ck <- do.call(rbind, compat$carrier_predictions)
  expect_true(all(ck$prediction_reason[ck$sample_index %in% c(1, 2, 4)] == "non_strict_phase_policy"))
  expect_true(all(ck$prediction_reason[ck$sample_index == 3] == "missing_call"))
  dbExecute(con, "CREATE TABLE elig_raw AS SELECT 1::UBIGINT event_index,0::UINTEGER seq_region,103::UBIGINT AS position,'G' AS reference,['A'] AS alternates,0::UINTEGER transcript_index,0::UINTEGER sample_index,'1|0' AS gt")
  rawres <- rduckvep_haplotypes(con, "SELECT * FROM elig_raw", "elig",
    phase_policy = "vep_compat", input_mode = "source_records")
  expect_true(all(rawres$prediction_status == "unsupported_context"))
  expect_true(all(rawres$prediction_reason == "non_strict_phase_policy"))
  expect_true(all(mapply(function(a, b) nrow(a) == nrow(b), rawres$contributors,
    rawres$contributor_provenance)))
  expect_true(any(rawres$edit_count == 0L))

  # Malformed identities and exhausted budgets fail explicitly instead of truncating.
  expect_error(rduckvep_haplotypes(con,
    "SELECT * REPLACE (NULL AS alt_index) FROM elig_calls WHERE sample_index=1", "elig"),
    "required input column 6 is NULL")
  expect_error(rduckvep_haplotypes(con, "SELECT * FROM elig_calls WHERE sample_index=1", "elig",
    max_leaf_events = 1), "max_leaf_events=1 exhausted")
})

# DuckDB vector boundaries and table-backed output: 2,100 transcripts with two SNVs each give
# 4,200 input rows and 2,800 output rows. Kind i %% 3: 0 cis pair, 1 unphased heterozygous
# pair (incomplete_input), 2 trans pair (two leaves).
local({
  con <- rduckvep_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,(100+30*i)::UBIGINT transcript_start,",
    "(117+30*i)::UBIGINT transcript_end,1::TINYINT strand,i::UINTEGER gene_index,3::UBIGINT transcript_flags,",
    "(100+30*i)::UBIGINT cds_start,(117+30*i)::UBIGINT cds_end,'ATGGCTGCTGAAGGTTAA'::BLOB cds_sequence,",
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2100) t(i)")
  exons <- paste("SELECT i::UINTEGER transcript_index,(100+30*i)::UBIGINT exon_start,(117+30*i)::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,18::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM range(2100) t(i)")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('elig_wide',",
    "'SELECT 0::UINTEGER seq_region',", dbQuoteString(con, tx), ",", dbQuoteString(con, exons), ")"))$loaded)
  calls <- paste("SELECT 2*i+j+1 AS event_index,0 AS seq_region,103+30*i+j AS position,",
    "CASE WHEN j=0 THEN 'G' ELSE 'C' END AS reference,CASE WHEN j=0 THEN 'A' ELSE 'T' END AS alternate,",
    "1 AS alt_index,i AS transcript_index,0 AS sample_index,",
    "CASE WHEN i%3=1 THEN [0,1] WHEN i%3=2 AND j=1 THEN [0,1] ELSE [1,0] END AS alleles,",
    "CASE WHEN i%3=1 THEN [false,false] ELSE [false,true] END AS phase_before,NULL::BIGINT AS phase_set",
    "FROM range(2100) t(i) CROSS JOIN range(2) e(j)")
  result <- rduckvep_haplotypes(con, calls, "elig_wide")
  expect_equal(nrow(result), 2800L)
  kind <- result$transcript_index %% 3L
  expect_equal(as.vector(table(kind)), c(700L, 700L, 1400L))
  expect_true(all(result$prediction_status[kind != 1L] == "predicted"))
  expect_true(all(is.na(result$haplotype_impact[kind == 1L])))
  expect_true(all(result$haplotype_impact[kind != 1L] == "MODERATE"))
  expect_true(all(vapply(result$haplotype_consequences[kind != 1L], identical, NA, "missense_variant")))
  expect_true(all(result$prediction_status[kind == 1L] == "incomplete_input"))
  expect_true(all(result$prediction_reason[kind == 1L] == "unphased_heterozygous"))
  for (i in seq_len(nrow(result))) {
    ids <- 2L * result$transcript_index[i] + 1:2
    expect_true(all(result$contributor_provenance[[i]]$event_index %in% ids))
    expect_identical(sort(result$contributor_provenance[[i]]$event_index),
      sort(result$contributors[[i]]$event_index))
    expect_equal(nrow(result$carrier_predictions[[i]]), result$carrier_count[i])
  }
  expect_equal(sum(vapply(result$contributor_provenance, nrow, 0L)), 4200L)
  expect_equal(sum(vapply(result$normalized_edits, nrow, 0L)), 2800L)

  rduckvep_haplotypes(con, calls, "elig_wide", table_name = "elig_wide_out")
  counts <- dbGetQuery(con, paste("SELECT count(*) n,sum(len(contributor_provenance)) prov,",
    "sum(len(carrier_predictions)) carriers,count(*) FILTER (WHERE prediction_status='incomplete_input') bad",
    "FROM elig_wide_out"))
  expect_equal(unlist(counts, use.names = FALSE), c(2800, 4200, 3500, 700))
})
