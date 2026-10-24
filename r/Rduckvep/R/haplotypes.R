#' Replay Phased Transcript Haplotypes
#'
#' Consume a SELECT query of flat event-by-transcript-by-sample calls through
#' the bundled native replay stream. DuckDB derives complete phase-set domains
#' and sorts the calls. Each output row is one occupied shared path, with CDS,
#' protein, carrier keys and all contributing events. Incomplete calls and
#' projection/edit failures retain provenance without inventing sequence.
#'
#' `nominal_length_diff` is the signed sum of projected replacement ALT-minus-REF
#' lengths before clipping at the current CDS end. It is zero for reference-only
#' replay and NA without a CDS. Known and conditional paths retain this value;
#' overlapping replacements can produce a different rebuilt CDS length change.
#'
#' `coding_blocks` groups physical edits that share an alternate codon or
#' displace then restore the reading frame. Each block gives one-based reference
#' `cds_start`, transcript-oriented `reference`/`alternate` spans, zero-based
#' `alt_start0` in the rebuilt CDS, `length_change`, `sequence_flags` and
#' `event_indices`: one source event ID per physical edit in ascending CDS order.
#' An ID can repeat for several differing islands, within or across blocks;
#' the list length is the physical edit count.
#' Join IDs to `contributors` for raw alleles and input provenance.
#' Spans include retained bases between edits; they are not aligned
#' differences or HGVS normalization. An insertion has empty reference; a deletion
#' has empty alternate. Unknown sequences have NULL blocks, not an empty known
#' result. Blocks reuse `max_leaf_edits` capacity and the per-call workspace limit.
#'
#' `cds_differences` instead contains aligned differing runs: zero-based
#' `ref_start0`, `alt_start0`, and `alignment_start0`, with borrowed sequence spans
#' materialized as `reference` and `alternate`. An empty span denotes a gap.
#' The reference uses replay's uppercase DNA spelling; model bytes stay immutable.
#' Runs join only when adjacent columns have the same gap/non-gap type on both
#' sides. Indel-bearing paths use VEP's pure-Perl global-alignment score and
#' tie order; substitution-only paths compare corresponding positions. Repeated
#' sequence can therefore place a difference away from its contributing event.
#' Differences are not HGVS normalization or event-provenance reassignment.
#' Unknown sequences have NULL differences. `max_alignment_cells` bounds the
#' exact traceback band and `max_leaf_differences` bounds output runs; exhaustion
#' is an error, not approximate alignment or discarded differences.
#'
#' `protein_differences` has the same span fields and alignment mode, with
#' positions in amino acids. Its reference follows Ensembl's start-methionine,
#' terminal-stop and curated single-residue peptide-edit rules. Haplosaurus then
#' appends `*` only for an exact uppercase TAA/TAG/TGA raw-CDS suffix, including
#' nonstandard-table and partial-CDS cases. The alternate is the displayed
#' first-stop prefix without reference peptide edits. Unknown paths and reference
#' CDS shorter than one codon have NULL protein differences; identical known
#' proteins have an empty list. Both difference axes reuse the same native
#' scratch and are separately subject to the alignment-cell and run limits.
#'
#' `prediction_status`, `prediction_reason`, `haplotype_consequences`,
#' `haplotype_impact`, `nmd_prediction` and `nmd_exceptions` classify the whole
#' edited product. The rules and the authority each is checked against are in the
#' haplotype contract (`design/duckvep_haplotype_contract.md`); in short:
#'
#' * Records compose at their own reference coordinates. Records that agree compose
#'   once and keep their identities; two that claim one base or insertion site
#'   differently are `edit_conflict`.
#' * `strict` reads GT and PS as written. An unphased heterozygous call needs no phase
#'   when it is the sample's only heterozygous call on the transcript. A missing allele
#'   or undefined phase is `incomplete_input`. `vep_compat` reads GT slots in order, as
#'   VEP's haplotype mode does, and then classifies like `strict`.
#' * The edited CDS is translated with the transcript's genetic code, with curated
#'   translations and RNA editing as part of the reference. A lost stop reads through
#'   the stored 3' flank. `start_lost` stands alone; otherwise the terms are
#'   `stop_gained`, `stop_lost`, `frameshift_variant` while the frame stays displaced,
#'   `protein_altering_variant` or `synonymous_variant` for a restored frame, and
#'   missense, in-frame, stop-retained and start-retained terms. Changes after the
#'   first stop are contributors with role `post_stop`.
#' * A path with no whole product to rebuild (an edit crossing a CDS end or an exon
#'   boundary, an essential splice-site record, an in-frame stop in the reference CDS)
#'   is `predicted` with reason `per_record_terms`: the union of its records'
#'   per-variant terms, without NMD.
#' * `nmd_prediction` applies rule `ejc50` to a new premature stop: `trigger` when the
#'   stop lies more than 50 bases before the last exon junction, otherwise `escape`.
#'   `nmd_exceptions` adds `start_proximal` and `long_exon` without changing it, and
#'   `nmd_stop_position`, `nmd_junction_position` and `nmd_contributors` give the
#'   evidence.
#'
#' `carrier_predictions` repeats the prediction per carrier, so one ineligible carrier
#' (a triploid call, say) does not hide the others; the row-level fields need every
#' carrier to be predicted. `contributor_provenance` lists every source of the path
#' (applied, omitted, shadowed, unapplied and post-stop) with its original operands,
#' carrier evidence and role. Sequence deduplication never drops a contributor.
#'
#' With `hgvs = TRUE`, a multi-record path whose records are proven cis gets compound
#' `hgvsc` and `hgvsp` names for the joint allele, and a single-record path gets the
#' per-variant names. A status other than `ok` keeps NULL text with its reason. HGVS
#' buffers count toward `workspace_limit`; exhaustion is an error, never truncated text.
#'
#' `stop_in_displaced_frame` reports whether any of the first translated stop
#' codon's three bases overlaps a frame-displaced span of the rebuilt CDS.
#' Displacement starts at a frame-changing edit and ends after the alternate
#' bases of the restoring edit, or continues downstream if unrestored.
#' It is FALSE for no stop or a stop after frame restoration, and NA when
#' sequence is unavailable. It is a sequence fact, not a combined SO consequence
#' or a claim that a restored DNA frame rescues the protein.
#'
#' Input must contain one row per
#' `event_index`, `transcript_index`, `sample_index`, with columns `seq_region`,
#' `position`, `reference`, `alternate`, `alt_index`, `alleles`, `phase_before`
#' and nullable `phase_set`. Event indices identify individual ALT events; retain
#' their source-record mapping. Transcript ordinals belong to the named model.
#' Candidate selection is explicit in this input relation.
#'
#' With `input_mode = "source_records"`, the required columns are `event_index`,
#' `seq_region`, `position`, `reference`, `alternates` (a character list),
#' `transcript_index`, `sample_index`, and `gt` (original VCF text).
#' Here `event_index` identifies the whole source record. This mode requires
#' `vep_compat`: the pinned file profile consumes two slots and ignores PS.
#' Missing calls and undefined slots can yield `conditional` sequence with
#' evidence bit 8; this is not known phase or biological rescue. Contributor
#' `alt_index` is 0 for REF, a positive source ALT ordinal, or NA for an undefined
#' slot's full-REF deletion. Source ALT strings must be nonempty and nonmissing.
#' Projection failures still withhold sequence. Raw spelling must be retained at
#' ingestion; it cannot be reconstructed losslessly from decoded GT arrays.
#'
#' Preparation reads committed objects on the registry's retained connection;
#' caller-local temporary objects and uncommitted changes are not visible. One
#' preparation may run per registry at a time. Nested or concurrent preparation
#' returns a busy error; completed scans use independent native workspaces.
#'
#' @param con An open DuckDB connection with DuckVEP loaded.
#' @param table_name Optional output table name.
#' @param overwrite Whether to replace an existing output table.
#' @param calls_query One nonempty SELECT query supplying the call relation.
#' @param model_name Name of an already loaded DuckVEP model.
#' @param phase_policy Strict GT/PS interpretation, or `vep_compat` for the
#'   called-slot order of the pinned executable VEP release.
#'   Decoded missing calls remain incomplete; source-record input uses the
#'   pinned raw parser and explicitly conditional missing-slot interpretation.
#' @param input_mode `alt_events` for decoded per-ALT calls, or `source_records`
#'   for raw GT and complete source ALT lists under `vep_compat`.
#' @param hgvs Whether to request bounded protein HGVS for supported completed paths.
#' @param ... Named positive integer workspace capacities accepted by
#'   `duckvep_haplotypes`, such as `max_active_events`, `max_active_carriers`,
#'   `max_sequence_bases`, `max_ploidy`, `max_phase_sets`, and `workspace_limit`.
#' @return A data frame, or invisible `TRUE` when creating `table_name`.
#' @export
rduckvep_haplotypes <- function(con, calls_query, model_name,
                               phase_policy = c("strict", "vep_compat"),
                               ..., input_mode = c("alt_events", "source_records"),
                               hgvs = FALSE, table_name = NULL, overwrite = FALSE) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  .duckvep_check_table_target(con, table_name, overwrite)
  for (name in c("calls_query", "model_name")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(name, " must be one nonempty string", call. = FALSE)
    }
  }
  phase_policy <- match.arg(phase_policy)
  input_mode <- match.arg(input_mode)
  if (input_mode == "source_records" && phase_policy != "vep_compat") {
    stop("source_records requires phase_policy='vep_compat'", call. = FALSE)
  }
  if (!is.logical(hgvs) || length(hgvs) != 1L || is.na(hgvs)) {
    stop("hgvs must be TRUE or FALSE", call. = FALSE)
  }
  limits <- list(...)
  if (length(limits) && (is.null(names(limits)) || anyDuplicated(names(limits)) ||
      any(!grepl("^[a-z][a-z_]*$", names(limits))))) {
    stop("capacities must have unique SQL parameter names", call. = FALSE)
  }
  params <- list(phase_policy = sql_quote_string(con, phase_policy),
                 input_mode = sql_quote_string(con, input_mode), hgvs = if (hgvs) "TRUE" else "FALSE")
  for (name in names(limits)) {
    value <- limits[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 || value > 2^53 || value != floor(value)) {
      stop(name, " must be one positive, exactly representable integer", call. = FALSE)
    }
    params[[name]] <- format(value, scientific = FALSE, trim = TRUE)
  }
  query <- paste0("SELECT * FROM duckvep_haplotypes(", sql_quote_string(con, calls_query),
                  ",", sql_quote_string(con, model_name), build_param_str(params), ")")
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  .duckvep_create_table(con, table_name, query, overwrite)
  invisible(TRUE)
}

#' Find the Transcripts Whose Coding Sequence Events Overlap
#'
#' The fast discovery route into [rduckvep_haplotypes()]: one native lookup in the
#' resident transcript interval index and exon layout instead of annotating every
#' record. Each row of `events_query` is one VCF record and one ALT allele with its
#' shared anchor base (`reference`, `alternate`, one-based `position`); the event is
#' normalized exactly as the annotation builder normalizes it, so the pairs returned
#' are exactly those to which [rduckvep_annotate()] assigns the CDS region bit: an
#' indel's anchor base alone never creates an overlap, an insertion is placed on the
#' base the annotation builder places it on (and is not a CDS event at the outer CDS
#' edges), and an intron of at most 13 bases inside the CDS counts as coding. Symbolic
#' and breakend ALTs, missing ALTs and alleles identical to REF are not events. Records
#' with no coding overlap, the large majority of a genome, produce no rows, so a VCF
#' reduces to its coding records before any per-record work. The result is the
#' `transcript_index` column of a calls relation for [rduckvep_haplotypes()]: join it
#' with the records' genotypes.
#'
#' @param con An open DuckDB connection with DuckVEP loaded.
#' @param events_query One nonempty SELECT query supplying `event_index`,
#'   `seq_region`, `position`, `reference` and `alternate` (any integer types for the
#'   first three).
#' @param model_name Name of an already loaded DuckVEP model without wrapped
#'   circular objects.
#' @param table_name Optional output table name.
#' @param overwrite Whether to replace an existing output table.
#' @return A data frame of `event_index` and `transcript_index` (one row per pair,
#'   ascending by transcript within an event), or invisible `TRUE` when creating
#'   `table_name`.
#' @export
rduckvep_coding_transcripts <- function(con, events_query, model_name,
                                        table_name = NULL, overwrite = FALSE) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  .duckvep_check_table_target(con, table_name, overwrite)
  for (name in c("events_query", "model_name")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(name, " must be one nonempty string", call. = FALSE)
    }
  }
  query <- paste0("SELECT event_index, unnest(duckvep_coding_transcripts(",
                  sql_quote_string(con, model_name),
                  ", seq_region, \"position\", reference, alternate)) AS transcript_index FROM (",
                  events_query, ") AS events")
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  .duckvep_create_table(con, table_name, query, overwrite)
  invisible(TRUE)
}

#' Read the Coding Calls of a VCF or BCF File
#'
#' The fused native route into [rduckvep_haplotypes()]: one pass over a (bgzipped)
#' VCF or BCF file with the bundled HTSlib. Each record is mapped to the model's
#' `seq_region` by `seq_region_name`, each ALT allele goes through the same discovery
#' as [rduckvep_coding_transcripts()], and only a record that touches coding sequence
#' has its FORMAT parsed and its genotypes decoded. The other records, about 99\% of a
#' whole genome, cost one line read and a few interval lookups.
#'
#' The result is the calls relation that [rduckvep_haplotypes()] consumes: one row per
#' ALT allele, transcript and sample with `event_index` (record ordinal times 64 plus
#' the ALT ordinal minus one, where the ordinal counts every data record of the file),
#' `seq_region`, `position`, `reference`, `alternate`, `alt_index`, `transcript_index`,
#' `sample_index` (the sample column, from 0), `alleles` (missing alleles are NA),
#' `phase_before` (a leading `FALSE`, then one flag per separator, `|` being `TRUE`)
#' and `phase_set` (the integer `PS` of the sample; labels such as `PATMAT` and a
#' missing `PS` give NA). Records on contigs the model does not know, symbolic ALTs and
#' records with no ALT make no rows. A file without samples, without a `FORMAT/GT`
#' header line, a coding record without `GT`, more than 64 ALT alleles in a coding
#' record, an invalid `POS` and a model loaded without `seq_region_name` are errors.
#'
#' @param con An open DuckDB connection with DuckVEP loaded.
#' @param model_name Name of an already loaded DuckVEP model (loaded with
#'   `seq_region_name` in its regions query) without wrapped circular objects.
#' @param path Path of the VCF, bgzipped VCF or BCF file.
#' @param table_name Optional output table name.
#' @param overwrite Whether to replace an existing output table.
#' @return A data frame of the calls, or invisible `TRUE` when creating `table_name`.
#'   To feed [rduckvep_haplotypes()], create the table and pass
#'   `"SELECT * FROM <table_name>"`.
#' @export
rduckvep_coding_calls <- function(con, model_name, path, table_name = NULL, overwrite = FALSE) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  .duckvep_check_table_target(con, table_name, overwrite)
  for (name in c("model_name", "path")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop(name, " must be one nonempty string", call. = FALSE)
    }
  }
  query <- paste0("SELECT * FROM duckvep_coding_calls(", sql_quote_string(con, model_name), ", ",
                  sql_quote_string(con, path.expand(path)), ")")
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  .duckvep_create_table(con, table_name, query, overwrite)
  invisible(TRUE)
}
