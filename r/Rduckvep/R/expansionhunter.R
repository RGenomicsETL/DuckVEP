#' Prepare a single-repeat ExpansionHunter v5 VCF allele
#'
#' The supported contract is ExpansionHunter v5.0.0 INFO/REF, RL, RU, END,
#' symbolic `<STRn>` ALT, and FORMAT/GT, SO with either CN/CI or REPCN/REPCI. The caller supplies
#' the literal reference interval from POS+1 through END (excluding the VCF
#' padding base). An ALT is accepted only if a called, spanning allele has a
#' point confidence interval matching its symbolic copy number. Other records
#' retain their source and a stable preparation reason; no symbolic sequence is
#' treated as an observed literal allele.
#'
#' @param info VCF INFO text.
#' @param format VCF FORMAT keys.
#' @param sample VCF sample values.
#' @param ref VCF REF padding base.
#' @param alt Comma-separated VCF ALT field.
#' @param alt_index One-based selected ALT index.
#' @param reference_sequence Literal reference bases in the repeat interval.
#' @return A list with status, reason, source, reference_components,
#'   alternate_components, and sequence_exact. Components have columns
#'   `unit` and `count`, in reference order, for `duckvep_repeat_alleles`.
#' @export
rduckvep_prepare_expansionhunter <- function(info, format, sample, ref, alt,
                                             reference_sequence, alt_index = 1L) {
  source <- list(info = info, format = format, sample = sample, ref = ref,
                 alt = alt, reference_sequence = reference_sequence,
                 alt_index = alt_index)
  result <- function(status, reason, reference = NULL, alternate = NULL) {
    list(status = status, reason = reason, source = source,
         reference_components = reference, alternate_components = alternate,
         sequence_exact = identical(status, "ok"))
  }
  scalar <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  if (!all(vapply(source[seq_len(6L)], scalar, logical(1))))
    return(result("incomplete", "missing_field"))
  alts <- strsplit(alt, ",", fixed = TRUE)[[1L]]
  if (length(alt_index) != 1L || !is.numeric(alt_index) ||
      is.na(alt_index) || !is.finite(alt_index) ||
      alt_index < 1L || alt_index > length(alts) ||
      alt_index != floor(alt_index))
    return(result("invalid", "alt_index"))
  selected <- alts[[alt_index]]
  if (!grepl("^[ACGT]$", ref) || !grepl("^[ACGT]+$", reference_sequence))
    return(result("invalid", "literal_reference"))
  pairs <- strsplit(info, ";", fixed = TRUE)[[1L]]
  keys <- sub("=.*$", "", pairs)
  if (anyDuplicated(keys)) return(result("invalid", "duplicate_info"))
  values <- sub("^[^=]*=", "", pairs)
  names(values) <- keys
  required <- c("REF", "RL", "RU", "END")
  if (!all(required %in% keys)) return(result("incomplete", "missing_info"))
  unit <- values[["RU"]]
  number <- function(x) {
    if (!grepl("^(0|[1-9][0-9]*)$", x)) return(NA_real_)
    as.numeric(x)
  }
  reference_count <- number(values[["REF"]])
  reference_length <- number(values[["RL"]])
  if (!grepl("^[ACGT]+$", unit) || !is.finite(reference_count) ||
      !is.finite(reference_length) || !is.finite(number(values[["END"]])))
    return(result("invalid", "metadata_syntax"))
  if (reference_length > 5000)
    return(result("summary_only", "allele_capacity"))
  if (reference_count * nchar(unit) != reference_length ||
      reference_length != nchar(reference_sequence) ||
      !identical(strrep(unit, reference_count), reference_sequence))
    return(result("invalid", "reference_mismatch"))
  reference <- data.frame(unit = unit, count = reference_count)
  if (!grepl("^<STR(0|[1-9][0-9]*)>$", selected))
    return(result("summary_only", "unsupported_alt", reference))
  alternate_count <- number(sub("^<STR([0-9]+)>$", "\\1", selected))
  if (!is.finite(alternate_count))
    return(result("invalid", "metadata_syntax", reference))
  fields <- strsplit(format, ":", fixed = TRUE)[[1L]]
  entries <- strsplit(sample, ":", fixed = TRUE)[[1L]]
  if (length(fields) != length(entries) || anyDuplicated(fields))
    return(result("invalid", "format_shape", reference))
  names(entries) <- fields
  count_field <- if (all(c("CN", "CI") %in% fields)) c("CN", "CI") else
    if (all(c("REPCN", "REPCI") %in% fields)) c("REPCN", "REPCI") else character()
  if (!all(c("GT", "SO") %in% fields) || length(count_field) == 0L)
    return(result("incomplete", "missing_format", reference))
  if (all(c("CN", "CI", "REPCN", "REPCI") %in% fields))
    return(result("invalid", "ambiguous_count_fields", reference))
  gt <- strsplit(entries[["GT"]], "[/|]")[[1L]]
  copies <- strsplit(entries[[count_field[1L]]], "/", fixed = TRUE)[[1L]]
  ranges <- strsplit(entries[[count_field[2L]]], "/", fixed = TRUE)[[1L]]
  support <- strsplit(entries[["SO"]], "/", fixed = TRUE)[[1L]]
  if (length(gt) != length(copies) || length(gt) != length(ranges) ||
      length(gt) != length(support))
    return(result("invalid", "format_shape", reference))
  if (!all(gt %in% c(".", as.character(seq.int(0L, length(alts))))))
    return(result("invalid", "genotype_index", reference))
  called <- which(gt == as.character(alt_index))
  if (length(called) == 0L)
    return(result("summary_only", "uncalled_alt", reference))
  if (!all(vapply(copies[called], function(x) is.finite(number(x)) &&
                  number(x) == alternate_count, TRUE)))
    return(result("invalid", "count_mismatch", reference))
  if (!any(support[called] == "SPANNING" &
           ranges[called] == paste0(alternate_count, "-", alternate_count)))
    return(result("summary_only", "estimated_count", reference))
  if (alternate_count * nchar(unit) > 5000)
    return(result("summary_only", "allele_capacity", reference))
  result("ok", "exact", reference, data.frame(unit = unit, count = alternate_count))
}
