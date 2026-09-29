#' Prepare VEP 116 VCF structural geometry and insertion provenance
#'
#' For symbolic alleles the VCF padding base is outside the affected interval.
#' Confidence bounds are VEP 116's inner/outer coordinates; annotation uses
#' nominal start/end. The supported symbolic contract has a one-base A/C/G/T/N
#' padding REF, an integral END, and coordinates within signed 32-bit range. A literal anchored insertion takes the small-variant path
#' and carries its inserted bases; INFO/SEQ on a symbolic INS is provenance only.
#'
#' @param pos One-based VCF POS.
#' @param ref VCF REF.
#' @param alt One VCF ALT.
#' @param info VCF INFO text.
#' @return A one-row data frame with status, mode, nominal and confidence
#'   coordinates, inserted_sequence, and source_sequence. Join by the caller's
#'   event identifier to annotations produced by `duckvep_annotate_sql`.
#' @export
rduckvep_prepare_sv_geometry <- function(pos, ref, alt, info) {
  if (length(pos) != 1L || !is.numeric(pos) || is.na(pos) ||
      !is.finite(pos) || pos < 1 || pos != floor(pos) || pos > 2147483647 ||
      length(ref) != 1L || is.na(ref) || length(alt) != 1L || is.na(alt) ||
      length(info) != 1L || is.na(info))
    stop("invalid VCF coordinate or allele", call. = FALSE)
  pairs <- strsplit(info, ";", fixed = TRUE)[[1L]]
  keys <- sub("=.*$", "", pairs)
  if (anyDuplicated(keys)) stop("duplicate INFO key", call. = FALSE)
  values <- sub("^[^=]*=", "", pairs)
  names(values) <- keys
  scalar <- function(key) if (key %in% keys) values[[key]] else NA_character_
  sequence <- scalar("SEQ")
  symbolic <- grepl("^<[^>]+>$", alt)
  literal <- grepl("^[ACGT]+$", ref) && grepl("^[ACGT]+$", alt) &&
    nchar(alt) > nchar(ref) && startsWith(alt, ref)
  mode <- if (symbolic) "structural" else if (literal) "literal_insertion" else "unsupported"
  start <- pos + nchar(ref)
  end_text <- scalar("END")
  end <- if (symbolic && !is.na(end_text)) {
    suppressWarnings(as.numeric(end_text))
  } else start - 1
  interval <- function(key, origin) {
    if (!key %in% keys) return(c(NA_real_, NA_real_))
    offsets <- strsplit(scalar(key), ",", fixed = TRUE)[[1L]]
    if (length(offsets) != 2L || !all(grepl("^[+-]?[0-9]+$", offsets)))
      return(c(NA_real_, NA_real_))
    origin + as.numeric(offsets)
  }
  cipos <- interval("CIPOS", start)
  ciend <- interval("CIEND", end)
  bounds_ok <- function(bounds, key) {
    !key %in% keys || (all(is.finite(bounds) & bounds >= 0 &
      bounds <= 2147483647) && bounds[1L] <= bounds[2L])
  }
  valid <- mode != "unsupported" &&
    (!symbolic || (grepl("^[ACGTN]$", ref) &&
      !is.na(end_text) && grepl("^[1-9][0-9]*$", end_text))) &&
    is.finite(end) && end >= start - 1 &&
    start <= 2147483647 && end <= 2147483647 &&
    bounds_ok(cipos, "CIPOS") && bounds_ok(ciend, "CIEND")
  data.frame(status = if (valid) "ok" else "unsupported_geometry", mode = mode,
    nominal_start = start, nominal_end = end,
    outer_start = cipos[1L], inner_start = cipos[2L],
    inner_end = ciend[1L], outer_end = ciend[2L],
    inserted_sequence = if (literal) substring(alt, nchar(ref) + 1L) else NA_character_,
    source_sequence = sequence)
}
