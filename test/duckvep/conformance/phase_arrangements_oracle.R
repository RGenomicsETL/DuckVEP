args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) > 0L) normalizePath(args[[1L]]) else normalizePath(".")
probe <- file.path(root, "test/duckvep/conformance/phase_arrangements_probe.c")
phase_source <- file.path(root, "src/kernel/src/duckvep_phase.c")
exe <- file.path(tempdir(), "duckvep_phase_arrangements_probe")
cc <- Sys.getenv("CC", unset = "cc")
status <- system2(cc, c("-std=c11", "-Wall", "-Wextra", "-Werror",
  "-I", file.path(root, "src/kernel/src"), probe, phase_source, "-o", exe))
if (status != 0L) stop("cannot compile phase arrangement probe", call. = FALSE)
observed <- strsplit(system2(exe, stdout = TRUE), "\t", fixed = TRUE)

make_sites <- function(ids, first, second, phased = rep(FALSE, length(ids)),
                       phase_set = rep(NA_real_, length(ids))) {
  data.frame(id = ids, first = first, second = second, phased = phased,
    phase_set = phase_set, stringsAsFactors = FALSE)
}

phase_units <- function(sites) {
  unit <- integer(nrow(sites))
  count <- 0L
  for (i in seq_len(nrow(sites))) {
    if (!sites$phased[[i]]) {
      count <- count + 1L
      unit[[i]] <- count
    } else {
      prior <- which(sites$phased[seq_len(i - 1L)] &
        ((is.na(sites$phase_set[seq_len(i - 1L)]) & is.na(sites$phase_set[[i]])) |
          (!is.na(sites$phase_set[seq_len(i - 1L)]) & !is.na(sites$phase_set[[i]]) &
            sites$phase_set[seq_len(i - 1L)] == sites$phase_set[[i]])))
      if (length(prior) > 0L) {
        unit[[i]] <- unit[[prior[[1L]]]]
      } else {
        count <- count + 1L
        unit[[i]] <- count
      }
    }
  }
  unit
}

exhaustive_arrangements <- function(sites) {
  unit <- phase_units(sites)
  unit_count <- max(unit)
  alternatives <- 2L ^ (unit_count - 1L)
  unlist(lapply(0L:(alternatives - 1L), function(value) {
    flips <- as.integer(intToBits(value))[seq_len(max(0L, unit_count - 1L))]
    lanes <- ifelse(sites$phased & sites$first == 1L, 1L, 2L)
    for (i in seq_len(nrow(sites))) {
      if (unit[[i]] > 1L && flips[[unit[[i]] - 1L]] == 1L) lanes[[i]] <- 3L - lanes[[i]]
    }
    paste(paste0(sites$id, ":", lanes), collapse = ",")
  }), use.names = FALSE)
}

cases <- list(
  unphased = make_sites(c(101L, 102L), c(0L, 1L), c(1L, 0L)),
  three = make_sites(c(111L, 112L, 113L), c(0L, 0L, 1L), c(1L, 1L, 0L)),
  block = make_sites(c(121L, 122L, 123L), c(0L, 1L, 1L), c(1L, 0L, 0L),
    c(TRUE, TRUE, FALSE), c(8, 8, NA)),
  separate = make_sites(c(131L, 132L), c(0L, 1L), c(1L, 0L), c(TRUE, TRUE), c(8, 9)),
  mixed_default = make_sites(c(141L, 142L, 143L), c(0L, 0L, 1L), c(1L, 1L, 0L),
    c(FALSE, TRUE, TRUE))
)
expected <- unlist(lapply(names(cases), function(name) {
  paste(name, seq_along(exhaustive_arrangements(cases[[name]])) - 1L,
    exhaustive_arrangements(cases[[name]]), sep = "\t")
}), use.names = FALSE)
actual <- vapply(observed, paste, character(1L), collapse = "\t")
if (!identical(actual, expected)) stop("kernel arrangements differ from independent exhaustive oracle", call. = FALSE)

codon_outcomes <- function(lanes) {
  edits <- list(`101` = list(position = 1L, base = "T"), `102` = list(position = 2L, base = "T"))
  # NCBI genetic code table 1: GAA=E, GTA=V, TAA=STOP, TTA=L.
  amino <- c(GAA = "E", GTA = "V", TAA = "*", TTA = "L")
  vapply(1:2, function(lane) {
    ids <- names(lanes)[lanes == lane]
    codon <- strsplit("GAA", "", fixed = TRUE)[[1L]]
    for (id in ids) codon[[edits[[id]]$position]] <- edits[[id]]$base
    paste0("lane", lane, ":", amino[[paste0(codon, collapse = "")]], "[",
      paste(ids, collapse = ","), "]")
  }, character(1L))
}
frame_outcomes <- function(lanes) {
  delta <- c(`201` = -1L, `202` = 1L)
  vapply(1:2, function(lane) {
    ids <- names(lanes)[lanes == lane]
    paste0("lane", lane, ":", sum(delta[ids]), "[", paste(ids, collapse = ","), "]")
  }, character(1L))
}

same_codon <- exhaustive_arrangements(make_sites(c(101L, 102L), c(0L, 0L), c(1L, 1L)))
same_codon_lanes <- lapply(strsplit(same_codon, ",", fixed = TRUE), function(row) {
  values <- as.integer(sub(".*:", "", row))
  names(values) <- sub(":.*", "", row)
  values
})
observed_codon <- vapply(same_codon_lanes, function(lanes) paste(codon_outcomes(lanes), collapse = ";"), character(1L))
expected_codon <- c("lane1:E[];lane2:L[101,102]", "lane1:V[102];lane2:*[101]")
if (!identical(observed_codon, expected_codon)) stop("same-codon cis/trans outcome or contributors changed", call. = FALSE)

frame <- exhaustive_arrangements(make_sites(c(201L, 202L), c(0L, 0L), c(1L, 1L)))
frame_lanes <- lapply(strsplit(frame, ",", fixed = TRUE), function(row) {
  values <- as.integer(sub(".*:", "", row))
  names(values) <- sub(":.*", "", row)
  values
})
observed_frame <- vapply(frame_lanes, function(lanes) paste(frame_outcomes(lanes), collapse = ";"), character(1L))
expected_frame <- c("lane1:0[];lane2:0[201,202]", "lane1:1[202];lane2:-1[201]")
if (!identical(observed_frame, expected_frame)) stop("compensating-frame cis/trans outcome or contributors changed", call. = FALSE)

message("phase arrangement oracle: ", length(actual), " kernel arrangements; cis/trans contributors verified")
