#!/usr/bin/env Rscript

# Publication-only controls: synthetic receipts exercise the real Rmd gate.
# No genome input, extension, oracle, cache builder or benchmark is executed.
main <- function() {
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  rmd <- readLines(file.path(root, "benchmarks/benchmark_duckvep_fastvep.Rmd"))
  begin <- which(rmd == "```{r complete-field-observations}")
  stopifnot(length(begin) == 1L)
  end <- which(seq_along(rmd) > begin & rmd == "```")[[1L]]
  gate <- parse(text = rmd[seq.int(begin + 1L, end - 1L)])
  registry <- read.delim(file.path(root, "r/duckhtsbench/inst/benchmark_registry.tsv"),
    check.names = FALSE, colClasses = "character"
  )
  cache_row <- registry[registry$id == "fastvep_ensembl116_cache", ]
  stopifnot(nrow(cache_row) == 1L)
  pairs <- strsplit(strsplit(cache_row$supplier_identity, ";", fixed = TRUE)[[1L]],
    "=",
    fixed = TRUE
  )
  stopifnot(all(lengths(pairs) == 2L))
  identity <- setNames(
    vapply(pairs, `[[`, character(1L), 2L),
    vapply(pairs, `[[`, character(1L), 1L)
  )
  stopifnot(
    identity[["source_commit"]] == "18177c26a0d1d2419fe43c3e8f6d4a0b5c4a3eb6",
    identity[["version"]] == "0.3.0", identity[["cache_format"]] == "FSTVEP05",
    identity[["preparation"]] == "full_gff_hgvs", identity[["transcripts"]] == "646577"
  )
  directory <- tempfile("fastvep-field-report-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  sha256 <- function(path) digest::digest(file = path, algo = "sha256")
  write_csv <- function(value, path) utils::write.csv(value, path, row.names = FALSE)
  read_csv <- function(path) utils::read.csv(path, colClasses = "character")
  write_fields <- function(values, path, tab = FALSE) {
    values <- data.frame(field = names(values), value = unname(values))
    if (tab) {
      utils::write.table(values, path, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      write_csv(values, path)
    }
  }
  revision <- strrep("a", 40L)
  configurations <- c(
    "duckvep_operational17", "duckvep_native_tab17",
    "duckvep_vep_csq", "fastvep_native_tab17", "fastvep_vep_csq"
  )
  matrix <- expand.grid(
    configuration = configurations, threads = c(1L, 4L),
    run = 1:3, stringsAsFactors = FALSE
  )
  labels <- with(matrix, paste(configuration, threads, run, sep = "_"))
  stopifnot(length(labels) == 30L, !anyDuplicated(labels))

  fixture <- function(path) {
    dir.create(path, recursive = TRUE)
    hashes <- setNames(
      vapply(seq_len(7L), function(i) strrep(as.character(i), 64L), character(1L)),
      c("input", "model", "fasta", "gff3", "cache", "fastvep", "fasta_index")
    )
    cache <- c(
      source_commit = identity[["source_commit"]],
      executable_version = paste("fastvep", identity[["version"]]),
      cache_format = identity[["cache_format"]], preparation = "full_gff_hgvs",
      transcript_count = identity[["transcripts"]],
      cache_sha256 = unname(hashes[["cache"]]), executable_sha256 = unname(hashes[["fastvep"]]),
      gff3_sha256 = unname(hashes[["gff3"]]), fasta_sha256 = unname(hashes[["fasta"]]),
      fasta_index_sha256 = unname(hashes[["fasta_index"]])
    )
    cache_path <- file.path(path, "fastvep_cache_receipt.tsv")
    write_fields(cache, cache_path, tab = TRUE)
    hashes <- c(hashes, cache_receipt = sha256(cache_path))
    write_csv(data.frame(
      artifact = names(hashes), path = paste0("synthetic/", names(hashes)),
      sha256 = unname(hashes)
    ), file.path(path, "inputs.csv"))
    metadata <- c(
      source_revision = revision, binding = "source_bound",
      input_id = "variantkey_giab_hg002_v421", input_records = "4048342",
      input_alt_alleles = "4096123", eligible_literal_alleles = "4095611",
      fastvep_source_revision = identity[["source_commit"]],
      fastvep_version = paste("fastvep", identity[["version"]]),
      supplementary_providers = "none", distance = "5000", output_filesystem = "synthetic"
    )
    write_fields(metadata, file.path(path, "metadata.csv"))
    for (i in seq_len(nrow(matrix))) {
      configuration <- matrix$configuration[[i]]
      group <- match(configuration, configurations)
      timing <- paste0(labels[[i]], ".time")
      writeLines("synthetic timing fixture", file.path(path, timing))
      observation <- data.frame(
        tool = sub("_.*$", "", configuration),
        output_contract = sub("^[^_]+_", "", configuration),
        threads = matrix$threads[[i]], run = matrix$run[[i]],
        row_count = 100L + group, bytes = 1000L + group,
        xor_hash = as.character(group), low32_sum = as.character(20L + group),
        high32_sum = as.character(40L + group), timing_file = timing,
        multiset_checked = TRUE, fingerprint_scope = "full_row"
      )
      write_csv(observation, file.path(path, paste0(labels[[i]], ".csv")))
    }
    write_csv(data.frame(
      source_revision = revision, completed = TRUE,
      binding = "source_bound", observations = 30L, repetitions = 3L,
      artifacts_sha256 = ""
    ), file.path(path, "completion.csv"))
    invisible(path)
  }
  seal <- function(path) {
    files <- setdiff(list.files(path), c("artifacts.csv", "completion.csv"))
    write_csv(
      data.frame(
        file = files,
        sha256 = unname(vapply(file.path(path, files), sha256, character(1L)))
      ),
      file.path(path, "artifacts.csv")
    )
    completion <- read_csv(file.path(path, "completion.csv"))
    completion$artifacts_sha256 <- sha256(file.path(path, "artifacts.csv"))
    write_csv(completion, file.path(path, "completion.csv"))
  }
  edit_csv <- function(path, file, edit) {
    target <- file.path(path, file)
    write_csv(edit(read_csv(target)), target)
  }
  edit_metadata <- function(path, field, value) {
    edit_csv(path, "metadata.csv", function(x) {
      stopifnot(sum(x$field == field) == 1L)
      x$value[x$field == field] <- value
      x
    })
  }
  edit_cache <- function(path, field, value) {
    target <- file.path(path, "fastvep_cache_receipt.tsv")
    cache <- read.delim(target, colClasses = "character")
    stopifnot(sum(cache$field == field) == 1L)
    cache$value[cache$field == field] <- value
    utils::write.table(cache, target, sep = "\t", quote = FALSE, row.names = FALSE)
    edit_csv(path, "inputs.csv", function(x) {
      x$sha256[x$artifact == "cache_receipt"] <- sha256(target)
      x
    })
  }
  evaluate <- function(path) {
    env <- new.env(parent = globalenv())
    env$root <- root
    env$data_dir <- dirname(path)
    # These denominators are independent literals, not read from mutated receipts.
    env$source_records <- 4048342
    env$source_alt_alleles <- 4096123
    timing_reads <- 0L
    env$read_gnu_time <- function(path, engine, threads, run) {
      timing_reads <<- timing_reads + 1L
      stopifnot(file.exists(path))
      data.frame(engine, threads, run,
        elapsed_seconds = 1, user_seconds = 0.8,
        system_seconds = 0.2, cpu_percent = 100, maximum_rss_kib = 1024
      )
    }
    output <- tryCatch(capture.output(result <- eval(gate, envir = env)),
      error = function(error) {
        attr(error, "timing_reads") <- timing_reads
        stop(error)
      })
    list(value = result, output = output, env = env, timing_reads = timing_reads)
  }
  rejected <- character()
  check <- function(label, mutate, reseal = TRUE) {
    path <- fixture(file.path(directory, label, "field_contracts"))
    seal(path)
    mutate(path)
    if (reseal) seal(path)
    error <- tryCatch(
      {
        evaluate(path)
        NULL
      },
      error = function(error) error
    )
    if (!inherits(error, "error")) stop("publication accepted mutation: ", label)
    if (!identical(attr(error, "timing_reads"), 0L)) {
      stop("publication read timing data before rejecting mutation: ", label)
    }
    rejected <<- c(rejected, label)
  }

  valid_path <- fixture(file.path(directory, "valid", "field_contracts"))
  seal(valid_path)
  valid <- evaluate(valid_path)
  stopifnot(
    inherits(valid$value, "knitr_kable"), valid$timing_reads == 30L,
    nrow(valid$env$complete_rows) == 30L, nrow(valid$env$medians) == 10L,
    all(valid$env$medians$elapsed_seconds == 1),
    any(grepl("Measured source:", valid$output, fixed = TRUE))
  )

  check("artifact_bytes_changed", function(path) {
    writeLines("corrupted timing", file.path(path, paste0(labels[[1L]], ".time")))
  }, reseal = FALSE)
  check("manifest_bytes_changed", function(path) {
    edit_csv(path, "artifacts.csv", function(x) {
      x$sha256[[1L]] <- strrep("0", 64L)
      x
    })
  }, reseal = FALSE)
  check("completion_diagnostic", function(path) {
    edit_csv(path, "completion.csv", function(x) {
      x$binding <- "diagnostic_unbound"
      x
    })
  })
  for (field in c(
    "binding", "input_id", "input_records", "input_alt_alleles",
    "eligible_literal_alleles", "fastvep_source_revision", "fastvep_version"
  )) {
    values <- c(
      binding = "diagnostic_unbound", input_id = "fastvep_cache_probe",
      input_records = "3", input_alt_alleles = "3", eligible_literal_alleles = "3",
      fastvep_source_revision = "7038e7c17708e7d2226149e78e0bb297bcc6d1d6",
      fastvep_version = "fastvep 0.2.0"
    )
    check(paste0("metadata_", field), function(path) edit_metadata(path, field, values[[field]]))
  }
  for (field in c(
    "source_commit", "executable_version", "cache_format", "preparation",
    "transcript_count", "cache_sha256", "executable_sha256", "gff3_sha256", "fasta_sha256",
    "fasta_index_sha256"
  )) {
    values <- c(
      source_commit = "7038e7c17708e7d2226149e78e0bb297bcc6d1d6",
      executable_version = "fastvep 0.2.0", cache_format = "FSTVEP04",
      preparation = "probe_only", transcript_count = "644427",
      cache_sha256 = strrep("0", 64L), executable_sha256 = strrep("0", 64L),
      gff3_sha256 = strrep("0", 64L), fasta_sha256 = strrep("0", 64L),
      fasta_index_sha256 = strrep("0", 64L)
    )
    check(paste0("cache_", field), function(path) edit_cache(path, field, values[[field]]))
  }
  check("missing_observation", function(path) {
    unlink(file.path(path, paste0(labels[[1L]], ".csv")))
  })
  check("missing_timing", function(path) {
    unlink(file.path(path, paste0(labels[[1L]], ".time")))
  })
  check("wrong_timing_observation", function(path) {
    edit_csv(path, paste0(labels[[1L]], ".csv"), function(x) {
      x$timing_file <- paste0(labels[[2L]], ".time")
      x
    })
  })
  check("changed_repeat_fingerprint", function(path) {
    edit_csv(path, paste0(labels[[1L]], ".csv"), function(x) {
      x$xor_hash <- "999"
      x
    })
  })
  check("duplicate_metadata_key", function(path) {
    edit_csv(path, "metadata.csv", function(x) rbind(x, x[x$field == "binding", ]))
  })
  check("duplicate_metadata_header", function(path) {
    edit_csv(path, "metadata.csv", function(x) {
      names(x) <- c("field", "field")
      x
    })
  })
  check("duplicate_input_key", function(path) {
    edit_csv(path, "inputs.csv", function(x) rbind(x, x[1L, ]))
  })
  cat("Complete-field report gate: valid matrix rendered;", length(rejected), "mutation controls rejected\n")
}

main()
