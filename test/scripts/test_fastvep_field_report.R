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
  setup_begin <- which(rmd == "```{r setup, include=FALSE}")
  stopifnot(length(setup_begin) == 1L)
  setup_end <- which(seq_along(rmd) > setup_begin & rmd == "```")[[1L]]
  setup <- parse(text = rmd[seq.int(setup_begin + 1L, setup_end - 1L)])
  timing_helpers <- new.env(parent = globalenv())
  for (name in c("time_value", "elapsed_seconds", "read_gnu_time")) {
    definition <- Filter(function(expr) is.call(expr) && identical(expr[[1L]], quote(`<-`)) &&
      identical(expr[[2L]], as.name(name)), setup)
    stopifnot(length(definition) == 1L)
    eval(definition[[1L]], envir = timing_helpers)
  }
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
      vapply(seq_len(8L), function(i) strrep(as.character(i), 64L), character(1L)),
      c("input", "model", "fasta", "gff3", "cache", "fastvep", "fasta_index", "source_map")
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
    writeLines("synthetic successful fresh Cargo build", file.path(path, "build.log"))
    build <- c(binding = "cargo_fresh_release_locked_offline", source_commit = identity[["source_commit"]],
      cargo_lock_sha256 = strrep("b", 64L), toolchain = "1.98.1", rustc = "fixture rustc",
      cargo = "fixture cargo", rustflags = "-C target-cpu=native", command = "fixture cargo build",
      executable_sha256 = hashes[["fastvep"]], log = "build.log",
      log_sha256 = sha256(file.path(path, "build.log")), exit_status = "0")
    write_fields(build, file.path(path, "fastvep_build.tsv"), tab = TRUE)
    hashes <- c(hashes, fastvep_build = sha256(file.path(path, "fastvep_build.tsv")),
      fastvep_build_log = build[["log_sha256"]])
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
      fastvep_binding = "cargo_fresh_release_locked_offline",
      supplementary_providers = "none", distance = "5000", output_filesystem = "synthetic"
    )
    write_fields(metadata, file.path(path, "metadata.csv"))
    coverage <- list()
    for (i in seq_len(nrow(matrix))) {
      configuration <- matrix$configuration[[i]]
      group <- match(configuration, configurations)
      timing <- paste0(labels[[i]], ".time")
      writeLines(c(
        "User time (seconds): 0.8", "System time (seconds): 0.2",
        "Percent of CPU this job got: 100%", "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00",
        "Maximum resident set size (kbytes): 1024", "File system outputs: 0", "Exit status: 0"
      ), file.path(path, timing))
      observation <- data.frame(
        tool = sub("_.*$", "", configuration),
        output_contract = sub("^[^_]+_", "", configuration),
        threads = matrix$threads[[i]], run = matrix$run[[i]],
        row_count = 4095611L + group, bytes = (4095611L + group) * 100,
        sha256 = strrep(as.character(group), 64L),
        xor_hash = as.character(group), low32_sum = as.character(20L + group),
        high32_sum = as.character(40L + group), timing_file = timing,
        multiset_checked = TRUE, fingerprint_scope = "full_row"
      )
      write_csv(observation, file.path(path, paste0(labels[[i]], ".csv")))
      coverage[[i]] <- data.frame(observation[c("tool", "output_contract", "threads", "run")],
        scope = "final_output", input_sha256 = hashes[["input"]],
        output_sha256 = observation$sha256, source_map_sha256 = hashes[["source_map"]],
        source_alleles = "4095611", covered_alleles = "4095611",
        missing_alleles = "0", unknown_alleles = "0", ambiguous_alleles = "0")
    }
    write_csv(do.call(rbind, coverage), file.path(path, "allele_coverage.csv"))
    write_csv(data.frame(
      source_revision = revision, completed = TRUE,
      binding = "source_bound", observations = 30L, repetitions = 3L,
      allele_coverage = "verified",
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
      timing_helpers$read_gnu_time(path, engine, threads, run)
    }
    output <- tryCatch(capture.output(result <- eval(gate, envir = env)),
      error = function(error) {
        attr(error, "timing_reads") <- timing_reads
        stop(error)
      })
    list(value = result, output = output, env = env, timing_reads = timing_reads)
  }
  rejected <- character()
  check <- function(label, mutate, reseal = TRUE, expected_timing_reads = 0L) {
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
    if (!identical(attr(error, "timing_reads"), expected_timing_reads)) {
      stop("publication rejected mutation after an unexpected number of timing reads: ", label)
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
    "eligible_literal_alleles", "fastvep_source_revision", "fastvep_version",
    "fastvep_binding"
  )) {
    values <- c(
      binding = "diagnostic_unbound", input_id = "fastvep_cache_probe",
      input_records = "3", input_alt_alleles = "3", eligible_literal_alleles = "3",
      fastvep_source_revision = "7038e7c17708e7d2226149e78e0bb297bcc6d1d6",
      fastvep_version = "fastvep 0.2.0", fastvep_binding = "diagnostic_binary_unbound"
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
  for (field in c("binding", "source_commit", "executable_sha256", "log_sha256", "exit_status")) {
    check(paste0("build_", field), function(path) {
      target <- file.path(path, "fastvep_build.tsv")
      build <- utils::read.delim(target, colClasses = "character")
      build$value[build$field == field] <- "invalid"
      utils::write.table(build, target, sep = "\t", quote = FALSE, row.names = FALSE)
    })
  }
  check("missing_build_receipt", function(path) unlink(file.path(path, "fastvep_build.tsv")))
  check("missing_build_log", function(path) unlink(file.path(path, "build.log")))
  check("changed_build_log", function(path) writeLines("stale build", file.path(path, "build.log")))
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
  check("unverified_source_coverage", function(path) {
    edit_csv(path, "completion.csv", function(x) {
      x$allele_coverage <- "unverified"
      x
    })
  })
  check("missing_source_coverage", function(path) {
    edit_csv(path, "completion.csv", function(x) {
      x$allele_coverage <- NULL
      x
    })
  })
  check("verified_flag_without_evidence", function(path) {
    unlink(file.path(path, "allele_coverage.csv"))
  })
  for (field in c("input_sha256", "output_sha256", "source_map_sha256", "source_alleles", "covered_alleles",
      "missing_alleles", "unknown_alleles", "ambiguous_alleles", "scope")) {
    check(paste0("coverage_", field), function(path) {
      edit_csv(path, "allele_coverage.csv", function(x) {
        x[[field]][[1L]] <- if (field == "scope") "before_projection" else "1"
        x
      })
    })
  }
  check("duplicate_coverage_observation", function(path) {
    edit_csv(path, "allele_coverage.csv", function(x) rbind(x, x[1L, ]))
  })
  check("fewer_output_rows_than_covered_alleles", function(path) {
    edit_csv(path, paste0(labels[[1L]], ".csv"), function(x) {
      x$row_count <- "1"
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
  for (field in names(read_csv(file.path(valid_path, paste0(labels[[1L]], ".csv"))))) {
    check(paste0("missing_observation_column_", field), function(path) {
      for (label in labels) {
        edit_csv(path, paste0(label, ".csv"), function(x) {
          x[[field]] <- NULL
          x
        })
      }
    })
  }
  for (value in list("", NA_character_)) {
    label <- if (is.na(value)) "missing" else "empty"
    check(paste0(label, "_observation_hash"), function(path) {
      edit_csv(path, paste0(labels[[1L]], ".csv"), function(x) {
        x$sha256 <- value
        x
      })
    })
  }
  check("duplicate_observation_hash_column", function(path) {
    edit_csv(path, paste0(labels[[1L]], ".csv"), function(x) {
      x$extra <- x$sha256
      names(x)[names(x) == "extra"] <- "sha256"
      x
    })
  })
  for (file in c("completion.csv", "metadata.csv", "inputs.csv", "allele_coverage.csv")) {
    fields <- names(read_csv(file.path(valid_path, file)))
    if (file == "inputs.csv") fields <- setdiff(fields, "path")
    for (field in fields) {
      check(paste0("missing_", file, "_", field), function(path) {
        edit_csv(path, file, function(x) {
          x[[field]] <- NULL
          x
        })
      }, reseal = file != "completion.csv")
    }
  }
  for (status in c("1", "137", "missing", "signal")) {
    check(paste0("failed_command_", status), function(path) {
      timing <- file.path(path, paste0(sort(labels)[[1L]], ".time"))
      lines <- readLines(timing)
      if (status == "signal") {
        lines <- c("Command terminated by signal 9", lines)
      } else {
        lines <- lines[!grepl("^Exit status:", lines)]
        if (status != "missing") lines <- c(lines, paste("Exit status:", status))
      }
      writeLines(lines, timing)
    }, expected_timing_reads = 1L)
  }
  check("retained_failed_command", function(path) {
    failed <- file.path(root, "benchmarks/data/duckvep_fastvep",
      "field_contracts_incomplete_0d2bdcb/duckvep_native_tab17_1_1.time")
    stopifnot(file.copy(failed, file.path(path, paste0(sort(labels)[[1L]], ".time")), overwrite = TRUE))
  }, expected_timing_reads = 1L)
  check_cli <- function(label, args, message) {
    output <- file.path(directory, paste0("cli-", label))
    result <- suppressWarnings(system2("Rscript", shQuote(c(
      file.path(root, "benchmarks/benchmark_duckvep_fastvep_run.R"), "--output", output, args
    )), stdout = TRUE, stderr = TRUE))
    stopifnot(!is.null(attr(result, "status")), attr(result, "status") != 0L,
      any(grepl(message, result, fixed = TRUE)), !file.exists(output))
  }
  check_cli("missing-build", character(), "published observations require --fastvep-build-receipt")
  check_cli("invalid-build", c("--diagnostic", "--fastvep-build-receipt", file.path(valid_path, "metadata.csv")),
    "invalid FastVEP build receipt schema")
  check_cli("mismatched-binary", c("--diagnostic", "--fastvep-build-receipt", file.path(valid_path, "fastvep_build.tsv"),
    "--fastvep", file.path(valid_path, "metadata.csv")), "FastVEP executable differs from its build receipt")
  cat("Complete-field report gate: valid matrix rendered;", length(rejected), "mutation controls rejected\n")
  cat("Timing runner: missing/malformed build receipts and mismatched executables rejected\n")
}

main()
