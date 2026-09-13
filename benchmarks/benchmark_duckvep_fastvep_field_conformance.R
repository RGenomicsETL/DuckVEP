#!/usr/bin/env Rscript

# Exact field comparisons over the existing eight projection fixture models.
# Original records, labelled comparison copies and every discrepancy are retained.
duckvep_fastvep_fixture_model <- function(con, gff, reference, output) {
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  execute <- function(sql) invisible(DBI::dbExecute(con, sql))
  if (file.exists(output)) stop("fixture model output already exists")
  execute(paste0(
    "CREATE OR REPLACE TEMP TABLE field_fixture_gff AS SELECT feature,
    attributes_map AS attributes FROM read_gff(", q(gff),
    ", attributes_map := TRUE, scan_mode := 'sequential')"
  ))
  index <- utils::read.delim(paste0(reference, ".fai"), header = FALSE)
  DBI::dbWriteTable(con, "field_fixture_fai", data.frame(
    name = index[[1L]],
    sequence_length = index[[2L]]
  ), temporary = TRUE, overwrite = TRUE)
  execute(paste0("ATTACH ", q(output), " AS field_fixture"))
  on.exit(execute("DETACH field_fixture"), add = TRUE)
  execute("CREATE TABLE field_fixture.model_regions AS
    SELECT seq_region, name AS seq_region_name, sequence_length::UBIGINT AS sequence_length
    FROM duckvep_sequence_regions JOIN field_fixture_fai USING(name)")
  # CDS feature IDs are not inferred protein accessions. Only an explicit,
  # unambiguous protein_id attribute supplies the cold translation identifier.
  execute("CREATE TABLE field_fixture.model_transcripts AS
    WITH proteins AS (
      SELECT regexp_replace(attributes['Parent'], '^transcript:', '') AS transcript_id,
        CASE WHEN count(DISTINCT attributes['protein_id']) = 1
          THEN min(attributes['protein_id']) END AS protein_id
      FROM field_fixture_gff WHERE feature = 'CDS' GROUP BY transcript_id
    ) SELECT p.*, n.transcript_id AS transcript_stable_id,
      regexp_replace(g.attributes['Parent'], '^gene:', '') AS gene_stable_id,
      try_cast(g.attributes['version'] AS BIGINT) AS transcript_version,
      g.attributes['biotype'] AS transcript_biotype,
      proteins.protein_id AS translation_stable_id, NULL::BIGINT AS translation_version,
      NULL::VARCHAR AS mane_select_refseq, NULL::VARCHAR AS mane_plus_clinical_refseq,
      []::STRUCT(mature_mirna_start UBIGINT, mature_mirna_end UBIGINT)[] AS mature_mirna_regions
    FROM projection_transcripts p JOIN duckvep_transcript_names n USING(transcript_index)
    JOIN field_fixture_gff g ON regexp_replace(g.attributes['ID'], '^transcript:', '') = n.transcript_id
    LEFT JOIN proteins ON proteins.transcript_id = n.transcript_id")
  counts <- DBI::dbGetQuery(con, "SELECT
    (SELECT count(*) FROM projection_transcripts) AS expected,
    count(*) AS actual FROM field_fixture.model_transcripts")
  stopifnot(counts$expected == counts$actual)
  invisible(output)
}

duckvep_fastvep_prepare_source <- function(con, raw, expected_records = NULL) {
  DBI::dbExecute(con, paste0("CREATE TEMP TABLE field_input AS
    SELECT row_number() OVER ()::UBIGINT AS record_index, 1::UBIGINT AS alt_index,
      ID AS Uploaded_variation, CHROM, POS, REF, ALT FROM ", raw))
  counts <- DBI::dbGetQuery(con, "SELECT count(*) AS records,
    count(*) FILTER (WHERE Uploaded_variation IS NULL OR Uploaded_variation IN ('', '.')
      OR ALT IS NULL OR ALT IN ('', '.') OR contains(ALT, ',')) AS invalid FROM field_input")
  duplicates <- DBI::dbGetQuery(con, "SELECT count(*) n FROM
    (SELECT Uploaded_variation FROM field_input GROUP BY ALL HAVING count(*) != 1)")$n
  if (counts$records == 0 || counts$invalid != 0 || duplicates != 0) {
    stop("source VCF requires nonempty biallelic records with unique nonmissing IDs")
  }
  if (!is.null(expected_records) && counts$records != expected_records) {
    stop("source VCF record count differs from the generated physical-record count")
  }
  counts$records
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--output", default = ""),
    optparse::make_option("--seed", type = "integer", default = 173L),
    optparse::make_option("--random-cases", dest = "random_cases", type = "integer", default = 1000L),
    optparse::make_option("--replay-input", dest = "replay_input", default = "",
      help = "retained unique-ID biallelic VCF; requires --case and bypasses generation and relabelling"),
    optparse::make_option("--case", default = "", help = "existing fixture case for --replay-input"),
    optparse::make_option("--extension", default = "build/release/duckhts.duckdb_extension"),
    optparse::make_option("--extension-receipt", dest = "extension_receipt", default = NULL),
    optparse::make_option("--vep-prefix", dest = "vep_prefix", default = Sys.getenv("VEP_PREFIX")),
    optparse::make_option("--fastvep", default = ".sync/fastVEP/target/release/fastvep"),
    optparse::make_option("--fastvep-build-receipt",
      dest = "fastvep_build_receipt", default = NULL,
      help = "fresh pinned-source build receipt; required for publishable evidence"
    )
  )))
  stopifnot(!is.na(opt$seed), !is.na(opt$random_cases), opt$random_cases >= 0L)
  replay <- nzchar(opt$replay_input)
  if (replay != nzchar(opt$case)) stop("--replay-input and --case must be supplied together")
  if (!nzchar(opt$vep_prefix)) stop("set VEP_PREFIX or provide --vep-prefix")
  if (!is.null(opt$extension_receipt) && is.null(opt$fastvep_build_receipt)) {
    stop("publishable evidence requires --fastvep-build-receipt for the selected executable artifact")
  }
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  imports <- c(
    "scripts/duckvep_evidence.R", "test/duckvep/conformance/projection_fixtures.R",
    "benchmarks/benchmark_duckvep_fastvep_fields.R", "benchmarks/benchmark_duckvep_fastvep_extract.R",
    "benchmarks/benchmark_duckvep_fastvep_compare.R",
    "r/duckhtsbench/R/duckvep.R", "r/duckhtsbench/R/fastvep.R"
  )
  for (script in imports) source(file.path(root, script), local = TRUE)
  if (replay && !opt$case %in% duckvep_projection_cases) stop("unknown replay fixture case: ", opt$case)
  cases <- if (replay) opt$case else duckvep_projection_cases
  if (replay) opt$replay_input <- normalizePath(opt$replay_input, mustWork = TRUE)
  sources <- file.path(root, c(
    imports, "benchmarks/benchmark_duckvep_fastvep_field_conformance.R",
    "benchmarks/benchmark_duckvep_fastvep_worker.R", "test/duckvep/conformance/generate_witnesses.R",
    "test/duckvep/conformance/minimal_model.sql"
  ))
  source_hashes <- vapply(sources, duckvep_evidence_sha256, character(1L))
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  fastvep <- normalizePath(opt$fastvep, mustWork = TRUE)
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  revision <- duckvep_evidence_revision(root)
  binding <- "diagnostic_unbound"
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_command(
      "git", c("-C", root, "ls-files", "--error-unmatch", sources),
      "publishable conformance requires checked-in driver sources"
    )
    duckvep_evidence_assert_checkout(root, revision, allowed_outputs = opt$output)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)$binding
  }
  if (!nzchar(opt$output)) {
    label <- if (replay) paste0("fastvep_fields_replay_", opt$case, "_") else paste0("fastvep_fields_seed", opt$seed, "_")
    opt$output <- tempfile(label,
      tmpdir = file.path(root, "test/duckvep/conformance/results")
    )
  }
  directory <- normalizePath(opt$output, mustWork = FALSE)
  if (dir.exists(directory) || !dir.create(directory, recursive = TRUE)) stop("output must be a new directory")
  message("Artifacts: ", directory)
  command <- function(exe, args) duckvep_evidence_command(exe, args, paste("failed:", exe))
  pins <- c(
    vep = "57ea5c52340acc1f156267f810ad162e26597082",
    variation = "2fb834b987ede3824e200197a838ce11e91aeb4b",
    fastvep = "18177c26a0d1d2419fe43c3e8f6d4a0b5c4a3eb6"
  )
  mirrors <- file.path(root, ".sync", c("ensembl-vep", "ensembl-variation", "fastVEP"))
  names(mirrors) <- names(pins)
  check_pins <- function() {
    for (name in names(pins)) {
      stopifnot(
        identical(command("git", c("-C", mirrors[[name]], "rev-parse", "HEAD")), pins[[name]]),
        !length(command("git", c("-C", mirrors[[name]], "status", "--porcelain", "--untracked-files=no")))
      )
    }
  }
  check_pins()
  environment <- command("micromamba", c("list", "-p", prefix, "--explicit"))
  lock <- file.path(root, "test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt")
  stopifnot(identical(
    duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines(lock))
  ))
  writeLines(environment, file.path(directory, "environment.txt"))
  extension_hash <- duckvep_evidence_sha256(extension)
  fastvep_hash <- duckvep_evidence_sha256(fastvep)
  fastvep_binding <- "diagnostic_binary_unbound"
  if (!is.null(opt$fastvep_build_receipt)) {
    build <- duckhts_bench_read_fastvep_build(opt$fastvep_build_receipt, pins[["fastvep"]], fastvep)
    fastvep_binding <- build[["binding"]]
    build_files <- c(opt$fastvep_build_receipt,
      file.path(dirname(opt$fastvep_build_receipt), build[["log"]]))
    if (!all(file.copy(build_files, file.path(directory, c("fastvep_build.tsv", build[["log"]]))))) {
      stop("could not retain FastVEP build provenance")
    }
  }
  inputs <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root, "duckvep-projection")
  if (replay) inputs <- c(inputs, replay_input = opt$replay_input)
  snapshot_inputs <- function() {
    states <- t(vapply(
      inputs, duckvep_evidence_file_state,
      c(size = "", mtime = "", ctime = "")
    ))
    data.frame(
      id = names(inputs), path = unname(inputs),
      sha256 = unname(vapply(inputs, duckvep_evidence_sha256, character(1L))),
      states, row.names = NULL
    )
  }
  inputs_before <- snapshot_inputs()
  utils::write.table(inputs_before, file.path(directory, "inputs_before.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  perl_library <- paste(c(
    file.path(mirrors[c("vep", "variation")], "modules"),
    file.path(prefix, "share/ensembl-vep-116.0-0")
  ), collapse = .Platform$path.sep)
  summaries <- errors <- commands <- list()
  denominators <- data.frame(
    case = cases,
    input_records = NA_real_, input_alleles = NA_real_, prepared = FALSE
  )
  run <- function(case, stage, exe, args) {
    commands[[length(commands) + 1L]] <<- data.frame(case, stage,
      executable = exe,
      arguments = jsonlite::toJSON(unname(args), auto_unbox = FALSE)
    )
    value <- command(exe, args)
    writeLines(value, file.path(directory, case, paste0(stage, ".log")))
    invisible(TRUE)
  }
  attempt <- function(case, stage, expression) {
    tryCatch(force(expression), error = function(error) {
      errors[[length(errors) + 1L]] <<- data.frame(case, stage, message = conditionMessage(error))
      message(case, "/", stage, ": ", conditionMessage(error))
      FALSE
    })
  }
  csq_fields <- duckvep_fastvep_fields("vep_csq")
  tab_fields <- duckvep_fastvep_fields("native_tab17")
  for (case in cases) {
    case_directory <- file.path(directory, case)
    dir.create(case_directory)
    con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
    q <- function(x) as.character(DBI::dbQuoteString(con, x))
    execute <- function(sql) invisible(DBI::dbExecute(con, sql))
    prepared <- attempt(case, "prepare", {
      execute(paste("LOAD", q(extension)))
      execute("SET threads=1")
      gff <- duckvep_projection_fixture(con, root, inputs, case, case_directory)
      model <- file.path(case_directory, "model.duckdb")
      duckvep_fastvep_fixture_model(con, gff, inputs[["projection_reference"]], model)
      generated <- file.path(case_directory, "generated.vcf")
      compressed <- replay && grepl("\\.(gz|bgz)$", opt$replay_input)
      vcf <- file.path(case_directory, if (compressed) "input.vcf.gz" else "input.vcf")
      input_count <- NULL
      if (replay) {
        if (!file.copy(opt$replay_input, vcf)) stop("could not retain replay input")
        stopifnot(identical(duckvep_evidence_sha256(vcf), duckvep_evidence_sha256(opt$replay_input)))
      } else {
        run(case, "generate", "Rscript", c(
          file.path(root, "test/duckvep/conformance/generate_witnesses.R"),
          "--gff", gff, "--fasta", inputs[["projection_reference"]], "--ext", extension,
          "--out", generated, "--random-cases", opt$random_cases, "--seed", opt$seed
        ))
        input_count <- duckvep_projection_label_records(generated, vcf)
      }
      raw <- duckvep_fastvep_vcf_relation(con, vcf, duckvep_fastvep_vcf_header(vcf))
      input_count <- duckvep_fastvep_prepare_source(con, raw, input_count)
      execute(paste0("COPY field_input TO ", q(file.path(case_directory, "source_keys.parquet")), " (FORMAT PARQUET)"))
      stopifnot(system2("bgzip", c("-c", shQuote(gff)), stdout = paste0(gff, ".gz")) == 0L)
      run(case, "index_gff", "tabix", c("-p", "gff", paste0(gff, ".gz")))
      TRUE
    })
    if (!isTRUE(prepared)) {
      DBI::dbDisconnect(con, shutdown = TRUE)
      next
    }
    denominators[denominators$case == case, c("input_records", "input_alleles", "prepared")] <-
      list(input_count, input_count, TRUE)
    vep_home <- file.path(case_directory, "vep_home")
    dir.create(vep_home)
    vep_output <- file.path(case_directory, "vep.vcf")
    vep_ok <- attempt(case, "vep", run(case, "vep", "micromamba", c(
      "run", "--clean-env",
      "--env", paste0("HOME=", vep_home), "--env", paste0("PERL5LIB=", perl_library), "-p", prefix,
      "perl", file.path(mirrors[["vep"]], "vep"), "-i", vcf, "--fasta", inputs[["projection_reference"]],
      "--gff", paste0(gff, ".gz"), "--dir", vep_home, "--vcf", "--hgvs", "--numbers",
      "--symbol", "--biotype", "--canonical", "--mane", "--tsl", "--appris", "--ccds", "--protein",
      "--show_ref_allele", "--uploaded_allele", "--fields", paste(setdiff(csq_fields, "Uploaded_variation"), collapse = ","),
      "--distance", "5000", "--no_stats", "--force_overwrite", "-o", vep_output
    )))
    imported <- list()
    for (contract in c("native_tab17", "vep_csq")) {
      fields <- duckvep_fastvep_fields(contract)
      for (tool in c("duckvep", "fastvep")) {
        table <- paste(tool, contract, sep = "_")
        output <- file.path(case_directory, paste0(table, if (tool == "fastvep" && contract == "vep_csq") ".vcf" else ".tsv"))
        ok <- attempt(case, table, {
          if (tool == "duckvep") {
            args <- c(
              file.path(root, "benchmarks/benchmark_duckvep_fastvep_worker.R"),
              "--extension", extension, "--model", model, "--input", vcf, "--output", output,
              "--output-contract", contract, "--gff3", gff
            )
            if (contract == "vep_csq") args <- c(args, "--fasta", inputs[["projection_reference"]])
            run(case, table, "Rscript", args)
          } else {
            run(case, table, "env", c(
              "RAYON_NUM_THREADS=1", fastvep, "annotate", "--input", vcf,
              "--output", output, "--gff3", gff, "--fasta", inputs[["projection_reference"]],
              "--transcript-cache", file.path(case_directory, "fastvep.cache"),
              "--output-format", if (contract == "vep_csq") "vcf" else "tab",
              "--hgvs", "--distance", "5000", "--no-progress"
            ))
          }
          if (tool == "fastvep" && contract == "vep_csq") {
            duckvep_fastvep_extract_csq(con, output, paste0(table, "_raw"), fields)
          } else {
            duckvep_fastvep_read_field_tab(con, output, paste0(table, "_raw"), fields)
          }
          execute(paste0("CREATE TEMP TABLE ", table, " AS SELECT s.record_index, s.alt_index, a.*
            FROM ", table, "_raw a LEFT JOIN field_input s USING(Uploaded_variation)"))
          TRUE
        })
        imported[[table]] <- isTRUE(ok)
      }
    }
    oracle <- attempt(case, "vep_extract", {
      if (!isTRUE(vep_ok)) stop("VEP execution unavailable")
      duckvep_fastvep_extract_csq(con, vep_output, "vep_raw", csq_fields)
      execute("CREATE TEMP TABLE vep AS SELECT s.record_index, s.alt_index, a.*
        FROM vep_raw a LEFT JOIN field_input s USING(Uploaded_variation)")
      TRUE
    })
    comparisons <- list(
      native_tab17 = c("duckvep_native_tab17", "fastvep_native_tab17"),
      duckvep_vep_csq = c("duckvep_vep_csq", "vep"), fastvep_vep_csq = c("fastvep_vep_csq", "vep")
    )
    for (name in names(comparisons)) {
      attempt(case, paste0("compare_", name), {
        tables <- comparisons[[name]]
        ready <- vapply(tables, function(table) if (table == "vep") isTRUE(oracle) else isTRUE(imported[[table]]), logical(1L))
        if (!all(ready)) stop("comparison input unavailable: ", paste(tables[!ready], collapse = ", "))
        contract <- if (name == "native_tab17") name else "vep_csq"
        result <- duckvep_fastvep_compare(
          con, tables[[1L]], tables[[2L]], "field_input",
          duckvep_fastvep_fields(contract), file.path(case_directory, paste0("comparison_", name))
        )
        summaries[[length(summaries) + 1L]] <- cbind(case,
          comparison = name,
          input_records = input_count, input_alleles = input_count, result
        )
        print(summaries[[length(summaries)]])
        TRUE
      })
    }
    DBI::dbDisconnect(con, shutdown = TRUE)
  }
  inputs_unchanged <- attempt(NA_character_, "input_snapshot", {
    inputs_after <- snapshot_inputs()
    utils::write.table(inputs_after, file.path(directory, "inputs_after.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    if (!identical(inputs_before, inputs_after)) stop("shared staged inputs changed during conformance")
    TRUE
  })
  changed_sources <- sources[source_hashes != vapply(sources, duckvep_evidence_sha256, character(1L))]
  if (length(changed_sources)) {
    errors[[length(errors) + 1L]] <- data.frame(
      case = NA_character_,
      stage = "source_changed", message = paste(changed_sources, collapse = "; ")
    )
  }
  summary <- if (length(summaries)) do.call(rbind, summaries) else data.frame()
  errors <- if (length(errors)) do.call(rbind, errors) else data.frame(case = character(), stage = character(), message = character())
  utils::write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
  utils::write.csv(errors, file.path(directory, "errors.csv"), row.names = FALSE)
  utils::write.csv(denominators, file.path(directory, "cases.csv"), row.names = FALSE)
  if (length(commands)) utils::write.csv(do.call(rbind, commands), file.path(directory, "commands.csv"), row.names = FALSE)
  receipt <- c(
    source_revision = revision, build_binding = binding, extension_sha256 = extension_hash,
    fastvep_sha256 = fastvep_hash, fastvep_binding = fastvep_binding,
    pins, mode = if (replay) "replay" else "generated",
    seed = if (replay) NA_integer_ else opt$seed, random_cases = if (replay) NA_integer_ else opt$random_cases,
    requested_cases = length(cases), completed_cases = sum(denominators$prepared),
    input_records = sum(denominators$input_records), input_alleles = sum(denominators$input_alleles),
    completed_comparisons = nrow(summary), shared_inputs_unchanged = isTRUE(inputs_unchanged),
    errors = nrow(errors)
  )
  utils::write.table(data.frame(field = names(receipt), value = unname(receipt)),
    file.path(directory, "receipt.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  artifacts <- c(
    inputs, lock, extension, fastvep, sources,
    list.files(directory, recursive = TRUE, full.names = TRUE)
  )
  utils::write.table(data.frame(path = artifacts, sha256 = vapply(artifacts, duckvep_evidence_sha256, character(1L))),
    file.path(directory, "artifacts.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  check_pins()
  stopifnot(
    identical(extension_hash, duckvep_evidence_sha256(extension)),
    identical(fastvep_hash, duckvep_evidence_sha256(fastvep))
  )
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root, revision, allowed_outputs = opt$output)
  if (nrow(errors) || nrow(summary) != length(cases) * 3L || !all(summary$passed)) quit(status = 1L)
}

if (sys.nframe() == 0L) main()
