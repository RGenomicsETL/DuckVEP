#!/usr/bin/env Rscript
# Seeded paired-allele expansion and evidence-status checks against base R.
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(optparse) })
op <- OptionParser()
op <- add_option(op, '--extension', default = 'build/release/duckhts.duckdb_extension')
op <- add_option(op, '--trials', type = 'integer', default = 100000L)
op <- add_option(op, '--seed', type = 'integer', default = 173L)
op <- add_option(op, '--out', default = '')
opt <- parse_args(op)
stopifnot(!is.na(opt$trials), opt$trials >= 49L, !is.na(opt$seed))

same_nullable <- function(actual, expected) {
  same <- (is.na(actual) & is.na(expected)) |
    (!is.na(actual) & !is.na(expected) & actual == expected)
  all(same)
}

check_comparisons <- function(x, ids) {
  fields <- c('reference', 'alternate', 'reference_length', 'alternate_length',
              'length_change', 'length_direction', 'status')
  required_columns <- c('scene', paste0('expected_', fields), paste0('actual_', fields))
  stopifnot(all(required_columns %in% names(x)))
  stopifnot(nrow(x) == length(ids), !anyDuplicated(x$scene),
            identical(as.integer(x$scene), as.integer(ids)))
  for (field in fields)
    stopifnot(same_nullable(x[[paste0('actual_', field)]],
                            x[[paste0('expected_', field)]]))
  TRUE
}

empty_components <- function() {
  data.frame(scene = integer(), allele = character(), ordinal = integer(),
             unit = character(), unit_count = double())
}

run <- function() {
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  out <- if (nzchar(opt$out)) opt$out else tempfile(paste0('repeat_alleles_seed', opt$seed, '_'),
    'test/duckvep/conformance/results')
  stopifnot(!dir.exists(out), dir.create(out, recursive = TRUE))
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(extension)))
  dbExecute(con, 'SET threads=1')
  set.seed(opt$seed)

  alphabet <- strsplit('ACGTRYSWKMBDHVNacgtryswkmbdhvn', '', fixed = TRUE)[[1L]]
  states <- c('exact', 'missing_count', 'missing_unit', 'fractional', 'empty',
              'missing_components', 'summary')
  cells <- as.vector(outer(states, states, paste, sep = '__'))
  coverage <- setNames(integer(length(cells)), cells)
  status_coverage <- setNames(integer(4L),
                              c('ok', 'summary_only', 'incomplete_input', 'nonintegral_count'))
  direction_coverage <- setNames(integer(3L), c('GAIN', 'LOSS', 'NEUTRAL'))
  precedence_coverage <- setNames(integer(4L), c('summary_over_incomplete',
    'summary_over_nonintegral', 'incomplete_over_nonintegral_ref',
    'incomplete_over_nonintegral_alt'))
  controls <- NULL
  capacity_checks <- list()

  make_axis <- function(scene, allele, state) {
    if (state %in% c('empty', 'missing_components'))
      return(list(has_components = state == 'empty', rows = empty_components(), sequence = ''))
    n <- sample.int(6L, 1L)
    units <- vapply(seq_len(n), function(j)
      paste0(sample(alphabet, sample.int(9L, 1L), replace = TRUE), collapse = ''), '')
    counts <- as.double(sample(c(0:3, 10, 31, 100), n, replace = TRUE))
    counts[1L] <- sample(c(1, 2, 3, 10, 31, 100), 1L)
    if (state == 'missing_count')
      counts[sample.int(n, 1L)] <- NA_real_
    if (state == 'missing_unit')
      units[sample.int(n, 1L)] <- NA_character_
    if (state == 'fractional') {
      at <- sample.int(n, 1L)
      counts[at] <- counts[at] + 0.5
    }
    if (state == 'summary')
      counts[sample.int(n, 1L)] <- 1e300
    rows <- data.frame(scene = rep(scene, n), allele = rep(allele, n), ordinal = seq_len(n),
                       unit = units, unit_count = counts)
    sequence <- if (state == 'exact') paste0(strrep(units, as.integer(counts)), collapse = '')
      else NA_character_
    list(has_components = TRUE, rows = rows, sequence = sequence)
  }

  expected_status <- function(reference_state, alternate_state, sequence_exact) {
    if (!sequence_exact)
      return('summary_only')
    incomplete <- c('missing_count', 'missing_unit', 'missing_components')
    if (reference_state %in% incomplete || alternate_state %in% incomplete)
      return('incomplete_input')
    if (reference_state == 'fractional' || alternate_state == 'fractional')
      return('nonintegral_count')
    'ok'
  }

  for (first in seq.int(1L, opt$trials, by = 1000L)) {
    ids <- seq.int(first, min(first + 999L, opt$trials))
    cases <- data.frame(scene = ids, reference_state = NA_character_,
      alternate_state = NA_character_, sequence_exact = TRUE, capacity = 5000,
      has_reference_components = TRUE, has_alternate_components = TRUE,
      expected_reference = NA_character_, expected_alternate = NA_character_,
      expected_reference_length = NA_real_, expected_alternate_length = NA_real_,
      expected_length_change = NA_real_, expected_length_direction = NA_character_,
      expected_status = NA_character_)
    components <- vector('list', length(ids) * 2L)
    for (i in seq_along(ids)) {
      cell <- (ids[i] - 1L) %% length(cells)
      reference_state <- states[cell %/% length(states) + 1L]
      alternate_state <- states[cell %% length(states) + 1L]
      cell_name <- paste(reference_state, alternate_state, sep = '__')
      coverage[cell_name] <- coverage[cell_name] + 1L
      sequence_exact <- reference_state != 'summary' && alternate_state != 'summary'
      reference <- make_axis(ids[i], 'reference', reference_state)
      alternate <- make_axis(ids[i], 'alternate', alternate_state)
      status <- expected_status(reference_state, alternate_state, sequence_exact)
      status_coverage[status] <- status_coverage[status] + 1L

      cases$reference_state[i] <- reference_state
      cases$alternate_state[i] <- alternate_state
      cases$sequence_exact[i] <- sequence_exact
      cases$has_reference_components[i] <- reference$has_components
      cases$has_alternate_components[i] <- alternate$has_components
      cases$expected_status[i] <- status
      if (status == 'ok') {
        cases$expected_reference[i] <- reference$sequence
        cases$expected_alternate[i] <- alternate$sequence
        cases$expected_reference_length[i] <- nchar(reference$sequence, type = 'bytes')
        cases$expected_alternate_length[i] <- nchar(alternate$sequence, type = 'bytes')
        cases$expected_length_change[i] <- cases$expected_alternate_length[i] -
          cases$expected_reference_length[i]
        cases$expected_length_direction[i] <- if (cases$expected_length_change[i] > 0) 'GAIN'
          else if (cases$expected_length_change[i] < 0) 'LOSS' else 'NEUTRAL'
        cases$capacity[i] <- max(cases$expected_reference_length[i],
                                 cases$expected_alternate_length[i]) + sample(0:2, 1L)
        direction_coverage[cases$expected_length_direction[i]] <-
          direction_coverage[cases$expected_length_direction[i]] + 1L
      } else if (!sequence_exact) {
        cases$capacity[i] <- 0
      }
      if (!sequence_exact && (reference_state %in% c('missing_count', 'missing_unit',
          'missing_components') || alternate_state %in% c('missing_count', 'missing_unit',
          'missing_components')))
        precedence_coverage['summary_over_incomplete'] <-
          precedence_coverage['summary_over_incomplete'] + 1L
      if (!sequence_exact && (reference_state == 'fractional' || alternate_state == 'fractional'))
        precedence_coverage['summary_over_nonintegral'] <-
          precedence_coverage['summary_over_nonintegral'] + 1L
      if (sequence_exact && reference_state %in% c('missing_count', 'missing_unit',
          'missing_components') && alternate_state == 'fractional')
        precedence_coverage['incomplete_over_nonintegral_ref'] <-
          precedence_coverage['incomplete_over_nonintegral_ref'] + 1L
      if (sequence_exact && alternate_state %in% c('missing_count', 'missing_unit',
          'missing_components') && reference_state == 'fractional')
        precedence_coverage['incomplete_over_nonintegral_alt'] <-
          precedence_coverage['incomplete_over_nonintegral_alt'] + 1L
      components[[2L * i - 1L]] <- reference$rows
      components[[2L * i]] <- alternate$rows
    }
    components <- do.call(rbind, components)
    dbWriteTable(con, 'repeat_cases', cases, overwrite = TRUE)
    dbWriteTable(con, 'repeat_components', components, overwrite = TRUE)
    query <- paste(
      'WITH refs AS (SELECT scene,list(struct_pack(unit:=unit,count:=unit_count)',
      "ORDER BY ordinal) parts FROM repeat_components WHERE allele='reference' GROUP BY scene),",
      'alts AS (SELECT scene,list(struct_pack(unit:=unit,count:=unit_count)',
      "ORDER BY ordinal) parts FROM repeat_components WHERE allele='alternate' GROUP BY scene),",
      'evaluated AS (SELECT c.*,duckvep_repeat_alleles(',
      'CASE WHEN c.has_reference_components THEN coalesce(r.parts,',
      '[]::STRUCT(unit VARCHAR,count DOUBLE)[]) ELSE NULL END,',
      'CASE WHEN c.has_alternate_components THEN coalesce(a.parts,',
      '[]::STRUCT(unit VARCHAR,count DOUBLE)[]) ELSE NULL END,',
      'c.sequence_exact,max_allele_bases:=c.capacity) result',
      'FROM repeat_cases c LEFT JOIN refs r USING(scene) LEFT JOIN alts a USING(scene))',
      'SELECT scene,expected_reference,expected_alternate,expected_reference_length,',
      'expected_alternate_length,expected_length_change,expected_length_direction,expected_status,',
      'result.reference AS actual_reference,result.alternate AS actual_alternate,',
      'result.reference_length AS actual_reference_length,',
      'result.alternate_length AS actual_alternate_length,',
      'result.length_change AS actual_length_change,',
      'result.length_direction AS actual_length_direction,',
      'result.status AS actual_status FROM evaluated ORDER BY scene')
    actual <- try(dbGetQuery(con, query), silent = TRUE)
    if (inherits(actual, 'try-error')) {
      saveRDS(list(seed = opt$seed, first = first, cases = cases, components = components,
        error = as.character(actual)), file.path(out, 'counterexample.rds'))
      stop(actual)
    }
    dbWriteTable(con, 'cases_all', cases, append = first != 1L)
    dbWriteTable(con, 'components_all', components, append = first != 1L)
    dbWriteTable(con, 'comparisons', actual, append = first != 1L)
    ok <- try(check_comparisons(actual, ids), silent = TRUE)
    if (inherits(ok, 'try-error')) {
      saveRDS(list(seed = opt$seed, first = first, cases = cases, components = components,
        comparisons = actual, error = as.character(ok)), file.path(out, 'counterexample.rds'))
      stop(ok)
    }

    if (is.null(controls)) {
      exact_row <- which(actual$expected_status == 'ok')[1L]
      duplicate <- actual
      duplicate$scene[2L] <- duplicate$scene[1L]
      reference <- actual
      reference$actual_reference[exact_row] <- paste0(reference$actual_reference[exact_row], 'A')
      alternate <- actual
      alternate$actual_alternate[exact_row] <- paste0(alternate$actual_alternate[exact_row], 'C')
      length <- actual
      length$actual_reference_length[exact_row] <- length$actual_reference_length[exact_row] + 1
      direction <- actual
      direction$actual_length_direction[exact_row] <-
        if (direction$actual_length_direction[exact_row] == 'GAIN') 'LOSS' else 'GAIN'
      status <- actual
      status$actual_status[exact_row] <- 'summary_only'
      controls <- vapply(list(drop = actual[-1L, ], duplicate = duplicate,
        reference = reference, alternate = alternate, length = length, direction = direction,
        status = status), function(x)
        inherits(try(check_comparisons(x, ids), silent = TRUE), 'try-error'), TRUE)
      stopifnot(all(controls))

      capacity_candidates <- list(
        reference = which(cases$expected_status == 'ok' &
          cases$expected_reference_length > cases$expected_alternate_length),
        alternate = which(cases$expected_status == 'ok' &
          cases$expected_alternate_length > cases$expected_reference_length))
      for (axis in names(capacity_candidates)) {
        stopifnot(length(capacity_candidates[[axis]]) > 0L)
        for (i in head(capacity_candidates[[axis]], 16L)) {
          required <- cases[[paste0('expected_', axis, '_length')]][i]
          limit <- required - 1
          dbExecute(con, sprintf('UPDATE repeat_cases SET capacity=%d WHERE scene=%d',
                                 limit, ids[i]))
          failure <- try(dbGetQuery(con, query), silent = TRUE)
          dbExecute(con, sprintf('UPDATE repeat_cases SET capacity=%d WHERE scene=%d',
                                 cases$capacity[i], ids[i]))
          passed <- inherits(failure, 'try-error') &&
            grepl(paste0('duckvep_repeat_alleles: ', axis, ' requires '),
                  as.character(failure), fixed = TRUE) &&
            grepl('exceeds max_allele_bases=', as.character(failure), fixed = TRUE)
          capacity_checks[[length(capacity_checks) + 1L]] <- data.frame(axis = axis,
            scene = ids[i], capacity = limit, required = required, passed = passed)
          if (!passed) {
            saveRDS(list(cases = cases, components = components, axis = axis, scene = ids[i],
              capacity = limit, result = failure),
              file.path(out, paste0(axis, '_capacity_counterexample.rds')))
            stop(axis, ' capacity exhaustion did not report the named limit')
          }
        }
      }
      recovered <- dbGetQuery(con, query)
      recovery_ok <- try(check_comparisons(recovered, ids), silent = TRUE)
      if (inherits(recovery_ok, 'try-error') || !identical(recovered, actual)) {
        saveRDS(list(cases = cases, components = components, before = actual, after = recovered),
                file.path(out, 'capacity_recovery_counterexample.rds'))
        stop('capacity recovery was not exact')
      }
    }
  }

  stopifnot(sum(coverage) == opt$trials, all(coverage >= opt$trials %/% length(cells)),
    all(status_coverage > 0L), all(direction_coverage > 0L),
    all(precedence_coverage > 0L),
    dbGetQuery(con, 'SELECT count(*) n FROM comparisons')$n == opt$trials,
    dbGetQuery(con, 'SELECT count(DISTINCT scene) n FROM comparisons')$n == opt$trials)
  capacity_checks <- do.call(rbind, capacity_checks)
  stopifnot(all(capacity_checks$passed),
            identical(sort(unique(capacity_checks$axis)), c('alternate', 'reference')))
  for (table in c('cases_all', 'components_all', 'comparisons'))
    dbExecute(con, paste('COPY', table, 'TO', q(file.path(out, paste0(table, '.parquet'))),
                        '(FORMAT PARQUET)'))
  write.csv(capacity_checks, file.path(out, 'capacity_checks.csv'), row.names = FALSE)
  sha <- function(path) digest::digest(file = path, algo = 'sha256', serialize = FALSE)
  sources <- c('test/duckvep/conformance/repeat_alleles_differential.R',
               'src/duckvep/duckvep_sql.c')
  for (path in sources) {
    destination <- file.path(out, 'source', path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(path, destination), sha(path) == sha(destination))
  }
  files <- list.files(out, recursive = TRUE, full.names = TRUE)
  jsonlite::write_json(list(scope = paste('paired ordered exact repeat alleles, signed length',
    'change, length direction and explicit status precedence'), seed = opt$seed,
    trials = opt$trials, threads = 1L, batch_size = 1000L,
    paired_strata = as.list(coverage), statuses = as.list(status_coverage),
    length_directions = as.list(direction_coverage), precedence = as.list(precedence_coverage),
    corruption_controls = as.list(controls),
    capacity_exhaustion_checks = as.list(table(capacity_checks$axis)), recovery_exact = TRUE,
    oracle = 'base R strrep and byte lengths over independently generated ordered alleles',
    R_version = R.version.string, RNG_kind = RNGkind(),
    duckdb_version = as.character(packageVersion('duckdb')),
    source_revision = system2('git', c('rev-parse', 'HEAD'), stdout = TRUE),
    source_dirty = length(system2('git', c('status', '--porcelain'), stdout = TRUE)) > 0L,
    build_binding = 'diagnostic_unbound', extension_sha256 = sha(extension),
    failed = 0L, skipped = 0L, failures_waived = 0L,
    interpretation = paste('Generated paired-component scenes and status strata, not',
      'independent biological samples or a VEP parser-conformance claim.'),
    sha256 = as.list(setNames(vapply(files, sha, ''), substring(files, nchar(out) + 2L)))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  cat(opt$trials, 'scenes passed; all 49 paired strata, seven corruption controls, and',
      'both capacity axes recovered exactly;', out, '\n')
}
run()
