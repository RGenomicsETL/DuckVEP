#!/usr/bin/env Rscript
# Seeded original alleles on transcripts whose reference peptide carries an
# Ensembl translation SeqEdit (initial_met, _selenocysteine, _stop_codon_rt,
# amino_acid_sub). The pinned oracle applies the edits the way cache-mode VEP
# does: only the reference TranscriptVariationAllele peptide is edited, and only
# when the allele's translation window overlaps the edited residue.
source('scripts/duckvep_evidence.R')
source('test/duckvep/conformance/ambiguous_codon_differential.R')

edit_families <- c('stop_codon_rt_x', 'stop_codon_rt_w', 'selenocysteine',
  'initial_met', 'amino_acid_sub')

edit_standard_code <- function() {
  bases <- c('T', 'C', 'A', 'G')
  residues <- strsplit('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', '')[[1L]]
  codons <- character()
  for (first in bases) for (second in bases) for (third in bases)
    codons <- c(codons, paste0(first, second, third))
  setNames(residues, codons)
}

edit_random_bases <- function(n) paste(sample(c('A', 'C', 'G', 'T'), n, replace = TRUE), collapse = '')

# One model: a single-exon forward CDS with exactly one declared edit, flanked
# by ten A bases on each side so VCF anchors and parser spans stay inside the
# contig. The edited codon is the stratum under test; its neighbours and the
# surrounding indel window are the unedited controls on the same transcript.
edit_model <- function(index, family, code_table) {
  sense <- names(code_table)[code_table != '*']
  body <- sample(5:11, 1L)
  codons <- c('ATG', sample(sense, body, replace = TRUE), sample(c('TAA', 'TAG', 'TGA'), 1L))
  edited <- sample(seq.int(2L, body + 1L), 1L)
  edit <- switch(family,
    stop_codon_rt_x = {
      codons[edited] <- sample(c('TAA', 'TAG', 'TGA'), 1L)
      list(code = '_stop_codon_rt', position1 = edited, alternate = 'X')
    },
    stop_codon_rt_w = {
      codons[edited] <- 'TGA'
      list(code = '_stop_codon_rt', position1 = edited, alternate = 'W')
    },
    selenocysteine = {
      codons[edited] <- 'TGA'
      list(code = '_selenocysteine', position1 = edited, alternate = 'U')
    },
    initial_met = {
      edited <- 1L
      codons[1L] <- sample(c('CTG', 'TTG', 'GTG', 'ATT', 'ATC', 'ACG'), 1L)
      list(code = 'initial_met', position1 = 1L, alternate = 'M')
    },
    amino_acid_sub = {
      residue <- unname(code_table[[codons[edited]]])
      alternate <- sample(setdiff(strsplit('ACDEFGHIKLMNPQRSTVWY', '')[[1L]], residue), 1L)
      list(code = 'amino_acid_sub', position1 = edited, alternate = alternate)
    })
  cds <- paste(codons, collapse = '')
  variants <- list()
  seen <- character()
  add <- function(position, reference, alternate, variant_family) {
    key <- paste(position, reference, alternate)
    # A random replacement can spell the reference again; VEP has no ALT then.
    if (key %in% seen || identical(reference, alternate)) return(invisible(NULL))
    seen <<- c(seen, key)
    variants[[length(variants) + 1L]] <<- list(position1 = position, reference = reference,
      alternate = alternate, variant_family = variant_family)
    invisible(NULL)
  }
  # Every SNV in the edited codon and in the adjacent codons.
  for (codon in max(1L, edited - 1L):min(length(codons), edited + 1L)) {
    for (offset in 0:2) {
      position <- (codon - 1L) * 3L + offset + 1L
      reference <- substr(cds, position, position)
      for (alternate in setdiff(c('A', 'C', 'G', 'T'), reference))
        add(position, reference, alternate, 'snv')
    }
  }
  # Twelve seeded indels anchored in a window around the edited codon.
  first <- max(1L, (edited - 1L) * 3L + 1L - 3L)
  last <- min(nchar(cds) - 1L, edited * 3L + 2L)
  for (draw in seq_len(12L)) {
    position <- sample(first:last, 1L)
    anchor <- substr(cds, position, position)
    shape <- sample(c('ins', 'del', 'delins'), 1L)
    if (shape == 'ins') {
      add(position, anchor, paste0(anchor, edit_random_bases(sample(1:4, 1L))), 'ins')
    } else if (shape == 'del') {
      deleted <- sample(1:6, 1L)
      if (position + deleted <= nchar(cds))
        add(position, substr(cds, position, position + deleted), anchor, 'del')
    } else {
      replaced <- sample(2:3, 1L)
      if (position + replaced <= nchar(cds))
        add(position, substr(cds, position, position + replaced),
          paste0(anchor, edit_random_bases(sample(1:4, 1L))), 'delins')
    }
  }
  list(id = paste0('model', index, '_', family), family = family, cds = cds, table = 1L,
    edit = edit, variants = variants)
}

edit_inputs <- function(models, seed) {
  set.seed(seed)
  code_table <- edit_standard_code()
  cases <- vector('list', models)
  events <- vector('list', models)
  event_count <- 0L
  for (i in seq_len(models)) {
    family <- edit_families[(i - 1L) %% length(edit_families) + 1L]
    model <- edit_model(i, family, code_table)
    ids <- event_count + seq_along(model$variants)
    variants <- lapply(seq_along(model$variants), function(j) {
      v <- model$variants[[j]]
      list(id = as.character(ids[[j]]), position1 = v$position1,
        reference = v$reference, alternate = v$alternate)
    })
    cases[[i]] <- list(id = model$id, cds = model$cds, table = model$table,
      edits = list(model$edit), variants = variants, variant_format = 'vcf',
      genomic_sequence = paste0(strrep('A', 10L), model$cds, strrep('A', 10L)),
      cds_start1 = 11L)
    edit_start <- (model$edit$position1 - 1L) * 3L + 1L
    events[[i]] <- data.frame(event_index = ids, seq_region = ids - 1L,
      transcript_index = ids - 1L,
      position = vapply(model$variants, `[[`, 0L, 'position1') + 10L,
      reference = vapply(model$variants, `[[`, '', 'reference'),
      alternate = vapply(model$variants, `[[`, '', 'alternate'),
      cds = model$cds, table = model$table, case_id = model$id, family = family,
      variant_family = vapply(model$variants, `[[`, '', 'variant_family'),
      edit_code = model$edit$code, edit_position = model$edit$position1,
      edit_amino = model$edit$alternate,
      # Bases between the uploaded REF span and the edited codon; zero when
      # they overlap. Insertions use the anchor base as their span.
      edit_distance = vapply(model$variants, function(v) {
        ref_start <- v$position1
        ref_end <- v$position1 + nchar(v$reference) - 1L
        max(0L, edit_start - ref_end, ref_start - (edit_start + 2L))
      }, 0L), stringsAsFactors = FALSE)
    event_count <- event_count + length(model$variants)
  }
  events <- do.call(rbind, events)
  stopifnot(!anyDuplicated(events$event_index), identical(events$event_index, seq_len(nrow(events))))
  list(events = events, cases = cases)
}

# One native transcript per original event, each carrying the model's edit.
edit_native <- function(con, events, expected, fasta, out) {
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbWriteTable(con, 'inputs', events)
  tx <- paste('SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,',
    '11::UBIGINT transcript_start,10+length(cds)::UBIGINT transcript_end,1::TINYINT strand,',
    'transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,',
    'transcript_start cds_start,transcript_end cds_end,cds::BLOB cds_sequence,',
    "\"table\"::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence",
    'FROM inputs ORDER BY transcript_index')
  ex <- paste('SELECT transcript_index::UINTEGER transcript_index,11::UBIGINT exon_start,',
    '10+length(cds)::UBIGINT exon_end,1::UBIGINT exon_cdna_start,length(cds)::UBIGINT exon_cdna_end,',
    '0::TINYINT phase,0::TINYINT end_phase FROM inputs ORDER BY transcript_index')
  regions <- paste('SELECT seq_region::UINTEGER seq_region,(length(cds)+20)::UBIGINT sequence_length,',
    'event_index::VARCHAR seq_region_name FROM inputs ORDER BY seq_region')
  edits <- paste('SELECT transcript_index::UINTEGER transcript_index,',
    'edit_position::UINTEGER protein_position,edit_amino::VARCHAR alternate_amino_acid,',
    'edit_code::VARCHAR edit_code FROM inputs ORDER BY transcript_index')
  stopifnot(DBI::dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('edits',",
    q(regions), ',', q(tx), ',', q(ex), ',reference_fasta:=', q(fasta),
    ',peptide_edit_query:=', q(edits), ')'))$loaded)
  DBI::dbExecute(con, paste('CREATE TABLE events AS SELECT event_index,seq_region,position,reference,alternate,',
    'NULL::UBIGINT end_position,NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,',
    'NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position FROM inputs ORDER BY seq_region,position'))
  comparisons <- list()
  for (threads in c(1L, 4L)) {
    DBI::dbExecute(con, paste('SET threads=', threads))
    label <- paste0('independent_', threads)
    message('Native route: ', label)
    DBI::dbExecute(con, paste0('CREATE TABLE ', label, " AS SELECT * FROM query(duckvep_annotate_sql('events',",
      "'edits', struct_pack(hgvs:=true,rich:=true,upstream_distance:=0,downstream_distance:=0)))"))
    DBI::dbExecute(con, paste('COPY', label, 'TO', q(file.path(out, paste0(label, '.parquet'))),
      '(FORMAT PARQUET)'))
    geometry <- DBI::dbGetQuery(con, paste('SELECT event_index,transcript_index FROM', label))
    stopifnot(!anyNA(geometry$event_index), !anyDuplicated(geometry$event_index),
      all(geometry$transcript_index == geometry$event_index - 1L))
    actual <- DBI::dbGetQuery(con, paste('SELECT a.event_index::INTEGER event_index,i.alternate allele,',
      'a.protein_hgvs hgvsp,(SELECT string_agg(t.consequence,\'&\' ORDER BY t.consequence)',
      'FROM duckvep_so_terms() t WHERE (a.consequence_mask & t.consequence_mask)<>0) so',
      'FROM', label, 'a LEFT JOIN inputs i USING(event_index)'))
    pairs <- codon_equal(actual, expected)
    facts <- DBI::dbGetQuery(con, paste('SELECT event_index::INTEGER event_index,duckvep_status,',
      'duckvep_reason,reference_amino_acid,alternate_amino_acid,protein_position FROM', label))
    comparisons[[label]] <- merge(pairs, facts, by = 'event_index', all.x = TRUE, sort = FALSE)
    DBI::dbRemoveTable(con, label)
  }
  stopifnot(DBI::dbGetQuery(con, "SELECT duckvep_model_drop('edits') dropped")$dropped)
  DBI::dbRemoveTable(con, 'events')
  DBI::dbRemoveTable(con, 'inputs')
  comparisons
}

edit_main <- function() {
  suppressPackageStartupMessages(library(DBI))
  opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option('--extension', default = 'build/release/duckvep.duckdb_extension'),
    optparse::make_option('--models', type = 'integer', default = 300L,
      help = 'seeded single-edit models; each carries about 40 original alleles'),
    optparse::make_option('--seed', type = 'integer', default = 173L),
    optparse::make_option('--evidence-out', dest = 'evidence_out', default = '',
      help = 'Retain a compact complete comparison bundle in a new directory'),
    optparse::make_option('--vep-prefix', dest = 'vep_prefix',
      default = Sys.getenv('VEP_PREFIX'))
  )))
  stopifnot(opt$models >= length(edit_families), !is.na(opt$seed))
  revision <- duckvep_evidence_revision('.')
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  sources <- c('test/duckvep/conformance/peptide_edit_differential.R',
    'test/duckvep/conformance/ambiguous_codon_differential.R',
    'test/duckvep/conformance/reference_translation_oracle.pl', 'scripts/duckvep_evidence.R',
    list.files('src/kernel', recursive = TRUE, full.names = TRUE, pattern = '\\.[ch]$'))
  source_hashes <- vapply(sources, duckvep_evidence_sha256, '')
  extension_hash <- duckvep_evidence_sha256(extension)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082',
    variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- setNames(normalizePath(paste0('.sync/ensembl-', names(pins))), names(pins))
  for (name in names(pins)) {
    stopifnot(identical(duckvep_evidence_command('git',
      c('-C', mirrors[[name]], 'rev-parse', 'HEAD'), 'oracle revision'), unname(pins[name])),
      !length(duckvep_evidence_command('git', c('-C', mirrors[[name]], 'status', '--porcelain'),
        'oracle checkout')))
  }
  if (!nzchar(opt$vep_prefix)) stop('--vep-prefix or VEP_PREFIX is required', call. = FALSE)
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  out <- tempfile(paste0('peptide_edit_seed', opt$seed, '_'), 'test/duckvep/conformance/results')
  stopifnot(dir.create(out))
  out <- normalizePath(out)
  message('Peptide-edit artifacts: ', out)
  command <- function(args, label) {
    status <- system2('micromamba', shQuote(args),
      stdout = file.path(out, paste0(label, '.stdout')),
      stderr = file.path(out, paste0(label, '.stderr')))
    stopifnot(status == 0L)
  }
  command(c('list', '--explicit', '-p', prefix), 'environment')
  stopifnot(identical(duckvep_evidence_explicit_packages(readLines(file.path(out, 'environment.stdout'))),
    duckvep_evidence_explicit_packages(readLines(
      'test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  modules <- c(file.path(prefix, 'share/ensembl-vep-116.0-0/Bio/EnsEMBL',
    c('Transcript.pm', 'Translation.pm', 'SeqEdit.pm', 'Slice.pm')),
    file.path(prefix, 'share/ensembl-vep-116.0-0/modules/Bio/EnsEMBL/VEP/AnnotationSource/Database/Transcript.pm'),
    file.path(mirrors[['vep']], 'modules/Bio/EnsEMBL/VEP', c('Config.pm', 'Parser.pm', 'Parser/VCF.pm')),
    file.path(mirrors[['variation']], 'modules/Bio/EnsEMBL/Variation',
      c('BaseTranscriptVariation.pm', 'TranscriptVariation.pm', 'TranscriptVariationAllele.pm',
        'VariationFeatureOverlapAllele.pm', 'Utils/VariationEffect.pm', 'Utils/FastaSequence.pm')))
  module_hashes <- vapply(modules, duckvep_evidence_sha256, '')
  inputs <- edit_inputs(opt$models, opt$seed)
  events <- inputs$events
  cases <- inputs$cases
  # Every event owns a native transcript; the model stores region ordinals in 16 bits.
  stopifnot(nrow(events) <= 56000L)
  saveRDS(events, file.path(out, 'events.rds'))
  input <- file.path(out, 'cases.jsonl')
  writeLines(vapply(cases, jsonlite::toJSON, '', auto_unbox = TRUE), input)
  libs <- paste(c(file.path(mirrors, 'modules'),
    file.path(prefix, 'share/ensembl-vep-116.0-0')), collapse = ':')
  message('VEP oracle: ', length(cases), ' models / ', nrow(events), ' original alleles')
  command(c('run', '--clean-env', '--env', paste0('PERL5LIB=', libs), '-p', prefix, 'perl',
    normalizePath('test/duckvep/conformance/reference_translation_oracle.pl'), input), 'oracle')
  oracle <- lapply(readLines(file.path(out, 'oracle.stdout')), jsonlite::fromJSON,
    simplifyVector = FALSE)
  stopifnot(identical(vapply(oracle, `[[`, '', 'id'), vapply(cases, `[[`, '', 'id')))
  # The oracle's reference protein must carry every declared edit.
  for (i in seq_along(cases)) {
    edit <- cases[[i]]$edits[[1L]]
    stopifnot(identical(substr(oracle[[i]]$reference, edit$position1, edit$position1), edit$alternate))
  }
  expected <- codon_expected(oracle, events, TRUE)
  stopifnot(setequal(expected$event_index, events$event_index))
  controls <- codon_controls(expected)
  write.csv(expected, file.path(out, 'expected.csv'), row.names = FALSE)
  write.csv(controls, file.path(out, 'controls.csv'), row.names = FALSE)
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(extension)))
  fasta <- file.path(out, 'reference.fa')
  writeLines(as.vector(rbind(paste0('>', events$event_index),
    paste0(strrep('A', 10L), events$cds, strrep('A', 10L)))), fasta)
  stopifnot(system2('samtools', c('faidx', shQuote(fasta)),
    stdout = file.path(out, 'faidx.stdout'), stderr = file.path(out, 'faidx.stderr')) == 0L)
  comparisons <- edit_native(con, events, expected, fasta, out)
  pairs <- do.call(rbind, lapply(names(comparisons), function(route)
    cbind(route = route, comparisons[[route]])))
  pairs <- merge(pairs, events[c('event_index', 'case_id', 'family', 'variant_family',
    'edit_code', 'edit_position', 'edit_amino', 'edit_distance', 'position', 'reference', 'cds')],
    by = 'event_index', all.x = TRUE, sort = FALSE)
  pairs$edit_overlap <- ifelse(is.na(pairs$edit_distance), 'unmatched',
    ifelse(pairs$edit_distance == 0L, 'edited_codon', 'control'))
  write.csv(pairs, file.path(out, 'pairs.csv'), row.names = FALSE)
  summary <- aggregate(cbind(pairs = rep(1L, nrow(pairs)), hgvsp_failures = !pairs$hgvsp_equal,
    so_failures = !pairs$so_equal), pairs[c('route', 'family', 'variant_family', 'edit_overlap')],
    sum, na.rm = TRUE)
  write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
  files <- list.files(out, full.names = TRUE)
  manifest <- list(source_revision = revision, source_binding = 'diagnostic_unbound',
    extension_path = extension, extension_sha256 = extension_hash, oracle_revisions = pins,
    seed = opt$seed, models = length(cases), source_alleles = nrow(events), oracle_pairs = nrow(expected),
    scope = 'independent_SO_and_HGVSp_original_VCF_alleles_on_single_SeqEdit_forward_single_exon_transcripts',
    sha256 = as.list(c(source_hashes, module_hashes, vapply(files, duckvep_evidence_sha256, ''))))
  jsonlite::write_json(manifest, file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  print(summary)
  stopifnot(identical(source_hashes, vapply(sources, duckvep_evidence_sha256, '')),
    identical(module_hashes, vapply(modules, duckvep_evidence_sha256, '')),
    identical(extension_hash, duckvep_evidence_sha256(extension)),
    identical(revision, duckvep_evidence_revision('.')))
  if (nzchar(opt$evidence_out)) {
    destination <- opt$evidence_out
    stopifnot(!dir.exists(destination), dir.create(destination))
    dbWriteTable(con, 'comparison_pairs', pairs)
    dbExecute(con, paste('COPY comparison_pairs TO', q(file.path(destination, 'pairs.parquet')),
      '(FORMAT PARQUET)'))
    for (name in c('cases.jsonl', 'oracle.stdout', 'oracle.stderr')) {
      compressed <- gzfile(file.path(destination, paste0(name, '.gz')), 'wt')
      writeLines(readLines(file.path(out, name)), compressed)
      close(compressed)
    }
    stopifnot(all(file.copy(file.path(out, c('summary.csv', 'controls.csv', 'environment.stdout')),
      destination)))
    manifest$local_receipt_sha256 <- duckvep_evidence_sha256(file.path(out, 'receipt.json'))
    manifest$source_sha256 <- as.list(c(source_hashes, module_hashes))
    published <- list.files(destination, full.names = TRUE)
    manifest$sha256 <- as.list(setNames(vapply(published, duckvep_evidence_sha256, ''), basename(published)))
    jsonlite::write_json(manifest, file.path(destination, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  }
  stopifnot(all(controls$rejected), all(pairs$hgvsp_equal), all(pairs$so_equal))
}
if (sys.nframe() == 0L) edit_main()
