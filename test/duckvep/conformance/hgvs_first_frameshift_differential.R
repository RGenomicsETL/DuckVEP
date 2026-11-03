#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 3L, args[3L] %in% c('baseline', 'fixed'))
  extension <- normalizePath(args[1L], mustWork = TRUE)
  out <- args[2L]
  stopifnot(!dir.exists(out), dir.create(out, recursive = TRUE))
  out <- normalizePath(out)
  prefix <- Sys.getenv('DUCKVEP_VEP_PREFIX')
  if (!nzchar(prefix)) {
    stop('Set DUCKVEP_VEP_PREFIX to the VEP 116 Conda environment')
  }
  vep <- file.path(prefix, 'share/ensembl-vep-116.0-0/vep')
  rc <- function(x) {
    paste0(
      rev(strsplit(chartr('ACGT', 'TGCA', x), '', fixed = TRUE)[[1L]]),
      collapse = ''
    )
  }
  codons <- c(
    AAA = 'K',
    CCT = 'P',
    CAA = 'Q',
    TAA = '*',
    TCA = 'S',
    ACA = 'T',
    ACC = 'T',
    ATA = 'I',
    GCC = 'A',
    AAT = 'N',
    CAC = 'H'
  )
  translate_fixture <- function(x) {
    starts <- seq.int(1L, nchar(x) - 2L, by = 3L)
    peptide <- codons[substring(x, starts, starts + 2L)]
    stopifnot(!anyNA(peptide))
    paste0(peptide, collapse = '')
  }
  cmd <- function(exe, args, label) {
    status <- system2(
      exe,
      shQuote(args),
      stdout = file.path(out, paste0(label, '.stdout')),
      stderr = file.path(out, paste0(label, '.stderr'))
    )
    stopifnot(status == 0L)
  }
  models <- expand.grid(
    strand = c(1L, -1L),
    site = c('first', 'later', 'shifted'),
    stringsAsFactors = FALSE
  )
  models$transcript_index <- models$seq_region <- seq_len(nrow(models)) - 1L
  models$chrom <- paste0('chr', seq_len(nrow(models)))
  models$cds <- ifelse(models$site == 'shifted', 'AAATCACAATAA', 'AAACCTCAATAA')
  models$post_cds <- ifelse(models$site == 'shifted', 'ATAA', 'AATAA')
  models$tx_start <- 1001L
  models$tx_end <- models$tx_start +
    nchar(models$cds) +
    nchar(models$post_cds) -
    1L
  models$sequence_length <- 2000L + nchar(models$cds) + nchar(models$post_cds)
  models$cds_start <- ifelse(
    models$strand == 1L,
    models$tx_start,
    models$tx_start + nchar(models$post_cds)
  )
  models$cds_end <- models$cds_start + 11L
  fasta <- gff <- records <- vector('list', nrow(models))
  for (i in seq_len(nrow(models))) {
    m <- models[i, ]
    transcript <- paste0(m$cds, m$post_cds)
    genome <- paste0(
      strrep('GTCGATCGTA', 100L),
      if (m$strand == 1L) transcript else rc(transcript),
      strrep('CGTAGCTAGC', 100L)
    )
    fasta[[i]] <- c(paste0('>', m$chrom), genome)
    id <- paste0('TX', i)
    attrs <- c(
      paste0('ID=gene:', id, ';biotype=protein_coding'),
      paste0(
        'ID=transcript:',
        id,
        ';Parent=gene:',
        id,
        ';biotype=protein_coding;tag=cds_start_NF'
      ),
      paste0('ID=exon:', id, ';Parent=transcript:', id),
      paste0('Parent=transcript:', id)
    )
    gff[[i]] <- paste(
      m$chrom,
      'fixture',
      c('gene', 'mRNA', 'exon', 'CDS'),
      c(rep(m$tx_start, 3L), m$cds_start),
      c(rep(m$tx_end, 3L), m$cds_end),
      '.',
      if (m$strand == 1L) '+' else '-',
      c('.', '.', '.', '0'),
      attrs,
      sep = '\t'
    )
    anchor <- if (m$site == 'later') 3L else 1L
    if (m$site == 'shifted') {
      ref <- substr(m$cds, anchor, anchor + 1L)
      alt <- substr(ref, 1L, 1L)
      expected_cds <- paste0(
        substr(m$cds, 1L, anchor),
        substring(m$cds, anchor + 2L)
      )
    } else {
      inserted <- if (m$site == 'first') 'C' else 'G'
      ref <- substr(m$cds, anchor, anchor)
      alt <- paste0(ref, inserted)
      expected_cds <- paste0(
        substr(m$cds, 1L, anchor),
        inserted,
        substring(m$cds, anchor + 1L)
      )
    }
    position <- if (m$strand == 1L) {
      m$cds_start + anchor - 1L
    } else {
      m$cds_end - anchor - nchar(ref) + 2L
    }
    if (m$strand == -1L) {
      ref <- rc(ref)
      alt <- rc(alt)
    }
    stopifnot(substr(genome, position, position + nchar(ref) - 1L) == ref)
    rebuilt <- paste0(
      substr(genome, 1L, position - 1L),
      alt,
      substring(genome, position + nchar(ref))
    )
    tx <- substr(rebuilt, m$tx_start, m$tx_end + nchar(alt) - nchar(ref))
    if (m$strand == -1L) {
      tx <- rc(tx)
    }
    stopifnot(tx == paste0(expected_cds, m$post_cds))
    ref_peptide <- if (m$site == 'shifted') 'KSQ*' else 'KPQ*'
    alt_peptide <- switch(
      m$site,
      first = 'TTSIK*',
      later = 'KASIK*',
      shifted = 'NHNK*'
    )
    stopifnot(
      translate_fixture(m$cds) == ref_peptide,
      translate_fixture(tx) == alt_peptide
    )
    records[[i]] <- data.frame(
      event_index = i,
      seq_region = m$seq_region,
      transcript_index = m$transcript_index,
      chrom = m$chrom,
      position = position,
      id = id,
      reference = ref,
      alternate = alt,
      expected_cds = expected_cds,
      ref_peptide = ref_peptide,
      alt_peptide = alt_peptide,
      strand = m$strand,
      site = m$site
    )
  }
  records <- do.call(rbind, records)
  writeLines(unlist(fasta), file.path(out, 'reference.fa'))
  writeLines(c('##gff-version 3', unlist(gff)), file.path(out, 'model.gff3'))
  writeLines(
    c(
      '##fileformat=VCFv4.2',
      paste0(
        '##contig=<ID=',
        models$chrom,
        ',length=',
        models$sequence_length,
        '>'
      ),
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
      with(
        records,
        paste(
          chrom,
          position,
          id,
          reference,
          alternate,
          '.',
          'PASS',
          '.',
          sep = '\t'
        )
      )
    ),
    file.path(out, 'input.vcf')
  )
  cmd('samtools', c('faidx', file.path(out, 'reference.fa')), 'faidx')
  cmd('bgzip', file.path(out, 'model.gff3'), 'bgzip')
  cmd('tabix', c('-p', 'gff', file.path(out, 'model.gff3.gz')), 'tabix')
  # The GFF reader does not transfer cds_start_NF to transcript attributes.
  # Install the fixture attribute after parsing, before VEP computes consequences.
  vep_lib <- dirname(vep)
  old_env <- Sys.getenv(c('PERL5LIB', 'PERL5OPT'), unset = NA_character_)
  on.exit(
    {
      Sys.unsetenv(names(old_env))
      present <- !is.na(old_env)
      if (any(present)) do.call(Sys.setenv, as.list(old_env[present]))
    },
    add = TRUE
  )
  Sys.setenv(
    PERL5LIB = paste(
      normalizePath('test/duckvep/conformance'),
      vep_lib,
      file.path(vep_lib, 'modules'),
      Sys.getenv('PERL5LIB'),
      sep = ':'
    ),
    PERL5OPT = '-MFirstFrameshiftFixture'
  )
  cmd(
    'micromamba',
    c(
      'run',
      '-p',
      prefix,
      'env',
      paste0('PERL5LIB=', Sys.getenv('PERL5LIB')),
      paste0('PERL5OPT=', Sys.getenv('PERL5OPT')),
      'perl',
      vep,
      '--gff',
      file.path(out, 'model.gff3.gz'),
      '--fasta',
      file.path(out, 'reference.fa'),
      '--format',
      'vcf',
      '--vcf',
      '--hgvs',
      '--no_stats',
      '--force_overwrite',
      '--input_file',
      file.path(out, 'input.vcf'),
      '--output_file',
      file.path(out, 'vep.vcf')
    ),
    'vep'
  )
  Sys.unsetenv(c('PERL5LIB', 'PERL5OPT'))
  lines <- readLines(file.path(out, 'vep.vcf'))
  stopifnot(any(startsWith(lines, '##VEP="v116.0" API="v116"')))
  columns <- strsplit(
    sub(
      '.*Format: ([^"]+).*',
      '\\1',
      lines[startsWith(lines, '##INFO=<ID=CSQ,')]
    ),
    '|',
    fixed = TRUE
  )[[1L]]
  rows <- strsplit(lines[!startsWith(lines, '#')], '\t', fixed = TRUE)
  stopifnot(
    length(rows) == nrow(records),
    identical(vapply(rows, `[`, '', 3L), records$id)
  )
  oracle <- vapply(
    rows,
    function(x) {
      fields <- strsplit(sub('^CSQ=', '', x[8L]), '|', fixed = TRUE)[[1L]]
      length(fields) <- length(columns)
      utils::URLdecode(fields[match('HGVSp', columns)])
    },
    ''
  )
  oracle <- sub('^.*:p\\.', 'p.', oracle)
  con <- dbConnect(duckdb(
    config = list(allow_unsigned_extensions = 'true'),
    shared_home = FALSE
  ))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(extension)))
  dbWriteTable(con, 'models', models)
  dbWriteTable(con, 'events', records)
  tx <- paste(
    'SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,',
    'tx_start::UBIGINT transcript_start,tx_end::UBIGINT transcript_end,strand::TINYINT strand,',
    'transcript_index::UINTEGER gene_index,19::UBIGINT transcript_flags,',
    'cds_start::UBIGINT cds_start,cds_end::UBIGINT cds_end,cds::BLOB cds_sequence,',
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,post_cds::BLOB post_cds_sequence",
    'FROM models ORDER BY transcript_index'
  )
  ex <- paste(
    'SELECT transcript_index::UINTEGER transcript_index,tx_start::UBIGINT exon_start,',
    'tx_end::UBIGINT exon_end,1::UBIGINT exon_cdna_start,',
    '(tx_end-tx_start+1)::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase',
    'FROM models ORDER BY transcript_index'
  )
  regions <- 'SELECT seq_region::UINTEGER seq_region,sequence_length::UBIGINT sequence_length,chrom::VARCHAR seq_region_name FROM models ORDER BY seq_region'
  stopifnot(
    dbGetQuery(
      con,
      paste0(
        "SELECT loaded FROM duckvep_model_load('fixture',",
        q(regions),
        ',',
        q(tx),
        ',',
        q(ex),
        ',reference_fasta:=',
        q(file.path(out, 'reference.fa')),
        ')'
      )
    )$loaded
  )
  dbExecute(
    con,
    paste(
      'CREATE TABLE input_events AS SELECT event_index,seq_region,position,reference,alternate,',
      'NULL::UBIGINT end_position,NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,',
      'NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position FROM events'
    )
  )
  actual <- dbGetQuery(
    con,
    "SELECT event_index,transcript_hgvs,protein_hgvs,hgvs_shift FROM query(duckvep_annotate_sql('input_events','fixture',struct_pack(hgvs:=true,upstream_distance:=0,downstream_distance:=0))) ORDER BY event_index"
  )
  stopifnot(
    nrow(actual) == nrow(records),
    identical(as.numeric(actual$event_index), as.numeric(records$event_index))
  )
  result <- data.frame(
    records,
    transcript_hgvs = actual$transcript_hgvs,
    hgvs_shift = actual$hgvs_shift,
    duckvep = actual$protein_hgvs,
    vep = oracle
  )
  write.csv(
    result,
    file.path(out, 'comparisons.csv'),
    row.names = FALSE,
    na = ''
  )
  print(result)
  expected <- ifelse(
    records$site == 'first',
    'p.Lys1ThrfsTer?',
    ifelse(records$site == 'later', 'p.Pro2AlafsTer5', 'p.Lys1AsnfsTer?')
  )
  stopifnot(identical(oracle, expected))
  stopifnot(all(actual$hgvs_shift[records$site == 'shifted'] > 0))
  if (args[3L] == 'fixed') {
    stopifnot(identical(actual$protein_hgvs, oracle))
  } else {
    baseline <- ifelse(
      records$site == 'first',
      'p.Lys1ThrfsTer6',
      ifelse(records$site == 'later', 'p.Pro2AlafsTer5', 'p.Lys1AsnfsTer5')
    )
    stopifnot(identical(actual$protein_hgvs, baseline))
  }
}
main()
