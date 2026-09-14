#!/usr/bin/env Rscript
# Small source-backed worker smoke; no oracle or whole-genome cache is required.
library(DBI)
library(duckdb)
local({
  root <- normalizePath(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE))
  source(file.path(root, "benchmarks/benchmark_duckvep_fastvep_fields.R"))
  for (contract in c("operational17", "native_tab17", "vep_csq")) {
    classification <- duckvep_fastvep_field_sources(contract)
    stopifnot(
      identical(classification$field, duckvep_fastvep_fields(contract)),
      !anyNA(classification$source), !anyDuplicated(classification$field)
    )
    if (contract != "operational17") stopifnot(!any(classification$source == "placeholder"))
  }
  con <- dbConnect(duckdb())
  escape <- function(field, value) {
    DBI::dbGetQuery(con, paste(
      "SELECT",
      duckvep_fastvep_csq_sql(field, as.character(DBI::dbQuoteString(con, value))), "AS value"
    ))$value
  }
  # These expected bytes follow pinned production OutputFactory, not FastVEP's
  # differently escaping test-only helper. Literal percent sequences stay literal.
  stopifnot(
    escape("SYMBOL", "A B,C;D|E") == "A%20B&C%3BD&E",
    escape("HGVSp", "ENSP:p.Val2=") == "ENSP:p.Val2%3D",
    escape("CCDS", "A\tB%3B") == "A_B%3B",
    escape("Allele", "-") == "-", escape("Codons", "-") == "",
    escape("Codons", NA_character_) == "", escape("Codons", "") == ""
  )
  dbDisconnect(con, shutdown = TRUE)
  extension <- normalizePath(file.path(root, "build/release/duckhts.duckdb_extension"))
  directory <- tempfile("fastvep-fields-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  model <- file.path(directory, "model.duckdb")
  con <- dbConnect(duckdb(), dbdir = model)
  sql <- paste(readLines(file.path(root, "test/duckvep/conformance/minimal_model.sql")), collapse = "\n")
  for (statement in strsplit(sql, ";", fixed = TRUE)[[1L]]) {
    if (nzchar(trimws(statement))) dbExecute(con, statement)
  }
  index <- read.delim(file.path(root, "test/data/duckvep/minimal.fa.fai"), header = FALSE)
  dbWriteTable(con, "fasta_index", data.frame(name = index[[1L]], sequence_length = index[[2L]]))
  dbExecute(con, "CREATE TABLE model_regions AS SELECT seq_region, name AS seq_region_name,
  sequence_length::UBIGINT AS sequence_length FROM duckvep_sequence_regions JOIN fasta_index USING(name)")
  dbExecute(con, "CREATE TABLE model_transcripts AS SELECT t.*,
  'DUCK1' AS gene_stable_id, 'DUCK1-201' AS transcript_stable_id,
  NULL::BIGINT AS transcript_version, 'protein_coding' AS transcript_biotype,
  'DUCK1-P' AS translation_stable_id, NULL::BIGINT AS translation_version,
  NULL::VARCHAR AS mane_select_refseq, NULL::VARCHAR AS mane_plus_clinical_refseq,
  list(struct_pack(exon_start := e.exon_start, exon_end := e.exon_end,
    exon_cdna_start := e.exon_cdna_start, exon_cdna_end := e.exon_cdna_end,
    phase := e.phase, end_phase := e.end_phase) ORDER BY e.exon_cdna_start) AS exons,
  []::STRUCT(mature_mirna_start UBIGINT, mature_mirna_end UBIGINT)[] AS mature_mirna_regions,
  []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] AS peptide_edits
  FROM duckvep_transcripts t JOIN duckvep_exons e USING(transcript_index) GROUP BY ALL")
  dbDisconnect(con, shutdown = TRUE)

  outputs <- list()
  for (contract in c("operational17", "native_tab17", "vep_csq")) {
    output <- file.path(directory, paste0(contract, ".tsv"))
    arguments <- c(
      file.path(root, "benchmarks/benchmark_duckvep_fastvep_worker.R"),
      "--extension", extension, "--model", model, "--input",
      file.path(root, "test/data/duckvep/minimal_bcsq.vcf"), "--output", output,
      "--output-contract", contract
    )
    if (contract != "operational17") {
      arguments <- c(
        arguments,
        "--gff3", file.path(root, "test/data/duckvep/minimal.gff3")
      )
    }
    if (contract == "vep_csq") {
      arguments <- c(
        arguments,
        "--fasta", file.path(root, "test/data/duckvep/minimal.fa")
      )
    }
    stopifnot(system2("Rscript", shQuote(arguments)) == 0L)
    outputs[[contract]] <- read.delim(output,
      check.names = FALSE, colClasses = "character",
      quote = "", na.strings = character()
    )
    stopifnot(identical(names(outputs[[contract]]), duckvep_fastvep_fields(contract)))
  }
  tab <- outputs$native_tab17
  csq <- outputs$vep_csq
  stopifnot(
    nrow(tab) == 3L, nrow(csq) == 3L,
    # The fixture has genomic 123..125 = GTA; 124 T>C changes GTA to GCA.
    tab$Codons[tab$Uploaded_variation == "duck_missense"] == "gTa/gCa",
    tab$FLAGS[tab$Uploaded_variation == "duck_missense"] == "-",
    tab$cDNA_position[tab$Uploaded_variation == "duck_intron"] == "",
    csq$cDNA_position[csq$Uploaded_variation == "duck_intron"] == "",
    csq$FLAGS[csq$Uploaded_variation == "duck_missense"] == "",
    tab$IMPACT[tab$Uploaded_variation == "duck_intergenic"] == "-",
    csq$SYMBOL[csq$Uploaded_variation == "duck_missense"] == "DUCK1",
    csq$EXON[csq$Uploaded_variation == "duck_missense"] == "1/2",
    all(csq$Existing_variation == ""),
    all(tab$Existing_variation == "-")
  )

  # Duplicate physical records and multiallelic rows keep original source ordinals;
  # a missing ID uses the complete physical-record allele list, not one ALT.
  source_vcf <- readLines(file.path(root, "test/data/duckvep/minimal_bcsq.vcf"))
  records <- source_vcf[!startsWith(source_vcf, "#")]
  multi <- strsplit(records[[1L]], "\t", fixed = TRUE)[[1L]]
  multi[[3L]] <- "."
  multi[[5L]] <- "C,G"
  input <- file.path(directory, "identity.vcf")
  writeLines(c(source_vcf, records[[1L]], paste(multi, collapse = "\t")), input)
  output <- file.path(directory, "identity.tsv")
  arguments <- c(
    file.path(root, "benchmarks/benchmark_duckvep_fastvep_worker.R"),
    "--extension", extension, "--model", model, "--input", input, "--output", output,
    "--output-contract", "native_tab17", "--threads", "4", "--include-identity", "--gff3",
    file.path(root, "test/data/duckvep/minimal.gff3")
  )
  stopifnot(system2("Rscript", shQuote(arguments)) == 0L)
  observed <- read.delim(output,
    check.names = FALSE, colClasses = "character",
    quote = "", na.strings = character()
  )
  stopifnot(
    nrow(observed) == 6L,
    identical(names(observed), c(duckvep_fastvep_identity_fields, duckvep_fastvep_fields("native_tab17"))),
    identical(sort(observed$record_index[observed$Uploaded_variation == "duck_missense"]), c("1", "4")),
    all(observed$Uploaded_variation[observed$record_index == "5"] == "chrDuck:124_T/C/G"),
    identical(sort(observed$alt_index[observed$record_index == "5"]), c("1", "2"))
  )

  # Pinned VEP116 retains the shared first-anchor-stripped REF for both ALTs.
  # Keep the complete witness here: per-ALT minimisation would yield GTA/TA
  # instead of TA/TA for this physical record.
  input <- file.path(directory, "multiallelic_reference.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chrDuck,length=260>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chrDuck\t123\tshared_ref\tGTA\tGCA,G\t.\tPASS\t."
  ), input)
  output <- file.path(directory, "multiallelic_reference.tsv")
  arguments <- c(
    file.path(root, "benchmarks/benchmark_duckvep_fastvep_worker.R"),
    "--extension", extension, "--model", model, "--input", input, "--output", output,
    "--output-contract", "vep_csq", "--include-identity",
    "--gff3", file.path(root, "test/data/duckvep/minimal.gff3"),
    "--fasta", file.path(root, "test/data/duckvep/minimal.fa")
  )
  stopifnot(system2("Rscript", shQuote(arguments)) == 0L)
  observed <- read.delim(output,
    check.names = FALSE, colClasses = "character",
    quote = "", na.strings = character()
  )
  stopifnot(
    nrow(observed) == 2L, all(observed$record_index == "1"),
    identical(sort(observed$alt_index), c("1", "2")),
    all(observed$REF_ALLELE == "TA"), all(observed$UPLOADED_ALLELE == "GTA/GCA&G")
  )

  # Ensembl gene features include ncRNA_gene and pseudogene, with ID=gene:...
  # independently of the transcript's own feature class. Only the label changes.
  source_gff <- readLines(file.path(root, "test/data/duckvep/minimal.gff3"))
  gene_row <- which(grepl("\tgene\t", source_gff, fixed = TRUE))
  stopifnot(length(gene_row) == 1L)
  for (feature in c("ncRNA_gene", "pseudogene", "miRNA_gene")) {
    gff <- source_gff
    fields <- strsplit(gff[[gene_row]], "\t", fixed = TRUE)[[1L]]
    stopifnot(startsWith(fields[[9L]], "ID=gene:"))
    fields[[3L]] <- feature
    gff[[gene_row]] <- paste(fields, collapse = "\t")
    path <- file.path(directory, paste0(feature, ".gff3"))
    writeLines(gff, path)
    arguments[[match("--gff3", arguments) + 1L]] <- path
    arguments[[match("--output", arguments) + 1L]] <- file.path(directory, paste0(feature, ".tsv"))
    stopifnot(system2("Rscript", shQuote(arguments)) == 0L)
    observed <- read.delim(arguments[[match("--output", arguments) + 1L]],
      check.names = FALSE, colClasses = "character", quote = "", na.strings = character()
    )
    stopifnot(nrow(observed) == 2L, all(observed$SYMBOL == "DUCK1"))
  }

  # Exercise the complete final projection, including every output field, with
  # a fixed source/transcript multiplicity and an explicit no-spill budget.
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true",
    threads = "1", memory_limit = "192MB", max_temp_directory_size = "0B")))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, paste("LOAD", dbQuoteString(con, extension)))
  dbExecute(con, "ATTACH ':memory:' AS duckvep_bench_model")
  dbExecute(con, "CREATE TABLE duckvep_bench_model.model_transcripts AS SELECT
    i::UINTEGER transcript_index, 1::UINTEGER seq_region, 1::UBIGINT transcript_start,
    12000::UBIGINT transcript_end, 1::TINYINT strand, 3::UBIGINT transcript_flags,
    1::UBIGINT cds_start, 12000::UBIGINT cds_end,
    ('ATG'||repeat('AAA',3998)||'TAA')::BLOB cds_sequence,
    'ACGT'::BLOB post_cds_sequence, 1::UTINYINT codon_table,
    [struct_pack(exon_start:=1::UBIGINT, exon_end:=12000::UBIGINT,
      exon_cdna_start:=1::UBIGINT, exon_cdna_end:=12000::UBIGINT,
      phase:=0::TINYINT, end_phase:=0::TINYINT)] exons,
    []::STRUCT(protein_position UINTEGER, alternate_amino_acid VARCHAR, edit_code VARCHAR)[] peptide_edits,
    'GENE' AS gene_stable_id, 'TX'||i::VARCHAR AS transcript_stable_id,
    1::BIGINT AS transcript_version, 'protein_coding' AS transcript_biotype,
    'P'||i::VARCHAR AS translation_stable_id, 1::BIGINT AS translation_version,
    NULL::VARCHAR AS mane_select_refseq, NULL::VARCHAR AS mane_plus_clinical_refseq
    FROM range(64) r(i)")
  dbExecute(con, "CREATE TEMP TABLE fastvep_events AS SELECT
    (i+1)::UBIGINT AS event_index, (i+1)::UBIGINT AS record_index, 1::BIGINT AS alt_index,
    1::UINTEGER AS seq_region, 'chr1' AS chrom, 4::UBIGINT AS position,
    'site_'||i::VARCHAR AS variant_id, 'A' AS reference, 'G' AS alternate, ['G'] AS alternates,
    'A' AS uploaded_reference, 'A' AS native_reference, ['G'] AS native_alternates,
    'chr1:4' AS native_location FROM range(4096) r(i)")
  dbExecute(con, "CREATE TEMP TABLE fastvep_annotations AS SELECT event_index, transcript_index,
    'missense_variant' AS consequence, (SELECT consequence_mask FROM duckvep_so_terms()
      WHERE consequence='missense_variant') AS consequence_mask, 'MODERATE' AS impact,
    'c.4A>G' AS transcript_hgvs, 'p.Lys2Glu' AS protein_hgvs, 0::BIGINT AS hgvs_shift
    FROM fastvep_events, duckvep_bench_model.model_transcripts")
  dbExecute(con, "CREATE TEMP TABLE fastvep_metadata AS SELECT transcript_index, 'GENE' AS symbol,
    TRUE AS canonical, NULL::VARCHAR AS tsl, NULL::VARCHAR AS appris, NULL::VARCHAR AS ccds
    FROM duckvep_bench_model.model_transcripts")
  for (contract in c("native_tab17", "vep_csq")) {
    query <- duckvep_fastvep_field_query(con, contract, include_identity = contract == "vep_csq")
    fields <- duckvep_fastvep_transport_fields(contract)
    full_hash <- paste0("hash(", paste(dbQuoteIdentifier(con, fields), collapse = ","), ")")
    result <- dbGetQuery(con, paste0("SELECT count(*) AS output_rows,",
      "count(DISTINCT Uploaded_variation) AS source_records, count(DISTINCT Feature) AS transcripts,",
      "count(DISTINCT (Uploaded_variation, Feature)) AS source_transcript_pairs,",
      "count(*) FILTER (WHERE Codons='Aaa/Gaa' AND Amino_acids='K/E'",
      " AND CDS_position='4' AND Protein_position='2') AS exact_rows,",
      "bit_xor(", full_hash, ")::VARCHAR AS fingerprint FROM (", query, ")"))
    stopifnot(result$output_rows == 262144, result$source_records == 4096,
      result$transcripts == 64, result$source_transcript_pairs == 262144,
      result$exact_rows == 262144, !is.na(result$fingerprint))
  }
  cat("FastVEP field projection smoke: OK\n")
})
