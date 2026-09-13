# Field contracts for provider-excluded transcript comparisons. VEP 116 owns
# biological presentation; FastVEP's pinned output.rs owns native tab spelling.
duckvep_fastvep_fields <- function(contract) {
  tab <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature",
    "Feature_type", "Consequence", "cDNA_position", "CDS_position",
    "Protein_position", "Amino_acids", "Codons", "Existing_variation",
    "IMPACT", "DISTANCE", "STRAND", "FLAGS")
  csq <- c("Uploaded_variation", "Allele", "Consequence", "IMPACT", "SYMBOL",
    "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc",
    "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "CANONICAL",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP",
    "REF_ALLELE", "UPLOADED_ALLELE", "HGVS_OFFSET")
  switch(contract, operational17 = tab, native_tab17 = tab, vep_csq = csq,
    stop("unknown output contract: ", contract, call. = FALSE))
}

duckvep_fastvep_identity_fields <- c("record_index", "alt_index")

duckvep_fastvep_tab_header <- function(path, fields) {
  connection <- file(path, "rt")
  on.exit(close(connection), add = TRUE)
  skip <- 0L
  repeat {
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(header)) stop("tab output has no header")
    if (!startsWith(header, "##")) break
    skip <- skip + 1L
  }
  columns <- strsplit(sub("^#", "", header), "\t", fixed = TRUE)[[1L]]
  if (!identical(columns, fields)) stop("tab output differs from declared field schema")
  skip
}

duckvep_fastvep_field_sources <- function(contract) {
  fields <- duckvep_fastvep_fields(contract)
  source <- stats::setNames(rep("typed SQL projection", length(fields)), fields)
  assign <- function(names, value) source[intersect(names, fields)] <<- value
  assign(c("Uploaded_variation", "Location", "REF_ALLELE", "UPLOADED_ALLELE"), "physical input record")
  assign(c("Consequence", "IMPACT", "HGVSc", "HGVSp", "HGVS_OFFSET"), "native annotation and SQL formatting")
  assign(c("Gene", "Feature", "Feature_type", "STRAND", "SYMBOL", "BIOTYPE", "CANONICAL",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP"), "cold source metadata")
  assign("Existing_variation", "provider excluded")
  if (contract == "native_tab17") assign("FLAGS", "cold source canonical status")
  if (contract == "native_tab17") assign("Allele", "physical input record")
  if (contract == "operational17") {
    assign(c("Codons", "DISTANCE", "FLAGS"), "placeholder")
    assign(c("cDNA_position", "CDS_position", "Protein_position", "Amino_acids"), "compact scalar display")
    assign("Allele", "raw input ALT")
  }
  data.frame(contract, field = fields, source = unname(source))
}

duckvep_fastvep_csq_sql <- function(field, value) {
  # VEP116 OutputFactory/VCF.pm:401-404; SYMBOL and HGVSp are pre-escaped
  # by OutputFactory.pm:1526/1757. This is not URL decoding or JSON escaping.
  if (field == "HGVSp") value <- paste0("replace(", value, ", '=', '%3D')")
  if (field == "SYMBOL") value <- paste0("regexp_replace(", value, ", '\\s', '%20', 'g')")
  if (field != "Allele") value <- paste0("nullif(", value, ", '-')")
  paste0("replace(regexp_replace(replace(replace(coalesce(", value,
    ", ''), ',', '&'), ';', '%3B'), '\\s+', '_', 'g'), '|', '&')")
}

duckvep_fastvep_prepare_fields <- function(con, input, contract, distance, gff3 = "") {
  stopifnot(contract %in% c("native_tab17", "vep_csq"), length(distance) == 1L,
    is.numeric(distance), is.finite(distance), distance == floor(distance),
    distance >= 0, distance <= 2^32 - 1)
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  execute <- function(sql) invisible(DBI::dbExecute(con, sql))
  execute("SET preserve_insertion_order = true")
  execute(paste0("CREATE TEMP TABLE fastvep_source AS SELECT
    row_number() OVER ()::UBIGINT AS record_index, CHROM AS chrom,
    POS::UBIGINT AS position, ID AS variant_id, REF AS reference, ALT AS alternates
    FROM read_bcf(", q(input), ", scan_mode := 'sequential', decompression_threads := 0)"))
  execute("SET preserve_insertion_order = false")
  # fastvep-io/src/vcf.rs removes one shared first base using every ALT of
  # the physical record; IDs and Location are not per-ALT normalized keys.
  execute("CREATE TEMP TABLE fastvep_source_spelling AS
    WITH anchored AS (
      SELECT *, len(list_filter(alternates, x -> len(x) != len(reference))) > 0
        AND len(list_filter(alternates, x -> starts_with(x, '<') AND ends_with(x, '>')
          AND x NOT IN ('<NON_REF>', '<*>'))) = 0
        AND len(list_filter(alternates, x -> NOT contains(x, '*'))) > 0
        AND len(list_filter(alternates, x -> NOT contains(x, '*')
          AND left(x, 1) != left(reference, 1))) = 0 AS strip_anchor
      FROM fastvep_source
    ), alleles AS (
      SELECT *, position + strip_anchor::UBIGINT AS native_start,
        position + len(reference) - 1 AS native_end,
        CASE WHEN strip_anchor THEN coalesce(nullif(substring(reference, 2), ''), '-')
          ELSE reference END AS native_reference,
        list_transform(alternates, x -> CASE WHEN strip_anchor AND NOT contains(x, '*')
          THEN coalesce(nullif(substring(x, 2), ''), '-') ELSE x END) AS native_alternates
      FROM anchored
    ) SELECT *, chrom || ':' || native_start::VARCHAR || CASE
      WHEN native_start = native_end THEN '' ELSE '-' || native_end::VARCHAR END AS native_location
    FROM alleles")
  execute("CREATE TEMP TABLE fastvep_events AS
    SELECT row_number() OVER (ORDER BY s.record_index, a.alt_index)::UBIGINT AS event_index,
      s.record_index, a.alt_index, r.seq_region, s.chrom, s.position, s.variant_id,
      upper(s.reference) AS reference, upper(a.alternate) AS alternate, s.alternates,
      s.reference AS uploaded_reference, s.native_reference, s.native_alternates, s.native_location,
      NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
      NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
      NULL::UBIGINT AS mate_position
    FROM fastvep_source_spelling s
    CROSS JOIN UNNEST(s.alternates) WITH ORDINALITY a(alternate, alt_index)
    LEFT JOIN duckvep_bench_regions r ON r.name = s.chrom
      OR r.name = regexp_replace(s.chrom, '^chr', '')
    WHERE regexp_full_match(s.reference, '[ACGTNacgtn]+')
      AND regexp_full_match(a.alternate, '[ACGTNacgtn]+')
      AND upper(s.reference) <> upper(a.alternate)
    ORDER BY seq_region, position, record_index, alt_index")
  invalid <- DBI::dbGetQuery(con, "SELECT count(*) n FROM (
    SELECT record_index, alt_index FROM fastvep_events
    GROUP BY record_index, alt_index HAVING count(*) != 1 OR count(seq_region) != 1)")$n
  if (invalid != 0) stop("literal input alleles need exactly one model region", call. = FALSE)
  execute("CREATE TEMP VIEW fastvep_ordered_events AS SELECT * FROM fastvep_events
    ORDER BY seq_region, position, record_index, alt_index")
  execute(paste0("CREATE TEMP TABLE fastvep_annotations AS SELECT
    event_index, transcript_index, consequence, consequence_mask, impact,
    transcript_hgvs, protein_hgvs, hgvs_shift FROM duckvep_annotate(
    'fastvep_ordered_events', 'fastvep_comparison', rich := TRUE, hgvs := ",
    if (contract == "vep_csq") "TRUE" else "FALSE", ", upstream_distance := ", distance,
    ", downstream_distance := ", distance, ")"))

  execute("CREATE TEMP TABLE fastvep_metadata AS SELECT transcript_index,
    NULL::VARCHAR AS symbol, NULL::BOOLEAN AS canonical, NULL::VARCHAR AS tsl,
    NULL::VARCHAR AS appris, NULL::VARCHAR AS ccds
    FROM duckvep_bench_model.model_transcripts")
  core <- DBI::dbGetQuery(con, "SELECT table_name FROM information_schema.tables
    WHERE table_catalog = 'duckvep_bench_model' AND table_schema = 'ensembl_core'")$table_name
  if (all(c("transcript_attrib", "attrib_type", "gene") %in% core)) {
    conflicts <- DBI::dbGetQuery(con, "SELECT count(*) n FROM (
      SELECT t.transcript_id, a.code
      FROM duckvep_bench_model.ensembl_core.transcript_attrib t
      JOIN duckvep_bench_model.ensembl_core.attrib_type a USING(attrib_type_id)
      WHERE a.code IN ('TSL', 'appris', 'ccds_transcript')
      GROUP BY t.transcript_id, a.code HAVING count(DISTINCT t.value) > 1)")$n
    if (conflicts != 0) stop("transcript display attributes have conflicting scalar values", call. = FALSE)
    execute("CREATE TEMP TABLE fastvep_attributes AS
      SELECT t.transcript_id,
        max(t.value) FILTER (WHERE a.code = 'TSL') AS tsl,
        max(t.value) FILTER (WHERE a.code = 'appris') AS appris,
        max(t.value) FILTER (WHERE a.code = 'ccds_transcript') AS ccds
      FROM duckvep_bench_model.ensembl_core.transcript_attrib t
      JOIN duckvep_bench_model.ensembl_core.attrib_type a USING(attrib_type_id)
      GROUP BY t.transcript_id")
    execute("UPDATE fastvep_metadata m SET canonical =
      t.source_transcript_id = g.canonical_transcript_id,
      tsl = a.tsl, appris = a.appris, ccds = a.ccds
      FROM duckvep_bench_model.model_transcripts t
      JOIN duckvep_bench_model.ensembl_core.gene g ON g.gene_id = t.source_gene_id
      LEFT JOIN fastvep_attributes a ON a.transcript_id = t.source_transcript_id
      WHERE m.transcript_index = t.transcript_index")
  }
  if (nzchar(gff3)) {
    execute(paste0("CREATE TEMP TABLE fastvep_gff AS SELECT feature, attributes_map AS attributes
      FROM read_gff(", q(gff3), ", attributes_map := TRUE, scan_mode := 'sequential')
      WHERE feature NOT IN ('exon', 'CDS', 'chromosome', 'biological_region',
        'five_prime_UTR', 'three_prime_UTR')"))
    execute("UPDATE fastvep_metadata m SET symbol = g.attributes['Name']
      FROM duckvep_bench_model.model_transcripts t JOIN fastvep_gff g
      ON (g.feature IN ('gene', 'pseudogene') OR ends_with(g.feature, '_gene'))
        AND regexp_replace(g.attributes['ID'], '^gene:', '') = t.gene_stable_id
      WHERE m.transcript_index = t.transcript_index")
    execute("UPDATE fastvep_metadata m SET
      canonical = coalesce(m.canonical, list_contains(string_split(g.attributes['tag'], ','), 'Ensembl_canonical'), false),
      tsl = coalesce(m.tsl, g.attributes['transcript_support_level']),
      ccds = coalesce(m.ccds, g.attributes['ccdsid'])
      FROM duckvep_bench_model.model_transcripts t JOIN fastvep_gff g
      ON regexp_replace(g.attributes['ID'], '^transcript:', '') = t.transcript_stable_id
      WHERE m.transcript_index = t.transcript_index")
  }
  if (DBI::dbGetQuery(con, "SELECT count(*) n FROM fastvep_metadata WHERE canonical IS NULL")$n != 0) {
    stop("canonical status requires Ensembl core gene metadata or matching GFF transcripts", call. = FALSE)
  }
}

duckvep_fastvep_field_query <- function(con, contract, include_identity = FALSE) {
  stopifnot(contract %in% c("native_tab17", "vep_csq"))
  is_csq <- contract == "vep_csq"
  range <- function(first, last, absent = "''") paste0("CASE WHEN ", first,
    " IS NULL AND ", last, " IS NULL THEN ", absent,
    " WHEN ", first, " = ", last, " THEN ", first,
    "::VARCHAR ELSE coalesce(", first, "::VARCHAR, '?') || '-' || coalesce(",
    last, "::VARCHAR, '?') END")
  pair <- function(ref, alt) paste0("CASE WHEN ", ref, " IS NULL AND ", alt,
    " IS NULL THEN NULL ELSE coalesce(", ref, ", '') || '/' || coalesce(", alt, ", '') END")
  aa <- pair("p.reference_amino_acids", "p.alternate_amino_acids")
  if (is_csq) aa <- paste0("CASE WHEN p.reference_amino_acids = p.alternate_amino_acids
    THEN p.reference_amino_acids ELSE ", aa, " END")
  accession <- function(stable, version) paste0(stable,
    " || CASE WHEN ", version, " IS NULL THEN '' ELSE '.' || ", version, "::VARCHAR END")
  expressions <- c(
    Uploaded_variation = "coalesce(e.variant_id, e.native_location || '_' ||
      e.native_reference || '/' || array_to_string(e.native_alternates, '/'))",
    Location = "e.native_location",
    Allele = if (is_csq) "p.output_allele" else "e.native_alternates[e.alt_index]",
    Gene = "t.gene_stable_id", Feature = "t.transcript_stable_id",
    Feature_type = "CASE WHEN a.transcript_index IS NOT NULL THEN 'Transcript' END",
    Consequence = if (is_csq) "a.consequence" else "replace(a.consequence, '&', ',')",
    cDNA_position = range("p.cdna_start", "p.cdna_end"),
    CDS_position = range("p.cds_start", "p.cds_end"),
    Protein_position = range("p.protein_start", "p.protein_end"),
    Amino_acids = aa, Codons = pair("p.reference_codons", "p.alternate_codons"),
    Existing_variation = "NULL::VARCHAR", IMPACT = "a.impact",
    DISTANCE = "p.transcript_distance::VARCHAR", STRAND = "t.strand::VARCHAR",
    FLAGS = if (is_csq) "concat_ws('&', CASE WHEN p.cds_start_nf THEN 'cds_start_NF' END,
      CASE WHEN p.cds_end_nf THEN 'cds_end_NF' END)" else "CASE WHEN m.canonical THEN 'canonical' END",
    SYMBOL = "m.symbol", BIOTYPE = "t.transcript_biotype",
    EXON = paste0("CASE WHEN p.exon_first IS NOT NULL THEN (", range("p.exon_first", "p.exon_last"), ") || '/' || p.exon_total END"),
    INTRON = paste0("CASE WHEN p.intron_first IS NOT NULL THEN (", range("p.intron_first", "p.intron_last"), ") || '/' || p.intron_total END"),
    HGVSc = paste0("CASE WHEN a.transcript_hgvs IS NOT NULL THEN (",
      accession("t.transcript_stable_id", "t.transcript_version"), ") || ':' || a.transcript_hgvs END"),
    HGVSp = paste0("CASE WHEN a.protein_hgvs IS NOT NULL THEN (",
      accession("t.translation_stable_id", "t.translation_version"), ") || ':' || a.protein_hgvs END"),
    CANONICAL = "CASE WHEN m.canonical THEN 'YES' END",
    MANE_SELECT = "t.mane_select_refseq", MANE_PLUS_CLINICAL = "t.mane_plus_clinical_refseq",
    TSL = "regexp_extract(m.tsl, '^(?:tsl)?([0-9]+)', 1)",
    APPRIS = "replace(replace(m.appris, 'principal', 'P'), 'alternative', 'A')", CCDS = "m.ccds",
    ENSP = "t.translation_stable_id",
    # Parser.pm::minimise_alleles retains the shared physical REF for multiple
    # ALTs and minimises unequal-length biallelic alleles. OutputFactory reads
    # REF_ALLELE from that feature, not the individual transcript allele.
    REF_ALLELE = "CASE WHEN length(e.alternates) > 1 THEN e.native_reference
      WHEN length(e.reference) = length(e.alternate) THEN e.reference
      ELSE coalesce(nullif(substring(e.reference, geometry.reference_difference_offset + 1,
        geometry.reference_difference_length), ''), '-') END",
    UPLOADED_ALLELE = "e.uploaded_reference || '/' || array_to_string(e.alternates, ',')",
    # VEP OutputFactory multiplies the transcript-oriented shift by strand.
    HGVS_OFFSET = "CASE WHEN a.transcript_hgvs IS NOT NULL OR a.protein_hgvs IS NOT NULL
      THEN (nullif(a.hgvs_shift, 0)::BIGINT * t.strand)::VARCHAR END")
  fields <- duckvep_fastvep_fields(contract)
  values <- vapply(fields, function(field) {
    value <- paste0("(", expressions[[field]], ")")
    if (is_csq && field != "Uploaded_variation") {
      value <- duckvep_fastvep_csq_sql(field, value)
    } else if (!is_csq) {
      if (field %in% c("cDNA_position", "CDS_position", "Protein_position")) {
        value <- paste0("CASE WHEN a.transcript_index IS NULL THEN '-' ELSE ", value, " END")
      } else if (field %in% c("IMPACT", "DISTANCE", "STRAND", "FLAGS")) {
        value <- paste0("CASE WHEN a.transcript_index IS NULL THEN '-' ELSE coalesce(", value, ", '-') END")
      } else value <- paste0("coalesce(", value, ", '-')")
    }
    paste(value, "AS", DBI::dbQuoteIdentifier(con, field))
  }, character(1L))
  if (include_identity) values <- c("e.record_index", "e.alt_index", values)
  paste0("SELECT ", paste(values, collapse = ",\n"), "
    FROM fastvep_annotations a
    JOIN fastvep_events e USING(event_index)
    JOIN duckvep_transcript_projection('fastvep_events', 'fastvep_annotations',
      'duckvep_bench_model.model_transcripts') p
      ON p.event_index = a.event_index AND p.transcript_index IS NOT DISTINCT FROM a.transcript_index
    LEFT JOIN duckvep_bench_model.model_transcripts t ON t.transcript_index = a.transcript_index
    LEFT JOIN fastvep_metadata m ON m.transcript_index = a.transcript_index
    CROSS JOIN LATERAL (SELECT duckvep_allele_geometry(e.position, e.reference, e.alternate) AS geometry)")
}
