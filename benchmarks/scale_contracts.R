library(DBI)
library(duckdb)

args <- commandArgs(trailingOnly = TRUE)
record <- "--record" %in% args
model <- Sys.getenv("DUCKVEP_SCALE_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
extension <- Sys.getenv("DUCKVEP_SCALE_EXTENSION", "build/release/extension/duckvep/duckvep.duckdb_extension")
root <- normalizePath(".")
output <- file.path(root, "benchmarks/data/scale_contracts")
source_file <- file.path(root, "benchmarks/data/duckvep_fastvep/fields_replay_9bf888e/worker-*.parquet")
metadata <- file.path(output, "canonical-metadata.parquet")
if (!file.exists(model) || !file.exists(metadata)) stop("model or canonical metadata is missing")
if (!file.exists(extension)) stop("release extension is missing")
sha256 <- function(path) {
  result <- system2("sha256sum", path, stdout = TRUE)
  if (length(result) != 1L || !grepl("^[[:xdigit:]]{64} ", result)) stop("sha256sum failed")
  substr(result, 1L, 64L)
}
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
q <- function(x) as.character(dbQuoteString(con, x))
run <- function(sql) invisible(dbExecute(con, sql))
get <- function(sql) dbGetQuery(con, sql)
run(paste0("LOAD ", q(normalizePath(extension))))
run("INSTALL json")
run("LOAD json")
run("SET threads=1")
run("SET memory_limit='8GB'")
run(paste0("SET temp_directory=", q(file.path(tempdir(), "scale-spill"))))
run("SET preserve_insertion_order=false")
run(paste0("ATTACH ", q(normalizePath(model)), " AS duckvep_bench_model (READ_ONLY)"))
source("benchmarks/benchmark_duckvep_fastvep_fields.R", local = TRUE)

# The sorted length-prefixed JSON stream is versioned separately from DuckDB hash().
# JSON object keys follow the declared SELECT order; NULL fields are retained.
receipt <- function(name, query, fields, partitions = 16L) {
  run(paste0("CREATE OR REPLACE TEMP TABLE scale_rows AS SELECT a.*, hash(a) AS row_hash, ",
    "to_json(a) AS row_payload FROM (", query, ") a"))
  schema <- get("DESCRIBE SELECT * EXCLUDE(row_hash, row_payload) FROM scale_rows")
  if (!identical(schema$column_name, fields)) stop(name, ": schema differs from contract: ",
    paste(schema$column_name, collapse = ","), " expected ", paste(fields, collapse = ","))
  profile <- get("SELECT count(*)::VARCHAR AS rows, sum(row_hash::HUGEINT)::VARCHAR AS hash_sum, bit_xor(row_hash)::VARCHAR AS hash_xor FROM scale_rows")
  result <- data.frame(contract = name, partition = -1L, rows = profile$rows,
    hash_sum = profile$hash_sum, hash_xor = profile$hash_xor,
    sha256 = NA_character_)
  for (partition in seq_len(partitions) - 1L) {
    path <- tempfile("scale-partition-")
    sql <- paste0("COPY (SELECT octet_length(encode(row_payload))::VARCHAR || ':' || row_payload AS line ",
      "FROM scale_rows WHERE (event_index % ", partitions, ") = ", partition,
      " ORDER BY row_payload) TO ", q(path), " (FORMAT CSV, HEADER FALSE, QUOTE '')")
    run(sql)
    n <- get(paste0("SELECT count(*)::VARCHAR AS n FROM scale_rows WHERE (event_index % ",
      partitions, ") = ", partition))$n
    result <- rbind(result, data.frame(contract = name, partition = partition,
      rows = n, hash_sum = NA_character_, hash_xor = NA_character_, sha256 = sha256(path)))
    unlink(path)
  }
  schema$contract <- name
  list(receipt = result, schema = schema[, c("contract", "column_name", "column_type")])
}

model_sha <- sha256(model)
metadata_sha <- sha256(metadata)
source_tree <- system2("git", c("rev-parse", "HEAD:src"), stdout = TRUE)
if (length(source_tree) != 1L || !grepl("^[[:xdigit:]]{40}$", source_tree)) {
  stop("cannot resolve native source tree")
}
fields_sha <- sha256("benchmarks/benchmark_duckvep_fastvep_fields.R")
literal_gate <- "seq_region IS NOT NULL AND
 regexp_full_match(reference, '[ACGTNacgtn]+') AND
 regexp_full_match(alternate, '[ACGTNacgtn]+') AND
 upper(reference) <> upper(alternate)"
fixture <- get(paste0("SELECT ", literal_gate, " AS eligible FROM (VALUES
 (1, 'A', 'C'), (1, 'a', 'c'), (1, 'A', 'N'), (1, 'A', 'A'),
 (1, 'A', '<DEL>'), (1, 'A', '*'), (NULL, 'A', 'C'), (1, 'A', 'R'),
 (1, '-', 'C')) v(seq_region, reference, alternate)"))
if (!identical(as.logical(fixture$eligible), c(TRUE, TRUE, TRUE, rep(FALSE, 6L)))) {
  stop("literal source eligibility differs from contract")
}
input_summary <- get(paste0("SELECT count(*) AS physical_records,
 count(*) FILTER (WHERE ", literal_gate, ") AS eligible_literal,
 count(*) FILTER (WHERE NOT (", literal_gate, ")) AS unsupported,
 count(*) - count(DISTINCT (seq_region, position, reference, alternate)) AS duplicate_records
 FROM duckvep_bench_model.bench_variants"))
run("CREATE TEMP TABLE events AS SELECT row_number() OVER (ORDER BY seq_region, position, reference, alternate)::UBIGINT AS event_index, *,
 NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change,
 NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position
 FROM duckvep_bench_model.bench_variants ORDER BY seq_region, position, reference, alternate")
run("CREATE TEMP VIEW ordered_events AS SELECT * FROM events ORDER BY seq_region, position, event_index")
run("CREATE TEMP TABLE fastvep_events AS WITH anchored AS (
 SELECT e.*, r.name AS chrom, length(reference) != length(alternate) AND left(reference,1)=left(alternate,1) AS strip_anchor
 FROM events e JOIN duckvep_bench_model.duckvep_sequence_regions r USING(seq_region)
 ) SELECT * EXCLUDE(strip_anchor), event_index AS record_index, 1::BIGINT AS alt_index,
 NULL::VARCHAR AS variant_id, [alternate] AS alternates, reference AS uploaded_reference,
 CASE WHEN strip_anchor THEN coalesce(nullif(substr(reference,2),''),'-') ELSE reference END AS native_reference,
 [CASE WHEN strip_anchor THEN coalesce(nullif(substr(alternate,2),''),'-') ELSE alternate END] AS native_alternates,
 chrom || ':' || (position + strip_anchor::UBIGINT)::VARCHAR || CASE WHEN position + strip_anchor::UBIGINT = position + length(reference) - 1 THEN ''
 ELSE '-' || (position + length(reference) - 1)::VARCHAR END AS native_location FROM anchored")
run("CREATE TEMP VIEW fastvep_ordered_events AS SELECT * FROM fastvep_events ORDER BY seq_region, position, record_index, alt_index")
run(paste0("CREATE TEMP TABLE fastvep_metadata AS SELECT transcript_index, NULL::VARCHAR AS symbol,
 canonical, NULL::VARCHAR AS tsl, NULL::VARCHAR AS appris, NULL::VARCHAR AS ccds FROM read_parquet(", q(metadata), ")"))
if (get("SELECT count(*) AS n FROM fastvep_metadata WHERE canonical IS NULL")$n != 0) stop("missing canonical status")
run("SELECT * FROM duckvep_model_load('grch38',
 'SELECT seq_region, sequence_length FROM duckvep_bench_model.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM duckvep_bench_model.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM duckvep_bench_model.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM duckvep_bench_model.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM duckvep_bench_model.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 transcript_coverage_complete := TRUE)")
compact_fields <- c("event_index", "transcript_index", "gene_index", "consequence_mask", "region_mask",
  "impact_code", "status_code", "reason_code", "cdna_position", "cds_position",
  "protein_position", "reference_amino_acid_code", "alternate_amino_acid_code",
  "nmd_prediction_code", "nmd_escape_reasons", "regulation_feature_index", "overlap_object_code")
compact <- paste0("SELECT ", paste(compact_fields, collapse = ", "),
  " FROM duckvep_annotate('ordered_events', 'grch38')")
complete_fields <- c("event_index", "record_index", "alt_index", duckvep_fastvep_fields("native_tab17"))
complete <- duckvep_fastvep_field_query(con, "native_tab17", include_identity = TRUE)
complete <- sub("SELECT p.record_index", "SELECT p.event_index, p.record_index", complete, fixed = TRUE)
complete <- sub("'fastvep_comparison'", "'grch38'", complete, fixed = TRUE)
receipts <- list(receipt("compact", compact, compact_fields),
  receipt("complete17", complete, complete_fields))
run("SELECT duckvep_model_drop('grch38')")
run("SELECT * FROM duckvep_model_load('grch38',
 'SELECT seq_region, sequence_length FROM duckvep_bench_model.duckvep_sequence_regions ORDER BY seq_region',
 'SELECT * FROM duckvep_bench_model.duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
 'SELECT * FROM duckvep_bench_model.duckvep_exons ORDER BY transcript_index, exon_cdna_start',
 mature_mirna_query := 'SELECT * FROM duckvep_bench_model.duckvep_mature_mirna ORDER BY transcript_index, mature_mirna_start',
 peptide_edit_query := 'SELECT * FROM duckvep_bench_model.duckvep_peptide_edits ORDER BY transcript_index, protein_position',
 interval_feature_query := 'SELECT regulation_feature_index, seq_region, feature_start, feature_end, feature_kind FROM duckvep_bench_model.duckvep_regulation_features ORDER BY seq_region, feature_start, regulation_feature_index',
 transcript_coverage_complete := TRUE)")
receipts[[3L]] <- receipt("compact_regulation", compact, compact_fields)
run("SELECT duckvep_model_drop('grch38')")

# The historical replay is an immutable disagreement inventory, not a new VEP run.
replay_fields <- c("event_index", "task_index", "profile", "original_record_index", "input_id",
  "input_sha256", "model_sha256", "status", "passed", "lane", "record_index",
  "alt_index", "Feature", "field", "actual", "expected")
replay <- paste0("SELECT row_number() OVER (ORDER BY task_index, lane, record_index, alt_index, Feature, field, actual, expected)::UBIGINT AS event_index, ",
  paste(replay_fields[-1L], collapse = ", "), " FROM read_parquet(", q(source_file), ")")
replay_stats <- get(paste0("SELECT count(*) AS cells, count(DISTINCT task_index) AS tasks, ",
  "count(*) FILTER (WHERE passed IS DISTINCT FROM true) AS failed, ",
  "count(*) FILTER (WHERE actual IS DISTINCT FROM expected) AS disagreements FROM (", replay, ")"))
if (replay_stats$cells != 41942 || replay_stats$tasks != 10100 || replay_stats$failed != 0 ||
    replay_stats$disagreements != 41942) stop("historical retained inventory differs")
receipts[[4L]] <- receipt("retained_failures", replay, replay_fields)
results <- do.call(rbind, lapply(receipts, `[[`, "receipt"))
schemas <- do.call(rbind, lapply(receipts, `[[`, "schema"))
manifest <- data.frame(object = c("model", "canonical_metadata", "native_source_tree", "field_contract_source", "duckdb_version"),
  value = c(model_sha, metadata_sha, source_tree, fields_sha,
    get("SELECT version() AS version")$version),
  stringsAsFactors = FALSE)
check <- function(value, file) {
  if (record) {
    if (file.exists(file)) stop("baseline already exists: ", file)
    write.table(value, file, sep = "\t", quote = TRUE, row.names = FALSE, na = "-")
  } else {
    expected <- read.delim(file, colClasses = "character", check.names = FALSE, na.strings = "-")
    actual <- data.frame(lapply(value, as.character), check.names = FALSE)
    if (!identical(expected, actual)) stop("contract drift: ", file)
  }
}
check(manifest, file.path(output, "inputs.tsv"))
check(input_summary, file.path(output, "eligibility.tsv"))
check(schemas, file.path(output, "schemas.tsv"))
check(results, file.path(output, "receipts.tsv"))
dbDisconnect(con, shutdown = TRUE)
cat("verified: compact, compact_regulation, complete17, 41942 retained cells; model ", model_sha, "\n", sep = "")
