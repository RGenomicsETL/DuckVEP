library(DBI)
library(duckdb)

local({
model <- Sys.getenv("DUCKVEP_GRCH38_MODEL")
dispositions <- Sys.getenv("DUCKVEP_REFERENCE_DISPOSITIONS")
if (!nzchar(model) || !nzchar(dispositions)) {
  stop("Set DUCKVEP_GRCH38_MODEL and DUCKVEP_REFERENCE_DISPOSITIONS to the pinned model and source-dispositions CSV")
}
extension <- Sys.getenv("DUCKVEP_EXTENSION", "build/release/duckvep.duckdb_extension")
con <- dbConnect(duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
quote <- function(x) as.character(dbQuoteString(con, x))
dbExecute(con, paste("LOAD", quote(normalizePath(extension))))
dbExecute(con, paste("ATTACH", quote(normalizePath(model)), "AS source (READ_ONLY)"))
name <- dbGetQuery(con, "SELECT transcript_index, transcript_id FROM source.duckvep_transcript_names WHERE transcript_id='ENST00000336066'")
stopifnot(nrow(name) == 1L)
index <- name$transcript_index[[1L]]
provenance <- read.csv(dispositions)
provenance <- provenance[provenance$transcript_index == index, ]
stopifnot(nrow(provenance) == 1L, provenance$transcript_stable_id == name$transcript_id,
  provenance$transcript_version == 8L, provenance$translation_stable_id == "ENSP00000477229",
  provenance$translation_version == 2L)
tx <- dbGetQuery(con, paste0("SELECT seq_region, strand, transcript_flags, cds_start, cds_end, codon_table, cast(cds_sequence AS VARCHAR) cds FROM source.duckvep_transcripts WHERE transcript_index=", index))
stopifnot(nrow(tx) == 1L, tx$codon_table == 1L, tx$strand == 1L,
  tx$seq_region == 1L, tx$cds_start == 3889431, substr(tx$cds, 1L, 3L) == "CTG",
  nchar(tx$cds) == provenance$cds_length, substr(tx$cds, 1L, 15L) == provenance$cds_first15,
  bitwAnd(tx$transcript_flags, 16L + 32L) == 0L)
region <- dbGetQuery(con, "SELECT * FROM source.duckvep_sequence_regions WHERE seq_region=1")
stopifnot(nrow(region) == 1L, region$name == "20", region$sequence_length == 64444167)
exons <- dbGetQuery(con, paste0("SELECT * FROM source.duckvep_exons WHERE transcript_index=", index,
  " ORDER BY exon_cdna_start"))
stopifnot(nrow(exons) == 7L, exons$exon_start[[1L]] <= tx$cds_start,
  exons$exon_end[[1L]] >= tx$cds_start)
queries <- c("SELECT 0::UINTEGER seq_region,sequence_length FROM source.duckvep_sequence_regions WHERE seq_region=1",
  paste0("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,* EXCLUDE(transcript_index,seq_region) FROM source.duckvep_transcripts WHERE transcript_index=", index),
  paste0("SELECT 0::UINTEGER transcript_index,* EXCLUDE(transcript_index) FROM source.duckvep_exons WHERE transcript_index=", index))
loaded <- dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('named_start', ",
  paste(quote(queries), collapse = ","), ")"))
stopifnot(identical(loaded$loaded, TRUE))
call <- paste0("SELECT 1::BIGINT event_index,0::UINTEGER seq_region,",
  "3889431::BIGINT AS position,'C' AS reference,'T' AS alternate,1::INTEGER AS alt_index,",
  "0::UINTEGER transcript_index,1::INTEGER sample_index,",
  "[1,0]::INTEGER[] alleles,[false,true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set")
result <- dbGetQuery(con, paste0("SELECT transcript_index, cds, protein, c.prediction_status,",
  " c.prediction_reason, c.haplotype_consequences FROM ",
  "(SELECT transcript_index, cast(cds AS VARCHAR) cds, protein, unnest(carrier_predictions) c ",
  "FROM duckvep_haplotypes(", quote(call), ",'named_start')) ORDER BY transcript_index"))
print(provenance[c("transcript_stable_id", "transcript_version", "translation_stable_id", "translation_version")])
print(region)
print(tx[c("seq_region", "strand", "transcript_flags", "cds_start", "cds_end", "codon_table")])
print(result[c("transcript_index", "prediction_status", "prediction_reason", "haplotype_consequences")])
print(substr(result$cds, 1L, 12L))
stopifnot(nrow(result) == 1L, result$prediction_status == "predicted",
  result$prediction_reason == "supported_domain", substr(result$cds, 1L, 3L) == "TTG",
  identical(unname(result$haplotype_consequences[[1L]]),
    c("start_retained_variant", "synonymous_variant")))
})
