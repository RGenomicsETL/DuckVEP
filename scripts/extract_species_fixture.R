#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: extract_species_fixture.R PF_MODEL TETRA_MODEL EXTENSION_COPY PF_FASTA TETRA_FASTA", call. = FALSE)
}
root <- normalizePath(".")
out <- file.path(root, "test/data/duckvep/species")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
con <- dbConnect(duckdb(shared_home = FALSE,
                        config = list(allow_unsigned_extensions = "true")),
                 dbdir = ":memory:")
on.exit(dbDisconnect(con, shutdown = TRUE))
quote <- function(x) as.character(dbQuoteString(con, x))
extension <- normalizePath(args[[3L]], mustWork = TRUE)
if (grepl("/build/release/", extension, fixed = TRUE)) stop("use an immutable extension copy")
dbExecute(con, paste("LOAD", quote(extension)))

extract <- function(name, region_sql, transcript_sql) {
  dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE fixture_regions AS", region_sql))
  dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE fixture_transcripts AS", transcript_sql))
  stopifnot(dbGetQuery(con, "SELECT count(*) n FROM fixture_regions")$n > 0L,
            dbGetQuery(con, "SELECT count(*) n FROM fixture_transcripts")$n > 0L)
  for (table in c("regions", "transcripts")) {
    file <- file.path(out, paste0(name, "_", table, ".parquet"))
    dbExecute(con, paste0("COPY fixture_", table, " TO ", quote(file),
                          " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  }
}
for (i in 1:2) {
  name <- c("plasmodium", "tetrahymena")[[i]]
  dbExecute(con, paste("ATTACH", quote(normalizePath(args[[i]], mustWork = TRUE)),
                       "AS source (READ_ONLY)"))
  ids <- if (i == 1L) "'PF3D7_MIT01400.1'" else "'EAR80522', 'EAR80553'"
  dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE picked AS SELECT * FROM source.model_transcripts WHERE transcript_stable_id IN (", ids, ")"))
  extract(name,
    "SELECT (row_number() OVER (ORDER BY seq_region) - 1)::UINTEGER AS seq_region, r.* EXCLUDE (seq_region) FROM source.model_regions r WHERE seq_region IN (SELECT seq_region FROM picked) ORDER BY seq_region",
    "SELECT (row_number() OVER (ORDER BY t.seq_region, t.transcript_start, t.transcript_index) - 1)::UINTEGER AS transcript_index, r.seq_region AS seq_region, t.* EXCLUDE (transcript_index, seq_region) FROM picked t JOIN fixture_regions r ON t.source_seq_region_id = r.source_seq_region_id ORDER BY seq_region, transcript_start")
  fasta <- normalizePath(args[[i + 3L]], mustWork = TRUE)
  regions <- dbGetQuery(con, "SELECT seq_region_name FROM fixture_regions ORDER BY seq_region")$seq_region_name
  excerpt <- file.path(out, paste0(name, ".fa"))
  stopifnot(system2("samtools", c("faidx", shQuote(fasta), regions), stdout = excerpt) == 0L,
            system2("samtools", c("faidx", shQuote(excerpt))) == 0L)
  dbExecute(con, "DETACH source")
}

# GRCh37's release-116 core excerpt is already pinned under test/data/duckvep.
dbExecute(con, "CREATE SCHEMA human_core")
base <- file.path(root, "test/data/duckvep/ensembl_core/grch37")
for (table in c("attrib_type", "coord_system", "seq_region", "seq_region_attrib",
                "gene", "transcript", "transcript_attrib", "translation",
                "translation_attrib", "exon", "exon_transcript")) {
  dbExecute(con, paste0("CREATE VIEW human_core.", table, " AS FROM read_parquet(",
                        quote(file.path(base, paste0(table, ".parquet"))), ")"))
}
dbExecute(con, paste0("CREATE VIEW human_reference AS FROM read_parquet(",
                      quote(file.path(base, "reference_chunks.parquet")), ")"))
dbExecute(con, "CREATE TEMP TABLE human_regions AS FROM query(duckvep_ensembl_regions_sql('human_core', 'human_reference', 'GRCh37'))")
dbExecute(con, "CREATE TEMP TABLE human_transcripts AS FROM query(duckvep_ensembl_transcripts_sql('human_core', 'human_reference', 'GRCh37'))")
dbExecute(con, "CREATE OR REPLACE TEMP TABLE picked AS SELECT * FROM human_transcripts WHERE codon_table = 1 ORDER BY transcript_index LIMIT 1")
extract("human_grch37",
  "SELECT (row_number() OVER (ORDER BY seq_region) - 1)::UINTEGER AS seq_region, r.* EXCLUDE (seq_region) FROM human_regions r WHERE seq_region IN (SELECT seq_region FROM picked) ORDER BY seq_region",
  "SELECT (row_number() OVER (ORDER BY t.seq_region, t.transcript_start, t.transcript_index) - 1)::UINTEGER AS transcript_index, r.seq_region AS seq_region, t.* EXCLUDE (transcript_index, seq_region) FROM picked t JOIN fixture_regions r ON t.source_seq_region_id = r.source_seq_region_id ORDER BY seq_region, transcript_start")
package_out <- file.path(root, "r/Rduckvep/inst/extdata/species")
dir.create(package_out, recursive = TRUE, showWarnings = FALSE)
files <- list.files(out, pattern = "[.](parquet|fa|fai)$", full.names = TRUE)
stopifnot(all(file.copy(files, package_out, overwrite = TRUE)))
print(dbGetQuery(con, "SELECT transcript_stable_id, seq_region, transcript_index, codon_table FROM fixture_transcripts"))
