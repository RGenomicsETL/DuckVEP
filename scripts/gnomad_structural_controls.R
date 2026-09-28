library(DBI)
library(duckdb)
library(Rsamtools)
library(GenomicRanges)
root <- "/root/duckvep/data/gnomad-v4.1"
out <- file.path(root, "structural-controls.parquet")
if (file.exists(out)) stop("structural control artifact already exists")
reference <- "/root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
index <- read.delim(paste0(reference, ".fai"), header = FALSE)
lengths <- setNames(index$V2, index$V1)
contigs <- c(as.character(1:22), "X", "Y")
events <- do.call(rbind, lapply(c("DEL", "DUP", "INV", "CNV", "BND"), function(kind) {
  n <- seq_len(2000L)
  chrom <- contigs[(n - 1L) %% length(contigs) + 1L]
  mate <- contigs[(n + 5L) %% length(contigs) + 1L]
  edge <- n %% 20L == 0L
  position <- ifelse(edge, ifelse(n %% 40L == 0L, lengths[chrom] - 2000, 2),
    10000 + (n * 7919) %% (lengths[chrom] - 100000))
  mate_position <- ifelse(edge, ifelse(n %% 40L == 0L, lengths[mate] - 3000, 100 + n %% 200),
    12000 + (n * 10301) %% (lengths[mate] - 100000))
  span <- ifelse(n %% 10L == 0L, 100000L, 100 + n %% 10000L)
  data.frame(kind, event_id = sprintf("%s-%04d", kind, n), chrom, position,
    mate, mate_position, end_position = ifelse(kind == "BND", NA_real_,
      pmin(lengths[chrom], position + span)), orientation = (n - 1L) %% 4L,
    uncertain = n %% 13L == 0L)
}))
expand <- rep(seq_len(nrow(events)), ifelse(events$kind == "BND", 2L, 1L))
rows <- events[expand, ]
rows$side <- sequence(ifelse(events$kind == "BND", 2L, 1L))
bnd <- rows$kind == "BND"
rows$mate_chrom <- ifelse(rows$side == 1L, rows$mate, rows$chrom)
rows$mate_pos <- ifelse(rows$side == 1L, rows$mate_position, rows$position)
rows$chrom[bnd & rows$side == 2L] <- rows$mate[bnd & rows$side == 2L]
rows$position[bnd & rows$side == 2L] <- rows$mate_position[bnd & rows$side == 2L]
fa <- FaFile(reference)
open(fa)
bases <- as.character(scanFa(fa, param = GRanges(rows$chrom, IRanges(rows$position, width = 1L))))
if (any(!bases %in% c("A", "C", "G", "T", "N"))) stop("non-DNA reference in control")
remote <- paste0(rows$mate_chrom, ":", rows$mate_pos)
orientation <- ifelse(rows$side == 2L, (rows$orientation + 2L) %% 4L, rows$orientation)
alt <- ifelse(orientation == 0L, paste0(bases, "[", remote, "["),
  ifelse(orientation == 1L, paste0(bases, "]", remote, "]"),
    ifelse(orientation == 2L, paste0("]", remote, "]", bases),
      paste0("[", remote, "[", bases))))
alt[!bnd] <- paste0("<", rows$kind[!bnd], ">")
data <- data.frame(source_object = "synthetic:structural-controls-v1",
  source_md5_base64 = NA_character_, record_index = seq_len(nrow(rows)),
  alt_index = 1L, chrom = rows$chrom,
  seq_region = NA_integer_, position = rows$position,
  id = paste0(rows$event_id, ifelse(bnd, paste0("-", rows$side), "")),
  reference = bases, alternate = alt, filter = "PASS",
  status = ifelse(bnd, "breakend", "symbolic"),
  end_position = rows$end_position, svtype = rows$kind,
  svlen = ifelse(bnd, NA_character_, as.character(rows$end_position - rows$position)),
  cipos = ifelse(rows$uncertain, "-5,5", NA_character_),
  ciend = ifelse(rows$uncertain & !bnd, "-5,5", NA_character_),
  mateid = ifelse(bnd, paste0(rows$event_id, "-", 3L - rows$side), NA_character_),
  event = rows$event_id)
con <- dbConnect(duckdb(shared_home = FALSE))
dbExecute(con, "ATTACH '/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb' AS model (READ_ONLY)")
regions <- dbGetQuery(con, "SELECT name, seq_region FROM model.duckvep_sequence_regions")
data$seq_region <- unname(setNames(as.integer(regions$seq_region), regions$name)[data$chrom])
if (anyNA(data$seq_region)) stop("control contig is absent from model")
dbWriteTable(con, "controls", data, temporary = TRUE)
dbExecute(con, paste0("COPY controls TO ", as.character(dbQuoteString(con, out)),
  " (FORMAT PARQUET, COMPRESSION ZSTD)"))
checks <- dbGetQuery(con, "SELECT count(*) AS physical_records, count(DISTINCT event) AS events,
 count(*) FILTER (WHERE status = 'breakend') AS breakends,
 count(DISTINCT event) FILTER (WHERE svtype = 'BND') AS bnd_events FROM controls")
if (!identical(as.integer(unlist(checks)), c(12000L, 10000L, 4000L, 2000L))) stop("control counts differ")
dbDisconnect(con, shutdown = TRUE)
close(fa)
write.table(checks, file.path(root, "structural-controls-counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(out, "\n")
