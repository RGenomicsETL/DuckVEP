library(DBI)
library(duckdb)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: stage_gnomad_v4.1.R lane object md5 output.partial")
lane <- args[[1L]]
object <- args[[2L]]
source_md5 <- args[[3L]]
out <- args[[4L]]
if (!lane %in% c("genomes", "exomes", "sv") || !dir.exists(out)) stop("invalid staging request")
con <- dbConnect(duckdb(shared_home = FALSE))
execute <- function(sql) invisible(dbExecute(con, sql))
quote <- function(x) as.character(dbQuoteString(con, x))
execute("SET threads=2")
execute("SET memory_limit='2GB'")
model <- Sys.getenv("DUCKVEP_SCALE_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
pin <- read.delim("benchmarks/data/scale_contracts/inputs.tsv", colClasses = "character")
expected_model <- pin$value[pin$object == "model"]
observed_model <- substr(system2("sha256sum", model, stdout = TRUE), 1L, 64L)
if (length(expected_model) != 1L || expected_model != observed_model) stop("GRCh38 model receipt differs")
execute(paste0("ATTACH ", quote(normalizePath(model)), " AS model (READ_ONLY)"))
regions <- dbGetQuery(con, "SELECT name, seq_region FROM model.duckvep_sequence_regions")
region_codes <- setNames(as.integer(regions$seq_region), regions$name)
if (anyDuplicated(names(region_codes)) || anyNA(region_codes)) stop("invalid model region map")
execute(paste0("SET temp_directory=", quote(file.path(out, "spill"))))
execute("CREATE TEMP TABLE alleles (
 source_object VARCHAR, source_md5_base64 VARCHAR, record_index BIGINT, alt_index INTEGER,
 chrom VARCHAR, seq_region INTEGER, position BIGINT, id VARCHAR, reference VARCHAR,
 alternate VARCHAR, filter VARCHAR, status VARCHAR, end_position BIGINT, svtype VARCHAR,
 svlen VARCHAR, cipos VARCHAR, ciend VARCHAR, mateid VARCHAR, event VARCHAR)")
info_value <- function(info, key) {
  hit <- regexpr(paste0("(^|;)", key, "="), info)
  start <- as.integer(hit) + attr(hit, "match.length")
  value <- rep(NA_character_, length(info))
  found <- hit > 0L
  if (any(found)) value[found] <- sub(";.*$", "", substring(info[found], start[found]))
  value
}
sha256 <- function(path) substr(system2("sha256sum", path, stdout = TRUE), 1L, 64L)
record <- 0L
alts <- 0L
part <- 0L
parts <- list()
counts <- integer(5L)
flush <- function() {
  n <- dbGetQuery(con, "SELECT count(*) AS n FROM alleles")$n
  if (n == 0L) return(invisible(NULL))
  part <<- part + 1L
  path <- sprintf("part-%06d.parquet", part)
  execute(paste0("COPY alleles TO ", quote(file.path(out, path)),
    " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  parts[[part]] <<- data.frame(part = path, rows = n,
    bytes = file.size(file.path(out, path)), sha256 = sha256(file.path(out, path)))
  existing <- suppressWarnings(as.numeric(Sys.getenv("DUCKVEP_STAGED_BYTES", "0")))
  if (!is.finite(existing) || existing + sum(vapply(parts, function(x) x$bytes, numeric(1L))) > 25e9) {
    stop("25 GB staged-output limit reached")
  }
  space <- tail(system2("df", c("-PB1", out), stdout = TRUE), 1L)
  free <- as.numeric(strsplit(trimws(space), "[[:space:]]+")[[1L]][4L])
  if (!is.finite(free) || free < 40e9) stop("shared disk headroom below 40 GB")
  execute("DELETE FROM alleles")
  invisible(NULL)
}
input <- file("stdin", "r")
repeat {
  lines <- readLines(input, n = 5000L, warn = FALSE)
  if (length(lines) == 0L) break
  cols <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(cols) < 8L)) stop("truncated VCF record")
  field <- function(i) vapply(cols, `[[`, character(1L), i)
  chrom <- field(1L)
  position <- as.numeric(field(2L))
  id <- field(3L)
  reference <- field(4L)
  raw_alt <- field(5L)
  filter <- field(7L)
  info <- field(8L)
  alternates <- strsplit(raw_alt, ",", fixed = TRUE)
  size <- lengths(alternates)
  if (any(!is.finite(position) | position < 1 | size < 1L)) stop("invalid physical record")
  index <- rep(seq_along(lines), size)
  alt <- unlist(alternates, use.names = FALSE)
  seq_region <- unname(region_codes[sub("^chr", "", chrom)])
  ref <- reference[index]
  literal <- grepl("^[ACGTNacgtn]+$", ref) & grepl("^[ACGTNacgtn]+$", alt)
  status <- ifelse(literal & toupper(ref) != toupper(alt) & !is.na(seq_region[index]), "literal",
    ifelse(alt == "*", "star", ifelse(grepl("^<[^>]+>$", alt), "symbolic",
      ifelse(grepl("[", alt, fixed = TRUE) | grepl("]", alt, fixed = TRUE),
        "breakend", "unsupported"))))
  if (lane == "sv") {
    extra <- lapply(c("END", "SVTYPE", "SVLEN", "CIPOS", "CIEND", "MATEID", "EVENT"),
      function(key) info_value(info, key)[index])
  } else extra <- rep(list(rep(NA_character_, length(alt))), 7L)
  end <- suppressWarnings(as.numeric(extra[[1L]]))
  end[!is.finite(end) | end < 1] <- NA_real_
  data <- data.frame(source_object = object, source_md5_base64 = source_md5,
    record_index = record + index, alt_index = sequence(size), chrom = chrom[index],
    seq_region = seq_region[index], position = position[index], id = id[index],
    reference = ref, alternate = alt, filter = filter[index], status = status,
    end_position = end, svtype = extra[[2L]], svlen = extra[[3L]], cipos = extra[[4L]],
    ciend = extra[[5L]], mateid = extra[[6L]], event = extra[[7L]])
  dbWriteTable(con, "alleles", data, append = TRUE, temporary = TRUE)
  record <- record + length(lines)
  alts <- alts + length(alt)
  counts <- counts + tabulate(match(status, c("literal", "star", "symbolic", "breakend", "unsupported")), 5L)
  if (dbGetQuery(con, "SELECT count(*) AS n FROM alleles")$n >= 100000L) flush()
}
flush()
close(input)
if (record == 0L || sum(counts) != alts || sum(vapply(parts, function(x) as.numeric(x$rows), numeric(1L))) != alts) {
  stop("record/ALT conservation failed")
}
write.table(do.call(rbind, parts), file.path(out, "parts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(records = record, alts = alts, literal = counts[1L],
  star = counts[2L], symbolic = counts[3L], breakend = counts[4L], unsupported = counts[5L]),
  file.path(out, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
dbDisconnect(con, shutdown = TRUE)
cat("staged", record, "physical records and", alts, "ALTs in", part, "parts\n")
