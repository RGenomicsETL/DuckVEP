library(DBI)
library(duckdb)
root <- "/root/duckvep/data/gnomad-v4.1"
manifest <- read.delim("benchmarks/data/gnomad-v4.1-manifest.tsv", header = FALSE,
  col.names = c("object", "bytes", "generation", "md5", "crc32c"),
  colClasses = "character")
objects <- manifest[grepl("(\\.vcf\\.bgz|\\.vcf\\.gz)$", manifest$object), ]
if (nrow(objects) != 49L || anyDuplicated(objects$object)) stop("frozen manifest does not have 49 VCF objects")
sha256 <- function(path) substr(system2("sha256sum", path, stdout = TRUE), 1L, 64L)
model <- Sys.getenv("DUCKVEP_SCALE_MODEL", "/root/duckvep/data/models/homo_sapiens_116_GRCh38_final.duckdb")
model_sha <- sha256(model)
pin <- read.delim("benchmarks/data/scale_contracts/inputs.tsv", colClasses = "character")
if (!identical(pin$value[pin$object == "model"], model_sha)) stop("model receipt differs")
con <- dbConnect(duckdb(shared_home = FALSE))
dbExecute(con, paste0("ATTACH ", as.character(dbQuoteString(con, normalizePath(model))),
  " AS model (READ_ONLY)"))
q <- function(x) as.character(dbQuoteString(con, x))
shards <- list.dirs(root, recursive = FALSE, full.names = TRUE)
shards <- shards[file.exists(file.path(shards, "source.tsv"))]
receipts <- lapply(shards, function(shard) {
  source <- read.delim(file.path(shard, "source.tsv"), colClasses = "character")
  if (nrow(source) != 1L) stop("invalid source receipt in ", shard)
  ref <- objects[objects$object == source$object, ]
  if (nrow(ref) != 1L || ref$bytes != source$bytes_read ||
      ref$generation != source$generation || ref$md5 != source$md5_base64 ||
      ref$crc32c != source$crc32c_base64) stop("source object not pinned: ", shard)
  counts <- read.delim(file.path(shard, "counts.tsv"))
  parts <- read.delim(file.path(shard, "parts.tsv"), colClasses = "character")
  if (nrow(counts) != 1L || !nrow(parts) || anyDuplicated(parts$part) ||
      any(grepl("/|\\.\\.", parts$part))) stop("invalid staging counts or paths")
  paths <- file.path(shard, parts$part)
  if (any(!file.exists(paths)) || any(vapply(paths, sha256, character(1L)) != parts$sha256) ||
      any(as.numeric(file.size(paths)) != as.numeric(parts$bytes))) stop("Parquet checksum mismatch")
  actual <- dbGetQuery(con, paste0("SELECT count(*) AS alts, count(DISTINCT record_index) AS records,
    count(*) FILTER (WHERE status='literal') AS literal,
    count(*) FILTER (WHERE status='star') AS star,
    count(*) FILTER (WHERE status='symbolic') AS symbolic,
    count(*) FILTER (WHERE status='breakend') AS breakend,
    count(*) FILTER (WHERE status='unsupported') AS unsupported FROM read_parquet([",
    paste(vapply(paths, q, character(1L)), collapse = ","), "])"))
  if (any(unlist(actual) != unlist(counts[c("alts", "records", "literal", "star", "symbolic", "breakend", "unsupported")])) ||
      sum(as.numeric(parts$rows)) != counts$alts) stop("record/ALT conservation mismatch")
  mapping <- dbGetQuery(con, paste0("SELECT count(*) AS n FROM read_parquet([",
    paste(vapply(paths, q, character(1L)), collapse = ","),
    "]) a LEFT JOIN model.duckvep_sequence_regions r ON r.name = regexp_replace(a.chrom, '^chr', '')
    WHERE a.seq_region IS DISTINCT FROM r.seq_region OR r.seq_region IS NULL"))$n
  if (mapping != 0) stop("Parquet model region index differs from pinned model")
  filter_counts <- dbGetQuery(con, paste0("SELECT filter, count(*) AS alts FROM read_parquet([",
    paste(vapply(paths, q, character(1L)), collapse = ","),
    "]) GROUP BY filter ORDER BY filter"))
  if (sum(filter_counts$alts) != counts$alts) stop("FILTER conservation mismatch")
  filter_file <- file.path(shard, "filters.tsv")
  write.table(filter_counts, filter_file, sep = "\t", row.names = FALSE, quote = FALSE)
  list(receipt = data.frame(object = source$object, generation = source$generation,
    bytes_read = source$bytes_read, source_md5_base64 = source$md5_base64,
    model_sha256 = model_sha,
    records = counts$records, alts = counts$alts, literal = counts$literal,
    star = counts$star, symbolic = counts$symbolic, breakend = counts$breakend,
    unsupported = counts$unsupported, parquet_bytes = sum(file.size(paths)),
    parts_sha256 = sha256(file.path(shard, "parts.tsv")),
    filters_sha256 = sha256(filter_file)),
    parts = data.frame(object = source$object, parts))
})
if (length(receipts) == 0L) stop("no completed shards")
part_inventory <- do.call(rbind, lapply(receipts, `[[`, "parts"))
receipts <- do.call(rbind, lapply(receipts, `[[`, "receipt"))
receipts <- receipts[order(receipts$object), ]
part_inventory <- part_inventory[order(part_inventory$object, part_inventory$part), ]
part_inventory$rows <- format(as.numeric(part_inventory$rows), scientific = FALSE, trim = TRUE)
part_inventory$bytes <- format(as.numeric(part_inventory$bytes), scientific = FALSE, trim = TRUE)
receipts$parquet_bytes <- format(as.numeric(receipts$parquet_bytes), scientific = FALSE, trim = TRUE)
if (anyDuplicated(receipts$object)) stop("duplicate completed source")
missing <- objects[!objects$object %in% receipts$object, c("object", "bytes", "generation", "md5")]
write.table(receipts, "benchmarks/data/gnomad-v4.1-staging.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(part_inventory, "benchmarks/data/gnomad-v4.1-parts.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(missing, "benchmarks/data/gnomad-v4.1-remaining.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
geometry <- file.path(root, "sv-sites", "geometry.parquet")
if (file.exists(geometry)) {
  geometry_counts <- dbGetQuery(con, paste0("SELECT geometry_status, count(*) AS alts
    FROM read_parquet(", q(geometry), ") GROUP BY geometry_status ORDER BY geometry_status"))
  sv <- receipts[grepl("/genome_sv/", receipts$object), ]
  if (nrow(sv) != 1L || sum(geometry_counts$alts) != sv$alts) stop("invalid SV geometry")
  geometry_counts$parquet_sha256 <- sha256(geometry)
  write.table(geometry_counts, "benchmarks/data/gnomad-v4.1-sv-geometry.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)
}
panel_dirs <- list.dirs(file.path(root, "panels"), recursive = FALSE, full.names = TRUE)
if (length(panel_dirs) > 0L) {
  source_counts <- vapply(panel_dirs, function(dir) {
    path <- file.path(dir, "sources.tsv")
    if (!file.exists(path)) return(0L)
    nrow(read.delim(path))
  }, integer(1L))
  selected <- unlist(lapply(c("genomes-v1-", "exomes-v1-"), function(prefix) {
    choices <- which(startsWith(basename(panel_dirs), prefix))
    if (length(choices) == 0L) return(character())
    panel_dirs[choices[which.max(source_counts[choices])]]
  }), use.names = FALSE)
  if (length(selected) > 0L) {
    status <- lapply(selected, function(dir) {
      source_sha <- sha256(file.path(dir, "sources.tsv"))
      if (startsWith(basename(dir), "genomes")) {
        rows <- read.delim(file.path(dir, "panels.tsv"), colClasses = "character")
        for (j in which(rows$status == "complete")) {
          panel <- file.path(dir, paste0("panel-", rows$quota[j], ".parquet"))
          if (!file.exists(panel) || sha256(panel) != rows$sha256[j]) {
            stop("genome panel checksum differs")
          }
        }
        data.frame(lane = "genomes", required = rows$quota, actual = rows$rows,
          status = rows$status, sha256 = rows$sha256, sources_sha256 = source_sha)
      } else {
        rows <- read.delim(file.path(dir, "panel.tsv"), colClasses = "character")
        if (rows$status == "complete") {
          panel <- file.path(dir, "exome-2m.parquet")
          if (!file.exists(panel) || sha256(panel) != rows$sha256) {
            stop("exome panel checksum differs")
          }
        }
        if (file.exists(file.path(dir, "availability.tsv"))) {
          availability <- read.delim(file.path(dir, "availability.tsv"))
          write.table(availability, "benchmarks/data/gnomad-v4.1-exome-availability.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
          unlink("benchmarks/data/gnomad-v4.1-exome-bounds.tsv")
        } else {
          bounds <- read.delim(file.path(dir, "bounds.tsv"))
          write.table(bounds, "benchmarks/data/gnomad-v4.1-exome-bounds.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
          unlink("benchmarks/data/gnomad-v4.1-exome-availability.tsv")
        }
        data.frame(lane = "exomes", required = "2000000", actual = rows$rows,
          status = rows$status, sha256 = rows$sha256, sources_sha256 = source_sha)
      }
    })
    write.table(do.call(rbind, status), "benchmarks/data/gnomad-v4.1-panels.tsv",
      sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
control <- file.path(root, "structural-controls.parquet")
if (file.exists(control)) {
  summary <- dbGetQuery(con, paste0("SELECT svtype, count(DISTINCT event) AS events,
    count(*) AS physical_records, count(*) FILTER (WHERE cipos IS NOT NULL) AS imprecise
    FROM read_parquet(", q(control), ") GROUP BY svtype ORDER BY svtype"))
  if (!identical(as.integer(summary$events), rep(2000L, 5L)) ||
      sum(summary$physical_records) != 12000L) stop("invalid structural controls")
  mapping <- dbGetQuery(con, paste0("SELECT count(*) AS n FROM read_parquet(", q(control),
    ") a LEFT JOIN model.duckvep_sequence_regions r ON r.name = a.chrom
    WHERE a.seq_region IS DISTINCT FROM r.seq_region OR r.seq_region IS NULL"))$n
  if (mapping != 0) stop("structural control region index differs from model")
  pairs <- dbGetQuery(con, paste0("SELECT count(*) AS n FROM read_parquet(", q(control), ") a
    WHERE a.status = 'breakend' AND NOT EXISTS (
      SELECT 1 FROM read_parquet(", q(control), ") b
      WHERE b.id = a.mateid AND b.mateid = a.id AND b.event = a.event
        AND b.record_index <> a.record_index)"))$n
  if (pairs != 0) stop("structural breakends are not paired")
  orientations <- dbGetQuery(con, paste0("SELECT CASE
    WHEN starts_with(alternate, '[') THEN 'left_open'
    WHEN starts_with(alternate, ']') THEN 'left_close'
    WHEN contains(alternate, '[') THEN 'right_open'
    ELSE 'right_close' END AS orientation, count(*) AS events
    FROM read_parquet(", q(control), ") WHERE status = 'breakend' AND ends_with(id, '-1')
    GROUP BY orientation ORDER BY orientation"))
  if (nrow(orientations) != 4L || any(orientations$events != 500L)) {
    stop("structural control orientations differ")
  }
  write.table(orientations, "benchmarks/data/gnomad-v4.1-bnd-orientations.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)
  summary$parquet_sha256 <- sha256(control)
  summary$reference_sha256 <- sha256("/root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
  write.table(summary, "benchmarks/data/gnomad-v4.1-controls.tsv", sep = "\t",
    row.names = FALSE, quote = FALSE)
}
dbDisconnect(con, shutdown = TRUE)
cat(nrow(receipts), "completed objects;", nrow(missing), "remaining\n")
