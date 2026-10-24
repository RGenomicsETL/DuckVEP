# Evidence shown in README.Rmd, computed from the committed ledgers and benchmark receipts.
#
# Ledgers record the DuckHTS source revision that produced each run. DuckVEP's
# history was extracted with git filter-repo, so revisions are translated through
# design/duckhts-commit-map.txt before choosing, per corpus, the newest run whose
# revision is an ancestor of HEAD (the same rule as benchmarks/duckvep_conformance.Rmd).

readme_evidence <- function(root = ".") {
  path <- function(...) file.path(root, ...)

  map_lines <- readLines(path("design", "duckhts-commit-map.txt"))
  map_lines <- map_lines[!startsWith(map_lines, "#") & !startsWith(map_lines, "old ")]
  map_fields <- do.call(rbind, strsplit(map_lines, " +"))
  commit_map <- setNames(map_fields[, 2L], map_fields[, 1L])
  commit_map[commit_map == strrep("0", 40L)] <- NA_character_
  ancestry <- system2("git", c("-C", shQuote(root), "rev-list", "--topo-order", "HEAD"),
    stdout = TRUE)
  translate <- function(revision) unname(commit_map[revision])

  latest <- function(rows) {
    mapped <- translate(rows$source_revision)
    tested <- ancestry[ancestry %in% mapped]
    if (length(tested) > 0L) {
      return(rows[!is.na(mapped) & mapped == tested[[1L]], , drop = FALSE])
    }
    rows[rows$source_revision == rows$source_revision[[nrow(rows)]], , drop = FALSE]
  }

  history <- read.csv(path("test", "duckvep", "conformance", "data", "conformance_history.csv"),
    colClasses = c(source_revision = "character"))
  history <- history[history$engine == "duckvep" & history$stratum_kind == "all", ]
  corpora <- data.frame(
    corpus = c("final_dbsnp157_windows_seed29", "final_giab_grch38_seed71",
      "final_clinvar_coding_seed113", "final_clinvar_crosschrom_seed17",
      "final_grch37_cache_seed37", "plasmodium-falciparum-vep63-seed11663",
      "breakend_multichrom_seed31_isolated", "regulation_giab_chr21_seed1",
      "regulation_sv_chr21_seed17"),
    label = c("dbSNP 157 windows", "GIAB HG002 small variants", "ClinVar coding",
      "ClinVar, all chromosomes", "GRCh37 cache corpus", "*P. falciparum* (Ensembl Genomes 63)",
      "Paired breakends, multi-chromosome", "GIAB + regulatory and motif features",
      "Exact structural variants + regulation"),
    assembly = c("GRCh38", "GRCh38", "GRCh38", "GRCh38", "GRCh37", "GCA_000002765v3", "GRCh38",
      "GRCh38", "GRCh38"),
    stringsAsFactors = FALSE)
  conformance <- do.call(rbind, lapply(seq_len(nrow(corpora)), function(i) {
    rows <- history[history$corpus == corpora$corpus[[i]], , drop = FALSE]
    if (!nrow(rows)) return(NULL)
    run <- latest(rows)
    data.frame(corpora[i, c("label", "assembly")], oracle = paste("VEP", run$oracle_version[[1L]]),
      pairs = sum(run$n), exact = sum(run$exact_agree), unresolved = sum(run$unresolved),
      differences = sum(run$n) - sum(run$exact_agree) - sum(run$unresolved),
      stringsAsFactors = FALSE)
  }))

  hgvs <- read.csv(path("test", "duckvep", "conformance", "data", "hgvs_history.csv"),
    colClasses = c(source_revision = "character", extension_sha256 = "character"))
  hgvs <- hgvs[hgvs$extension_build_binding == "htslib_distclean_make_release" &
    grepl("^[0-9a-f]{64}$", hgvs$extension_sha256), ]
  hgvs_runs <- do.call(rbind, lapply(split(hgvs, hgvs$corpus), latest))
  hgvs_summary <- do.call(rbind, lapply(split(hgvs_runs, hgvs_runs$metric), function(x) {
    total <- sum(x$n)
    exact <- sum(x$n[x$comparison %in% c("match", "both_absent")])
    data.frame(field = toupper(sub("^hgvs", "HGVS", x$metric[[1L]])),
      pairs = total, exact = exact, discordant = total - exact,
      clinvar_pairs = sum(x$n[x$corpus == "clinvar_chr21_hgvs_seed113"]),
      stringsAsFactors = FALSE)
  }))
  hgvs_summary$field <- sub("^HGVSC$", "HGVSc", sub("^HGVSP$", "HGVSp", hgvs_summary$field))

  # The public SQL builder on full GIAB HG002: the newest measured revision that is
  # part of this checkout's history, every output mode and thread count.
  public <- read.csv(path("benchmarks", "data", "duckvep_throughput.csv"),
    colClasses = c(source_revision = "character"))
  public <- public[public$workload == "ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation", ]
  in_history <- vapply(public$source_revision, function(rev) {
    system2("git", c("merge-base", "--is-ancestor", rev, "HEAD"), stdout = FALSE, stderr = FALSE) == 0L
  }, logical(1))
  public <- public[in_history, ]
  stopifnot(nrow(public) > 0L)
  when <- vapply(public$source_revision, function(rev) {
    as.numeric(system2("git", c("show", "-s", "--format=%ct", rev), stdout = TRUE))
  }, numeric(1))
  public <- public[when == max(when), ]
  # Repeated measurements of one configuration: keep the last recorded.
  public <- public[!duplicated(public[c("output_mode", "threads")], fromLast = TRUE), ]

  list(conformance = conformance, hgvs = hgvs_summary, public = public)
}

# DuckDB box tables printed into a log, as data frames (first table whose header has `first`).
read_box <- function(path, first) {
  lines <- readLines(path, warn = FALSE)
  cells <- function(x) trimws(strsplit(sub("^│", "", sub("│$", "", x)), "│", fixed = TRUE)[[1L]])
  rows <- which(startsWith(lines, "│"))
  for (i in rows) {
    header <- cells(lines[[i]])
    if (!identical(header[[1L]], first)) next
    values <- list(); j <- i + 2L
    while (j <= length(lines) && !startsWith(lines[[j]], "└")) {
      if (startsWith(lines[[j]], "│")) values[[length(values) + 1L]] <- cells(lines[[j]])
      j <- j + 1L
    }
    frame <- as.data.frame(do.call(rbind, values), stringsAsFactors = FALSE)
    names(frame) <- header
    return(type.convert(frame, as.is = TRUE))
  }
  stop("no table starting with ", first, " in ", path)
}

# Newest benchmark directory (by its receipt's modification time) among `pattern` under `dir`.
newest_dir <- function(dir, pattern, receipt) {
  dirs <- list.files(dir, pattern = pattern, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, receipt))]
  stopifnot(length(dirs) > 0L)
  dirs[which.max(file.mtime(file.path(dirs, receipt)))]
}

speed_evidence <- function(root = ".") {
  path <- function(...) file.path(root, ...)
  vep_rs_dir <- newest_dir(path("benchmarks", "data", "vep_rs"), "_paired$", "timings.tsv")
  timings <- read.delim(file.path(vep_rs_dir, "timings.tsv"))
  adjudication <- file.path(vep_rs_dir, "adjudication.txt")
  tuples <- read_box(adjudication, "duckvep_tuples")
  verdicts <- read_box(adjudication, "verdict")
  vep_rs <- list(
    duckvep_s = median(timings$wall_s[timings$tool == "duckvep"]),
    veprs_s = median(timings$wall_s[timings$tool == "vep-rs"]),
    duckvep_rss = median(timings$maxrss_kib[timings$tool == "duckvep"]),
    veprs_rss = median(timings$maxrss_kib[timings$tool == "vep-rs"]),
    runs = sum(timings$tool == "duckvep"), tuples = tuples,
    vep_agrees_duckvep = sum(verdicts$n[verdicts$verdict == "VEP agrees with DuckVEP"]),
    disagreements = tuples$only_duckvep, revision = basename(vep_rs_dir))
  load <- read.csv(path("benchmarks", "data", "duckvep_model_load.csv"),
    colClasses = c(source_revision = "character"))
  load <- load[load$source_revision == load$source_revision[[nrow(load)]], ]
  model_load <- list(
    relation_1 = median(load$seconds[load$method == "relation_load" & load$threads == 1L]),
    relation_8 = median(load$seconds[load$method == "relation_load" & load$threads == 8L]),
    restore = median(load$seconds[load$method == "snapshot_restore"]),
    relation_rss = median(load$max_rss_kib[load$method == "relation_load" & load$threads == 1L]),
    restore_rss = median(load$max_rss_kib[load$method == "snapshot_restore"]))
  hap_dir <- newest_dir(path("benchmarks", "data", "haplotype_scale"), "^hg002_", "process.tsv")
  runs <- read.delim(file.path(hap_dir, "process.tsv"))
  ok <- runs$exit == 0L
  csq_s <- median(runs$wall_s[ok & grepl("^csq_", runs$label)])
  cli_s <- median(runs$wall_s[ok & grepl("^duckvep_F_", runs$label)])
  r_s <- median(runs$wall_s[ok & grepl("^duckvep_R_F_", runs$label)])
  haplotypes <- list(csq_s = csq_s, cli_s = cli_s, r_s = r_s, cli_speedup = csq_s / cli_s,
    r_speedup = csq_s / r_s, dir = basename(hap_dir))
  list(vep_rs = vep_rs, model_load = model_load, haplotypes = haplotypes)
}
