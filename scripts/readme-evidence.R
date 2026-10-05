# Evidence shown in README.Rmd, computed from the committed ledgers.
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

  throughput <- read.csv(path("benchmarks", "data", "duckvep_throughput.csv"),
    colClasses = c(source_revision = "character"))
  throughput <- throughput[throughput$workload == "ensembl116_grch38_giab_sites_hash40" &
    throughput$output_mode == "rich" & throughput$threads == 1L, ]
  throughput <- throughput[order(throughput$run_date), ]
  throughput <- throughput[nrow(throughput), ]

  # The public SQL builder on full GIAB, one core, compact: the newest measured
  # revision that is part of this checkout's history.
  public <- read.csv(path("benchmarks", "data", "duckvep_throughput.csv"),
    colClasses = c(source_revision = "character"))
  public <- public[public$workload == "ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation" &
    public$output_mode == "compact" & public$threads == 1L, ]
  in_history <- vapply(public$source_revision, function(rev) {
    system2("git", c("merge-base", "--is-ancestor", rev, "HEAD"), stdout = FALSE, stderr = FALSE) == 0L
  }, logical(1))
  public <- public[in_history, ]
  stopifnot(nrow(public) > 0L)
  when <- vapply(public$source_revision, function(rev) {
    as.numeric(system2("git", c("show", "-s", "--format=%ct", rev), stdout = TRUE))
  }, numeric(1))
  public <- public[when == max(when), ]
  public <- public[nrow(public), ]

  haplotype_runs <- read.delim(path("benchmarks", "data", "haplotype_scale", "perf34_slice1", "process.tsv"))
  successful <- haplotype_runs$exit == 0L
  csq_s <- median(haplotype_runs$wall_s[successful & grepl("^csq_", haplotype_runs$label)])
  cli_s <- median(haplotype_runs$wall_s[successful & grepl("^duckvep_F_", haplotype_runs$label)])
  r_s <- median(haplotype_runs$wall_s[successful & grepl("^duckvep_R_F_", haplotype_runs$label)])
  haplotype <- list(csq_s = csq_s, cli_s = cli_s, cli_speedup = csq_s / cli_s,
    r_s = r_s, r_speedup = csq_s / r_s)

  list(conformance = conformance, hgvs = hgvs_summary,
    throughput = throughput, haplotype = haplotype, public = public,
    giab_alt_alleles = 4096123L)
}
