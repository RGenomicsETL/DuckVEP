#!/usr/bin/env Rscript
# Build the two-level project site into _site/:
#   index.html        the README (the extension landing page)
#   docs/*.html       function reference, contracts and evidence reports
#   Rduckvep/         the R package's pkgdown site

for (package in c("pkgdown", "litedown", "xml2")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Building the site requires ", package, call. = FALSE)
  }
}

repository <- "https://github.com/RGenomicsETL/DuckVEP"
package_dir <- "r/Rduckvep"

unlink("_site", recursive = TRUE)
dir.create("_site")
file.create(file.path("_site", ".nojekyll"))
site_root <- normalizePath("_site", mustWork = TRUE)

# pkgdown refuses a non-pkgdown destination, so build into an empty directory.
destination <- file.path(site_root, "Rduckvep")
dir.create(destination)
pkgdown::build_site(pkg = package_dir, new_process = FALSE, install = FALSE,
  preview = FALSE, override = list(destination = destination))

landing_css <- normalizePath("scripts/site/landing.css", winslash = "/", mustWork = TRUE)
dir.create(file.path(site_root, "assets"))
stopifnot(file.copy(landing_css, file.path(site_root, "assets", "landing.css")))
landing_header <- normalizePath("scripts/site/landing-header.html", winslash = "/", mustWork = TRUE)
docs_header <- normalizePath("scripts/site/docs-header.html", winslash = "/", mustWork = TRUE)

metadata <- function(title, header, css, toc = TRUE) {
  c("---", paste0("title: ", title), "output:", "  html:", "    options:",
    paste0("      toc: ", tolower(toc)), "      embed_resources: false", "    meta:",
    paste0("      css: [\"@default@1.14.69\", \"@article@1.14.69\", ",
      "\"@site@1.14.69\", \"", css, "\"]"),
    paste0("      include_before: \"", header, "\""), "---")
}

# Documents published as site pages; every other repository link goes to GitHub.
documents <- c(
  functions = "docs/functions.md",
  errata = "ERRATA.md",
  `v2-host` = "docs/v2-host.md",
  design = "design/duckvep.md",
  `haplotype-contract` = "design/duckvep_haplotype_contract.md",
  `structural-hgvs` = "docs/structural-identity-hgvs.md",
  `corpus-workflow` = "design/duckvep_corpus_workflow.md",
  `mane-grch37` = "design/duckvep_mane_grch37.md",
  conformance = "benchmarks/duckvep_conformance.md",
  haplotypes = "benchmarks/duckvep_haplotypes.md",
  `haplotypes-indel` = "benchmarks/duckvep_haplotypes_indel.md",
  `haplotype-scale` = "benchmarks/data/haplotype_scale/README.md",
  projection = "benchmarks/benchmark_duckvep_projection.md",
  throughput = "benchmarks/duckvep_throughput.md",
  throughput_evidence = "benchmarks/duckvep_throughput_evidence.md",
  vep_rs = "benchmarks/benchmark_duckvep_vep_rs.md",
  `scale-runner` = "docs/scale-runner.md"
)
titles <- c(functions = "SQL function reference", errata = "Compatibility and errata",
  `v2-host` = "DuckDB hosts and release profile", design = "Design and implementation contract",
  `haplotype-contract` = "The duckvep-coding haplotype contract",
  `structural-hgvs` = "Breakend identity, fusion evidence and structural HGVS",
  `corpus-workflow` = "Corpus workflow", `mane-grch37` = "MANE v1.5 mapped to GRCh37",
  conformance = "Conformance against Ensembl VEP", haplotypes = "Haplotypes",
  `haplotypes-indel` = "Haplotypes with indels", `haplotype-scale` = "Haplotype scale qualification",
  projection = "Transcript projection", throughput = "Throughput",
  throughput_evidence = "Throughput methods and receipts", vep_rs = "DuckVEP and vep-rs",
  `scale-runner` = "Scale runner")
# Sections of the docs index; every published document belongs to exactly one.
sections <- list(
  Reference = c("functions", "errata", "v2-host"),
  Contracts = c("design", "haplotype-contract", "structural-hgvs", "corpus-workflow", "mane-grch37"),
  Evidence = c("conformance", "haplotypes", "haplotypes-indel", "haplotype-scale", "projection"),
  Performance = c("throughput", "throughput_evidence", "vep_rs", "scale-runner"))
stopifnot(setequal(names(titles), names(documents)),
  setequal(unlist(sections), names(documents)), !anyDuplicated(unlist(sections)))

published <- c(README.md = "index.html",
  setNames(paste0("docs/", names(documents), ".html"), documents))

revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
stopifnot(length(revision) == 1L, grepl("^[0-9a-f]{40}$", revision))

# Point relative links at site pages when published here, otherwise at the source
# file on GitHub at this revision. `from` is the source file's directory.
page_ids <- function(path) {
  xml2::xml_attr(xml2::xml_find_all(xml2::read_html(path), "//*[@id]"), "id")
}

# litedown prefixes heading ids with "sec:" and may fold punctuation such as "_"
# into "-"; GitHub-style anchors (`#duckvep_coding_calls`) keep the heading text.
site_anchor <- function(anchor, ids) {
  candidates <- c(anchor, paste0("sec:", anchor),
    paste0("sec:", gsub("^-+|-+$", "", gsub("[^[:alnum:]]+", "-", tolower(anchor)))))
  found <- candidates[candidates %in% ids]
  if (length(found) > 0L) found[[1L]] else anchor
}

relink <- function(path, from, site_prefix) {
  document <- xml2::read_html(path)
  ids <- xml2::xml_attr(xml2::xml_find_all(document, "//*[@id]"), "id")
  for (node in xml2::xml_find_all(document, "//a[@href][not(ancestor::nav)] | //img[@src]")) {
    attribute <- if (xml2::xml_name(node) == "img") "src" else "href"
    href <- xml2::xml_attr(node, attribute)
    if (startsWith(href, "#")) {
      xml2::xml_set_attr(node, attribute, paste0("#", site_anchor(substring(href, 2L), ids)))
      next
    }
    if (grepl("^([a-z]+:|//)", href)) {
      next
    }
    target <- sub("[?#].*$", "", href)
    fragment <- if (grepl("#", href, fixed = TRUE)) sub("^[^#]*", "", href) else ""
    source <- normalizePath(file.path(from, target), winslash = "/", mustWork = FALSE)
    relative <- sub(paste0("^", normalizePath(".", winslash = "/"), "/"), "", source)
    page <- unname(published[relative])
    if (!is.na(page)) {
      if (nzchar(fragment)) {
        target_ids <- page_ids(file.path(site_root, page))
        fragment <- paste0("#", site_anchor(substring(fragment, 2L), target_ids))
      }
      xml2::xml_set_attr(node, attribute, paste0(site_prefix, page, fragment))
    } else if (attribute == "src" && file.exists(relative)) {
      copied <- file.path(site_root, "assets", relative)
      dir.create(dirname(copied), recursive = TRUE, showWarnings = FALSE)
      file.copy(relative, copied, overwrite = TRUE)
      xml2::xml_set_attr(node, attribute, paste0(site_prefix, "assets/", relative))
    } else if (file.exists(relative) || dir.exists(relative)) {
      kind <- if (dir.exists(relative)) "tree" else "blob"
      xml2::xml_set_attr(node, attribute,
        paste0(repository, "/", kind, "/", revision, "/", relative, fragment))
    } else {
      stop(path, ": link has no source file: ", href, call. = FALSE)
    }
  }
  xml2::write_html(document, path)
}

readme <- readLines("README.md", warn = FALSE, encoding = "UTF-8")
readme <- readme[readme != "# DuckVEP"]
landing <- file.path(site_root, "index.html")
litedown::mark(text = c(metadata("DuckVEP", landing_header, "assets/landing.css", toc = FALSE), readme), output = landing)

docs_root <- file.path(site_root, "docs")
dir.create(docs_root)
for (page in names(documents)) {
  source <- documents[[page]]
  markdown <- readLines(source, warn = FALSE, encoding = "UTF-8")
  # The page title comes from the metadata; drop the document's own first H1
  # (ATX `# Title` or setext `Title` over `====`).
  atx <- which(startsWith(markdown, "# "))
  setext <- which(grepl("^=+$", markdown)) - 1L
  first <- suppressWarnings(min(c(atx, setext[setext >= 1L])))
  if (is.finite(first)) {
    markdown <- markdown[-c(first, if (first %in% setext) first + 1L)]
  }
  output <- file.path(docs_root, paste0(page, ".html"))
  # The function reference carries its own index table; a 60-entry TOC would bury it.
  litedown::mark(text = c(metadata(titles[[page]], docs_header, "../assets/landing.css",
    toc = page != "functions"), markdown), output = output)
}
index <- unlist(lapply(names(sections), function(section) {
  pages <- sections[[section]]
  c(paste("##", section), "", sprintf("- [%s](%s.html)", titles[pages], pages), "")
}))
litedown::mark(text = c(metadata("Documentation", docs_header, "../assets/landing.css", toc = FALSE),
  sprintf("Reference, contracts and evidence for DuckVEP, built from revision [`%s`](%s/tree/%s).",
    substr(revision, 1L, 10L), repository, revision),
  "The [R package reference](../Rduckvep/) is published separately.", "", index),
  output = file.path(docs_root, "index.html"))

# Relink only once every page exists, so fragments resolve against real ids.
relink(landing, ".", "")
for (page in names(documents)) {
  relink(file.path(docs_root, paste0(page, ".html")), dirname(documents[[page]]), "../")
}

message("Site pages: ", length(list.files(site_root, pattern = "\\.html$", recursive = TRUE)))
