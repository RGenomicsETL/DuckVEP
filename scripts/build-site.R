#!/usr/bin/env Rscript
# Build the two-level project site into _site/:
#   index.html        the README (the extension landing page)
#   docs/*.html       design contract and evidence reports
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
landing_header <- normalizePath("scripts/site/landing-header.html", winslash = "/", mustWork = TRUE)
docs_header <- normalizePath("scripts/site/docs-header.html", winslash = "/", mustWork = TRUE)

metadata <- function(title, header) {
  c("---", paste0("title: ", title), "output:", "  html:", "    options:",
    "      toc: true", "    meta:",
    paste0("      css: [\"@default@1.14.69\", \"@article@1.14.69\", ",
      "\"@site@1.14.69\", \"", landing_css, "\"]"),
    paste0("      include_before: \"", header, "\""), "---")
}

# Documents published as site pages; every other repository link goes to GitHub.
documents <- c(
  errata = "ERRATA.md",
  design = "design/duckvep.md",
  conformance = "benchmarks/duckvep_conformance.md",
  throughput = "benchmarks/duckvep_throughput.md",
  fastvep = "benchmarks/benchmark_duckvep_fastvep.md",
  haplotypes = "benchmarks/duckvep_haplotypes.md",
  `haplotypes-indel` = "benchmarks/duckvep_haplotypes_indel.md",
  projection = "benchmarks/benchmark_duckvep_projection.md",
  `corpus-workflow` = "design/duckvep_corpus_workflow.md",
  `mane-grch37` = "design/duckvep_mane_grch37.md"
)
titles <- c(errata = "Compatibility and errata", design = "Design and implementation contract", conformance = "Conformance against Ensembl VEP",
  throughput = "Throughput", fastvep = "DuckVEP and FastVEP", haplotypes = "Haplotypes",
  `haplotypes-indel` = "Haplotypes with indels", projection = "Transcript projection",
  `corpus-workflow` = "Corpus workflow", `mane-grch37` = "MANE v1.5 mapped to GRCh37")

revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
stopifnot(length(revision) == 1L, grepl("^[0-9a-f]{40}$", revision))

# Point relative links at site pages when published here, otherwise at the source
# file on GitHub at this revision. `from` is the source file's directory.
page_ids <- function(path) {
  xml2::xml_attr(xml2::xml_find_all(xml2::read_html(path), "//*[@id]"), "id")
}

relink <- function(path, from, page_prefix) {
  document <- xml2::read_html(path)
  ids <- xml2::xml_attr(xml2::xml_find_all(document, "//*[@id]"), "id")
  for (node in xml2::xml_find_all(document, "//a[@href][not(ancestor::nav)] | //img[@src]")) {
    attribute <- if (xml2::xml_name(node) == "img") "src" else "href"
    href <- xml2::xml_attr(node, attribute)
    if (startsWith(href, "#")) {
      # litedown prefixes heading ids with "sec:"; GitHub-style links do not.
      anchor <- substring(href, 2L)
      if (!anchor %in% ids && paste0("sec:", anchor) %in% ids) {
        xml2::xml_set_attr(node, attribute, paste0("#sec:", anchor))
      }
      next
    }
    if (grepl("^([a-z]+:|//)", href)) {
      next
    }
    target <- sub("[?#].*$", "", href)
    fragment <- if (grepl("#", href, fixed = TRUE)) sub("^[^#]*", "", href) else ""
    source <- normalizePath(file.path(from, target), winslash = "/", mustWork = FALSE)
    relative <- sub(paste0("^", normalizePath(".", winslash = "/"), "/"), "", source)
    page <- names(documents)[documents == relative]
    if (length(page) == 1L) {
      if (nzchar(fragment)) {
        target_ids <- page_ids(file.path(site_root, "docs", paste0(page, ".html")))
        anchor <- substring(fragment, 2L)
        if (!anchor %in% target_ids && paste0("sec:", anchor) %in% target_ids) {
          fragment <- paste0("#sec:", anchor)
        }
      }
      xml2::xml_set_attr(node, attribute, paste0(page_prefix, page, ".html", fragment))
    } else if (attribute == "src" && file.exists(relative)) {
      copied <- file.path(site_root, "assets", relative)
      dir.create(dirname(copied), recursive = TRUE, showWarnings = FALSE)
      file.copy(relative, copied, overwrite = TRUE)
      xml2::xml_set_attr(node, attribute, paste0(if (page_prefix == "") "../" else "", "assets/", relative))
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
litedown::mark(text = c(metadata("DuckVEP", landing_header), readme), output = landing)

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
  litedown::mark(text = c(metadata(titles[[page]], docs_header), markdown), output = output)
}
index <- c("# Reports and design", "",
  sprintf("- [%s](%s.html)", titles[names(documents)], names(documents)))
litedown::mark(text = c(metadata("Reports and design", docs_header), index),
  output = file.path(docs_root, "index.html"))

# Relink only once every page exists, so fragments resolve against real ids.
relink(landing, ".", "docs/")
for (page in names(documents)) {
  relink(file.path(docs_root, paste0(page, ".html")), dirname(documents[[page]]), "")
}

message("Site pages: ", length(list.files(site_root, pattern = "\\.html$", recursive = TRUE)))
