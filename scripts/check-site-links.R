#!/usr/bin/env Rscript
# Every local link and image on the landing and report pages, and every navigation
# link on the package pages, must resolve inside _site/ (including #anchors).

site_root <- normalizePath("_site", mustWork = TRUE)
pages <- list.files(site_root, pattern = "\\.html$", recursive = TRUE)
required <- c("index.html", "Rduckvep/index.html", "docs/index.html", "docs/conformance.html")
missing <- setdiff(required, pages)
if (length(missing)) stop("Missing site pages: ", paste(missing, collapse = ", "), call. = FALSE)

checked <- 0L
for (page in pages) {
  document <- xml2::read_html(file.path(site_root, page))
  xpath <- if (startsWith(page, "Rduckvep/")) "//nav//a[@href]" else "//a[@href] | //img[@src]"
  for (node in xml2::xml_find_all(document, xpath)) {
    href <- xml2::xml_attr(node, if (xml2::xml_name(node) == "img") "src" else "href")
    if (grepl("^([a-z]+:|//)", href)) next
    path <- sub("[?#].*$", "", href)
    target <- if (nzchar(path)) file.path(site_root, dirname(page), path) else file.path(site_root, page)
    target <- utils::URLdecode(target)
    if (dir.exists(target)) target <- file.path(target, "index.html")
    if (!file.exists(target)) stop(page, ": broken link ", href, call. = FALSE)
    if (grepl("#", href, fixed = TRUE) && grepl("\\.html$", target)) {
      fragment <- utils::URLdecode(sub("^.*#", "", href))
      ids <- xml2::xml_attr(xml2::xml_find_all(xml2::read_html(target), "//*[@id]"), "id")
      if (nzchar(fragment) && !fragment %in% ids) stop(page, ": missing anchor ", href, call. = FALSE)
    }
    checked <- checked + 1L
  }
}
cat("Checked", checked, "local links across", length(pages), "HTML pages\n")
