#!/usr/bin/env Rscript

# Generates everything the documentation site renders from the package itself,
# so the two cannot drift:
#
#   vignettes/*.Rmd  ->  site/src/pages/*.md
#   man/*.Rd         ->  site/src/pages/reference/*.md and its index
#
# The vignettes are the source of record for the narrative pages. Editing a
# page under site/src/pages is pointless: the next build overwrites it.
#
# Run from the package root, and before `npm run build`.

root <- normalizePath(".")
pages <- file.path(root, "site", "src", "pages")
refdir <- file.path(pages, "reference")
dir.create(refdir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

# ---- narrative pages ------------------------------------------------------

field <- function(yaml, key) {
  hit <- grep(paste0("^", key, ":"), yaml, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  v <- trimws(sub(paste0("^", key, ":"), "", hit[1]))
  gsub('^"|"$', "", v)
}

escape_yaml <- function(x) gsub('"', '\\\\"', x)

vignettes <- sort(list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$",
                             full.names = TRUE))
if (length(vignettes) == 0) stop("No vignettes found to build the site from.")

for (path in vignettes) {
  lines <- readLines(path, warn = FALSE)
  ends <- which(lines == "---")
  if (length(ends) < 2) stop("Vignette without a YAML header: ", path)
  yaml <- lines[(ends[1] + 1):(ends[2] - 1)]
  body <- lines[(ends[2] + 1):length(lines)]
  while (length(body) && !nzchar(body[1])) body <- body[-1]

  slug <- sub("[.]Rmd$", "", basename(path))
  title <- field(yaml, "title") %||% slug
  description <- field(yaml, "description") %||% ""
  order <- suppressWarnings(as.integer(field(yaml, "site_order")))
  wide <- identical(field(yaml, "site_wide"), "true")

  header <- c(
    "---",
    "layout: ../layouts/Base.astro",
    sprintf('title: "%s"', escape_yaml(title)),
    sprintf('description: "%s"', escape_yaml(description)),
    # The order travels with the page, so the navigation is read from the pages
    # that exist rather than from a file generated beside them.
    sprintf("site_order: %d", if (is.na(order)) 99L else order),
    if (wide) "wide: true",
    "---",
    "",
    paste0("# ", title),
    "",
    sprintf('<p class="lede">%s</p>', description),
    ""
  )
  writeLines(c(header, body), file.path(pages, paste0(slug, ".md")))
  message("page  ", slug)
}

# ---- function reference ---------------------------------------------------

rd_files <- sort(list.files(file.path(root, "man"), pattern = "[.]Rd$",
                            full.names = TRUE))

topic <- function(path) {
  rd <- tools::parse_Rd(path)
  tags <- vapply(rd, function(x) attr(x, "Rd_tag") %||% "", character(1))
  get1 <- function(tag) {
    i <- which(tags == tag)
    if (length(i) == 0) return(NA_character_)
    trimws(paste(unlist(rd[[i[1]]]), collapse = ""))
  }
  keywords <- unlist(lapply(which(tags == "\\keyword"), function(i)
    trimws(paste(unlist(rd[[i]]), collapse = ""))))
  list(name = get1("\\name"), title = get1("\\title"),
       internal = "internal" %in% keywords, rd = rd)
}

# Cross references inside generated help point at an installed tree. Those that
# name a topic in this package become links within the site; the rest keep
# their text and lose the link, since there is nothing here to point at.
relink <- function(html, known) {
  html <- gsub('<a href="[^"]*/([A-Za-z0-9._]+)\\.html"', '<a href="./\\1"', html)
  html <- gsub('<a href="[^"]*html/([A-Za-z0-9._]+)\\.html"', '<a href="./\\1"', html)
  m <- gregexpr('<a href="\\./([A-Za-z0-9._]+)"[^>]*>(.*?)</a>', html, perl = TRUE)
  regmatches(html, m) <- lapply(regmatches(html, m), function(links) {
    vapply(links, function(a) {
      target <- sub('.*<a href="\\./([A-Za-z0-9._]+)".*', "\\1", a)
      if (target %in% known) a else sub("<a [^>]*>|</a>", "", a)
    }, character(1))
  })
  html
}

topics <- lapply(rd_files, topic)
topics <- Filter(function(t) !t$internal && !is.na(t$name), topics)
known <- vapply(topics, function(t) t$name, character(1))

for (t in topics) {
  tmp <- tempfile(fileext = ".html")
  # Supplying a package name makes Rd2HTML look for an installed copy solely
  # to print its version. The site is generated from this source tree, so an
  # installed copy is neither needed nor desirable in a clean CI checkout.
  tools::Rd2HTML(t$rd, out = tmp, package = "")
  x <- readLines(tmp, warn = FALSE)
  i <- grep("<h2", x)[1]
  j <- grep("</main>", x)[1]
  if (is.na(i) || is.na(j)) next
  html <- paste(x[(i + 1):(j - 1)], collapse = "\n")
  html <- relink(html, known)
  # the generated page repeats the topic name in a footer link block
  html <- sub("<hr />\\s*<div style[^>]*>\\[Package.*$", "", html)

  writeLines(c(
    "---",
    "layout: ../../layouts/Base.astro",
    sprintf('title: "%s"', escape_yaml(t$name)),
    sprintf('description: "%s"', escape_yaml(t$title)),
    "wide: true",
    "---",
    "",
    paste0("# `", t$name, "`"),
    "",
    sprintf('<p class="lede">%s</p>', t$title),
    "",
    '<div class="rd">',
    html,
    "</div>"
  ), file.path(refdir, paste0(t$name, ".md")))
}
message("reference topics: ", length(topics))

titles <- vapply(topics, function(t) t$title, character(1))
ord <- order(known)
writeLines(c(
  "---",
  "layout: ../../layouts/Base.astro",
  'title: "Function reference"',
  'description: "Every exported function, its arguments and what it returns."',
  "wide: true",
  "---",
  "",
  "# Function reference",
  "",
  '<p class="lede">Every exported function, its arguments and what it returns. Generated from the package sources.</p>',
  "",
  "| Function | Purpose |",
  "|---|---|",
  sprintf("| [`%s`](./%s) | %s |", known[ord], known[ord], titles[ord])
), file.path(refdir, "index.md"))
message("done")
