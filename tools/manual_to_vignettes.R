#!/usr/bin/env Rscript

# Writes the narrative vignettes from the manual, which is the source of record
# for them:
#
#   knowledge/MANUAL.md  the problem and the response  ->  why.Rmd
#                        how the estimator works      ->  how-it-works.Rmd
#                        the input data               ->  input-format.Rmd
#                        the outputs                  ->  outputs.Rmd
#                        the mathematics, the conditions and the sources
#                                                     ->  methods.Rmd
#
# vignettes/tutorial.Rmd and vignettes/simulation.Rmd are written by hand and
# are never touched here: they teach a workflow rather than describing the
# method, and have no counterpart in the manual.
#
# The manual sits outside the repository, so this runs when the manual changes
# and the vignettes it produces are committed. `site/sync.R` then builds the
# site from the committed vignettes, which is what continuous integration runs.
#
# Usage, from the package root:
#   Rscript tools/manual_to_vignettes.R [path/to/MANUAL.md]

args <- commandArgs(trailingOnly = TRUE)
manual <- if (length(args) > 0) args[[1]] else file.path("..", "knowledge", "MANUAL.md")
if (!file.exists(manual)) {
  stop("Manual not found at '", manual, "'. Pass its path as the first argument.")
}

lines <- readLines(manual, warn = FALSE)

# Each section runs from its heading to the next top-level heading.
starts <- grep("^## ", lines)
heading <- function(pattern) {
  i <- grep(pattern, lines[starts])
  if (length(i) != 1) {
    stop("Expected exactly one '", pattern, "' heading in the manual; found ",
         length(i), ".")
  }
  starts[i]
}
section <- function(from, to) {
  stop_at <- if (is.na(to)) length(lines) else to - 1L
  lines[(from + 1L):stop_at]
}

i_problem <- heading("^## The problem and the response")
i_how <- heading("^## How the estimator works")
i_input <- heading("^## The input data")
i_outputs <- heading("^## The outputs and how to read them")
i_maths <- heading("^## Statistics, mathematics and implementation")
i_cond <- heading("^## Conditions of applicability")

# The manual nests one level deeper than a vignette, which carries its title in
# its own header, so its subsections move up a level.
promote <- function(x) sub("^### ", "## ", x)
trim <- function(x) {
  while (length(x) && !nzchar(x[1])) x <- x[-1]
  while (length(x) && !nzchar(x[length(x)])) x <- x[-length(x)]
  x
}

spec <- list(
  list(slug = "why", order = 1L, wide = FALSE,
       title = "Why CHORALE exists",
       description = paste("What integration across disjoint cohorts requires,",
                           "what existing methods provide, and what is missing."),
       body = section(i_problem, i_how)),
  list(slug = "how-it-works", order = 2L, wide = FALSE,
       title = "How the estimator works",
       description = "The six stages of the pipeline in plain terms.",
       body = section(i_how, i_input)),
  list(slug = "input-format", order = 3L, wide = TRUE,
       title = "The input data",
       description = paste("The container the tables become, the scale each",
                           "matrix is on, and what every column does."),
       body = section(i_input, i_outputs)),
  list(slug = "outputs", order = 5L, wide = TRUE,
       title = "The outputs and how to read them",
       description = paste("Every table the pipeline writes, what each column",
                           "means, and how to interpret it."),
       body = section(i_outputs, i_maths)),
  # The conditions and the sources belong with the mathematics they qualify.
  list(slug = "methods", order = 6L, wide = TRUE,
       title = "Statistics, mathematics and implementation",
       description = paste("Every quantity the pipeline computes, the equation",
                           "behind it, where it is implemented, and the",
                           "literature it rests on."),
       body = c(trim(section(i_maths, i_cond)), "", lines[i_cond:length(lines)]))
)

for (s in spec) {
  body <- trim(promote(s$body))
  yaml <- c(
    "---",
    sprintf('title: "%s"', s$title),
    sprintf('description: "%s"', s$description),
    sprintf("site_order: %d", s$order),
    sprintf("site_wide: %s", tolower(s$wide)),
    "output:",
    "  rmarkdown::html_vignette:",
    "    math_method: mathml",
    "vignette: >",
    sprintf("  %%\\VignetteIndexEntry{%s}", s$title),
    "  %\\VignetteEngine{knitr::rmarkdown}",
    "  %\\VignetteEncoding{UTF-8}",
    "---",
    ""
  )
  path <- file.path("vignettes", paste0(s$slug, ".Rmd"))
  writeLines(c(yaml, body), path)
  message(sprintf("%-16s %4d lines", basename(path), length(body)))
}

if (!file.exists(file.path("vignettes", "tutorial.Rmd"))) {
  warning("vignettes/tutorial.Rmd is absent; it is written by hand, not generated.")
} else {
  message("tutorial.Rmd     left untouched, written by hand")
}
