#' Coordinated variation that no concept explains
#'
#' The vocabulary bounds what the concept channel can say. Biology outside it
#' does not disappear: it is carried by the free dimensions, and a free
#' dimension that moves with the phenotype is a candidate the vocabulary has no
#' name for. Reporting those separately is what keeps the bound on the
#' vocabulary visible instead of leaving it to be inferred from an absence.
#'
#' Whether such a dimension exceeds noise is a question the same permutations
#' answer. The phenotype is permuted within strata of the design, the effect on
#' each free dimension is recomputed from the scores already fitted, and the
#' largest statistic over every free dimension of the collection forms the null
#' the observed values are read against. A dimension is therefore called only
#' when it beats the strongest thing the search could have turned up by chance.
#'
#' A free dimension has no name, so what it is made of is reported as the
#' features that load on it most heavily. That is a description, not an
#' identification.
#'
#' The concept channel is fitted first and the free dimensions are what remains
#' orthogonal to it, so a phenotype-linked direction that shares its direction
#' with a concept is carried by that concept and not seen here. What this table
#' reports is variation the vocabulary does not reach at all, which is the
#' question it is for.
#'
#' @param fit A `chorale_concept_fit`.
#' @param n_permutations Permutations calibrating the free dimensions. `NULL`
#'   uses the count the fit was calibrated with.
#' @param n_top Features named per dimension.
#' @param seed Integer seed.
#'
#' @returns A data frame with one row per free dimension per modality, carrying
#'   the share of variance it holds, how reproducibly it was recovered, its
#'   adjusted phenotype effect, the permutation p-values and the features that
#'   load on it most heavily.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
#'                            n_permutations = 99, n_init = 2)
#' chorale_free_dimensions(fit, n_permutations = 99)
chorale_free_dimensions <- function(fit, n_permutations = NULL, n_top = 10L,
                                    seed = 1L) {
  if (!inherits(fit, "chorale_concept_fit")) {
    rlang::abort("`fit` must be a chorale_concept_fit object.")
  }
  spec <- fit$evidence$spec
  terms <- chorale_profile_terms(spec$levels)
  phenotype_terms <- terms[startsWith(terms, paste0(spec$phenotype, "="))]
  n_perm <- n_permutations %||% fit$evidence$n_permutations
  designs <- fit$designs

  observed <- chorale_free_effects(fit, spec, phenotype_terms, designs)
  if (nrow(observed) == 0) return(chorale_empty_free_table())

  strata <- lapply(designs, function(d) {
    chorale_permutation_strata(d, spec$phenotype)
  })
  null_max <- rep(NA_real_, n_perm)
  null_by_dimension <- matrix(NA_real_, nrow = n_perm, ncol = nrow(observed))
  for (b in seq_len(n_perm)) {
    permuted <- Map(function(d, s, i) {
      chorale_permute_within(d, spec$phenotype, s, seed + 1000L * b + i)
    }, designs, strata, seq_along(designs))
    e <- chorale_free_effects(fit, spec, phenotype_terms, permuted)
    v <- abs(e$z[match(observed$key, e$key)])
    v[!is.finite(v)] <- 0
    null_by_dimension[b, ] <- v
    null_max[b] <- max(v, na.rm = TRUE)
  }

  statistic <- abs(observed$z)
  observed$p_value <- vapply(seq_along(statistic), function(i) {
    (1 + sum(null_by_dimension[, i] >= statistic[i], na.rm = TRUE)) /
      (1 + n_perm)
  }, numeric(1))
  observed$p_family <- vapply(statistic, function(s) {
    (1 + sum(null_max >= s, na.rm = TRUE)) / (1 + n_perm)
  }, numeric(1))
  observed$outside_vocabulary <- observed$p_family <= fit$evidence$alpha

  observed$top_features <- vapply(seq_len(nrow(observed)), function(i) {
    l <- fit$encoding$encodings[[observed$modality[i]]]$free_loadings
    if (!observed$dimension[i] %in% colnames(l)) return(NA_character_)
    v <- l[, observed$dimension[i]]
    paste(names(sort(abs(v), decreasing = TRUE))[seq_len(min(n_top, length(v)))],
          collapse = ", ")
  }, character(1))

  observed$key <- NULL
  observed <- observed[order(-statistic), , drop = FALSE]
  rownames(observed) <- NULL
  observed
}

#' @keywords internal
#' @noRd
chorale_empty_free_table <- function() {
  data.frame(modality = character(), dimension = character(),
             variance_share = numeric(), reproducibility = numeric(),
             term = character(), effect = numeric(), se = numeric(),
             z = numeric(), p_value = numeric(), p_family = numeric(),
             outside_vocabulary = logical(), top_features = character(),
             stringsAsFactors = FALSE)
}

#' Adjusted phenotype effect on every free dimension
#' @keywords internal
#' @noRd
chorale_free_effects <- function(fit, spec, phenotype_terms, designs) {
  rows <- list()
  for (m in fit$modalities) {
    e <- fit$encoding$encodings[[m]]
    scores <- e$free_scores
    if (ncol(scores) == 0) next
    profile <- chorale_adjusted_profile(scores, designs[[m]], spec)
    keep <- intersect(phenotype_terms, colnames(profile$effects))
    if (length(keep) == 0) next
    stability <- if (is.data.frame(e$stability) &&
                     "subspace" %in% colnames(e$stability)) {
      chorale_finite_mean(e$stability$subspace)
    } else {
      NA_real_
    }
    share <- fit$encoding$variance$free_share[
      fit$encoding$variance$modality == m]
    for (term in keep) {
      rows[[length(rows) + 1L]] <- data.frame(
        key = paste(m, colnames(scores), term, sep = "|"),
        modality = m,
        dimension = colnames(scores),
        # The share is the whole free channel's; it is divided evenly rather
        # than attributed, since the components are not ordered by variance.
        variance_share = round(share / ncol(scores), 5),
        reproducibility = round(stability, 3),
        term = term,
        effect = round(unname(profile$effects[, term]), 4),
        se = round(unname(profile$se[, term]), 4),
        z = round(unname(profile$z[, term]), 4),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(chorale_empty_free_table())
  do.call(rbind, rows)
}

#' @rdname chorale_added_value
#' @export
chorale_added_value.chorale_concept_fit <- function(fit, ...) {
  j <- fit$evidence$joint
  if (is.null(j) || nrow(j) == 0) return(data.frame())
  out <- data.frame(
    concept = j$concept,
    term = j$term,
    n_modalities = j$n_modalities,
    modalities = j$modalities,
    joint_z = j$joint_z,
    best_single_z = j$best_single_z,
    # On the scale the statistics are already on, so a positive margin means the
    # concept is stronger across the modalities than in any one of them.
    margin = round(abs(j$joint_z) - abs(j$best_single_z), 3),
    sign_agreement = j$sign_agreement,
    # One modality cannot need more than one modality, whatever its statistic.
    needs_multiple = j$n_modalities > 1L &
      abs(j$joint_z) > abs(j$best_single_z),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$margin), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @rdname chorale_report
#' @param n_top Features named per free dimension.
#' @export
chorale_report.chorale_concept_fit <- function(fit, null = NULL, path,
                                               n_top = 10L, ...) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  written <- character()

  concepts <- chorale_concept_table(fit)
  written <- c(written, chorale_write(concepts, path, "concepts.tsv"))
  written <- c(written, chorale_write(fit$evidence$per_modality, path,
                                      "concept_effects.tsv"))
  written <- c(written, chorale_write(fit$concepts$coverage, path,
                                      "concept_coverage.tsv"))
  written <- c(written, chorale_write(fit$encoding$variance, path,
                                      "variance.tsv"))
  written <- c(written, chorale_write(fit$encoding$concept_variance, path,
                                      "concept_variance.tsv"))

  free <- chorale_free_dimensions(fit, n_top = n_top)
  written <- c(written, chorale_write(free, path, "free_dimensions.tsv"))

  added <- chorale_added_value(fit)
  written <- c(written, chorale_write(added, path, "added_value.tsv"))

  # A run without controls writes no control file, rather than an empty one a
  # reader would have to open to discover is empty.
  controls <- if (!is.null(null)) null$controls else data.frame()
  if (nrow(controls) > 0) {
    written <- c(written, chorale_write(controls, path, "controls.tsv"))
  }

  for (m in fit$modalities) {
    e <- fit$encoding$encodings[[m]]
    written <- c(written, chorale_write(
      data.frame(sample_id = rownames(e$concept_scores),
                 as.data.frame(e$concept_scores), check.names = FALSE),
      path, paste0("concept_scores_", m, ".tsv")))
    if (ncol(e$free_scores) > 0) {
      written <- c(written, chorale_write(
        data.frame(sample_id = rownames(e$free_scores),
                   as.data.frame(e$free_scores), check.names = FALSE),
        path, paste0("free_scores_", m, ".tsv")))
    }
  }

  written <- c(written, chorale_write_concept_html(fit, concepts, free, added,
                                                   controls, path))
  invisible(written)
}

#' One row per concept, with everything the report claims about it
#' @keywords internal
#' @noRd
chorale_concept_table <- function(fit) {
  j <- fit$evidence$joint
  if (is.null(j) || nrow(j) == 0) return(data.frame())
  pm <- fit$evidence$per_modality
  coverage <- fit$concepts$summary

  wide <- lapply(fit$modalities, function(m) {
    d <- pm[pm$modality == m, , drop = FALSE]
    key <- paste(d$concept, d$term, sep = "|")
    stats::setNames(
      list(round(d$effect[match(paste(j$concept, j$term, sep = "|"), key)], 4),
           round(d$z[match(paste(j$concept, j$term, sep = "|"), key)], 4)),
      paste0(c("effect_", "z_"), m))
  })
  wide <- do.call(c, wide)

  out <- data.frame(
    concept = j$concept,
    term = j$term,
    n_modalities = j$n_modalities,
    modalities = j$modalities,
    in_all_modalities = coverage$in_all_modalities[
      match(j$concept, coverage$concept)],
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(wide, stringsAsFactors = FALSE))
  added <- chorale_added_value(fit)
  out$joint_effect <- j$joint_effect
  out$joint_z <- j$joint_z
  out$sign_agreement <- j$sign_agreement
  out$heterogeneity_p <- j$heterogeneity_p
  out$max_jaccard <- j$max_jaccard
  out$attributed_z <- j$attributed_z
  out$p_value <- j$p_value
  out$p_family <- j$p_family
  out$q_value <- j$q_value
  out$significant <- j$significant
  out$beats_best_single <- added$needs_multiple[
    match(paste(j$concept, j$term), paste(added$concept, added$term))]
  out
}

#' A page a reader can take the result from
#' @keywords internal
#' @noRd
chorale_write_concept_html <- function(fit, concepts, free, added, controls,
                                       path) {
  file <- file.path(path, "report.html")
  n_supported <- if (nrow(concepts) > 0) sum(concepts$significant) else 0L
  head <- utils::head(concepts, 25)
  parts <- c(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<title>chorale concept report</title>",
    "<style>body{font-family:system-ui,sans-serif;margin:2rem;max-width:70rem}",
    "table{border-collapse:collapse;margin:1rem 0}",
    "th,td{border:1px solid #999;padding:.25rem .5rem;font-size:.85rem}",
    "caption{text-align:left;font-weight:600;margin-bottom:.25rem}",
    "p{max-width:44rem}</style></head><body>",
    "<h1>What the collection shares</h1>",
    paste0("<p>The modalities were measured on different individuals and share ",
           "no sample. They meet in a vocabulary of ",
           length(fit$concepts$vocabulary), " named concepts, of which ",
           fit$concepts$n_in_all, " are expressible in every modality. Of the ",
           nrow(concepts), " tested, ", n_supported,
           " separate cases from controls at a false discovery rate of ",
           fit$evidence$alpha, ".</p>"),
    chorale_html_table(
      head, "Concepts",
      paste0("One row per concept. The effect columns give the adjusted ",
             "case-control contrast in each modality; the joint statistic ",
             "combines them by inverse-variance weighting, so effects in ",
             "opposite directions cancel. Sign agreement says how far the ",
             "modalities agree. The attributed statistic is what survives ",
             "once overlapping concepts are regressed out.")),
    chorale_html_table(
      utils::head(free, 15), "Free dimensions",
      paste0("Coordinated variation no concept explains. A dimension called ",
             "here is a candidate outside the vocabulary, described by the ",
             "features that load on it and not identified by them.")),
    chorale_html_table(
      utils::head(added, 15), "Added value",
      paste0("Whether a concept's evidence across the modalities exceeds the ",
             "best any one of them carries. A concept expressed in one ",
             "modality cannot need more than one.")),
    if (nrow(controls) > 0) {
      chorale_html_table(
        controls, "Controls",
        paste0("Each control with the smallest p-value it can attain, so a ",
               "control that cannot reach the threshold it is read against is ",
               "visible as one."))
    } else {
      ""
    },
    "</body></html>")
  writeLines(parts, file)
  file
}
