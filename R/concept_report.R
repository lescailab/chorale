#' Emit interpretable outputs from a concept fit
#'
#' Writes the concept table, the free-dimension table and the added-value table,
#' with the coverage and variance behind them, a page that leads with what the
#' collection shares, and the fit itself so an analysis can be resumed without
#' re-reading anything.
#'
#' @param fit A `chorale_concept_fit`.
#' @param null A `chorale_null` object, as returned by [chorale_null()]. Where it
#'   is absent no control table is written.
#' @param path Directory to write into. Created if absent.
#' @param family_evidence Optional fitted [chorale_family_evidence()] object.
#' @param joint_state Optional fitted [chorale_joint_state()] object.
#' @param joint_evidence Optional fitted [chorale_joint_evidence()] object.
#' @param joint_transfer Optional fitted [chorale_joint_transfer()] object.
#'   Reporting never computes these analyses: only objects supplied here are
#'   written, so component counts, nuisance covariates, family thresholds and
#'   permutation counts remain explicit analysis decisions.
#' @param ... Passed to the method.
#' @returns Invisibly, a character vector of the files written.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
#'                            n_permutations = 99, n_init = 2)
#' basename(chorale_report(fit, path = withr::local_tempdir()))
chorale_report <- function(fit, null = NULL, path, ...) {
  UseMethod("chorale_report")
}

#' @export
chorale_report.default <- function(fit, null = NULL, path, ...) {
  rlang::abort("`fit` must be a chorale_concept_fit object.")
}

#' Whether a concept needed more than one modality
#'
#' The claim that integration produced a result is only worth making if the
#' result is stronger than what any one modality could have produced alone. The
#' combined statistic is therefore reported beside the best any single modality
#' reached, and a concept only one modality expresses cannot need more than one
#' whatever its statistic.
#'
#' @param fit A `chorale_concept_fit`.
#' @param ... Passed to the method.
#' @returns A data frame with one row per concept carrying the combined
#'   statistic, the best single-modality statistic, the margin between them and
#'   whether more than one modality was needed.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
#'                            n_permutations = 99)
#' chorale_added_value(fit)
chorale_added_value <- function(fit, ...) {
  UseMethod("chorale_added_value")
}

#' @export
chorale_added_value.default <- function(fit, ...) {
  rlang::abort("`fit` must be a chorale_concept_fit object.")
}

#' Coordinated variation that no concept explains
#'
#' The vocabulary bounds what the concept channel can say. Biology outside it
#' does not disappear: it is carried by the free dimensions, and a free
#' dimension that moves with the phenotype is a candidate the vocabulary has no
#' name for. Reporting those separately is what keeps the bound on the
#' vocabulary visible instead of leaving it to be inferred from an absence.
#'
#' Whether such a dimension exceeds noise is a question the same null answers.
#' The design is held fixed and the part of each free score the adjusting
#' covariates do not explain is exchanged, and the largest statistic over every
#' free dimension of the collection forms the null the observed values are read
#' against. A dimension is therefore called only
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

  blocks <- lapply(designs, chorale_exchangeability_blocks,
                   columns = fit$control$exchangeability_blocks)
  null_max <- rep(NA_real_, n_perm)
  null_by_dimension <- matrix(NA_real_, nrow = n_perm, ncol = nrow(observed))
  for (b in seq_len(n_perm)) {
    # The design is held fixed and the response is rebuilt, so the phenotype
    # keeps the relation to the covariates it has in the data.
    permuted_scores <- Map(function(m, d, block, i) {
      scores <- fit$encoding$encodings[[m]]$free_scores
      # A modality whose vocabulary spans its sample space has no free
      # dimensions, so there is nothing to rebuild and nothing to exchange.
      if (is.null(scores) || ncol(scores) == 0) return(scores)
      chorale_freedman_lane_scores(scores, d, spec, blocks = block,
                                   seed = seed + 1000L * b + i)
    }, fit$modalities, designs, blocks, seq_along(designs))
    e <- chorale_free_effects(fit, spec, phenotype_terms, designs,
                              scores_by_modality = permuted_scores)
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

#' What each free dimension reconstructs of its modality
#'
#' Independent component analysis returns its components in no particular order,
#' which is sometimes taken to mean they contribute equally. They do not: the
#' fitted loadings can differ in norm by a wide margin, so dividing the channel's
#' share evenly among them would report invented numbers.
#'
#' Each component is a rank-one term, so what it reconstructs is the product of
#' the squared norms of its score and its loading. The components are
#' uncorrelated after whitening, so those contributions should sum to the share
#' the channel carries as a whole. Where they do not, the cross-terms are
#' material and no per-dimension attribution is reported rather than an
#' allocation being forced.
#'
#' @keywords internal
#' @noRd
chorale_component_variance <- function(encoding, channel_share,
                                       tolerance = 0.05) {
  scores <- encoding$free_scores
  loadings <- encoding$free_loadings
  total <- sum(encoding$analysis_matrix^2)
  if (ncol(scores) == 0 || !is.finite(total) || total <= 0) {
    return(rep(NA_real_, ncol(scores)))
  }
  # The reconstruction of one component is the outer product of its score and
  # its loading, whose squared Frobenius norm is the product of theirs.
  share <- vapply(seq_len(ncol(scores)), function(j) {
    sum(scores[, j]^2) * sum(loadings[, j]^2) / total
  }, numeric(1))
  if (is.finite(channel_share) &&
      abs(sum(share) - channel_share) > tolerance * max(channel_share, 1e-8)) {
    return(rep(NA_real_, ncol(scores)))
  }
  round(share, 5)
}

#' Adjusted phenotype effect on every free dimension
#' @keywords internal
#' @noRd
chorale_free_effects <- function(fit, spec, phenotype_terms, designs,
                                 scores_by_modality = NULL) {
  rows <- list()
  for (m in fit$modalities) {
    e <- fit$encoding$encodings[[m]]
    scores <- if (is.null(scores_by_modality)) {
      e$free_scores
    } else {
      scores_by_modality[[m]]
    }
    if (is.null(scores) || ncol(scores) == 0) next
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
    component_share <- chorale_component_variance(e, share)
    for (term in keep) {
      rows[[length(rows) + 1L]] <- data.frame(
        key = paste(m, colnames(scores), term, sep = "|"),
        modality = m,
        dimension = colnames(scores),
        variance_share = component_share,
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
chorale_report.chorale_concept_fit <- function(
    fit, null = NULL, path, family_evidence = NULL, joint_state = NULL,
    joint_evidence = NULL, joint_transfer = NULL, n_top = 10L, ...) {
  chorale_validate_report_additions(
    fit, family_evidence, joint_state, joint_evidence, joint_transfer)
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

  if (!is.null(family_evidence)) {
    membership <- chorale_report_family_membership(family_evidence)
    written <- c(written, chorale_write(membership, path,
                                        "concept_families.tsv"))
    written <- c(written, chorale_write(family_evidence$families, path,
                                        "family_evidence.tsv"))
  }

  if (!is.null(joint_state)) {
    scores <- data.frame(
      sample_id = rownames(joint_state$scores),
      modality = joint_state$modality,
      as.data.frame(joint_state$scores, check.names = FALSE),
      check.names = FALSE, stringsAsFactors = FALSE)
    loadings <- data.frame(
      concept = rownames(joint_state$loadings),
      as.data.frame(joint_state$loadings, check.names = FALSE),
      check.names = FALSE, stringsAsFactors = FALSE)
    nuisance <- data.frame(
      modality = names(joint_state$nuisance),
      covariates = vapply(joint_state$nuisance, function(x) {
        if (length(x)) paste(x, collapse = "; ") else ""
      }, character(1)),
      stringsAsFactors = FALSE)
    written <- c(written, chorale_write(scores, path, "joint_scores.tsv"))
    written <- c(written, chorale_write(loadings, path,
                                        "joint_loadings.tsv"))
    written <- c(written, chorale_write(joint_state$variance, path,
                                        "joint_variance.tsv"))
    written <- c(written, chorale_write(joint_state$coverage, path,
                                        "joint_coverage.tsv"))
    written <- c(written, chorale_write(nuisance, path,
                                        "joint_nuisance.tsv"))
  }

  if (!is.null(joint_evidence)) {
    written <- c(written, chorale_write(joint_evidence$components, path,
                                        "joint_components.tsv"))
    written <- c(written, chorale_write(joint_evidence$per_modality, path,
                                        "joint_component_effects.tsv"))
  }

  if (!is.null(joint_transfer)) {
    written <- c(written, chorale_write(joint_transfer$transfer, path,
                                        "joint_transfer.tsv"))
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

  # The fit itself, so an analysis can be resumed without re-reading anything
  # and without the tables having to be lossless.
  fit_file <- file.path(path, "fit.rds")
  saveRDS(fit, fit_file)
  written <- c(written, fit_file)

  written <- c(written, chorale_write_concept_html(
    fit, concepts, free, added, controls, path,
    family_evidence = family_evidence, joint_state = joint_state,
    joint_evidence = joint_evidence, joint_transfer = joint_transfer))
  invisible(written)
}

#' Validate optional analyses supplied to the report
#' @keywords internal
#' @noRd
chorale_validate_report_additions <- function(
    fit, family_evidence, joint_state, joint_evidence, joint_transfer) {
  expected <- list(
    family_evidence = "chorale_family_evidence",
    joint_state = "chorale_joint_state",
    joint_evidence = "chorale_joint_evidence",
    joint_transfer = "chorale_joint_transfer")
  supplied <- list(family_evidence = family_evidence, joint_state = joint_state,
                   joint_evidence = joint_evidence,
                   joint_transfer = joint_transfer)
  for (nm in names(expected)) {
    if (!is.null(supplied[[nm]]) && !inherits(supplied[[nm]], expected[[nm]])) {
      rlang::abort(paste0("`", nm, "` must be a ", expected[[nm]],
                          " object."),
                   class = "chorale_invalid_report_analysis")
    }
  }

  vocabulary <- fit$concepts$vocabulary
  if (!is.null(family_evidence)) {
    membership <- chorale_report_family_membership(family_evidence)
    unknown <- setdiff(membership$concept, vocabulary)
    if (length(unknown)) {
      rlang::abort("`family_evidence` was fitted on a different vocabulary.",
                   class = "chorale_incompatible_report_analysis")
    }
  }
  if (!is.null(joint_state)) {
    unknown <- setdiff(rownames(joint_state$loadings), vocabulary)
    if (length(unknown)) {
      rlang::abort("`joint_state` was fitted on a different vocabulary.",
                   class = "chorale_incompatible_report_analysis")
    }
  }
  if (!is.null(joint_state) && !is.null(joint_evidence)) {
    unknown <- setdiff(unique(joint_evidence$components$component),
                       colnames(joint_state$scores))
    if (length(unknown)) {
      rlang::abort("`joint_evidence` does not describe `joint_state`.",
                   class = "chorale_incompatible_report_analysis")
    }
  }
  invisible(TRUE)
}

#' Recover one row per concept-family assignment for reporting
#' @keywords internal
#' @noRd
chorale_report_family_membership <- function(family_evidence) {
  membership <- family_evidence$membership
  if (!is.null(membership)) {
    return(membership[, c("concept", "family", "family_size"), drop = FALSE])
  }
  rows <- lapply(seq_len(nrow(family_evidence$families)), function(i) {
    r <- family_evidence$families[i, , drop = FALSE]
    concepts <- trimws(strsplit(r$members, ";", fixed = TRUE)[[1]])
    data.frame(concept = concepts, family = r$family,
               family_size = length(concepts), stringsAsFactors = FALSE)
  })
  if (!length(rows)) {
    return(data.frame(concept = character(), family = character(),
                      family_size = integer(), stringsAsFactors = FALSE))
  }
  unique(do.call(rbind, rows))
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
chorale_write_concept_html <- function(
    fit, concepts, free, added, controls, path, family_evidence = NULL,
    joint_state = NULL, joint_evidence = NULL, joint_transfer = NULL) {
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
           " pass the configured phenotype test at a false discovery rate of ",
           fit$evidence$alpha, ".</p>"),
    chorale_html_table(
      head, "Concepts",
      paste0("One row per concept. The effect columns give the adjusted ",
             "phenotype contrast in each modality; the joint statistic ",
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
    "<h2>Additional analyses</h2>",
    chorale_report_analysis_status(family_evidence, joint_state,
                                   joint_evidence, joint_transfer),
    if (!is.null(family_evidence)) {
      chorale_html_table(
        utils::head(family_evidence$families, 25), "Concept families",
        paste0("Related concepts tested as a directional group. The family ",
               "p-value is calibrated for that family; p_family compares it ",
               "with the strongest family in each permutation."))
    } else "",
    if (!is.null(joint_state)) {
      chorale_html_table(
        joint_state$variance, "Joint state",
        paste0("A low-rank representation of the stacked, within-modality ",
               "standardised concept scores. Component sign and order are ",
               "arbitrary; a component is not an identified mechanism."))
    } else "",
    if (!is.null(joint_evidence)) {
      chorale_html_table(
        utils::head(joint_evidence$components, 25), "Joint component evidence",
        paste0("Adjusted phenotype effects on the joint component scores. ",
               "Read magnitude and calibrated error rates; the sign depends ",
               "on the arbitrary orientation of the component."))
    } else "",
    if (!is.null(joint_transfer)) {
      chorale_html_table(
        utils::head(joint_transfer$transfer, 25), "Joint transfer",
        paste0("Each direction was fitted without the held-out modality. ",
               "p_value is component-specific, not family-wise adjusted; ",
               "loading_agreement records orientation against the reference."))
    } else "",
    "</body></html>")
  writeLines(parts, file)
  file
}

#' State which optional analyses were supplied to a report
#' @keywords internal
#' @noRd
chorale_report_analysis_status <- function(family_evidence, joint_state,
                                           joint_evidence, joint_transfer) {
  supplied <- list(
    "concept-family evidence" = family_evidence,
    "joint state" = joint_state,
    "joint component evidence" = joint_evidence,
    "leave-one-modality-out transfer" = joint_transfer)
  status <- data.frame(
    analysis = names(supplied),
    status = vapply(supplied, function(x) {
      if (is.null(x)) "not supplied" else "supplied"
    }, character(1)),
    stringsAsFactors = FALSE)
  chorale_html_table(
    status, "Analysis status",
    paste0("Not supplied means the analysis was not included in this report; ",
           "it is not a negative or nonsignificant result."))
}


#' @keywords internal
#' @noRd
chorale_write <- function(x, path, file) {
  f <- file.path(path, file)
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

#' Escape text for HTML
#' @keywords internal
#' @noRd
chorale_esc <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' A table with its legend
#' @keywords internal
#' @noRd
chorale_html_table <- function(d, caption, legend) {
  head_block <- paste0("<h3>", chorale_esc(caption), "</h3>",
                       "<p class='legend'>", legend, "</p>")
  if (is.null(d) || nrow(d) == 0) {
    return(paste0(head_block, "<p class='empty'>No rows.</p>"))
  }
  hdr <- paste0("<tr>", paste0("<th>", chorale_esc(colnames(d)), "</th>",
                               collapse = ""), "</tr>")
  body <- apply(d, 1, function(r) {
    paste0("<tr>", paste0("<td>", chorale_esc(r), "</td>", collapse = ""), "</tr>")
  })
  paste0(head_block, "<div class='scroll'><table>", hdr,
         paste(body, collapse = ""), "</table></div>")
}
