#' Whether a jointly estimated component separates cases from controls
#'
#' A component of [chorale_joint_state()] is one direction in concept space,
#' carried by every sample in the collection. This asks whether a sample's
#' position along it moves with the phenotype, adjusting for the covariates the
#' modalities share and for the modality itself.
#'
#' The modality term matters. Without it a component that merely orders the
#' modalities would appear to separate cases from controls whenever the
#' modalities differ in case-control composition. With it, only variation
#' between samples measured the same way contributes, which is the variation the
#' stacking was arranged to keep.
#'
#' Calibration follows [chorale_concept_evidence()]: the design is held fixed,
#' the reduced model without the phenotype is fitted, and its residuals are
#' exchanged. Exchange happens **within modality**, so a permuted collection has
#' the same modality composition as the observed one and the null describes the
#' phenotype coefficient rather than the modality difference.
#'
#' @section What is reported and what is not required:
#' Each component's effect is also read within each modality on its own, and
#' those per-modality effects are reported beside the joint one. They are a
#' description, not a condition: a component carried by every modality and a
#' component carried by one are different findings, and the reader is given the
#' numbers to tell them apart rather than a rule that discards one of them.
#'
#' @param state A `chorale_joint_state` object.
#' @param n_permutations Permutations calibrating the components.
#' @param control A [chorale_control()] object.
#' @param seed Integer seed.
#' @param anchor The design covariate whose effect is tested. `NULL` tests the
#'   phenotype.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_joint_evidence` with `components`, one
#'   row per component per phenotype term carrying the effect, its standard
#'   error, the permutation p-values and the error control; and
#'   `per_modality`, the same effect estimated within each modality alone.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' ids <- sprintf("feature_%05d", seq_len(80))
#' sim$modalities <- lapply(sim$modalities, function(m) {
#'   rownames(m) <- ids
#'   m
#' })
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' cc <- chorale_concepts(containers, list(one = ids[1:20], two = ids[21:45]),
#'                        min_features = 5)
#' enc <- chorale_encode(containers, cc, n_free = 1, n_init = 2)
#' state <- chorale_joint_state(enc, n_components = 1)
#' chorale_joint_evidence(state, n_permutations = 99)
chorale_joint_evidence <- function(state, n_permutations = 999L,
                                   control = chorale_control(), seed = 1L,
                                   anchor = NULL, ...) {
  if (!inherits(state, "chorale_joint_state")) {
    rlang::abort("`state` must be a chorale_joint_state object.")
  }
  control <- chorale_merge_control(control, list(...))
  if (state$n_components == 0) {
    return(structure(list(components = chorale_empty_joint_components(),
                          per_modality = chorale_empty_joint_per_modality(),
                          anchor = NA_character_, n_permutations = 0L,
                          alpha = control$alpha),
                     class = "chorale_joint_evidence"))
  }

  spec <- chorale_joint_signature(state, control)
  anchor <- anchor %||% spec$phenotype
  anchor_terms <- chorale_anchor_terms(spec, anchor)
  if (length(anchor_terms) == 0) {
    rlang::abort(paste0("`", anchor,
                        "` contributes no term to compare components on."),
                 class = "chorale_unusable_anchor")
  }

  observed <- chorale_joint_effects(state, spec, anchor_terms)
  statistic <- abs(observed$z)

  # Exchange happens inside a modality, so every permuted collection has the
  # composition the study produced.
  blocks <- paste(state$modality,
                  chorale_exchangeability_blocks(
                    state$design, control$exchangeability_blocks),
                  sep = "|")
  null_by_component <- matrix(NA_real_, nrow = n_permutations,
                              ncol = length(statistic))
  null_max <- rep(NA_real_, n_permutations)
  for (b in seq_len(n_permutations)) {
    permuted <- chorale_freedman_lane_scores(
      state$scores, state$design, spec, anchor = anchor, blocks = blocks,
      seed = seed + 1000L * b)
    e <- chorale_joint_effects(state, spec, anchor_terms, scores = permuted)
    v <- abs(e$z)
    v[!is.finite(v)] <- 0
    null_by_component[b, ] <- v
    null_max[b] <- max(v, na.rm = TRUE)
  }

  components <- observed$table
  components$p_value <- vapply(seq_along(statistic), function(i) {
    (1 + sum(null_by_component[, i] >= statistic[i], na.rm = TRUE)) /
      (1 + n_permutations)
  }, numeric(1))
  components$p_family <- vapply(statistic, function(s) {
    (1 + sum(null_max >= s, na.rm = TRUE)) / (1 + n_permutations)
  }, numeric(1))
  components$q_value <- chorale_permutation_fdr(statistic, null_by_component)
  components$significant <- components$q_value <= control$alpha
  components$family_significant <- components$p_family <= control$alpha

  structure(
    list(components = components,
         per_modality = chorale_joint_per_modality(state, spec, anchor_terms),
         spec = spec, anchor = anchor,
         null = null_max,
         n_permutations = as.integer(n_permutations),
         alpha = control$alpha,
         smallest_attainable_p = 1 / (n_permutations + 1)),
    class = "chorale_joint_evidence")
}

#' The design signature of a stacked collection
#'
#' The covariates are those the modalities share, resolved exactly as
#' [chorale_concept_evidence()] resolves them, with the modality itself added as
#' a further term so that a component ordering the modalities is not read as a
#' phenotype effect.
#'
#' @keywords internal
#' @noRd
chorale_joint_signature <- function(state, control) {
  designs <- split(state$design, state$modality)
  spec <- chorale_resolve_signature(
    designs,
    phenotype_column = control$phenotype_column,
    phenotype_reference = control$phenotype_reference,
    profile_covariates = control$profile_covariates)
  modalities <- sort(unique(state$modality))
  if (length(modalities) > 1) {
    spec$covariates <- c(spec$covariates, "modality")
    spec$levels[["modality"]] <- modalities
    spec$secondary <- c(spec$secondary, "modality")
  }
  spec
}

#' @keywords internal
#' @noRd
chorale_joint_effects <- function(state, spec, anchor_terms, scores = NULL) {
  scores <- scores %||% state$scores
  design <- state$design
  design$sample_id <- rownames(state$scores)
  profile <- chorale_adjusted_profile(scores, design, spec)
  keep <- intersect(anchor_terms, colnames(profile$effects))
  rows <- list()
  for (term in keep) {
    rows[[term]] <- data.frame(
      component = rownames(profile$effects),
      term = term,
      share = state$variance$share[match(rownames(profile$effects),
                                         state$variance$component)],
      effect = unname(profile$effects[, term]),
      se = unname(profile$se[, term]),
      z = unname(profile$z[, term]),
      stringsAsFactors = FALSE)
  }
  table <- if (length(rows)) do.call(rbind, rows) else
    chorale_empty_joint_components()
  rownames(table) <- NULL
  list(table = table, z = table$z)
}

#' The same effect read within each modality alone
#'
#' Reported so that a component the whole collection carries can be told from
#' one a single modality carries. Nothing is filtered on it.
#'
#' @keywords internal
#' @noRd
chorale_joint_per_modality <- function(state, spec, anchor_terms) {
  local_spec <- spec
  local_spec$covariates <- setdiff(spec$covariates, "modality")
  local_spec$secondary <- setdiff(spec$secondary, "modality")
  local_spec$levels[["modality"]] <- NULL

  rows <- list()
  for (m in sort(unique(state$modality))) {
    at <- state$modality == m
    scores <- state$scores[at, , drop = FALSE]
    design <- state$design[at, , drop = FALSE]
    design$sample_id <- rownames(scores)
    profile <- try(chorale_adjusted_profile(scores, design, local_spec),
                   silent = TRUE)
    if (inherits(profile, "try-error")) next
    for (term in intersect(anchor_terms, colnames(profile$effects))) {
      rows[[length(rows) + 1L]] <- data.frame(
        component = rownames(profile$effects), modality = m, term = term,
        n_samples = sum(at),
        effect = unname(profile$effects[, term]),
        se = unname(profile$se[, term]),
        z = unname(profile$z[, term]),
        stringsAsFactors = FALSE)
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    chorale_empty_joint_per_modality()
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
chorale_empty_joint_components <- function() {
  data.frame(component = character(), term = character(), share = numeric(),
             effect = numeric(), se = numeric(), z = numeric(),
             stringsAsFactors = FALSE)
}

#' @keywords internal
#' @noRd
chorale_empty_joint_per_modality <- function() {
  data.frame(component = character(), modality = character(),
             term = character(), n_samples = integer(), effect = numeric(),
             se = numeric(), z = numeric(), stringsAsFactors = FALSE)
}

#' @export
print.chorale_joint_evidence <- function(x, ...) {
  cat("<chorale_joint_evidence>\n")
  if (nrow(x$components) == 0) {
    cat("  no component to test\n")
    return(invisible(x))
  }
  cat("  anchor: ", x$anchor, "\n", sep = "")
  print(utils::head(x$components[order(-abs(x$components$z)), ], 10))
  invisible(x)
}

#' Whether a jointly estimated direction predicts phenotype in a modality that
#' did not contribute to it
#'
#' The control a permutation cannot supply. Each modality is held out in turn,
#' the joint state is estimated from the remaining ones, and the held-out
#' modality's samples are projected onto the concept loadings that came out.
#' Whether those projections separate cases from controls is then tested in the
#' held-out modality alone.
#'
#' The loadings carry nothing from the held-out modality: not its samples, not
#' its batches, not its assay, not its preparation. Structure that a modality's
#' own nuisance produced can inflate a component fitted with that modality in,
#' but it cannot produce a direction over concepts that predicts phenotype in a
#' cohort sharing none of it. Transfer is therefore evidence about the biology
#' and not only about the arithmetic.
#'
#' It is also the positive claim rather than only its defence: a direction
#' estimated without a cohort, which orders that cohort by disease state, is
#' what it means for the collection to carry something no member carries.
#'
#' @param encoding A `chorale_encode` object.
#' @param n_components Components to fit on each training subset.
#' @param n_permutations Permutations calibrating the held-out test.
#' @param nuisance Layer-local nuisance covariates, as in
#'   [chorale_joint_state()].
#' @param control A [chorale_control()] object.
#' @param seed Integer seed.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_joint_transfer` with `transfer`, one row
#'   per held-out modality per component per phenotype term, carrying the effect
#'   in the held-out modality and its permutation p-value.
#' @export
chorale_joint_transfer <- function(encoding, n_components = 2L,
                                   n_permutations = 999L, nuisance = NULL,
                                   control = chorale_control(), seed = 1L,
                                   ...) {
  if (!inherits(encoding, "chorale_encode")) {
    rlang::abort("`encoding` must be a chorale_encode object.")
  }
  control <- chorale_merge_control(control, list(...))
  modalities <- encoding$modalities
  if (length(modalities) < 3) {
    rlang::abort(paste0("Holding a modality out needs at least three: with two",
                        " the training set is a single modality and the fit is",
                        " no longer joint."),
                 class = "chorale_transfer_needs_three")
  }

  # Each held-out fit is its own decomposition, so its second component is not
  # the second component of any other fit. A reference fitted on the whole
  # collection supplies one set of names, and each held-out component is matched
  # to it by how far their concept loadings agree. The reference contributes
  # nothing to the statistic: the direction tested is the one fitted without the
  # held-out modality, and the p-value is on the magnitude of its effect, so
  # neither the name nor the orientation can change what is concluded.
  reference <- chorale_joint_state(encoding, n_components = n_components,
                                   nuisance = nuisance, control = control,
                                   seed = seed)

  rows <- list()
  for (held in modalities) {
    trained_on <- setdiff(modalities, held)
    trained <- chorale_subset_encoding(encoding, trained_on)
    # Subsetting a list by names it does not carry yields entries named NA, so
    # the modalities to keep are intersected rather than indexed directly: a
    # collection where only some modalities declare nuisance covariates is the
    # ordinary case, not an error.
    state <- chorale_joint_state(
      trained, n_components = n_components,
      nuisance = nuisance[intersect(names(nuisance), trained_on)],
      control = control, seed = seed)
    if (state$n_components == 0) next

    projected <- chorale_project_modality(encoding, held, state$loadings,
                                          nuisance[[held]],
                                          control$phenotype_column)
    if (is.null(projected)) next

    design <- encoding$designs[[held]]
    design <- design[match(rownames(projected), design$sample_id), ,
                     drop = FALSE]
    spec <- chorale_resolve_signature(
      list(design, design),
      phenotype_column = control$phenotype_column,
      phenotype_reference = control$phenotype_reference,
      profile_covariates = control$profile_covariates)
    anchor_terms <- chorale_anchor_terms(spec, spec$phenotype)
    profile <- chorale_adjusted_profile(projected, design, spec)

    blocks <- chorale_exchangeability_blocks(design,
                                             control$exchangeability_blocks)
    null_z <- matrix(NA_real_, nrow = n_permutations, ncol = ncol(projected))
    for (b in seq_len(n_permutations)) {
      permuted <- chorale_freedman_lane_scores(
        projected, design, spec, anchor = spec$phenotype, blocks = blocks,
        seed = seed + 1000L * b)
      p <- chorale_adjusted_profile(permuted, design, spec)
      term <- intersect(anchor_terms, colnames(p$z))[1]
      null_z[b, ] <- abs(p$z[, term])
    }

    matched <- chorale_match_components(state$loadings, reference$loadings)

    for (term in intersect(anchor_terms, colnames(profile$effects))) {
      statistic <- abs(profile$z[, term])
      rows[[length(rows) + 1L]] <- data.frame(
        held_out = held,
        trained_on = paste(trained_on, collapse = ", "),
        component = matched$reference[match(rownames(profile$effects),
                                            matched$fitted)],
        fitted_component = rownames(profile$effects),
        loading_agreement = matched$agreement[match(rownames(profile$effects),
                                                    matched$fitted)],
        term = term,
        n_samples = nrow(projected),
        n_concepts = attr(projected, "n_concepts"),
        effect = unname(profile$effects[, term]),
        se = unname(profile$se[, term]),
        z = unname(profile$z[, term]),
        p_value = vapply(seq_along(statistic), function(i) {
          (1 + sum(null_z[, i] >= statistic[i], na.rm = TRUE)) /
            (1 + n_permutations)
        }, numeric(1)),
        stringsAsFactors = FALSE)
    }
  }

  transfer <- if (length(rows)) do.call(rbind, rows) else
    data.frame(held_out = character(), trained_on = character(),
               component = character(), fitted_component = character(),
               loading_agreement = numeric(), term = character(),
               n_samples = integer(), n_concepts = integer(),
               effect = numeric(), se = numeric(), z = numeric(),
               p_value = numeric(), stringsAsFactors = FALSE)
  rownames(transfer) <- NULL
  structure(list(transfer = transfer, n_permutations = as.integer(n_permutations),
                 alpha = control$alpha),
            class = "chorale_joint_transfer")
}

#' Project one modality's scores onto concept loadings fitted without it
#'
#' The projection uses only the concepts the held-out modality can be scored on,
#' so a loading over a concept it does not measure contributes nothing rather
#' than standing in for a measurement it never made.
#'
#' @keywords internal
#' @noRd
chorale_project_modality <- function(encoding, modality, loadings, nuisance,
                                     phenotype_column) {
  scores <- encoding$encodings[[modality]]$concept_scores
  if (is.null(scores) || ncol(scores) == 0) return(NULL)
  design <- encoding$designs[[modality]]
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  if (length(nuisance) > 0) {
    if (phenotype_column %in% nuisance) {
      rlang::abort(paste0("`nuisance` names the phenotype for modality `",
                          modality, "`."),
                   class = "chorale_nuisance_is_phenotype")
    }
    scores <- chorale_residualise_scores(scores, design, nuisance, modality)
  }
  scores <- scale(scores)
  scores[!is.finite(scores)] <- 0

  common <- intersect(colnames(scores), rownames(loadings))
  if (length(common) < ncol(loadings) + 1L) return(NULL)
  u <- loadings[common, , drop = FALSE]
  projected <- scores[, common, drop = FALSE] %*% u %*%
    solve(crossprod(u) + 1e-8 * diag(ncol(u)))
  colnames(projected) <- colnames(loadings)
  rownames(projected) <- rownames(scores)
  attr(projected, "n_concepts") <- length(common)
  projected
}

#' Match the components of one fit to those of another
#'
#' Two decompositions of overlapping data recover the same directions in a
#' different order and with arbitrary sign, since neither is identified by the
#' method. Components are therefore paired on how far their concept loadings
#' agree, over the concepts both carry, and the assignment is the one maximising
#' total agreement rather than a greedy pass, so one strong pairing cannot force
#' a poor one on what is left.
#'
#' The agreement is returned signed and unrounded. A pairing at a low absolute
#' correlation is reported rather than suppressed: it says the held-out fit
#' recovered something the reference did not, which is a result about the
#' collection and not a failure of matching.
#'
#' @param fitted,reference Loading matrices, concepts by components.
#' @returns A data frame with `fitted`, `reference` and `agreement`.
#' @keywords internal
#' @noRd
chorale_match_components <- function(fitted, reference) {
  common <- intersect(rownames(fitted), rownames(reference))
  blank <- data.frame(fitted = colnames(fitted),
                      reference = colnames(fitted),
                      agreement = NA_real_, stringsAsFactors = FALSE)
  if (length(common) < 3 || ncol(fitted) == 0 || ncol(reference) == 0) {
    return(blank)
  }
  agreement <- suppressWarnings(
    stats::cor(fitted[common, , drop = FALSE],
               reference[common, , drop = FALSE]))
  agreement[!is.finite(agreement)] <- 0
  if (all(agreement == 0)) return(blank)

  # solve_LSAP minimises, and the pairing wanted is the one maximising
  # agreement whichever way each component happens to point.
  cost <- max(abs(agreement)) - abs(agreement)
  square <- matrix(max(cost), nrow = nrow(cost),
                   ncol = max(ncol(cost), nrow(cost)))
  square[, seq_len(ncol(cost))] <- cost
  assignment <- clue::solve_LSAP(square)

  taken <- as.integer(assignment)[seq_len(nrow(cost))]
  data.frame(
    fitted = rownames(agreement),
    reference = ifelse(taken <= ncol(agreement),
                       colnames(agreement)[pmin(taken, ncol(agreement))],
                       rownames(agreement)),
    agreement = vapply(seq_along(taken), function(i) {
      if (taken[i] > ncol(agreement)) return(NA_real_)
      agreement[i, taken[i]]
    }, numeric(1)),
    stringsAsFactors = FALSE)
}

#' Restrict an encoding to a subset of its modalities
#' @keywords internal
#' @noRd
chorale_subset_encoding <- function(encoding, modalities) {
  out <- encoding
  out$modalities <- modalities
  out$encodings <- encoding$encodings[modalities]
  out$designs <- encoding$designs[modalities]
  if (!is.null(encoding$concepts$membership)) {
    out$concepts$membership <- encoding$concepts$membership[modalities]
  }
  out
}

#' @export
print.chorale_joint_transfer <- function(x, ...) {
  cat("<chorale_joint_transfer>\n")
  if (nrow(x$transfer) == 0) {
    cat("  nothing transferred\n")
    return(invisible(x))
  }
  print(x$transfer[order(x$transfer$p_value), ])
  invisible(x)
}
