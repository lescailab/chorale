#' Whether each concept separates cases from controls
#'
#' A concept is the same concept in every modality, so nothing has to be matched
#' across modalities before it can be tested. What is tested is whether the
#' concept's score moves with the phenotype: in each modality the score is
#' regressed on the design, and the phenotype contrast is read off adjusted for
#' the covariates that modality shares with the others. A concept is supported
#' when those effects agree in sign and are large relative to their standard
#' errors in the modalities that express it.
#'
#' The modalities are combined by inverse-variance weighting, so a modality that
#' measures the concept precisely counts for more than one that measures it
#' loosely, and effects in opposite directions cancel rather than accumulate.
#' Agreement is reported beside the combined statistic, because a large combined
#' statistic carried by one modality against the others is a different finding
#' from the same statistic carried by all of them.
#'
#' Significance is calibrated by permuting the phenotype within strata of
#' otherwise-alike samples. The permutation acts on the design and leaves the
#' assay untouched, so the concept scores are computed once and reused: nothing
#' is refitted, and the number of permutations is set by the resolution the
#' report needs rather than by what a refit can afford. Two levels of error
#' control follow from the same permutations. The family-wise p-value compares
#' each concept with the largest statistic anywhere in the vocabulary, which is
#' exact under the dependence curated sets have because they overlap. The
#' q-value is the false discovery rate the same null implies at that threshold.
#'
#' Overlapping concepts are not independent claims. Each supported concept is
#' therefore also reported after the concepts it overlaps have been regressed
#' out of its score, which says whether the signal belongs to the concept or to
#' its neighbourhood.
#'
#' @param encoding A `chorale_encode` object.
#' @param n_permutations Permutations calibrating the vocabulary. The smallest
#'   attainable p-value is `1 / (n_permutations + 1)`.
#' @param min_overlap Jaccard overlap above which two concepts count as
#'   neighbours for the attribution check.
#' @param control A [chorale_control()] object.
#' @param seed Integer seed.
#' @param anchor The design covariate whose effect is tested. `NULL`, the
#'   default, tests the phenotype. Naming another covariate tests the vocabulary
#'   against that contrast instead, which is what
#'   [chorale_concept_specificity()] compares the phenotype with.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_concept_evidence` with `per_modality`,
#'   one row per concept per modality carrying the adjusted effect and its
#'   standard error; `joint`, one row per concept carrying the combined
#'   statistic, how far the modalities agree, the permutation p-values, the
#'   error control and the overlap-corrected attribution; and `spec`, the design
#'   terms the effects were adjusted for.
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
#' chorale_concept_evidence(enc, n_permutations = 99)
chorale_concept_evidence <- function(encoding, n_permutations = 999L,
                                     min_overlap = 0.1,
                                     control = chorale_control(),
                                     seed = 1L, anchor = NULL, ...) {
  if (!inherits(encoding, "chorale_encode")) {
    rlang::abort("`encoding` must be a chorale_encode object.")
  }
  control <- chorale_merge_control(control, list(...))
  modalities <- encoding$modalities
  designs <- encoding$designs

  spec <- chorale_resolve_signature(
    designs,
    phenotype_column = control$phenotype_column,
    phenotype_reference = control$phenotype_reference,
    profile_covariates = control$profile_covariates)
  anchor <- anchor %||% spec$phenotype
  anchor_terms <- chorale_anchor_terms(spec, anchor)
  if (length(anchor_terms) == 0) {
    rlang::abort(paste0("`", anchor,
                        "` contributes no term to compare concepts on."),
                 class = "chorale_unusable_anchor")
  }

  observed <- chorale_concept_effects(encoding, spec, anchor_terms,
                                      min_overlap = min_overlap)
  if (nrow(observed$joint) == 0) {
    return(structure(
      list(per_modality = observed$per_modality, joint = observed$joint,
           spec = spec, anchor = anchor, null = numeric(), n_permutations = 0L,
           alpha = control$alpha),
      class = "chorale_concept_evidence"))
  }

  # The design is held fixed and the score is what moves, so the relation
  # between the phenotype and every covariate the model adjusts for is exactly
  # the relation the data have. The encoder is never refitted: only the
  # response of the already-computed concept scores is rebuilt.
  null_max <- rep(NA_real_, n_permutations)
  null_by_concept <- matrix(NA_real_, nrow = n_permutations,
                            ncol = nrow(observed$joint))
  null_signed_by_concept <- null_by_concept
  blocks <- lapply(designs, chorale_exchangeability_blocks,
                   columns = control$exchangeability_blocks)
  for (b in seq_len(n_permutations)) {
    permuted_scores <- Map(function(m, d, block, i) {
      chorale_freedman_lane_scores(
        encoding$encodings[[m]]$concept_scores, d, spec, anchor = anchor,
        blocks = block, seed = seed + 1000L * b + i)
    }, encoding$modalities, designs, blocks, seq_along(designs))
    e <- chorale_concept_effects(encoding, spec, anchor_terms,
                                 designs = designs,
                                 scores_by_modality = permuted_scores,
                                 attribution = FALSE)
    signed <- e$joint$joint_z[match(observed$joint$key, e$joint$key)]
    signed[!is.finite(signed)] <- 0
    # The signed statistic is kept as well as its magnitude. A statistic summed
    # over a group of concepts is directional, and the sum of magnitudes with
    # signs drawn afterwards would break the dependence between overlapping
    # concepts that this null exists to carry.
    null_signed_by_concept[b, ] <- signed
    v <- abs(signed)
    null_by_concept[b, ] <- v
    null_max[b] <- max(v, na.rm = TRUE)
  }

  joint <- observed$joint
  statistic <- abs(joint$joint_z)
  joint$p_value <- vapply(seq_along(statistic), function(i) {
    (1 + sum(null_by_concept[, i] >= statistic[i], na.rm = TRUE)) /
      (1 + n_permutations)
  }, numeric(1))
  joint$p_family <- vapply(statistic, function(s) {
    (1 + sum(null_max >= s, na.rm = TRUE)) / (1 + n_permutations)
  }, numeric(1))
  joint$q_value <- chorale_permutation_fdr(statistic, null_by_concept)
  joint$significant <- joint$q_value <= control$alpha
  joint$family_significant <- joint$p_family <= control$alpha
  joint <- joint[order(-statistic), , drop = FALSE]
  rownames(joint) <- NULL

  structure(
    list(
      per_modality = observed$per_modality,
      joint = joint,
      spec = spec,
      anchor = anchor,
      null = null_max,
      # Kept so that a statistic over a group of concepts can be calibrated
      # against the same resamples, which carry the dependence overlapping
      # concepts have and a fresh set of permutations would not.
      null_by_concept = `colnames<-`(null_by_concept, observed$joint$key),
      null_signed_by_concept = `colnames<-`(null_signed_by_concept,
                                            observed$joint$key),
      n_permutations = as.integer(n_permutations),
      alpha = control$alpha,
      min_overlap = min_overlap,
      smallest_attainable_p = 1 / (n_permutations + 1)
    ),
    class = "chorale_concept_evidence"
  )
}

#' Adjusted phenotype effects on every concept, in every modality
#'
#' @param encoding A `chorale_encode` object.
#' @param spec The resolved design signature.
#' @param phenotype_terms The phenotype contrasts to read off.
#' @param designs Designs to use, defaulting to the encoding's own.
#' @param attribution Whether to recompute each effect with overlapping
#'   concepts regressed out. Skipped under permutation, where only the combined
#'   statistic is needed.
#' @param scores_by_modality Optional replacement concept scores, one matrix per
#'   modality. The null supplies its rebuilt responses here, which is what lets
#'   it hold the design fixed and never refit the encoder.
#' @keywords internal
#' @noRd
chorale_concept_effects <- function(encoding, spec, phenotype_terms,
                                    designs = NULL, attribution = TRUE,
                                    min_overlap = 0.1,
                                    scores_by_modality = NULL) {
  designs <- designs %||% encoding$designs
  rows <- list()
  for (m in encoding$modalities) {
    scores <- if (is.null(scores_by_modality)) {
      encoding$encodings[[m]]$concept_scores
    } else {
      scores_by_modality[[m]]
    }
    if (is.null(scores) || ncol(scores) == 0) next
    profile <- chorale_adjusted_profile(scores, designs[[m]], spec)
    keep <- intersect(phenotype_terms, colnames(profile$effects))
    if (length(keep) == 0) next

    adjusted <- NULL
    if (attribution) {
      residualised <- chorale_neighbour_residual(
        scores, encoding$concepts$membership[[m]], min_overlap)
      adjusted <- chorale_adjusted_profile(residualised$scores, designs[[m]],
                                           spec)
    }

    for (term in keep) {
      rows[[length(rows) + 1L]] <- data.frame(
        concept = colnames(scores),
        modality = m,
        term = term,
        n_features = as.integer(colSums(
          encoding$concepts$membership[[m]][, colnames(scores),
                                            drop = FALSE] > 0)),
        effect = unname(profile$effects[, term]),
        se = unname(profile$se[, term]),
        z = unname(profile$z[, term]),
        estimable = unname(profile$estimable[, term]),
        n_neighbours = if (is.null(adjusted)) NA_integer_ else
          residualised$n_neighbours,
        max_jaccard = if (is.null(adjusted)) NA_real_ else
          residualised$max_jaccard,
        attributed_z = if (is.null(adjusted)) NA_real_ else
          unname(adjusted$z[, term]),
        stringsAsFactors = FALSE
      )
    }
  }
  per_modality <- if (length(rows)) do.call(rbind, rows) else
    chorale_empty_per_modality()
  per_modality <- per_modality[per_modality$estimable, , drop = FALSE]
  rownames(per_modality) <- NULL

  list(per_modality = per_modality,
       joint = chorale_combine_modalities(per_modality))
}

#' @keywords internal
#' @noRd
chorale_empty_per_modality <- function() {
  data.frame(concept = character(), modality = character(), term = character(),
             n_features = integer(), effect = numeric(), se = numeric(),
             z = numeric(), estimable = logical(), n_neighbours = integer(),
             max_jaccard = numeric(), attributed_z = numeric(),
             stringsAsFactors = FALSE)
}

#' Combine one concept's effects across the modalities that express it
#'
#' Inverse-variance weighting is what puts a precisely measured modality ahead
#' of a loosely measured one without either being chosen by hand. Effects in
#' opposite directions cancel, so a concept the modalities disagree about cannot
#' reach a large combined statistic by accumulating disagreement.
#'
#' @keywords internal
#' @noRd
chorale_combine_modalities <- function(per_modality) {
  if (nrow(per_modality) == 0) {
    return(data.frame(key = character(), concept = character(),
                      term = character(), n_modalities = integer(),
                      modalities = character(), joint_effect = numeric(),
                      joint_se = numeric(), joint_z = numeric(),
                      sign_agreement = numeric(), heterogeneity_p = numeric(),
                      best_single_z = numeric(), attributed_z = numeric(),
                      max_jaccard = numeric(), stringsAsFactors = FALSE))
  }
  key <- paste(per_modality$concept, per_modality$term, sep = "|")
  parts <- split(per_modality, key)
  rows <- lapply(parts, function(d) {
    w <- 1 / d$se^2
    w[!is.finite(w)] <- 0
    total <- sum(w)
    effect <- if (total > 0) sum(w * d$effect) / total else NA_real_
    se <- if (total > 0) sqrt(1 / total) else NA_real_
    z <- if (is.finite(se) && se > 0) effect / se else NA_real_
    # Cochran's Q, the weighted spread of the effects about their combination:
    # it says whether the modalities are describing one quantity.
    q <- if (nrow(d) > 1 && total > 0) sum(w * (d$effect - effect)^2) else NA_real_
    het <- if (is.finite(q)) stats::pchisq(q, df = nrow(d) - 1,
                                           lower.tail = FALSE) else NA_real_
    attributed <- if (all(is.na(d$attributed_z))) NA_real_ else {
      wa <- w[is.finite(d$attributed_z)]
      za <- d$attributed_z[is.finite(d$attributed_z)]
      if (length(za) == 0 || sum(wa) == 0) NA_real_ else
        sum(sqrt(wa) * za) / sqrt(sum(wa))
    }
    data.frame(
      key = paste(d$concept[1], d$term[1], sep = "|"),
      concept = d$concept[1],
      term = d$term[1],
      n_modalities = nrow(d),
      modalities = paste(sort(d$modality), collapse = ", "),
      joint_effect = round(effect, 4),
      joint_se = round(se, 4),
      joint_z = round(z, 4),
      sign_agreement = round(mean(sign(d$effect) == sign(effect)), 3),
      heterogeneity_p = round(het, 4),
      best_single_z = round(max(abs(d$z), na.rm = TRUE), 4),
      attributed_z = round(attributed, 4),
      max_jaccard = round(suppressWarnings(max(d$max_jaccard, na.rm = TRUE)), 3),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$max_jaccard[!is.finite(out$max_jaccard)] <- NA_real_
  rownames(out) <- NULL
  out
}

#' Each concept's score with its overlapping neighbours regressed out
#'
#' Curated concepts share features, so a signal in one appears in every concept
#' that overlaps it. Removing the neighbours before the effect is re-estimated
#' says whether the concept carries the signal itself or is reporting the
#' neighbourhood it sits in.
#'
#' @keywords internal
#' @noRd
chorale_neighbour_residual <- function(scores, membership, min_overlap = 0.1) {
  concepts <- colnames(scores)
  n <- length(concepts)
  out <- scores
  n_neighbours <- rep(0L, n)
  max_jaccard <- rep(0, n)
  if (n < 2) {
    return(list(scores = out, n_neighbours = n_neighbours,
                max_jaccard = max_jaccard))
  }
  members <- lapply(concepts, function(cn) {
    rownames(membership)[membership[, cn] > 0]
  })
  for (i in seq_len(n)) {
    j <- vapply(seq_len(n), function(k) {
      if (k == i) return(0)
      inter <- length(intersect(members[[i]], members[[k]]))
      if (inter == 0) return(0)
      inter / length(union(members[[i]], members[[k]]))
    }, numeric(1))
    neighbours <- which(j >= min_overlap)
    n_neighbours[i] <- length(neighbours)
    max_jaccard[i] <- max(j)
    if (length(neighbours) == 0) next
    z <- scores[, neighbours, drop = FALSE]
    fit <- stats::lm.fit(cbind(1, z), scores[, i])
    out[, i] <- fit$residuals
  }
  list(scores = out, n_neighbours = n_neighbours, max_jaccard = max_jaccard)
}

#' Blocks resampling must stay inside
#'
#' Exchangeability is a property of the study rather than something a table can
#' be read for. A repeated measure, a litter or a processing batch makes samples
#' exchangeable only within itself, and only the person who ran the study knows
#' which column records that. Blocks are therefore taken from the columns
#' declared in [chorale_control()], and where none are declared every sample is
#' exchangeable with every other.
#'
#' Inferring blocks from whatever the design happens to carry is what this
#' replaces. A continuous covariate read that way produces one sample per block,
#' after which resampling is the identity or covariates are dropped until it is
#' not, and neither is a statement about the study.
#'
#' @keywords internal
#' @noRd
chorale_exchangeability_blocks <- function(design, columns = NULL) {
  if (is.null(columns) || length(columns) == 0) {
    return(rep("all", nrow(design)))
  }
  missing <- setdiff(columns, colnames(design))
  if (length(missing) > 0) {
    rlang::abort(paste0(
      "`exchangeability_blocks` names column(s) a design does not carry: ",
      paste(missing, collapse = ", "), "."))
  }
  parts <- lapply(columns, function(cv) {
    v <- as.character(design[[cv]])
    v[is.na(v)] <- "unknown"
    v
  })
  do.call(paste, c(parts, sep = "|"))
}

#' The terms a covariate contributes to the model
#'
#' A categorical covariate contributes one contrast per level beyond its
#' reference; a continuous one contributes a single standardised slope named for
#' the covariate itself. Reading the terms from the resolved signature is what
#' lets a continuous covariate be tested at all: turning one into a categorical
#' contrast would give a level per sample and leave no residual degrees of
#' freedom, and the model would fail rather than report a small effect.
#'
#' @keywords internal
#' @noRd
chorale_anchor_terms <- function(spec, anchor) {
  if (!anchor %in% names(spec$levels)) return(character())
  lv <- spec$levels[[anchor]]
  if (length(lv) == 1L && is.na(lv)) return(anchor)
  paste0(anchor, "=", lv[-1L])
}

#' Rebuild a score under the null that the phenotype has no adjusted effect
#'
#' The quantity under test is the phenotype coefficient of a model that also
#' carries the covariates the modalities share. Permuting the phenotype label
#' would break its relation to those covariates as well, so the permuted designs
#' would be designs the study could not have produced and the null would not
#' correspond to the coefficient being tested.
#'
#' Freedman and Lane (1983) permute what is left of the response once the
#' nuisance covariates have explained what they can. The reduced model, carrying
#' the intercept and every covariate except the phenotype, is fitted; its
#' residuals are permuted; and the permuted residuals are added back to the
#' reduced fit. The design matrix never changes, so the phenotype keeps exactly
#' the relation to the covariates it has in the data, and only the part of the
#' score the covariates do not explain is exchanged.
#'
#' @param scores A samples-by-concepts matrix.
#' @param design The design for those samples.
#' @param spec The resolved design signature.
#' @param blocks Resampling blocks, one label per sample.
#' @param seed Integer seed.
#' @param anchor The covariate under test, which is the one the reduced model
#'   leaves out.
#'
#' @returns A matrix of the same shape as `scores`.
#' @keywords internal
#' @noRd
chorale_freedman_lane_scores <- function(scores, design, spec, blocks, seed,
                                        anchor = spec$phenotype) {
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  reduced <- setdiff(spec$covariates, anchor)
  mm <- chorale_signature_matrix(
    design, list(covariates = spec$covariates, levels = spec$levels),
    include = reduced)
  x <- mm$x

  ok <- stats::complete.cases(x) & stats::complete.cases(scores)
  out <- scores
  if (sum(ok) <= ncol(x) || qr(x[ok, , drop = FALSE])$rank < ncol(x)) {
    # With no estimable reduced model there is nothing for the covariates to
    # explain, so the whole score is exchangeable.
    fitted <- matrix(0, nrow = nrow(scores), ncol = ncol(scores))
    residual <- scores
  } else {
    fit <- stats::lm.fit(x[ok, , drop = FALSE], scores[ok, , drop = FALSE])
    fitted <- matrix(NA_real_, nrow = nrow(scores), ncol = ncol(scores))
    residual <- matrix(NA_real_, nrow = nrow(scores), ncol = ncol(scores))
    fitted[ok, ] <- as.matrix(scores[ok, , drop = FALSE]) -
      as.matrix(fit$residuals)
    residual[ok, ] <- as.matrix(fit$residuals)
    # A sample the reduced model could not use keeps its own value, so it
    # contributes nothing to the exchange rather than dropping out of the fit.
    fitted[!ok, ] <- scores[!ok, , drop = FALSE]
    residual[!ok, ] <- 0
  }

  set.seed(seed)
  order_new <- seq_len(nrow(scores))
  for (lv in unique(blocks)) {
    at <- which(blocks == lv)
    if (length(at) > 1) order_new[at] <- sample(at)
  }
  out[] <- fitted + residual[order_new, , drop = FALSE]
  out
}

#' False discovery rate the permutation null implies at each threshold
#'
#' The expected number of concepts a null run puts above a threshold, over the
#' number the data put there. Because the permutations carry the same overlap
#' between concepts as the data, the dependence curated sets have is in the
#' estimate rather than assumed away.
#'
#' @keywords internal
#' @noRd
chorale_permutation_fdr <- function(statistic, null_by_concept) {
  n_perm <- nrow(null_by_concept)
  q <- vapply(statistic, function(s) {
    if (!is.finite(s)) return(NA_real_)
    expected_false <- sum(null_by_concept >= s, na.rm = TRUE) / n_perm
    observed <- sum(statistic >= s, na.rm = TRUE)
    if (observed == 0) return(1)
    min(1, expected_false / observed)
  }, numeric(1))
  # A concept cannot have a smaller false discovery rate than one the data
  # rank above it, so each rate is the smallest any threshold that would also
  # reject this concept attains.
  ord <- order(statistic)
  q[ord] <- cummin(q[ord])
  round(q, 4)
}

#' @export
print.chorale_concept_evidence <- function(x, ...) {
  cat("<chorale_concept_evidence>\n")
  cat("  concepts tested:", nrow(x$joint), "\n")
  cat("  permutations:", x$n_permutations,
      sprintf("(smallest attainable p-value %.3g)\n",
              1 / (x$n_permutations + 1)))
  if (nrow(x$joint) > 0) {
    cat("  supported at q <=", x$alpha, ":", sum(x$joint$significant), "\n")
    cat("  family-wise supported:", sum(x$joint$family_significant), "\n")
    top <- utils::head(x$joint, 5)
    for (i in seq_len(nrow(top))) {
      cat(sprintf("    %-40s z = %6.2f in %d modalities, q = %.3g\n",
                  substr(top$concept[i], 1, 40), top$joint_z[i],
                  top$n_modalities[i], top$q_value[i]))
    }
  }
  invisible(x)
}

#' Compare the phenotype with what a nuisance covariate reaches
#'
#' A vocabulary that answers to any strong contrast is not telling you about the
#' disease. Each nuisance covariate is tested in the phenotype's place, over the
#' same vocabulary and against the same kind of null, so the phenotype result can
#' be read against what the procedure recovers when something else drives it. A
#' vocabulary whose evidence under the phenotype is no stronger than under sex or
#' age is reporting cohort structure.
#'
#' The covariate is tested on the terms it actually contributes to the model: a
#' categorical covariate on its contrasts, a continuous one on its standardised
#' slope. Recoding a covariate into the phenotype column instead would turn a
#' continuous covariate into a contrast with one level per sample, which leaves
#' no residual degrees of freedom, and the resulting failure to estimate would
#' read as evidence that the vocabulary is specific to the phenotype.
#'
#' Anchors are compared on the family-wise p-value rather than on the largest
#' statistic. A covariate with four levels contributes three contrasts and
#' therefore three times as many chances to reach a large maximum, so raw
#' statistics favour it for a reason that has nothing to do with the biology.
#' The family-wise p-value is taken against the maximum over that covariate's
#' own family under permutation, which is what puts anchors of different sizes on
#' one scale.
#'
#' The comparison is inexpensive because the encoding does not depend on the
#' design. Concept scores are computed once and every anchor reuses them, so an
#' anchor costs the permutations and nothing else.
#'
#' @param fit A `chorale_concept_fit`.
#' @param covariates Nuisance covariates to test in the phenotype's place.
#'   `NULL` uses every covariate the modalities share, which is what the
#'   phenotype effect was adjusted for.
#' @param n_permutations Permutations calibrating each anchor. `NULL` uses the
#'   count the fit was calibrated with.
#' @param seed Integer seed.
#'
#' @returns A data frame with one row per anchor, carrying the number of terms
#'   it contributes, the strongest statistic it reaches over the vocabulary, the
#'   family-wise p-value of that statistic, the same two quantities for the
#'   phenotype, and whether the anchor reaches at least as far as the phenotype
#'   does.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
#'                            n_permutations = 99)
#' chorale_concept_specificity(fit, covariates = "sex", n_permutations = 99)
chorale_concept_specificity <- function(fit, covariates = NULL,
                                        n_permutations = NULL, seed = 1L) {
  if (!inherits(fit, "chorale_concept_fit")) {
    rlang::abort("`fit` must be a chorale_concept_fit object.")
  }
  n_perm <- n_permutations %||% fit$evidence$n_permutations
  spec <- fit$evidence$spec
  phenotype <- spec$phenotype
  covariates <- covariates %||% setdiff(spec$covariates, phenotype)

  # The phenotype is re-read against the same permutation count, so the two
  # sides of the comparison are calibrated the same way.
  reference <- if (identical(n_perm, fit$evidence$n_permutations)) {
    fit$evidence
  } else {
    chorale_concept_evidence(fit$encoding, n_permutations = n_perm,
                             control = fit$control, seed = seed)
  }
  phenotype_row <- chorale_anchor_summary(reference)

  if (length(covariates) == 0) {
    return(chorale_empty_specificity())
  }

  rows <- lapply(covariates, function(cv) {
    empty <- function(reason) {
      data.frame(anchor = cv, n_terms = NA_integer_, statistic = NA_real_,
                 p_value = NA_real_,
                 phenotype_statistic = phenotype_row$statistic,
                 phenotype_p = phenotype_row$p_value,
                 reaches_phenotype = NA, reason = reason,
                 stringsAsFactors = FALSE)
    }
    terms <- chorale_anchor_terms(spec, cv)
    if (length(terms) == 0) {
      return(empty(paste0(
        "not a covariate the modalities share on comparable terms")))
    }
    ev <- try(chorale_concept_evidence(fit$encoding, n_permutations = n_perm,
                                       control = fit$control, seed = seed,
                                       anchor = cv), silent = TRUE)
    if (inherits(ev, "try-error")) {
      return(empty("the contrast could not be estimated"))
    }
    got <- chorale_anchor_summary(ev)
    if (!is.finite(got$statistic)) {
      return(empty("no concept could be tested on this contrast"))
    }
    data.frame(
      anchor = cv,
      n_terms = length(terms),
      statistic = got$statistic,
      p_value = got$p_value,
      phenotype_statistic = phenotype_row$statistic,
      phenotype_p = phenotype_row$p_value,
      # An anchor reaches the phenotype when its calibrated evidence is at
      # least as strong, which is a comparison of p-values rather than of
      # statistics drawn from families of different size.
      reaches_phenotype = isTRUE(got$p_value <= phenotype_row$p_value),
      reason = NA_character_,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$p_value, -out$statistic), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' The strongest evidence an anchor reaches, and how calibrated it is
#' @keywords internal
#' @noRd
chorale_anchor_summary <- function(evidence) {
  j <- evidence$joint
  if (is.null(j) || nrow(j) == 0) {
    return(list(statistic = NA_real_, p_value = NA_real_))
  }
  at <- which.max(abs(j$joint_z))
  list(statistic = round(abs(j$joint_z[at]), 4),
       p_value = j$p_family[at])
}

#' @keywords internal
#' @noRd
chorale_empty_specificity <- function() {
  data.frame(anchor = character(), n_terms = integer(), statistic = numeric(),
             p_value = numeric(), phenotype_statistic = numeric(),
             phenotype_p = numeric(), reaches_phenotype = logical(),
             reason = character(), stringsAsFactors = FALSE)
}
