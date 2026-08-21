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
                                     seed = 1L, ...) {
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
  terms <- chorale_profile_terms(spec$levels)
  phenotype_terms <- terms[startsWith(terms,
                                      paste0(spec$phenotype, "="))]
  if (length(phenotype_terms) == 0) {
    rlang::abort("The phenotype contributes no term to compare concepts on.")
  }

  observed <- chorale_concept_effects(encoding, spec, phenotype_terms,
                                      min_overlap = min_overlap)
  if (nrow(observed$joint) == 0) {
    return(structure(
      list(per_modality = observed$per_modality, joint = observed$joint,
           spec = spec, null = numeric(), n_permutations = 0L,
           alpha = control$alpha),
      class = "chorale_concept_evidence"))
  }

  # The permutation acts on the design only, so the scores computed once above
  # are the scores every permutation uses.
  null_max <- rep(NA_real_, n_permutations)
  null_by_concept <- matrix(NA_real_, nrow = n_permutations,
                            ncol = nrow(observed$joint))
  strata <- lapply(designs, function(d) {
    chorale_permutation_strata(d, spec$phenotype)
  })
  for (b in seq_len(n_permutations)) {
    permuted <- Map(function(d, s, i) {
      chorale_permute_within(d, spec$phenotype, s, seed + 1000L * b + i)
    }, designs, strata, seq_along(designs))
    e <- chorale_concept_effects(encoding, spec, phenotype_terms,
                                 designs = permuted, attribution = FALSE)
    v <- abs(e$joint$joint_z[match(observed$joint$key, e$joint$key)])
    v[!is.finite(v)] <- 0
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
      null = null_max,
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
#' @param designs Designs to use, defaulting to the encoding's own. A
#'   permutation supplies its own here, which is what lets the null reuse the
#'   scores rather than recompute them.
#' @param attribution Whether to recompute each effect with overlapping
#'   concepts regressed out. Skipped under permutation, where only the combined
#'   statistic is needed.
#' @keywords internal
#' @noRd
chorale_concept_effects <- function(encoding, spec, phenotype_terms,
                                    designs = NULL, attribution = TRUE,
                                    min_overlap = 0.1) {
  designs <- designs %||% encoding$designs
  rows <- list()
  for (m in encoding$modalities) {
    scores <- encoding$encodings[[m]]$concept_scores
    if (ncol(scores) == 0) next
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

#' Blocks of otherwise-alike samples the phenotype may be permuted within
#'
#' Permuting the phenotype across the whole cohort would break its relation to
#' the data and to the rest of the design at the same time, so the null would
#' contain designs the study could not have produced. Permuting within blocks of
#' samples alike on every other covariate leaves the design intact and breaks
#' only what is under test. Blocks are built from whatever the design carries,
#' since naming the covariates in advance would suit one study and exclude
#' another, and a covariate so finely divided that no block holds two samples is
#' dropped rather than making the permutation the identity.
#'
#' @keywords internal
#' @noRd
chorale_permutation_strata <- function(design, phenotype_column = "phenotype") {
  block_cols <- setdiff(chorale_candidate_covariates(list(design)),
                        phenotype_column)
  parts <- lapply(block_cols, function(cv) {
    v <- as.character(design[[cv]])
    v[is.na(v)] <- "unknown"
    v
  })
  if (length(parts) == 0) return(rep("all", nrow(design)))
  repeat {
    key <- do.call(paste, c(parts, sep = "|"))
    if (max(table(key)) > 1 || length(parts) == 0) break
    parts <- parts[-length(parts)]
  }
  if (length(parts) == 0) rep("all", nrow(design)) else
    do.call(paste, c(parts, sep = "|"))
}

#' Permute one column within strata
#' @keywords internal
#' @noRd
chorale_permute_within <- function(design, column, strata, seed) {
  set.seed(seed)
  for (lv in unique(strata)) {
    idx <- which(strata == lv)
    if (length(idx) > 1) design[[column]][idx] <- sample(design[[column]][idx])
  }
  design
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
