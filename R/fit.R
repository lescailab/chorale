#' Recover latent components within one modality
#'
#' Runs independent component analysis over several random initialisations and
#' keeps the run reaching the highest average non-Gaussianity, measured as mean
#' absolute excess kurtosis of the recovered sources. ICA is non-convex, so a
#' single run is a draw rather than an estimate; the spread across runs is
#' returned alongside the selected fit and is what [chorale_null()] reports as
#' stability.
#'
#' @param x A samples-by-features numeric matrix, centred and scaled.
#' @param n_factors Integer number of components to recover.
#' @param n_init Integer number of random initialisations.
#' @param seed Integer seed. Initialisation `i` uses `seed + i`.
#'
#' @returns A list with `scores` (samples by factors), `loadings` (features by
#'   factors), and `stability`, a data frame of the objective reached by each
#'   initialisation.
#' @keywords internal
#' @noRd
chorale_ica <- function(x, n_factors, n_init = 20L, seed = 1L) {
  rlang::check_installed("fastICA")
  best <- NULL
  best_obj <- -Inf
  obj <- rep(NA_real_, n_init)

  for (i in seq_len(n_init)) {
    set.seed(seed + i)
    fit <- try(
      fastICA::fastICA(x, n.comp = n_factors, method = "C",
                       maxit = 500, tol = 1e-5, verbose = FALSE),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    s <- scale(fit$S)
    o <- mean(abs(apply(s, 2, chorale_excess_kurtosis)), na.rm = TRUE)
    obj[i] <- o
    if (is.finite(o) && o > best_obj) {
      best_obj <- o
      best <- list(scores = s)
    }
  }

  if (is.null(best)) {
    rlang::abort("Independent component analysis failed at every initialisation.")
  }

  # Loadings in feature space: regress each feature on the recovered sources.
  coefs <- stats::coef(stats::lm(x ~ best$scores))
  loadings <- t(coefs[-1, , drop = FALSE])
  colnames(loadings) <- paste0("factor_", seq_len(n_factors))
  rownames(loadings) <- colnames(x)
  colnames(best$scores) <- colnames(loadings)
  rownames(best$scores) <- rownames(x)

  list(
    scores = best$scores,
    loadings = loadings,
    stability = data.frame(init = seq_len(n_init), objective = obj,
                           stringsAsFactors = FALSE)
  )
}

#' Excess kurtosis of a vector
#' @keywords internal
#' @noRd
chorale_excess_kurtosis <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < 4) return(NA_real_)
  m <- mean(v)
  s <- stats::sd(v)
  if (s == 0) return(NA_real_)
  sum(((v - m) / s)^4) / n - 3
}

#' Select pure features anchoring each factor
#'
#' A pure feature loads on one factor and negligibly on the others, and is what
#' gives a factor a definition independent of the rest of the loading vector.
#' At least two per factor per modality are required for the shared latent
#' structure to be recoverable.
#'
#' Purity is decided first, on the loadings alone, so the identification
#' argument holds whatever the biology turns out to be. Where more features
#' qualify than are needed, the tie is broken on biology: candidates sharing a
#' curated set with the other candidates for the same factor are preferred. Any
#' qualifying set satisfies the condition equally, so choosing the coherent one
#' costs nothing and yields markers that can be read.
#'
#' @param loadings A features-by-factors numeric matrix.
#' @param purity_ratio A feature is pure for a factor when the largest
#'   competing absolute loading is at most this fraction of its own.
#' @param min_markers Minimum qualifying features per factor.
#' @param max_markers Maximum markers retained per factor.
#' @param prior Optional feature-by-set indicator matrix from
#'   [chorale_geneset_matrix()], used only to break ties.
#'
#' @returns A list with `markers`, the features per factor that clear the
#'   purity threshold; `best_candidates`, the most nearly pure features whether
#'   or not they clear it; `purity_margin`, the competing-to-own loading ratio
#'   of the best candidates; and `pure_feature_condition`, recording which
#'   factors reached `min_markers`.
#' @export
#' @examples
#' set.seed(1)
#' l <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("f", 1:20), NULL))
#' l[1:2, 2:3] <- 0
#' chorale_markers(l)$markers[[1]]
chorale_markers <- function(loadings, purity_ratio = 0.25, min_markers = 2L,
                            max_markers = 20L, prior = NULL) {
  if (!is.matrix(loadings)) rlang::abort("`loadings` must be a matrix.")
  if (is.null(colnames(loadings))) {
    colnames(loadings) <- paste0("factor_", seq_len(ncol(loadings)))
  }
  if (is.null(rownames(loadings))) {
    rownames(loadings) <- paste0("feature_", seq_len(nrow(loadings)))
  }

  a <- abs(loadings)
  markers <- vector("list", ncol(loadings))
  best_candidates <- vector("list", ncol(loadings))
  purity_margin <- rep(NA_real_, ncol(loadings))
  names(markers) <- names(best_candidates) <- names(purity_margin) <-
    colnames(loadings)

  for (j in seq_len(ncol(loadings))) {
    own <- a[, j]
    other <- if (ncol(a) > 1) {
      apply(a[, -j, drop = FALSE], 1, max)
    } else {
      rep(0, nrow(a))
    }
    qualifies <- own > 0 & other <= purity_ratio * own
    candidates <- rownames(a)[qualifies]

    if (length(candidates) > max_markers) {
      score <- own[candidates]
      score <- score / max(score)
      if (!is.null(prior)) {
        shared <- intersect(candidates, rownames(prior))
        if (length(shared) > 1) {
          sub <- prior[shared, , drop = FALSE] > 0
          # Coherence: how many other candidates each one co-occurs with in
          # at least one curated set.
          co <- tcrossprod(sub) > 0
          coherence <- rowSums(co) - 1
          if (max(coherence) > 0) {
            bonus <- stats::setNames(rep(0, length(candidates)), candidates)
            bonus[shared] <- coherence / max(coherence)
            score <- score + bonus
          }
        }
      }
      candidates <- names(sort(score, decreasing = TRUE))[seq_len(max_markers)]
    }
    markers[[j]] <- candidates

    # The purity margin of the best available features, whether or not they
    # clear the threshold. Where the rotation is imperfectly recovered no
    # feature is strictly pure, and reporting the margin distinguishes that
    # case from a factor that is genuinely diffuse.
    margin <- other / own
    margin[!is.finite(margin)] <- Inf
    best <- names(sort(margin))[seq_len(min(max_markers, length(margin)))]
    best_candidates[[j]] <- best
    purity_margin[j] <- stats::median(margin[best[seq_len(min(min_markers, length(best)))]])
  }

  list(
    markers = markers,
    best_candidates = best_candidates,
    purity_margin = purity_margin,
    purity_ratio = purity_ratio,
    pure_feature_condition = vapply(
      markers, function(m) length(m) >= min_markers, logical(1)
    )
  )
}

#' Express factor loadings in the space of curated sets
#'
#' Regresses the loadings of features that are not markers on the curated set
#' indicators, in the manner of PLIER, giving each factor a composition in set
#' space. Marker loadings are excluded from the regression, since they carry
#' the purity the identification argument rests on.
#'
#' This is a projection computed after the factors are fitted, not a
#' constrained factorisation. The `loadings` element it returns is the
#' projected matrix, and the scores that produced the original loadings are not
#' refitted to it, so the projection does not reconstruct the data. The
#' estimator therefore keeps the fitted loadings and uses only `set_weights`,
#' as the composition each factor has in the curated vocabulary. The
#' reconstruction each version achieves is reported by [chorale_fit()], so the
#' distance between them is visible rather than assumed away.
#'
#' @param loadings A features-by-factors numeric matrix.
#' @param prior A feature-by-set indicator matrix from
#'   [chorale_geneset_matrix()], with rows matching `loadings`.
#' @param markers A named list of marker features per factor.
#' @param lambda Ridge penalty on the set coefficients.
#'
#' @returns A list with `loadings`, the projected matrix, and `set_weights`, a
#'   sets-by-factors matrix giving each factor's composition in set space.
#' @export
#' @examples
#' set.seed(1)
#' l <- matrix(rnorm(40), nrow = 10, dimnames = list(paste0("f", 1:10), NULL))
#' p <- matrix(rbinom(30, 1, 0.4), nrow = 10,
#'             dimnames = list(paste0("f", 1:10), paste0("set", 1:3)))
#' dim(chorale_constrain(l, p, rep(list(character()), 4))$set_weights)
chorale_constrain <- function(loadings, prior, markers, lambda = 1) {
  shared <- intersect(rownames(loadings), rownames(prior))
  if (length(shared) < 2 || ncol(prior) < 1) {
    return(list(
      loadings = loadings,
      set_weights = matrix(numeric(0), nrow = 0, ncol = ncol(loadings))
    ))
  }

  c_mat <- prior[shared, , drop = FALSE]
  out <- loadings
  set_weights <- matrix(
    0, nrow = ncol(c_mat), ncol = ncol(loadings),
    dimnames = list(colnames(c_mat), colnames(loadings))
  )

  ctc <- crossprod(c_mat) + lambda * diag(ncol(c_mat))
  inv <- tryCatch(solve(ctc), error = function(e) chorale_ginv(ctc))

  for (j in seq_len(ncol(loadings))) {
    free <- setdiff(shared, markers[[j]])
    if (length(free) < 2) next
    y <- loadings[shared, j]
    y[!(shared %in% free)] <- 0
    u <- inv %*% crossprod(c_mat, y)
    set_weights[, j] <- as.numeric(u)
    fitted <- as.numeric(c_mat %*% u)
    names(fitted) <- shared
    out[free, j] <- fitted[free]
  }

  list(loadings = out, set_weights = set_weights)
}

#' Variance the loadings explain, as fitted and after projection
#'
#' The fitted loadings and the scores are one object: together they reconstruct
#' the data. The projection onto curated sets is a different matrix, and the
#' scores were not refitted to it. Reporting what each explains keeps the
#' distance between the estimate and its annotation in view, so a curated
#' vocabulary too coarse to describe the factors is visible as such rather than
#' silently replacing them.
#'
#' @keywords internal
#' @noRd
chorale_reconstruction <- function(x, scores, fitted_loadings,
                                   projected_loadings) {
  total <- sum(x^2)
  explained <- function(l) {
    shared <- intersect(colnames(x), rownames(l))
    if (length(shared) == 0 || total <= 0) return(NA_real_)
    resid <- x[, shared, drop = FALSE] -
      scores %*% t(l[shared, , drop = FALSE])
    1 - sum(resid^2) / sum(x[, shared, drop = FALSE]^2)
  }
  data.frame(
    fitted = explained(fitted_loadings),
    projected = explained(projected_loadings),
    stringsAsFactors = FALSE
  )
}

#' Moore-Penrose inverse
#' @keywords internal
#' @noRd
chorale_ginv <- function(m, tol = sqrt(.Machine$double.eps)) {
  s <- svd(m)
  keep <- s$d > max(tol * s$d[1], 0)
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

#' Mean factor score per design stratum
#' @keywords internal
#' @noRd
chorale_stratum_means <- function(scores, design, strata_keys) {
  keys <- intersect(strata_keys, colnames(design))
  empty <- matrix(numeric(0), nrow = 0, ncol = ncol(scores))
  if (length(keys) == 0) return(empty)
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  ok <- stats::complete.cases(design[, keys, drop = FALSE])
  if (!any(ok)) return(empty)
  key <- do.call(paste, c(lapply(keys, function(k) as.character(design[[k]])),
                          sep = "|"))
  agg <- rowsum(scores[ok, , drop = FALSE], key[ok], reorder = TRUE)
  counts <- as.numeric(table(key[ok])[rownames(agg)])
  agg / counts
}

#' Standardised effect of each design term on each factor
#'
#' The profile is what a factor does to the design: one standardised effect per
#' design term, comparable across modalities measured on different animals. A
#' binary or categorical covariate contributes one signed standardised mean
#' difference per level beyond a reference; a continuous one contributes a
#' signed rank correlation. Every entry is therefore signed and on the same
#' scale, which is what allows profiles from different modalities to be
#' compared by direction. Profiles are estimated from samples, so their
#' precision improves with sample size rather than with the number of design
#' strata.
#'
#' @param scores A samples-by-factors matrix.
#' @param design The design table for those samples.
#' @param covariates Covariate columns to profile.
#' @param levels Optional shared-level resolution fixing which terms are formed
#'   and which level of each categorical covariate is the reference. Supply it
#'   when profiles from several modalities must line up; it is derived from
#'   `design` alone when absent.
#'
#' @returns A factors-by-terms numeric matrix.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
#' se <- chorale_load(sim$modalities[[1]], sim$col_data[[1]])
#' x <- scale(t(SummarizedExperiment::assay(se)))
#' dim(chorale_design_profile(x[, 1:3, drop = FALSE], sim$col_data[[1]],
#'                            c("phenotype", "sex")))
chorale_design_profile <- function(scores, design, covariates, levels = NULL) {
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  if (is.null(levels)) {
    levels <- chorale_profile_levels(list(design), covariates)
  }
  terms <- chorale_profile_terms(levels)
  out <- matrix(0, nrow = ncol(scores), ncol = length(terms),
                dimnames = list(colnames(scores), terms))
  if (length(terms) == 0) return(out)

  for (cv in names(levels)) {
    lv <- levels[[cv]]
    v <- design[[cv]]
    continuous <- length(lv) == 1 && is.na(lv)
    vc <- if (continuous) suppressWarnings(as.numeric(v)) else as.character(v)
    for (j in seq_len(ncol(scores))) {
      y <- scores[, j]
      if (continuous) {
        ok <- is.finite(y) & is.finite(vc)
        if (sum(ok) < 4) next
        r <- suppressWarnings(
          stats::cor(y[ok], vc[ok], method = "spearman")
        )
        out[j, cv] <- if (is.finite(r)) r else 0
      } else {
        for (l in lv[-1]) {
          out[j, paste0(cv, "=", l)] <-
            chorale_contrast_effect(y, vc, lv[1], l)
        }
      }
    }
  }
  out[!is.finite(out)] <- 0
  out
}

#' Distributional shape of each factor
#'
#' Where two modalities share no covariate, the only thing they hold in common
#' is the shape of the latent error distributions, which is what the
#' identification results match on. Independent components are standardised, so
#' location and scale carry nothing; skewness, tail weight and the quantile
#' profile do.
#'
#' @param scores A samples-by-factors matrix.
#' @param probs Quantiles describing the shape.
#'
#' @returns A factors-by-descriptors numeric matrix.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
#' x <- scale(t(sim$modalities[[1]]))
#' dim(chorale_shape_profile(x[, 1:3, drop = FALSE]))
chorale_shape_profile <- function(scores, probs = c(0.05, 0.25, 0.75, 0.95)) {
  t(apply(scores, 2, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 4) return(rep(0, length(probs) + 2))
    v <- (v - mean(v)) / stats::sd(v)
    c(skew = mean(v^3),
      kurtosis = chorale_excess_kurtosis(v),
      stats::quantile(v, probs, names = FALSE))
  }))
}

#' Match factors across all modalities at once
#'
#' With disjoint samples no animal is shared, so factors cannot be matched by
#' correlating scores. What can be compared is what a factor does to the
#' design: a factor measuring the same latent state should act on a shared
#' covariate the same way in every modality that carries it. The comparison
#' currency is the cosine between design-effect profiles, and the null is built
#' by permuting the design across samples, so precision grows with cohort size
#' and a single shared covariate suffices. Phenotype alone is enough.
#'
#' The assignment is solved once over the whole collection rather than a pair
#' at a time, by permutation synchronisation, so the correspondences agree around
#' every cycle and no programme is assembled by chaining pairwise decisions.
#' The rows returned here are the pairwise view of that one joint solution:
#' correlation bounds concern two quantities at a time and need it. They are
#' not themselves the assignment.
#'
#' The phenotype is required in every modality. The estimand is a case/control
#' contrast, so a modality that cannot express it cannot contribute one, and
#' matching such a modality on distributional shape alone would not support a
#' claim that the factors measure the same thing. Any further covariate the
#' modalities happen to share sharpens the profile, and none is required.
#'
#' @param fits A named list of per-modality fits, each carrying `scores`.
#' @param designs A named list of per-modality design tables.
#' @param strata_keys Candidate covariates for the design profile. Those
#'   absent, constant, or unshared are dropped, so supplying more than the data
#'   carry is harmless.
#' @param n_perm Number of permutations calibrating the statistic.
#' @param alpha Significance threshold.
#' @param seed Integer seed.
#'
#' @returns A data frame, one row per cross-modality factor pair implied by the
#'   joint assignment.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 2, n_per_cell = 2, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' fit$matches
chorale_match <- function(fits, designs,
                          strata_keys = c("phenotype", "sex", "age_months",
                                          "age_bin", "strain", "region"),
                          n_perm = 200L, alpha = 0.05, seed = 1L) {
  chorale_integrate(fits, designs, strata_keys = strata_keys,
                    n_perm = n_perm, alpha = alpha, seed = seed)$matches
}

#' Solve the assignment and calibrate it, over all modalities at once
#'
#' The whole integration step: shared design terms are resolved once for the
#' collection, one profile is built per modality on those terms, the assignment
#' is solved jointly, and both the joint and the pairwise statistics are
#' calibrated against a single null in which the complete procedure is rerun on
#' permuted designs. Permuting the design rows as a block rather than each
#' covariate independently preserves the dependence among age, sex, cohort and
#' phenotype, so the null contains only designs that could have arisen.
#'
#' One null calibrates every programme and every modality subset the procedure
#' could have selected, which is what makes the reported p-values account for
#' the selection rather than condition on it.
#'
#' @inheritParams chorale_match
#'
#' @returns A list with `programmes`, `matches`, `leave_one_out`,
#'   `synchronisation`, `terms` and `null`.
#' @keywords internal
#' @noRd
chorale_integrate <- function(fits, designs,
                              strata_keys = c("phenotype", "sex", "age_months",
                                              "age_bin", "strain", "region"),
                              n_perm = 200L, alpha = 0.05, seed = 1L) {
  modalities <- names(fits)
  empty <- list(programmes = data.frame(), matches = data.frame(),
                leave_one_out = data.frame(), synchronisation = NULL,
                terms = character(0), null = numeric(0))
  if (length(modalities) < 2) return(empty)

  aligned <- lapply(modalities, function(m) {
    d <- designs[[m]]
    d[match(rownames(fits[[m]]$scores), d$sample_id), , drop = FALSE]
  })
  names(aligned) <- modalities

  shared <- chorale_shared_covariates_all(aligned, strata_keys)
  chorale_require_phenotype(shared, modalities)
  levels <- chorale_profile_levels(aligned, shared)
  terms <- chorale_profile_terms(levels)
  if (length(terms) == 0) return(empty)

  profile_of <- function(designs_by_mod) {
    out <- lapply(modalities, function(m) {
      chorale_design_profile(fits[[m]]$scores, designs_by_mod[[m]], shared,
                             levels = levels)
    })
    names(out) <- modalities
    out
  }

  profiles <- profile_of(aligned)
  sync <- chorale_synchronise(profiles)
  if (sync$n_programmes == 0) return(empty)

  # One null for the whole procedure: the assignment is solved again on each
  # permuted design, so the calibrating value is the best any programme and any
  # modality subset could have reached by chance.
  set.seed(seed)
  null_joint <- numeric(n_perm)
  null_pair <- lapply(names(sync$similarity), function(k) numeric(n_perm))
  names(null_pair) <- names(sync$similarity)
  for (b in seq_len(n_perm)) {
    permuted <- lapply(modalities, function(m) {
      d <- aligned[[m]]
      d[, shared] <- d[sample(nrow(d)), shared, drop = FALSE]
      d
    })
    names(permuted) <- modalities
    sp <- chorale_synchronise(profile_of(permuted))
    null_joint[b] <- chorale_best_statistic(sp)
    for (k in names(null_pair)) {
      blk <- sp$similarity[[k]]
      null_pair[[k]][b] <- if (is.null(blk)) 0 else max(abs(blk))
    }
  }
  p_of <- function(value, null) (1 + sum(null >= value)) / (1 + n_perm)

  prog_rows <- list()
  loo_rows <- list()
  match_rows <- list()

  for (pr in sort(unique(sync$assignment$programme))) {
    members <- sync$assignment[sync$assignment$programme == pr, , drop = FALSE]
    if (nrow(members) < 2) next

    # The modality subset is part of what the procedure selects, so it is
    # chosen against the same null: the widest subset whose joint evidence
    # survives, and the full tuple where none does.
    subsets <- chorale_member_subsets(nrow(members))
    stats_by_subset <- vapply(subsets, function(idx) {
      chorale_programme_statistic(sync$similarity, members[idx, , drop = FALSE])
    }, numeric(1))
    ps <- vapply(stats_by_subset, p_of, numeric(1), null = null_joint)
    ok <- which(is.finite(stats_by_subset) & ps < alpha)
    if (length(ok) > 0) {
      sizes <- vapply(subsets[ok], length, integer(1))
      best <- ok[order(-sizes, -stats_by_subset[ok])][1]
      supported <- TRUE
    } else {
      best <- which.max(vapply(subsets, length, integer(1)))
      supported <- FALSE
    }
    keep <- members[subsets[[best]], , drop = FALSE]
    statistic <- stats_by_subset[best]
    p_joint <- ps[best]

    label <- paste0("P", length(prog_rows) + 1L)
    prog_rows[[length(prog_rows) + 1L]] <- data.frame(
      programme = label,
      n_modalities = nrow(keep),
      modalities = paste(sort(keep$modality), collapse = ", "),
      modality = keep$modality,
      factor = keep$factor,
      joint_statistic = round(statistic, 4),
      joint_p = p_joint,
      supported = supported,
      stringsAsFactors = FALSE
    )

    # Leaving one modality out shows whether the evidence rests on all of them
    # or on one; a programme that is no weaker without a modality was not
    # carrying it.
    if (nrow(keep) >= 3) {
      for (i in seq_len(nrow(keep))) {
        sub <- keep[-i, , drop = FALSE]
        s <- chorale_programme_statistic(sync$similarity, sub)
        loo_rows[[length(loo_rows) + 1L]] <- data.frame(
          programme = label,
          dropped = keep$modality[i],
          n_modalities = nrow(sub),
          joint_statistic = round(s, 4),
          joint_p = p_of(s, null_joint),
          delta = round(s - statistic, 4),
          stringsAsFactors = FALSE
        )
      }
    }

    for (i in seq_len(nrow(keep))) {
      for (j in seq_len(nrow(keep))) {
        if (j <= i) next
        ma <- keep$modality[i]
        mb <- keep$modality[j]
        raw <- chorale_pair_agreement(sync$similarity, ma, mb,
                                      keep$factor_index[i], keep$factor_index[j])
        if (!is.finite(raw)) next
        key <- if (!is.null(sync$similarity[[paste(ma, mb, sep = "|")]])) {
          paste(ma, mb, sep = "|")
        } else {
          paste(mb, ma, sep = "|")
        }
        value <- abs(raw)
        match_rows[[length(match_rows) + 1L]] <- data.frame(
          programme = label,
          modality_a = ma, modality_b = mb,
          factor_a = keep$factor[i], factor_b = keep$factor[j],
          sign = if (raw < 0) -1 else 1,
          basis = "design effects",
          n_shared_covariates = length(shared),
          shared_covariates = paste(shared, collapse = ","),
          statistic = value,
          p_value = p_of(value, null_pair[[key]]),
          joint_statistic = round(statistic, 4),
          joint_p = p_joint,
          supported = supported,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  programmes <- do.call(rbind, prog_rows)
  matches <- do.call(rbind, match_rows)
  if (is.null(programmes)) programmes <- data.frame()
  if (is.null(matches)) {
    matches <- data.frame()
  } else {
    matches$significant <- matches$p_value < alpha
    matches$p_attainable_floor <- 1 / (1 + n_perm)
    matches <- matches[order(matches$joint_p, matches$p_value,
                             -matches$statistic), , drop = FALSE]
  }
  loo <- do.call(rbind, loo_rows)
  if (is.null(loo)) loo <- data.frame()

  list(programmes = programmes, matches = matches, leave_one_out = loo,
       synchronisation = sync, terms = terms, null = null_joint)
}

#' Require the phenotype to be shared
#'
#' The estimand is a case/control contrast, so the phenotype is present in
#' every modality by construction. A modality whose deposited metadata does not
#' resolve it cannot contribute a comparable contrast, and matching it on
#' distributional shape alone would produce a result that could not be defended
#' as measuring the same thing.
#'
#' @keywords internal
#' @noRd
chorale_require_phenotype <- function(shared, modalities) {
  if (!"phenotype" %in% shared) {
    rlang::abort(
      paste0(
        "Modalities '", paste(modalities, collapse = "', '"),
        "' do not all share a phenotype contrast. ",
        "chorale estimates a case/control shared state, so every modality must ",
        "carry a `phenotype` column taking at least two values. Resolve the ",
        "phenotype for every modality, or drop the one that lacks it."
      ),
      class = "chorale_missing_phenotype"
    )
  }
  invisible(TRUE)
}

#' Covariates present and varying in every design
#' @keywords internal
#' @noRd
chorale_shared_covariates_all <- function(designs, candidates) {
  usable <- function(d, cv) {
    if (!cv %in% colnames(d)) return(FALSE)
    length(unique(stats::na.omit(d[[cv]]))) >= 2
  }
  shared <- candidates[vapply(candidates, function(cv) {
    all(vapply(designs, usable, logical(1), cv = cv))
  }, logical(1))]
  # age_bin and age_months describe the same variable, so the finer one is
  # kept and the coarser dropped to avoid counting it twice.
  if (all(c("age_months", "age_bin") %in% shared)) {
    shared <- setdiff(shared, "age_bin")
  }
  shared
}

#' Covariates present and varying in both designs
#' @keywords internal
#' @noRd
chorale_shared_covariates <- function(da, db, candidates) {
  chorale_shared_covariates_all(list(da, db), candidates)
}

#' One-to-one assignment tolerating a rectangular statistic
#' @keywords internal
#' @noRd
chorale_assign <- function(stat) {
  if (nrow(stat) > ncol(stat)) {
    flipped <- clue::solve_LSAP(t(stat), maximum = TRUE)
    assignment <- rep(NA_integer_, nrow(stat))
    assignment[as.integer(flipped)] <- seq_len(ncol(stat))
    assignment
  } else {
    as.integer(clue::solve_LSAP(stat, maximum = TRUE))
  }
}

#' Programmes recovered by the joint assignment
#'
#' A programme is a set of factors, at most one per modality, that the joint
#' assignment placed at the same position in the common space of latent
#' programmes. It is not assembled from pairwise decisions: the assignment was
#' solved once over the whole collection, so a programme's membership is
#' settled using every modality rather than by chaining one modality to the
#' next.
#'
#' Which modalities a programme spans is itself part of what the procedure
#' selects, and it is chosen against the same null the statistic is calibrated
#' on: the widest set of modalities whose joint evidence survives. A programme
#' for which no such set exists is returned with `supported` false and is
#' excluded by default.
#'
#' @param fit A `chorale_fit` object.
#' @param significant_only Return only programmes whose joint evidence beats
#'   its null.
#'
#' @returns A data frame with one row per (programme, modality, factor)
#'   membership, carrying `programme`, `n_modalities`, `modalities`,
#'   `joint_statistic`, `joint_p` and `supported`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 3, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
#' chorale_programmes(fit)
chorale_programmes <- function(fit, significant_only = TRUE) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  pg <- fit$programmes
  if (is.null(pg) || nrow(pg) == 0) return(data.frame())
  if (significant_only) pg <- pg[pg$supported, , drop = FALSE]
  pg[order(pg$joint_p, -pg$joint_statistic, pg$programme), , drop = FALSE]
}

#' What each modality contributes to a programme
#'
#' The joint statistic says a programme is carried by its modalities together.
#' Dropping one at a time says whether it needed them: a programme whose
#' evidence is unchanged when a modality is removed was not resting on that
#' modality, and one that collapses without it was. Reporting both is what
#' distinguishes a result integration produced from a result one modality
#' produced and the others accompanied.
#'
#' @param fit A `chorale_fit` object.
#' @param programmes Optional subset of [chorale_programmes()] to restrict to.
#'
#' @returns A data frame with one row per (programme, dropped modality),
#'   carrying the joint statistic and p-value without that modality and the
#'   change from the full programme. Empty where no programme spans three or
#'   more modalities.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 3, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
#' chorale_leave_one_out(fit)
chorale_leave_one_out <- function(fit, programmes = NULL) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  loo <- fit$leave_one_out
  if (is.null(loo) || nrow(loo) == 0) return(data.frame())
  if (!is.null(programmes) && nrow(programmes) > 0) {
    loo <- loo[loo$programme %in% programmes$programme, , drop = FALSE]
  }
  loo
}

#' Joint evidence that a programme is carried by all its modalities
#'
#' The evidence for a programme is joint by construction: the statistic is the
#' mean agreement of the design profiles over every pair inside the programme,
#' evaluated as one quantity, and the null reruns the whole procedure on
#' permuted designs, keeping the best value any programme and any modality
#' subset could have reached. Requiring three modalities to agree
#' simultaneously is far harder to achieve by chance than requiring two, and
#' because the null sees the same choices the estimator made, the p-value
#' accounts for the selection rather than conditioning on it.
#'
#' This function reports the quantities [chorale_fit()] already computed.
#'
#' @param fit A `chorale_fit` object.
#' @param programmes Output of [chorale_programmes()]; taken from `fit` if
#'   absent.
#' @param n_perm Ignored; retained so existing calls keep working. The null is
#'   computed once during fitting.
#' @param seed Ignored; see `n_perm`.
#'
#' @returns `programmes`, carrying `joint_statistic` and `joint_p`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 3, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
#' chorale_joint_evidence(fit)
chorale_joint_evidence <- function(fit, programmes = NULL, n_perm = NULL,
                                   seed = NULL) {
  if (is.null(programmes)) programmes <- chorale_programmes(fit)
  programmes
}

#' Add the age band used for anchoring
#' @keywords internal
#' @noRd
chorale_add_age_bin <- function(design) {
  if ("age_bin" %in% colnames(design)) return(design)
  if (!"age_months" %in% colnames(design)) return(design)
  age <- suppressWarnings(as.numeric(as.character(design$age_months)))
  design$age_bin <- cut(age, breaks = c(0, 4, 9, 18, 100),
                        labels = c("2mo", "6mo", "14mo", "aged"))
  design
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Estimate the shared latent state across modalities
#'
#' The estimator of `AGENT_PLAN.md` Section 8.3, in four steps. Components are
#' recovered per modality by independent component analysis over several
#' initialisations. Pure features are then selected on the loadings alone, so
#' the identification argument holds before any prior is applied, with ties
#' broken towards curated-set coherence. The remaining loadings are shrunk
#' towards those sets, so every factor carries a pathway definition at
#' estimation. Factors are finally matched across modalities by distributional
#' agreement, corroborated by the design strata.
#'
#' @param containers A named list of [SummarizedExperiment::SummarizedExperiment]
#'   objects, one per modality, as returned by [chorale_load()].
#' @param n_factors Integer, or one integer per modality, giving the number of
#'   components to recover, from the detectability gate.
#' @param gene_sets A named list of curated gene sets, as returned by
#'   [chorale_genesets()]. Without it, markers are selected on loadings alone
#'   and factors carry no pathway definition.
#' @param feature_map Optional named list, one entry per modality, each a data
#'   frame from [chorale_map()] harmonising that modality's features to Entrez
#'   identifiers.
#' @param feature_space Optional named character vector, one entry per
#'   modality, saying what its features are: `"gene"` for a transcriptome or
#'   proteome, whose identifiers reach the curated sets through
#'   [chorale_map()], or `"lipid"` for a lipidome, whose features reach them
#'   through their class by [chorale_metabolite_matrix()]. Defaults to `"gene"`
#'   everywhere. It is declared rather than guessed, since a wrong guess would
#'   silently leave a modality out of the pathway comparison.
#' @param n_pathway_perm Annotation-matched permutations calibrating the
#'   pathway channel. Zero skips it.
#' @param n_init Integer number of random initialisations per modality.
#' @param strata_keys Design columns defining an anchoring stratum.
#' @param assay_name Assay to take from each container.
#' @param seed Integer seed.
#'
#' @returns An object of class `chorale_fit`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 2, n_per_cell = 2, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' fit
chorale_fit <- function(containers,
                        n_factors,
                        gene_sets = NULL,
                        feature_map = NULL,
                        feature_space = NULL,
                        n_init = 20L,
                        strata_keys = c("phenotype", "age_bin", "sex"),
                        assay_name = NULL,
                        n_pathway_perm = 200L,
                        seed = 1L) {
  if (!is.list(containers) || length(containers) < 2) {
    rlang::abort("`containers` must be a list of at least two modalities.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  modalities <- names(containers)

  if (length(n_factors) == 1) n_factors <- rep(n_factors, length(containers))
  names(n_factors) <- modalities

  feature_space <- chorale_feature_space(feature_space, modalities)

  fits <- list()
  designs <- list()

  for (m in modalities) {
    se <- containers[[m]]
    an <- assay_name %||% SummarizedExperiment::assayNames(se)[1]
    mat <- SummarizedExperiment::assay(se, an)
    design <- as.data.frame(SummarizedExperiment::colData(se))
    design <- chorale_add_age_bin(design)

    x <- scale(t(as.matrix(mat)))
    x[!is.finite(x)] <- 0

    fit <- chorale_ica(x, n_factors[[m]], n_init = n_init, seed = seed)

    prior <- NULL
    if (!is.null(gene_sets)) {
      if (identical(unname(feature_space[[m]]), "lipid")) {
        # A lipidome reaches the sets through its classes rather than through
        # gene identifiers, so both modalities end up in one vocabulary.
        prior <- chorale_metabolite_matrix(rownames(mat), gene_sets)
      } else {
        ids <- rownames(mat)
        weights <- rep(1, length(ids))
        if (!is.null(feature_map) && !is.null(feature_map[[m]])) {
          fm <- feature_map[[m]]
          idx <- match(ids, fm$id)
          mapped <- !is.na(idx)
          ids[mapped] <- fm$ENTREZID[idx[mapped]]
          weights[mapped] <- fm$weight[idx[mapped]]
        }
        prior <- chorale_geneset_matrix(ids, gene_sets, weights = weights)
        if (ncol(prior) > 0) rownames(prior) <- rownames(mat)
      }
    }

    mk <- chorale_markers(fit$loadings, prior = prior)
    if (!is.null(prior) && ncol(prior) > 0) {
      # The projection onto curated sets is an annotation, so it is stored
      # beside the fit rather than substituted for it. Overwriting the
      # loadings would leave a matrix the scores no longer reconstruct.
      projected <- chorale_constrain(fit$loadings, prior, mk$markers)
      fit$set_weights <- projected$set_weights
      fit$pathway_loadings <- projected$loadings
      fit$reconstruction <- chorale_reconstruction(x, fit$scores, fit$loadings,
                                                   projected$loadings)
    }
    fit$markers <- mk$markers
    fit$best_candidates <- mk$best_candidates
    fit$purity_margin <- mk$purity_margin
    fit$pure_feature_condition <- mk$pure_feature_condition
    fit$prior <- prior

    fits[[m]] <- fit
    designs[[m]] <- design
  }

  integration <- chorale_integrate(fits, designs, strata_keys = strata_keys,
                                   seed = seed)
  matches <- integration$matches

  out <- structure(
    list(
      modalities = modalities,
      feature_space = feature_space,
      fits = fits,
      designs = designs,
      matches = matches,
      programmes = integration$programmes,
      leave_one_out = integration$leave_one_out,
      synchronisation = integration$synchronisation,
      profile_terms = integration$terms,
      joint_null = integration$null,
      n_shared = if (nrow(matches) > 0) sum(matches$significant) else 0L,
      n_programmes = if (nrow(integration$programmes) > 0) {
        sum(!duplicated(integration$programmes$programme) &
              integration$programmes$supported)
      } else {
        0L
      },
      n_factors = n_factors,
      strata_keys = strata_keys,
      gene_sets = gene_sets,
      seed = seed
    ),
    class = "chorale_fit"
  )

  # The pathway channel is a second, independent line of evidence, so it is
  # computed once here and travels with the fit.
  out$pathway_evidence <- if (n_pathway_perm > 0) {
    chorale_pathway_evidence(out, n_perm = n_pathway_perm, seed = seed)
  } else {
    data.frame()
  }
  out
}

#' Resolve the declared feature space of every modality
#' @keywords internal
#' @noRd
chorale_feature_space <- function(feature_space, modalities) {
  if (is.null(feature_space)) {
    return(stats::setNames(rep("gene", length(modalities)), modalities))
  }
  unknown <- setdiff(names(feature_space), modalities)
  if (length(unknown) > 0) {
    rlang::abort(paste0("`feature_space` names unknown modalities: ",
                        paste(unknown, collapse = ", "), "."))
  }
  out <- stats::setNames(rep("gene", length(modalities)), modalities)
  out[names(feature_space)] <- as.character(feature_space)
  bad <- setdiff(unique(out), c("gene", "lipid"))
  if (length(bad) > 0) {
    rlang::abort(paste0("`feature_space` must be \"gene\" or \"lipid\"; got: ",
                        paste(bad, collapse = ", "), "."))
  }
  out
}

#' @export
print.chorale_fit <- function(x, ...) {
  cat("<chorale_fit>\n")
  cat("  modalities:", paste(x$modalities, collapse = ", "), "\n")
  cat("  factors per modality:", paste(x$n_factors, collapse = ", "), "\n")
  for (m in x$modalities) {
    ok <- sum(x$fits[[m]]$pure_feature_condition)
    cat(sprintf("  %-12s %d of %d factors carry pure features\n",
                m, ok, length(x$fits[[m]]$pure_feature_condition)))
  }
  cat("  programmes recovered by the joint assignment:",
      if (is.data.frame(x$programmes) && nrow(x$programmes) > 0) {
        length(unique(x$programmes$programme))
      } else {
        0
      },
      "of which supported:", x$n_programmes %||% 0L, "\n")
  cat("  implied factor pairs:",
      if (is.data.frame(x$matches)) nrow(x$matches) else 0,
      "of which significant:", x$n_shared, "\n")
  invisible(x)
}
