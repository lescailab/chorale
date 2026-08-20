#' Recover latent components within one modality
#'
#' Runs independent component analysis over several random initialisations and
#' keeps the medoid run: the valid run with the greatest mean matched-factor
#' agreement with every other valid run. This retains an actual ICA solution
#' while choosing the most reproducible one. ICA is non-convex, so the spread
#' across runs is returned alongside the selected fit.
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
chorale_ica <- function(x, n_factors, n_init = 20L, seed = 1L,
                        consensus = TRUE) {
  rlang::check_installed("fastICA")
  best <- NULL
  best_obj <- -Inf
  obj <- rep(NA_real_, n_init)
  runs <- vector("list", n_init)

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
    runs[[i]] <- s
    if (is.finite(o) && o > best_obj) {
      best_obj <- o
      best <- list(scores = s)
    }
  }

  if (is.null(best)) {
    rlang::abort("Independent component analysis failed at every initialisation.")
  }

  # Averaging aligned runs does not generally produce another ICA solution.
  # Select the stability medoid instead: an observed run whose factors are most
  # reproducible across starts. `consensus` is retained for API compatibility;
  # both settings now return a valid representative run.
  valid <- which(!vapply(runs, is.null, logical(1)))
  if (length(valid) > 1L) {
    agreement <- matrix(NA_real_, length(valid), length(valid))
    diag(agreement) <- 1
    for (a in seq_along(valid)) for (b in seq_along(valid)) {
      if (b <= a) next
      cor_mat <- abs(suppressWarnings(stats::cor(runs[[valid[a]]],
                                                   runs[[valid[b]]])))
      cor_mat[!is.finite(cor_mat)] <- 0
      assignment <- clue::solve_LSAP(cor_mat, maximum = TRUE)
      value <- mean(cor_mat[cbind(seq_len(nrow(cor_mat)),
                                  as.integer(assignment))])
      agreement[a, b] <- agreement[b, a] <- value
    }
    medoid <- valid[which.max(rowMeans(agreement, na.rm = TRUE))]
    best$scores <- runs[[medoid]]
    best$selected_init <- medoid
  } else {
    best$selected_init <- valid[1]
  }

  # Loadings in feature space: regress each feature on the recovered sources.
  coefs <- stats::coef(stats::lm(x ~ best$scores))
  loadings <- t(coefs[-1, , drop = FALSE])
  colnames(loadings) <- paste0("factor_", seq_len(n_factors))
  rownames(loadings) <- colnames(x)
  colnames(best$scores) <- colnames(loadings)
  rownames(best$scores) <- rownames(x)

  # Stability is a statement about the recovered factors, not the objective:
  # each other run's factors are matched to the selected run's, and the mean
  # matched correlation says how reproducibly the same subspace is recovered. A
  # stable objective can coexist with rotated or permuted factors.
  subspace <- chorale_subspace_stability(best$scores, runs)

  list(
    scores = best$scores,
    loadings = loadings,
    stability = data.frame(init = seq_len(n_init), objective = obj,
                           subspace = subspace,
                           selected = seq_len(n_init) == best$selected_init,
                           stringsAsFactors = FALSE),
    selected_init = best$selected_init
  )
}

#' Match every run's factors to a reference and align their signs
#'
#' A factor recovered at one initialisation is the same factor recovered at
#' another up to order and sign, which are not identified. Matching one to one
#' by absolute correlation and flipping the sign to agree puts the runs in a
#' common frame, which is what allows them to be averaged.
#'
#' @keywords internal
#' @noRd
chorale_align_runs <- function(reference, runs) {
  out <- list()
  for (s in runs) {
    if (is.null(s) || !identical(dim(s), dim(reference))) next
    a <- suppressWarnings(stats::cor(reference, s))
    a[!is.finite(a)] <- 0
    assign <- as.integer(clue::solve_LSAP(abs(a), maximum = TRUE))
    matched <- s[, assign, drop = FALSE]
    signs <- sign(a[cbind(seq_len(nrow(a)), assign)])
    signs[signs == 0] <- 1
    out[[length(out) + 1L]] <- sweep(matched, 2, signs, `*`)
  }
  out
}

#' Reproducibility of the recovered factors across initialisations
#'
#' For each other run, its factors are matched one-to-one to the selected run's
#' by absolute correlation, and the mean matched correlation is that run's
#' agreement. A value near one means the same factors were recovered whatever
#' the start; a low value means the recovery is a draw. `NA` for the selected
#' run itself and for runs that failed.
#'
#' @keywords internal
#' @noRd
chorale_subspace_stability <- function(reference, runs) {
  vapply(runs, function(s) {
    if (is.null(s) || !identical(dim(s), dim(reference))) return(NA_real_)
    a <- abs(suppressWarnings(stats::cor(reference, s)))
    a[!is.finite(a)] <- 0
    if (nrow(a) != ncol(a)) return(NA_real_)
    assign <- clue::solve_LSAP(a, maximum = TRUE)
    mean(a[cbind(seq_len(nrow(a)), as.integer(assign))])
  }, numeric(1))
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
#' A pure feature loads strongly on one factor and weakly on the others. This is
#' a loading-concentration diagnostic; CHORALE does not use it as proof that a
#' latent factor is uniquely recovered.
#'
#' Purity is calculated first, on the loadings alone. Where more features
#' qualify than are retained, the tie is broken on biology: candidates sharing a
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
#' space. Marker loadings are excluded so the same features are not used both
#' to define the loading-purity diagnostic and to fit the set projection.
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

#' Adjusted standardised effect of each design term on each factor
#'
#' The profile is estimated from one multivariable model per factor. A binary or
#' categorical covariate contributes adjusted contrasts beyond a declared
#' reference; a continuous covariate contributes an adjusted slope per standard
#' deviation. Factor scores are standardised before this function is called, so
#' the coefficients share a factor-score scale. Standard errors, covariance and
#' estimability are retained as attributes for uncertainty-aware matching.
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
  spec <- list(covariates = names(levels), levels = levels,
               phenotype = names(levels)[1],
               secondary = names(levels)[-1])
  profile <- chorale_adjusted_profile(scores, design, spec)
  out <- profile$effects
  attr(out, "standard_errors") <- profile$se
  attr(out, "covariance") <- profile$covariance
  attr(out, "estimable") <- profile$estimable
  attr(out, "term_covariate") <- profile$term_covariate
  out
}

#' Distributional shape of each factor
#'
#' Describes skewness, tail weight and selected quantiles of each standardised
#' factor. This is a diagnostic only: current matching requires phenotype and
#' does not use distributional shape to identify or pair factors.
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
#' design: compatible factors should have compatible adjusted phenotype
#' effects. Covariance-aware quadratic losses define phenotype-compatible
#' alternatives; other eligible shared covariates can refine only those
#' alternatives. The null uses reduced-model residual permutations, so
#' phenotype alone is sufficient and nuisance covariates cannot carry a match.
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
#' @param strata_keys Candidate covariates for the design profile. `NULL`, the
#'   default, considers every covariate the designs carry, so a design the
#'   package has not seen anchors without naming its columns in advance. Those
#'   absent, constant, or unshared are dropped, so supplying more than the data
#'   carry is harmless.
#' @param n_perm Number of permutations calibrating the statistic.
#' @param alpha Significance threshold.
#' @param seed Integer seed.
#' @param phenotype_column Mandatory phenotype column.
#' @param phenotype_reference Reference level for the phenotype contrast.
#' @param profile_covariates Optional secondary shared covariates. `NULL`
#'   discovers all eligible shared covariates.
#' @param exchangeability_blocks Optional columns restricting permutations.
#'
#' @returns A data frame, one row per cross-modality factor pair implied by the
#'   joint assignment.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 2, n_per_cell = 2, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' fit$matches
chorale_match <- function(fits, designs,
                          strata_keys = NULL,
                          n_perm = 200L, alpha = 0.05, seed = 1L,
                          phenotype_column = "phenotype",
                          phenotype_reference = "control",
                          profile_covariates = NULL,
                          exchangeability_blocks = NULL) {
  chorale_integrate(
    fits, designs, strata_keys = strata_keys,
    n_perm = n_perm, alpha = alpha, seed = seed,
    phenotype_column = phenotype_column,
    phenotype_reference = phenotype_reference,
    profile_covariates = profile_covariates,
    exchangeability_blocks = exchangeability_blocks)$matches
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
                              strata_keys = NULL,
                              n_perm = 200L, alpha = 0.05, seed = 1L,
                              phenotype_column = "phenotype",
                              phenotype_reference = "control",
                              profile_covariates = NULL,
                              exchangeability_blocks = NULL,
                              ambiguity_level = 0.95,
                              n_ambiguity_boot = 0L,
                              n_cores = 1L) {
  modalities <- names(fits)
  empty <- list(programmes = data.frame(), matches = data.frame(),
                leave_one_out = data.frame(), synchronisation = NULL,
                terms = character(0), null = numeric(0),
                strata_keys = character(0), signature = NULL,
                excluded_covariates = data.frame(), profiles = list())
  if (length(modalities) < 2) return(empty)

  aligned <- lapply(modalities, function(m) {
    d <- designs[[m]]
    d[match(rownames(fits[[m]]$scores), d$sample_id), , drop = FALSE]
  })
  names(aligned) <- modalities

  # `strata_keys` used to select both matching covariates and bound strata. It
  # remains a compatibility alias for profile covariates, while new fits store
  # the two concepts separately.
  candidates <- profile_covariates %||% strata_keys
  spec <- chorale_resolve_signature(
    aligned, phenotype_column = phenotype_column,
    phenotype_reference = phenotype_reference,
    profile_covariates = candidates)

  profile_of <- function(scores_by_mod) {
    out <- lapply(modalities, function(m) {
      chorale_adjusted_profile(scores_by_mod[[m]], aligned[[m]], spec)
    })
    names(out) <- modalities
    out
  }

  observed_scores <- lapply(fits, `[[`, "scores")
  profiles <- profile_of(observed_scores)
  blocks <- chorale_hierarchical_blocks(profiles, spec, ambiguity_level)
  phenotype_terms <- blocks$phenotype_terms
  task_seeds <- chorale_resampling_seeds(
    seed, c(n_ambiguity_boot, n_perm, n_perm))
  ambiguity_seeds <- task_seeds[[1L]]
  secondary_seeds <- task_seeds[[2L]]
  phenotype_seeds <- task_seeds[[3L]]
  if (n_ambiguity_boot > 0L) {
    boot_blocks <- chorale_deterministic_lapply(seq_len(n_ambiguity_boot), function(b) {
      set.seed(ambiguity_seeds[b])
      boot_profiles <- lapply(modalities, function(m) {
        at <- chorale_bootstrap_rows(aligned[[m]], exchangeability_blocks)
        d <- aligned[[m]][at, , drop = FALSE]
        s <- observed_scores[[m]][at, , drop = FALSE]
        d$sample_id <- rownames(s) <- paste0("bootstrap_", seq_len(nrow(d)))
        chorale_adjusted_profile(s, d, spec)
      })
      names(boot_profiles) <- modalities
      chorale_hierarchical_blocks(boot_profiles, spec, ambiguity_level)
    }, n_cores)
    blocks <- chorale_bootstrap_candidates(blocks, boot_blocks,
                                           ambiguity_level)
  }
  factor_names <- lapply(observed_scores, colnames)
  secondary_null <- numeric(n_perm)
  if (length(spec$secondary)) {
    secondary_null <- unlist(chorale_deterministic_lapply(seq_len(n_perm), function(b) {
      set.seed(secondary_seeds[b])
      max(vapply(spec$secondary, function(cv) {
        conditional_scores <- lapply(modalities, function(m) {
          chorale_covariate_null_scores(
            observed_scores[[m]], aligned[[m]], spec, cv,
            exchangeability_blocks = exchangeability_blocks)
        })
        names(conditional_scores) <- modalities
        chorale_secondary_max(profile_of(conditional_scores), spec, cv)
      }, numeric(1)))
    }, n_cores), use.names = FALSE)
    secondary_cutoff <- unname(stats::quantile(secondary_null, 1 - alpha,
                                                type = 8))
    blocks <- chorale_gate_secondary(blocks, profiles, spec, secondary_cutoff)
  } else {
    secondary_cutoff <- Inf
  }
  sync <- chorale_synchronise_affinity(blocks$final, factor_names)
  if (sync$n_programmes == 0) return(empty)

  # Freedman--Lane residual permutation removes phenotype while retaining the
  # fitted contribution of every secondary shared covariate. The ICA
  # decomposition is deliberately fixed: this calibrates the supervised
  # matching/search step without paying for irrelevant unsupervised refits.
  null_draws <- chorale_deterministic_lapply(seq_len(n_perm), function(b) {
    set.seed(phenotype_seeds[b])
    null_scores <- lapply(modalities, function(m) {
      chorale_phenotype_null_scores(
        observed_scores[[m]], aligned[[m]], spec,
        exchangeability_blocks = exchangeability_blocks)
    })
    names(null_scores) <- modalities
    null_profiles <- profile_of(null_scores)
    bp <- chorale_hierarchical_blocks(null_profiles, spec, ambiguity_level)
    if (length(spec$secondary)) {
      bp <- chorale_gate_secondary(bp, null_profiles, spec, secondary_cutoff)
    }
    sp <- chorale_synchronise_affinity(bp$final, factor_names)
    # Assignment follows the full hierarchy, but phenotype support is scored
    # on the primary phenotype statistic, on the same scale as the observation.
    sp$similarity <- bp$primary
    pair <- vapply(names(blocks$primary), function(k) {
      blk <- bp$primary[[k]]
      if (is.null(blk)) 0 else max(abs(blk))
    }, numeric(1))
    list(joint = chorale_best_statistic(sp), pair = pair)
  }, n_cores)
  null_joint <- vapply(null_draws, `[[`, numeric(1), "joint")
  null_pair <- lapply(names(blocks$primary), function(k) {
    vapply(null_draws, function(x) unname(x$pair[k]), numeric(1))
  })
  names(null_pair) <- names(blocks$primary)
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
      chorale_programme_statistic(blocks$primary,
                                  members[idx, , drop = FALSE])
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
    # Loading purity is retained as a diagnostic and optional user-requested
    # filter. It is not evidence of unique latent-state recovery.
    pure <- vapply(seq_len(nrow(keep)), function(r) {
      cond <- fits[[keep$modality[r]]]$pure_feature_condition
      isTRUE(unname(cond[keep$factor[r]]))
    }, logical(1))
    pair_diagnostics <- list()
    if (nrow(keep) >= 2L) for (i in seq_len(nrow(keep))) {
      for (j in seq_len(nrow(keep))) {
        if (j <= i) next
        ma <- keep$modality[i]
        mb <- keep$modality[j]
        direct <- paste(ma, mb, sep = "|")
        reverse <- paste(mb, ma, sep = "|")
        key <- if (!is.null(blocks$diagnostics[[direct]])) direct else reverse
        dg <- blocks$diagnostics[[key]]
        ia <- keep$factor_index[i]
        ib <- keep$factor_index[j]
        if (key == reverse) {
          tmp <- ia; ia <- ib; ib <- tmp
        }
        candidates_at <- which(dg$candidate[ia, ])
        selected_secondary <- dg$secondary %||% blocks$secondary[[key]][ia, ib]
        sec_row <- blocks$secondary[[key]][ia, candidates_at, drop = TRUE]
        alternatives <- sec_row[seq_along(sec_row) != match(ib, candidates_at)]
        secondary_margin <- if (any(is.finite(alternatives))) {
          blocks$secondary[[key]][ia, ib] - max(alternatives, na.rm = TRUE)
        } else if (length(candidates_at) == 1L) {
          Inf
        } else {
          NA_real_
        }
        secondary_contributions <- vapply(spec$secondary, function(cv) {
          ta <- which(profiles[[ma]]$term_covariate == cv)
          tb <- which(profiles[[mb]]$term_covariate == cv)
          if (!length(ta) || length(ta) != length(tb)) return(NA_real_)
          block <- chorale_secondary_block_affinity(
            profiles[[ma]], profiles[[mb]], cv, dg$orientation[ia, ib])
          block[keep$factor_index[i], keep$factor_index[j]]
        }, numeric(1))
        pair_diagnostics[[length(pair_diagnostics) + 1L]] <- list(
          key = key, ia = ia, ib = ib, candidates = candidates_at,
          candidate_names = colnames(dg$candidate)[candidates_at],
          compatible = isTRUE(dg$candidate[ia, ib]),
          secondary_margin = secondary_margin,
          orientation = dg$orientation[ia, ib],
          phenotype_signal = dg$signal[ia, ib],
          phenotype_loss = dg$loss[ia, ib],
          phenotype_effect_a = paste(signif(profiles[[ma]]$effects[
            keep$factor_index[i], phenotype_terms], 5),
                                     collapse = ";"),
          phenotype_effect_b = paste(signif(profiles[[mb]]$effects[
            keep$factor_index[j], phenotype_terms] *
                                               dg$orientation[ia, ib], 5),
                                     collapse = ";"),
          phenotype_contrast = paste(phenotype_terms, collapse = ";"),
          secondary_contributions = if (length(secondary_contributions)) {
            paste0(names(secondary_contributions), "=",
                   signif(secondary_contributions, 4), collapse = ";")
          } else "")
        pair_diagnostics[[length(pair_diagnostics)]]$phenotype_margin_lower <-
          if (!is.null(dg$phenotype_margin_lower))
            dg$phenotype_margin_lower[ia, ib] else NA_real_
        pair_diagnostics[[length(pair_diagnostics)]]$phenotype_margin_upper <-
          if (!is.null(dg$phenotype_margin_upper))
            dg$phenotype_margin_upper[ia, ib] else NA_real_
      }
    }
    compatible <- all(vapply(pair_diagnostics, `[[`, logical(1), "compatible"))
    separated <- all(vapply(pair_diagnostics, function(x) {
      length(x$candidates) == 1L ||
        (is.finite(x$secondary_margin) && x$secondary_margin > 0)
    }, logical(1)))
    status <- if (!compatible) "incompatible" else if (!supported) {
      "phenotype_unsupported"
    } else if (!separated) "ambiguous" else "resolved"
    secondary_evidence <- if (length(spec$secondary) == 0L) {
      "uninformative"
    } else if (any(vapply(pair_diagnostics, function(x) {
      is.finite(x$secondary_margin) && x$secondary_margin < 0
    }, logical(1)))) {
      "conflicted"
    } else if (any(vapply(pair_diagnostics, function(x) length(x$candidates) > 1L,
                          logical(1)))) {
      "selected among phenotype-compatible candidates"
    } else {
      "confirmed"
    }

    prog_rows[[length(prog_rows) + 1L]] <- data.frame(
      programme = label,
      n_modalities = nrow(keep),
      modalities = paste(sort(keep$modality), collapse = ", "),
      modality = keep$modality,
      factor = keep$factor,
      joint_statistic = round(statistic, 4),
      joint_p = p_joint,
      supported = supported,
      resolution_status = status,
      secondary_evidence = secondary_evidence,
      phenotype_column = spec$phenotype,
      phenotype_reference = spec$levels[[spec$phenotype]][1],
      pure_features = pure,
      all_pure = all(pure),
      stringsAsFactors = FALSE
    )

    # Leaving one modality out shows whether the evidence rests on all of them
    # or on one; a programme that is no weaker without a modality was not
    # carrying it.
    if (nrow(keep) >= 3) {
      for (i in seq_len(nrow(keep))) {
        sub <- keep[-i, , drop = FALSE]
        s <- chorale_programme_statistic(blocks$primary, sub)
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

    pair_at <- 0L
    for (i in seq_len(nrow(keep))) {
      for (j in seq_len(nrow(keep))) {
        if (j <= i) next
        pair_at <- pair_at + 1L
        pd <- pair_diagnostics[[pair_at]]
        ma <- keep$modality[i]
        mb <- keep$modality[j]
        key <- pd$key
        value <- blocks$primary[[key]][pd$ia, pd$ib]
        match_rows[[length(match_rows) + 1L]] <- data.frame(
          programme = label,
          modality_a = ma, modality_b = mb,
          factor_a = keep$factor[i], factor_b = keep$factor[j],
          sign = pd$orientation,
          basis = "phenotype-led adjusted design effects",
          n_shared_covariates = length(spec$covariates),
          shared_covariates = paste(spec$covariates, collapse = ","),
          phenotype_column = spec$phenotype,
          phenotype_reference = spec$levels[[spec$phenotype]][1],
          phenotype_statistic = value,
          phenotype_loss = pd$phenotype_loss,
          phenotype_signal = pd$phenotype_signal,
          adjusted_phenotype_effect_a = pd$phenotype_effect_a,
          adjusted_phenotype_effect_b = pd$phenotype_effect_b,
          phenotype_contrast = pd$phenotype_contrast,
          secondary_contributions = pd$secondary_contributions,
          secondary_p = if (nzchar(pd$secondary_contributions)) {
            observed_secondary <- suppressWarnings(max(abs(as.numeric(sub(
              ".*=", "", strsplit(pd$secondary_contributions, ";", fixed = TRUE)[[1]]
            ))), na.rm = TRUE))
            if (is.finite(observed_secondary)) p_of(observed_secondary, secondary_null)
            else NA_real_
          } else NA_real_,
          phenotype_margin_lower = pd$phenotype_margin_lower,
          phenotype_margin_upper = pd$phenotype_margin_upper,
          candidate_set = paste(pd$candidate_names, collapse = ","),
          n_candidates = length(pd$candidates),
          secondary_margin = pd$secondary_margin,
          resolution_status = status,
          secondary_evidence = secondary_evidence,
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
    matches$significant <- matches$p_value < alpha &
      matches$resolution_status == "resolved"
    matches$p_attainable_floor <- 1 / (1 + n_perm)
    matches <- matches[order(matches$joint_p, matches$p_value,
                             -matches$statistic), , drop = FALSE]
  }
  loo <- do.call(rbind, loo_rows)
  if (is.null(loo)) loo <- data.frame()

  list(programmes = programmes, matches = matches, leave_one_out = loo,
       synchronisation = sync, terms = spec$terms, null = null_joint,
       strata_keys = spec$covariates, signature = spec,
       excluded_covariates = spec$excluded, profiles = profiles,
       affinity = blocks)
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

#' Every covariate a set of designs could anchor on
#'
#' The identifiers of samples and modalities describe the record rather than the
#' subject, and a per-sample identifier cannot contrast anything, so they are
#' excluded. Everything else is a candidate, which is what lets a design the
#' package has never seen anchor a comparison without naming its columns in
#' advance.
#'
#' @keywords internal
#' @noRd
chorale_candidate_covariates <- function(designs) {
  bookkeeping <- c("sample_id", "modality", "sample", "id", "run", "file",
                   "filename", "replicate")
  all_cols <- unique(unlist(lapply(designs, colnames)))
  candidates <- setdiff(all_cols, bookkeeping)
  # A column holding a distinct value for nearly every sample identifies the
  # sample rather than grouping it, so it cannot serve as a contrast.
  keep <- vapply(candidates, function(cv) {
    any(vapply(designs, function(d) {
      if (!cv %in% colnames(d)) return(FALSE)
      v <- stats::na.omit(d[[cv]])
      length(v) > 0 && (is.numeric(v) ||
        length(unique(v)) < max(2L, floor(0.9 * length(v))))
    }, logical(1)))
  }, logical(1))
  candidates[keep]
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
#' @param require_pure_features Apply the optional loading-purity filter. Off by
#'   default because simulations have not established a discovery threshold for
#'   this diagnostic.
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
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' chorale_programmes(fit)
chorale_programmes <- function(fit, significant_only = TRUE,
                              require_pure_features = FALSE) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  pg <- fit$programmes
  if (is.null(pg) || nrow(pg) == 0) return(data.frame())
  if (significant_only) pg <- pg[pg$supported, , drop = FALSE]
  if (require_pure_features && "all_pure" %in% colnames(pg)) {
    # This is an explicit diagnostic filter, not a default identification gate.
    pg <- pg[pg$all_pure, , drop = FALSE]
  }
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
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
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
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' chorale_joint_evidence(fit)
chorale_joint_evidence <- function(fit, programmes = NULL, n_perm = NULL,
                                   seed = NULL) {
  if (is.null(programmes)) programmes <- chorale_programmes(fit)
  programmes
}

#' Add the age band used for anchoring
#' @keywords internal
#' @noRd
chorale_add_age_bin <- function(design, n_bins = 3L) {
  if ("age_bin" %in% colnames(design)) return(design)
  age_col <- intersect(c("age", "age_months", "age_years", "age_days"),
                       colnames(design))
  if (length(age_col) == 0) return(design)
  age <- suppressWarnings(as.numeric(as.character(design[[age_col[1]]])))
  if (all(is.na(age))) return(design)

  distinct <- unique(stats::na.omit(age))
  if (length(distinct) <= n_bins) {
    # Few enough distinct ages that each is its own band; imposing cut points
    # on them would merge groups the design deliberately separates.
    design$age_bin <- ifelse(is.na(age), NA_character_, as.character(age))
    return(design)
  }
  # Bands are quantiles of the ages observed, so they follow the cohort rather
  # than an assumption about the organism's lifespan. Labels carry the range
  # each band covers, so a band means something without knowing the units.
  breaks <- unique(stats::quantile(age, probs = seq(0, 1, length.out = n_bins + 1L),
                                   na.rm = TRUE))
  if (length(breaks) < 3) return(design)
  design$age_bin <- as.character(cut(age, breaks = breaks,
                                     include.lowest = TRUE, dig.lab = 4))
  design
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Match phenotype-responsive factors across modalities
#'
#' Components are recovered per modality by independent component analysis over
#' several initialisations. Loading-purity diagnostics are then calculated,
#' with ties broken towards curated-set coherence. The remaining loadings are
#' expressed in the curated vocabulary and stored beside the fit as a
#' description of it, with the variance each explains reported.
#'
#' Factors are then matched across modalities on their design-effect profiles,
#' and the assignment is solved once over the whole collection rather than a
#' pair at a time, so the correspondences agree around every cycle. Each
#' resulting programme is scored as one object against a null that reruns the
#' complete procedure on permuted designs. Where curated sets are supplied the
#' programmes are corroborated a second time on biology, against an
#' score-residual null, and the two channels are reported separately.
#'
#' @param containers A named list of [SummarizedExperiment::SummarizedExperiment]
#'   objects, one per modality, as returned by [chorale_load()].
#' @param n_factors Integer, or one integer per modality, giving the number of
#'   components to recover, from the detectability gate. `"auto"` sets it per
#'   modality by parallel analysis, through [chorale_n_factors()].
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
#' @param transform Per-modality measurement-model transform, one of `"auto"`,
#'   `"none"`, `"log"` or `"vst"`, given as a single value or a named vector.
#'   See [chorale_transform()]. Defaults to `"auto"`, which reads the transform
#'   from each matrix.
#' @param strata_keys Design columns defining an anchoring stratum. `NULL`, the
#'   default, uses every covariate the designs share. Deprecated: use
#'   `profile_covariates` for matching and `bound_strata` for bounds.
#' @param profile_covariates Optional covariates allowed to refine the mandatory
#'   phenotype-led match. `NULL` discovers every eligible shared covariate.
#' @param bound_strata Covariates defining strata for coupling bounds. `NULL`
#'   uses the eligible shared covariates.
#' @param assay_name Assay to take from each container.
#' @param control Every decision the run takes, as [chorale_control()]: the
#'   threshold a programme must beat, the permutation counts, what counts as a
#'   pure feature, and the rest. It travels with the fit, so a result records
#'   what decided it.
#' @param seed Integer seed.
#' @param ... Named settings overriding `control`, so changing one decision does
#'   not mean restating the others.
#'
#' @returns An object of class `chorale_fit`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 2, n_per_cell = 2, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' fit
chorale_fit <- function(containers,
                        n_factors,
                        gene_sets = NULL,
                        feature_map = NULL,
                        feature_space = NULL,
                        transform = "auto",
                        strata_keys = NULL,
                        profile_covariates = NULL,
                        bound_strata = NULL,
                        assay_name = NULL,
                        control = chorale_control(),
                        seed = 1L,
                        ...) {
  if (!is.list(containers) || length(containers) < 2) {
    rlang::abort("`containers` must be a list of at least two modalities.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  modalities <- names(containers)

  # Anything named directly overrides the control object, so a caller changing
  # one decision need not rebuild the whole set, and the set that actually
  # applied is what travels with the fit.
  control <- chorale_merge_control(control, list(...))
  if (!is.null(strata_keys)) {
    lifecycle::deprecate_warn(
      "0.1.0", "chorale_fit(strata_keys)",
      details = paste0("Use `profile_covariates` for matching and ",
                       "`bound_strata` for coupling bounds."))
    profile_covariates <- profile_covariates %||% strata_keys
    bound_strata <- bound_strata %||% strata_keys
  }
  profile_covariates <- profile_covariates %||% control$profile_covariates
  bound_strata <- bound_strata %||% control$bound_strata

  # "auto" defers the count to parallel analysis on each modality, which is the
  # only parameter the estimator cannot infer from its own objective.
  auto_factors <- identical(n_factors, "auto")
  if (!auto_factors) {
    if (length(n_factors) == 1) n_factors <- rep(n_factors, length(containers))
    names(n_factors) <- modalities
  } else {
    n_factors <- stats::setNames(rep(NA_integer_, length(modalities)), modalities)
  }

  feature_space <- chorale_feature_space(feature_space, modalities)
  transform_of <- chorale_transform_spec(transform, modalities)

  fits <- list()
  designs <- list()

  for (m in modalities) {
    se <- containers[[m]]
    an <- assay_name %||% SummarizedExperiment::assayNames(se)[1]
    mat <- SummarizedExperiment::assay(se, an)
    design <- as.data.frame(SummarizedExperiment::colData(se))
    design <- chorale_add_age_bin(design)

    tf <- chorale_transform(mat, transform = transform_of[[m]])
    x <- scale(t(tf$matrix))
    x[!is.finite(x)] <- 0

    if (auto_factors) {
      n_factors[[m]] <- chorale_n_factors(
        x, quantile = control$n_factors_quantile,
        max_factors = control$max_factors, seed = seed)
    }
    fit <- chorale_ica(x, n_factors[[m]], n_init = control$n_init, seed = seed,
                       consensus = control$consensus)

    prior <- NULL
    if (!is.null(gene_sets)) {
      if (identical(unname(feature_space[[m]]), "lipid")) {
        # A lipidome reaches the sets through its classes rather than through
        # gene identifiers, so both modalities end up in one vocabulary.
        prior <- chorale_metabolite_matrix(
          rownames(mat), gene_sets,
          min_compounds = control$min_lipid_compounds,
          min_specificity = control$min_lipid_specificity,
          min_features = control$min_set_features)
      } else {
        ids <- rownames(mat)
        mapping <- NULL
        if (!is.null(feature_map) && !is.null(feature_map[[m]])) {
          mapping <- feature_map[[m]]
        }
        prior <- chorale_geneset_matrix(ids, gene_sets, mapping = mapping,
                                        min_features = control$min_set_features)
        if (ncol(prior) > 0) rownames(prior) <- rownames(mat)
      }
    }

    mk <- chorale_markers(fit$loadings, purity_ratio = control$purity_ratio,
                          min_markers = control$min_markers,
                          max_markers = control$max_markers, prior = prior)
    if (!is.null(prior) && ncol(prior) > 0) {
      # The projection onto curated sets is an annotation, so it is stored
      # beside the fit rather than substituted for it. Overwriting the
      # loadings would leave a matrix the scores no longer reconstruct.
      projected <- chorale_constrain(fit$loadings, prior, mk$markers,
                                     lambda = control$lambda)
      fit$set_weights <- projected$set_weights
      fit$pathway_loadings <- projected$loadings
      fit$reconstruction <- chorale_reconstruction(x, fit$scores, fit$loadings,
                                                   projected$loadings)
    }
    fit$transform <- tf$applied
    # Retained on the standardised analysis scale so correlation-preserving
    # pathway nulls can recompute loadings without refitting ICA.
    fit$analysis_matrix <- x
    fit$markers <- mk$markers
    fit$best_candidates <- mk$best_candidates
    fit$purity_margin <- mk$purity_margin
    fit$pure_feature_condition <- mk$pure_feature_condition
    fit$prior <- prior

    fits[[m]] <- fit
    designs[[m]] <- design
  }

  integration <- chorale_integrate(
                                   fits, designs,
                                   strata_keys = NULL,
                                   profile_covariates = profile_covariates,
                                   n_perm = control$n_perm,
                                   alpha = control$phenotype_alpha, seed = seed,
                                   phenotype_column = control$phenotype_column,
                                   phenotype_reference = control$phenotype_reference,
                                   exchangeability_blocks =
                                     control$exchangeability_blocks,
                                   ambiguity_level = control$ambiguity_level,
                                   n_ambiguity_boot = control$n_ambiguity_boot,
                                   n_cores = control$n_cores)
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
      design_profiles = integration$profiles,
      signature = integration$signature,
      excluded_covariates = integration$excluded_covariates,
      joint_null = integration$null,
      n_shared = if (nrow(matches) > 0) sum(matches$significant) else 0L,
      n_programmes = if (nrow(integration$programmes) > 0) {
        sum(!duplicated(integration$programmes$programme) &
              integration$programmes$supported)
      } else {
        0L
      },
      n_factors = n_factors,
      # The covariates actually resolved, so everything downstream conditions on
      # what anchored the comparison rather than on the request that produced it.
      strata_keys = integration$strata_keys,
      profile_covariates = integration$strata_keys,
      bound_strata = bound_strata %||% integration$strata_keys,
      gene_sets = gene_sets,
      control = control,
      phenotype_column = control$phenotype_column,
      phenotype_reference = control$phenotype_reference,
      seed = seed
    ),
    class = "chorale_fit"
  )

  # The pathway channel is separate but shares fitted factors with the design
  # channel. It is computed once here and travels with the fit.
  out$pathway_evidence <- if (control$n_pathway_perm > 0) {
    chorale_pathway_evidence(out, n_perm = control$n_pathway_perm,
                             alpha = control$alpha, seed = seed)
  } else {
    data.frame()
  }
  out
}

#' Resolve the transform requested for every modality
#' @keywords internal
#' @noRd
chorale_transform_spec <- function(transform, modalities) {
  valid <- c("auto", "none", "log", "vst")
  out <- stats::setNames(rep("auto", length(modalities)), modalities)
  if (is.null(transform)) return(out)
  if (is.null(names(transform))) {
    if (length(transform) == 1) {
      out[] <- as.character(transform)
    } else if (length(transform) == length(modalities)) {
      out[] <- as.character(transform)
    } else {
      rlang::abort("`transform` must be one value, one per modality, or named.")
    }
  } else {
    unknown <- setdiff(names(transform), modalities)
    if (length(unknown) > 0) {
      rlang::abort(paste0("`transform` names unknown modalities: ",
                          paste(unknown, collapse = ", "), "."))
    }
    out[names(transform)] <- as.character(transform)
  }
  bad <- setdiff(unique(out), valid)
  if (length(bad) > 0) {
    rlang::abort(paste0("`transform` must be one of ",
                        paste(valid, collapse = ", "), "; got: ",
                        paste(bad, collapse = ", "), "."))
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

#' Put a modality on the scale its measurement model implies
#'
#' Generic centring and scaling treats every assay as though its errors were
#' additive and its variance constant, which none of these are. Sequencing
#' counts have variance growing with the mean, so a highly expressed gene
#' dominates a factor for no biological reason. Label-free proteomic and
#' lipidomic intensities span orders of magnitude and are conventionally read
#' on the log scale. Applying the right transform first is what makes the
#' subsequent centring meaningful.
#'
#' `"auto"` reads the transform from the matrix rather than from a promise
#' about it: values that are non-negative, whole, and spread over several orders
#' of magnitude are counts and take the variance-stabilising transform; values
#' that are non-negative and heavily right-skewed are intensities and take the
#' log; anything already on a symmetric scale is left alone. The choice is
#' returned so it appears in the record rather than happening silently.
#'
#' Batch and study are not removed here. When shared and estimable, they enter
#' as adjusted secondary covariates; removing them first would erase information
#' the fitted multivariable model is meant to separate from phenotype.
#'
#' @param mat A features-by-samples numeric matrix.
#' @param transform One of `"auto"`, `"none"`, `"log"` or `"vst"`.
#'
#' @returns A list with `matrix`, the transformed features-by-samples matrix,
#'   and `applied`, the transform used.
#' @export
#' @examples
#' counts <- matrix(rpois(200, lambda = 50), nrow = 20)
#' chorale_transform(counts)$applied
chorale_transform <- function(mat, transform = c("auto", "none", "log", "vst")) {
  transform <- match.arg(transform)
  m <- as.matrix(mat)
  storage.mode(m) <- "double"
  finite <- m[is.finite(m)]

  if (transform == "auto") {
    if (length(finite) == 0) {
      transform <- "none"
    } else if (min(finite) < 0) {
      # Negative values are already on a symmetric scale, typically a log ratio.
      transform <- "none"
    } else {
      whole <- all(abs(finite - round(finite)) < 1e-8)
      spread <- stats::quantile(finite, 0.99) / max(stats::quantile(finite, 0.5), 1)
      skew <- mean(((finite - mean(finite)) / max(stats::sd(finite), 1e-9))^3)
      transform <- if (whole && spread > 5) {
        "vst"
      } else if (skew > 1 && spread > 5) {
        "log"
      } else {
        "none"
      }
    }
  }

  out <- switch(
    transform,
    none = m,
    # log1p keeps zeros finite, which a plain log would not.
    log = log1p(pmax(m, 0)),
    # Anscombe's variance-stabilising transform for counts, under which a
    # Poisson variance no longer grows with the mean.
    vst = 2 * sqrt(pmax(m, 0) + 3 / 8)
  )
  list(matrix = out, applied = transform)
}

#' How many components a modality can support
#'
#' The number of factors is the one parameter the estimator cannot infer from
#' its own objective, and setting it too high manufactures axes that are noise.
#' Parallel analysis answers it from the data: each feature is permuted
#' independently, which destroys the covariance while leaving every marginal
#' intact, and a component counts only where its eigenvalue exceeds what that
#' permuted null reaches. The threshold is therefore calibrated against the
#' modality's own marginals rather than against an assumed noise level.
#'
#' The count is a property of the cohort rather than of the feature space: a
#' modality with few samples supports few components whatever the number of
#' features measured.
#'
#' @param x A samples-by-features numeric matrix, centred and scaled.
#' @param n_perm Permutations forming the null.
#' @param quantile Quantile of the null eigenvalues a component must exceed.
#' @param max_factors Upper bound on the count returned.
#' @param seed Integer seed.
#'
#' @returns An integer count, at least two.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, seed = 1)
#' x <- scale(t(sim$modalities[[1]]))
#' chorale_n_factors(x, n_perm = 20)
chorale_n_factors <- function(x, n_perm = 100L, quantile = 0.95,
                              max_factors = 20L, seed = 1L) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  n <- nrow(x)
  observed <- svd(x, nu = 0, nv = 0)$d^2 / max(n - 1L, 1L)

  set.seed(seed)
  null <- matrix(NA_real_, nrow = n_perm, ncol = length(observed))
  for (b in seq_len(n_perm)) {
    xb <- apply(x, 2, sample)
    d <- svd(xb, nu = 0, nv = 0)$d^2 / max(n - 1L, 1L)
    null[b, seq_along(d)] <- d
  }
  threshold <- apply(null, 2, stats::quantile, probs = quantile, na.rm = TRUE)
  k <- sum(observed > threshold[seq_along(observed)], na.rm = TRUE)

  # A modality cannot support more components than a fifth of its samples
  # without the recovered axes describing individual animals.
  as.integer(max(2L, min(k, max_factors, floor(n / 5))))
}
