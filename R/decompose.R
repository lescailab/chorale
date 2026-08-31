# Primitives the estimator is built from: the decomposition, the measurement
# scale each modality is read on, and the small numerical helpers both of those
# need. Nothing here knows what the components are later used for.

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
      # Three properties decide the transform, and each is deliberately coarse:
      # the question is which measurement family a matrix belongs to, not a
      # precise characterisation of it. `whole` tolerates the rounding a count
      # matrix picks up from being written and read as a double. `spread` is the
      # 99th percentile over the median or over 1, whichever is larger, so a
      # matrix whose median is zero does not divide by it; above 5 the matrix
      # covers more than one order of magnitude, which centring alone would not
      # handle. `skew` is the third standardised moment, whose sign and rough
      # size separate a right-skewed intensity distribution from a symmetric
      # one; the 1e-9 floor only guards a constant matrix. `spread` gates both
      # transforms, so a narrow count matrix is left alone rather than
      # transformed for the sake of being whole.
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

#' The matrix a modality is read on
#'
#' The transform its measurement model implies, then centring and scaling by
#' feature, then the substitution of a standardised feature's mean for anything
#' that is not finite. Order matters: standardising first and substituting
#' afterwards leaves the mean and variance of a feature decided by the values
#' that were measured, where substituting first would let the filler move both.
#'
#' Every part of the package that reads an assay reads it through this function,
#' so a diagnostic and the estimator it describes cannot drift onto different
#' matrices.
#'
#' @param mat A features-by-samples numeric matrix.
#' @param transform One of `"auto"`, `"none"`, `"log"` or `"vst"`.
#'
#' @returns A list with `matrix`, the samples-by-features analysis matrix, and
#'   `applied`, the transform used.
#' @keywords internal
#' @noRd
chorale_analysis_matrix <- function(mat, transform = "auto") {
  tf <- chorale_transform(mat, transform = transform)
  x <- scale(t(tf$matrix))
  x[!is.finite(x)] <- 0
  list(matrix = x, applied = tf$applied)
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
#' A matrix whose eigenvalues never reach the permuted null returns zero. That
#' is the answer parallel analysis gave, and reporting it is what keeps a
#' diagnostic from being read on components the data did not support. A caller
#' that must run a decomposition whatever the count raises the floor itself,
#' where the choice is visible.
#'
#' @param x A samples-by-features numeric matrix, centred and scaled.
#' @param n_perm Permutations forming the null.
#' @param quantile Quantile of the null eigenvalues a component must exceed.
#' @param max_factors Optional further upper bound on the count returned. The
#'   count is capped at one component per five samples whether or not this is
#'   supplied, so a modality's own size sets the ceiling and no default number
#'   stands between the data and the answer.
#' @param seed Integer seed.
#'
#' @returns An integer count, zero where no component clears the null.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, seed = 1)
#' x <- scale(t(sim$modalities[[1]]))
#' chorale_n_factors(x, n_perm = 20)
chorale_n_factors <- function(x, n_perm = 100L, quantile = 0.95,
                              max_factors = NULL, seed = 1L) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  n <- nrow(x)
  observed <- svd(x, nu = 0, nv = 0)$d^2 / max(n - 1L, 1L)

  set.seed(seed)
  null <- matrix(NA_real_, nrow = n_perm, ncol = length(observed))
  for (b in seq_len(n_perm)) {
    # Each feature is permuted independently, which is what makes this the
    # parallel-analysis null rather than a resampling of whole samples:
    # permuting rows would leave the feature-feature covariance intact and there
    # would be nothing for an eigenvalue to stand above.
    xb <- apply(x, 2, sample)
    d <- svd(xb, nu = 0, nv = 0)$d^2 / max(n - 1L, 1L)
    null[b, seq_along(d)] <- d
  }
  threshold <- apply(null, 2, stats::quantile, probs = quantile, na.rm = TRUE)
  k <- sum(observed > threshold[seq_along(observed)], na.rm = TRUE)

  # A modality cannot support more components than a fifth of its samples
  # without the recovered axes describing individual samples.
  bound <- c(k, floor(n / 5))
  if (!is.null(max_factors) && !is.na(max_factors)) {
    bound <- c(bound, as.integer(max_factors))
  }
  as.integer(max(0L, min(bound)))
}

#' The same matrix, given to the factorisation in the smaller basis
#'
#' `fastICA` transposes its input and forms a feature-by-feature covariance
#' before it whitens. With many more features than samples that matrix is mostly
#' structural zeros: a modality measured on `n` samples carries at most `n`
#' independent directions however many features it holds, so a covariance over
#' tens of thousands of features built from a few hundred samples spends its
#' cost decomposing an object of rank at most `n`.
#'
#' Passing the samples-by-components representation instead leaves the whitened
#' matrix the factorisation actually iterates on unchanged. Writing the input as
#' `U D V'`, whitening on the feature side yields `sqrt(n) U'` truncated to the
#' requested components; whitening `U D` yields the same `sqrt(n) U'`, because
#' its covariance is already diagonal. The two paths therefore reach the same
#' sources from the same initialisation, and the loadings are regressed back
#' onto the original features either way.
#'
#' Column centring is preserved rather than assumed: `fastICA` centres its
#' input, and a column-centred matrix has column-centred scores, so the
#' reduction is applied to the centred matrix.
#'
#' @param x A samples-by-features numeric matrix.
#'
#' @returns `x` where it has no more features than samples, and otherwise a
#'   samples-by-`min(n, p)` matrix spanning the same row space.
#' @keywords internal
#' @noRd
chorale_ica_basis <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  if (p <= n) return(x)
  centred <- scale(x, center = TRUE, scale = FALSE)
  centred[!is.finite(centred)] <- 0
  s <- svd(centred, nu = min(n, p), nv = 0)
  z <- s$u %*% diag(s$d, nrow = length(s$d))
  rownames(z) <- rownames(x)
  z
}

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
  # The factorisation is run in the smaller basis; the loadings below are
  # regressed onto the original features, so nothing downstream sees it.
  basis <- chorale_ica_basis(x)
  best <- NULL
  best_obj <- -Inf
  obj <- rep(NA_real_, n_init)
  runs <- vector("list", n_init)

  for (i in seq_len(n_init)) {
    set.seed(seed + i)
    # `method = "C"` rather than fastICA's default `"R"`, because the medoid
    # selection below needs twenty initialisations and the R implementation is
    # too slow to run that many. The tolerance is tightened from fastICA's 1e-4
    # and the cap raised from its 200 to match: a looser tolerance leaves runs
    # stopping at different points, which the agreement matrix below would
    # report as unstable components when the instability is the optimiser's.
    fit <- try(
      fastICA::fastICA(basis, n.comp = n_factors, method = "C",
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
      # Absolute correlation, then a one-to-one assignment: ICA identifies
      # neither the sign nor the order of its components, so two runs recovering
      # the same subspace can look unrelated until their columns are paired.
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
  # This is what returns the fit to the original features after the
  # factorisation ran in the reduced basis, and it is also why the basis change
  # is invisible downstream. Only the slopes are kept: the intercept is a
  # feature's mean, which the centring has already removed from `x`.
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

#' Moore-Penrose inverse
#' @keywords internal
#' @noRd
chorale_ginv <- function(m, tol = sqrt(.Machine$double.eps)) {
  s <- svd(m)
  # Singular values are cut relative to the largest, as MASS::ginv does, so the
  # cut does not depend on the scale of the input. It is a rank decision, and it
  # is what makes the inverse defined for the rank-deficient Gram matrix that
  # overlapping concepts produce.
  keep <- s$d > max(tol * s$d[1], 0)
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
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
  # sample rather than grouping it, so it cannot serve as a contrast. The
  # threshold is nine tenths rather than all of them because a table can carry a
  # handful of repeated identifiers by accident. A numeric column is exempt: it
  # contributes a single slope however many values it realises, so having one
  # per sample costs it nothing.
  #
  # A covariate is kept where any one modality can group on it. Whether every
  # modality can is decided later, in chorale_resolve_signature(), which reports
  # the reason a covariate was excluded; deciding it here would drop the
  # covariate silently.
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
