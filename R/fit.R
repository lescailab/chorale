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

#' Shrink non-anchor loadings towards curated sets
#'
#' Expresses the loadings of features that are not markers as a combination of
#' curated gene sets, in the manner of PLIER, so each factor carries a pathway
#' definition from estimation rather than from later annotation. Marker
#' loadings are left untouched, since they carry the purity the identification
#' argument rests on and a dense prior would erode it.
#'
#' @param loadings A features-by-factors numeric matrix.
#' @param prior A feature-by-set indicator matrix from
#'   [chorale_geneset_matrix()], with rows matching `loadings`.
#' @param markers A named list of marker features per factor.
#' @param lambda Ridge penalty on the set coefficients.
#'
#' @returns A list with `loadings`, the constrained matrix, and `set_weights`,
#'   a sets-by-factors matrix giving each factor's composition in set space.
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

#' Standardised effect of each shared covariate on each factor
#'
#' The profile is what a factor does to the design: one standardised effect per
#' covariate, comparable across modalities measured on different animals. A
#' binary covariate contributes a standardised mean difference, a continuous
#' one a rank correlation. Profiles are the matching currency because they are
#' estimated from samples, so their precision improves with sample size rather
#' than with the number of design strata.
#'
#' @param scores A samples-by-factors matrix.
#' @param design The design table for those samples.
#' @param covariates Covariate columns to profile.
#'
#' @returns A factors-by-covariates numeric matrix.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
#' se <- chorale_load(sim$modalities[[1]], sim$col_data[[1]])
#' x <- scale(t(SummarizedExperiment::assay(se)))
#' dim(chorale_design_profile(x[, 1:3, drop = FALSE], sim$col_data[[1]],
#'                            c("phenotype", "sex")))
chorale_design_profile <- function(scores, design, covariates) {
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  out <- matrix(0, nrow = ncol(scores), ncol = length(covariates),
                dimnames = list(colnames(scores), covariates))
  for (cv in covariates) {
    v <- design[[cv]]
    for (j in seq_len(ncol(scores))) {
      y <- scores[, j]
      ok <- is.finite(y) & !is.na(v)
      if (sum(ok) < 4) next
      out[j, cv] <- chorale_effect(y[ok], v[ok])$effect
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

#' Match factors across modalities
#'
#' With disjoint samples no animal is shared, so factors cannot be matched by
#' correlating scores. Two things can still be compared.
#'
#' Where the modalities share at least one design covariate, a factor
#' measuring the same latent state should act on that covariate the same way in
#' both. The statistic is the inner product of the two design-effect profiles,
#' and the null is built by permuting the covariate labels across samples
#' within each modality. Permuting samples rather than strata matters: the null
#' space is then the sample permutations, so precision improves as the cohorts
#' grow, and a single shared covariate suffices. Phenotype alone is enough.
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
#' @returns A data frame, one row per assigned cross-modality factor pair.
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
  modalities <- names(fits)
  rows <- list()

  for (i in seq_along(modalities)) {
    for (j in seq_along(modalities)) {
      if (j <= i) next
      a <- modalities[i]
      b <- modalities[j]
      sa <- fits[[a]]$scores
      sb <- fits[[b]]$scores
      da <- designs[[a]][match(rownames(sa), designs[[a]]$sample_id), , drop = FALSE]
      db <- designs[[b]][match(rownames(sb), designs[[b]]$sample_id), , drop = FALSE]

      shared <- chorale_shared_covariates(da, db, strata_keys)
      chorale_require_phenotype(shared, a, b)
      basis <- "design effects"
      pa <- chorale_design_profile(sa, da, shared)
      pb <- chorale_design_profile(sb, db, shared)

      stat <- abs(pa %*% t(pb))
      stat[!is.finite(stat)] <- 0
      assignment <- chorale_assign(stat)

      set.seed(seed)
      null_best <- numeric(n_perm)
      for (b_i in seq_len(n_perm)) {
        # Permute the covariates across samples, which is the null of a factor
        # unrelated to the design. The permutation space is the samples, so it
        # does not shrink when few strata are populated.
        da_p <- da
        db_p <- db
        for (cv in shared) {
          da_p[[cv]] <- sample(da_p[[cv]])
          db_p[[cv]] <- sample(db_p[[cv]])
        }
        qa <- chorale_design_profile(sa, da_p, shared)
        qb <- chorale_design_profile(sb, db_p, shared)
        g <- abs(qa %*% t(qb))
        g[!is.finite(g)] <- 0
        null_best[b_i] <- max(g)
      }

      for (ca in seq_len(ncol(sa))) {
        cb <- assignment[ca]
        if (is.na(cb)) next
        value <- stat[ca, cb]
        p <- (1 + sum(null_best >= value)) / (1 + n_perm)
        direction <- sum(pa[ca, ] * pb[cb, ])
        rows[[length(rows) + 1]] <- data.frame(
          modality_a = a, modality_b = b,
          factor_a = colnames(sa)[ca], factor_b = colnames(sb)[cb],
          sign = if (direction < 0) -1 else 1,
          basis = basis,
          n_shared_covariates = length(shared),
          shared_covariates = paste(shared, collapse = ","),
          statistic = as.numeric(value),
          p_value = p,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out$significant <- out$p_value < alpha
  out$p_attainable_floor <- 1 / (1 + n_perm)
  out[order(out$p_value, -out$statistic), , drop = FALSE]
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
chorale_require_phenotype <- function(shared, a, b) {
  if (!"phenotype" %in% shared) {
    rlang::abort(
      paste0(
        "Modalities '", a, "' and '", b, "' do not share a phenotype contrast. ",
        "chorale estimates a case/control shared state, so every modality must ",
        "carry a `phenotype` column taking at least two values. Resolve the ",
        "phenotype for both modalities, or drop the one that lacks it."
      ),
      class = "chorale_missing_phenotype"
    )
  }
  invisible(TRUE)
}

#' Covariates present and varying in both designs
#' @keywords internal
#' @noRd
chorale_shared_covariates <- function(da, db, candidates) {
  usable <- function(d, cv) {
    if (!cv %in% colnames(d)) return(FALSE)
    v <- d[[cv]]
    length(unique(stats::na.omit(v))) >= 2
  }
  shared <- candidates[vapply(candidates, function(cv) {
    usable(da, cv) && usable(db, cv)
  }, logical(1))]
  # age_bin and age_months describe the same variable, so the finer one is
  # kept and the coarser dropped to avoid counting it twice.
  if (all(c("age_months", "age_bin") %in% shared)) {
    shared <- setdiff(shared, "age_bin")
  }
  shared
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

#' Consolidate pairwise matches into programmes spanning every modality
#'
#' Matching is evaluated a pair of modalities at a time, because a design
#' profile is compared between two modalities. A latent programme, though, is
#' not a pair: one factor may match partners in several modalities, and those
#' partners are then the same programme seen in each. Factors are therefore
#' linked into connected components, so a programme carries every modality that
#' measures it rather than appearing once per pair.
#'
#' Three modalities agreeing on one programme is stronger evidence than three
#' separate agreements, which is the sense in which false-discovery control
#' improves with the number of modalities.
#'
#' @param fit A `chorale_fit` object.
#' @param significant_only Link only matches that beat their permutation null.
#'
#' @returns A data frame with one row per (programme, modality, factor)
#'   membership, carrying `programme`, `n_modalities` and `modalities`.
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
  m <- fit$matches
  if (nrow(m) == 0) return(data.frame())
  if (significant_only && "significant" %in% colnames(m)) {
    m <- m[m$significant, , drop = FALSE]
  }
  if (nrow(m) == 0) return(data.frame())

  node_a <- paste(m$modality_a, m$factor_a, sep = "|")
  node_b <- paste(m$modality_b, m$factor_b, sep = "|")
  nodes <- unique(c(node_a, node_b))

  # Components by label propagation: a factor and every partner it matches,
  # transitively, form one programme.
  label <- stats::setNames(seq_along(nodes), nodes)
  repeat {
    changed <- FALSE
    for (r in seq_len(nrow(m))) {
      a <- node_a[r]
      b <- node_b[r]
      lo <- min(label[[a]], label[[b]])
      if (label[[a]] != lo || label[[b]] != lo) {
        label[[a]] <- lo
        label[[b]] <- lo
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  root <- unname(label[nodes])

  # Order programmes by the strongest evidence they carry.
  comp <- unique(root)
  best <- vapply(comp, function(rt) {
    members <- nodes[root == rt]
    rows <- m[node_a %in% members | node_b %in% members, , drop = FALSE]
    min(rows$p_value)
  }, numeric(1))
  ord <- comp[order(best)]

  out <- list()
  for (i in seq_along(ord)) {
    members <- nodes[root == ord[i]]
    parts <- do.call(rbind, strsplit(members, "|", fixed = TRUE))
    mods <- unique(parts[, 1])
    rows <- m[node_a %in% members | node_b %in% members, , drop = FALSE]
    out[[i]] <- data.frame(
      programme = paste0("P", i),
      n_modalities = length(mods),
      modalities = paste(sort(mods), collapse = ", "),
      modality = parts[, 1],
      factor = parts[, 2],
      best_match_p = signif(min(rows$p_value), 3),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Joint evidence that a programme is carried by all its modalities
#'
#' Consolidating pairwise matches groups the factors correctly but leaves the
#' evidence pairwise, and taking the strongest pairwise p-value as a
#' programme's would be anti-conservative and would forfeit the one statistical
#' advantage of measuring many modalities.
#'
#' The statistic here is joint: the mean agreement of the design profiles over
#' every pair within the programme, evaluated as one quantity. The null is
#' joint too. Covariate labels are permuted in every contributing modality at
#' once, and the null keeps the best value attainable by any combination of
#' factors, one per modality, over the same modality set. Requiring three
#' modalities to agree simultaneously is far harder to achieve by chance than
#' requiring two, so the p-value tightens as modalities are added, which is the
#' finite-sample behaviour the modality-count bound describes.
#'
#' @param fit A `chorale_fit` object.
#' @param programmes Output of [chorale_programmes()]; recomputed if absent.
#' @param n_perm Number of joint permutations.
#' @param seed Integer seed.
#'
#' @returns `programmes` with `joint_statistic` and `joint_p` added.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 3, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
#' chorale_joint_evidence(fit, n_perm = 20)
chorale_joint_evidence <- function(fit, programmes = NULL, n_perm = 200L,
                                   seed = 1L) {
  if (is.null(programmes)) programmes <- chorale_programmes(fit)
  if (nrow(programmes) == 0) return(programmes)

  profiles <- function(designs_by_mod, mods, covs) {
    lapply(mods, function(m) {
      chorale_design_profile(fit$fits[[m]]$scores, designs_by_mod[[m]], covs)
    })
  }
  # Mean agreement over every pair inside one candidate programme.
  joint_stat <- function(prof, pick) {
    mods <- seq_along(prof)
    if (length(mods) < 2) return(0)
    vals <- c()
    for (i in mods) {
      for (j in mods) {
        if (j <= i) next
        vals <- c(vals, abs(sum(prof[[i]][pick[i], ] * prof[[j]][pick[j], ])))
      }
    }
    mean(vals)
  }

  out <- programmes
  out$joint_statistic <- NA_real_
  out$joint_p <- NA_real_

  for (pr in unique(programmes$programme)) {
    d <- programmes[programmes$programme == pr, , drop = FALSE]
    mods <- d$modality
    if (length(mods) < 2) next
    covs <- Reduce(intersect, lapply(mods, function(m) {
      chorale_shared_covariates(fit$designs[[m]], fit$designs[[m]],
                                fit$strata_keys)
    }))
    covs <- Reduce(intersect, lapply(mods, function(m) {
      keep <- covs[vapply(covs, function(cv) {
        cv %in% colnames(fit$designs[[m]]) &&
          length(unique(stats::na.omit(fit$designs[[m]][[cv]]))) >= 2
      }, logical(1))]
      keep
    }))
    if (length(covs) == 0) next

    prof <- profiles(fit$designs, mods, covs)
    pick <- vapply(seq_along(mods), function(i) {
      match(d$factor[i], rownames(prof[[i]]))
    }, numeric(1))
    if (any(is.na(pick))) next
    observed <- joint_stat(prof, pick)

    grid <- expand.grid(lapply(prof, function(pm) seq_len(nrow(pm))))
    set.seed(seed)
    null <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      permuted <- lapply(mods, function(m) {
        dm <- fit$designs[[m]]
        for (cv in covs) dm[[cv]] <- sample(dm[[cv]])
        dm
      })
      names(permuted) <- mods
      pp <- profiles(permuted, mods, covs)
      null[b] <- max(apply(grid, 1, function(g) joint_stat(pp, as.integer(g))))
    }
    p <- (1 + sum(null >= observed)) / (1 + n_perm)
    out$joint_statistic[out$programme == pr] <- round(observed, 4)
    out$joint_p[out$programme == pr] <- p
  }
  out
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
                        n_init = 20L,
                        strata_keys = c("phenotype", "age_bin", "sex"),
                        assay_name = NULL,
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

    mk <- chorale_markers(fit$loadings, prior = prior)
    if (!is.null(prior) && ncol(prior) > 0) {
      constrained <- chorale_constrain(fit$loadings, prior, mk$markers)
      fit$loadings <- constrained$loadings
      fit$set_weights <- constrained$set_weights
    }
    fit$markers <- mk$markers
    fit$best_candidates <- mk$best_candidates
    fit$purity_margin <- mk$purity_margin
    fit$pure_feature_condition <- mk$pure_feature_condition
    fit$prior <- prior

    fits[[m]] <- fit
    designs[[m]] <- design
  }

  matches <- chorale_match(fits, designs, strata_keys = strata_keys,
                           seed = seed)

  structure(
    list(
      modalities = modalities,
      fits = fits,
      designs = designs,
      matches = matches,
      n_shared = if (nrow(matches) > 0) sum(matches$significant) else 0L,
      n_factors = n_factors,
      strata_keys = strata_keys,
      gene_sets = gene_sets,
      seed = seed
    ),
    class = "chorale_fit"
  )
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
  cat("  assigned factor pairs:",
      if (is.data.frame(x$matches)) nrow(x$matches) else 0,
      "of which significant:", x$n_shared, "\n")
  invisible(x)
}
