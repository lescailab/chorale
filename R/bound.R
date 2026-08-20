#' Bound the cross-modality coupling
#'
#' No animal is measured in more than one modality, so the joint distribution
#' of a factor pair is not identified: the data fix each marginal and say
#' nothing about how they are coupled. What is available is the Frechet class,
#' the set of joints consistent with those marginals, and its bounds are the
#' honest report. For a correlation the extremes are attained by the
#' comonotone and countermonotone couplings, giving the widest interval the
#' data cannot exclude.
#'
#' Conditioning on the design narrows the interval, because the coupling is
#' then only free within a stratum rather than across the whole sample. The
#' difference between the unconditional and conditional widths is what the
#' anchors buy, and reporting both is the point: where the bounds stay wide,
#' the width is the result.
#'
#' @param fit A `chorale_fit` object, as returned by [chorale_fit()].
#' @param n_grid Number of quantiles used to represent each marginal.
#' @param include_ambiguous Include correspondences that phenotype-led matching
#'   did not resolve. The default bounds only resolved correspondences.
#'
#' @returns An object of class `chorale_bound`, wrapping a data frame with one
#'   row per matched factor pair: the Frechet bounds on their correlation
#'   without anchors, the bounds within design strata, and the narrowing the
#'   anchors achieve.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' chorale_bound(fit)
chorale_bound <- function(fit, n_grid = 200L, include_ambiguous = FALSE) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  matches <- fit$matches
  if (!include_ambiguous && "resolution_status" %in% names(matches)) {
    matches <- matches[matches$resolution_status == "resolved", , drop = FALSE]
  }
  if (nrow(matches) == 0) {
    return(structure(list(bounds = data.frame()), class = "chorale_bound"))
  }

  rows <- list()
  strata <- fit$bound_strata %||% fit$strata_keys
  for (i in seq_len(nrow(matches))) {
    r <- matches[i, ]
    va <- fit$fits[[r$modality_a]]$scores[, r$factor_a]
    vb <- fit$fits[[r$modality_b]]$scores[, r$factor_b] * r$sign

    da <- fit$designs[[r$modality_a]]
    db <- fit$designs[[r$modality_b]]
    ka <- chorale_stratum_key(names(va), da, strata)
    kb <- chorale_stratum_key(names(vb), db, strata)
    common <- sort(intersect(unique(stats::na.omit(ka)),
                             unique(stats::na.omit(kb))))
    weights <- if (length(common) >= 2) {
      chorale_target_weights(ka, kb, common)
    } else {
      NULL
    }

    if (is.null(weights)) {
      # No design distribution is common to both, so there is nothing to
      # condition on and the two intervals coincide.
      unconditional <- chorale_frechet_correlation(va, vb, n_grid)
      anchored <- unconditional
      n_used <- 0L
    } else {
      # Both marginals are reweighted to the same target design distribution
      # before either interval is computed, so the widths are comparable and
      # the narrowing is what conditioning on the design actually buys.
      probs <- seq(0.5 / n_grid, 1 - 0.5 / n_grid, length.out = n_grid)
      unconditional <- chorale_frechet_from_quantiles(
        chorale_weighted_quantile(va, weights$a, probs),
        chorale_weighted_quantile(vb, weights$b, probs)
      )
      anchored <- chorale_anchored_correlation(va, vb, ka, kb, common, n_grid,
                                               weights)
      n_used <- length(common)
    }
    lower <- anchored$lower
    upper <- anchored$upper

    rows[[i]] <- data.frame(
      modality_a = r$modality_a, modality_b = r$modality_b,
      factor_a = r$factor_a, factor_b = r$factor_b,
      resolution_status = if ("resolution_status" %in% names(r)) {
        r$resolution_status
      } else {
        "legacy"
      },
      lower_no_anchor = unconditional$lower,
      upper_no_anchor = unconditional$upper,
      width_no_anchor = unconditional$upper - unconditional$lower,
      lower_anchored = lower,
      upper_anchored = upper,
      width_anchored = upper - lower,
      narrowing = (unconditional$upper - unconditional$lower) - (upper - lower),
      n_strata_used = n_used,
      stringsAsFactors = FALSE
    )
  }

  structure(list(bounds = do.call(rbind, rows)), class = "chorale_bound")
}

#' Quantiles of a weighted empirical distribution
#'
#' Reweighting a sample to a target design distribution is what puts two
#' modalities on the same population before their coupling is bounded. The
#' quantile function of the reweighted marginal is obtained by inverting its
#' weighted empirical distribution.
#'
#' @keywords internal
#' @noRd
chorale_weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w[o])
  # Midpoint convention, so a point mass is represented at its centre rather
  # than at either edge.
  cw <- cumsum(w) - w / 2
  stats::approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

#' Frechet bounds on the correlation of two reweighted marginals
#'
#' Both marginals are represented on the same quantile grid and both sets of
#' moments are taken from that representation, so the cross-moment in the
#' numerator and the standard deviations in the denominator describe the same
#' target population. The extremes are the comonotone and countermonotone
#' couplings.
#'
#' @keywords internal
#' @noRd
chorale_frechet_from_quantiles <- function(qa, qb) {
  qa <- qa[is.finite(qa)]
  qb <- qb[is.finite(qb)]
  if (length(qa) < 3 || length(qb) < 3 || length(qa) != length(qb)) {
    return(list(lower = -1, upper = 1))
  }
  mu_a <- mean(qa)
  mu_b <- mean(qb)
  sd_a <- sqrt(mean((qa - mu_a)^2))
  sd_b <- sqrt(mean((qb - mu_b)^2))
  if (!is.finite(sd_a) || !is.finite(sd_b) || sd_a == 0 || sd_b == 0) {
    return(list(lower = -1, upper = 1))
  }
  hi <- (mean(qa * qb) - mu_a * mu_b) / (sd_a * sd_b)
  lo <- (mean(qa * rev(qb)) - mu_a * mu_b) / (sd_a * sd_b)
  chorale_clamp_bounds(min(lo, hi), max(lo, hi))
}

#' Order and clamp an identified set for a correlation
#'
#' A correlation lies in `[-1, 1]`, so each endpoint is clamped on its own and
#' the pair is ordered. Clamping one endpoint using the other, as an earlier
#' version did, can leave an interval that no correlation satisfies.
#'
#' @keywords internal
#' @noRd
chorale_clamp_bounds <- function(lower, upper) {
  if (!is.finite(lower)) lower <- -1
  if (!is.finite(upper)) upper <- 1
  lo <- min(lower, upper)
  hi <- max(lower, upper)
  list(lower = max(-1, min(1, lo)), upper = max(-1, min(1, hi)))
}

#' Weights putting both modalities on one target design distribution
#'
#' The two modalities realise the design in different proportions, so before
#' their coupling can be bounded they have to be reweighted to one declared
#' target. The target is the average of the two stratum distributions over the
#' strata both populate, and every sample in a stratum carries the same weight.
#'
#' @keywords internal
#' @noRd
chorale_target_weights <- function(ka, kb, common) {
  na <- table(factor(ka[ka %in% common], levels = common))
  nb <- table(factor(kb[kb %in% common], levels = common))
  if (any(na == 0) || any(nb == 0)) return(NULL)
  w <- (as.numeric(na) / sum(na) + as.numeric(nb) / sum(nb)) / 2
  w <- w / sum(w)
  names(w) <- common
  list(
    stratum = w,
    a = ifelse(ka %in% common, w[match(ka, common)] / as.numeric(na)[match(ka, common)], 0),
    b = ifelse(kb %in% common, w[match(kb, common)] / as.numeric(nb)[match(kb, common)], 0)
  )
}

#' Correlation bounds conditional on the design strata
#'
#' Both marginals are first reweighted to one target design distribution, so
#' the moments in the numerator and the denominator describe the same
#' population. Under that target the covariance splits into a between-stratum
#' term, fixed by the stratum means both modalities report, and a
#' within-stratum term that the disjoint samples leave free. Only the free part
#' is extremised, which is what the anchors buy, and the total variance is
#' taken under the same target, so the ratio is a correlation of the target
#' population and cannot leave `[-1, 1]`.
#'
#' @keywords internal
#' @noRd
chorale_anchored_correlation <- function(va, vb, ka, kb, common, n_grid,
                                         weights) {
  w <- weights$stratum
  probs <- seq(0.5 / n_grid, 1 - 0.5 / n_grid, length.out = n_grid)

  mu_a <- vapply(common, function(k) mean(va[which(ka == k)]), numeric(1))
  mu_b <- vapply(common, function(k) mean(vb[which(kb == k)]), numeric(1))
  if (!all(is.finite(mu_a)) || !all(is.finite(mu_b))) {
    return(list(lower = -1, upper = 1))
  }
  grand_a <- sum(w * mu_a)
  grand_b <- sum(w * mu_b)
  between <- sum(w * (mu_a - grand_a) * (mu_b - grand_b))

  var_a <- sum(w * (mu_a - grand_a)^2)
  var_b <- sum(w * (mu_b - grand_b)^2)
  within_lo <- 0
  within_hi <- 0
  for (i in seq_along(common)) {
    k <- common[i]
    ra <- va[which(ka == k)] - mu_a[i]
    rb <- vb[which(kb == k)] - mu_b[i]
    qa <- if (length(ra) >= 2) {
      stats::quantile(ra, probs, names = FALSE, type = 7)
    } else {
      rep(0, n_grid)
    }
    qb <- if (length(rb) >= 2) {
      stats::quantile(rb, probs, names = FALSE, type = 7)
    } else {
      rep(0, n_grid)
    }
    # Moments of the within-stratum residuals are read off the same quantile
    # grid the extremal cross-moments use, so the Cauchy-Schwarz inequality
    # that keeps the ratio inside [-1, 1] holds exactly.
    qa <- qa - mean(qa)
    qb <- qb - mean(qb)
    var_a <- var_a + w[i] * mean(qa^2)
    var_b <- var_b + w[i] * mean(qb^2)
    hi <- mean(qa * qb)
    lo <- mean(qa * rev(qb))
    within_hi <- within_hi + w[i] * max(hi, lo)
    within_lo <- within_lo + w[i] * min(hi, lo)
  }

  denom <- sqrt(var_a) * sqrt(var_b)
  if (!is.finite(denom) || denom <= 0) return(list(lower = -1, upper = 1))
  chorale_clamp_bounds((between + within_lo) / denom,
                       (between + within_hi) / denom)
}

#' Frechet bounds on the correlation of two unweighted marginals
#'
#' Used where the modalities populate no design stratum in common, so there is
#' no target distribution to reweight to and each marginal is taken as it was
#' observed.
#'
#' @keywords internal
#' @noRd
chorale_frechet_correlation <- function(a, b, n_grid = 200L) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3) {
    return(list(lower = -1, upper = 1))
  }
  probs <- seq(0.5 / n_grid, 1 - 0.5 / n_grid, length.out = n_grid)
  chorale_frechet_from_quantiles(
    stats::quantile(a, probs, names = FALSE, type = 7),
    stats::quantile(b, probs, names = FALSE, type = 7)
  )
}

#' Stratum key for a set of samples
#' @keywords internal
#' @noRd
chorale_stratum_key <- function(sample_ids, design, strata_keys) {
  keys <- intersect(strata_keys, colnames(design))
  if (length(keys) == 0) return(rep(NA_character_, length(sample_ids)))
  idx <- match(sample_ids, design$sample_id)
  parts <- lapply(keys, function(k) as.character(design[[k]][idx]))
  key <- do.call(paste, c(parts, sep = "|"))
  key[!stats::complete.cases(do.call(cbind, parts))] <- NA_character_
  key
}

#' @export
print.chorale_bound <- function(x, ...) {
  cat("<chorale_bound>\n")
  if (nrow(x$bounds) == 0) {
    cat("  no matched factor pairs to bound\n")
    return(invisible(x))
  }
  cat("  matched pairs bounded:", nrow(x$bounds), "\n")
  cat("  median width without anchors:",
      round(stats::median(x$bounds$width_no_anchor), 3), "\n")
  cat("  median width with anchors:   ",
      round(stats::median(x$bounds$width_anchored), 3), "\n")
  invisible(x)
}

#' Sampling uncertainty in the identified set
#'
#' The interval [chorale_bound()] reports is a plug-in estimate: it is the set
#' of correlations consistent with the marginals as observed, and says nothing
#' about how those marginals would move under resampling. A programme built on
#' forty animals and one built on four hundred can carry the same plug-in
#' interval and very different confidence in it.
#'
#' Samples are resampled within design stratum, so the design distribution the
#' bounds are conditioned on is preserved, and the interval is recomputed on
#' each resample. The reported sensitivity envelope runs from the lower percentile of the
#' lower endpoints to the upper percentile of the upper endpoints, which is the
#' outer envelope of the identified sets the data support rather than an
#' interval for a point.
#'
#' Factor estimation is not resampled by default, so the region carries sampling
#' error in the marginals and not the uncertainty of the factorisation itself.
#' Setting `refit` reruns the whole estimator on each resample, which does carry
#' it, at a cost of one full fit per replicate.
#'
#' @param fit A `chorale_fit` object.
#' @param containers The modality containers the fit was built from. Required
#'   when `refit` is set.
#' @param n_boot Number of resamples.
#' @param level Central proportion of bootstrap endpoints retained in the
#'   sensitivity envelope. This is not asserted to be confidence coverage.
#' @param n_grid Quantiles representing each marginal.
#' @param refit Refit the estimator on every resample, so the region carries
#'   the uncertainty of the factorisation as well as of the marginals.
#' @param n_init Initialisations per refit, when `refit` is set.
#' @param seed Integer seed.
#'
#' @returns A data frame with one row per matched factor pair, carrying the
#'   plug-in endpoints and the bootstrap sensitivity envelope around them.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' chorale_bound_uncertainty(fit, n_boot = 20)
chorale_bound_uncertainty <- function(fit, containers = NULL, n_boot = 200L,
                                      level = 0.95, n_grid = 200L,
                                      refit = FALSE, n_init = 5L, seed = 1L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  if (refit && is.null(containers)) {
    rlang::abort("`containers` is required when `refit` is set.")
  }
  base <- chorale_bound(fit, n_grid = n_grid)$bounds
  if (nrow(base) == 0) return(data.frame())

  key <- paste(base$modality_a, base$factor_a, base$modality_b, base$factor_b)
  lo <- matrix(NA_real_, nrow = n_boot, ncol = nrow(base),
               dimnames = list(NULL, key))
  hi <- lo

  # Resampling within stratum keeps each modality's design distribution, so the
  # conditioning the bounds rest on is not itself resampled away.
  strata_of <- function(m) {
    s <- fit$fits[[m]]$scores
    chorale_stratum_key(rownames(s), fit$designs[[m]],
                        fit$bound_strata %||% fit$strata_keys)
  }

  set.seed(seed)
  for (b in seq_len(n_boot)) {
    boot <- fit
    resampled_containers <- if (refit) {
      stats::setNames(vector("list", length(fit$modalities)), fit$modalities)
    } else {
      NULL
    }
    ok <- TRUE
    for (m in fit$modalities) {
      s <- fit$fits[[m]]$scores
      k <- strata_of(m)
      idx <- unlist(lapply(split(seq_len(nrow(s)), k), function(i) {
        if (length(i) == 0) return(integer(0))
        sample(i, length(i), replace = TRUE)
      }), use.names = FALSE)
      if (length(idx) < 4) {
        ok <- FALSE
        break
      }
      resampled <- s[idx, , drop = FALSE]
      # Names must stay unique for the design join, and each resampled row
      # keeps the design row it was drawn from.
      d <- fit$designs[[m]]
      d <- d[match(rownames(s)[idx], d$sample_id), , drop = FALSE]
      new_id <- make.unique(as.character(d$sample_id))
      rownames(resampled) <- d$sample_id <- new_id
      boot$fits[[m]]$scores <- resampled
      boot$designs[[m]] <- d
      if (refit) {
        se <- containers[[m]]
        at <- match(rownames(s)[idx], colnames(se))
        if (anyNA(at)) {
          ok <- FALSE
          break
        }
        se <- se[, at, drop = FALSE]
        colnames(se) <- new_id
        cd <- as.data.frame(SummarizedExperiment::colData(se))
        cd$sample_id <- new_id
        SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)
        resampled_containers[[m]] <- se
      }
    }
    if (!ok) next
    if (refit) {
      rf <- try(chorale_fit(resampled_containers, n_factors = fit$n_factors,
                            n_init = n_init,
                            profile_covariates = fit$profile_covariates,
                            bound_strata = fit$bound_strata,
                            phenotype_column = fit$phenotype_column,
                            phenotype_reference = fit$phenotype_reference,
                            n_perm = 19L, n_ambiguity_boot = 0L,
                            n_pathway_perm = 0L, seed = seed + b), silent = TRUE)
      if (inherits(rf, "try-error")) next
      rf <- chorale_align_refit_factors(rf, fit)
      # Bounds quantify the uncertainty of the correspondence selected in the
      # original fit. Re-running selection inside each bootstrap would mix
      # assignment uncertainty into endpoint sensitivity and often change the
      # target. After alignment, reuse the original labelled correspondence.
      rf$matches <- fit$matches
      boot <- rf
    }
    bb <- try(chorale_bound(boot, n_grid = n_grid)$bounds, silent = TRUE)
    if (inherits(bb, "try-error") || nrow(bb) == 0) next
    bkey <- paste(bb$modality_a, bb$factor_a, bb$modality_b, bb$factor_b)
    j <- match(key, bkey)
    lo[b, ] <- bb$lower_anchored[j]
    hi[b, ] <- bb$upper_anchored[j]
  }

  a <- (1 - level) / 2
  region_lower <- apply(lo, 2, stats::quantile, probs = a, na.rm = TRUE)
  region_upper <- apply(hi, 2, stats::quantile, probs = 1 - a, na.rm = TRUE)
  n_used <- apply(lo, 2, function(v) sum(is.finite(v)))

  data.frame(
    base[, c("modality_a", "modality_b", "factor_a", "factor_b")],
    lower_anchored = base$lower_anchored,
    upper_anchored = base$upper_anchored,
    region_lower = round(pmax(-1, region_lower), 4),
    region_upper = round(pmin(1, region_upper), 4),
    region_width = round(pmin(1, region_upper) - pmax(-1, region_lower), 4),
    n_resamples = as.integer(n_used),
    sensitivity_level = level,
    stringsAsFactors = FALSE
  )
}

#' Align refitted factor names to the original loading directions
#' @keywords internal
#' @noRd
chorale_align_refit_factors <- function(refit, original) {
  for (m in intersect(refit$modalities, original$modalities)) {
    a <- original$fits[[m]]$loadings
    b <- refit$fits[[m]]$loadings
    shared <- intersect(rownames(a), rownames(b))
    if (length(shared) < 3L || ncol(a) != ncol(b)) next
    agreement <- abs(suppressWarnings(stats::cor(a[shared, , drop = FALSE],
                                                  b[shared, , drop = FALSE])))
    agreement[!is.finite(agreement)] <- 0
    assignment <- as.integer(clue::solve_LSAP(agreement, maximum = TRUE))
    old <- colnames(b)
    new <- colnames(a)[match(seq_len(ncol(b)), assignment)]
    new[is.na(new)] <- old[is.na(new)]
    colnames(refit$fits[[m]]$scores) <- new
    colnames(refit$fits[[m]]$loadings) <- new
    if (nrow(refit$matches)) {
      at_a <- refit$matches$modality_a == m
      at_b <- refit$matches$modality_b == m
      refit$matches$factor_a[at_a] <- new[match(refit$matches$factor_a[at_a], old)]
      refit$matches$factor_b[at_b] <- new[match(refit$matches$factor_b[at_b], old)]
    }
  }
  refit
}
