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
chorale_bound <- function(fit, n_grid = 200L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  if (nrow(fit$matches) == 0) {
    return(structure(list(bounds = data.frame()), class = "chorale_bound"))
  }

  rows <- list()
  for (i in seq_len(nrow(fit$matches))) {
    r <- fit$matches[i, ]
    va <- fit$fits[[r$modality_a]]$scores[, r$factor_a]
    vb <- fit$fits[[r$modality_b]]$scores[, r$factor_b] * r$sign

    unconditional <- chorale_frechet_correlation(va, vb, n_grid)

    da <- fit$designs[[r$modality_a]]
    db <- fit$designs[[r$modality_b]]
    ka <- chorale_stratum_key(names(va), da, fit$strata_keys)
    kb <- chorale_stratum_key(names(vb), db, fit$strata_keys)
    common <- intersect(unique(stats::na.omit(ka)), unique(stats::na.omit(kb)))

    if (length(common) >= 2) {
      # Conditioning on the design splits the covariance into a part that is
      # identified and a part that is not. Stratum means are observed in both
      # modalities, so the between-stratum covariance is fixed; only the
      # within-stratum coupling is free, and only that part is bounded.
      anchored <- chorale_anchored_correlation(va, vb, ka, kb, common, n_grid)
      lower <- anchored$lower
      upper <- anchored$upper
    } else {
      lower <- unconditional$lower
      upper <- unconditional$upper
    }

    rows[[i]] <- data.frame(
      modality_a = r$modality_a, modality_b = r$modality_b,
      factor_a = r$factor_a, factor_b = r$factor_b,
      lower_no_anchor = unconditional$lower,
      upper_no_anchor = unconditional$upper,
      width_no_anchor = unconditional$upper - unconditional$lower,
      lower_anchored = lower,
      upper_anchored = upper,
      width_anchored = upper - lower,
      narrowing = (unconditional$upper - unconditional$lower) - (upper - lower),
      n_strata_used = length(common),
      stringsAsFactors = FALSE
    )
  }

  structure(list(bounds = do.call(rbind, rows)), class = "chorale_bound")
}

#' Correlation bounds conditional on the design strata
#'
#' The covariance decomposes into a between-stratum term, fixed by the stratum
#' means both modalities report, and a within-stratum term that the disjoint
#' samples leave free. Bounding only the free part is what the anchors buy, and
#' it can never widen the interval, since the between-stratum contribution is
#' pinned rather than extremised.
#'
#' @keywords internal
#' @noRd
chorale_anchored_correlation <- function(va, vb, ka, kb, common, n_grid) {
  sd_a <- stats::sd(va)
  sd_b <- stats::sd(vb)
  if (!is.finite(sd_a) || !is.finite(sd_b) || sd_a == 0 || sd_b == 0) {
    return(list(lower = -1, upper = 1))
  }

  in_a <- ka %in% common
  in_b <- kb %in% common
  # Stratum weights are shared by construction, so the two modalities are
  # given the same design distribution before their couplings are bounded.
  wa <- table(ka[in_a])[common]
  wb <- table(kb[in_b])[common]
  w <- (as.numeric(wa) / sum(wa) + as.numeric(wb) / sum(wb)) / 2

  mu_a <- vapply(common, function(k) mean(va[which(ka == k)]), numeric(1))
  mu_b <- vapply(common, function(k) mean(vb[which(kb == k)]), numeric(1))
  grand_a <- sum(w * mu_a)
  grand_b <- sum(w * mu_b)
  between <- sum(w * (mu_a - grand_a) * (mu_b - grand_b))

  within_lo <- 0
  within_hi <- 0
  for (i in seq_along(common)) {
    k <- common[i]
    ra <- va[which(ka == k)] - mu_a[i]
    rb <- vb[which(kb == k)] - mu_b[i]
    if (length(ra) < 2 || length(rb) < 2) next
    probs <- seq(0.5 / n_grid, 1 - 0.5 / n_grid, length.out = n_grid)
    qa <- stats::quantile(ra, probs, names = FALSE, type = 7)
    qb <- stats::quantile(rb, probs, names = FALSE, type = 7)
    hi <- mean(qa * qb)
    lo <- mean(qa * rev(qb))
    within_hi <- within_hi + w[i] * max(hi, lo)
    within_lo <- within_lo + w[i] * min(hi, lo)
  }

  lower <- (between + within_lo) / (sd_a * sd_b)
  upper <- (between + within_hi) / (sd_a * sd_b)
  list(lower = max(-1, min(lower, upper)), upper = min(1, max(lower, upper)))
}

#' Frechet bounds on the correlation of two marginals
#'
#' The extremes of the Frechet class are the comonotone and countermonotone
#' couplings, obtained by pairing the two samples in the same and in opposite
#' rank order on a common quantile grid.
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
  qa <- stats::quantile(a, probs, names = FALSE, type = 7)
  qb <- stats::quantile(b, probs, names = FALSE, type = 7)
  upper <- suppressWarnings(stats::cor(qa, qb))
  lower <- suppressWarnings(stats::cor(qa, rev(qb)))
  if (!is.finite(upper)) upper <- 1
  if (!is.finite(lower)) lower <- -1
  list(lower = min(lower, upper), upper = max(lower, upper))
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
