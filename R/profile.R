#' Summarise the properties of a real matrix a simulation must reproduce
#'
#' The estimator reads a normalised feature-by-sample matrix, so the properties
#' a simulation has to get right are the properties of that matrix rather than
#' of the measurement that produced it. Six of them are what the estimator can
#' respond to: the per-feature marginal distribution, the mean-variance
#' relation, missingness as a function of abundance, the feature-feature
#' correlation structure, non-Gaussianity relevant to ICA behaviour, and the
#' margins of the design. This function measures all six and returns them as one
#' object, which [chorale_simulate()] consumes to give simulated data the shape
#' of real data without reproducing the measurement model that produced it.
#'
#' The object holds summaries, not data. No sample identifier enters it and no
#' value of any covariate the caller withholds. Feature identifiers are kept,
#' because a planted pathway is defined on them.
#'
#' @param assay A feature-by-sample numeric matrix, on the scale the estimator
#'   would read it. Missing values are recorded rather than imputed.
#' @param col_data Optional design table for the same samples, one row per
#'   column of `assay`, from which the design margins are taken.
#' @param covariates Design columns to profile. `NULL`, the default, takes every
#'   column except `sample_id`. Name the columns explicitly to keep identifying
#'   fields out of the profile.
#' @param n_quantiles Number of points in the per-feature quantile grid. The
#'   grid is what a simulated feature is mapped onto, so it sets how finely a
#'   marginal is reproduced.
#' @param n_eigen Number of leading eigenvalues retained from the correlation
#'   structure.
#' @param layer Optional name recorded with the profile.
#'
#' @returns An object of class `chorale_data_profile`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 40, seed = 1)
#' chorale_data_profile(sim$modalities[[1]], sim$col_data[[1]],
#'                      covariates = c("phenotype", "sex"))
chorale_data_profile <- function(assay, col_data = NULL, covariates = NULL,
                                 n_quantiles = 25L, n_eigen = 50L,
                                 layer = NULL) {
  m <- as.matrix(assay)
  storage.mode(m) <- "double"
  if (nrow(m) < 2 || ncol(m) < 2) {
    rlang::abort("`assay` needs at least two features and two samples.")
  }
  if (n_quantiles < 3) {
    rlang::abort("`n_quantiles` must be at least 3.")
  }
  if (is.null(rownames(m))) {
    rownames(m) <- sprintf("feature_%05d", seq_len(nrow(m)))
  }

  probs <- seq(0, 1, length.out = as.integer(n_quantiles))

  # The quantile grid is the marginal. It is taken over observed values only,
  # so a feature's missingness is a separate property rather than a spike at
  # the bottom of its distribution.
  qgrid <- matrix(NA_real_, nrow = nrow(m), ncol = length(probs),
                  dimnames = list(rownames(m), NULL))
  abundance <- rep(NA_real_, nrow(m))
  spread <- rep(NA_real_, nrow(m))
  mu <- rep(NA_real_, nrow(m))
  sdv <- rep(NA_real_, nrow(m))
  kurt <- rep(NA_real_, nrow(m))
  missing <- rep(0, nrow(m))

  for (i in seq_len(nrow(m))) {
    v <- m[i, ]
    ok <- is.finite(v)
    missing[i] <- 1 - mean(ok)
    v <- v[ok]
    if (length(v) < 2) next
    qgrid[i, ] <- stats::quantile(v, probs, names = FALSE, type = 7)
    abundance[i] <- stats::median(v)
    spread[i] <- stats::IQR(v)
    mu[i] <- mean(v)
    sdv[i] <- stats::sd(v)
    if (is.finite(sdv[i]) && sdv[i] > 0) {
      kurt[i] <- chorale_excess_kurtosis((v - mu[i]) / sdv[i])
    }
  }

  usable <- is.finite(abundance)
  if (sum(usable) < 2) {
    rlang::abort("Fewer than two features carry enough observed values to profile.")
  }

  feature <- data.frame(
    feature_id = rownames(m),
    abundance = abundance,
    spread = spread,
    mean = mu,
    sd = sdv,
    kurtosis = kurt,
    missing = missing,
    stringsAsFactors = FALSE
  )

  bin <- chorale_abundance_bins(abundance)
  missing_curve <- chorale_bin_summary(bin, missing, "missing")
  mean_variance <- chorale_bin_summary(bin, sdv^2, "variance")
  mean_variance$mean <- chorale_bin_summary(bin, mu, "mean")$mean

  # The correlation structure enters as its eigenspectrum: how much of the
  # variance the leading directions carry is what an estimator recovering
  # components responds to, and it is comparable across feature spaces of
  # different size once normalised.
  x <- t(m[usable, , drop = FALSE])
  x <- scale(x)
  x[!is.finite(x)] <- 0
  d <- svd(x, nu = 0, nv = 0)$d
  ev <- d^2 / max(nrow(x) - 1L, 1L)
  ev <- ev / sum(ev)
  keep <- seq_len(min(as.integer(n_eigen), length(ev)))
  component_kurtosis <- rep(NA_real_, length(keep))
  if (length(keep) > 0) {
    pcs <- svd(x, nu = length(keep), nv = 0)$u
    for (k in seq_along(keep)) {
      component_kurtosis[k] <- chorale_excess_kurtosis(scale(pcs[, k]))
    }
  }

  design_cells <- NULL
  if (!is.null(col_data)) {
    design_cells <- chorale_design_cells(col_data, covariates)
  }

  structure(
    list(
      layer = layer %||% NA_character_,
      n_features = nrow(m),
      n_samples = ncol(m),
      probs = probs,
      quantiles = qgrid,
      feature = feature,
      missing_curve = missing_curve,
      mean_variance = mean_variance,
      eigenvalues = ev[keep],
      component_kurtosis = component_kurtosis,
      design_cells = design_cells,
      transform = chorale_transform(m)$applied
    ),
    class = "chorale_data_profile"
  )
}

#' Decile bins of feature abundance
#' @keywords internal
#' @noRd
chorale_abundance_bins <- function(abundance) {
  ok <- is.finite(abundance)
  breaks <- stats::quantile(abundance[ok], probs = seq(0, 1, by = 0.1),
                            names = FALSE, na.rm = TRUE)
  breaks <- unique(breaks)
  if (length(breaks) < 2) return(rep(1L, length(abundance)))
  b <- cut(abundance, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  b[is.na(b)] <- 1L
  as.integer(b)
}

#' Mean of a per-feature quantity within each abundance bin
#' @keywords internal
#' @noRd
chorale_bin_summary <- function(bin, value, name) {
  agg <- stats::aggregate(list(v = value), by = list(bin = bin),
                          FUN = function(z) mean(z[is.finite(z)]))
  out <- data.frame(bin = agg$bin, stringsAsFactors = FALSE)
  out[[name]] <- agg$v
  out
}

#' Observed design cells and their occupancy
#'
#' The margins of a real design, as the cells that were populated and how often.
#' Cells no animal occupies are absent, so a simulation drawn from this table
#' leaves them empty rather than filling a balanced grid.
#'
#' @keywords internal
#' @noRd
chorale_design_cells <- function(col_data, covariates = NULL) {
  d <- as.data.frame(col_data)
  d <- chorale_blank_to_na(d)
  keep <- covariates %||% setdiff(names(d), "sample_id")
  keep <- intersect(keep, names(d))
  if (length(keep) == 0) return(NULL)
  d <- d[, keep, drop = FALSE]
  # A covariate constant across the cohort contributes no contrast, and one that
  # is absent for every sample cannot be drawn from.
  varying <- vapply(d, function(v) length(unique(v[!is.na(v)])) > 0, logical(1))
  d <- d[, varying, drop = FALSE]
  if (ncol(d) == 0) return(NULL)
  for (nm in names(d)) {
    if (is.factor(d[[nm]])) d[[nm]] <- as.character(d[[nm]])
  }
  key <- do.call(paste, c(unname(as.list(d)), sep = "\r"))
  first <- !duplicated(key)
  cells <- d[first, , drop = FALSE]
  cells$n <- as.integer(table(key)[key[first]])
  rownames(cells) <- NULL
  cells
}

#' @export
print.chorale_data_profile <- function(x, ...) {
  cat("<chorale_data_profile>\n")
  cat("  layer:      ", x$layer, "\n", sep = "")
  cat("  dimensions: ", x$n_features, " features x ", x$n_samples,
      " samples\n", sep = "")
  cat("  transform:  ", x$transform, "\n", sep = "")
  cat("  missing:    ", format(round(100 * mean(x$feature$missing), 2),
                               nsmall = 2), "%\n", sep = "")
  cat("  kurtosis:   median ",
      format(round(stats::median(x$feature$kurtosis, na.rm = TRUE), 2),
             nsmall = 2), " over features\n", sep = "")
  cat("  leading eigenvalue: ",
      format(round(x$eigenvalues[1], 3), nsmall = 3),
      " of the total variance\n", sep = "")
  if (!is.null(x$design_cells)) {
    cat("  design:     ", nrow(x$design_cells), " populated cells over ",
        paste(setdiff(names(x$design_cells), "n"), collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

#' Compare two matrices on the properties a simulation must reproduce
#'
#' Places a simulated matrix beside the real one its profile was taken from,
#' on the six properties [chorale_data_profile()] measures. The comparison is
#' what a claim of realism is bounded to: agreement here says the simulated
#' matrix is indistinguishable from the real one on what the estimator reads,
#' and says nothing about the measurement process.
#'
#' @param simulated A `chorale_data_profile` of the simulated matrix.
#' @param reference A `chorale_data_profile` of the matrix the simulation was
#'   built to resemble.
#'
#' @returns A data frame with one row per property, carrying its value in each
#'   profile and the discrepancy between them.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 40, seed = 1)
#' a <- chorale_data_profile(sim$modalities[[1]])
#' b <- chorale_data_profile(sim$modalities[[2]])
#' chorale_profile_agreement(a, b)
chorale_profile_agreement <- function(simulated, reference) {
  if (!inherits(simulated, "chorale_data_profile") ||
      !inherits(reference, "chorale_data_profile")) {
    rlang::abort("Both arguments must be `chorale_data_profile` objects.")
  }

  qs <- c(0.1, 0.5, 0.9)
  quantile_of <- function(p) {
    v <- as.numeric(p$quantiles)
    stats::quantile(v[is.finite(v)], qs, names = FALSE)
  }
  qa <- quantile_of(simulated)
  qb <- quantile_of(reference)

  # The mean-variance relation is compared as its slope on the log scale, which
  # is what distinguishes a count-like matrix from an intensity-like one.
  slope_of <- function(p) {
    d <- p$mean_variance
    ok <- is.finite(d$mean) & is.finite(d$variance) & d$variance > 0
    if (sum(ok) < 3) return(NA_real_)
    x <- log(pmax(d$mean[ok] - min(d$mean[ok]) + 1, 1e-8))
    unname(stats::coef(stats::lm(log(d$variance[ok]) ~ x))[2])
  }

  # Missingness is compared where it is informative: the difference between the
  # least and most abundant features is the abundance dependence itself.
  missing_slope <- function(p) {
    d <- p$missing_curve
    if (nrow(d) < 2) return(NA_real_)
    d$missing[nrow(d)] - d$missing[1]
  }

  n_ev <- min(length(simulated$eigenvalues), length(reference$eigenvalues))
  ev_gap <- if (n_ev >= 1) {
    max(abs(cumsum(simulated$eigenvalues[seq_len(n_ev)]) -
              cumsum(reference$eigenvalues[seq_len(n_ev)])))
  } else {
    NA_real_
  }

  row <- function(property, a, b) {
    data.frame(property = property, simulated = as.numeric(a),
               reference = as.numeric(b), stringsAsFactors = FALSE)
  }
  out <- rbind(
    row("marginal_q10", qa[1], qb[1]),
    row("marginal_median", qa[2], qb[2]),
    row("marginal_q90", qa[3], qb[3]),
    row("mean_variance_slope", slope_of(simulated), slope_of(reference)),
    row("missing_fraction", mean(simulated$feature$missing),
        mean(reference$feature$missing)),
    row("missing_abundance_slope", missing_slope(simulated),
        missing_slope(reference)),
    row("eigenvalue_1", simulated$eigenvalues[1], reference$eigenvalues[1]),
    row("feature_kurtosis",
        stats::median(simulated$feature$kurtosis, na.rm = TRUE),
        stats::median(reference$feature$kurtosis, na.rm = TRUE)),
    row("component_kurtosis",
        stats::median(simulated$component_kurtosis, na.rm = TRUE),
        stats::median(reference$component_kurtosis, na.rm = TRUE))
  )
  out$simulated <- round(out$simulated, 4)
  out$reference <- round(out$reference, 4)
  out$discrepancy <- round(abs(out$simulated - out$reference), 4)
  # The cumulative eigenspectrum is a curve rather than a number, so its
  # discrepancy is the largest gap anywhere along it.
  out <- rbind(out, data.frame(
    property = "eigenspectrum_max_gap",
    simulated = NA_real_, reference = NA_real_,
    discrepancy = round(ev_gap, 4), stringsAsFactors = FALSE
  ))
  out
}
