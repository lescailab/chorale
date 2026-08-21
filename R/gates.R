#' Whether a collection meets the conditions the method rests on
#'
#' This function reports whether phenotype and shared-design effects can be
#' estimated, whether the design overlaps and has full rank, how many factors
#' are detectable, and how stable ICA is. Distributional summaries are retained
#' as diagnostics rather than treated as proof of cross-modal recovery.
#'
#' \describe{
#'   \item{non-Gaussianity diagnostic}{Are a modality's recovered components further from
#'     normal than components recovered the same way from data Gaussian by
#'     construction? Independent component analysis maximises non-Gaussianity,
#'     so the question is only meaningful against that calibration.}
#'   \item{modality-difference diagnostic}{How different are the component
#'     distributions across modalities? This is descriptive and does not gate
#'     matching.}
#'   \item{detectability}{How many components stand above what the same matrix
#'     produces with its covariance destroyed, and above the spiked-covariance
#'     threshold implied by its shape?}
#'   \item{anchor richness}{How many design strata are populated, and how many
#'     of those are shared across the modalities? The design is what stands in
#'     for matched individuals: no animal appears in two modalities, so the
#'     comparison is anchored on the strata the modalities have in common.}
#' }
#'
#' Distributional diagnostics are evaluated in R on components recovered by
#' the same estimator the free dimensions are recovered with. The estimator does not
#' claim unique latent-state recovery: its permutation inference is conditional
#' on stable, detectable fitted factors. Consequently cross-modality shape is
#' not an identification condition for the reported phenotype correspondence,
#' and non-Gaussianity is shown as an ICA diagnostic rather than a theorem gate.
#'
#' @param containers A named list of `SummarizedExperiment` objects, as
#'   [chorale_load()] returns, or of feature-by-sample matrices.
#' @param designs A named list of design tables. Required only where
#'   `containers` holds bare matrices; taken from the containers otherwise.
#' @param control A [chorale_control()] object. `alpha`, `n_init`,
#'   `n_factors_quantile` and `max_factors` are read from it.
#' @param n_factors Optional named integer vector fixing the component count
#'   per modality. Taken from the detectability condition where absent.
#' @param n_surrogate Integer number of Gaussian surrogates the non-Gaussianity
#'   condition is calibrated against.
#' @param n_perm Integer number of permutations for the detectability
#'   condition.
#' @param seed Integer seed.
#'
#' @returns An object of class `chorale_gates`: a list of one data frame per
#'   condition, plus `n_factors`, the component count each modality supports.
#'
#' @examples
#' \dontrun{
#' sim <- chorale_simulate(n_modalities = 2, n_samples = 60, n_features = 200)
#' containers <- lapply(sim$modalities, chorale_load, design = sim$design)
#' chorale_gates(containers)
#' }
#' @export
chorale_gates <- function(containers, designs = NULL,
                          control = chorale_control(),
                          n_factors = NULL, n_surrogate = 100L,
                          n_perm = 200L, seed = 1L) {
  if (!is.list(containers) || length(containers) < 1) {
    rlang::abort("`containers` must be a named list of modalities.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }

  assays <- lapply(containers, function(x) {
    if (inherits(x, "SummarizedExperiment")) SummarizedExperiment::assay(x) else as.matrix(x)
  })
  if (is.null(designs)) {
    designs <- lapply(containers, function(x) {
      if (!inherits(x, "SummarizedExperiment")) {
        rlang::abort("`designs` is required where `containers` holds matrices.")
      }
      as.data.frame(SummarizedExperiment::colData(x))
    })
  }

  # Samples by features, centred and scaled, which is what both the estimator
  # and the surrogate construction expect.
  xs <- lapply(assays, function(a) {
    x <- t(as.matrix(a))
    x[!is.finite(x)] <- 0
    x <- scale(x)
    x[, apply(x, 2, function(col) all(is.finite(col))), drop = FALSE]
  })

  detect <- chorale_gate_detectability(xs, n_perm = n_perm,
                                       quantile = control$n_factors_quantile,
                                       max_factors = control$max_factors,
                                       seed = seed)
  if (is.null(n_factors)) {
    n_factors <- stats::setNames(detect$n_factors, detect$modality)
  }

  ica_fn <- chorale_gate_ica(n_init = control$n_init, consensus = control$consensus)
  ng <- lapply(names(xs), function(m) {
    chorale_gate_nongaussianity(
      xs[[m]], as.integer(n_factors[[m]]), m, ica_fn,
      seed = as.integer(seed), n_surrogate = as.integer(n_surrogate),
      alpha = control$alpha)
  })
  names(ng) <- names(xs)

  gate_fits <- lapply(names(xs), function(m) {
    chorale_ica(xs[[m]], as.integer(n_factors[[m]]),
                n_init = control$n_init, seed = as.integer(seed),
                consensus = control$consensus)
  })
  names(gate_fits) <- names(xs)
  sources <- lapply(gate_fits, function(x) unname(as.matrix(x$scores)))
  names(sources) <- names(xs)
  difference <- chorale_gate_modality_difference(sources, alpha = control$alpha)

  design_estimability <- chorale_gate_design(designs, control)
  factor_stability <- do.call(rbind, lapply(names(gate_fits), function(m) {
    s <- gate_fits[[m]]$stability
    data.frame(
      modality = m,
      selected_initialisation = gate_fits[[m]]$selected_init,
      mean_subspace_agreement = mean(s$subspace, na.rm = TRUE),
      minimum_subspace_agreement = min(s$subspace, na.rm = TRUE),
      n_failed = sum(!is.finite(s$objective)),
      stringsAsFactors = FALSE)
  }))

  out <- list(
    non_gaussianity = do.call(rbind, lapply(ng, `[[`, 1L)),
    non_gaussianity_components = do.call(rbind, lapply(ng, `[[`, 2L)),
    modality_difference = difference,
    design_estimability = design_estimability,
    factor_stability = factor_stability,
    detectability = detect,
    anchor_richness = chorale_gate_anchors(designs),
    n_factors = n_factors
  )
  class(out) <- "chorale_gates"
  out
}

#' Phenotype estimability and joint design rank
#' @keywords internal
#' @noRd
chorale_gate_design <- function(designs, control = chorale_control()) {
  designs <- lapply(designs, as.data.frame)
  spec <- try(chorale_resolve_signature(
    designs, phenotype_column = control$phenotype_column,
    phenotype_reference = control$phenotype_reference,
    profile_covariates = control$profile_covariates), silent = TRUE)
  if (inherits(spec, "try-error")) {
    return(data.frame(modality = names(designs), phenotype_estimable = FALSE,
                      full_rank = FALSE, n_complete = NA_integer_,
                      reason = as.character(spec), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(names(designs), function(m) {
    x <- chorale_signature_matrix(designs[[m]], spec)$x
    ok <- stats::complete.cases(x)
    data.frame(
      modality = m,
      phenotype_estimable = spec$phenotype %in% spec$covariates,
      full_rank = sum(ok) > ncol(x) && qr(x[ok, , drop = FALSE])$rank == ncol(x),
      n_complete = sum(ok), reason = "", stringsAsFactors = FALSE)
  }))
}

#' Anderson--Darling statistic against a fitted normal distribution
#' @keywords internal
#' @noRd
chorale_ad_normal <- function(x) {
  x <- sort(as.numeric(x[is.finite(x)]))
  n <- length(x)
  if (n < 5L || stats::sd(x) == 0) return(NA_real_)
  z <- (x - mean(x)) / stats::sd(x)
  p <- pmin(1 - 1e-12, pmax(1e-12, stats::pnorm(z)))
  i <- seq_len(n)
  -n - mean((2 * i - 1) * (log(p) + log(1 - rev(p))))
}

#' R-native non-Gaussianity diagnostic calibrated by Gaussian surrogates
#' @keywords internal
#' @noRd
chorale_gate_nongaussianity <- function(x, k, modality, estimator,
                                        seed = 1L, n_surrogate = 100L,
                                        alpha = 0.05) {
  observed <- estimator(x, k, seed)
  obs <- apply(observed, 2, chorale_ad_normal)
  set.seed(seed)
  surrogate <- numeric(n_surrogate)
  for (b in seq_len(n_surrogate)) {
    # Preserve location, scale and matrix dimensions while removing
    # non-Gaussian marginal structure.
    gx <- matrix(stats::rnorm(length(x)), nrow = nrow(x), ncol = ncol(x))
    gs <- estimator(scale(gx), k, seed + b)
    surrogate[b] <- stats::median(apply(gs, 2, chorale_ad_normal), na.rm = TRUE)
  }
  value <- stats::median(obs, na.rm = TRUE)
  p <- (1 + sum(surrogate >= value)) / (1 + n_surrogate)
  summary <- data.frame(
    modality = modality, median_A2_observed = value,
    median_A2_surrogate = stats::median(surrogate, na.rm = TRUE),
    p_value = p,
    verdict = if (p < alpha) "non-Gaussian pattern detected" else
      "non-Gaussian pattern not detected",
    role = "diagnostic; does not gate matching", gates_matching = FALSE,
    stringsAsFactors = FALSE)
  detail <- data.frame(modality = modality, component = seq_along(obs),
                       A2 = obs, stringsAsFactors = FALSE)
  list(summary, detail)
}

#' R-native cross-modality distribution diagnostic
#' @keywords internal
#' @noRd
chorale_gate_modality_difference <- function(sources, alpha = 0.05) {
  mods <- names(sources)
  rows <- list()
  for (i in seq_along(mods)) for (j in seq_along(mods)) {
    if (j <= i) next
    a <- as.numeric(scale(sources[[i]]))
    b <- as.numeric(scale(sources[[j]]))
    ks <- suppressWarnings(stats::ks.test(a, b, exact = FALSE))
    component_p <- numeric()
    for (ca in seq_len(ncol(sources[[i]]))) {
      for (cb in seq_len(ncol(sources[[j]]))) {
        component_p <- c(component_p, suppressWarnings(stats::ks.test(
          sources[[i]][, ca], sources[[j]][, cb], exact = FALSE)$p.value))
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      pair = paste(mods[i], mods[j], sep = " vs "),
      pooled_KS_D = unname(ks$statistic), pooled_KS_p = ks$p.value,
      pct_pairs_indistinguishable = mean(component_p >= alpha) * 100,
      verdict = if (ks$p.value < alpha) "distribution difference detected" else
        "distribution difference not detected",
      role = "diagnostic; does not gate matching", gates_matching = FALSE,
      stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}


#' The component estimator, wrapped for the gate module to call
#'
#' Returns a function of `(x, k, seed)` giving standardised sources, so the
#' components the conditions are judged on are recovered by the same estimator
#' that will consume them.
#' @keywords internal
#' @noRd
chorale_gate_ica <- function(n_init = 20L, consensus = TRUE) {
  function(x, k, seed) {
    x <- as.matrix(x)
    fit <- chorale_ica(x, n_factors = as.integer(k), n_init = n_init,
                       seed = as.integer(seed), consensus = consensus)
    unname(as.matrix(fit$scores))
  }
}


#' The detectability condition
#'
#' Two readings of the same question. Parallel analysis compares each
#' eigenvalue with what the matrix produces once each feature is permuted
#' independently, destroying covariance while preserving every marginal. The
#' spiked-covariance threshold asks the same of the shape of the matrix alone,
#' at `sigma^2 (1 + sqrt(gamma))` with `sigma^2` taken as the median bulk
#' eigenvalue.
#' @keywords internal
#' @noRd
chorale_gate_detectability <- function(xs, n_perm = 200L, quantile = 0.95,
                                       max_factors = 20L, seed = 1L) {
  rows <- lapply(names(xs), function(m) {
    x <- xs[[m]]
    n <- nrow(x)
    p <- ncol(x)
    ev <- svd(x, nu = 0, nv = 0)$d^2 / max(n - 1L, 1L)
    sigma2 <- stats::median(ev)
    gamma <- p / n
    threshold <- sigma2 * (1 + sqrt(gamma))
    k <- chorale_n_factors(x, n_perm = n_perm, quantile = quantile,
                           max_factors = max_factors, seed = seed)
    data.frame(
      modality = m,
      n_samples = n,
      n_features = p,
      gamma_p_over_n = round(gamma, 3),
      n_spikes_above_threshold = sum(ev > threshold),
      n_factors = k,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' The anchor-richness condition
#'
#' Reports, at each level of coarsening, how many design strata each modality
#' populates and how many of those every modality shares. Nothing beyond
#' `sample_id` and `phenotype` is required of a design, so the coarsening is
#' built from whatever covariates the designs actually carry.
#' @keywords internal
#' @noRd
chorale_gate_anchors <- function(designs) {
  designs <- lapply(designs, as.data.frame)
  shared_cols <- Reduce(intersect, lapply(designs, colnames))
  candidates <- setdiff(shared_cols, c("sample_id", "modality", "cohort"))
  varies <- vapply(candidates, function(cv) {
    all(vapply(designs, function(d) {
      length(unique(stats::na.omit(d[[cv]]))) >= 2
    }, logical(1)))
  }, logical(1))
  candidates <- candidates[varies]
  if (!"phenotype" %in% candidates) {
    rlang::abort("No shared covariate varies in every modality; nothing can be anchored on.")
  }
  # Coarsen from every shared covariate down to the phenotype alone, so a
  # collection that cannot support the finest anchoring is still reported at
  # the level it can support.
  others <- setdiff(candidates, "phenotype")
  levels_of <- lapply(seq(length(others), 0), function(k) {
    c("phenotype", utils::head(others, k))
  })

  rows <- list()
  for (keys in levels_of) {
    label <- paste(keys, collapse = " x ")
    per_mod <- lapply(designs, function(d) {
      sub <- d[stats::complete.cases(d[, keys, drop = FALSE]), keys, drop = FALSE]
      if (nrow(sub) == 0) return(character())
      unique(do.call(paste, c(sub, sep = "|")))
    })
    shared <- Reduce(intersect, per_mod)
    for (m in names(designs)) {
      rows[[length(rows) + 1L]] <- data.frame(
        coarsening = label, modality = m,
        n_strata = length(per_mod[[m]]),
        n_strata_shared = length(shared),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


#' Turn the gate module's records into a data frame
#'
#' Gate implementations return named R records. A record missing a field is
#' filled rather than dropped, so one modality cannot silently shorten a table.
#' @keywords internal
#' @noRd
chorale_records_to_df <- function(records) {
  if (is.null(records) || length(records) == 0) return(data.frame())
  if (!is.null(names(records)) && !is.list(records[[1]])) records <- list(records)
  records <- unname(records)
  cols <- unique(unlist(lapply(records, names)))
  out <- lapply(cols, function(cn) {
    vals <- lapply(records, function(r) {
      v <- r[[cn]]
      if (is.null(v) || length(v) == 0) NA else v[[1]]
    })
    if (all(vapply(vals, function(v) is.numeric(v) || is.na(v), logical(1)))) {
      as.numeric(unlist(vals))
    } else {
      as.character(unlist(lapply(vals, as.character)))
    }
  })
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE)
}


#' @export
print.chorale_gates <- function(x, ...) {
  cat("<chorale_gates>\n")
  cat("\nnon-Gaussianity diagnostic\n")
  print(x$non_gaussianity[, c("modality", "median_A2_observed",
                              "median_A2_surrogate", "p_value", "verdict")],
        row.names = FALSE)
  cat("\nmodality-distribution diagnostic\n")
  print(x$modality_difference[, c("pair", "pooled_KS_D", "pooled_KS_p",
                                  "pct_pairs_indistinguishable", "verdict")],
        row.names = FALSE)
  cat("\ndetectability\n")
  print(x$detectability, row.names = FALSE)
  cat("\nphenotype estimability and design rank\n")
  print(x$design_estimability, row.names = FALSE)
  cat("\nfactor stability\n")
  print(x$factor_stability, row.names = FALSE)
  cat("\nanchor richness\n")
  print(utils::head(x$anchor_richness, 12), row.names = FALSE)
  invisible(x)
}
