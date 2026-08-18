#' Permutation calibration and stability diagnostics
#'
#' Every reported result is accompanied by the control that could have
#' falsified it. Three run here.
#'
#' Permuting the phenotype label within cohort and stratum breaks its relation
#' to the data while leaving the design intact, so the anchor agreement of a
#' genuine shared factor should exceed what permuted labels produce.
#'
#' Permuting the modality label pools the samples and reassigns them, which
#' destroys the modality-specific structure the identification argument
#' consumes. Matches that survive it were never resting on that structure.
#'
#' Stability across initialisations is reported because independent component
#' analysis is non-convex: a factor recovered at one initialisation and not at
#' another is a draw, not an estimate.
#'
#' @param fit A `chorale_fit` object, as returned by [chorale_fit()].
#' @param containers The modality containers the fit was built from, needed to
#'   refit under permutation.
#' @param n_permutations Number of label permutations.
#' @param n_init Initialisations per refit.
#' @param seed Integer seed.
#'
#' @returns An object of class `chorale_null` holding the permutation nulls,
#'   the modality-shuffle null, and per-factor stability.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' chorale_null(fit, containers, n_permutations = 3, n_init = 2)
chorale_null <- function(fit, containers, n_permutations = 100L,
                         n_init = 5L, seed = 1L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }

  significant <- if (nrow(fit$matches) > 0) {
    fit$matches[fit$matches$significant, , drop = FALSE]
  } else {
    fit$matches
  }
  observed <- if (nrow(significant) > 0) {
    max(significant$anchor_agreement)
  } else if (nrow(fit$matches) > 0) {
    max(fit$matches$anchor_agreement)
  } else {
    NA_real_
  }
  n_observed <- nrow(significant)

  # Phenotype permutation, stratified so the design is preserved.
  phenotype_null <- numeric(n_permutations)
  for (i in seq_len(n_permutations)) {
    permuted <- containers
    for (m in names(permuted)) {
      cd <- SummarizedExperiment::colData(permuted[[m]])
      d <- chorale_add_age_bin(as.data.frame(cd))
      parts <- list(
        as.character(d[["cohort"]] %||% rep("all", nrow(d))),
        as.character(d[["age_bin"]] %||% rep("all", nrow(d))),
        as.character(d[["sex"]] %||% rep("all", nrow(d)))
      )
      parts <- lapply(parts, function(v) {
        v[is.na(v)] <- "unknown"
        v
      })
      strata <- do.call(paste, c(parts, sep = "|"))
      set.seed(seed + i)
      # Permute the phenotype within each stratum, so the design is preserved
      # and only its relation to the data is broken.
      for (lv in unique(strata)) {
        idx <- which(strata == lv)
        if (length(idx) > 1) d$phenotype[idx] <- sample(d$phenotype[idx])
      }
      SummarizedExperiment::colData(permuted[[m]]) <- S4Vectors::DataFrame(d)
    }
    refit <- chorale_fit(permuted, n_factors = fit$n_factors,
                         n_init = n_init, strata_keys = fit$strata_keys,
                         seed = seed + i)
    phenotype_null[i] <- if (nrow(refit$matches) > 0) {
      max(refit$matches$anchor_agreement)
    } else {
      0
    }
  }

  # Modality shuffle: pool the samples and reassign them across modalities,
  # which is only possible where the modalities share a feature space.
  modality_null <- chorale_modality_shuffle(containers, fit, n_init, seed)

  stability <- do.call(rbind, lapply(fit$modalities, function(m) {
    s <- fit$fits[[m]]$stability
    data.frame(
      modality = m,
      n_init = nrow(s),
      objective_median = stats::median(s$objective, na.rm = TRUE),
      objective_sd = stats::sd(s$objective, na.rm = TRUE),
      objective_cv = stats::sd(s$objective, na.rm = TRUE) /
        abs(stats::median(s$objective, na.rm = TRUE)),
      n_failed = sum(is.na(s$objective)),
      stringsAsFactors = FALSE
    )
  }))

  p_phenotype <- if (is.na(observed)) {
    NA_real_
  } else {
    (1 + sum(phenotype_null >= observed)) / (1 + n_permutations)
  }

  structure(
    list(
      observed_agreement = observed,
      n_observed_matches = n_observed,
      phenotype_null = phenotype_null,
      p_phenotype = p_phenotype,
      modality_null = modality_null,
      stability = stability,
      n_permutations = n_permutations
    ),
    class = "chorale_null"
  )
}

#' Refit after reassigning samples across modalities
#' @keywords internal
#' @noRd
chorale_modality_shuffle <- function(containers, fit, n_init, seed) {
  features <- lapply(containers, function(se) rownames(SummarizedExperiment::assay(se)))
  common <- Reduce(intersect, features)
  if (length(common) < 10) {
    # Disjoint feature spaces cannot be pooled, so the shuffle is undefined
    # rather than passed.
    return(list(applicable = FALSE, agreement = NA_real_,
                reason = "modalities share fewer than ten features"))
  }

  mats <- lapply(containers, function(se) {
    SummarizedExperiment::assay(se)[common, , drop = FALSE]
  })
  designs <- lapply(containers, function(se) {
    chorale_add_age_bin(as.data.frame(SummarizedExperiment::colData(se)))
  })
  pooled <- do.call(cbind, mats)
  pooled_design <- do.call(rbind, lapply(names(designs), function(m) {
    d <- designs[[m]]
    d[, intersect(colnames(d), colnames(designs[[1]])), drop = FALSE]
  }))
  colnames(pooled) <- pooled_design$sample_id <-
    make.unique(as.character(pooled_design$sample_id))

  set.seed(seed)
  assignment <- sample(rep(names(containers), length.out = ncol(pooled)))
  shuffled <- lapply(names(containers), function(m) {
    keep <- which(assignment == m)
    chorale_load(pooled[, keep, drop = FALSE], pooled_design[keep, , drop = FALSE])
  })
  names(shuffled) <- names(containers)

  refit <- try(
    chorale_fit(shuffled, n_factors = fit$n_factors, n_init = n_init,
                strata_keys = fit$strata_keys, seed = seed),
    silent = TRUE
  )
  if (inherits(refit, "try-error")) {
    return(list(applicable = FALSE, agreement = NA_real_,
                reason = "refit failed under modality shuffle"))
  }
  list(
    applicable = TRUE,
    agreement = if (nrow(refit$matches) > 0) max(refit$matches$anchor_agreement) else 0,
    n_matches = nrow(refit$matches),
    reason = NA_character_
  )
}

#' @export
print.chorale_null <- function(x, ...) {
  cat("<chorale_null>\n")
  cat("  observed matches:", x$n_observed_matches, "\n")
  cat("  strongest observed agreement:", round(x$observed_agreement, 3), "\n")
  cat("  phenotype permutation p-value:", signif(x$p_phenotype, 3),
      sprintf("(%d permutations)\n", x$n_permutations))
  if (isTRUE(x$modality_null$applicable)) {
    cat("  modality shuffle agreement:", round(x$modality_null$agreement, 3), "\n")
  } else {
    cat("  modality shuffle: not applicable,", x$modality_null$reason, "\n")
  }
  cat("  initialisation stability (coefficient of variation):\n")
  for (i in seq_len(nrow(x$stability))) {
    cat(sprintf("    %-12s %.3f\n", x$stability$modality[i],
                x$stability$objective_cv[i]))
  }
  invisible(x)
}
