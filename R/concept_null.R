#' Permutation calibration and stability diagnostics
#'
#' Every reported result is accompanied by the check that would have caught it
#' if it were an artefact. Three run here.
#'
#' The phenotype null holds the design fixed and exchanges the part of each
#' concept score the adjusting covariates do not explain, so the phenotype keeps
#' the relation to those covariates that the data give it. The encoder is never
#' refitted.
#'
#' The modality shuffle pools the samples and deals them out again. It is a
#' description of how much of the result survives that treatment rather than a
#' test, because pooling restricts the collection to the features every modality
#' shares and the reallocation preserves neither each modality's sample count
#' nor its design composition.
#'
#' Stability across initialisations is reported because independent component
#' analysis is non-convex: a free dimension recovered at one initialisation and
#' not at another is a draw, not an estimate.
#'
#' @param fit A `chorale_concept_fit`.
#' @param containers The modality containers the fit was built from, needed only
#'   by the modality shuffle.
#' @param n_permutations Number of permutations. `NULL` reuses the permutations
#'   the fit already paid for.
#' @param n_init Retained for compatibility; the encoder is not refitted.
#' @param seed Integer seed.
#' @param ... Passed to the method.
#'
#' @returns An object of class `chorale_concept_null` holding the phenotype
#'   null, the modality-shuffle description, per-modality stability, and one row
#'   per control with the smallest p-value it can attain.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
#'                            n_permutations = 99, n_init = 2)
#' chorale_null(fit, fx$containers, n_shuffles = 3)
chorale_null <- function(fit, containers = NULL, n_permutations = 100L,
                         n_init = 5L, seed = 1L, ...) {
  UseMethod("chorale_null")
}

#' @export
chorale_null.default <- function(fit, containers = NULL, n_permutations = 100L,
                                 n_init = 5L, seed = 1L, ...) {
  rlang::abort("`fit` must be a chorale_concept_fit object.")
}

#' @export

#' @rdname chorale_null
#' @param n_shuffles Reassignments of samples across modalities forming the
#'   modality-shuffle null. Zero skips that control.
#' @export
chorale_null.chorale_concept_fit <- function(fit, containers = NULL,
                                             n_permutations = NULL,
                                             n_init = 5L, seed = 1L,
                                             n_shuffles = 20L, ...) {
  evidence <- fit$evidence
  if (!is.null(n_permutations)) {
    # Recalibrating costs nothing but the permutations themselves: the concept
    # scores are already computed and the permutation acts on the design.
    evidence <- chorale_concept_evidence(
      fit$encoding, n_permutations = n_permutations,
      control = fit$control, seed = seed)
  }
  observed <- chorale_best_concept(evidence)
  n_perm <- evidence$n_permutations

  p_phenotype <- if (!is.finite(observed) || n_perm == 0L) {
    NA_real_
  } else {
    (1 + sum(evidence$null >= observed, na.rm = TRUE)) / (1 + n_perm)
  }

  modality_null <- if (n_shuffles > 0 && !is.null(containers)) {
    chorale_concept_modality_shuffle(containers, fit, n_shuffles, seed)
  } else {
    list(applicable = FALSE, statistic = NA_real_, p_value = NA_real_,
         n_shuffles = 0L,
         reason = if (is.null(containers)) {
           "the containers the fit was built from were not supplied"
         } else {
           "no shuffles were requested"
         })
  }

  stability <- do.call(rbind, lapply(fit$modalities, function(m) {
    s <- fit$encoding$encodings[[m]]$stability
    sub <- if (is.data.frame(s) && "subspace" %in% colnames(s)) {
      s$subspace
    } else {
      NA_real_
    }
    data.frame(
      modality = m,
      n_free = fit$encoding$encodings[[m]]$n_free,
      n_init = if (is.data.frame(s)) nrow(s) else 0L,
      subspace_agreement = round(chorale_finite_mean(sub), 3),
      subspace_min = round(chorale_finite_min(sub), 3),
      stringsAsFactors = FALSE
    )
  }))

  structure(
    list(
      observed_statistic = observed,
      n_supported = if (nrow(evidence$joint) > 0) {
        sum(evidence$joint$significant)
      } else {
        0L
      },
      phenotype_null = evidence$null,
      p_phenotype = p_phenotype,
      modality_null = modality_null,
      stability = stability,
      n_permutations = n_perm,
      alpha = evidence$alpha,
      controls = chorale_concept_control_table(observed, p_phenotype, n_perm,
                                               modality_null, stability)
    ),
    class = c("chorale_concept_null", "chorale_null")
  )
}

#' The strongest statistic anywhere in the vocabulary
#' @keywords internal
#' @noRd
chorale_best_concept <- function(evidence) {
  j <- evidence$joint
  if (is.null(j) || nrow(j) == 0) return(0)
  v <- abs(j$joint_z[is.finite(j$joint_z)])
  if (length(v) == 0) 0 else max(v)
}

#' @keywords internal
#' @noRd
chorale_finite_mean <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) NA_real_ else mean(v)
}

#' Re-encode after reassigning samples across modalities
#'
#' Pooling the samples and dealing them out again asks whether a concept is a
#' property of the biology rather than of the assay it was measured with. It is
#' reported as a description and not as a test, because the shuffled collections
#' are not comparable with the fit they would have to be read against.
#'
#' Three differences make them incomparable. Pooling keeps only the features the
#' modalities share, so the shuffled fit is scored on a smaller feature space and
#' often on a smaller vocabulary than the observed one. The reallocation does not
#' preserve each modality's original sample count, nor the phenotype and
#' covariate composition within it, so the precision of the regressions changes.
#' And the observed statistic it would be compared against was computed on the
#' full assays. A p-value from that comparison would not mean that the observed
#' result survived a shuffle of the assay labels.
#'
#' An inferential form of this control needs the observed statistic refitted on
#' exactly the common-feature subset, feature spaces that can be pooled, and
#' modality labels permuted within declared design cells while each modality
#' keeps its original count in every cell.
#'
#' @keywords internal
#' @noRd
chorale_concept_modality_shuffle <- function(containers, fit, n_shuffles, seed) {
  features <- lapply(containers, function(se) {
    rownames(SummarizedExperiment::assay(se))
  })
  # Pooling requires one feature space, so the collection is restricted to the
  # features every modality measures. This is the first of the three reasons the
  # result is descriptive: below ten shared features there is nothing to pool
  # and the control is reported as inapplicable rather than run on a handful.
  common <- Reduce(intersect, features)
  if (length(common) < 10) {
    return(list(applicable = FALSE, statistic = NA_real_, p_value = NA_real_,
                n_shuffles = 0L,
                reason = "modalities share fewer than ten features"))
  }

  mats <- lapply(containers, function(se) {
    SummarizedExperiment::assay(se)[common, , drop = FALSE]
  })
  designs <- lapply(containers, function(se) {
    as.data.frame(SummarizedExperiment::colData(se))
  })
  shared_columns <- Reduce(intersect, lapply(designs, colnames))
  pooled <- do.call(cbind, mats)
  pooled_design <- do.call(rbind, lapply(designs, function(d) {
    d[, shared_columns, drop = FALSE]
  }))
  # Two modalities may reuse an identifier, and the pooled table would then
  # carry two rows under one name, which chorale_load() refuses.
  colnames(pooled) <- pooled_design$sample_id <-
    make.unique(as.character(pooled_design$sample_id))

  observed <- chorale_best_concept(fit$evidence)
  null <- rep(NA_real_, n_shuffles)
  for (b in seq_len(n_shuffles)) {
    set.seed(seed + b)
    # Labels are dealt out evenly and without regard to the design, which is the
    # second and third reasons the result is descriptive: a shuffled modality
    # keeps neither its original sample count nor its phenotype composition, so
    # the precision of the regressions it supports is not the precision the
    # observed fit had.
    assignment <- sample(rep(names(containers), length.out = ncol(pooled)))
    shuffled <- lapply(names(containers), function(m) {
      keep <- which(assignment == m)
      chorale_load(pooled[, keep, drop = FALSE],
                   pooled_design[keep, , drop = FALSE])
    })
    names(shuffled) <- names(containers)
    refit <- try(chorale_concept_fit(
      shuffled, fit$concepts$sets,
      n_free = 0, n_permutations = 1L,
      min_features = fit$concepts$min_features,
      min_modalities = fit$concepts$min_modalities,
      control = fit$control, seed = seed + b), silent = TRUE)
    if (!inherits(refit, "try-error")) {
      null[b] <- chorale_best_concept(refit$evidence)
    }
  }
  null <- null[is.finite(null)]
  if (length(null) == 0) {
    return(list(applicable = FALSE, statistic = NA_real_, p_value = NA_real_,
                n_shuffles = 0L,
                reason = "every modality shuffle failed to re-encode"))
  }
  list(
    applicable = FALSE,
    statistic = round(stats::median(null), 4),
    observed = observed,
    p_value = NA_real_,
    n_shuffles = length(null),
    reason = paste0(
      "descriptive only: pooling restricts the feature space and modality ",
      "labels are not exchangeable under the fitted design")
  )
}

#' One row per control, with the smallest p-value it can attain
#'
#' A control that cannot reach the threshold it is read against says nothing,
#' and the reader cannot tell that from the p-value alone. Stating the smallest
#' attainable value beside each result is what makes an uninformative control
#' visible as one.
#'
#' @keywords internal
#' @noRd
chorale_concept_control_table <- function(observed, p_phenotype, n_perm,
                                          modality_null, stability) {
  rows <- list(
    data.frame(
      control = "phenotype permutation",
      applicable = TRUE,
      observed = round(observed, 4),
      null_value = NA_real_,
      p_value = signif(p_phenotype, 3),
      smallest_attainable_p = signif(1 / (n_perm + 1), 3),
      n_resamples = n_perm,
      reason = NA_character_,
      stringsAsFactors = FALSE),
    data.frame(
      control = "modality shuffle",
      applicable = isTRUE(modality_null$applicable),
      observed = round(observed, 4),
      null_value = modality_null$statistic,
      p_value = signif(modality_null$p_value, 3),
      smallest_attainable_p = NA_real_,
      n_resamples = modality_null$n_shuffles %||% 0L,
      reason = modality_null$reason %||% NA_character_,
      stringsAsFactors = FALSE)
  )
  # Stability is a diagnostic rather than a test, so it carries no p-value and
  # reports the weakest component instead.
  rows[[3]] <- data.frame(
    control = "initialisation stability",
    applicable = any(is.finite(stability$subspace_agreement)),
    observed = round(chorale_finite_min(stability$subspace_min), 3),
    null_value = NA_real_,
    p_value = NA_real_,
    smallest_attainable_p = NA_real_,
    n_resamples = sum(stability$n_init),
    reason = if (any(is.finite(stability$subspace_agreement))) {
      NA_character_
    } else {
      "no free dimensions were fitted"
    },
    stringsAsFactors = FALSE)
  do.call(rbind, rows)
}

#' @export
print.chorale_concept_null <- function(x, ...) {
  cat("<chorale_concept_null>\n")
  cat("  concepts supported at q <=", x$alpha, ":", x$n_supported, "\n")
  cat("  strongest concept statistic:", round(x$observed_statistic, 3), "\n")
  for (i in seq_len(nrow(x$controls))) {
    r <- x$controls[i, ]
    if (!r$applicable && is.finite(r$null_value)) {
      cat(sprintf("  %-24s median %.3f across %d shuffles; %s\n", r$control,
                  r$null_value, r$n_resamples, r$reason))
    } else if (!r$applicable) {
      cat(sprintf("  %-24s not applicable: %s\n", r$control, r$reason))
    } else if (is.na(r$p_value)) {
      cat(sprintf("  %-24s weakest component %.3f\n", r$control, r$observed))
    } else {
      cat(sprintf("  %-24s p = %s (smallest attainable %s, %d resamples)\n",
                  r$control, format(r$p_value), format(r$smallest_attainable_p),
                  r$n_resamples))
    }
  }
  invisible(x)
}
