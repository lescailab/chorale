#' False-discovery control across the levels of the search
#'
#' A run tests many things at once: every factor against the design, every
#' programme against its joint null, and every programme again on curated
#' biology. Reporting the raw p-values of the survivors would control nothing,
#' since the survivors were selected. The permutation nulls already account for
#' the selection *within* a level, by keeping the best value the procedure could
#' reach; what remains is the multiplicity *across* the objects reported at each
#' level, which this controls.
#'
#' Levels are corrected separately rather than pooled, because they answer
#' different questions on different objects and pooling them would penalise a
#' programme for the number of factors that happened to be fitted. Within a
#' level the Benjamini-Hochberg procedure gives the false discovery rate
#' (Benjamini and Hochberg 1995), which is the appropriate target when the aim
#' is to rank candidates for follow-up rather than to protect a single
#' confirmatory claim.
#'
#' @param fit A `chorale_fit` object.
#' @param alpha The false-discovery rate to control at.
#' @param associations A precomputed table of factor-covariate associations, as
#'   `chorale_report()` builds internally. It costs a permutation test per
#'   factor and covariate, so a caller that already holds it passes it rather
#'   than paying for it again.
#'
#' @returns A data frame with one row per tested object, carrying `level`
#'   (`factor`, `programme` or `pathway`), the object's identifier, its
#'   `p_value`, the `q_value` within its level, and whether it is `significant`
#'   at `alpha`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' chorale_fdr(fit)
chorale_fdr <- function(fit, alpha = 0.05, associations = NULL) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  rows <- list()

  # Level one: each factor's association with the phenotype, within modality.
  # The table costs a permutation test per factor and covariate, so a caller
  # that already has it passes it rather than paying for it again.
  assoc <- associations %||% try(chorale_association_table(fit), silent = TRUE)
  if (!inherits(assoc, "try-error") && is.data.frame(assoc) && nrow(assoc) > 0) {
    a <- assoc[assoc$covariate == "phenotype", , drop = FALSE]
    if (nrow(a) > 0) {
      rows[[length(rows) + 1L]] <- data.frame(
        level = "factor",
        object = paste(a$modality, a$factor),
        p_value = a$p_permutation,
        stringsAsFactors = FALSE
      )
    }
  }

  # Level two: each programme's joint evidence.
  pg <- fit$programmes
  if (!is.null(pg) && nrow(pg) > 0) {
    u <- pg[!duplicated(pg$programme), , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      level = "programme",
      object = u$programme,
      p_value = u$joint_p,
      stringsAsFactors = FALSE
    )
  }

  # Level three: each programme's corroboration on curated biology.
  pe <- fit$pathway_evidence
  if (!is.null(pe) && nrow(pe) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(
      level = "pathway",
      object = pe$programme,
      p_value = pe$pathway_p,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  out <- out[is.finite(out$p_value), , drop = FALSE]
  if (nrow(out) == 0) return(data.frame())

  out$q_value <- NA_real_
  for (lv in unique(out$level)) {
    i <- out$level == lv
    out$q_value[i] <- stats::p.adjust(out$p_value[i], method = "BH")
  }
  out$significant <- out$q_value < alpha
  out[order(out$level, out$q_value, out$p_value), , drop = FALSE]
}

#' Whether a programme needed more than one modality
#'
#' The claim that integration produced a result is only worth making if the
#' result is stronger than what any one modality could have produced alone. Two
#' comparisons decide it, and a programme should carry both.
#'
#' The **best single modality** is the strongest phenotype association any one
#' member factor achieves on its own. A programme whose joint evidence does not
#' exceed it was visible in one modality and merely accompanied by the others.
#'
#' **Leave-one-modality-out** removes each member in turn. Where the programme
#' spans two modalities the remainder is a single modality and no joint
#' statistic exists, so the comparison rests on the single-modality column
#' alone.
#'
#' @param fit A `chorale_fit` object.
#' @param programmes Output of [chorale_programmes()]; taken from `fit` if
#'   absent.
#' @param associations A precomputed table of factor-covariate associations.
#'
#' @returns A data frame with one row per programme carrying `joint_p`, the
#'   `best_single_p` among its members, the `margin` between them on the
#'   negative log scale, the worst leave-one-out p-value, and `needs_multiple`,
#'   true when the joint evidence beats every single modality and no single
#'   removal leaves the programme as strong.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 3, n_features = 150,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 4, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
#' chorale_added_value(fit)
chorale_added_value <- function(fit, programmes = NULL, associations = NULL) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  if (is.null(programmes)) programmes <- chorale_programmes(fit)
  if (nrow(programmes) == 0) return(data.frame())

  assoc <- associations %||% try(chorale_association_table(fit), silent = TRUE)
  single_p <- function(mod, fac) {
    if (inherits(assoc, "try-error") || !is.data.frame(assoc) || nrow(assoc) == 0) {
      return(NA_real_)
    }
    r <- assoc[assoc$modality == mod & assoc$factor == fac &
                 assoc$covariate == "phenotype", , drop = FALSE]
    if (nrow(r) == 0) return(NA_real_)
    as.numeric(r$p_permutation[1])
  }

  loo <- fit$leave_one_out
  rows <- lapply(unique(programmes$programme), function(pr) {
    d <- programmes[programmes$programme == pr, , drop = FALSE]
    singles <- vapply(seq_len(nrow(d)), function(i)
      single_p(d$modality[i], d$factor[i]), numeric(1))
    best <- suppressWarnings(min(singles, na.rm = TRUE))
    if (!is.finite(best)) best <- NA_real_
    joint <- d$joint_p[1]
    l <- if (!is.null(loo) && nrow(loo) > 0) {
      loo[loo$programme == pr, , drop = FALSE]
    } else {
      NULL
    }
    worst_loo <- if (!is.null(l) && nrow(l) > 0) max(l$joint_p, na.rm = TRUE) else NA_real_
    # A permutation p-value saturates at its floor, so two results at the floor
    # cannot be ordered by it. Whether a modality was needed is therefore read
    # from the change in the statistic when it is removed, which does not
    # saturate; a programme needs its modalities when removing any one of them
    # weakens it.
    worst_delta <- if (!is.null(l) && nrow(l) > 0) max(l$delta, na.rm = TRUE) else NA_real_
    data.frame(
      programme = pr,
      n_modalities = d$n_modalities[1],
      joint_p = joint,
      best_single_p = best,
      # On the negative log scale, so a positive margin means the programme is
      # stronger than anything one modality carried.
      margin = round(-log10(joint) + log10(best), 3),
      worst_leave_one_out_p = worst_loo,
      worst_leave_one_out_delta = worst_delta,
      needs_multiple = isTRUE(joint <= best) &&
        (is.na(worst_delta) || worst_delta < 0),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$margin), , drop = FALSE]
}

#' Whether a programme is specific to the disease phenotype
#'
#' A programme that answers to any strong contrast is not a statement about the
#' disease. Each nuisance covariate is put in the phenotype's place and the
#' whole procedure rerun, so a programme can be compared against what the same
#' pipeline recovers when a covariate that is not the disease drives the
#' matching. A programme whose evidence under the phenotype is no stronger than
#' under a nuisance anchor is reporting cohort structure.
#'
#' @param fit A `chorale_fit` object.
#' @param containers The modality containers the fit was built from.
#' @param covariates Nuisance covariates to substitute for the phenotype.
#' @param n_init Initialisations per refit.
#' @param seed Integer seed.
#'
#' @returns A data frame with one row per substituted covariate, carrying the
#'   strongest joint evidence the pipeline reaches with that covariate in the
#'   phenotype's place, beside the observed value under the true phenotype.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' chorale_specificity(fit, containers, covariates = "sex", n_init = 2)
chorale_specificity <- function(fit, containers, covariates = c("sex", "batch"),
                                n_init = 5L, seed = 1L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  observed <- chorale_best_joint(fit)
  rows <- list()
  for (cv in covariates) {
    swapped <- containers
    usable <- TRUE
    for (m in names(swapped)) {
      d <- as.data.frame(SummarizedExperiment::colData(swapped[[m]]))
      if (!cv %in% colnames(d) ||
          length(unique(stats::na.omit(d[[cv]]))) < 2) {
        usable <- FALSE
        break
      }
      # The nuisance covariate takes the phenotype's place, so the matching is
      # anchored on it and on nothing else that distinguishes the design.
      d$phenotype <- as.character(d[[cv]])
      SummarizedExperiment::colData(swapped[[m]]) <- S4Vectors::DataFrame(d)
    }
    if (!usable) {
      rows[[length(rows) + 1L]] <- data.frame(
        anchor = cv, joint_statistic = NA_real_, observed_phenotype = observed,
        stronger_than_phenotype = NA, reason = "covariate absent or constant",
        stringsAsFactors = FALSE
      )
      next
    }
    refit <- try(
      chorale_fit(swapped, n_factors = fit$n_factors, n_init = n_init,
                  strata_keys = fit$strata_keys, n_pathway_perm = 0L,
                  seed = seed),
      silent = TRUE
    )
    value <- if (inherits(refit, "try-error")) NA_real_ else chorale_best_joint(refit)
    rows[[length(rows) + 1L]] <- data.frame(
      anchor = cv,
      joint_statistic = round(value, 4),
      observed_phenotype = round(observed, 4),
      stronger_than_phenotype = isTRUE(value >= observed),
      reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
