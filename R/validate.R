#' Align recovered factors to the planted truth within one modality
#'
#' Independent component analysis returns factors in an arbitrary order and
#' sign, so a recovered factor is identified with a planted one by the absolute
#' correlation of their scores. Each recovered factor is labelled with the
#' planted factor it matches best, `shared_k` or `private_j`, which is what lets
#' a recovered assignment be scored against the answer.
#'
#' @param fit A `chorale_fit` object.
#' @param sim The `chorale_simulate()` output the fit was built from.
#'
#' @returns A data frame with one row per recovered factor, carrying the
#'   modality, the recovered factor name, the planted label it matches, whether
#'   that label is shared, and the correlation.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 3)
#' chorale_align_truth(fit, sim)
chorale_align_truth <- function(fit, sim) {
  rows <- list()
  for (mi in seq_along(fit$modalities)) {
    m <- fit$modalities[mi]
    recovered <- fit$fits[[m]]$scores
    planted <- cbind(sim$truth$scores[[mi]]$shared,
                     sim$truth$scores[[mi]]$private)
    labels <- c(paste0("shared_", seq_len(ncol(sim$truth$scores[[mi]]$shared))),
                paste0("private_", seq_len(ncol(sim$truth$scores[[mi]]$private))))
    cormat <- abs(suppressWarnings(stats::cor(recovered, planted)))
    cormat[!is.finite(cormat)] <- 0
    for (f in seq_len(ncol(recovered))) {
      best <- which.max(cormat[f, ])
      rows[[length(rows) + 1L]] <- data.frame(
        modality = m,
        factor = colnames(recovered)[f],
        planted = labels[best],
        shared = grepl("^shared_", labels[best]),
        correlation = round(cormat[f, best], 3),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Score one fit against the structure that was planted
#'
#' Reduces a fit on simulated data to the quantities a validation cares about:
#' whether the programmes it recovered correspond to the shared factors that
#' were planted, and whether the assignments it made are between factors that
#' are genuinely the same planted factor.
#'
#' A programme is counted correct when its members all align to one shared
#' planted factor. A recovered assignment between two factors is a true match
#' when both align to the same shared factor, and a false match otherwise;
#' the false-match rate among supported assignments is the quantity the
#' modality-count argument is about.
#'
#' @param fit A `chorale_fit` object.
#' @param sim The `chorale_simulate()` output the fit was built from.
#'
#' @returns A one-row data frame with `programmes_supported`,
#'   `programmes_correct`, `assignment_accuracy`, `false_match_rate` and
#'   `shared_recovered`, the fraction of planted shared factors recovered as a
#'   supported programme spanning at least two modalities.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 3, n_features = 150,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 4, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 3)
#' chorale_score_recovery(fit, sim)
chorale_score_recovery <- function(fit, sim) {
  align <- chorale_align_truth(fit, sim)
  key <- stats::setNames(align$planted, paste(align$modality, align$factor))
  is_shared <- stats::setNames(align$shared, paste(align$modality, align$factor))

  m <- fit$matches
  supported <- if (is.data.frame(m) && nrow(m) > 0) m[m$significant, , drop = FALSE] else m
  n_assign <- if (is.data.frame(supported)) nrow(supported) else 0L
  true_match <- 0L
  if (n_assign > 0) {
    a <- key[paste(supported$modality_a, supported$factor_a)]
    b <- key[paste(supported$modality_b, supported$factor_b)]
    sa <- is_shared[paste(supported$modality_a, supported$factor_a)]
    true_match <- sum(a == b & sa, na.rm = TRUE)
  }

  pg <- chorale_programmes(fit, significant_only = TRUE)
  n_prog <- if (nrow(pg) > 0) length(unique(pg$programme)) else 0L
  correct <- 0L
  recovered_shared <- character(0)
  if (n_prog > 0) {
    for (pr in unique(pg$programme)) {
      d <- pg[pg$programme == pr, , drop = FALSE]
      planted <- key[paste(d$modality, d$factor)]
      shared <- is_shared[paste(d$modality, d$factor)]
      if (nrow(d) >= 2 && all(shared) && length(unique(planted)) == 1) {
        correct <- correct + 1L
        recovered_shared <- c(recovered_shared, unique(planted))
      }
    }
  }

  data.frame(
    programmes_supported = n_prog,
    programmes_correct = correct,
    assignment_accuracy = if (n_assign > 0) round(true_match / n_assign, 3) else NA_real_,
    false_match_rate = if (n_assign > 0) round(1 - true_match / n_assign, 3) else NA_real_,
    shared_recovered = round(length(unique(recovered_shared)) /
                               sim$truth$n_shared_factors, 3),
    stringsAsFactors = FALSE
  )
}

#' Run the estimator across a grid of simulated regimes
#'
#' The validation matrix. Each row of `grid` is a regime the estimator is run
#' on, and the recovery, false-match rate and, where relevant, interval
#' coverage are measured against the planted truth. A regime that satisfies the
#' identification conditions should recover its shared programmes; one that
#' violates them, or that plants distinct programmes sharing a phenotype
#' response, is where the estimator's limits are read.
#'
#' @param grid A data frame of regimes. Recognised columns are `label`,
#'   `n_modalities`, `n_features`, `n_shared_factors`, `n_private_factors`,
#'   `n_strains`, `n_per_cell`, `effect_size`, `imbalance` and `n_init`; any
#'   absent column takes a default. A `same_response` column, when `TRUE`,
#'   plants two shared factors with the same phenotype response, the adversarial
#'   case in which distinct programmes should not be merged.
#' @param n_rep Replicates per regime, over which the metrics are averaged.
#' @param seed Integer seed; replicate `r` of regime `i` uses `seed + 100 * i + r`.
#'
#' @returns `grid` with the mean recovery metrics joined on.
#' @export
#' @examples
#' grid <- data.frame(label = c("clean", "null"), effect_size = c(3, 0))
#' chorale_validate(grid, n_rep = 1)
chorale_validate <- function(grid, n_rep = 3L, seed = 1L) {
  default <- list(n_modalities = 3L, n_features = 150L, n_shared_factors = 2L,
                  n_private_factors = 1L, n_strains = 4L, n_per_cell = 4L,
                  effect_size = 3, imbalance = 0, n_init = 5L,
                  same_response = FALSE)
  get_col <- function(row, name) if (name %in% names(grid)) grid[[name]][row] else default[[name]]

  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    metrics <- vector("list", n_rep)
    for (r in seq_len(n_rep)) {
      ns <- as.integer(get_col(i, "n_shared_factors"))
      sig <- NULL
      if (isTRUE(get_col(i, "same_response"))) {
        # Two shared factors responding to the phenotype the same way: distinct
        # programmes the estimator must keep apart rather than merge.
        sig <- matrix(0, nrow = ns, ncol = 3,
                      dimnames = list(NULL, c("phenotype", "age", "sex")))
        sig[, 1] <- 1
        if (ns >= 2) sig[2, ] <- c(1, 0, 0)
      }
      sim <- chorale_simulate(
        n_modalities = as.integer(get_col(i, "n_modalities")),
        n_features = as.integer(get_col(i, "n_features")),
        n_shared_factors = ns,
        n_private_factors = as.integer(get_col(i, "n_private_factors")),
        n_strains = as.integer(get_col(i, "n_strains")),
        n_per_cell = as.integer(get_col(i, "n_per_cell")),
        effect_size = get_col(i, "effect_size"),
        imbalance = get_col(i, "imbalance"),
        signature = sig,
        seed = seed + 100L * i + r
      )
      containers <- Map(chorale_load, sim$modalities, sim$col_data)
      nf <- ns + as.integer(get_col(i, "n_private_factors"))
      fit <- try(chorale_fit(containers, n_factors = nf,
                             n_init = as.integer(get_col(i, "n_init")),
                             seed = seed + 100L * i + r), silent = TRUE)
      if (inherits(fit, "try-error")) next
      metrics[[r]] <- chorale_score_recovery(fit, sim)
    }
    metrics <- do.call(rbind, metrics)
    agg <- if (is.null(metrics)) {
      data.frame(programmes_supported = NA_real_, programmes_correct = NA_real_,
                 assignment_accuracy = NA_real_, false_match_rate = NA_real_,
                 shared_recovered = NA_real_)
    } else {
      data.frame(
        programmes_supported = mean(metrics$programmes_supported, na.rm = TRUE),
        programmes_correct = mean(metrics$programmes_correct, na.rm = TRUE),
        assignment_accuracy = round(mean(metrics$assignment_accuracy, na.rm = TRUE), 3),
        false_match_rate = round(mean(metrics$false_match_rate, na.rm = TRUE), 3),
        shared_recovered = round(mean(metrics$shared_recovered, na.rm = TRUE), 3)
      )
    }
    out[[i]] <- cbind(grid[i, , drop = FALSE], agg, n_rep = n_rep)
  }
  do.call(rbind, out)
}

#' Calibration of the joint null under no shared structure
#'
#' The p-value of a genuine null must be uniform, or the reported significance
#' does not mean what it says. Data with no design-linked shared structure are
#' simulated repeatedly, the estimator is run in full on each, and the joint
#' p-value of its strongest programme is collected. A Kolmogorov-Smirnov test of
#' those p-values against the uniform distribution is the calibration: a small
#' p-value there is evidence the procedure is miscalibrated, and the false
#' positive rate at a threshold reports how often a null run would be called a
#' discovery.
#'
#' @param n_sim Number of null datasets.
#' @param alpha Threshold at which the false positive rate is reported.
#' @param n_perm Permutations calibrating each fit.
#' @param n_init Initialisations per fit.
#' @param ... Passed to [chorale_simulate()], except `effect_size`, which is
#'   fixed at zero.
#' @param seed Integer seed.
#'
#' @returns A list with the collected `p_values`, the `ks_p` of their
#'   uniformity test, and the `false_positive_rate` at `alpha`.
#' @export
#' @examples
#' chorale_null_calibration(n_sim = 10, n_perm = 99, n_init = 3,
#'                          n_features = 120, n_per_cell = 3)
chorale_null_calibration <- function(n_sim = 100L, alpha = 0.05,
                                     n_perm = 200L, n_init = 5L, ...,
                                     seed = 1L) {
  pvals <- rep(NA_real_, n_sim)
  for (s in seq_len(n_sim)) {
    sim <- chorale_simulate(effect_size = 0, ..., seed = seed + s)
    containers <- Map(chorale_load, sim$modalities, sim$col_data)
    nf <- sim$truth$n_shared_factors + sim$truth$n_private_factors
    fit <- try(chorale_fit(containers, n_factors = nf, n_init = n_init,
                           n_pathway_perm = 0L, seed = seed + s), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit$programmes) ||
        nrow(fit$programmes) == 0) next
    pvals[s] <- min(fit$programmes$joint_p, na.rm = TRUE)
  }
  p <- pvals[is.finite(pvals)]
  ks <- if (length(p) >= 5) {
    suppressWarnings(stats::ks.test(p, "punif")$p.value)
  } else {
    NA_real_
  }
  list(
    p_values = p,
    ks_p = ks,
    false_positive_rate = if (length(p) > 0) mean(p <= alpha) else NA_real_,
    n_evaluated = length(p)
  )
}
