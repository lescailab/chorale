#' Recover a correspondence that is known, with the pairing withheld
#'
#' Every other check asks whether the estimator's output beats a null. This asks
#' the harder question: on data where the true cross-modality correspondence is
#' known, does the estimator find it? Two modalities measured on the *same*
#' individuals carry that answer. The pairing is discarded, the estimator is run
#' as though the samples were disjoint, and the recovered correspondence is then
#' compared with the pairing that was withheld.
#'
#' The truth is a correspondence between factors, not between samples: a factor
#' recovered in one modality and a factor recovered in the other are the same
#' programme when their scores agree across the animals both were measured on.
#' That agreement is computable only because the pairing exists, and it is
#' withheld from the estimator throughout.
#'
#' Three quantities frame the result, in the manner the plan sets out. The
#' **paired benchmark** is the correspondence obtained by using the pairing,
#' which is the best any method could do. **Random alignment** is the
#' correspondence obtained by pairing factors at random, which is the worst.
#' The estimator's recovery is reported between them, and a result near the
#' lower bound is a failure however small its p-value.
#'
#' @param paired_a,paired_b Feature-by-sample matrices for two modalities
#'   measured on the same individuals, with the same sample identifiers in the
#'   same order.
#' @param design A design table for those individuals, carrying `sample_id` and
#'   the covariates the estimator anchors on.
#' @param n_factors Factors per modality, or `"auto"`.
#' @param n_init Initialisations per fit.
#' @param n_random Random alignments forming the lower bound.
#' @param seed Integer seed.
#'
#' @returns A list with `truth`, the factor correspondence the pairing implies;
#'   `recovered`, the correspondence the estimator found without it; and
#'   `summary`, a one-row data frame carrying the recovered agreement against
#'   the paired upper bound and the random lower bound, and the fraction of
#'   factors whose partner was recovered correctly.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' # Give the two modalities the same animals, so a pairing exists to withhold.
#' a <- sim$modalities[[1]]
#' b <- sim$modalities[[2]][, seq_len(ncol(a))]
#' colnames(b) <- colnames(a)
#' chorale_destroy_pairing(a, b, sim$col_data[[1]], n_factors = 3, n_init = 2,
#'                         n_random = 20)$summary
chorale_destroy_pairing <- function(paired_a, paired_b, design,
                                    n_factors = "auto", n_init = 5L,
                                    n_random = 200L, seed = 1L) {
  shared <- intersect(colnames(paired_a), colnames(paired_b))
  shared <- intersect(shared, as.character(design$sample_id))
  if (length(shared) < 20) {
    rlang::abort(paste0(
      "The two modalities share ", length(shared), " samples. A pairing that ",
      "can be withheld needs at least twenty."
    ))
  }
  a <- paired_a[, shared, drop = FALSE]
  b <- paired_b[, shared, drop = FALSE]
  d <- design[match(shared, design$sample_id), , drop = FALSE]

  # The estimator never sees that the samples correspond: each modality is given
  # its own identifiers, as a disjoint pair of cohorts would be.
  da <- d
  db <- d
  da$sample_id <- paste0("a_", seq_along(shared))
  db$sample_id <- paste0("b_", seq_along(shared))
  colnames(a) <- da$sample_id
  colnames(b) <- db$sample_id
  containers <- list(A = chorale_load(a, da), B = chorale_load(b, db))

  fit <- chorale_fit(containers, n_factors = n_factors, n_init = n_init,
                     n_pathway_perm = 0L, seed = seed)

  # The withheld truth: factors correspond when their scores agree across the
  # animals both were measured on.
  sa <- fit$fits$A$scores
  sb <- fit$fits$B$scores
  agreement <- abs(suppressWarnings(stats::cor(sa, sb)))
  agreement[!is.finite(agreement)] <- 0
  truth_assign <- chorale_assign(agreement)
  truth <- data.frame(
    factor_a = colnames(sa),
    factor_b = colnames(sb)[truth_assign],
    agreement = round(agreement[cbind(seq_len(nrow(agreement)),
                                      truth_assign)], 4),
    stringsAsFactors = FALSE
  )
  paired_bound <- mean(truth$agreement, na.rm = TRUE)

  # What the phenotype-led estimator supported, knowing nothing of the pairing.
  m <- fit$matches
  recovered <- if (is.data.frame(m) && nrow(m) > 0) {
    eligible <- m$supported & m$resolution_status %in% c("resolved", "ambiguous")
    m <- m[eligible, , drop = FALSE]
    data.frame(factor_a = m$factor_a, factor_b = m$factor_b,
               significant = m$significant, stringsAsFactors = FALSE)
  } else {
    data.frame(factor_a = character(0), factor_b = character(0),
               significant = logical(0), stringsAsFactors = FALSE)
  }
  recovered_agreement <- if (nrow(recovered) > 0) {
    mean(vapply(seq_len(nrow(recovered)), function(i) {
      agreement[recovered$factor_a[i], recovered$factor_b[i]]
    }, numeric(1)), na.rm = TRUE)
  } else {
    NA_real_
  }
  correct <- if (nrow(recovered) > 0) {
    key <- stats::setNames(truth$factor_b, truth$factor_a)
    mean(recovered$factor_b == key[recovered$factor_a], na.rm = TRUE)
  } else {
    NA_real_
  }

  # The lower bound: what pairing factors at random achieves on the same data.
  set.seed(seed)
  random <- vapply(seq_len(n_random), function(i) {
    j <- sample(seq_len(ncol(agreement)), nrow(agreement), replace = TRUE)
    mean(agreement[cbind(seq_len(nrow(agreement)), j)])
  }, numeric(1))
  random_bound <- mean(random)

  # Where the recovery sits between the bounds: one is the paired fit, zero is
  # random alignment.
  placement <- if (is.finite(recovered_agreement) &&
                   paired_bound > random_bound) {
    (recovered_agreement - random_bound) / (paired_bound - random_bound)
  } else {
    NA_real_
  }

  summary <- data.frame(
    n_samples = length(shared),
    n_factors_a = ncol(sa),
    n_factors_b = ncol(sb),
    paired_upper_bound = round(paired_bound, 4),
    recovered_agreement = round(recovered_agreement, 4),
    random_lower_bound = round(random_bound, 4),
    placement_between_bounds = round(placement, 3),
    fraction_partner_correct = round(correct, 3),
    n_recovered_pairs = nrow(recovered),
    stringsAsFactors = FALSE
  )
  list(truth = truth, recovered = recovered, summary = summary, fit = fit)
}
