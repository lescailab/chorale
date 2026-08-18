#' Simulate disjoint multi-modal data with known shared latent structure
#'
#' Generates synthetic data for `n_modalities` modalities measured on disjoint
#' sample sets, following the model of `MATHEMATICAL_FOUNDATION.md` Section 8:
#' each sample's features are a linear mixture of a shared latent state and a
#' modality-private latent state, both non-Gaussian, plus additive Gaussian
#' measurement noise. Samples carry a genetic-panel design (strain, phenotype,
#' age, sex) that supplies the anchor structure of Section 6, and every
#' shared factor is given at least two pure marker features per modality
#' (Section 6, D5), so recovery can be checked against a known answer.
#'
#' The shared latent components are drawn from a Laplace distribution and the
#' modality-private components from a distinct, non-symmetric distribution
#' per modality (shifted chi-squared with a modality-specific degrees of
#' freedom), so that the non-Gaussianity (D7) and pairwise-difference (D4)
#' identification conditions hold by construction.
#'
#' @param n_modalities Integer number of modalities to simulate.
#' @param n_features Integer vector of length `n_modalities`, features per
#'   modality. Recycled if length 1.
#' @param n_shared_factors Integer number of shared latent factors.
#' @param n_private_factors Integer number of modality-private latent
#'   factors, common to every modality.
#' @param n_strains Integer number of strains in the genetic panel.
#' @param n_per_cell Integer number of animals per strain by phenotype by age
#'   by sex cell, per modality. Samples are drawn independently per modality,
#'   so no animal appears in more than one modality.
#' @param noise_sd Standard deviation of the additive Gaussian measurement
#'   noise.
#' @param seed Integer random seed, for reproducibility.
#'
#' @returns A list with components:
#'   \describe{
#'     \item{modalities}{Named list of feature-by-sample matrices, one per
#'       modality.}
#'     \item{col_data}{Named list of per-modality sample metadata data
#'       frames, with columns `sample_id`, `cohort`, `modality`, `strain`,
#'       `phenotype`, `age_months`, `sex`, `region`, `batch`.}
#'     \item{truth}{List with the shared and private loading matrices per
#'       modality, the shared and private factor scores per modality, and the
#'       indices of the pure marker features per shared factor per modality.}
#'   }
#' @export
#' @examples
#' sim <- chorale_simulate(
#'   n_modalities = 2, n_features = 40, n_shared_factors = 3,
#'   n_private_factors = 2, n_strains = 4, n_per_cell = 2, seed = 1
#' )
#' names(sim$modalities)
#' dim(sim$modalities[[1]])
chorale_simulate <- function(n_modalities = 3,
                              n_features = 500,
                              n_shared_factors = 5,
                              n_private_factors = 3,
                              n_strains = 8,
                              n_per_cell = 3,
                              noise_sd = 0.1,
                              seed = 1) {
  if (n_modalities < 2) {
    rlang::abort("`n_modalities` must be at least 2.")
  }
  if (length(n_features) == 1) {
    n_features <- rep(n_features, n_modalities)
  }
  if (length(n_features) != n_modalities) {
    rlang::abort("`n_features` must have length 1 or `n_modalities`.")
  }
  min_features_needed <- 2L * n_shared_factors + n_private_factors
  if (any(n_features < min_features_needed)) {
    rlang::abort(paste(
      "Each modality needs at least", min_features_needed, "features:",
      "two pure markers per shared factor plus the private factors."
    ))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  })
  set.seed(seed)

  design <- expand.grid(
    strain = paste0("BXD", seq_len(n_strains)),
    phenotype = c("Ntg", "5XFAD"),
    age_months = c(6, 14),
    sex = c("F", "M"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  n_total_factors <- n_shared_factors + n_private_factors

  # Strain-level random effect on the shared factors, so strain is
  # informative of the shared state (anchor structure, Section 6).
  strain_effect <- matrix(
    stats::rnorm(n_strains * n_shared_factors, sd = 0.8),
    nrow = n_strains, ncol = n_shared_factors,
    dimnames = list(unique(design$strain), NULL)
  )

  modalities <- vector("list", n_modalities)
  col_data <- vector("list", n_modalities)
  truth_loadings <- vector("list", n_modalities)
  truth_scores <- vector("list", n_modalities)
  truth_markers <- vector("list", n_modalities)
  names(modalities) <- names(col_data) <- names(truth_loadings) <-
    names(truth_scores) <- names(truth_markers) <- paste0("modality_", seq_len(n_modalities))

  for (m in seq_len(n_modalities)) {
    cells <- design[rep(seq_len(nrow(design)), each = n_per_cell), ]
    rownames(cells) <- NULL
    n_m <- nrow(cells)
    p_m <- n_features[m]

    sample_id <- sprintf("modality%d_sample%04d", m, seq_len(n_m))

    # Shared latent scores: Laplace (non-Gaussian), shifted by the
    # sample's strain effect.
    shared_scores <- matrix(
      stats::rexp(n_m * n_shared_factors) - stats::rexp(n_m * n_shared_factors),
      nrow = n_m, ncol = n_shared_factors
    )
    shared_scores <- shared_scores + strain_effect[cells$strain, , drop = FALSE]

    # Private latent scores: shifted chi-squared, degrees of freedom set by
    # modality index, so the private-component distribution is
    # non-symmetric and pairwise different across modalities (D4).
    df_m <- 2 + m
    private_scores <- matrix(
      stats::rchisq(n_m * n_private_factors, df = df_m) - df_m,
      nrow = n_m, ncol = n_private_factors
    )

    scores <- cbind(shared_scores, private_scores)

    # Mixing matrix: random full column-rank loadings, with two pure
    # marker features per shared factor (Section 6, D5).
    loadings <- matrix(
      stats::rnorm(p_m * n_total_factors, sd = 1),
      nrow = p_m, ncol = n_total_factors
    )
    markers <- vector("list", n_shared_factors)
    for (k in seq_len(n_shared_factors)) {
      marker_idx <- c(2L * k - 1L, 2L * k)
      loadings[marker_idx, -k] <- 0
      loadings[marker_idx, k] <- stats::runif(2, min = 0.5, max = 1.5)
      markers[[k]] <- marker_idx
    }
    names(markers) <- paste0("factor_", seq_len(n_shared_factors))

    noise <- matrix(stats::rnorm(n_m * p_m, sd = noise_sd), nrow = n_m, ncol = p_m)
    x <- scores %*% t(loadings) + noise

    feature_id <- sprintf("modality%d_feature%05d", m, seq_len(p_m))
    x <- t(x)
    dimnames(x) <- list(feature_id, sample_id)

    modalities[[m]] <- x
    col_data[[m]] <- data.frame(
      sample_id = sample_id,
      cohort = "simulated",
      modality = names(modalities)[m],
      strain = cells$strain,
      phenotype = cells$phenotype,
      age_months = cells$age_months,
      sex = cells$sex,
      region = "simulated",
      batch = "batch1",
      stringsAsFactors = FALSE
    )
    truth_loadings[[m]] <- list(
      shared = loadings[, seq_len(n_shared_factors), drop = FALSE],
      private = loadings[, n_shared_factors + seq_len(n_private_factors), drop = FALSE]
    )
    truth_scores[[m]] <- list(shared = shared_scores, private = private_scores)
    truth_markers[[m]] <- markers
  }

  list(
    modalities = modalities,
    col_data = col_data,
    truth = list(
      loadings = truth_loadings,
      scores = truth_scores,
      markers = truth_markers,
      n_shared_factors = n_shared_factors,
      n_private_factors = n_private_factors
    )
  )
}
