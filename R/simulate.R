#' Simulate disjoint multi-modal data with known shared latent structure
#'
#' Generates synthetic data for `n_modalities` modalities measured on disjoint
#' sample sets. Each sample's features are a linear mixture of latent factors
#' plus additive Gaussian measurement noise, and the latent factors are drawn
#' so that the conditions the identification results rest on hold by
#' construction: every factor is non-Gaussian, no two factors share a
#' distribution, and each factor's distribution is non-symmetric (Sturma et al.
#' 2023). A validation of the estimator needs the data to satisfy the
#' assumptions it is judged against, which the earlier symmetric,
#' identically-distributed design did not.
#'
#' A shared factor is what the estimator should recover as one cross-modality
#' programme. Its correspondence across modalities is carried by its
#' **design signature**, one weight per design contrast, which is the same in
#' every modality; its remaining, design-independent variation is drawn
#' independently per modality, since the animals are disjoint. Matching two
#' modalities on their design response is therefore what the truth rewards, and
#' the signature is the answer a recovered assignment is scored against. Each
#' shared factor is given a distinct dominant contrast, so the factors' design
#' responses are uncorrelated and the recovered components stay separable.
#' Private factors carry no signature and no cross-modality partner.
#'
#' @param n_modalities Integer number of modalities to simulate.
#' @param n_features Integer vector of length `n_modalities`, features per
#'   modality. Recycled if length 1.
#' @param n_shared_factors Integer number of shared latent factors.
#' @param n_private_factors Integer number of modality-private latent factors,
#'   common to every modality.
#' @param n_strains Integer number of strains in the genetic panel.
#' @param n_per_cell Integer number of animals per design cell, per modality.
#'   Samples are drawn independently per modality, so no animal appears in more
#'   than one modality.
#' @param noise_sd Standard deviation of the additive Gaussian measurement
#'   noise.
#' @param effect_size Magnitude of the design signatures carried by the shared
#'   factors, in units of the latent scale. Zero detaches the shared factors
#'   from the design, which is the complete null the matching is calibrated
#'   against.
#' @param signature Optional `n_shared_factors` by three matrix of design
#'   weights over the phenotype, age and sex contrasts, shared across
#'   modalities. When absent each shared factor is given a distinct dominant
#'   contrast. Supplying a signature in which two factors respond to the
#'   phenotype the same way is how a same-response, distinct-programme
#'   adversarial case is built.
#' @param confounder Optional named list adding a nuisance covariate correlated
#'   with the phenotype, with elements `name`, `rho` (the correlation) and
#'   `loading` (its weight on the shared factors), for testing robustness to
#'   confounding.
#' @param imbalance Optional numeric in `[0, 1)`. When positive, design cells
#'   are populated unevenly across modalities, for testing behaviour under the
#'   unequal mixtures disjoint cohorts produce.
#' @param seed Integer random seed.
#'
#' @returns A list with `modalities` (feature-by-sample matrices),
#'   `col_data` (per-modality design tables) and `truth`, which carries the
#'   shared and private loadings and scores per modality, the marker indices,
#'   the design `signature`, and the factor counts, for scoring recovery.
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
                             effect_size = 1,
                             signature = NULL,
                             confounder = NULL,
                             imbalance = 0,
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
  if (imbalance < 0 || imbalance >= 1) {
    rlang::abort("`imbalance` must lie in [0, 1).")
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
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

  # Each shared factor is given a distinct dominant contrast, so their design
  # responses are uncorrelated and the recovered components remain separable.
  # The first factor is phenotype-dominant, the factor the case/control anchor
  # should recover.
  if (is.null(signature)) {
    signature <- matrix(0, nrow = n_shared_factors, ncol = 3,
                        dimnames = list(NULL, c("phenotype", "age", "sex")))
    for (k in seq_len(n_shared_factors)) {
      dom <- ((k - 1L) %% 3L) + 1L
      signature[k, dom] <- 1
    }
  }
  signature <- signature * effect_size

  n_total_factors <- n_shared_factors + n_private_factors

  modalities <- vector("list", n_modalities)
  col_data <- vector("list", n_modalities)
  truth_loadings <- vector("list", n_modalities)
  truth_scores <- vector("list", n_modalities)
  truth_markers <- vector("list", n_modalities)
  names(modalities) <- names(col_data) <- names(truth_loadings) <-
    names(truth_scores) <- names(truth_markers) <-
    paste0("modality_", seq_len(n_modalities))

  # A latent source that is non-Gaussian and non-symmetric, with a shape that
  # sets its skewness. Distinct shapes across factors and modalities make the
  # sources pairwise different, standardised to zero mean and unit variance so
  # scale carries nothing.
  skewed <- function(n, shape, sign = 1) {
    v <- stats::rgamma(n, shape = shape)
    v <- (v - shape) / sqrt(shape)
    sign * v
  }
  # Shapes are held below four so every source stays clearly non-Gaussian;
  # a gamma of large shape approaches the normal and cannot be recovered.
  shape_grid <- seq(0.8, 3.5, length.out = max(n_total_factors * n_modalities, 2))
  shape_at <- function(m, k) shape_grid[((m - 1L) * n_total_factors + k - 1L) %% length(shape_grid) + 1L]

  for (m in seq_len(n_modalities)) {
    cells <- design[rep(seq_len(nrow(design)), each = n_per_cell), ]
    if (imbalance > 0) {
      # Thin each cell independently per modality, so the modalities realise the
      # design in different proportions.
      keep <- stats::runif(nrow(cells)) > imbalance * stats::runif(1)
      if (sum(keep) >= 2 * n_total_factors) cells <- cells[keep, , drop = FALSE]
    }
    rownames(cells) <- NULL
    n_m <- nrow(cells)
    p_m <- n_features[m]

    sample_id <- sprintf("modality%d_sample%04d", m, seq_len(n_m))

    contrasts <- cbind(
      phenotype = ifelse(cells$phenotype == "5XFAD", 0.5, -0.5),
      age = ifelse(cells$age_months == 14, 0.5, -0.5),
      sex = ifelse(cells$sex == "M", 0.5, -0.5)
    )

    # Shared factor k: an independent non-symmetric base, plus the design
    # response its shared signature prescribes. The signature is the same in
    # every modality, which is the correspondence the estimator must recover.
    shared_scores <- matrix(0, nrow = n_m, ncol = n_shared_factors)
    for (k in seq_len(n_shared_factors)) {
      base <- skewed(n_m, shape_at(m, k), sign = if (k %% 2 == 0) -1 else 1)
      shared_scores[, k] <- base + as.numeric(contrasts %*% signature[k, ])
    }

    # Private factor: non-symmetric, distinct shape, no design response and no
    # cross-modality partner.
    private_scores <- matrix(0, nrow = n_m, ncol = n_private_factors)
    for (j in seq_len(n_private_factors)) {
      private_scores[, j] <- skewed(n_m, shape_at(m, n_shared_factors + j),
                                    sign = if (j %% 2 == 0) 1 else -1)
    }

    batch <- rep("batch1", n_m)
    if (!is.null(confounder)) {
      # A nuisance covariate correlated with the phenotype, loading on the
      # shared factors, so a matched programme can be probed for whether it
      # tracks the phenotype or the confounder.
      z <- ifelse(cells$phenotype == "5XFAD", 1, -1) * confounder$rho +
        stats::rnorm(n_m) * sqrt(1 - confounder$rho^2)
      shared_scores <- shared_scores + outer(z, rep(confounder$loading, n_shared_factors))
      batch <- ifelse(z > stats::median(z), "batchB", "batchA")
    }

    scores <- cbind(shared_scores, private_scores)

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
    x <- t(scores %*% t(loadings) + noise)
    dimnames(x) <- list(sprintf("modality%d_feature%05d", m, seq_len(p_m)), sample_id)

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
      batch = batch,
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
      signature = signature,
      n_shared_factors = n_shared_factors,
      n_private_factors = n_private_factors
    )
  )
}
