#' Simulate disjoint multi-modal data with known shared latent structure
#'
#' Generates synthetic data for `n_modalities` modalities measured on disjoint
#' sample sets. Each sample's features are a linear mixture of latent factors
#' plus additive Gaussian measurement noise. Latent factors are non-Gaussian,
#' have different shapes and are non-symmetric (Sturma et al. 2023), which
#' provides a demanding ICA simulation. These choices are not presented as
#' identification conditions for CHORALE's phenotype correspondence.
#'
#' A shared factor is what a component-recovery benchmark should recover across
#' modalities. Its correspondence across modalities is carried by its
#' **design signature**, one weight per design contrast, which is the same in
#' every modality; its remaining, design-independent variation is drawn
#' independently per modality, since the individuals are disjoint. Matching two
#' modalities on their design response is therefore what the truth rewards, and
#' the signature is the answer a recovered assignment is scored against. Each
#' shared factor is given a distinct dominant contrast, so the factors' design
#' responses are uncorrelated and the recovered components stay separable.
#' Private factors carry no signature and no cross-modality partner.
#'
#' Supplying `profile` replaces the parametric observation model with the
#' marginals of real data. The latent structure is unchanged, and each simulated
#' feature is mapped onto the empirical distribution of a real feature through
#' its quantile grid, carrying that feature's missingness with it. The design is
#' drawn from the cells the real cohort populated, so a cell no sample occupied
#' stays empty. Realism is therefore inherited from a
#' [chorale_data_profile()] rather than asserted by a distributional choice, and
#' the map is monotone per feature, so the sign and the ordering of a planted
#' design response survive it.
#'
#' @param n_modalities Integer number of modalities to simulate.
#' @param n_features Integer vector of length `n_modalities`, features per
#'   modality. Recycled if length 1. `NULL` is allowed only with `profile`, and
#'   takes every feature of each profile, keeping its identifier, which is what
#'   planting a named pathway requires.
#' @param n_shared_factors Integer number of shared latent factors.
#' @param n_private_factors Integer number of modality-private latent factors,
#'   common to every modality.
#' @param n_strains Integer number of strains in the genetic panel. Ignored
#'   where `profile` supplies the design.
#' @param n_per_cell Integer number of samples per design cell, per modality.
#'   Samples are drawn independently per modality, so no individual appears in more
#'   than one modality. Ignored where `profile` supplies the design.
#' @param noise_sd Standard deviation of the additive Gaussian measurement
#'   noise.
#' @param effect_size Magnitude of the design signatures carried by the shared
#'   factors, in units of the latent scale. Zero detaches the shared factors
#'   from the design, which is the complete null the recovery benchmark is calibrated
#'   against.
#' @param signature Optional matrix of design weights, one row per shared factor
#'   and one column per design term, shared across modalities. The terms are the
#'   phenotype, age and sex contrasts by default, and the contrasts the profiles
#'   share where `profile` is given; [chorale_signature_terms()] reports them.
#'   When absent each shared factor is given a distinct dominant contrast.
#'   Supplying a signature in which two factors respond to the phenotype the
#'   same way is how a same-response, distinct-component adversarial case is
#'   built.
#' @param confounder Optional named list adding a nuisance covariate correlated
#'   with the phenotype, with elements `name`, `rho` (the correlation) and
#'   `loading` (its weight on the shared factors), for testing robustness to
#'   confounding.
#' @param imbalance Optional numeric in `[0, 1)`. When positive, design cells
#'   are populated unevenly across modalities, for testing behaviour under the
#'   unequal mixtures disjoint cohorts produce.
#' @param profile Optional [chorale_data_profile()], or a list of one per
#'   modality, supplying the marginals, the missingness and the design margins
#'   of real data.
#' @param loadings Optional list, one entry per modality, each a
#'   features-by-factors matrix used in place of the drawn loadings. This is
#'   how [chorale_plant()] places a named pathway on the features that belong
#'   to it.
#' @param background Whether to add correlated variation reproducing the tail of
#'   the profile's correlation eigenspectrum. Without it a simulated matrix
#'   carries only the planted factors and independent noise, and its covariance
#'   is far lower-dimensional than a real one. It has no effect without
#'   `profile`.
#' @param max_background Largest number of background directions drawn.
#' @param seed Integer random seed.
#'
#' @returns A list with `modalities` (feature-by-sample matrices),
#'   `col_data` (per-modality design tables) and `truth`, which carries the
#'   shared and private loadings and scores per modality, the marker indices,
#'   the design `signature`, its `terms`, and the factor counts, for scoring
#'   recovery.
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
                             profile = NULL,
                             loadings = NULL,
                             background = TRUE,
                             max_background = 50L,
                             seed = 1) {
  if (n_modalities < 2) {
    rlang::abort("`n_modalities` must be at least 2.")
  }
  profiles <- chorale_profile_list(profile, n_modalities)
  if (is.null(n_features)) {
    if (is.null(profiles)) {
      rlang::abort("`n_features` is required unless `profile` is supplied.")
    }
    n_features <- vapply(profiles, function(p) as.integer(p$n_features), integer(1))
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
  if (!is.null(loadings) && length(loadings) != n_modalities) {
    rlang::abort("`loadings` must have one entry per modality.")
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
    strain = paste0("group_", seq_len(n_strains)),
    phenotype = c("control", "case"),
    age_months = c(6, 14),
    sex = c("F", "M"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  terms <- chorale_signature_terms(profiles)
  n_terms <- length(terms)

  # Each shared factor has a distinct dominant contrast. The first is
  # phenotype-responsive; factors driven only by a secondary covariate provide
  # the negative case that hierarchical matching must leave unsupported.
  if (is.null(signature)) {
    signature <- matrix(0, nrow = n_shared_factors, ncol = n_terms,
                        dimnames = list(NULL, terms))
    for (k in seq_len(n_shared_factors)) {
      dom <- ((k - 1L) %% n_terms) + 1L
      signature[k, dom] <- 1
    }
  }
  if (ncol(signature) != n_terms) {
    rlang::abort(paste0(
      "`signature` must have one column per design term: ",
      paste(terms, collapse = ", "), "."
    ))
  }
  colnames(signature) <- terms
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
    cells <- if (is.null(profiles)) {
      design[rep(seq_len(nrow(design)), each = n_per_cell), ]
    } else {
      chorale_draw_cells(profiles[[m]])
    }
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

    contrasts <- chorale_sim_contrasts(cells, terms)

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
      # shared factors, so a matched component can be probed for whether it
      # tracks the phenotype or the confounder.
      z <- ifelse(cells$phenotype == "case", 1, -1) * confounder$rho +
        stats::rnorm(n_m) * sqrt(1 - confounder$rho^2)
      shared_scores <- shared_scores + outer(z, rep(confounder$loading, n_shared_factors))
      batch <- ifelse(z > stats::median(z), "batchB", "batchA")
    }

    scores <- cbind(shared_scores, private_scores)

    if (is.null(loadings)) {
      built <- chorale_draw_loadings(p_m, n_shared_factors, n_private_factors)
    } else {
      built <- chorale_supplied_loadings(loadings[[m]], p_m, n_shared_factors,
                                         n_private_factors)
    }
    load_m <- built$loadings
    markers <- built$markers

    signal <- scores %*% t(load_m)
    noise <- matrix(stats::rnorm(n_m * p_m, sd = noise_sd), nrow = n_m, ncol = p_m)

    feature_id <- sprintf("modality%d_feature%05d", m, seq_len(p_m))
    if (is.null(profiles)) {
      x <- t(signal + noise)
    } else {
      pr <- profiles[[m]]
      donor <- chorale_marginal_donors(
        t(signal + noise), pr,
        identity_map = is.null(loadings) && p_m == pr$n_features)
      shape <- function(mat) chorale_apply_marginals(t(mat), pr, donor)
      if (isTRUE(background)) {
        bg <- chorale_background_directions(n_m, p_m, pr, max_background)
        if (!is.null(bg)) {
          signal <- signal + chorale_calibrate_background(
            signal + noise, bg, shape, pr$eigenvalues[1])
        }
      }
      x <- shape(signal + noise)
      feature_id <- chorale_donor_ids(pr, donor)
    }
    if (!is.null(loadings) && !is.null(rownames(loadings[[m]]))) {
      feature_id <- rownames(loadings[[m]])
    }
    dimnames(x) <- list(feature_id, sample_id)

    modalities[[m]] <- x
    col_data[[m]] <- chorale_sim_col_data(cells, sample_id,
                                          names(modalities)[m], batch)
    truth_loadings[[m]] <- list(
      shared = load_m[, seq_len(n_shared_factors), drop = FALSE],
      private = load_m[, n_shared_factors + seq_len(n_private_factors), drop = FALSE]
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
      terms = terms,
      n_shared_factors = n_shared_factors,
      n_private_factors = n_private_factors
    )
  )
}

#' Design terms a simulated signature is written over
#'
#' The columns of a design signature. Without profiles these are the phenotype,
#' age and sex contrasts the parametric design carries. With profiles they are
#' the contrasts every profile's design shares, since a signature is the same in
#' every modality and can only be written over terms every modality holds.
#'
#' @param profile A [chorale_data_profile()], a list of them, or `NULL`.
#'
#' @returns A character vector of design term names.
#' @export
#' @examples
#' chorale_signature_terms(NULL)
chorale_signature_terms <- function(profile) {
  profiles <- if (is.null(profile)) {
    NULL
  } else if (inherits(profile, "chorale_data_profile")) {
    list(profile)
  } else {
    profile
  }
  if (is.null(profiles)) return(c("phenotype", "age", "sex"))

  per <- lapply(profiles, function(p) {
    cells <- p$design_cells
    if (is.null(cells)) {
      rlang::abort("A profile used for simulation must carry design cells.")
    }
    nm <- setdiff(names(cells), "n")
    nm[vapply(nm, function(v) length(unique(cells[[v]][!is.na(cells[[v]])])) >= 2,
              logical(1))]
  })
  shared <- Reduce(intersect, per)
  if (!"phenotype" %in% shared) {
    rlang::abort(paste0(
      "Every profile must carry a varying `phenotype`; shared varying terms: ",
      if (length(shared) == 0) "none" else paste(shared, collapse = ", "), "."
    ))
  }
  # The phenotype leads, since the estimand is the case/control contrast, and
  # the rest follow in a fixed order so a signature means the same thing in
  # every modality.
  c("phenotype", sort(setdiff(shared, "phenotype")))
}

#' Resolve a profile argument to one profile per modality
#' @keywords internal
#' @noRd
chorale_profile_list <- function(profile, n_modalities) {
  if (is.null(profile)) return(NULL)
  if (inherits(profile, "chorale_data_profile")) {
    return(rep(list(profile), n_modalities))
  }
  if (!is.list(profile) ||
      !all(vapply(profile, inherits, logical(1), "chorale_data_profile"))) {
    rlang::abort("`profile` must be a `chorale_data_profile` or a list of them.")
  }
  if (length(profile) != n_modalities) {
    rlang::abort("`profile` must have one entry per modality.")
  }
  profile
}

#' Draw a design from the cells a real cohort populated
#' @keywords internal
#' @noRd
chorale_draw_cells <- function(profile) {
  cells <- profile$design_cells
  if (is.null(cells)) {
    rlang::abort("A profile used for simulation must carry design cells.")
  }
  covariates <- setdiff(names(cells), "n")
  idx <- rep(seq_len(nrow(cells)), times = cells$n)
  cells[idx, covariates, drop = FALSE]
}

#' Signed contrasts for the design terms a signature is written over
#'
#' Every entry is on one scale, so a signature weight means the same thing
#' whichever term it multiplies. A two-level covariate contributes plus or minus
#' a half, a covariate with more levels contributes them evenly spaced across
#' that range, and a numeric covariate is rescaled onto it.
#'
#' @keywords internal
#' @noRd
chorale_sim_contrasts <- function(cells, terms) {
  out <- matrix(0, nrow = nrow(cells), ncol = length(terms),
                dimnames = list(NULL, terms))
  for (j in seq_along(terms)) {
    tm <- terms[j]
    v <- if (tm %in% names(cells)) {
      cells[[tm]]
    } else if (tm == "age" && "age_months" %in% names(cells)) {
      cells$age_months
    } else {
      next
    }
    out[, j] <- chorale_signed_contrast(v, control_first = tm == "phenotype")
  }
  out
}

#' One covariate as a signed contrast in `[-0.5, 0.5]`
#'
#' The reference level is the one that sorts first, so a contrast means the
#' same thing in every modality. For the phenotype the reference is the control
#' arm wherever the label registry recognises it, which puts cases on the
#' positive side of every signature.
#'
#' @keywords internal
#' @noRd
chorale_signed_contrast <- function(v, control_first = FALSE) {
  if (is.numeric(v)) {
    ok <- is.finite(v)
    if (!any(ok)) return(rep(0, length(v)))
    rng <- range(v[ok])
    out <- if (diff(rng) == 0) rep(0, length(v)) else (v - rng[1]) / diff(rng) - 0.5
    out[!ok] <- 0
    return(as.numeric(out))
  }
  ch <- as.character(v)
  lev <- sort(unique(ch[!is.na(ch)]))
  if (length(lev) < 2) return(rep(0, length(ch)))
  if (control_first) {
    map <- chorale_label_registry()$phenotype
    is_control <- map[match(tolower(trimws(lev)), names(map))] %in% "control"
    is_control[is.na(is_control)] <- FALSE
    if (any(is_control) && !all(is_control)) lev <- c(lev[is_control], lev[!is_control])
  }
  pos <- seq(-0.5, 0.5, length.out = length(lev))
  out <- pos[match(ch, lev)]
  out[is.na(out)] <- 0
  as.numeric(out)
}

#' Loadings with two planted pure features per shared factor
#' @keywords internal
#' @noRd
chorale_draw_loadings <- function(p_m, n_shared_factors, n_private_factors) {
  n_total_factors <- n_shared_factors + n_private_factors
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
  list(loadings = loadings, markers = markers)
}

#' Check and adopt a supplied loading matrix
#' @keywords internal
#' @noRd
chorale_supplied_loadings <- function(l, p_m, n_shared_factors, n_private_factors) {
  n_total_factors <- n_shared_factors + n_private_factors
  l <- as.matrix(l)
  if (nrow(l) != p_m || ncol(l) != n_total_factors) {
    rlang::abort(paste0(
      "A supplied loading matrix must be ", p_m, " by ", n_total_factors, "."
    ))
  }
  # A pure feature loads on one factor and negligibly on the others, which is a
  # property of the matrix rather than of how it was built, so it is read off
  # rather than assumed.
  markers <- vector("list", n_shared_factors)
  for (k in seq_len(n_shared_factors)) {
    own <- abs(l[, k])
    other <- if (ncol(l) > 1) {
      apply(abs(l[, -k, drop = FALSE]), 1, max)
    } else {
      rep(0, nrow(l))
    }
    pure <- which(own > 0 & other <= 0.25 * own)
    markers[[k]] <- utils::head(pure[order(-own[pure])], 20L)
  }
  names(markers) <- paste0("factor_", seq_len(n_shared_factors))
  list(loadings = l, markers = markers)
}

#' Choose the real feature each simulated feature takes its marginal from
#'
#' With an identity map a simulated feature keeps the identity of the real one
#' it was built on, which is what a planted pathway needs. Otherwise donors are
#' drawn and assigned by rank of spread, so a simulated feature carrying more
#' signal is given a marginal with more room to carry it, and the mean-variance
#' relation of the real data stays attached to the structure planted in the
#' simulated one.
#'
#' @keywords internal
#' @noRd
chorale_marginal_donors <- function(x, profile, identity_map = FALSE) {
  p_m <- nrow(x)
  q <- profile$quantiles
  usable <- which(stats::complete.cases(q))
  if (length(usable) < 2) {
    rlang::abort("The profile carries fewer than two usable feature marginals.")
  }

  if (identity_map && p_m == nrow(q)) {
    donor <- seq_len(p_m)
    replace_with <- !donor %in% usable
    if (any(replace_with)) {
      donor[replace_with] <- sample(usable, sum(replace_with), replace = TRUE)
    }
    return(donor)
  }

  donor <- sample(usable, p_m, replace = p_m > length(usable))
  sim_spread <- apply(x, 1, stats::sd)
  donor[order(order(sim_spread))] <- donor[order(profile$feature$spread[donor])]
  donor
}

#' Identifiers the donors carry
#' @keywords internal
#' @noRd
chorale_donor_ids <- function(profile, donor) {
  make.unique(as.character(profile$feature$feature_id[donor]))
}

#' Map simulated features onto the marginals of real ones
#'
#' Each simulated feature is replaced by its own ranks read through a real
#' feature's quantile grid, so the values that come out have the distribution
#' the real feature had. The map is monotone, so a planted design response keeps
#' its sign and its ordering. The donor's missingness is applied to its least
#' abundant values, which is how missingness depends on abundance in the data
#' the grids come from.
#'
#' @keywords internal
#' @noRd
chorale_apply_marginals <- function(x, profile, donor) {
  p_m <- nrow(x)
  n_m <- ncol(x)
  q <- profile$quantiles
  probs <- profile$probs
  out <- matrix(NA_real_, nrow = p_m, ncol = n_m)
  for (i in seq_len(p_m)) {
    d <- donor[i]
    r <- rank(x[i, ], ties.method = "first") / (n_m + 1)
    mapped <- stats::approx(probs, q[d, ], xout = r, rule = 2)$y
    k <- floor(n_m * profile$feature$missing[d])
    if (is.finite(k) && k > 0) {
      mapped[order(mapped)[seq_len(min(k, n_m - 2L))]] <- NA_real_
    }
    out[i, ] <- mapped
  }
  out
}

#' Random directions whose variances follow a real eigenspectrum
#'
#' A matrix carrying only its planted factors and independent noise has a
#' covariance of far lower dimension than a real one, and an estimator choosing
#' how many components a modality supports responds to exactly that. These
#' directions restore it. Their scores are Gaussian, so they add covariance
#' without adding anything a component analysis driven by non-Gaussianity could
#' recover as a source, and they are orthonormal in feature space, so they
#' occupy it evenly rather than concentrating on the features a component
#' already uses.
#'
#' @keywords internal
#' @noRd
chorale_background_directions <- function(n, p, profile, max_background = 50L) {
  ev <- profile$eigenvalues
  ev <- ev[is.finite(ev) & ev > 0]
  if (length(ev) < 2) return(NULL)
  b <- min(length(ev), as.integer(max_background), n - 2L, p)
  if (b < 2L) return(NULL)

  s <- scale(matrix(stats::rnorm(n * b), nrow = n, ncol = b))
  s <- sweep(s, 2, sqrt(ev[seq_len(b)] / ev[1]), "*")
  l <- qr.Q(qr(matrix(stats::rnorm(p * b), nrow = p, ncol = b)))
  s %*% t(l)
}

#' How much background reproduces the share the leading direction carries
#'
#' The amplitude is read off a grid rather than solved in closed form, for two
#' reasons. The share the leading direction carries falls and then rises as the
#' background grows, since a little of it flattens the planted spectrum and a
#' lot of it imposes its own, so there is no single crossing to solve for. And
#' the quantity that has to match is measured on the matrix the estimator will
#' read, after the marginals have been imposed, which no expression in the
#' latent variances predicts.
#'
#' @keywords internal
#' @noRd
chorale_calibrate_background <- function(base, bg, shape, target) {
  zero <- matrix(0, nrow = nrow(bg), ncol = ncol(bg))
  if (!is.finite(target)) return(zero)
  # The directions are put on the scale of the signal first, so the amplitudes
  # tried mean the same thing whatever the feature count and the factor count
  # happen to be.
  scale_base <- stats::sd(base)
  scale_bg <- stats::sd(bg)
  if (!is.finite(scale_bg) || scale_bg <= 0) return(zero)
  bg <- bg * (scale_base / scale_bg)
  leading <- function(cc) {
    mat <- shape(base + cc * bg)
    z <- scale(t(mat))
    z[!is.finite(z)] <- 0
    e <- eigen(tcrossprod(z), symmetric = TRUE, only.values = TRUE)$values
    e <- pmax(e, 0)
    if (sum(e) <= 0) return(NA_real_)
    e[1] / sum(e)
  }
  grid <- c(0, 10^seq(-1, 2.5, length.out = 15))
  achieved <- vapply(grid, leading, numeric(1))
  ok <- is.finite(achieved)
  if (!any(ok)) return(zero)
  best <- grid[ok][which.min(abs(achieved[ok] - target))]
  if (best == 0) return(zero)
  bg * best
}

#' Assemble the per-modality design table
#' @keywords internal
#' @noRd
chorale_sim_col_data <- function(cells, sample_id, modality, batch) {
  out <- data.frame(sample_id = sample_id, stringsAsFactors = FALSE)
  if (!"cohort" %in% names(cells)) out$cohort <- "simulated"
  out$modality <- modality
  for (nm in setdiff(names(cells), c("sample_id", "modality", "batch"))) {
    out[[nm]] <- cells[[nm]]
  }
  if (!"region" %in% names(out)) out$region <- "simulated"
  # A confounder writes a batch, so it is set here rather than taken from the
  # profile, whose batch labels belong to the cohort it was measured on.
  out$batch <- batch
  out
}
