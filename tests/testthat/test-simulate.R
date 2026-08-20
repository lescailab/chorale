test_that("chorale_simulate produces the requested dimensions", {
  sim <- chorale_simulate(
    n_modalities = 3, n_features = c(30, 40, 50), n_shared_factors = 3,
    n_private_factors = 2, n_strains = 4, n_per_cell = 2, seed = 1
  )

  expect_length(sim$modalities, 3)
  expect_equal(nrow(sim$modalities[[1]]), 30)
  expect_equal(nrow(sim$modalities[[2]]), 40)
  expect_equal(nrow(sim$modalities[[3]]), 50)

  n_expected_samples <- 4 * 2 * 2 * 2 * 2 # strains x phenotype x age x sex x n_per_cell
  expect_equal(ncol(sim$modalities[[1]]), n_expected_samples)
  expect_equal(nrow(sim$col_data[[1]]), n_expected_samples)
})

test_that("chorale_simulate samples are disjoint across modalities", {
  sim <- chorale_simulate(
    n_modalities = 2, n_features = 30, n_shared_factors = 2,
    n_private_factors = 1, n_strains = 3, n_per_cell = 1, seed = 2
  )
  ids <- unlist(lapply(sim$col_data, `[[`, "sample_id"))
  expect_equal(length(ids), length(unique(ids)))
})

test_that("chorale_simulate col_data carries the required design columns", {
  sim <- chorale_simulate(
    n_modalities = 2, n_features = 30, n_shared_factors = 2,
    n_private_factors = 1, n_strains = 3, n_per_cell = 1, seed = 3
  )
  required <- c(
    "sample_id", "cohort", "modality", "strain", "phenotype",
    "age_months", "sex", "region", "batch"
  )
  expect_true(all(required %in% names(sim$col_data[[1]])))
})

test_that("chorale_simulate marker features are pure within each modality", {
  sim <- chorale_simulate(
    n_modalities = 2, n_features = 30, n_shared_factors = 2,
    n_private_factors = 1, n_strains = 3, n_per_cell = 1, seed = 4
  )
  for (m in seq_along(sim$modalities)) {
    shared_loadings <- sim$truth$loadings[[m]]$shared
    for (k in seq_along(sim$truth$markers[[m]])) {
      marker_idx <- sim$truth$markers[[m]][[k]]
      other_factors <- setdiff(seq_len(ncol(shared_loadings)), k)
      expect_true(all(shared_loadings[marker_idx, other_factors] == 0))
      expect_true(all(shared_loadings[marker_idx, k] != 0))
    }
  }
})

test_that("chorale_simulate is reproducible given the same seed", {
  sim1 <- chorale_simulate(
    n_modalities = 2, n_features = 30, n_shared_factors = 2,
    n_private_factors = 1, n_strains = 3, n_per_cell = 1, seed = 5
  )
  sim2 <- chorale_simulate(
    n_modalities = 2, n_features = 30, n_shared_factors = 2,
    n_private_factors = 1, n_strains = 3, n_per_cell = 1, seed = 5
  )
  expect_equal(sim1$modalities, sim2$modalities)
})

test_that("chorale_simulate rejects too few features for the requested factors", {
  expect_error(
    chorale_simulate(
      n_modalities = 2, n_features = 2, n_shared_factors = 3,
      n_private_factors = 2, n_strains = 2, n_per_cell = 1, seed = 1
    )
  )
})

test_that("chorale_simulate requires at least two modalities", {
  expect_error(
    chorale_simulate(n_modalities = 1, n_features = 30)
  )
})

test_that("planted sources satisfy the identification conditions", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 150,
                          n_shared_factors = 3, n_private_factors = 2,
                          n_strains = 4, n_per_cell = 4, effect_size = 3,
                          seed = 1)
  for (mi in seq_along(sim$modalities)) {
    s <- cbind(sim$truth$scores[[mi]]$shared, sim$truth$scores[[mi]]$private)
    kurt <- apply(s, 2, function(v) mean(scale(v)^4) - 3)
    # Non-Gaussian: every factor departs from the normal in kurtosis. A shared
    # factor carrying a strong two-point design response is bimodal rather than
    # skewed, which is non-Gaussian in the other direction, so kurtosis is the
    # condition that holds for all of them.
    expect_true(all(abs(kurt) > 0.1))
    # Non-symmetric: the private factors are the pure planted sources, with no
    # design response to symmetrise them, so their skewness stays away from zero.
    priv <- sim$truth$scores[[mi]]$private
    skew_priv <- apply(priv, 2, function(v) mean(scale(v)^3))
    expect_true(all(abs(skew_priv) > 0.2))
    # Design responses are near-orthogonal, so the shared factors stay
    # separable for independent component analysis.
    shared <- sim$truth$scores[[mi]]$shared
    cc <- abs(stats::cor(shared)[upper.tri(diag(ncol(shared)))])
    expect_lt(max(cc), 0.3)
  }
})

test_that("a shared factor carries the same design signature in every modality", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4, effect_size = 4,
                          seed = 2)
  profile <- function(mi, k) {
    s <- sim$truth$scores[[mi]]$shared[, k]
    d <- sim$col_data[[mi]]
    c(pheno = mean(s[d$phenotype == "case"]) - mean(s[d$phenotype == "control"]),
      age = mean(s[d$age_months == 14]) - mean(s[d$age_months == 6]))
  }
  # The first shared factor is phenotype-dominant in all three modalities, which
  # is the correspondence the estimator must recover.
  p1 <- vapply(seq_along(sim$modalities), profile, numeric(2), k = 1)
  expect_true(all(p1["pheno", ] > p1["age", ]))
})

test_that("the confounder and imbalance options change the data as intended", {
  conf <- chorale_simulate(n_modalities = 2, n_features = 80,
                           n_shared_factors = 2, n_private_factors = 1,
                           n_strains = 4, n_per_cell = 3, effect_size = 3,
                           confounder = list(name = "batch", rho = 0.8,
                                             loading = 1.5), seed = 1)
  expect_true("batchB" %in% conf$col_data[[1]]$batch)

  imb <- chorale_simulate(n_modalities = 2, n_features = 80,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4, effect_size = 3,
                          imbalance = 0.5, seed = 1)
  # Thinning cells unevenly leaves the two modalities on different sample counts.
  expect_true(ncol(imb$modalities[[1]]) != ncol(imb$modalities[[2]]) ||
                nrow(imb$col_data[[1]]) < 256)
})

test_that("a profile gives the simulation the marginals of the data it came from", {
  base <- chorale_simulate(n_modalities = 2, n_features = 40,
                           n_shared_factors = 2, n_private_factors = 1,
                           n_strains = 3, n_per_cell = 2, seed = 11)
  real <- exp(base$modalities[[1]])
  p <- chorale_data_profile(real, base$col_data[[1]],
                            covariates = c("phenotype", "sex"))
  sim <- chorale_simulate(n_modalities = 2, n_features = 40,
                          n_shared_factors = 2, n_private_factors = 1,
                          profile = p, seed = 12)
  expect_true(all(sim$modalities[[1]] > 0, na.rm = TRUE))
  # Every value comes from a real feature's own quantile grid, so none can fall
  # outside the range the real matrix covered.
  expect_gte(min(sim$modalities[[1]], na.rm = TRUE), min(real))
  expect_lte(max(sim$modalities[[1]], na.rm = TRUE), max(real))
  expect_equal(stats::median(sim$modalities[[1]], na.rm = TRUE),
               stats::median(real), tolerance = 0.1)
})

test_that("a profile supplies the design, and unpopulated cells stay empty", {
  d <- data.frame(sample_id = paste0("s", 1:12),
                  phenotype = rep(c("case", "control"), each = 6),
                  sex = c(rep("F", 6), rep("M", 6)),
                  stringsAsFactors = FALSE)
  x <- matrix(stats::rnorm(12 * 30), nrow = 30, dimnames = list(NULL, d$sample_id))
  p <- chorale_data_profile(x, d, covariates = c("phenotype", "sex"))
  sim <- chorale_simulate(n_modalities = 2, n_features = 30,
                          n_shared_factors = 1, n_private_factors = 1,
                          profile = p, seed = 13)
  cells <- unique(paste(sim$col_data[[1]]$phenotype, sim$col_data[[1]]$sex))
  expect_setequal(cells, c("case F", "control M"))
  expect_equal(ncol(sim$modalities[[1]]), 12)
})

test_that("a profile fixes the design terms a signature is written over", {
  d <- data.frame(sample_id = paste0("s", 1:8),
                  phenotype = rep(c("case", "control"), 4),
                  sex = rep(c("F", "M"), each = 4),
                  stringsAsFactors = FALSE)
  x <- matrix(stats::rnorm(8 * 30), nrow = 30, dimnames = list(NULL, d$sample_id))
  p <- chorale_data_profile(x, d, covariates = c("phenotype", "sex"))
  expect_equal(chorale_signature_terms(p), c("phenotype", "sex"))
  expect_equal(chorale_signature_terms(NULL), c("phenotype", "age", "sex"))
  expect_error(
    chorale_simulate(n_modalities = 2, n_features = 30, n_shared_factors = 1,
                     n_private_factors = 1, profile = p,
                     signature = matrix(1, nrow = 1, ncol = 3)),
    "one column per design term"
  )
})

test_that("a design carrying no varying phenotype cannot be simulated from", {
  d <- data.frame(sample_id = paste0("s", 1:6), phenotype = "case",
                  sex = rep(c("F", "M"), 3), stringsAsFactors = FALSE)
  x <- matrix(stats::rnorm(6 * 20), nrow = 20, dimnames = list(NULL, d$sample_id))
  p <- chorale_data_profile(x, d, covariates = c("phenotype", "sex"))
  expect_error(chorale_signature_terms(p), "varying `phenotype`")
})

test_that("the background reproduces the share the leading direction carries", {
  set.seed(31)
  n <- 60
  common <- stats::rnorm(n)
  real <- t(outer(common, stats::rnorm(80)) + matrix(stats::rnorm(n * 80, sd = 0.3), n, 80))
  colnames(real) <- paste0("s", seq_len(n))
  d <- data.frame(sample_id = colnames(real),
                  phenotype = rep(c("case", "control"), length.out = n),
                  sex = rep(c("F", "M"), each = n / 2), stringsAsFactors = FALSE)
  p <- chorale_data_profile(real, d, covariates = c("phenotype", "sex"))

  leading <- function(mat) {
    z <- scale(t(mat))
    z[!is.finite(z)] <- 0
    e <- pmax(eigen(tcrossprod(z), symmetric = TRUE, only.values = TRUE)$values, 0)
    e[1] / sum(e)
  }
  with_bg <- chorale_simulate(n_modalities = 2, n_features = 80,
                              n_shared_factors = 2, n_private_factors = 1,
                              profile = p, background = TRUE, seed = 32)
  without <- chorale_simulate(n_modalities = 2, n_features = 80,
                              n_shared_factors = 2, n_private_factors = 1,
                              profile = p, background = FALSE, seed = 32)
  target <- p$eigenvalues[1]
  expect_lt(abs(leading(with_bg$modalities[[1]]) - target),
            abs(leading(without$modalities[[1]]) - target))
})

test_that("supplied loadings are used and their pure features read off", {
  p <- chorale_data_profile(
    matrix(stats::rnorm(30 * 12), nrow = 30,
           dimnames = list(paste0("g", 1:30), paste0("s", 1:12))),
    data.frame(sample_id = paste0("s", 1:12),
               phenotype = rep(c("case", "control"), 6),
               sex = rep(c("F", "M"), each = 6), stringsAsFactors = FALSE),
    covariates = c("phenotype", "sex"))
  l <- matrix(0.01, nrow = 30, ncol = 2, dimnames = list(paste0("g", 1:30), NULL))
  l[1:5, 1] <- 1
  l[1:5, 2] <- 0
  sim <- chorale_simulate(n_modalities = 2, n_features = NULL,
                          n_shared_factors = 1, n_private_factors = 1,
                          profile = p, loadings = list(l, l), seed = 41)
  expect_equal(rownames(sim$modalities[[1]]), paste0("g", 1:30))
  expect_true(all(sim$truth$markers[[1]][[1]] %in% 1:5))
})
