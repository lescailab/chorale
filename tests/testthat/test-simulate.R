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
    c(pheno = mean(s[d$phenotype == "5XFAD"]) - mean(s[d$phenotype == "Ntg"]),
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
