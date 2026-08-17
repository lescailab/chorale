test_that("chorale_simulate produces the requested dimensions", {
  sim <- chorale_simulate(
    n_modalities = 3, n_features = c(30, 40, 50), n_shared_factors = 3,
    n_private_factors = 2, n_strains = 4, n_per_cell = 2, seed = 1
  )

  expect_length(sim$modalities, 3)
  expect_equal(nrow(sim$modalities[[1]]), 30)
  expect_equal(nrow(sim$modalities[[2]]), 40)
  expect_equal(nrow(sim$modalities[[3]]), 50)

  n_expected_samples <- 4 * 2 * 2 * 2 * 2 # strains x genotype x age x sex x n_per_cell
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
    "sample_id", "cohort", "modality", "strain", "genotype",
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
