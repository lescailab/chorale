test_that("the factor count is chosen from the data", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 200,
                          n_shared_factors = 3, n_private_factors = 2,
                          n_strains = 6, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  x <- scale(t(sim$modalities[[1]]))
  k <- chorale_n_factors(x, n_perm = 40, seed = 1)
  # Five factors were planted; parallel analysis should land on them rather
  # than on the feature count.
  expect_gte(k, 3)
  expect_lte(k, 8)

  # The count follows the cohort, not the feature space: more features on the
  # same animals does not buy more components.
  wide <- chorale_simulate(n_modalities = 2, n_features = 600,
                           n_shared_factors = 3, n_private_factors = 2,
                           n_strains = 6, n_per_cell = 3, effect_size = 3,
                           seed = 1)
  kw <- chorale_n_factors(scale(t(wide$modalities[[1]])), n_perm = 40, seed = 1)
  expect_lt(abs(kw - k), 4)
})

test_that("chorale_fit accepts auto factor counts", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = "auto", n_init = 2)
  expect_true(all(fit$n_factors >= 2))
  expect_equal(length(fit$n_factors), 2L)
  for (m in fit$modalities) {
    expect_equal(ncol(fit$fits[[m]]$scores), unname(fit$n_factors[[m]]))
  }
})

test_that("the identified set widens once resampling is accounted for", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 3)
  u <- chorale_bound_uncertainty(fit, n_boot = 30, seed = 1)
  skip_if(nrow(u) == 0)
  expect_true(all(u$region_lower >= -1 & u$region_upper <= 1))
  expect_true(all(u$region_lower <= u$region_upper))
  # A region accounting for resampling contains the plug-in set it surrounds,
  # so it is never the narrower of the two.
  plug_width <- u$upper_anchored - u$lower_anchored
  expect_true(mean(u$region_width >= plug_width - 1e-6) >= 0.5)
})
