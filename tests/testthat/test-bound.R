make_fit <- function(effect_size = 3, seed = 1) {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 5, n_per_cell = 3,
                          effect_size = effect_size, seed = seed)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  list(containers = containers,
       fit = chorale_fit(containers, n_factors = c(3, 3), n_init = 2, seed = seed))
}

test_that("bounds are returned for every matched pair", {
  f <- make_fit()
  b <- chorale_bound(f$fit)
  expect_s3_class(b, "chorale_bound")
  expect_equal(nrow(b$bounds), nrow(f$fit$matches))
  expect_output(print(b), "chorale_bound")
})

test_that("bounds contain the achievable range and are ordered", {
  f <- make_fit()
  b <- chorale_bound(f$fit)
  skip_if(nrow(b$bounds) == 0)
  expect_true(all(b$bounds$lower_no_anchor <= b$bounds$upper_no_anchor))
  expect_true(all(b$bounds$lower_anchored <= b$bounds$upper_anchored))
  expect_true(all(b$bounds$lower_no_anchor >= -1.0001))
  expect_true(all(b$bounds$upper_no_anchor <= 1.0001))
})

test_that("anchors do not widen the identified set", {
  f <- make_fit()
  b <- chorale_bound(f$fit)
  skip_if(nrow(b$bounds) == 0)
  # Conditioning on the design can only remove couplings, never add them.
  expect_true(all(b$bounds$width_anchored <= b$bounds$width_no_anchor + 1e-6))
  expect_true(all(b$bounds$narrowing >= -1e-6))
})

test_that("the unconditional bound is near the full range", {
  # With disjoint samples the joint is not identified, so the honest report
  # is an interval close to the whole of [-1, 1].
  f <- make_fit()
  b <- chorale_bound(f$fit)
  skip_if(nrow(b$bounds) == 0)
  expect_gt(median(b$bounds$width_no_anchor), 1.5)
})

test_that("a fit with no matches yields an empty bound table", {
  f <- make_fit(effect_size = 0, seed = 7)
  b <- chorale_bound(f$fit)
  expect_equal(nrow(b$bounds), nrow(f$fit$matches))
  expect_s3_class(b, "chorale_bound")
})

test_that("chorale_bound rejects a non-fit", {
  expect_error(chorale_bound(list()), class = "rlang_error")
})
