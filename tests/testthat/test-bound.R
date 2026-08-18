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

test_that("imbalanced strata still give a valid identified set", {
  # The two modalities realise the design in very different proportions, which
  # is the case that produced endpoints outside [-1, 1] before both marginals
  # were reweighted to one target.
  set.seed(3)
  ka <- rep(c("s1", "s2"), times = c(60, 6))
  kb <- rep(c("s1", "s2"), times = c(6, 60))
  va <- stats::rnorm(length(ka)) + ifelse(ka == "s1", -2, 3)
  vb <- stats::rnorm(length(kb)) * 4 + ifelse(kb == "s1", -5, 6)
  names(va) <- paste0("a", seq_along(va))
  names(vb) <- paste0("b", seq_along(vb))

  w <- chorale:::chorale_target_weights(ka, kb, c("s1", "s2"))
  out <- chorale:::chorale_anchored_correlation(va, vb, ka, kb, c("s1", "s2"),
                                                200L, w)
  expect_true(is.finite(out$lower) && is.finite(out$upper))
  expect_gte(out$lower, -1)
  expect_lte(out$upper, 1)
  expect_lte(out$lower, out$upper)

  # Conditioning on the design can only narrow the range, never widen it, so
  # the two intervals must be computed on the same target population.
  probs <- seq(0.5 / 200, 1 - 0.5 / 200, length.out = 200)
  unc <- chorale:::chorale_frechet_from_quantiles(
    chorale:::chorale_weighted_quantile(va, w$a, probs),
    chorale:::chorale_weighted_quantile(vb, w$b, probs)
  )
  expect_lte(out$upper - out$lower, unc$upper - unc$lower + 1e-8)
})

test_that("a known coupling lies inside the reported set", {
  set.seed(4)
  ka <- kb <- rep(c("s1", "s2"), each = 40)
  u <- stats::rnorm(80)
  va <- u + ifelse(ka == "s1", -1, 1)
  vb <- 0.7 * u + stats::rnorm(80, sd = 0.7) + ifelse(kb == "s1", -1, 1)
  names(va) <- paste0("a", 1:80)
  names(vb) <- paste0("b", 1:80)
  truth <- stats::cor(va, vb)

  w <- chorale:::chorale_target_weights(ka, kb, c("s1", "s2"))
  out <- chorale:::chorale_anchored_correlation(va, vb, ka, kb, c("s1", "s2"),
                                                200L, w)
  expect_gte(truth, out$lower - 1e-6)
  expect_lte(truth, out$upper + 1e-6)
})

test_that("endpoints are clamped independently", {
  # Clamping one endpoint using the other can leave a set no correlation
  # satisfies, which is the defect this guards.
  expect_equal(chorale:::chorale_clamp_bounds(1.4, 0.9),
               list(lower = 0.9, upper = 1))
  expect_equal(chorale:::chorale_clamp_bounds(-3, 2),
               list(lower = -1, upper = 1))
  expect_equal(chorale:::chorale_clamp_bounds(NA, 0.5),
               list(lower = -1, upper = 0.5))
})
