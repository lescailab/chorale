test_that("the transform is read from the matrix", {
  set.seed(1)
  # Sequencing counts: whole, non-negative, spanning orders of magnitude.
  mu <- 2^runif(400, 1, 14)
  counts <- matrix(stats::rnbinom(400 * 30, mu = mu, size = 5), nrow = 400)
  expect_equal(chorale_transform(counts)$applied, "vst")

  # Label-free intensities: continuous, non-negative, heavily right-skewed.
  intens <- matrix(stats::rlnorm(400 * 30, meanlog = 5, sdlog = 2), nrow = 400)
  expect_equal(chorale_transform(intens)$applied, "log")

  # Already on a log scale, so nothing to do.
  logged <- matrix(stats::rnorm(400 * 30, mean = 20, sd = 2), nrow = 400)
  expect_equal(chorale_transform(logged)$applied, "none")

  # Negative values are already symmetric, typically a ratio.
  ratios <- matrix(stats::rnorm(400 * 30), nrow = 400)
  expect_equal(chorale_transform(ratios)$applied, "none")
})

test_that("each transform does what it says", {
  m <- matrix(c(0, 1, 10, 100), nrow = 2)
  expect_equal(chorale_transform(m, "none")$matrix, m)
  expect_equal(chorale_transform(m, "log")$matrix, log1p(m))
  expect_equal(chorale_transform(m, "vst")$matrix, 2 * sqrt(m + 3 / 8))
  # Zeros stay finite, which a plain logarithm would not allow.
  expect_true(all(is.finite(chorale_transform(m, "log")$matrix)))
})

test_that("the variance-stabilising transform removes the mean-variance trend", {
  set.seed(2)
  mu <- rep(c(5, 50, 500, 5000), each = 50)
  counts <- t(vapply(mu, function(u) stats::rpois(60, u), numeric(60)))
  raw_trend <- stats::cor(rowMeans(counts), apply(counts, 1, stats::var))
  vst <- chorale_transform(counts, "vst")$matrix
  vst_trend <- stats::cor(rowMeans(vst), apply(vst, 1, stats::var))
  # Under the transform a feature's variance no longer follows its mean.
  expect_gt(raw_trend, 0.8)
  expect_lt(abs(vst_trend), raw_trend)
})

test_that("chorale_fit records the transform it applied", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
  applied <- vapply(fit$modalities, function(m) fit$fits[[m]]$transform,
                    character(1))
  expect_length(applied, 2L)
  expect_true(all(applied %in% c("auto", "none", "log", "vst")))

  # A per-modality request is honoured, and an unknown one is refused.
  fit2 <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
                      transform = c(modality_1 = "log"))
  expect_equal(unname(fit2$fits$modality_1$transform), "log")
  expect_error(chorale_fit(containers, n_factors = 2, n_init = 2,
                           transform = "quantile"), class = "rlang_error")
  expect_error(chorale_fit(containers, n_factors = 2, n_init = 2,
                           transform = c(nope = "log")), class = "rlang_error")
})

test_that("the medoid is the factor run the starts agree on", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 150,
                          n_shared_factors = 3, n_private_factors = 2,
                          n_strains = 5, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  x <- scale(t(sim$modalities[[1]]))
  truth <- cbind(sim$truth$scores[[1]]$shared, sim$truth$scores[[1]]$private)
  agree <- function(ic) {
    a <- abs(stats::cor(ic$scores, truth))
    m <- clue::solve_LSAP(a, maximum = TRUE)
    mean(a[cbind(seq_len(nrow(a)), as.integer(m))])
  }
  single <- chorale_ica(x, 5, n_init = 8, seed = 1, consensus = FALSE)
  cons <- chorale_ica(x, 5, n_init = 8, seed = 1, consensus = TRUE)
  expect_equal(dim(cons$scores), dim(single$scores))
  # The stable representative recovers the planted factors at least as well as
  # the legacy switch, which now maps to the same valid ICA medoid.
  expect_gt(agree(cons), agree(single) - 0.05)
  # Scores stay standardised, which the downstream profiles assume.
  expect_true(all(abs(colMeans(cons$scores)) < 1e-8))
})

test_that("run alignment resolves the order and sign a factor is recovered in", {
  set.seed(3)
  ref <- scale(matrix(stats::rnorm(120), ncol = 3))
  # The same factors, permuted and sign-flipped, which is all that separates
  # two initialisations of the same solution.
  other <- ref[, c(3, 1, 2)] %*% diag(c(-1, 1, -1))
  aligned <- chorale_align_runs(ref, list(other))
  expect_length(aligned, 1L)
  expect_equal(as.numeric(diag(stats::cor(ref, aligned[[1]]))), rep(1, 3),
               tolerance = 1e-8)
})
