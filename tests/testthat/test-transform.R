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


test_that("the transform each modality was read on travels with the encoding", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
  ids <- sprintf("feature_%05d", seq_len(60))
  sim$modalities <- lapply(sim$modalities, function(m) {
    rownames(m) <- ids
    m
  })
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  names(containers) <- c("A", "B")
  cc <- chorale_concepts(containers, list(one = ids[1:20], two = ids[21:45]),
                         min_features = 5)

  enc <- chorale_encode(containers, cc, n_free = 0)
  expect_true(all(vapply(enc$encodings, function(e) nzchar(e$transform),
                         logical(1))))

  # A modality whose scale is declared is read on the declared one.
  named <- chorale_encode(containers, cc, n_free = 0,
                          transform = c(A = "none", B = "log"))
  expect_equal(named$encodings$A$transform, "none")
  expect_equal(named$encodings$B$transform, "log")
  expect_error(chorale_encode(containers, cc, transform = c(Z = "log")),
               "unknown modalities")
  expect_error(chorale_encode(containers, cc, transform = "nonsense"),
               "must be one of")
})
