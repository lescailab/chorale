test_that("the reduced basis recovers the same components where they are determined", {
  skip_if_not_installed("fastICA")
  set.seed(1)
  n <- 60L
  p <- 1200L
  sources <- matrix(stats::rt(n * 3L, df = 3), n, 3L)
  x <- scale(sources %*% matrix(stats::rnorm(3L * p), 3L, p) +
               matrix(stats::rnorm(n * p), n, p))
  x[!is.finite(x)] <- 0

  # Neither sign nor order is identified by independent component analysis, so
  # the components are matched one to one before anything is read off. At the
  # count the data actually support the same components come back; above it the
  # optimiser is choosing among directions that are not determined, which is
  # what `chorale_select_factors()` exists to detect.
  for (k in c(1L, 3L)) {
    set.seed(2)
    full <- fastICA::fastICA(x, n.comp = k, method = "C", maxit = 500,
                             tol = 1e-5)$S
    set.seed(2)
    reduced <- fastICA::fastICA(chorale_ica_basis(x), n.comp = k, method = "C",
                                maxit = 500, tol = 1e-5)$S
    agreement <- abs(stats::cor(scale(full), scale(reduced)))
    matched <- agreement[cbind(seq_len(k),
                               as.integer(clue::solve_LSAP(agreement,
                                                           maximum = TRUE)))]
    expect_gt(min(matched), 0.99)
  }
})

test_that("the reduced basis spans the row space and is left alone when narrow", {
  set.seed(3)
  wide <- matrix(stats::rnorm(20L * 300L), 20L, 300L)
  z <- chorale_ica_basis(wide)
  expect_equal(nrow(z), 20L)
  expect_lte(ncol(z), 20L)
  # Same geometry between samples: distances and inner products are preserved.
  expect_equal(tcrossprod(z), tcrossprod(scale(wide, scale = FALSE)),
               tolerance = 1e-8, ignore_attr = TRUE)

  narrow <- matrix(stats::rnorm(40L * 10L), 40L, 10L)
  expect_identical(chorale_ica_basis(narrow), narrow)
})

test_that("chorale_ica returns loadings on the original features", {
  skip_if_not_installed("fastICA")
  set.seed(4)
  n <- 50L
  p <- 800L
  x <- scale(matrix(stats::rt(n * 2L, df = 3), n, 2L) %*%
               matrix(stats::rnorm(2L * p), 2L, p) +
               matrix(stats::rnorm(n * p), n, p))
  colnames(x) <- paste0("feature_", seq_len(p))
  rownames(x) <- paste0("sample_", seq_len(n))
  fit <- chorale_ica(x, n_factors = 2L, n_init = 2L, seed = 1L)
  expect_equal(nrow(fit$loadings), p)
  expect_equal(rownames(fit$loadings), colnames(x))
  expect_equal(rownames(fit$scores), rownames(x))
})
