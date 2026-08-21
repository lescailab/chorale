select_matrix <- function(n_factors = 2, n_samples = 60, n_features = 40,
                          noise = 2, seed = 1) {
  set.seed(seed)
  # Sources with heavy tails, so independent component analysis has something
  # to recover; anything beyond `n_factors` is noise with no reproducible
  # structure to find.
  s <- matrix(stats::rt(n_samples * n_factors, df = 3),
              nrow = n_samples, ncol = n_factors)
  l <- matrix(stats::rnorm(n_features * n_factors), nrow = n_factors)
  x <- s %*% l + matrix(stats::rnorm(n_samples * n_features, sd = noise),
                        nrow = n_samples)
  scale(x)
}

test_that("the table carries one row per candidate count", {
  x <- select_matrix()
  sel <- chorale_select_factors(x, max_factors = 4, n_init = 3, n_subsample = 3)

  expect_equal(sel$n_factors, 1:4)
  expect_true(all(c("init_weakest", "init_typical", "subsample_weakest",
                    "subsample_typical", "weakest", "admissible") %in%
                    names(sel)))
  expect_true(all(sel$weakest >= 0 & sel$weakest <= 1, na.rm = TRUE))
  expect_equal(attr(sel, "threshold"), 0.75)
})

test_that("the selected count stops where components stop reproducing", {
  x <- select_matrix(n_factors = 2)
  sel <- chorale_select_factors(x, max_factors = 5, n_init = 4, n_subsample = 4,
                                threshold = 0.8)
  selected <- attr(sel, "selected")

  expect_equal(selected, 2L)
  # Every count at or below the selection is admissible; the first count above
  # it is not, which is what stopped the selection there.
  expect_true(all(sel$admissible[sel$n_factors <= selected]))
  expect_false(sel$admissible[sel$n_factors == selected + 1L])
})

test_that("a stricter threshold cannot select a larger count", {
  x <- select_matrix(n_factors = 2)
  loose <- attr(chorale_select_factors(x, max_factors = 4, n_init = 3,
                                       n_subsample = 3, threshold = 0.6),
                "selected")
  strict <- attr(chorale_select_factors(x, max_factors = 4, n_init = 3,
                                        n_subsample = 3, threshold = 0.98),
                 "selected")
  expect_lte(strict, loose)
})

test_that("the selected count survives a change of seed", {
  x <- select_matrix(n_factors = 2)
  a <- chorale_select_factors(x, max_factors = 4, n_init = 4, n_subsample = 4,
                              seed = 1)
  b <- chorale_select_factors(x, max_factors = 4, n_init = 4, n_subsample = 4,
                              seed = 202)
  expect_equal(attr(a, "selected"), attr(b, "selected"))
})

test_that("selection reads the assay alone", {
  x <- select_matrix()
  # There is no design argument to pass, so the count cannot depend on one:
  # the same matrix under any labelling gives the same table.
  a <- chorale_select_factors(x, max_factors = 3, n_init = 3, n_subsample = 3)
  rownames(x) <- paste0("relabelled_", seq_len(nrow(x)))
  b <- chorale_select_factors(x, max_factors = 3, n_init = 3, n_subsample = 3)
  expect_equal(a, b)
  expect_false(any(c("design", "phenotype", "col_data") %in%
                     names(formals(chorale_select_factors))))
})

test_that("the subsample criterion can be switched off and is then absent", {
  x <- select_matrix()
  sel <- chorale_select_factors(x, max_factors = 3, n_init = 3, n_subsample = 0)
  expect_true(all(is.na(sel$subsample_weakest)))
  expect_true(all(is.finite(sel$init_weakest)))
})

test_that("impossible requests are refused rather than guessed at", {
  x <- select_matrix()
  expect_error(chorale_select_factors(x, threshold = 0), "threshold")
  expect_error(chorale_select_factors(x, subsample_fraction = 1),
               "subsample_fraction")
  expect_error(chorale_select_factors(x[1:3, ]), "at least four samples")
})

test_that("the count chosen for the free dimensions travels with the encoding", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_init = 2,
                        n_select_init = 2, n_subsample = 2, max_factors = 3L)
  for (m in names(fx$containers)) {
    sel <- enc$encodings[[m]]$selection
    expect_s3_class(sel, "data.frame")
    expect_equal(attr(sel, "selected"),
                 enc$variance$n_free[enc$variance$modality == m])
  }
})
