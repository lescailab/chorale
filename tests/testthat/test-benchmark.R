paired_data <- function(seed = 1, effect_size = 3, n_features = 150) {
  sim <- chorale_simulate(n_modalities = 2, n_features = n_features,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4,
                          effect_size = effect_size, seed = seed)
  a <- sim$modalities[[1]]
  b <- sim$modalities[[2]][, seq_len(ncol(a))]
  colnames(b) <- colnames(a)
  list(a = a, b = b, design = sim$col_data[[1]])
}

test_that("the estimator recovers a correspondence it was not shown", {
  p <- paired_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, n_factors = 3, n_init = 5,
                               n_random = 100, seed = 1)
  s <- r$summary
  expect_equal(s$n_samples, ncol(p$a))
  # The recovery must sit above random alignment; placing at or below it is a
  # failure whatever the p-values say.
  expect_gt(s$recovered_agreement, s$random_lower_bound)
  expect_gt(s$placement_between_bounds, 0.5)
  # On clean data the withheld partner should be recovered for most factors
  # whose phenotype evidence is supported. Phenotype-uninformative factors are
  # outside the production estimand and are not forced into a match.
  expect_gt(s$fraction_partner_correct, 0.5)
})

test_that("the pairing is withheld from the estimator", {
  p <- paired_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, n_factors = 3, n_init = 3,
                               n_random = 20, seed = 1)
  # The two modalities are presented under their own identifiers, so nothing in
  # the fit could reveal that a sample in one is a sample in the other.
  ids_a <- rownames(r$fit$fits$A$scores)
  ids_b <- rownames(r$fit$fits$B$scores)
  expect_length(intersect(ids_a, ids_b), 0)
  expect_true(all(grepl("^a_", ids_a)))
  expect_true(all(grepl("^b_", ids_b)))
})

test_that("a benchmark without enough shared samples is refused", {
  p <- paired_data()
  small_a <- p$a[, 1:10, drop = FALSE]
  small_b <- p$b[, 1:10, drop = FALSE]
  expect_error(
    chorale_destroy_pairing(small_a, small_b, p$design, n_factors = 3,
                            n_init = 2, n_random = 5),
    "at least twenty"
  )
})

test_that("noise degrades recovery towards the random bound", {
  clean <- paired_data(effect_size = 3)
  weak <- paired_data(effect_size = 0)
  rc <- chorale_destroy_pairing(clean$a, clean$b, clean$design, n_factors = 3,
                               n_init = 5, n_random = 50, seed = 1)
  rw <- chorale_destroy_pairing(weak$a, weak$b, weak$design, n_factors = 3,
                               n_init = 5, n_random = 50, seed = 1)
  # With no design signal there is nothing for the estimator to match on, so it
  # cannot place as well as it does when the signal is present.
  if (!is.finite(rw$summary$placement_between_bounds)) {
    expect_equal(rw$summary$n_recovered_pairs, 0L)
  } else {
    expect_gte(rc$summary$placement_between_bounds,
               rw$summary$placement_between_bounds)
  }
})
