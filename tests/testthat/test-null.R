make_fit <- function(effect_size = 3, seed = 1) {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 5, n_per_cell = 3,
                          effect_size = effect_size, seed = seed)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  list(containers = containers,
       fit = chorale_fit(containers, n_factors = c(3, 3), n_init = 2, seed = seed))
}

test_that("the controls run and report all three", {
  f <- make_fit()
  n <- chorale_null(f$fit, f$containers, n_permutations = 3, n_init = 2)
  expect_s3_class(n, "chorale_null")
  expect_length(n$phenotype_null, 3)
  expect_true(all(c("modality", "objective_cv") %in% colnames(n$stability)))
  expect_output(print(n), "chorale_null")
})

test_that("the permutation p-value is a valid probability", {
  f <- make_fit()
  n <- chorale_null(f$fit, f$containers, n_permutations = 5, n_init = 2)
  skip_if(is.na(n$p_phenotype))
  expect_gte(n$p_phenotype, 1 / 6)
  expect_lte(n$p_phenotype, 1)
})

test_that("the modality shuffle reports inapplicability rather than passing", {
  # The simulated modalities carry their own feature namespaces, so pooling
  # them is undefined and must be declared, not silently treated as a pass.
  f <- make_fit()
  n <- chorale_null(f$fit, f$containers, n_permutations = 2, n_init = 2)
  expect_false(isTRUE(n$modality_null$applicable))
  expect_true(nzchar(n$modality_null$reason))
})

test_that("the modality shuffle runs where a feature space is shared", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  shared_names <- paste0("gene", seq_len(nrow(sim$modalities[[1]])))
  mats <- lapply(sim$modalities, function(m) {
    rownames(m) <- shared_names
    m
  })
  containers <- Map(chorale_load, mats, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
  n <- chorale_null(fit, containers, n_permutations = 2, n_init = 2)
  expect_true(n$modality_null$applicable)
})

test_that("stability records every initialisation", {
  f <- make_fit()
  n <- chorale_null(f$fit, f$containers, n_permutations = 2, n_init = 2)
  expect_true(all(n$stability$n_init == 2))
  expect_true(all(n$stability$n_failed == 0))
})
