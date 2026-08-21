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
  expect_true(all(c("modality", "subspace_agreement", "subspace_min") %in%
                    colnames(n$stability)))
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
  # A shared feature space is what lets the shuffle run at all. What it then
  # reports is a description, since the pooled collections are not comparable
  # with the statistic they would be read against.
  expect_true(is.finite(n$modality_null$agreement))
  expect_gt(n$modality_null$n_shuffles, 0)
  expect_true(is.na(n$modality_null$p_value))
})

test_that("stability records every initialisation", {
  f <- make_fit()
  n <- chorale_null(f$fit, f$containers, n_permutations = 2, n_init = 2)
  expect_true(all(n$stability$n_init == 2))
  expect_true(all(n$stability$n_failed == 0))
})

test_that("factor stability measures the factors, not the objective", {
  # A stable objective can coexist with rotated or permuted factors, so
  # stability is reported as how well the recovered factors match across starts.
  f <- make_fit()
  n <- chorale_null(f$fit, f$containers, n_permutations = 2, n_init = 4)
  expect_true(all(n$stability$subspace_agreement >= -1 &
                    n$stability$subspace_agreement <= 1))
  # On clean data with a strong effect the same factors should recover whatever
  # the start.
  expect_gt(min(n$stability$subspace_min), 0.5)
})

test_that("the modality shuffle describes the shuffles and reports no p-value", {
  # Two modalities on one shared feature space, so the shuffle is defined and
  # returns a calibrated p-value rather than a single draw.
  set.seed(1)
  ids <- paste0("g", 1:60)
  mk <- function(tag) {
    m <- matrix(stats::rnorm(60 * 40), nrow = 60, dimnames = list(ids, NULL))
    d <- data.frame(sample_id = paste0(tag, 1:40),
                    cohort = "synthetic", modality = tag, strain = "group_1",
                    phenotype = rep(c("control", "case"), each = 20),
                    age_months = 6, sex = rep(c("F", "M"), 20),
                    region = "sim", batch = "b1",
                    stringsAsFactors = FALSE)
    colnames(m) <- d$sample_id
    chorale_load(m, d)
  }
  containers <- list(a = mk("a"), b = mk("b"))
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
  n <- chorale_null(fit, containers, n_permutations = 2, n_init = 2)
  skip_if(!is.finite(n$modality_null$agreement))
  # The shuffled collections run on the common features alone and do not keep
  # each modality's sample composition, so they describe rather than test.
  expect_false(n$modality_null$applicable)
  expect_true(is.na(n$modality_null$p_value))
  expect_true(is.finite(n$modality_null$agreement))
  expect_gt(n$modality_null$n_shuffles, 1)
  expect_match(n$modality_null$reason, "descriptive only")
})
