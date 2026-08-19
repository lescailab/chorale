make_profile <- function(seed = 1, ...) {
  sim <- chorale_simulate(n_modalities = 2, n_features = 40,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 3, n_per_cell = 2, seed = seed)
  chorale_data_profile(sim$modalities[[1]], sim$col_data[[1]],
                       covariates = c("phenotype", "sex"), ...)
}

test_that("a profile records the properties a simulation must reproduce", {
  p <- make_profile()
  expect_s3_class(p, "chorale_data_profile")
  expect_equal(p$n_features, 40)
  expect_equal(nrow(p$quantiles), 40)
  expect_equal(ncol(p$quantiles), length(p$probs))
  expect_equal(nrow(p$feature), 40)
  expect_true(all(c("abundance", "spread", "missing", "kurtosis") %in%
                    names(p$feature)))
  expect_true(length(p$eigenvalues) >= 2)
  expect_true(abs(sum(p$eigenvalues)) <= 1 + 1e-8)
})

test_that("a profile holds no sample identifier and no covariate it was not given", {
  p <- make_profile()
  flat <- unlist(lapply(p[setdiff(names(p), "feature")], as.character))
  expect_false(any(grepl("^modality1_sample", flat)))
  expect_equal(sort(setdiff(names(p$design_cells), "n")),
               c("phenotype", "sex"))
})

test_that("design cells record only the combinations that were populated", {
  d <- data.frame(sample_id = paste0("s", 1:6),
                  phenotype = c("case", "case", "case", "control", "control", "control"),
                  sex = c("F", "F", "M", "F", "F", "F"),
                  stringsAsFactors = FALSE)
  x <- matrix(rnorm(60), nrow = 10, dimnames = list(NULL, d$sample_id))
  p <- chorale_data_profile(x, d, covariates = c("phenotype", "sex"))
  expect_equal(nrow(p$design_cells), 3)
  expect_equal(sum(p$design_cells$n), 6)
})

test_that("missingness is recorded rather than imputed", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 3, n_per_cell = 2, seed = 4)
  x <- sim$modalities[[1]]
  x[1, 1:5] <- NA
  p <- chorale_data_profile(x, sim$col_data[[1]], covariates = "phenotype")
  expect_equal(p$feature$missing[1], 5 / ncol(x))
  expect_true(all(is.finite(p$quantiles[1, ])))
})

test_that("a profile needs at least two features and two samples", {
  expect_error(chorale_data_profile(matrix(1:4, nrow = 1)), "at least two")
})

test_that("agreement puts a simulation beside the profile it was drawn from", {
  p <- make_profile()
  sim <- chorale_simulate(n_modalities = 2, n_features = 40,
                          n_shared_factors = 2, n_private_factors = 1,
                          profile = p, seed = 5)
  a <- chorale_profile_agreement(chorale_data_profile(sim$modalities[[1]]), p)
  expect_true(all(c("property", "simulated", "reference", "discrepancy") %in%
                    names(a)))
  quantiles <- a[a$property %in% c("marginal_q10", "marginal_median",
                                   "marginal_q90"), ]
  # The marginals are imposed feature by feature, so they agree closely.
  expect_true(all(quantiles$discrepancy < 0.5))
})

test_that("agreement refuses anything that is not a profile", {
  p <- make_profile()
  expect_error(chorale_profile_agreement(p, list()), "chorale_data_profile")
})
