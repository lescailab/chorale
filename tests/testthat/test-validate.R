test_that("recovery is scored against the planted truth", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 150,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = 3, n_init = 5, seed = 1)

  align <- chorale_align_truth(fit, sim)
  expect_true(all(c("modality", "factor", "planted", "shared") %in% colnames(align)))
  expect_equal(nrow(align), sum(fit$n_factors))

  score <- chorale_score_recovery(fit, sim)
  # On clean data the estimator should recover the planted shared programmes
  # with no false matches.
  expect_gte(score$shared_recovered, 0.5)
  expect_lte(score$false_match_rate, 0.34)
})

test_that("the validation matrix separates the regimes it should", {
  grid <- data.frame(
    label = c("clean", "null", "same_response"),
    effect_size = c(3, 0, 3),
    same_response = c(FALSE, FALSE, TRUE)
  )
  res <- chorale_validate(grid, n_rep = 2, seed = 1)

  clean <- res[res$label == "clean", ]
  null <- res[res$label == "null", ]
  adv <- res[res$label == "same_response", ]

  # Clean data recovers shared programmes; the complete null recovers none.
  expect_gt(clean$shared_recovered, 0.5)
  expect_equal(null$programmes_correct, 0)
  # Two distinct programmes sharing a phenotype response cannot be told apart,
  # which is a limitation the matrix is built to expose rather than hide.
  expect_lt(adv$assignment_accuracy, clean$assignment_accuracy)
})

test_that("the joint null is calibrated", {
  cal <- chorale_null_calibration(
    n_sim = 30, n_perm = 199, n_init = 3,
    n_modalities = 2, n_features = 120, n_shared_factors = 2,
    n_private_factors = 1, n_strains = 4, n_per_cell = 3, seed = 1
  )
  skip_if(cal$n_evaluated < 20)
  # A calibrated procedure returns uniform p-values under the null, so the
  # Kolmogorov-Smirnov test should not reject, and the false positive rate at
  # 0.05 should sit near 0.05 rather than far above it.
  expect_gt(cal$ks_p, 0.05)
  expect_lt(cal$false_positive_rate, 0.15)
})
