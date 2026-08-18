sim_fit3 <- function(seed = 1, effect_size = 3) {
  sim <- chorale_simulate(n_modalities = 3, n_features = 150,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4,
                          effect_size = effect_size, seed = seed)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  list(sim = sim, containers = containers,
       fit = chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 3,
                         seed = seed))
}

test_that("false discovery is controlled within each level of the search", {
  f <- sim_fit3()
  q <- chorale_fdr(f$fit)
  skip_if(nrow(q) == 0)
  expect_true(all(c("level", "object", "p_value", "q_value", "significant")
                  %in% colnames(q)))
  expect_true(all(q$level %in% c("factor", "programme", "pathway")))
  # A q-value is never below its p-value, and both are probabilities.
  expect_true(all(q$q_value >= q$p_value - 1e-9))
  expect_true(all(q$q_value <= 1 & q$q_value >= 0))
  # Levels are corrected separately, so a level's q-values depend only on that
  # level's p-values.
  for (lv in unique(q$level)) {
    i <- q$level == lv
    expect_equal(q$q_value[i], stats::p.adjust(q$p_value[i], method = "BH"))
  }
})

test_that("added value compares a programme with its best single modality", {
  f <- sim_fit3()
  av <- chorale_added_value(f$fit)
  skip_if(nrow(av) == 0)
  expect_true(all(c("joint_p", "best_single_p", "margin", "needs_multiple")
                  %in% colnames(av)))
  # The margin is the gap on the negative log scale, so it is positive exactly
  # when the programme is stronger than anything one modality carried.
  pos <- is.finite(av$margin) & is.finite(av$best_single_p)
  expect_equal(av$margin[pos] > 0, av$joint_p[pos] < av$best_single_p[pos])
})

test_that("a programme resting on one modality is not called multi-modal", {
  f <- sim_fit3()
  fit <- f$fit
  # Removing a modality leaves the programme stronger, so it did not need it.
  if (!is.null(fit$leave_one_out) && nrow(fit$leave_one_out) > 0) {
    fit$leave_one_out$delta <- abs(fit$leave_one_out$delta) + 1
    av <- chorale_added_value(fit)
    skip_if(nrow(av) == 0)
    expect_false(any(av$needs_multiple))
  }
})

test_that("specificity reruns the pipeline with a nuisance anchor", {
  f <- sim_fit3()
  sp <- chorale_specificity(f$fit, f$containers, covariates = "sex", n_init = 2)
  expect_equal(nrow(sp), 1L)
  expect_true(all(c("anchor", "joint_statistic", "observed_phenotype",
                    "stronger_than_phenotype") %in% colnames(sp)))
  # A covariate that is absent or constant is declared, not silently passed.
  sp2 <- chorale_specificity(f$fit, f$containers, covariates = "not_a_column",
                             n_init = 2)
  expect_true(is.na(sp2$stronger_than_phenotype))
  expect_match(sp2$reason, "absent or constant")
})

test_that("cohort overlap and common support describe the shared population", {
  f <- sim_fit3()
  ov <- chorale_cohort_overlap(f$fit)
  skip_if(nrow(ov) == 0)
  expect_true(all(ov$max_total_variation >= 0 & ov$max_total_variation <= 1))
  # The simulated cohorts realise the design identically, so they overlap
  # completely and no claim is restricted.
  expect_true(all(ov$comparable))
  cs <- chorale_common_support(f$fit)
  expect_gt(length(cs$cells), 0)
  expect_true(all(cs$share > 0))
})

test_that("cohorts that do not overlap are reported as not comparable", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 100,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  # One cohort holds only females, the other only males: no common support on
  # sex, and a claim about sex cannot be made from the pair.
  cd <- sim$col_data
  cd[[1]]$sex <- "F"
  cd[[2]]$sex <- "M"
  containers <- Map(chorale_load, sim$modalities, cd)
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
  ov <- chorale_cohort_overlap(fit)
  expect_false("sex" %in% ov$covariate)
})
