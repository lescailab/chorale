sim_fit <- function(effect_size = 3, n_modalities = 2, seed = 1, n_init = 3,
                    n_per_cell = 3, n_strains = 6) {
  sim <- chorale_simulate(
    n_modalities = n_modalities, n_features = 150, n_shared_factors = 3,
    n_private_factors = 2, n_strains = n_strains, n_per_cell = n_per_cell,
    effect_size = effect_size, seed = seed
  )
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  list(sim = sim,
       fit = chorale_fit(containers, n_factors = rep(5, n_modalities),
                         n_init = n_init, seed = seed))
}

test_that("chorale_fit returns a fit over every modality", {
  out <- sim_fit()
  expect_s3_class(out$fit, "chorale_fit")
  expect_length(out$fit$fits, 2)
  for (m in out$fit$modalities) {
    expect_equal(ncol(out$fit$fits[[m]]$scores), 5)
    expect_equal(ncol(out$fit$fits[[m]]$loadings), 5)
    expect_equal(nrow(out$fit$fits[[m]]$stability), 3)
  }
  expect_output(print(out$fit), "chorale_fit")
})

test_that("recovered components track the planted latent scores", {
  sim <- chorale_simulate(
    n_modalities = 2, n_features = 200, n_shared_factors = 3,
    n_private_factors = 2, n_strains = 8, n_per_cell = 4, seed = 1
  )
  x <- scale(t(sim$modalities[[1]]))
  ica <- chorale_ica(x, 5, n_init = 3)
  truth <- cbind(sim$truth$scores[[1]]$shared, sim$truth$scores[[1]]$private)
  agreement <- abs(stats::cor(ica$scores, truth))
  assignment <- clue::solve_LSAP(agreement, maximum = TRUE)
  matched <- agreement[cbind(seq_len(nrow(agreement)), assignment)]
  # A rotation this close to the truth is what makes the pure-feature
  # condition checkable on the recovered loadings.
  expect_gt(mean(matched), 0.9)
})

test_that("matching finds shared factors when the design carries them", {
  out <- sim_fit(effect_size = 3, n_modalities = 3, n_per_cell = 4,
                 n_strains = 8)
  sig <- out$fit$matches[out$fit$matches$significant, , drop = FALSE]
  expect_gt(nrow(sig), 0)
  expect_true(all(sig$p_value < 0.05))
  expect_true(all(abs(sig$anchor_agreement) > 0.5))
})

test_that("matching is calibrated when the shared state ignores the design", {
  # With no design effect there is nothing for an anchor to corroborate, so
  # matches should be rare rather than guaranteed.
  found <- vapply(1:3, function(s) {
    m <- sim_fit(effect_size = 0, seed = s, n_per_cell = 4,
                 n_strains = 8)$fit$matches
    if (nrow(m) == 0) 0 else sum(m$significant)
  }, numeric(1))
  expect_lt(mean(found), 2)
})

test_that("matching declines without enough shared strata", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 80, seed = 1,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 2, n_per_cell = 2)
  fits <- list(a = list(scores = matrix(rnorm(20), 10, 2)),
               b = list(scores = matrix(rnorm(20), 10, 2)))
  designs <- list(
    a = data.frame(sample_id = paste0("s", 1:10), phenotype = "WT",
                   age_bin = "6mo", sex = "F"),
    b = data.frame(sample_id = paste0("t", 1:10), phenotype = "WT",
                   age_bin = "6mo", sex = "F")
  )
  rownames(fits$a$scores) <- designs$a$sample_id
  rownames(fits$b$scores) <- designs$b$sample_id
  expect_equal(nrow(chorale_match(fits, designs)), 0)
})

test_that("pure features are found on loadings built to be pure", {
  set.seed(1)
  l <- matrix(rnorm(90), nrow = 30,
              dimnames = list(paste0("f", 1:30), paste0("factor_", 1:3)))
  l[1:3, 2:3] <- 0
  l[4:6, c(1, 3)] <- 0
  mk <- chorale_markers(l)
  expect_true(all(paste0("f", 1:3) %in% mk$markers$factor_1))
  expect_true(all(paste0("f", 4:6) %in% mk$markers$factor_2))
  expect_true(mk$pure_feature_condition[["factor_1"]])
})

test_that("the purity margin is reported even where no feature is pure", {
  set.seed(1)
  l <- matrix(rnorm(90), nrow = 30,
              dimnames = list(paste0("f", 1:30), paste0("factor_", 1:3)))
  mk <- chorale_markers(l, purity_ratio = 0.001)
  expect_true(all(lengths(mk$markers) == 0))
  expect_true(all(lengths(mk$best_candidates) > 0))
  expect_true(all(is.finite(mk$purity_margin)))
  expect_false(any(mk$pure_feature_condition))
})

test_that("the pathway constraint leaves marker loadings untouched", {
  set.seed(1)
  l <- matrix(rnorm(40), nrow = 10,
              dimnames = list(paste0("f", 1:10), paste0("factor_", 1:4)))
  p <- matrix(rbinom(30, 1, 0.5), nrow = 10,
              dimnames = list(paste0("f", 1:10), paste0("set", 1:3)))
  markers <- list(c("f1", "f2"), character(), character(), character())
  out <- chorale_constrain(l, p, markers)
  expect_equal(out$loadings[c("f1", "f2"), 1], l[c("f1", "f2"), 1])
  # Non-marker loadings are shrunk towards the sets, so they move.
  expect_false(isTRUE(all.equal(out$loadings[3:10, 1], l[3:10, 1])))
  expect_equal(dim(out$set_weights), c(3, 4))
})

test_that("chorale_fit rejects a single modality", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  expect_error(chorale_fit(containers[1], n_factors = 2), class = "rlang_error")
})
