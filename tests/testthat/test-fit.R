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
  # condition checkable on the recovered loadings. The planted sources are
  # non-symmetric, as the identification results require and unlike the earlier
  # symmetric draws, which is a harder recovery than a symmetric source; random
  # alignment would sit near 0.1.
  expect_gt(mean(matched), 0.85)
})

test_that("matching finds shared factors when the design carries them", {
  out <- sim_fit(effect_size = 3, n_modalities = 3, n_per_cell = 4,
                 n_strains = 8)
  sig <- out$fit$matches[out$fit$matches$significant, , drop = FALSE]
  expect_gt(nrow(sig), 0)
  expect_true(all(sig$p_value < 0.05))
  expect_true(all(sig$statistic > 0))
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

test_that("matching works from the phenotype alone", {
  # Most deposited datasets share little beyond the case/control label, so a
  # single covariate must suffice.
  sim <- chorale_simulate(n_modalities = 2, n_features = 150,
                          n_shared_factors = 3, n_private_factors = 2,
                          n_strains = 4, n_per_cell = 4, effect_size = 2,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(5, 5), n_init = 3,
                     strata_keys = "phenotype", seed = 1)
  expect_equal(unique(fit$matches$n_shared_covariates), 1L)
  expect_gt(sum(fit$matches$significant), 0)
})

test_that("power grows with sample size rather than with strata count", {
  weak <- function(n_per_cell) {
    sim <- chorale_simulate(n_modalities = 2, n_features = 150,
                            n_shared_factors = 3, n_private_factors = 2,
                            n_strains = 4, n_per_cell = n_per_cell,
                            effect_size = 0.6, seed = 1)
    containers <- Map(chorale_load, sim$modalities, sim$col_data)
    m <- chorale_fit(containers, n_factors = c(5, 5), n_init = 3,
                     strata_keys = "phenotype", seed = 1)$matches
    min(m$p_value)
  }
  expect_lte(weak(8), weak(1))
})

test_that("a modality without a phenotype contrast is refused", {
  # The estimand is a case/control contrast, so a modality that cannot express
  # one is excluded rather than matched on distribution shape alone.
  sim <- chorale_simulate(n_modalities = 2, n_features = 150,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 3, n_per_cell = 3, effect_size = 2,
                          seed = 1)
  blank <- lapply(sim$col_data, function(d) {
    d$phenotype <- "unknown"
    d
  })
  containers <- Map(chorale_load, sim$modalities, blank)
  expect_error(chorale_fit(containers, n_factors = c(3, 3), n_init = 2),
               class = "chorale_missing_phenotype")
})

test_that("a constant phenotype is refused", {
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
  # A constant phenotype carries no contrast, so it is refused.
  expect_error(chorale_match(fits, designs, n_perm = 20),
               class = "chorale_missing_phenotype")
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

test_that("programme evidence is joint, not the best pairwise link", {
  # The plan requires false-discovery control that improves with the number of
  # modalities, which the strongest pairwise p-value cannot deliver.
  sim <- chorale_simulate(n_modalities = 3, n_features = 150,
                          n_shared_factors = 3, n_private_factors = 2,
                          n_strains = 6, n_per_cell = 4, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = rep(5, 3), n_init = 3, seed = 1)
  pg <- chorale_joint_evidence(fit)
  skip_if(nrow(pg) == 0)
  expect_true(all(c("joint_statistic", "joint_p") %in% colnames(pg)))
  ok <- !is.na(pg$joint_p)
  skip_if(!any(ok))
  floor <- 1 / (1 + length(fit$joint_null))
  expect_true(all(pg$joint_p[ok] >= floor))
  expect_true(all(pg$joint_p[ok] <= 1))
  # A programme spanning more modalities carries a larger joint statistic,
  # since agreement is required simultaneously rather than pair by pair.
  u <- unique(pg[ok, c("programme", "n_modalities", "joint_statistic")])
  if (length(unique(u$n_modalities)) > 1) {
    expect_gt(
      stats::median(u$joint_statistic[u$n_modalities == max(u$n_modalities)]),
      stats::median(u$joint_statistic[u$n_modalities == min(u$n_modalities)])
    )
  }
})

test_that("the assignment is solved once, so it agrees around every cycle", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
  pg <- fit$programmes
  skip_if(nrow(pg) == 0)

  # A factor belongs to at most one programme, and a programme carries at most
  # one factor per modality. Chaining pairwise decisions cannot guarantee this.
  key <- paste(pg$modality, pg$factor)
  expect_equal(anyDuplicated(key), 0L)
  expect_equal(anyDuplicated(pg[, c("programme", "modality")]), 0L)

  # Every implied pair is a pair of members of the same programme, so no link
  # crosses two programmes.
  m <- fit$matches
  skip_if(nrow(m) == 0)
  member <- stats::setNames(pg$programme, paste(pg$modality, pg$factor))
  expect_equal(unname(member[paste(m$modality_a, m$factor_a)]), m$programme)
  expect_equal(unname(member[paste(m$modality_b, m$factor_b)]), m$programme)
})

test_that("the joint assignment does not depend on the order of modalities", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  a <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
  b <- chorale_fit(containers[c(3, 1, 2)], n_factors = c(3, 3, 3), n_init = 2)

  grouping <- function(fit) {
    pg <- fit$programmes[fit$programmes$supported, , drop = FALSE]
    if (nrow(pg) == 0) return(character(0))
    sort(vapply(split(paste(pg$modality, pg$factor), pg$programme),
                function(v) paste(sort(v), collapse = " + "), character(1),
                USE.NAMES = FALSE))
  }
  expect_equal(grouping(a), grouping(b))
})

test_that("a design profile expands a multilevel covariate into signed terms", {
  scores <- matrix(c(1:9, 9:1), ncol = 2,
                   dimnames = list(paste0("s", 1:9), c("factor_1", "factor_2")))
  design <- data.frame(sample_id = paste0("s", 1:9),
                       region = rep(c("cortex", "hippocampus", "striatum"), 3),
                       stringsAsFactors = FALSE)
  p <- chorale_design_profile(scores, design, "region")
  # One signed contrast per level beyond the reference, never an unsigned
  # summary that cannot be compared in direction with the other terms.
  expect_equal(colnames(p), c("region=hippocampus", "region=striatum"))
  expect_true(all(is.finite(p)))
  expect_lt(p[1, "region=hippocampus"] * p[2, "region=hippocampus"], 0)
})

test_that("leaving a modality out reports what it contributed", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 2)
  loo <- chorale_leave_one_out(fit)
  skip_if(nrow(loo) == 0)
  expect_true(all(c("programme", "dropped", "joint_statistic", "joint_p",
                    "delta") %in% colnames(loo)))
  expect_true(all(loo$n_modalities == 2))
  expect_true(all(loo$joint_p >= 1 / (1 + length(fit$joint_null))))
})

test_that("pure features can be required as a gate on programmes", {
  sim <- chorale_simulate(n_modalities = 3, n_features = 150,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 3)

  pg <- chorale_programmes(fit)
  expect_true(all(c("pure_features", "all_pure") %in% colnames(pg)))

  # Strip the pure-feature condition from every factor. The gate must then
  # return nothing, while the ungated call is unaffected.
  for (m in fit$modalities) {
    fit$fits[[m]]$pure_feature_condition[] <- FALSE
  }
  fit$programmes$all_pure <- FALSE
  fit$programmes$pure_features <- FALSE
  expect_gt(nrow(chorale_programmes(fit, require_pure_features = FALSE)), 0)
  expect_equal(nrow(chorale_programmes(fit, require_pure_features = TRUE)), 0)
})
