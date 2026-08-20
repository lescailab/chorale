sim_collection <- function(n_features = 120L, seed = 1L) {
  sim <- chorale_simulate(n_modalities = 2L, n_features = n_features,
                          n_shared_factors = 2L, n_private_factors = 1L,
                          n_strains = 4L, n_per_cell = 2L, seed = seed)
  Map(function(assay, cd) chorale_load(assay, col_data = cd),
      sim$modalities, sim$col_data)
}

test_that("detectability reports a usable factor count", {
  containers <- sim_collection()
  xs <- lapply(containers, function(x) scale(t(SummarizedExperiment::assay(x))))
  d <- chorale_gate_detectability(xs, n_perm = 10L, max_factors = 6L, seed = 1L)
  expect_setequal(d$modality, names(containers))
  expect_true(all(d$n_factors >= 2L))
  expect_true(all(d$n_factors <= 6L))
})

test_that("anchor richness coarsens down to phenotype", {
  containers <- sim_collection()
  designs <- lapply(containers, function(x) {
    as.data.frame(SummarizedExperiment::colData(x))
  })
  a <- chorale_gate_anchors(designs)
  expect_true(all(grepl("phenotype", a$coarsening)))
  expect_true("phenotype" %in% a$coarsening)
  expect_true(all(a$n_strata_shared <= a$n_strata))
})

test_that("R-native normality statistic responds to heavy tails", {
  set.seed(1)
  gaussian <- chorale_ad_normal(stats::rnorm(500))
  heavy <- chorale_ad_normal(stats::rt(500, df = 2))
  expect_gt(heavy, gaussian)
})

test_that("modality distribution comparison is a diagnostic", {
  set.seed(1)
  s <- scale(matrix(stats::rt(80 * 3, df = 3), nrow = 80))
  d <- chorale_gate_modality_difference(list(a = s, b = s))
  expect_equal(nrow(d), 1L)
  expect_identical(d$role, "diagnostic; does not gate matching")
  expect_identical(d$verdict, "distribution difference not detected")
  expect_false(d$gates_matching)
})

test_that("component estimator is deterministic", {
  skip_if_not_installed("fastICA")
  set.seed(1)
  x <- scale(matrix(stats::rt(60 * 40, df = 3), nrow = 60))
  ica_fn <- chorale_gate_ica(n_init = 2L, consensus = FALSE)
  expect_equal(ica_fn(x, 3L, 1L), ica_fn(x, 3L, 1L))
})

test_that("production gates report design rank and factor stability", {
  containers <- sim_collection(n_features = 80L)
  g <- chorale_gates(containers, control = chorale_control(n_init = 2L),
                     n_surrogate = 2L, n_perm = 5L)
  expect_true(all(g$design_estimability$phenotype_estimable))
  expect_true(all(g$design_estimability$full_rank))
  expect_equal(nrow(g$factor_stability), length(containers))
  expect_true(all(is.finite(g$factor_stability$mean_subspace_agreement)))
})
