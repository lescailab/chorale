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

count_collection <- function(n_features = 60L, n_samples = 24L, seed = 1L) {
  set.seed(seed)
  mk <- function(tag) {
    mu <- rep(2^stats::runif(n_features, 1, 12), times = n_samples)
    a <- matrix(stats::rpois(n_features * n_samples, lambda = mu),
                nrow = n_features,
                dimnames = list(paste0("g", seq_len(n_features)),
                                paste0(tag, seq_len(n_samples))))
    d <- data.frame(sample_id = colnames(a),
                    phenotype = rep(c("control", "case"), length.out = n_samples),
                    stringsAsFactors = FALSE)
    chorale_load(a, d)
  }
  list(a = mk("a"), b = mk("b"))
}

test_that("the gates read a modality on the scale the encoder consumes", {
  skip_if_not_installed("fastICA")
  containers <- count_collection()
  g <- chorale_gates(containers, control = chorale_control(n_init = 2L),
                     n_surrogate = 2L, n_perm = 5L)
  # Counts spanning orders of magnitude take the variance-stabilising
  # transform, which is what the encoder chooses for the same matrices.
  expect_equal(g$transform$modality, names(containers))
  expect_true(all(g$transform$transform == "vst"))

  sets <- list(one = paste0("g", 1:20), two = paste0("g", 21:45))
  cc <- chorale_concepts(containers, sets, min_features = 5L)
  enc <- chorale_encode(containers, cc, n_free = 1L, n_init = 2L)
  encoder_scale <- vapply(enc$encodings, `[[`, character(1), "transform")
  expect_equal(g$transform$transform, unname(encoder_scale[g$transform$modality]))

  # The matrix each diagnostic is read on is the matrix the encoder builds.
  for (m in names(containers)) {
    gate_x <- chorale_analysis_matrix(
      SummarizedExperiment::assay(containers[[m]]))$matrix
    expect_equal(unname(gate_x),
                 unname(enc$encodings[[m]]$analysis_matrix))
  }
})

test_that("a missing value reaches the gates as it reaches the encoder", {
  skip_if_not_installed("fastICA")
  containers <- count_collection(n_features = 40L, n_samples = 20L, seed = 3L)
  # A feature measured in only part of the cohort, and one measured nowhere.
  for (m in names(containers)) {
    a <- SummarizedExperiment::assay(containers[[m]])
    a[1L, 1:5] <- NA_real_
    a[2L, ] <- NA_real_
    SummarizedExperiment::assay(containers[[m]]) <- a
  }

  sets <- list(one = paste0("g", 1:20), two = paste0("g", 21:35))
  cc <- chorale_concepts(containers, sets, min_features = 5L)
  enc <- chorale_encode(containers, cc, n_free = 1L, n_init = 2L)

  for (m in names(containers)) {
    gate_x <- chorale_analysis_matrix(
      SummarizedExperiment::assay(containers[[m]]))$matrix
    expect_equal(unname(gate_x), unname(enc$encodings[[m]]$analysis_matrix))
    # Substitution happens after standardisation, so the filled entry is the
    # standardised mean and the measured values keep their own mean and spread.
    expect_true(all(is.finite(gate_x)))
    expect_equal(unname(gate_x[, 2L]), rep(0, nrow(gate_x)))
    measured <- gate_x[6:nrow(gate_x), 1L]
    expect_equal(mean(measured), mean(scale(
      SummarizedExperiment::assay(containers[[m]])[1L, 6:ncol(
        SummarizedExperiment::assay(containers[[m]]))])[, 1]),
      tolerance = 1e-8)
  }

  expect_silent(chorale_gates(containers, control = chorale_control(n_init = 2L),
                              n_surrogate = 2L, n_perm = 5L))
})

test_that("gates refuse designs that name different modalities", {
  containers <- count_collection(n_features = 30L, n_samples = 12L)
  designs <- lapply(containers, function(x) {
    as.data.frame(SummarizedExperiment::colData(x))
  })
  names(designs) <- c("a", "wrong")
  expect_error(
    chorale_gates(containers, designs = designs,
                  control = chorale_control(n_init = 2L),
                  n_surrogate = 2L, n_perm = 5L),
    "names different modalities")

  # An unnamed list of the right length takes the collection's names in order.
  bare <- unname(designs)
  expect_error(
    chorale_gates(containers, designs = bare[1],
                  control = chorale_control(n_init = 2L),
                  n_surrogate = 2L, n_perm = 5L),
    "must be named")
})

test_that("the gate transform can be overridden per modality", {
  skip_if_not_installed("fastICA")
  containers <- count_collection()
  g <- chorale_gates(containers, control = chorale_control(n_init = 2L),
                     transform = c(a = "log", b = "none"),
                     n_surrogate = 2L, n_perm = 5L)
  expect_equal(g$transform$transform, c("log", "none"))
})
