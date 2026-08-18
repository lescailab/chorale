test_that("chorale_load builds a valid container from simulated data", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  se <- chorale_load(sim$modalities[[1]], sim$col_data[[1]])

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), dim(sim$modalities[[1]]))
  expect_true(all(chorale_required_col_data() %in% colnames(SummarizedExperiment::colData(se))))
})

test_that("chorale_load rejects col_data missing required columns", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  bad_col_data <- sim$col_data[[1]][, setdiff(names(sim$col_data[[1]]), "strain")]
  expect_error(chorale_load(sim$modalities[[1]], bad_col_data), class = "rlang_error")
})

test_that("chorale_load rejects mismatched sample order", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  shuffled <- sim$col_data[[1]][rev(seq_len(nrow(sim$col_data[[1]]))), ]
  expect_error(chorale_load(sim$modalities[[1]], shuffled), class = "rlang_error")
})

test_that("chorale_load rejects a non-matrix assay", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  expect_error(
    chorale_load(as.data.frame(sim$modalities[[1]]), sim$col_data[[1]]),
    class = "rlang_error"
  )
})

test_that("blank and placeholder strings are read as missing, not as levels", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3,
                          n_per_cell = 1, seed = 1)
  cd <- sim$col_data[[1]]
  # A design read from a delimited file carries "" where a value was absent.
  cd$phenotype[3:nrow(cd)] <- ""
  cd$sex[1:2] <- "unknown"
  se <- chorale_load(sim$modalities[[1]], cd)
  out <- SummarizedExperiment::colData(se)
  expect_true(all(is.na(out$phenotype[3:nrow(cd)])))
  expect_true(all(is.na(out$sex[1:2])))
  # One real level plus blanks is not a two-level contrast.
  expect_equal(length(unique(stats::na.omit(out$phenotype))), 1L)
})
