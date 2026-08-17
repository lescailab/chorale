test_that("chorale_fit signals not-yet-implemented pending Stage 3 gates", {
  expect_error(
    chorale_fit(containers = list(), n_factors = 5, gene_sets = list()),
    class = "chorale_not_implemented"
  )
})
