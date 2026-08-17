test_that("chorale_map resolves unambiguous gene symbols", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  out <- chorale_map(c("Bdnf", "Trem2"), from = "SYMBOL")
  expect_true(all(c("Bdnf", "Trem2") %in% out$id))
  expect_true(all(out$weight > 0 & out$weight <= 1))
})

test_that("chorale_map drops unmatched identifiers rather than erroring", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  out <- chorale_map(c("Bdnf", "not_a_real_gene_xyz"), from = "SYMBOL")
  expect_false("not_a_real_gene_xyz" %in% out$id)
})

test_that("chorale_map assigns fractional weights that sum to 1 per identifier", {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  out <- chorale_map(c("Bdnf", "Trem2"), from = "SYMBOL")
  totals <- stats::aggregate(weight ~ id, data = out, FUN = sum)
  expect_true(all(abs(totals$weight - 1) < 1e-8))
})

test_that("chorale_map rejects non-character input", {
  expect_error(chorale_map(1:3), class = "rlang_error")
})
