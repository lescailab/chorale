test_that("chorale_bound signals not-yet-implemented pending chorale_fit()", {
  expect_error(
    chorale_bound(fit = NULL),
    class = "chorale_not_implemented"
  )
})
