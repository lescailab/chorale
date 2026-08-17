test_that("chorale_null signals not-yet-implemented pending chorale_fit()", {
  expect_error(
    chorale_null(fit = NULL),
    class = "chorale_not_implemented"
  )
})
