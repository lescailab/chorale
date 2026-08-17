test_that("chorale_report signals not-yet-implemented pending upstream fits", {
  expect_error(
    chorale_report(fit = NULL, bound = NULL, null = NULL, path = withr::local_tempdir()),
    class = "chorale_not_implemented"
  )
})
