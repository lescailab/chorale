# How a covariate is compared across modalities, and the invariant that a
# covariate carried in a signature always has terms to build.

test_that("a covariate numeric in every modality is continuous whatever its spread", {
  # Two distinct values in one modality and many in another is an ordinary
  # two-timepoint design. It must not be sent down the categorical path.
  vals <- list(a = c(6, 6.2, 15.1, 15.2, 20), b = c(6, 14, 6, 14))
  expect_identical(chorale_covariate_kind(vals), "continuous")

  designs <- list(
    a = data.frame(sample_id = paste0("a", 1:5),
                   phenotype = c("case", "control", "case", "control", "case"),
                   age_months = c(6, 6.2, 15.1, 15.2, 20)),
    b = data.frame(sample_id = paste0("b", 1:4),
                   phenotype = c("case", "control", "case", "control"),
                   age_months = c(6, 14, 6, 14))
  )
  levels <- chorale_profile_levels(designs, c("phenotype", "age_months"))
  expect_true("age_months" %in% names(levels))
  expect_true(is.na(levels$age_months))
})

test_that("a numeric covariate is never matched as text", {
  # The same age written 6 and 6.0 is the same age. Matching levels as strings
  # would leave these two modalities with no shared value.
  vals <- list(a = c(6, 6, 14), b = c(6.0, 14.0, 14.0))
  expect_identical(chorale_covariate_kind(vals), "continuous")
})

test_that("a covariate categorical everywhere keeps its shared levels", {
  vals <- list(a = c("F", "M", "F"), b = c("M", "F", "M"))
  expect_identical(chorale_covariate_kind(vals), "categorical")
  designs <- list(
    a = data.frame(sample_id = paste0("a", 1:3), sex = c("F", "M", "F")),
    b = data.frame(sample_id = paste0("b", 1:3), sex = c("M", "F", "M"))
  )
  expect_setequal(chorale_profile_levels(designs, "sex")$sex, c("F", "M"))
})

test_that("a covariate numeric in one modality and categorical in another is refused", {
  expect_true(is.na(chorale_covariate_kind(list(a = c(1, 2, 3), b = c("x", "y")))))
  expect_true(is.na(chorale_covariate_kind(list())))
})

test_that("resolving a signature never keeps a covariate without terms", {
  designs <- list(
    a = data.frame(sample_id = paste0("a", 1:20),
                   phenotype = rep(c("case", "control"), 10),
                   age_months = seq(6, 20, length.out = 20),
                   site = rep(c("p", "q"), 10)),
    b = data.frame(sample_id = paste0("b", 1:20),
                   phenotype = rep(c("case", "control"), 10),
                   age_months = rep(c(6, 14), each = 10),
                   site = rep(c("r", "s"), 10))
  )
  spec <- chorale_resolve_signature(designs)

  # `site` shares no level, so it cannot contribute terms and must not be kept.
  expect_false("site" %in% spec$covariates)
  expect_true("site" %in% spec$excluded$covariate)
  # Whatever survives must have terms to build, which is what the model matrix
  # assumes when it makes a block per covariate.
  expect_true(all(spec$covariates %in% names(spec$levels)))
  expect_silent(chorale_signature_matrix(designs$a, spec))
  expect_silent(chorale_signature_matrix(designs$b, spec))
})

test_that("a two-timepoint design fits rather than failing on its model matrix", {
  # This is the shape that previously reached the model matrix with no levels
  # and failed inside `colnames<-`.
  designs <- list(
    a = data.frame(sample_id = paste0("a", 1:24),
                   phenotype = rep(c("case", "control"), 12),
                   age_months = seq(6, 20, length.out = 24)),
    b = data.frame(sample_id = paste0("b", 1:24),
                   phenotype = rep(c("case", "control"), 12),
                   age_months = rep(c(6, 14), each = 12))
  )
  spec <- chorale_resolve_signature(designs)
  expect_true("age_months" %in% spec$covariates)
  for (d in designs) {
    mm <- chorale_signature_matrix(d, spec)
    expect_equal(nrow(mm$x), nrow(d))
    # One term for the phenotype contrast and one for the continuous covariate,
    # plus the intercept.
    expect_equal(ncol(mm$x), 3L)
  }
})

test_that("a covariate carried without terms is refused with a clear message", {
  design <- data.frame(sample_id = paste0("a", 1:4), phenotype = rep(c("case", "control"), 2))
  broken <- list(covariates = c("phenotype", "ghost"),
                 levels = list(phenotype = c("control", "case")))
  expect_error(chorale_signature_matrix(design, broken),
               class = "chorale_unresolved_covariate")
})
