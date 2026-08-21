concept_null_fit <- function(seed = 1, ...) {
  fx <- chorale_concept_example(seed = seed)
  list(fx = fx,
       fit = chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
                                 n_permutations = 199, n_init = 2, ...))
}

test_that("the controls reuse the concept scores rather than refitting them", {
  h <- concept_null_fit()
  before <- h$fit$encoding$encodings$A$concept_scores
  nl <- chorale_null(h$fit, h$fx$containers, n_shuffles = 0)

  expect_s3_class(nl, "chorale_concept_null")
  expect_s3_class(nl, "chorale_null")
  expect_identical(before, h$fit$encoding$encodings$A$concept_scores)
  # No count was asked for, so the permutations the fit already paid for are
  # the permutations the controls read.
  expect_identical(nl$phenotype_null, h$fit$evidence$null)
  expect_equal(nl$n_permutations, 199L)
})

test_that("each control states the smallest p-value it can attain", {
  h <- concept_null_fit()
  nl <- chorale_null(h$fit, h$fx$containers, n_shuffles = 5)
  ct <- nl$controls

  expect_equal(nrow(ct), 3L)
  expect_setequal(ct$control, c("phenotype permutation", "modality shuffle",
                                "initialisation stability"))
  phenotype <- ct[ct$control == "phenotype permutation", ]
  expect_equal(phenotype$smallest_attainable_p, signif(1 / 200, 3))
  expect_gte(phenotype$p_value, phenotype$smallest_attainable_p)
  shuffle <- ct[ct$control == "modality shuffle", ]
  expect_equal(shuffle$n_resamples, 5L)
  expect_equal(shuffle$smallest_attainable_p, signif(1 / 6, 3))
})

test_that("the permutation count is a choice about resolution", {
  h <- concept_null_fit()
  coarse <- chorale_null(h$fit, n_permutations = 49L, n_shuffles = 0)
  fine <- chorale_null(h$fit, n_permutations = 999L, n_shuffles = 0)

  expect_equal(coarse$observed_statistic, fine$observed_statistic)
  expect_equal(coarse$controls$smallest_attainable_p[1], signif(1 / 50, 3))
  expect_equal(fine$controls$smallest_attainable_p[1], signif(1 / 1000, 3))
})

test_that("the phenotype permutation places the planted concept above the null", {
  h <- concept_null_fit()
  nl <- chorale_null(h$fit, n_shuffles = 0)
  expect_lt(nl$p_phenotype, 0.01)
  expect_gt(nl$observed_statistic, max(nl$phenotype_null))
})

test_that("the modality shuffle reports inapplicability rather than passing", {
  h <- concept_null_fit(seed = 2)
  disjoint <- h$fx$containers
  b <- SummarizedExperiment::assay(disjoint$B)
  rownames(b) <- paste0("unshared_", seq_len(nrow(b)))
  disjoint$B <- chorale_load(
    b, as.data.frame(SummarizedExperiment::colData(disjoint$B)))

  nl <- chorale_null(h$fit, disjoint, n_shuffles = 3)
  shuffle <- nl$controls[nl$controls$control == "modality shuffle", ]
  expect_false(shuffle$applicable)
  expect_match(shuffle$reason, "fewer than ten features")
  expect_true(is.na(shuffle$p_value))
})

test_that("without the containers the shuffle says why it did not run", {
  h <- concept_null_fit()
  nl <- chorale_null(h$fit)
  shuffle <- nl$controls[nl$controls$control == "modality shuffle", ]
  expect_false(shuffle$applicable)
  expect_match(shuffle$reason, "were not supplied")
})

test_that("stability is reported per modality for the free dimensions", {
  h <- concept_null_fit()
  nl <- chorale_null(h$fit, n_shuffles = 0)
  expect_setequal(nl$stability$modality, c("A", "B"))
  expect_true(all(nl$stability$n_free == 1L))
  expect_true(all(nl$stability$subspace_agreement >= 0 &
                    nl$stability$subspace_agreement <= 1))
})

test_that("an object that is not a fit is refused by name", {
  expect_error(chorale_null(list(a = 1)), "chorale_concept_fit")
})

test_that("printing lists every control", {
  h <- concept_null_fit()
  nl <- chorale_null(h$fit, n_shuffles = 0)
  expect_output(print(nl), "phenotype permutation")
  expect_output(print(nl), "initialisation stability")
})
