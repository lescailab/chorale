report_fit <- function(seed = 1, n_free = 2, ...) {
  fx <- chorale_concept_example(n_concepts = 4L, seed = seed, ...)
  list(fx = fx,
       fit = chorale_concept_fit(fx$containers, fx$sets, n_free = n_free,
                                 n_permutations = 199, n_init = 2))
}

test_that("the report writes the concept, free-dimension and added-value tables", {
  h <- report_fit()
  path <- withr::local_tempdir()
  written <- chorale_report(h$fit, null = chorale_null(h$fit, n_shuffles = 0),
                            path = path)

  expect_true(all(c("concepts.tsv", "free_dimensions.tsv", "added_value.tsv",
                    "concept_coverage.tsv", "controls.tsv", "report.html") %in%
                    basename(written)))
  expect_true(all(file.exists(written)))
  concepts <- utils::read.delim(file.path(path, "concepts.tsv"))
  expect_equal(nrow(concepts), length(h$fx$sets))
})

test_that("the concept table carries the effect in each modality and the joint one", {
  h <- report_fit()
  path <- withr::local_tempdir()
  chorale_report(h$fit, path = path)
  concepts <- utils::read.delim(file.path(path, "concepts.tsv"))

  expect_true(all(c("effect_A", "z_A", "effect_B", "z_B", "joint_z",
                    "sign_agreement", "attributed_z", "in_all_modalities",
                    "q_value", "beats_best_single") %in% names(concepts)))
  planted <- concepts[concepts$concept == h$fx$planted, ]
  expect_true(planted$significant)
  expect_true(planted$in_all_modalities)
  expect_gt(planted$effect_A, 0)
  expect_gt(planted$effect_B, 0)
})

test_that("free dimensions are reported separately from the concepts", {
  h <- report_fit()
  free <- chorale_free_dimensions(h$fit, n_permutations = 99)

  expect_equal(nrow(free), 4L)
  expect_setequal(free$modality, c("A", "B"))
  expect_true(all(c("variance_share", "reproducibility", "p_family",
                    "outside_vocabulary", "top_features") %in% names(free)))
  # Each dimension carries what it reconstructs, which the fitted loadings make
  # unequal; the shares of one modality sum to that modality's free share.
  expect_false(isTRUE(all.equal(free$variance_share[free$modality == "A"][1],
                               free$variance_share[free$modality == "A"][2])))
  for (m in c("A", "B")) {
    expect_equal(sum(free$variance_share[free$modality == m]),
                 h$fit$encoding$variance$free_share[
                   h$fit$encoding$variance$modality == m],
                 tolerance = 1e-3)
  }
  # The features naming a dimension are a description of it, and there are as
  # many as were asked for.
  expect_equal(lengths(strsplit(free$top_features, ", ")), rep(10L, 4L))
  # Nothing outside the vocabulary moves with the phenotype in this collection.
  expect_false(any(free$outside_vocabulary))
})

test_that("a phenotype-linked direction outside the vocabulary is called", {
  # No concept separates cases from controls here, so what does is outside the
  # vocabulary entirely rather than sharing a direction with something in it.
  fx <- chorale_concept_example(n_concepts = 3L, effect = 0, seed = 4)
  # Shift the features no concept covers, in cases only: coordinated variation
  # the vocabulary has no name for.
  uncovered <- setdiff(
    rownames(SummarizedExperiment::assay(fx$containers$A)),
    unlist(fx$sets, use.names = FALSE))
  for (m in c("A", "B")) {
    x <- SummarizedExperiment::assay(fx$containers[[m]])
    d <- as.data.frame(SummarizedExperiment::colData(fx$containers[[m]]))
    x[uncovered, d$phenotype == "case"] <-
      x[uncovered, d$phenotype == "case"] + 3
    fx$containers[[m]] <- chorale_load(x, d)
  }
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
                             n_permutations = 199, n_init = 2)
  free <- chorale_free_dimensions(fit, n_permutations = 199)

  expect_true(any(free$outside_vocabulary))
  expect_lt(min(free$p_family), 0.05)
})

test_that("no per-dimension share is reported where the cross-terms matter", {
  h <- report_fit()
  encoding <- h$fit$encoding$encodings$A
  # A channel share the components cannot account for: the allocation is
  # refused rather than forced.
  expect_true(all(is.na(chorale_component_variance(encoding, 0.9))))
  expect_true(all(is.finite(chorale_component_variance(
    encoding, h$fit$encoding$variance$free_share[1]))))
})

test_that("added value says whether more than one modality was needed", {
  h <- report_fit()
  added <- chorale_added_value(h$fit)
  expect_true(all(c("joint_z", "best_single_z", "margin", "needs_multiple") %in%
                    names(added)))
  planted <- added[added$concept == h$fx$planted, ]
  expect_true(planted$needs_multiple)
  expect_gt(planted$margin, 0)
})

test_that("a concept only one modality expresses cannot need more than one", {
  fx <- chorale_concept_example(seed = 5)
  b <- SummarizedExperiment::assay(fx$containers$B)
  rownames(b)[rownames(b) %in% fx$sets$planted] <-
    paste0("unshared_", seq_along(fx$sets$planted))
  fx$containers$B <- chorale_load(
    b, as.data.frame(SummarizedExperiment::colData(fx$containers$B)))
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                             n_permutations = 99)
  added <- chorale_added_value(fit)
  planted <- added[added$concept == "planted", ]
  expect_equal(planted$n_modalities, 1L)
  expect_false(planted$needs_multiple)
})

test_that("the report runs without controls, and refuses what is not a fit", {
  h <- report_fit()
  path <- withr::local_tempdir()
  written <- chorale_report(h$fit, path = path)
  expect_false("controls.tsv" %in% basename(written))
  expect_false(file.exists(file.path(path, "controls.tsv")))
  expect_error(chorale_report(list(), path = path), "chorale_concept_fit")
  expect_error(chorale_free_dimensions(list()), "chorale_concept_fit")
})

test_that("a fit with no free dimensions writes an empty free-dimension table", {
  h <- report_fit(n_free = 0)
  path <- withr::local_tempdir()
  chorale_report(h$fit, path = path)
  free <- utils::read.delim(file.path(path, "free_dimensions.tsv"))
  expect_equal(nrow(free), 0L)
  expect_false(file.exists(file.path(path, "free_scores_A.tsv")))
})
