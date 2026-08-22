test_that("every modality is scored on the concepts it expresses", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = 2, n_init = 2)

  expect_s3_class(enc, "chorale_encode")
  for (m in names(fx$containers)) {
    s <- enc$encodings[[m]]$concept_scores
    expect_equal(colnames(s), colnames(fx$concepts$membership[[m]]))
    expect_equal(nrow(s), ncol(SummarizedExperiment::assay(fx$containers[[m]])))
    # Scores are standardised, so an effect on one concept is comparable with
    # an effect on another.
    expect_equal(unname(round(apply(s, 2, stats::sd), 6)), rep(1, ncol(s)))
  }
  expect_setequal(chorale_scored_concepts(enc), names(fx$sets))
})

test_that("a concept score is the average level of the features carrying it", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = 0)
  x <- enc$encodings$A$analysis_matrix
  expected <- rowMeans(x[, fx$sets$one, drop = FALSE])
  observed <- enc$encodings$A$concept_scores[, "one"]
  expect_gt(stats::cor(expected, observed), 0.999)
})

test_that("scoring never reads the phenotype", {
  fx <- encode_fixture()
  scrambled <- fx$containers
  for (m in names(scrambled)) {
    d <- as.data.frame(SummarizedExperiment::colData(scrambled[[m]]))
    set.seed(7)
    d$phenotype <- sample(d$phenotype)
    SummarizedExperiment::colData(scrambled[[m]]) <- S4Vectors::DataFrame(d)
  }
  a <- chorale_encode(fx$containers, fx$concepts, n_free = 2, n_init = 2)
  b <- chorale_encode(scrambled, fx$concepts, n_free = 2, n_init = 2)
  expect_equal(a$encodings$A$concept_scores, b$encodings$A$concept_scores)
  expect_equal(a$encodings$A$free_scores, b$encodings$A$free_scores)
})

test_that("free dimensions carry what the concepts do not, and not what they do", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = 3, n_init = 2)
  v <- enc$variance

  expect_true(all(v$n_free == 3L))
  expect_true(all(v$concept_share > 0 & v$concept_share < 1))
  expect_true(all(v$free_share > 0))
  # The three shares partition the variance of the analysis matrix.
  expect_equal(v$concept_share + v$free_share + v$residual_share,
               rep(1, nrow(v)), tolerance = 1e-3)
  # The free dimensions come from the residual, so they hold no concept signal.
  expect_true(all(v$max_free_concept_correlation < 1e-6))
})

test_that("the free-dimension count is honoured per modality", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = c(A = 1, B = 4),
                        n_init = 2)
  expect_equal(enc$variance$n_free[enc$variance$modality == "A"], 1L)
  expect_equal(enc$variance$n_free[enc$variance$modality == "B"], 4L)
  expect_error(chorale_encode(fx$containers, fx$concepts, n_free = c(Z = 2)),
               "unknown modalities")
  expect_error(chorale_encode(fx$containers, fx$concepts, n_free = -1),
               "non-negative")
})

test_that("no free dimensions leaves the concept channel alone", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = 0)
  expect_equal(enc$variance$n_free, c(0L, 0L))
  expect_equal(enc$variance$free_share, c(0, 0))
  expect_equal(ncol(enc$encodings$A$free_scores), 0L)
})

test_that("the automatic count does not exceed what parallel analysis supports", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_init = 2, max_factors = 4L)
  expect_true(all(enc$variance$n_free <= enc$variance$free_ceiling))
  expect_true(all(enc$variance$n_free <= 4L))
})

test_that("each concept's own share of the variance is reported", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = 2, n_init = 2)
  cv <- enc$concept_variance
  expect_setequal(unique(cv$modality), names(fx$containers))
  expect_true(all(cv$variance_share >= 0 & cv$variance_share <= 1))
})

test_that("a collection the concepts were not built on is refused", {
  fx <- encode_fixture()
  expect_error(chorale_encode(fx$containers["A"], fx$concepts),
               "different set of modalities")
  expect_error(chorale_encode(fx$containers, list(a = 1)),
               "chorale_concepts object")
})

test_that("printing states what each channel carries", {
  fx <- encode_fixture()
  enc <- chorale_encode(fx$containers, fx$concepts, n_free = 2, n_init = 2)
  expect_output(print(enc), "concepts scored in every modality: 3")
})

test_that("a vocabulary that spans the sample space leaves no free dimensions", {
  # More concepts than samples: their scores span the whole sample space, the
  # projection reproduces the assay, and the residual is numerical dust.
  fx <- encode_fixture(n_features = 300, seed = 3)
  ids <- fx$ids
  set.seed(1)
  n_concepts <- nrow(SummarizedExperiment::colData(fx$containers[[1]])) + 40L
  many <- stats::setNames(
    lapply(seq_len(n_concepts), function(i) ids[sample(length(ids), 40)]),
    paste0("concept_", seq_len(n_concepts)))
  concepts <- chorale_concepts(fx$containers, many, min_features = 5)
  enc <- chorale_encode(fx$containers, concepts, n_free = "auto", n_init = 2)
  expect_true(all(enc$variance$n_free == 0L))
  expect_true(all(enc$variance$free_ceiling == 0L))
  expect_true(all(enc$variance$concept_share > 0.999))
})
