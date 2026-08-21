test_that("a planted concept is supported and the others are not", {
  fx <- chorale_concept_example(seed = 1)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
                             n_permutations = 199, n_init = 2)
  j <- fit$evidence$joint

  expect_s3_class(fit$evidence, "chorale_concept_evidence")
  expect_equal(j$concept[1], fx$planted)
  expect_true(j$significant[j$concept == fx$planted])
  expect_true(j$family_significant[j$concept == fx$planted])
  expect_false(any(j$significant[j$concept != fx$planted]))
  # The modalities agree on the direction, which is what makes the concept a
  # shared finding rather than one modality carrying the other.
  expect_equal(j$sign_agreement[j$concept == fx$planted], 1)
  expect_equal(j$n_modalities[j$concept == fx$planted], 2L)
})

test_that("effects are reported per concept per modality with their errors", {
  fx <- chorale_concept_example(seed = 2)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                             n_permutations = 99)
  pm <- fit$evidence$per_modality

  expect_setequal(unique(pm$modality), c("A", "B"))
  expect_setequal(unique(pm$concept), names(fx$sets))
  expect_true(all(pm$term == "phenotype=case"))
  expect_true(all(pm$se > 0))
  expect_equal(round(pm$z, 6), round(pm$effect / pm$se, 6))
  expect_true(all(pm$effect[pm$concept == fx$planted] > 0))
})

test_that("the null permutes the design and never refits the encoder", {
  fx <- chorale_concept_example(seed = 3)
  enc <- chorale_encode(fx$containers,
                        chorale_concepts(fx$containers, fx$sets,
                                         min_features = 5),
                        n_free = 1, n_init = 2)
  before <- enc$encodings$A$concept_scores
  ev <- chorale_concept_evidence(enc, n_permutations = 99)
  expect_identical(before, enc$encodings$A$concept_scores)
  expect_equal(ev$smallest_attainable_p, 1 / 100)
  expect_length(ev$null, 99L)
})

test_that("the permutation count sets the resolution and nothing else", {
  fx <- chorale_concept_example(seed = 4)
  cc <- chorale_concepts(fx$containers, fx$sets, min_features = 5)
  enc <- chorale_encode(fx$containers, cc, n_free = 0)
  coarse <- chorale_concept_evidence(enc, n_permutations = 49)
  fine <- chorale_concept_evidence(enc, n_permutations = 499)

  # The statistic is a property of the data, so it does not move with the
  # permutation count; only the smallest reportable p-value does.
  expect_equal(coarse$joint$joint_z, fine$joint$joint_z)
  expect_equal(min(coarse$joint$p_value), 1 / 50)
  expect_equal(min(fine$joint$p_value), 1 / 500)
})

test_that("the vocabulary null is calibrated on data with no signal", {
  # The stage condition: with no concept separating cases from controls, the
  # share of runs calling any concept must sit near the nominal level.
  alpha <- 0.1
  n_sim <- 40L
  called <- vapply(seq_len(n_sim), function(s) {
    fx <- chorale_concept_example(n_samples = 40L, n_features = 60L,
                                  effect = 0, seed = 100L + s)
    fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                               n_permutations = 199, min_features = 5,
                               alpha = alpha, seed = s)
    any(fit$evidence$joint$family_significant)
  }, logical(1))

  rate <- mean(called)
  # Forty runs give a standard error of about 0.05 at this level, so the test
  # is that the rate is near alpha rather than exactly it.
  expect_lt(rate, alpha + 0.15)
})

test_that("overlapping concepts are reported after their neighbours are removed", {
  fx <- chorale_concept_example(seed = 5)
  # A concept that is almost the planted one: on its own it inherits the
  # signal, and the attribution is what separates the two.
  sets <- c(fx$sets, list(shadow = fx$sets$planted[1:22]))
  fit <- chorale_concept_fit(fx$containers, sets, n_free = 0,
                             n_permutations = 199)
  j <- fit$evidence$joint
  shadow <- j[j$concept == "shadow", ]

  expect_gt(shadow$max_jaccard, 0.5)
  expect_gt(abs(shadow$joint_z), 5)
  # With the concept it overlaps regressed out, almost nothing of its own is
  # left, so the signal belongs to the neighbourhood rather than to it.
  expect_lt(abs(shadow$attributed_z), abs(shadow$joint_z) / 2)
})

test_that("a concept only one modality expresses is tested in that modality", {
  fx <- chorale_concept_example(seed = 6)
  b <- SummarizedExperiment::assay(fx$containers$B)
  rownames(b)[1:25] <- paste0("unshared_", 1:25)
  fx$containers$B <- chorale_load(
    b, as.data.frame(SummarizedExperiment::colData(fx$containers$B)))
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                             n_permutations = 99)
  planted <- fit$evidence$joint[fit$evidence$joint$concept == fx$planted, ]

  expect_equal(planted$n_modalities, 1L)
  expect_equal(planted$modalities, "A")
  expect_equal(planted$sign_agreement, 1)
  expect_true(is.na(planted$heterogeneity_p))
})

test_that("an empty vocabulary is refused rather than fitted", {
  fx <- chorale_concept_example(seed = 7)
  expect_error(
    chorale_concept_fit(fx$containers, list(nothing = paste0("absent_", 1:20))),
    class = "chorale_empty_vocabulary")
})

test_that("printing reports what was supported", {
  fx <- chorale_concept_example(seed = 8)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                             n_permutations = 99)
  expect_output(print(fit), "concepts in the vocabulary: 3")
  expect_output(print(fit$evidence), "concepts tested: 3")
})
