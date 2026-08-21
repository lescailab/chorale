test_that("the planted concepts are called and the rest stay quiet", {
  fx <- chorale_concept_example(n_concepts = 5L, n_planted = 2L,
                                n_features = 150L, seed = 1)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                             n_permutations = 199)
  scored <- chorale_score_concepts(fit, fx$planted)

  expect_equal(scored$summary$n_planted, 2L)
  expect_equal(scored$summary$recall, 1)
  # A ranking is threshold-free, so it says how far the planted concepts are
  # separated whatever cut is applied.
  expect_equal(scored$summary$rank_auc, 1)
  expect_true(all(scored$per_concept$planted[1:2]))

  # Holding the false discovery rate over the vocabulary permits an occasional
  # false call; holding the probability of any false call does not.
  strict <- chorale_score_concepts(fit, fx$planted, criterion = "family_wise")
  expect_equal(strict$summary$recall, 1)
  expect_equal(strict$summary$false_positive_rate, 0)
})

test_that("a planting with nothing in it is recovered as nothing", {
  fx <- chorale_concept_example(effect = 0, seed = 2)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
                             n_permutations = 199)
  scored <- chorale_score_concepts(fit, fx$planted)
  expect_equal(scored$summary$recall, 0)
  expect_equal(scored$summary$false_positive_rate, 0)
})

test_that("planting at the concept level keeps the sets that would be refused", {
  profiles <- toy_profiles()
  sets <- list(set_a = 1:20, set_b = 21:40, set_c = 41:60)
  membership <- toy_membership(profiles, sets)
  vocabulary <- list(set_a = as.character(1:20), set_b = as.character(21:40),
                     set_c = as.character(41:60))

  # At the factor level a planting set that coincides with a scoring set is
  # refused, because the pathway channel would recover it by construction.
  expect_error(
    chorale_plant(profiles, membership, vocabulary, vocabulary,
                  n_programmes = 2L, min_features = 5L),
    "max_jaccard")

  # At the concept level the two vocabularies are the same one on purpose: a
  # concept is planted and recovered by name.
  planted <- chorale_plant(profiles, membership, vocabulary, vocabulary,
                           n_programmes = 2L, min_features = 5L,
                           level = "concept")
  expect_equal(planted$level, "concept")
  expect_length(planted$concepts, 2L)
  expect_true(all(planted$concepts %in% names(vocabulary)))
})

test_that("a planting from chorale_plant names the concepts it planted", {
  profiles <- toy_profiles()
  sets <- list(set_a = 1:20, set_b = 21:40, set_c = 41:60)
  membership <- toy_membership(profiles, sets)
  vocabulary <- list(set_a = as.character(1:20), set_b = as.character(21:40),
                     set_c = as.character(41:60))
  planted <- chorale_plant(profiles, membership, vocabulary, vocabulary,
                           n_programmes = 2L, min_features = 5L,
                           level = "concept")
  expect_equal(chorale_planted_names(planted), planted$concepts)
  expect_error(chorale_score_concepts(list(), "a"), "chorale_concept_fit")
})

test_that("recovery falls away as the effect the concept carries falls away", {
  grid <- data.frame(label = c("strong", "absent"), effect = c(1, 0))
  curve <- chorale_validate_concepts(grid, n_rep = 2, n_permutations = 199)

  expect_equal(nrow(curve), 2L)
  expect_true(all(c("recall", "false_positive_rate", "rank_auc",
                    "n_evaluated") %in% names(curve)))
  expect_equal(curve$n_evaluated, c(2L, 2L))
  expect_gt(curve$recall[curve$label == "strong"],
            curve$recall[curve$label == "absent"])
})

test_that("the curve runs over sample size and vocabulary coverage", {
  grid <- data.frame(label = c("small", "large"), n_samples = c(24L, 80L),
                     effect = 0.4)
  by_size <- chorale_validate_concepts(grid, n_rep = 2, n_permutations = 199)
  expect_gte(by_size$rank_auc[by_size$label == "large"],
             by_size$rank_auc[by_size$label == "small"])

  covered <- data.frame(label = c("full", "partial"), coverage = c(1, 0.4),
                        n_features = 200L)
  by_coverage <- chorale_validate_concepts(covered, n_rep = 1,
                                           n_permutations = 199)
  expect_equal(by_coverage$n_evaluated, c(1L, 1L))
})

test_that("a regime the estimator cannot be run in says why", {
  grid <- data.frame(label = "too_sparse", coverage = 0.1, n_features = 90L)
  curve <- chorale_validate_concepts(grid, n_rep = 1, n_permutations = 99)
  expect_equal(curve$n_evaluated, 0L)
  expect_match(curve$reason, "cannot be scored on this vocabulary")
  expect_true(is.na(curve$recall))
})

test_that("an empty grid is refused", {
  expect_error(chorale_validate_concepts(data.frame()), "at least one regime")
})
