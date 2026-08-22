family_fixture <- function(n_features = 120, seed = 1) {
  sim <- chorale_simulate(n_modalities = 3, n_features = n_features,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4, effect_size = 3,
                          seed = seed)
  ids <- sprintf("feature_%05d", seq_len(n_features))
  sim$modalities <- lapply(sim$modalities, function(m) {
    rownames(m) <- ids
    m
  })
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  sets <- list(one = ids[1:30], two = ids[25:60], three = ids[55:90],
               four = ids[80:n_features])
  concepts <- chorale_concepts(containers, sets, min_features = 5)
  list(concepts = concepts,
       encoding = chorale_encode(containers, concepts, n_free = 1,
                                 n_init = 2))
}

test_that("families are formed from the vocabulary and not from the data", {
  fx <- family_fixture()
  families <- chorale_concept_families(fx$concepts, min_overlap = 0.1)

  expect_setequal(families$concept, names(fx$concepts$sets))
  expect_true(all(families$family_size >= 1))
  # Overlapping definitions land together; disjoint ones do not.
  expect_equal(families$family[families$concept == "one"],
               families$family[families$concept == "two"])
  expect_false(identical(families$family[families$concept == "one"],
                         families$family[families$concept == "four"]))
})

test_that("a concept in no family of the required size stands alone", {
  fx <- family_fixture()
  families <- chorale_concept_families(fx$concepts, min_overlap = 0.99)
  expect_equal(nrow(families), length(fx$concepts$sets))
  expect_true(all(families$family_size == 1))
  expect_setequal(families$family, families$concept)
})

test_that("min_overlap is checked", {
  fx <- family_fixture()
  expect_error(chorale_concept_families(fx$concepts, min_overlap = 0),
               "must lie in")
  expect_error(chorale_concept_families(fx$concepts, min_overlap = 2),
               "must lie in")
})

test_that("the family statistic is directional, so members that disagree cancel", {
  fx <- family_fixture()
  evidence <- chorale_concept_evidence(fx$encoding, n_permutations = 99)
  families <- chorale_concept_families(fx$concepts, min_overlap = 0.1)
  result <- chorale_family_evidence(evidence, families)

  expect_s3_class(result, "chorale_family_evidence")
  expect_true(all(result$families$p_value >= 1 / 100))
  expect_true(all(result$families$statistic >= 0))
  # A family whose members point opposite ways scores below its largest member.
  disagreeing <- result$families$mean_z
  at <- which.min(abs(disagreeing))
  expect_lt(result$families$statistic[at], result$families$best_single_z[at])
})

test_that("family evidence refuses an object carrying no per-concept null", {
  fx <- family_fixture()
  evidence <- chorale_concept_evidence(fx$encoding, n_permutations = 19)
  families <- chorale_concept_families(fx$concepts, min_overlap = 0.1)
  evidence$null_signed_by_concept <- NULL
  expect_error(chorale_family_evidence(evidence, families),
               class = "chorale_missing_concept_null")
})

test_that("the per-concept null is kept beside the maximum, signed and unsigned", {
  fx <- family_fixture()
  evidence <- chorale_concept_evidence(fx$encoding, n_permutations = 19)
  expect_false(is.null(evidence$null_by_concept))
  expect_false(is.null(evidence$null_signed_by_concept))
  expect_equal(nrow(evidence$null_signed_by_concept), 19L)
  expect_setequal(colnames(evidence$null_signed_by_concept), evidence$joint$key)
  # The magnitudes stored are those of the signed statistics, so a family sum
  # and a single-concept test are calibrated on one set of resamples.
  expect_equal(abs(evidence$null_signed_by_concept), evidence$null_by_concept)
  # Direction is carried rather than reconstructed, which is what keeps the
  # dependence between overlapping concepts in the null.
  expect_true(any(evidence$null_signed_by_concept < 0))
})
