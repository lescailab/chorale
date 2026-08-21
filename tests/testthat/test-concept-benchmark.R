test_that("the correspondence the vocabulary asserts is recovered", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               spaces = "concept", n_random = 100, n_init = 2)
  s <- r$concept_summary

  expect_equal(s$n_concepts, 3L)
  expect_gt(s$recovered_agreement, s$random_lower_bound)
  expect_gt(s$placement_between_bounds, 0.5)
  expect_equal(s$fraction_partner_correct, 1)
  expect_true(all(r$concept_recovery$partner_correct))
})

test_that("recovery is placed between the two bounds", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               spaces = "concept", n_random = 100, n_init = 2)
  s <- r$concept_summary
  # The identity the vocabulary asserts cannot beat the best alignment the
  # pairing admits, and the placement is where it sits between the two.
  expect_lte(s$recovered_agreement, s$paired_upper_bound)
  expect_gte(s$placement_between_bounds, 0)
  expect_lte(s$placement_between_bounds, 1)
})

test_that("the pairing is withheld from the encoder", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               spaces = "concept", n_random = 20, n_init = 2)
  ids_a <- rownames(r$concept_encoding$encodings$A$concept_scores)
  ids_b <- rownames(r$concept_encoding$encodings$B$concept_scores)
  expect_length(intersect(ids_a, ids_b), 0)
  expect_true(all(grepl("^a_", ids_a)))
  expect_true(all(grepl("^b_", ids_b)))
})

test_that("a concept with no shared state falls towards the random bound", {
  strong <- paired_concept_data(concept_effect = 1.5, seed = 2)
  none <- paired_concept_data(concept_effect = 0, seed = 2)
  rs <- chorale_destroy_pairing(strong$a, strong$b, strong$design,
                                sets = strong$sets, spaces = "concept",
                                n_random = 100, n_init = 2)
  rn <- chorale_destroy_pairing(none$a, none$b, none$design, sets = none$sets,
                                spaces = "concept", n_random = 100, n_init = 2)
  expect_gt(rs$concept_summary$recovered_agreement,
            rn$concept_summary$recovered_agreement)
})

test_that("both spaces can be scored in one run", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               n_factors = 3, n_init = 2, n_random = 20)
  expect_equal(nrow(r$summary), 1L)
  expect_equal(nrow(r$concept_summary), 1L)
  expect_s3_class(r$fit, "chorale_fit")
})

test_that("the factor space alone is what a run without a vocabulary scores", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, n_factors = 3, n_init = 2,
                               n_random = 20)
  expect_equal(nrow(r$summary), 1L)
  expect_null(r$concept_summary)
  expect_error(
    chorale_destroy_pairing(p$a, p$b, p$design, spaces = "concept"),
    "needs `sets`")
})

test_that("a vocabulary reaching neither modality is reported, not scored", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design,
                               sets = list(absent = paste0("no_such_", 1:20)),
                               spaces = "concept", n_random = 20)
  expect_true(is.na(r$concept_summary$placement_between_bounds))
  expect_match(r$concept_summary$reason, "fewer than two concepts")
})
