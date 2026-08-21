test_that("the correspondence the vocabulary asserts is recovered", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               n_random = 100, n_init = 2)
  s <- r$concept_summary

  expect_equal(s$n_concepts, 3L)
  expect_gt(s$recovered_agreement, s$random_baseline)
  expect_gt(s$placement_between_bounds, 0.5)
  expect_equal(s$fraction_partner_correct, 1)
  expect_true(all(r$concept_recovery$partner_correct))
})

test_that("a concept the modalities rank in opposite directions is not recovered", {
  p <- paired_concept_data()
  # The same concept, measured upside down in the second modality. Its score
  # ranks the same animals in the opposite order, which is not recovery.
  flipped <- p$b
  flipped[p$sets$planted, ] <- -flipped[p$sets$planted, ]
  upright <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets, n_random = 100,
                                     n_init = 2)
  inverted <- chorale_destroy_pairing(p$a, flipped, p$design, sets = p$sets, n_random = 100,
                                      n_init = 2)

  planted_up <- upright$concept_recovery[
    upright$concept_recovery$concept == "planted", ]
  planted_down <- inverted$concept_recovery[
    inverted$concept_recovery$concept == "planted", ]
  expect_gt(planted_up$self_agreement, 0.5)
  expect_lt(planted_down$self_agreement, -0.5)
  expect_lt(inverted$concept_summary$recovered_agreement,
            upright$concept_summary$recovered_agreement)
})

test_that("recovery is placed between the two bounds", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               n_random = 100, n_init = 2)
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
                               n_random = 20, n_init = 2)
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
                                sets = strong$sets,
                                n_random = 100, n_init = 2)
  rn <- chorale_destroy_pairing(none$a, none$b, none$design, sets = none$sets,
                                n_random = 100, n_init = 2)
  expect_gt(rs$concept_summary$recovered_agreement,
            rn$concept_summary$recovered_agreement)
})

test_that("a vocabulary that names nothing shared places below random", {
  p <- paired_concept_data(concept_effect = 1.5, seed = 9)
  # Every concept measured upside down in the second modality: the names line
  # up and the biology does not, which is worse than relabelling at random.
  flipped <- p$b
  flipped[unlist(p$sets, use.names = FALSE), ] <-
    -flipped[unlist(p$sets, use.names = FALSE), ]
  r <- chorale_destroy_pairing(p$a, flipped, p$design, sets = p$sets,
                               n_random = 100, n_init = 2)
  expect_lt(r$concept_summary$recovered_agreement,
            r$concept_summary$random_baseline)
  expect_lt(r$concept_summary$placement_between_bounds, 0)
})

test_that("the benchmark returns the concept result and nothing else", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design, sets = p$sets,
                               n_init = 2, n_random = 20)
  expect_equal(nrow(r$concept_summary), 1L)
  expect_setequal(names(r), c("concept_summary", "concept_recovery",
                              "concept_encoding"))
})

test_that("a run without a vocabulary is refused", {
  p <- paired_concept_data()
  # The vocabulary is what the two modalities are compared in, so there is no
  # benchmark to run without one.
  expect_error(chorale_destroy_pairing(p$a, p$b, p$design), "sets")
})

test_that("a vocabulary reaching neither modality is reported, not scored", {
  p <- paired_concept_data()
  r <- chorale_destroy_pairing(p$a, p$b, p$design,
                               sets = list(absent = paste0("no_such_", 1:20)),
                               n_random = 20)
  expect_true(is.na(r$concept_summary$placement_between_bounds))
  expect_match(r$concept_summary$reason, "fewer than two concepts")
})
