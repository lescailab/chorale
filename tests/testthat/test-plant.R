test_that("overlap reports the nearest scoring set and its Jaccard index", {
  a <- list(one = c("1", "2", "3", "4"), two = c("9", "10"))
  b <- list(x = c("1", "2", "3", "4"), y = c("5", "6"))
  ov <- chorale_set_overlap(a, b)
  expect_equal(ov$max_jaccard[ov$set == "one"], 1)
  expect_equal(ov$nearest[ov$set == "one"], "x")
  expect_equal(ov$max_jaccard[ov$set == "two"], 0)
})

test_that("planting places each set on the features that belong to it", {
  profiles <- toy_profiles()
  sets <- list(set_a = 1:20, set_b = 21:40, set_c = 41:60)
  membership <- toy_membership(profiles, sets)
  plant_sets <- list(set_a = as.character(1:20), set_b = as.character(21:40),
                     set_c = as.character(41:60))
  score_sets <- list(other = as.character(200:260))

  planted <- chorale_plant(profiles, membership, plant_sets, score_sets,
                           n_concepts = 2L, n_private_factors = 1L,
                           min_features = 5L, seed = 1L)

  expect_length(planted$sim$modalities, 2)
  expect_equal(names(planted$sim$modalities), c("A", "B"))
  expect_equal(nrow(planted$plant), 2 * 2)
  # A planted feature belongs to the set it was planted from.
  for (i in seq_len(nrow(planted$plant))) {
    s <- planted$plant$set[i]
    m <- planted$plant$modality[i]
    members <- rownames(membership[[m]])[membership[[m]][, s] > 0]
    expect_true(all(planted$plant$planted[[i]] %in% members))
    expect_true(all(planted$plant$markers[[i]] %in% planted$plant$planted[[i]]))
  }
})

test_that("a set overlapping the scoring vocabulary is refused", {
  profiles <- toy_profiles()
  sets <- list(set_a = 1:20, set_b = 21:40)
  membership <- toy_membership(profiles, sets)
  plant_sets <- list(set_a = as.character(1:20), set_b = as.character(21:40))
  score_sets <- list(mirror_a = as.character(1:20),
                     mirror_b = as.character(21:40))
  expect_error(
    chorale_plant(profiles, membership, plant_sets, score_sets,
                  n_concepts = 1L, min_features = 5L, max_jaccard = 0.4),
    "max_jaccard"
  )
})

test_that("planting refuses membership that does not match the profile", {
  profiles <- toy_profiles()
  sets <- list(set_a = 1:20, set_b = 21:40)
  membership <- toy_membership(profiles, sets)
  rownames(membership$A) <- rev(rownames(membership$A))
  expect_error(
    chorale_plant(profiles, membership, list(set_a = as.character(1:20)),
                  list(other = as.character(1:5)), n_concepts = 1L),
    "feature identifiers"
  )
})

test_that("the ranking statistic is a probability that truth outranks the rest", {
  expect_equal(chorale_rank_auc(c(3, 2, 1), c(TRUE, FALSE, FALSE)), 1)
  expect_equal(chorale_rank_auc(c(1, 2, 3), c(TRUE, FALSE, FALSE)), 0)
  expect_true(is.na(chorale_rank_auc(c(1, 2, 3), c(FALSE, FALSE, FALSE))))
})

test_that("planted signatures give every concept a different design response", {
  sig <- chorale_spread_signature(3, c("phenotype", "sex"))
  expect_equal(dim(sig), c(3, 2))
  unit <- sig / sqrt(rowSums(sig^2))
  cosines <- abs(tcrossprod(unit))
  diag(cosines) <- 0
  expect_true(max(cosines) < 0.95)
})
