joint_fixture <- function(n_modalities = 3, n_features = 120, seed = 1,
                          trim_third = FALSE) {
  sim <- chorale_simulate(n_modalities = n_modalities, n_features = n_features,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 4, effect_size = 3,
                          seed = seed)
  ids <- sprintf("feature_%05d", seq_len(n_features))
  sim$modalities <- lapply(sim$modalities, function(m) {
    rownames(m) <- ids
    m
  })
  if (trim_third && n_modalities >= 3) {
    sim$modalities[[3]] <- sim$modalities[[3]][seq_len(n_features / 2), ]
  }
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  sets <- list(one = ids[1:30], two = ids[25:60], three = ids[55:90],
               four = ids[80:n_features])
  concepts <- chorale_concepts(containers, sets, min_features = 5)
  list(containers = containers, concepts = concepts, ids = ids,
       encoding = chorale_encode(containers, concepts, n_free = 1,
                                 n_init = 2))
}

test_that("the stacked state carries every sample and the whole vocabulary", {
  fx <- joint_fixture()
  state <- chorale_joint_state(fx$encoding, n_components = 2)

  n_total <- sum(vapply(fx$encoding$encodings,
                        function(e) nrow(e$concept_scores), integer(1)))
  expect_equal(nrow(state$scores), n_total)
  expect_equal(ncol(state$scores), 2L)
  expect_equal(nrow(state$loadings), length(fx$concepts$vocabulary))
  expect_setequal(rownames(state$loadings), fx$concepts$vocabulary)
  expect_equal(sort(unique(state$modality)), sort(fx$encoding$modalities))
})

test_that("a concept a modality cannot express carries no weight", {
  fx <- joint_fixture(trim_third = TRUE)
  expressed <- fx$concepts$coverage
  expect_false(all(expressed$expressed))

  state <- chorale_joint_state(fx$encoding, n_components = 2)
  third <- fx$encoding$modalities[3]
  reported <- state$coverage$n_concepts[state$coverage$modality == third]
  expect_lt(reported, length(fx$concepts$vocabulary))
  # The narrowest modality does not bound the vocabulary the state is fitted in.
  expect_equal(nrow(state$loadings), length(fx$concepts$vocabulary))
  expect_true(state$converged)
})

test_that("a level shift in one modality does not move the components", {
  fx <- joint_fixture()
  base <- chorale_joint_state(fx$encoding, n_components = 2)

  shifted <- fx$encoding
  first <- shifted$modalities[1]
  shifted$encodings[[first]]$concept_scores <-
    shifted$encodings[[first]]$concept_scores + 25
  moved <- chorale_joint_state(shifted, n_components = 2)

  # Centring within modality is what makes this hold: only a sample's position
  # among samples measured the same way survives the stacking.
  agreement <- abs(stats::cor(base$scores[, 1], moved$scores[, 1]))
  expect_gt(agreement, 0.999)
})

test_that("layer-local nuisance covariates are removed before stacking", {
  fx <- joint_fixture()
  first <- fx$encoding$modalities[1]
  design <- fx$encoding$designs[[first]]
  n <- nrow(design)
  design$study <- rep(c("a", "b"), length.out = n)
  fx$encoding$designs[[first]] <- design

  planted <- fx$encoding
  offset <- ifelse(design$study == "a", 6, -6)
  planted$encodings[[first]]$concept_scores <-
    planted$encodings[[first]]$concept_scores + offset

  unadjusted <- chorale_joint_state(planted, n_components = 2)
  adjusted <- chorale_joint_state(planted, n_components = 2,
                                  nuisance = stats::setNames(list("study"),
                                                             first))
  at <- unadjusted$modality == first
  separation <- function(state) {
    abs(diff(tapply(state$scores[at, 1], design$study, mean)))
  }
  expect_gt(separation(unadjusted), separation(adjusted))
  expect_identical(adjusted$nuisance[[first]], "study")
})

test_that("nuisance may not name the phenotype", {
  fx <- joint_fixture()
  first <- fx$encoding$modalities[1]
  expect_error(
    chorale_joint_state(fx$encoding, n_components = 1,
                        nuisance = stats::setNames(list("phenotype"), first)),
    class = "chorale_nuisance_is_phenotype")
})

test_that("nuisance must name a modality the encoding carries", {
  fx <- joint_fixture()
  expect_error(
    chorale_joint_state(fx$encoding, n_components = 1,
                        nuisance = list(not_a_modality = "study")),
    "does not carry")
})

test_that("the evidence is calibrated on permutations within modality", {
  fx <- joint_fixture()
  state <- chorale_joint_state(fx$encoding, n_components = 2)
  evidence <- chorale_joint_evidence(state, n_permutations = 99)

  expect_s3_class(evidence, "chorale_joint_evidence")
  expect_equal(nrow(evidence$components), 2L)
  expect_true(all(evidence$components$p_value >= 1 / 100))
  expect_true(all(evidence$components$p_family >= evidence$components$p_value))
  expect_equal(evidence$smallest_attainable_p, 1 / 100)
})

test_that("the per-modality effects are reported and nothing is filtered on them", {
  fx <- joint_fixture()
  state <- chorale_joint_state(fx$encoding, n_components = 2)
  evidence <- chorale_joint_evidence(state, n_permutations = 49)

  expect_setequal(unique(evidence$per_modality$modality),
                  fx$encoding$modalities)
  # Every component appears for every modality, whatever its effect there.
  expect_equal(nrow(evidence$per_modality),
               2L * length(fx$encoding$modalities))
})

test_that("holding a modality out needs at least three", {
  fx <- joint_fixture(n_modalities = 2)
  expect_error(chorale_joint_transfer(fx$encoding, n_components = 1,
                                      n_permutations = 19),
               class = "chorale_transfer_needs_three")
})

test_that("a direction fitted without a modality is tested in that modality", {
  fx <- joint_fixture(trim_third = TRUE)
  transfer <- chorale_joint_transfer(fx$encoding, n_components = 2,
                                     n_permutations = 99)

  expect_s3_class(transfer, "chorale_joint_transfer")
  expect_setequal(unique(transfer$transfer$held_out), fx$encoding$modalities)
  expect_true(all(transfer$transfer$p_value >= 1 / 100))
  # The held-out modality is projected only on concepts it can be scored on.
  third <- fx$encoding$modalities[3]
  at <- transfer$transfer$held_out == third
  expect_lt(max(transfer$transfer$n_concepts[at]),
            length(fx$concepts$vocabulary))
  # A modality never contributes to the direction it is judged against.
  expect_false(any(grepl(third, transfer$transfer$trained_on[at],
                         fixed = TRUE)))
})

test_that("holding out works when only some modalities declare nuisance", {
  fx <- joint_fixture(trim_third = TRUE)
  first <- fx$encoding$modalities[1]
  design <- fx$encoding$designs[[first]]
  design$study <- rep(c("a", "b"), length.out = nrow(design))
  fx$encoding$designs[[first]] <- design

  # Only one of three modalities names a nuisance covariate, so each held-out
  # fit sees a nuisance list that does not mention two of its modalities.
  transfer <- chorale_joint_transfer(
    fx$encoding, n_components = 1, n_permutations = 19,
    nuisance = stats::setNames(list("study"), first))
  expect_setequal(unique(transfer$transfer$held_out), fx$encoding$modalities)
})

test_that("held-out components are named against one reference fit", {
  fx <- joint_fixture(trim_third = TRUE)
  transfer <- chorale_joint_transfer(fx$encoding, n_components = 2,
                                     n_permutations = 19)

  expect_true(all(c("component", "fitted_component", "loading_agreement") %in%
                    colnames(transfer$transfer)))
  # Every row carries the label of the reference component it matches, so two
  # rows with the same `component` describe the same direction.
  expect_true(all(!is.na(transfer$transfer$component)))
  # Each held-out fit contributes each reference label at most once per term.
  counts <- table(transfer$transfer$held_out, transfer$transfer$component)
  expect_true(all(counts <= 1))
  agreement <- transfer$transfer$loading_agreement
  expect_true(all(is.na(agreement) | abs(agreement) <= 1 + 1e-8))
})

test_that("matching pairs components whatever their order and sign", {
  loadings <- matrix(stats::rnorm(60), nrow = 20,
                     dimnames = list(paste0("c", 1:20), c("a", "b", "d")))
  # The same directions, reordered and with two of them pointing the other way.
  shuffled <- loadings[, c(3, 1, 2)] %*% diag(c(-1, 1, -1))
  colnames(shuffled) <- c("x", "y", "z")
  matched <- chorale:::chorale_match_components(shuffled, loadings)

  expect_equal(matched$reference[matched$fitted == "x"], "d")
  expect_equal(matched$reference[matched$fitted == "y"], "a")
  expect_equal(matched$reference[matched$fitted == "z"], "b")
  # Sign is recorded rather than hidden, and does not affect the pairing.
  expect_lt(matched$agreement[matched$fitted == "x"], 0)
  expect_gt(matched$agreement[matched$fitted == "y"], 0)
})

test_that("the state reports no component where the collection carries none", {
  fx <- joint_fixture()
  state <- chorale_joint_state(fx$encoding, n_components = 0)
  expect_equal(state$n_components, 0L)
  expect_equal(ncol(state$scores), 0L)
  evidence <- chorale_joint_evidence(state, n_permutations = 19)
  expect_equal(nrow(evidence$components), 0L)
})
