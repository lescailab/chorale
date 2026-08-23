report_fit <- function(seed = 1, n_free = 2, ...) {
  fx <- chorale_concept_example(n_concepts = 4L, seed = seed, ...)
  list(fx = fx,
       fit = chorale_concept_fit(fx$containers, fx$sets, n_free = n_free,
                                 n_permutations = 199, n_init = 2))
}

report_additions <- function(seed = 7) {
  fx <- chorale_concept_example(n_samples = 30, n_features = 120,
                                n_modalities = 3, n_concepts = 4,
                                seed = seed)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
                             n_permutations = 19, n_init = 2)
  membership <- chorale_concept_families(fit$concepts, min_overlap = 0.9)
  family_evidence <- chorale_family_evidence(fit$evidence, membership)
  joint_state <- chorale_joint_state(fit$encoding, n_components = 1)
  joint_evidence <- chorale_joint_evidence(joint_state, n_permutations = 19)
  joint_transfer <- chorale_joint_transfer(fit$encoding, n_components = 1,
                                           n_permutations = 19)
  list(fit = fit, family_evidence = family_evidence,
       joint_state = joint_state, joint_evidence = joint_evidence,
       joint_transfer = joint_transfer)
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

test_that("the report writes a full audit trail for supplied analyses", {
  h <- report_additions()
  path <- withr::local_tempdir()
  written <- chorale_report(
    h$fit, path = path, family_evidence = h$family_evidence,
    joint_state = h$joint_state, joint_evidence = h$joint_evidence,
    joint_transfer = h$joint_transfer)

  extra <- c(
    "concept_families.tsv", "family_evidence.tsv", "joint_scores.tsv",
    "joint_loadings.tsv", "joint_variance.tsv", "joint_coverage.tsv",
    "joint_nuisance.tsv", "joint_components.tsv",
    "joint_component_effects.tsv", "joint_transfer.tsv")
  expect_setequal(intersect(basename(written), extra), extra)
  expect_true(all(file.exists(file.path(path, extra))))

  scores <- utils::read.delim(file.path(path, "joint_scores.tsv"))
  loadings <- utils::read.delim(file.path(path, "joint_loadings.tsv"))
  transfer <- utils::read.delim(file.path(path, "joint_transfer.tsv"))
  expect_named(scores, c("sample_id", "modality", "joint_01"))
  expect_named(loadings, c("concept", "joint_01"))
  expect_true(all(c("component", "fitted_component", "loading_agreement",
                    "p_value") %in% names(transfer)))

  html <- paste(readLines(file.path(path, "report.html")), collapse = "\n")
  expect_match(html, "Concept families", fixed = TRUE)
  expect_match(html, "Joint component evidence", fixed = TRUE)
  expect_match(html, "component-specific", fixed = TRUE)
})

test_that("optional analyses are never computed or written implicitly", {
  h <- report_fit()
  path <- withr::local_tempdir()
  written <- chorale_report(h$fit, path = path)
  expect_false(any(startsWith(basename(written), "joint_")))
  expect_false("family_evidence.tsv" %in% basename(written))

  html <- paste(readLines(file.path(path, "report.html")), collapse = "\n")
  expect_match(html, "not supplied", fixed = TRUE)
  expect_match(html, "not a negative or nonsignificant result", fixed = TRUE)
})

test_that("the report validates optional analysis objects and vocabularies", {
  h <- report_additions()
  path <- withr::local_tempdir()
  expect_error(
    chorale_report(h$fit, path = path, joint_state = list()),
    class = "chorale_invalid_report_analysis")

  wrong <- h$family_evidence
  wrong$membership$concept[1] <- "not_in_this_vocabulary"
  expect_error(
    chorale_report(h$fit, path = path, family_evidence = wrong),
    class = "chorale_incompatible_report_analysis")
})

test_that("empty supplied joint results remain distinct from omitted analyses", {
  h <- report_additions()
  state <- chorale_joint_state(h$fit$encoding, n_components = 0)
  evidence <- chorale_joint_evidence(state, n_permutations = 19)
  transfer <- structure(
    list(transfer = h$joint_transfer$transfer[0, , drop = FALSE],
         n_permutations = 19L, alpha = 0.05),
    class = "chorale_joint_transfer")
  path <- withr::local_tempdir()
  chorale_report(h$fit, path = path, joint_state = state,
                 joint_evidence = evidence, joint_transfer = transfer)

  components <- utils::read.delim(file.path(path, "joint_components.tsv"))
  transferred <- utils::read.delim(file.path(path, "joint_transfer.tsv"))
  expect_equal(nrow(components), 0L)
  expect_equal(nrow(transferred), 0L)
  html <- paste(readLines(file.path(path, "report.html")), collapse = "\n")
  expect_match(html, "supplied", fixed = TRUE)
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

test_that("free dimensions report where one modality has none", {
  fx <- encode_fixture(n_features = 300, seed = 5)
  ids <- fx$ids
  # One modality is given fewer samples than the vocabulary has concepts, so
  # its concept scores span its sample space and it has no free dimensions,
  # while the other keeps a residual to decompose.
  small <- fx$containers[[2]][, seq_len(30)]
  containers <- list(A = fx$containers[[1]], B = small)
  set.seed(1)
  sets <- stats::setNames(
    lapply(seq_len(45), function(i) ids[sample(length(ids), 40)]),
    paste0("concept_", seq_len(45)))
  concepts <- chorale_concepts(containers, sets, min_features = 5)
  fit <- chorale_concept_fit(containers, concepts, n_free = "auto",
                             n_permutations = 19L, n_init = 2L)
  n_free <- fit$encoding$variance$n_free
  expect_true(any(n_free == 0L))
  expect_true(any(n_free > 0L))

  free <- chorale_free_dimensions(fit, n_permutations = 19L)
  expect_s3_class(free, "data.frame")
  expect_true(all(free$modality %in% fit$modalities))
  expect_false(any(free$modality %in%
                     fit$encoding$variance$modality[n_free == 0L]))
})
