joint_annotation_fixture <- function(seed = 7) {
  fx <- chorale_concept_example(n_samples = 30, n_features = 120,
                                n_modalities = 3, n_concepts = 4, seed = seed)
  fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
                             n_permutations = 19, n_init = 2)
  state <- chorale_joint_state(fit$encoding, n_components = 2)
  list(fit = fit, state = state,
       evidence = chorale_joint_evidence(state, n_permutations = 19),
       membership = chorale_concept_families(fit$concepts, min_overlap = 0.5))
}

test_that("each component is named by the concepts that load on it", {
  h <- joint_annotation_fixture()
  out <- chorale_joint_concepts(h$state, n_top = 2)

  expect_equal(out$component, colnames(h$state$loadings))
  expect_equal(out$share, h$state$variance$share)
  expect_true(all(out$n_concepts > 0))
  named <- unlist(strsplit(stats::na.omit(c(out$concepts_high,
                                            out$concepts_low)), "; "))
  named <- sub(" \\(.*", "", named)
  expect_true(all(named %in% rownames(h$state$loadings)))
  # Two concepts were asked for per direction, so no more can be named.
  expect_true(all(lengths(strsplit(stats::na.omit(out$concepts_high), "; ")) <= 2))
})

test_that("without the phenotype effect the orientation is declared arbitrary", {
  h <- joint_annotation_fixture()
  out <- chorale_joint_concepts(h$state, n_top = 2)
  expect_true(all(out$orientation == "arbitrary"))
  expect_false(any(c("leading_family", "leading_family_share") %in%
                     colnames(out)))
})

test_that("supplying the evidence orients each component by its effect", {
  h <- joint_annotation_fixture()
  plain <- chorale_joint_concepts(h$state, n_top = 3)
  oriented <- chorale_joint_concepts(h$state, evidence = h$evidence, n_top = 3)

  expect_true(all(grepl("^high = raised in ", oriented$orientation)))
  components <- h$evidence$components
  flipped <- components$component[components$effect < 0]
  # A component with a negative effect is turned over, so what it ranked high
  # becomes what it ranks low; one with a positive effect is left alone.
  for (k in seq_len(nrow(oriented))) {
    component <- oriented$component[k]
    at <- match(component, plain$component)
    if (component %in% flipped) {
      expect_equal(oriented$concepts_high[k], plain$concepts_low[at])
    } else {
      expect_equal(oriented$concepts_high[k], plain$concepts_high[at])
    }
  }
})

test_that("a family assignment names the group holding most of the mass", {
  h <- joint_annotation_fixture()
  out <- chorale_joint_concepts(h$state, families = h$membership, n_top = 2)

  expect_true(all(c("leading_family", "leading_family_share",
                    "leading_family_members") %in% colnames(out)))
  expect_true(all(out$leading_family %in% h$membership$family))
  expect_true(all(out$leading_family_share > 0 &
                    out$leading_family_share <= 1))
  members <- strsplit(out$leading_family_members[1], "; ")[[1]]
  expect_true(all(members %in%
                    h$membership$concept[h$membership$family ==
                                           out$leading_family[1]]))
})

test_that("a family table naming no concept of the state is refused", {
  h <- joint_annotation_fixture()
  other <- data.frame(concept = "not_a_concept", family = "family_0001",
                      stringsAsFactors = FALSE)
  expect_error(chorale_joint_concepts(h$state, families = other),
               "names no concept")
  expect_error(chorale_joint_concepts(h$fit), "chorale_joint_state")
})

test_that("the report writes the annotation and names concepts in the page", {
  h <- joint_annotation_fixture()
  family_evidence <- chorale_family_evidence(h$fit$evidence, h$membership)
  path <- withr::local_tempdir()
  written <- chorale_report(h$fit, path = path,
                            family_evidence = family_evidence,
                            joint_state = h$state,
                            joint_evidence = h$evidence)

  expect_true("joint_component_concepts.tsv" %in% basename(written))
  table <- utils::read.delim(file.path(path, "joint_component_concepts.tsv"))
  expect_equal(nrow(table), ncol(h$state$loadings))

  page <- paste(readLines(file.path(path, "report.html")), collapse = "\n")
  expect_true(grepl("What each component is made of", page, fixed = TRUE))
  expect_true(grepl("significant (family-wise)", page, fixed = TRUE))
  # The component is named beside the concepts that define it, so a reader
  # never meets a component identifier on its own.
  expect_true(grepl(rownames(h$state$loadings)[1], page, fixed = TRUE))
})
