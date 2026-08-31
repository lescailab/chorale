test_that("chorale_load builds a valid container from simulated data", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  se <- chorale_load(sim$modalities[[1]], sim$col_data[[1]])

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), dim(sim$modalities[[1]]))
  expect_true(all(chorale_required_col_data() %in% colnames(SummarizedExperiment::colData(se))))
})

test_that("chorale_load rejects col_data missing the contrast", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3,
                          n_per_cell = 1, seed = 1)
  # A covariate may be absent: it sharpens the comparison where present and is
  # dropped where not, so requiring it would exclude designs the method handles.
  without_strain <- sim$col_data[[1]][, setdiff(names(sim$col_data[[1]]), "strain")]
  expect_s4_class(chorale_load(sim$modalities[[1]], without_strain),
                  "SummarizedExperiment")
  # The contrast may not: without it there is no estimand.
  without_phenotype <- sim$col_data[[1]][, setdiff(names(sim$col_data[[1]]),
                                                   "phenotype")]
  expect_error(chorale_load(sim$modalities[[1]], without_phenotype),
               class = "rlang_error")
})

test_that("chorale_load rejects mismatched sample order", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  shuffled <- sim$col_data[[1]][rev(seq_len(nrow(sim$col_data[[1]]))), ]
  expect_error(chorale_load(sim$modalities[[1]], shuffled), class = "rlang_error")
})

test_that("chorale_load rejects a non-matrix assay", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3, n_per_cell = 1, seed = 1)
  expect_error(
    chorale_load(as.data.frame(sim$modalities[[1]]), sim$col_data[[1]]),
    class = "rlang_error"
  )
})

test_that("blank and placeholder strings are read as missing, not as levels", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 20, n_strains = 3,
                          n_per_cell = 1, seed = 1)
  cd <- sim$col_data[[1]]
  # A design read from a delimited file carries "" where a value was absent.
  cd$phenotype[3:nrow(cd)] <- ""
  cd$sex[1:2] <- "unknown"
  se <- chorale_load(sim$modalities[[1]], cd)
  out <- SummarizedExperiment::colData(se)
  expect_true(all(is.na(out$phenotype[3:nrow(cd)])))
  expect_true(all(is.na(out$sex[1:2])))
  # One real level plus blanks is not a two-level contrast.
  expect_equal(length(unique(stats::na.omit(out$phenotype))), 1L)
})

test_that("design labels are put in one vocabulary as the data are read", {
  d <- data.frame(
    sample_id = paste0("s", 1:8),
    phenotype = c("study_a", "study_b", "label_a", "label_b",
                  "group_a", "group_b",
                  "control", "case"),
    sex = c("male", "female", "F", "M", "Male", "FEMALE", "m", "f"),
    stringsAsFactors = FALSE
  )
  out <- chorale:::chorale_canonical_labels(chorale_blank_to_na(d))
  expect_equal(out$phenotype,
               c("study_a", "study_b", "label_a", "label_b",
                 "group_a", "group_b",
                 "control", "case"))
  expect_equal(out$sex, c("M", "F", "F", "M", "M", "F", "M", "F"))

  # A label outside the vocabulary is left alone rather than guessed at.
  odd <- data.frame(phenotype = "stage_iii", sex = "unknown",
                    stringsAsFactors = FALSE)
  kept <- chorale:::chorale_canonical_labels(odd)
  expect_equal(kept$phenotype, "stage_iii")
  expect_equal(kept$sex, "unknown")
})

test_that("the vocabulary carries no organism, model or study of its own", {
  reg <- chorale_label_registry()
  flat <- unlist(lapply(reg, names))
  # A tool carrying the labels of the data it was built on fits nothing else,
  # so names particular to one study must not ship with the package.
  study_specific <- c("study_case", "study_control", "cohort_code",
                      "dataset_identifier")
  expect_length(intersect(flat, study_specific), 0)
  expect_setequal(unname(reg$phenotype[c("control", "case")]),
                  c("control", "case"))
})

test_that("a study extends the vocabulary rather than editing the package", {
  reg <- chorale_label_registry()
  reg$phenotype <- c(reg$phenotype, "stage_iii" = "case", "stage_0" = "control")
  reg$treatment <- c("drug" = "treated", "placebo" = "untreated")
  d <- data.frame(sample_id = c("a", "b"), phenotype = c("stage_iii", "stage_0"),
                  treatment = c("drug", "placebo"), stringsAsFactors = FALSE)
  out <- chorale:::chorale_canonical_labels(d, reg)
  expect_equal(out$phenotype, c("case", "control"))
  expect_equal(out$treatment, c("treated", "untreated"))
})

test_that("only the columns an estimand needs are required", {
  m <- matrix(stats::rnorm(40), nrow = 10,
              dimnames = list(paste0("g", 1:10), paste0("s", 1:4)))
  minimal <- data.frame(sample_id = paste0("s", 1:4),
                        phenotype = c("case", "case", "control", "control"),
                        stringsAsFactors = FALSE)
  # A design carrying nothing but the identifier and the contrast is enough.
  expect_s4_class(chorale_load(m, minimal), "SummarizedExperiment")
  # Without the contrast there is no estimand, so it is refused.
  expect_error(chorale_load(m, minimal["sample_id"]), "phenotype")
})


test_that("chorale_check_design says what a collection can anchor on", {
  mk <- function(tag, ph, tis) {
    m <- matrix(stats::rnorm(60 * 20), 60,
                dimnames = list(paste0("g", 1:60), paste0(tag, 1:20)))
    d <- data.frame(sample_id = paste0(tag, 1:20), phenotype = rep(ph, each = 10),
                    tissue = tis, stringsAsFactors = FALSE)
    chorale_load(m, d)
  }
  # Tissue is constant within each cohort, so it carries no contrast and cannot
  # anchor, while the phenotype can.
  containers <- list(a = mk("a", c("control", "case"), "source_a"),
                     b = mk("b", c("control", "case"), "source_b"))
  chk <- chorale_check_design(containers)
  expect_true(chk$can_anchor[chk$covariate == "phenotype"])
  expect_false(chk$can_anchor[chk$covariate == "tissue"])
  expect_equal(attr(chk, "usable"), "phenotype")
})



test_that("chorale_check_disjoint is silent on disjoint modalities", {
  mk <- function(ids) {
    m <- matrix(stats::rnorm(20 * length(ids)), 20,
                dimnames = list(paste0("g", 1:20), ids))
    d <- data.frame(sample_id = ids,
                    phenotype = rep(c("control", "case"), length.out = length(ids)),
                    stringsAsFactors = FALSE)
    chorale_load(m, d)
  }
  containers <- list(a = mk(paste0("a", 1:6)), b = mk(paste0("b", 1:6)))
  collisions <- expect_silent(chorale_check_disjoint(containers))
  expect_equal(nrow(collisions), 0L)
  expect_silent(chorale_warn_shared_samples(containers))
})

test_that("chorale_check_disjoint names identifiers two modalities share", {
  mk <- function(ids) {
    m <- matrix(stats::rnorm(20 * length(ids)), 20,
                dimnames = list(paste0("g", 1:20), ids))
    d <- data.frame(sample_id = ids,
                    phenotype = rep(c("control", "case"), length.out = length(ids)),
                    stringsAsFactors = FALSE)
    chorale_load(m, d)
  }
  containers <- list(a = mk(c("s1", "s2", "s3", "s4")),
                     b = mk(c("s3", "s4", "b5", "b6")))
  collisions <- chorale_check_disjoint(containers)
  expect_equal(collisions$sample_id, c("s3", "s4"))
  expect_true(all(collisions$n_modalities == 2L))
  expect_true(all(collisions$modalities == "a, b"))
  expect_warning(chorale_warn_shared_samples(containers),
                 class = "chorale_shared_samples")
})

test_that("a modality carrying no identifier is reported as unchecked", {
  mk <- function(ids) {
    m <- matrix(stats::rnorm(20 * length(ids)), 20,
                dimnames = list(paste0("g", 1:20), ids))
    data.frame(sample_id = ids, stringsAsFactors = FALSE)
  }
  designs <- list(a = mk(paste0("a", 1:4)), b = data.frame(x = 1:4))
  expect_warning(collisions <- chorale_check_disjoint(designs),
                 class = "chorale_unchecked_disjointness")
  # An empty table is not a verdict on the modality that could not be checked.
  expect_equal(nrow(collisions), 0L)
  expect_equal(attr(collisions, "unchecked"), "b")

  # Where no modality carries one, there is nothing to check at all.
  expect_error(
    chorale_check_disjoint(list(a = data.frame(x = 1:3),
                                b = data.frame(x = 1:3))),
    "disjointness cannot be checked")
})

test_that("the collection-level check reaches the fit and the gates", {
  skip_if_not_installed("fastICA")
  fx <- chorale_concept_example(n_samples = 40L, n_features = 60L,
                                n_modalities = 2L, seed = 1L)
  shared <- fx$containers
  # Give the second modality one of the first modality's identifiers.
  cd <- SummarizedExperiment::colData(shared[[2]])
  first <- SummarizedExperiment::colData(shared[[1]])$sample_id[1]
  cd$sample_id[1] <- first
  SummarizedExperiment::colData(shared[[2]]) <- cd
  expect_warning(
    chorale_concept_fit(shared, fx$sets, n_free = 1L, n_permutations = 19L,
                        n_init = 2L),
    class = "chorale_shared_samples")
  expect_warning(
    chorale_gates(shared, control = chorale_control(n_init = 2L),
                  n_surrogate = 2L, n_perm = 5L),
    class = "chorale_shared_samples")
})
