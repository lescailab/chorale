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
    phenotype = c("WT", "transgenic", "Ntg", "tg", "wild-type", "TREATED",
                  "control", "case"),
    sex = c("male", "female", "F", "M", "Male", "FEMALE", "m", "f"),
    stringsAsFactors = FALSE
  )
  out <- chorale:::chorale_canonical_labels(chorale_blank_to_na(d))
  expect_equal(out$phenotype,
               c("control", "case", "control", "case", "control", "case",
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
  study_specific <- c("5xfad", "3xtg", "app/ps1", "ps19", "apoe4", "bxd",
                      "c57bl/6j", "tau", "alzheimer")
  expect_length(intersect(flat, study_specific), 0)
  expect_setequal(unname(reg$phenotype[c("wt", "case")]), c("control", "case"))
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

test_that("a design the package has never seen anchors without configuration", {
  mk <- function(tag, ph) {
    set.seed(nchar(tag))
    m <- matrix(stats::rnorm(80 * 24), 80,
                dimnames = list(paste0("f", 1:80), paste0(tag, 1:24)))
    d <- data.frame(sample_id = paste0(tag, 1:24),
                    phenotype = rep(ph, each = 12),
                    tissue_source = rep(c("liver", "kidney"), 12),
                    donor_age = rep(c(30, 60), 12), stringsAsFactors = FALSE)
    chorale_load(m, d)
  }
  containers <- list(x = mk("x", c("case", "control")),
                     y = mk("y", c("affected", "healthy")))
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
  anchored <- strsplit(unique(fit$matches$shared_covariates), ",")[[1]]
  expect_true("phenotype" %in% anchored)
  expect_true("tissue_source" %in% anchored)
  # The controls run on it too, so nothing downstream needs the column names
  # of one particular study.
  n <- chorale_null(fit, containers, n_permutations = 2, n_init = 2)
  expect_true(is.finite(n$p_phenotype))
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
  containers <- list(a = mk("a", c("WT", "case"), "brain"),
                     b = mk("b", c("control", "tg"), "liver"))
  chk <- chorale_check_design(containers)
  expect_true(chk$can_anchor[chk$covariate == "phenotype"])
  expect_false(chk$can_anchor[chk$covariate == "tissue"])
  expect_equal(attr(chk, "usable"), "phenotype")
})

test_that("age bands follow the cohort rather than an assumed lifespan", {
  # Ages in years, which no band expressed in months could describe.
  d <- data.frame(sample_id = paste0("s", 1:60),
                  age = rep(c(20, 45, 70), each = 20), stringsAsFactors = FALSE)
  out <- chorale:::chorale_add_age_bin(d)
  expect_true("age_bin" %in% colnames(out))
  expect_equal(length(unique(out$age_bin)), 3L)

  # Continuous ages are cut into quantile bands of the ages observed.
  d2 <- data.frame(sample_id = paste0("s", 1:90), age_months = runif(90, 1, 300))
  out2 <- chorale:::chorale_add_age_bin(d2)
  expect_equal(length(unique(out2$age_bin)), 3L)
})
