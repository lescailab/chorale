# The conditions of applicability, evaluated before a collection is fitted.
#
# The distributional conditions run in Python through reticulate, so every test
# touching them is skipped where the bound interpreter or its packages are
# absent. The conditions evaluated in R are tested unconditionally.

sim_collection <- function(n_features = 120L, seed = 1L) {
  sim <- chorale_simulate(n_modalities = 2L, n_features = n_features,
                          n_shared_factors = 2L, n_private_factors = 1L,
                          n_strains = 4L, n_per_cell = 2L, seed = seed)
  Map(function(assay, cd) chorale_load(assay, col_data = cd),
      sim$modalities, sim$col_data)
}

skip_without_python <- function() {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_available(initialize = TRUE),
              "no Python interpreter is bound")
  skip_if_not(reticulate::py_module_available("scipy"), "scipy is unavailable")
  skip_if_not(reticulate::py_module_available("pandas"), "pandas is unavailable")
}

test_that("the detectability condition reports a usable factor count", {
  containers <- sim_collection()
  xs <- lapply(containers, function(x) scale(t(SummarizedExperiment::assay(x))))
  d <- chorale_gate_detectability(xs, n_perm = 10L, max_factors = 6L, seed = 1L)

  expect_s3_class(d, "data.frame")
  expect_setequal(d$modality, names(containers))
  # The count is bounded below by two and above by the ceiling and by a fifth
  # of the samples, so it can never exceed what the modality can support.
  expect_true(all(d$n_factors >= 2L))
  expect_true(all(d$n_factors <= 6L))
  expect_true(all(d$n_factors <= floor(d$n_samples / 5)))
  expect_true(all(d$gamma_p_over_n > 0))
})

test_that("the anchor-richness condition coarsens down to the phenotype", {
  containers <- sim_collection()
  designs <- lapply(containers, function(x) {
    as.data.frame(SummarizedExperiment::colData(x))
  })
  a <- chorale_gate_anchors(designs)

  expect_s3_class(a, "data.frame")
  # Every coarsening includes the phenotype, and the coarsest is it alone.
  expect_true(all(grepl("phenotype", a$coarsening)))
  expect_true("phenotype" %in% a$coarsening)
  # Coarsening merges strata, so each level has no more strata than the one
  # above it, and the phenotype alone is the fewest.
  by_level <- stats::aggregate(n_strata ~ coarsening, data = a, FUN = max)
  by_level <- by_level[match(unique(a$coarsening), by_level$coarsening), ]
  expect_false(is.unsorted(rev(by_level$n_strata)))
  expect_equal(by_level$n_strata[by_level$coarsening == "phenotype"],
               min(by_level$n_strata))
  # Shared strata can never exceed what a modality populates.
  expect_true(all(a$n_strata_shared <= a$n_strata))
})

test_that("a design with no shared varying covariate is refused", {
  designs <- list(
    a = data.frame(sample_id = c("s1", "s2"), phenotype = c("case", "case")),
    b = data.frame(sample_id = c("s3", "s4"), phenotype = c("case", "case"))
  )
  expect_error(chorale_gate_anchors(designs), "anchored")
})

test_that("the component estimator handed to the gates is this package's own", {
  set.seed(1)
  x <- scale(matrix(stats::rt(60 * 40, df = 3), nrow = 60))
  ica_fn <- chorale_gate_ica(n_init = 2L, consensus = FALSE)
  s <- ica_fn(x, 3L, 1L)

  expect_true(is.matrix(s))
  expect_equal(nrow(s), nrow(x))
  expect_equal(ncol(s), 3L)
  # The same seed must give the same components, so a gate verdict is
  # reproducible from the arguments that produced it.
  expect_equal(s, ica_fn(x, 3L, 1L))
})

test_that("a matrix keeps its orientation crossing into Python and back", {
  skip_without_python()
  py <- chorale_gates_python()
  # Deliberately non-square, so a transpose cannot pass unnoticed, and with
  # known contents so the correspondence can be checked cell by cell.
  n <- 7L
  p <- 3L
  x <- matrix(seq_len(n * p), nrow = n, ncol = p)

  # The estimator is called once on the matrix and again on each surrogate, so
  # only the first call is the matrix that was sent.
  seen <- NULL
  recorder <- function(xx, k, seed) {
    if (is.null(seen)) seen <<- xx
    matrix(stats::rnorm(nrow(xx) * k), nrow = nrow(xx), ncol = k)
  }
  invisible(py$gate_nongaussianity(x, 2L, "m", recorder, seed = 1L,
                                   n_surrogate = 1L))

  # The estimator must see the same matrix R sent: same shape, same cells,
  # not its transpose.
  expect_equal(dim(seen), c(n, p))
  expect_equal(as.matrix(seen), x, ignore_attr = TRUE)
})

test_that("components are numbered from one, as R counts them", {
  skip_without_python()
  py <- chorale_gates_python()
  set.seed(1)
  s <- scale(matrix(stats::rt(60 * 4, df = 3), nrow = 60))
  rows <- py$component_stats(s)
  ids <- vapply(rows, function(r) as.integer(r$component), integer(1))

  # Python enumerates from zero; what crosses the boundary must not.
  expect_equal(sort(ids), 1:4)
  expect_false(any(ids == 0L))
})

test_that("the Anderson-Darling critical value is found by its level, not its position", {
  skip_without_python()
  py <- chorale_gates_python()
  set.seed(1)
  v <- matrix(stats::rnorm(200), ncol = 1)

  at5 <- py$component_stats(v, alpha = 0.05)[[1]]$anderson_critical_value
  at1 <- py$component_stats(v, alpha = 0.01)[[1]]$anderson_critical_value
  # A stricter level has a larger critical value; equality would mean the
  # level was ignored and a fixed position read instead.
  expect_gt(at1, at5)
  # A level the table does not carry must fail rather than fall back.
  expect_error(py$component_stats(v, alpha = 0.04), "per cent")
})

test_that("the gate module refuses to fit components itself", {
  skip_without_python()
  py <- chorale_gates_python()
  x <- scale(matrix(stats::rnorm(40 * 20), nrow = 40))
  # Passing an estimator of the wrong shape must fail loudly rather than be
  # silently reshaped, because the components decide every verdict.
  bad <- function(x, k, seed) matrix(0, nrow = 2L, ncol = k)
  expect_error(
    py$gate_nongaussianity(x, 2L, "m", bad, seed = 1L, n_surrogate = 1L),
    "samples by components"
  )
  # A transposed return is the same failure and must be refused too.
  transposed <- function(x, k, seed) matrix(0, nrow = k, ncol = nrow(x))
  expect_error(
    py$gate_nongaussianity(x, 2L, "m", transposed, seed = 1L, n_surrogate = 1L),
    "samples by components"
  )
})

test_that("non-Gaussian data passes the non-Gaussianity condition and Gaussian data does not", {
  skip_without_python()
  py <- chorale_gates_python()
  ica_fn <- chorale_gate_ica(n_init = 2L, consensus = FALSE)

  set.seed(1)
  heavy <- scale(matrix(stats::rt(80 * 30, df = 2), nrow = 80))
  gaussian <- scale(matrix(stats::rnorm(80 * 30), nrow = 80))

  # The gate returns a pair: the summary, then the per-component detail.
  h <- py$gate_nongaussianity(heavy, 3L, "heavy", ica_fn, seed = 1L,
                              n_surrogate = 20L)[[1]]
  g <- py$gate_nongaussianity(gaussian, 3L, "gaussian", ica_fn, seed = 1L,
                              n_surrogate = 20L)[[1]]

  # Calibration is the point: ICA maximises non-Gaussianity, so Gaussian data
  # still yields components far from normal in absolute terms. What separates
  # the two is the comparison with their own surrogates.
  expect_lt(h$p_value, g$p_value)
  expect_gt(h$median_A2_observed, h$median_A2_surrogate)
})

test_that("a modality cannot differ from a copy of itself", {
  skip_without_python()
  py <- chorale_gates_python()
  set.seed(1)
  s <- scale(matrix(stats::rt(80 * 3, df = 3), nrow = 80))
  d <- chorale_records_to_df(py$gate_modality_difference(list(a = s, b = s)))

  expect_equal(nrow(d), 1L)
  expect_identical(d$verdict, "fail")
  expect_gt(d$pooled_KS_p, 0.05)
})

test_that("chorale_gates runs end to end and reports every condition", {
  skip_without_python()
  containers <- sim_collection()
  g <- chorale_gates(containers, control = chorale_control(n_init = 2L),
                     n_surrogate = 5L, n_perm = 10L, seed = 1L)

  expect_s3_class(g, "chorale_gates")
  expect_named(g, c("non_gaussianity", "non_gaussianity_components",
                    "modality_difference", "detectability",
                    "anchor_richness", "n_factors"))
  expect_setequal(g$non_gaussianity$modality, names(containers))
  expect_true(all(g$non_gaussianity$verdict %in% c("pass", "fail")))
  expect_true(all(g$modality_difference$verdict %in% c("pass", "fail")))
  expect_output(print(g), "non-Gaussianity")
})
