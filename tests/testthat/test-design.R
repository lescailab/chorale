test_that("phenotype reference and adjusted effects are explicit", {
  set.seed(1)
  n <- 400L
  age <- stats::rnorm(n)
  phenotype <- ifelse(age + stats::rnorm(n, sd = 0.7) > 0, "case", "control")
  design <- data.frame(sample_id = paste0("s", seq_len(n)), phenotype, age)
  scores <- cbind(
    factor_1 = 1.5 * (phenotype == "case") + age + stats::rnorm(n),
    factor_2 = age + stats::rnorm(n))
  rownames(scores) <- design$sample_id
  spec <- chorale_resolve_signature(
    list(a = design, b = transform(design, sample_id = paste0("t", seq_len(n)))),
    phenotype_reference = "control")
  p <- chorale_adjusted_profile(scores, design, spec)
  expect_equal(spec$levels$phenotype[1], "control")
  expect_gt(p$effects["factor_1", "phenotype=case"], 0.5)
  expect_lt(abs(p$effects["factor_2", "phenotype=case"]), 0.8)
  expect_true(all(is.finite(p$se)))
})

test_that("phenotype is mandatory and categorical", {
  a <- data.frame(sample_id = paste0("a", 1:6), group = rep(c("x", "y"), 3))
  b <- data.frame(sample_id = paste0("b", 1:6), group = rep(c("x", "y"), 3))
  expect_error(chorale_resolve_signature(list(a = a, b = b)),
               class = "chorale_missing_phenotype")
})

test_that("incompatible phenotype levels fail and unusable covariates are reported", {
  a <- data.frame(sample_id = paste0("a", 1:8),
                  phenotype = rep(c("control", "case"), 4),
                  batch = c(rep(NA, 6), "one", "two"))
  b <- data.frame(sample_id = paste0("b", 1:8),
                  phenotype = rep(c("control", "case"), 4),
                  batch = rep(c("one", "two"), 4))
  spec <- chorale_resolve_signature(list(a = a, b = b))
  expect_equal(spec$excluded$reason[spec$excluded$covariate == "batch"],
               "excessive missingness")
  b$phenotype[1] <- "extra_level"
  expect_error(chorale_resolve_signature(list(a = a, b = b)),
               class = "chorale_incompatible_phenotype_levels")
})








test_that("one-to-many mappings contribute every mapped identifier", {
  sets <- list(set_a = c("1"), set_b = c("2"))
  mapping <- data.frame(id = c("feature_a", "feature_a"),
                        ENTREZID = c("1", "2"), weight = c(0.5, 0.5))
  out <- chorale_geneset_matrix("feature_a", sets, mapping = mapping,
                                min_features = 1)
  expect_equal(unname(out[1, ]), c(0.5, 0.5))
  expect_equal(attr(out, "mapping_provenance")$one_to_many, 1L)
})
