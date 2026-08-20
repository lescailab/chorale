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

test_that("secondary categorical variables receive equal block weight", {
  profile <- list(
    z = matrix(c(2, 3, 6, 9), nrow = 1,
               dimnames = list("factor_1",
                               c("phenotype=case", "site=b", "site=c", "age"))),
    term_covariate = c("phenotype=case" = "phenotype", "site=b" = "site",
                       "site=c" = "site", age = "age"))
  spec <- list(secondary = c("site", "age"))
  out <- chorale_secondary_signature(profile, spec)
  expect_equal(unname(out[1, 1:2]), c(3, 6) / sqrt(2))
  expect_equal(unname(out[1, 3]), 9)
})

test_that("secondary affinity cannot promote phenotype-incompatible factors", {
  effects_a <- matrix(c(2, 4, 1, 10), 2, byrow = TRUE,
                      dimnames = list(c("a1", "a2"),
                                      c("phenotype=case", "age")))
  effects_b <- matrix(c(2.1, -10, 8, 10), 2, byrow = TRUE,
                      dimnames = list(c("b1", "b2"),
                                      c("phenotype=case", "age")))
  make <- function(e) list(effects = e, se = matrix(0.2, 2, 2,
    dimnames = dimnames(e)), z = e / 0.2,
    term_covariate = c("phenotype=case" = "phenotype", age = "age"))
  spec <- list(phenotype = "phenotype", secondary = "age")
  blocks <- chorale_hierarchical_blocks(list(a = make(effects_a),
                                              b = make(effects_b)), spec,
                                        ambiguity_level = 0.95)
  expect_true(blocks$diagnostics[["a|b"]]$candidate[1, 1])
  expect_false(blocks$diagnostics[["a|b"]]$candidate[1, 2])
  expect_gt(blocks$final[["a|b"]][1, 1], blocks$final[["a|b"]][1, 2])
})

test_that("multi-level ambiguity cutoff is on the mean-loss scale", {
  effects <- matrix(c(2, 1, 0.5, -0.5), 2, byrow = TRUE,
                    dimnames = list(c("f1", "f2"),
                                    c("phenotype=b", "phenotype=c")))
  make <- function(e) list(
    effects = e, se = matrix(1, 2, 2, dimnames = dimnames(e)), z = e,
    covariance = list(diag(2), diag(2)),
    term_covariate = c("phenotype=b" = "phenotype",
                       "phenotype=c" = "phenotype"))
  spec <- list(phenotype = "phenotype", secondary = character())
  out <- chorale_hierarchical_blocks(list(a = make(effects), b = make(effects)),
                                      spec, ambiguity_level = 0.95)
  expect_equal(out$diagnostics[["a|b"]]$ambiguity_cutoff,
               stats::qchisq(0.95, df = 2) / 2)
})

test_that("secondary evidence refines rather than replaces phenotype affinity", {
  primary <- matrix(c(0.9, 0.1), 1, 2)
  secondary <- matrix(c(0, 0.9), 1, 2)
  candidate <- matrix(TRUE, 1, 2)
  out <- chorale_combine_affinity(primary, secondary, candidate)
  expect_gt(out[1, 1], out[1, 2])
  bounded <- primary / (1 + primary)
  expect_equal(out, 1 + bounded * (1 + secondary))

  candidate[1, 2] <- FALSE
  out <- chorale_combine_affinity(primary, secondary, candidate)
  expect_gt(out[1, 1], out[1, 2])
  expect_equal(out[1, 2], bounded[1, 2])
})

test_that("effect affinity is bounded and requires support plus compatibility", {
  make <- function(effect) list(
    effects = matrix(effect, 1, dimnames = list("f1", "phenotype=case")),
    se = matrix(1, 1, 1), covariance = list(matrix(1, 1, 1)),
    term_covariate = c("phenotype=case" = "phenotype"))
  strong_agree <- chorale_effect_block_affinity(
    make(4), make(4), 1, 1, 1, 1)$affinity
  strong_disagree <- chorale_effect_block_affinity(
    make(4), make(1), 1, 1, 1, 1, fixed_orientation = 1)$affinity
  weak_agree <- chorale_effect_block_affinity(
    make(0.1), make(0.1), 1, 1, 1, 1)$affinity
  expect_true(all(c(strong_agree, strong_disagree, weak_agree) >= 0 &
                    c(strong_agree, strong_disagree, weak_agree) <= 1))
  expect_gt(strong_agree, strong_disagree)
  expect_gt(strong_agree, weak_agree)
})

test_that("resampling streams have distinct deterministic seeds", {
  a <- chorale_resampling_seeds(7, c(10, 20, 30))
  b <- chorale_resampling_seeds(7, c(10, 20, 30))
  expect_identical(a, b)
  expect_length(unique(unlist(a)), 60)
})

test_that("bootstrap candidates retain phenotype-loss ties", {
  make <- function(loss) list(
    primary = list("a|b" = matrix(1, 1, 2)),
    secondary = list("a|b" = matrix(c(0.2, 0.8), 1, 2)),
    final = list("a|b" = matrix(1, 1, 2)),
    diagnostics = list("a|b" = list(
      loss = matrix(loss, 1, 2), candidate = matrix(TRUE, 1, 2))))
  observed <- make(c(1, 1.1))
  draws <- list(make(c(1.1, 1)), make(c(1, 1.2)), make(c(1.2, 1)))
  out <- chorale_bootstrap_candidates(observed, draws, 0.95)
  expect_true(all(out$diagnostics[["a|b"]]$candidate))
  expect_true(all(is.finite(out$diagnostics[["a|b"]]$phenotype_margin_lower)))
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
