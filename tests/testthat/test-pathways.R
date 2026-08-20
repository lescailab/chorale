test_that("the shipped metabolite table is well formed", {
  map <- chorale_metabolite_pathways()
  needed <- c("abbreviation", "lipidmaps_category", "chebi_id", "chebi_name",
              "pathway_id", "pathway_name", "msigdb_name", "n_compounds",
              "n_pathway_compounds", "specificity")
  expect_true(all(needed %in% colnames(map)))
  expect_gt(nrow(map), 100)
  expect_false(any(is.na(map$abbreviation)))
  expect_true(all(grepl("^CHEBI:[0-9]+$", map$chebi_id)))
  expect_true(all(grepl("^R-MMU-[0-9]+$", map$pathway_id)))
  # The pathway names are the gene-side vocabulary, which is what puts a lipid
  # and a gene in the same set rather than in two that have to be reconciled.
  expect_true(all(grepl("^REACTOME_", map$msigdb_name)))
  expect_true(all(map$specificity > 0 & map$specificity <= 1))
  expect_true(all(map$n_compounds <= map$n_pathway_compounds))
})

test_that("the table places the classes it should where it should", {
  map <- chorale_metabolite_pathways()
  of <- function(cls) map$msigdb_name[map$abbreviation == cls &
                                        (map$n_compounds >= 2 |
                                           map$specificity >= 0.05)]
  expect_true("REACTOME_CARNITINE_SHUTTLE" %in% of("CAR"))
  expect_true("REACTOME_SPHINGOLIPID_METABOLISM" %in% of("SM"))
  expect_true("REACTOME_SPHINGOLIPID_METABOLISM" %in% of("Cer"))
  # A class reached through one generic compound inside a large pathway is not
  # evidence about that pathway, and the specificity filter is what removes it.
  expect_false("REACTOME_BIOLOGICAL_OXIDATIONS" %in% of("CE"))
  expect_false("REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" %in% of("CAR"))
  # A single compound still counts where the pathway is small enough for one to
  # matter, which is how a class with one representative keeps its own pathway.
  expect_true("REACTOME_CHYLOMICRON_ASSEMBLY" %in% of("CE"))
})

test_that("a lipid shorthand name yields its class", {
  expect_equal(
    chorale_lipid_class(c("PC 34:1|[M+H]+", "SM 36:1;O2", "CAR 18:2")),
    c("PC", "SM", "CAR")
  )
  # A dialect the grammar rejects still resolves, because the class is the
  # leading token in every dialect the project has met.
  expect_equal(chorale_lipid_class("AHexCer 57:1;O3|AHexCer (O-16:0)18:1"),
               "AHexCer")
  expect_true(is.na(chorale_lipid_class("")))
})

test_that("a lipidome lands in the same sets as the genes", {
  sets <- stats::setNames(vector("list", 3),
                          c("REACTOME_SPHINGOLIPID_METABOLISM",
                            "REACTOME_CARNITINE_SHUTTLE",
                            "REACTOME_NOT_A_REAL_SET"))
  features <- c(paste0("SM ", 30:36, ":1;O2"), paste0("Cer ", 30:34, ":1;O2"),
                paste0("CAR ", 16:20, ":0"))
  m <- chorale_metabolite_matrix(features, sets, min_features = 1)
  expect_true("REACTOME_SPHINGOLIPID_METABOLISM" %in% colnames(m))
  expect_true("REACTOME_CARNITINE_SHUTTLE" %in% colnames(m))
  # A set the vocabulary does not contain cannot appear.
  expect_false("REACTOME_NOT_A_REAL_SET" %in% colnames(m))
  expect_equal(sum(m[, "REACTOME_CARNITINE_SHUTTLE"]), 5)
  expect_equal(nrow(m), length(features))
})

test_that("a pathway profile is signed and scaled by set size", {
  set.seed(1)
  ids <- paste0("f", 1:60)
  l <- matrix(stats::rnorm(120), nrow = 60,
              dimnames = list(ids, c("factor_1", "factor_2")))
  p <- matrix(0, nrow = 60, ncol = 2,
              dimnames = list(ids, c("up", "down")))
  p[1:20, "up"] <- 1
  p[21:40, "down"] <- 1
  l[1:20, "factor_1"] <- l[1:20, "factor_1"] + 3
  l[21:40, "factor_1"] <- l[21:40, "factor_1"] - 3

  prof <- chorale_pathway_profile(l, p)
  expect_equal(dim(prof), c(2L, 2L))
  expect_gt(prof["factor_1", "up"], 3)
  expect_lt(prof["factor_1", "down"], -3)
  # A factor with no relation to either set sits near zero on both.
  expect_lt(abs(prof["factor_2", "up"]), 3)
})

test_that("the pathway channel fires on agreement and not on disagreement", {
  build <- function(same) {
    ids <- paste0("f", 1:60)
    sets <- list(set_a = ids[1:20], set_b = ids[21:40])
    prior <- chorale_geneset_matrix(ids, sets, min_features = 5)
    mk <- function(shift_set) {
      set.seed(2)
      l <- matrix(stats::rnorm(180), nrow = 60,
                  dimnames = list(ids, paste0("factor_", 1:3)))
      l[shift_set, "factor_1"] <- l[shift_set, "factor_1"] + 4
      s <- matrix(stats::rnorm(120), nrow = 40,
                  dimnames = list(paste0("sample_", 1:40),
                                  paste0("factor_", 1:3)))
      x <- cbind(1, s) %*% rbind(rep(0, 60), t(l))
      colnames(x) <- ids
      list(loadings = l, prior = prior, scores = s, analysis_matrix = x)
    }
    designs <- list(
      a = data.frame(sample_id = paste0("sample_", 1:40)),
      b = data.frame(sample_id = paste0("sample_", 1:40))
    )
    fit <- structure(
      list(
        modalities = c("a", "b"),
        fits = list(a = mk(1:20), b = mk(if (same) 1:20 else 21:40)),
        designs = designs,
        programmes = data.frame(
          programme = c("P1", "P1"), modality = c("a", "b"),
          factor = c("factor_1", "factor_1"), n_modalities = 2L,
          modalities = "a, b", joint_statistic = 1, joint_p = 0.01,
          supported = TRUE, stringsAsFactors = FALSE
        )
      ),
      class = "chorale_fit"
    )
    chorale_pathway_evidence(fit, n_perm = 199L)
  }

  agree <- build(TRUE)
  disagree <- build(FALSE)
  expect_equal(nrow(agree), 1L)
  # Loading the same sets is corroboration; loading different ones is not, and
  # the statistic has to separate them or the channel says nothing.
  expect_gt(agree$pathway_statistic, disagree$pathway_statistic)
  expect_true(agree$pathway_supported)
  expect_true(all(agree$pathway_p >= 1 / 200))
})

test_that("evidence labels say which channels a programme rests on", {
  pg <- data.frame(programme = c("P1", "P2"), supported = c(TRUE, TRUE),
                   stringsAsFactors = FALSE)
  pe <- data.frame(programme = c("P1", "P2"),
                   pathway_statistic = c(9, 1), pathway_p = c(0.001, 0.6),
                   n_shared_sets = 5L, pathway_supported = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)
  out <- chorale_evidence_label(pg, pe)
  expect_equal(out$evidence, c("design and pathway", "design only"))

  none <- chorale_evidence_label(pg, data.frame())
  expect_equal(none$evidence, c("design only", "design only"))
  expect_false(any(none$pathway_supported))
})

test_that("the curated projection annotates the fit without replacing it", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  ids <- lapply(containers, function(se)
    rownames(SummarizedExperiment::assay(se)))
  span <- function(a, b) unlist(lapply(ids, function(v) v[a:b]))
  sets <- list(set_a = span(1, 40), set_b = span(30, 80), set_c = span(60, 120))
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
                     gene_sets = sets, n_pathway_perm = 0)

  f <- fit$fits[[1]]
  expect_false(is.null(f$pathway_loadings))
  # The stored loadings are the ones the scores were fitted with, so the fit
  # still reconstructs the data; the projection is kept beside them.
  expect_false(isTRUE(all.equal(f$loadings, f$pathway_loadings)))
  expect_gt(f$reconstruction$fitted, f$reconstruction$projected)
  expect_gt(f$reconstruction$fitted, 0.5)
})

test_that("an undeclared feature space is refused", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  expect_error(
    chorale_fit(containers, n_factors = 2, n_init = 2,
                feature_space = c(modality_1 = "protein")),
    class = "rlang_error"
  )
  expect_error(
    chorale_fit(containers, n_factors = 2, n_init = 2,
                feature_space = c(nonexistent = "lipid")),
    class = "rlang_error"
  )
})
