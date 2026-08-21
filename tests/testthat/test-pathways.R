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






test_that("a feature space the vocabulary cannot reach is declared, not guessed", {
  lipids <- c(paste0("SM ", 30:44, ":1;O2"), paste0("Cer ", 30:44, ":1;O2"))
  mk <- function(tag, ids) {
    m <- matrix(stats::rnorm(length(ids) * 20), length(ids),
                dimnames = list(ids, paste0(tag, 1:20)))
    chorale_load(m, data.frame(sample_id = paste0(tag, 1:20),
                               phenotype = rep(c("case", "control"), each = 10),
                               stringsAsFactors = FALSE))
  }
  genes <- paste0("g", 1:30)
  containers <- list(RNA = mk("r", genes), METAB = mk("m", lipids))
  sets <- stats::setNames(list(genes[1:20]),
                          "REACTOME_SPHINGOLIPID_METABOLISM")

  # Undeclared, the lipidome reaches nothing: its identifiers are not genes.
  silent <- chorale_concepts(containers, sets, min_features = 5)
  expect_equal(ncol(silent$membership$METAB), 0L)

  # Declared, it reaches the same concept through its classes.
  declared <- chorale_concepts(containers, sets,
                               feature_space = c(METAB = "lipid"),
                               min_features = 5)
  expect_equal(ncol(declared$membership$METAB), 1L)
  expect_error(chorale_concepts(containers, sets,
                                feature_space = c(METAB = "protein")),
               "must be")
})
