test_that("the registry names the default collections", {
  reg <- chorale_geneset_registry()
  expect_true(all(c("hallmark", "reactome", "cell_type", "kegg") %in% names(reg)))
  for (entry in reg) {
    expect_true(all(c("codes", "description") %in% names(entry)))
    expect_true(length(entry$codes) >= 1)
    expect_true("HS" %in% names(entry$codes))
  }
  # The three collections a factor is described in are registered under both
  # databases, so an analysis can be moved between them.
  for (nm in c("hallmark", "reactome", "cell_type")) {
    expect_true(all(c("MM", "HS") %in% names(reg[[nm]]$codes)))
  }
})

test_that("a collection registered under one database is refused for another", {
  # KEGG is distributed in the human collections and reaches other organisms
  # through orthologs, so asking for it mouse-native names the reason.
  expect_null(chorale_geneset_registry()$kegg$codes$MM)
  expect_error(chorale_genesets("kegg", db_species = "MM"), "no code registered")
})

test_that("the registry is exchangeable", {
  # A replacement registry changes the vocabulary without touching callers.
  alt <- list(only_hallmark = list(
    codes = list(MM = list(collection = "MH", subcollection = NULL)),
    description = "Hallmark alone"
  ))
  expect_error(
    chorale_genesets("reactome", registry = alt),
    class = "rlang_error"
  )
})

test_that("custom sets alone replace MSigDB entirely", {
  custom <- list(
    set_a = as.character(1:20),
    set_b = as.character(10:30)
  )
  sets <- chorale_genesets(collections = character(), custom = custom,
                           min_size = 5, max_size = 100)
  expect_named(sets, c("set_a", "set_b"))
  expect_identical(unname(attr(sets, "collection")), c("custom", "custom"))
})

test_that("sets outside the size window are dropped", {
  custom <- list(
    too_small = as.character(1:2),
    just_right = as.character(1:20),
    too_big = as.character(1:1000)
  )
  sets <- chorale_genesets(collections = character(), custom = custom,
                           min_size = 10, max_size = 500)
  expect_named(sets, "just_right")
})

test_that("an empty result is an error rather than a silent empty fit", {
  expect_error(
    chorale_genesets(collections = character(),
                     custom = list(a = as.character(1:3)),
                     min_size = 100),
    class = "rlang_error"
  )
})

test_that("the indicator matrix marks membership and drops thin sets", {
  sets <- list(kept = c("1", "2", "3"), thin = c("4"))
  m <- chorale_geneset_matrix(c("1", "2", "3", "4"), sets, min_features = 2)
  expect_equal(colnames(m), "kept")
  expect_equal(as.numeric(m[, "kept"]), c(1, 1, 1, 0))
})

test_that("fractional weights carry into the indicator matrix", {
  sets <- list(kept = c("1", "2"))
  m <- chorale_geneset_matrix(c("1", "2"), sets, weights = c(0.5, 1),
                              min_features = 1)
  expect_equal(as.numeric(m[, "kept"]), c(0.5, 1))
})

test_that("mismatched weights are rejected", {
  expect_error(
    chorale_geneset_matrix(c("1", "2"), list(a = "1"), weights = 1),
    class = "rlang_error"
  )
})

test_that("MSigDB collections retrieve for mouse and for human", {
  skip_if_not_installed("msigdbr")
  skip_on_cran()
  mouse <- chorale_genesets("hallmark", species = "Mus musculus",
                            db_species = "MM")
  expect_gt(length(mouse), 20)
  expect_true(all(vapply(mouse, is.character, logical(1))))
  # Swapping species and source database serves a human analysis unchanged.
  human <- chorale_genesets("hallmark", species = "Homo sapiens",
                            db_species = "HS")
  expect_gt(length(human), 20)
})
