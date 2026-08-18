fixture_path <- testthat::test_path("fixtures")

test_that("every modality fixture loads and is well formed", {
  skip_if_not(dir.exists(fixture_path))
  for (layer in c("RNA", "PROT", "METAB")) {
    fx <- chorale_fixture(layer, path = fixture_path)
    expect_true(is.matrix(fx$assay))
    expect_true(is.numeric(fx$assay))
    expect_gt(nrow(fx$assay), 0)
    expect_identical(colnames(fx$assay), fx$col_data$sample_id)
    expect_false(any(duplicated(rownames(fx$assay))))
  }
})

test_that("fixtures carry the required design metadata", {
  skip_if_not(dir.exists(fixture_path))
  for (layer in c("RNA", "PROT", "METAB")) {
    fx <- chorale_fixture(layer, path = fixture_path)
    expect_true(all(chorale_required_col_data() %in% colnames(fx$col_data)))
  }
})

test_that("fixtures build valid containers through chorale_load()", {
  skip_if_not(dir.exists(fixture_path))
  for (layer in c("RNA", "PROT", "METAB")) {
    fx <- chorale_fixture(layer, path = fixture_path)
    se <- chorale_load(fx$assay, fx$col_data)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(dim(se), dim(fx$assay))
  }
})

test_that("samples are disjoint across modalities", {
  skip_if_not(dir.exists(fixture_path))
  ids <- lapply(c("RNA", "PROT", "METAB"), function(l) {
    chorale_fixture(l, path = fixture_path)$col_data$sample_id
  })
  expect_length(intersect(ids[[1]], ids[[2]]), 0)
  expect_length(intersect(ids[[1]], ids[[3]]), 0)
  expect_length(intersect(ids[[2]], ids[[3]]), 0)
})

test_that("layers resolving covariates populate anchoring strata at least twice", {
  skip_if_not(dir.exists(fixture_path))
  for (layer in c("RNA", "METAB")) {
    cd <- chorale_fixture(layer, path = fixture_path)$col_data
    keys <- c("phenotype", "age_bin", "sex")
    complete <- cd[stats::complete.cases(cd[, keys, drop = FALSE]), , drop = FALSE]
    expect_gt(nrow(complete), 0)
    counts <- table(interaction(complete[, keys], drop = TRUE))
    expect_true(all(counts >= 2))
  }
})

test_that("the proteome fixture preserves a real missingness pattern", {
  skip_if_not(dir.exists(fixture_path))
  fx <- chorale_fixture("PROT", path = fixture_path)
  expect_true(is.numeric(fx$assay))
})

test_that("fixtures load quickly enough for continuous integration", {
  skip_if_not(dir.exists(fixture_path))
  elapsed <- system.time({
    for (layer in c("RNA", "PROT", "METAB")) {
      chorale_fixture(layer, path = fixture_path)
    }
  })[["elapsed"]]
  expect_lt(elapsed, 10)
})
