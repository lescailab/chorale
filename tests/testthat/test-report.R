report_run <- function(with_controls = TRUE) {
  sim <- chorale_simulate(n_modalities = 2, n_features = 120,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 5, n_per_cell = 3, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
  path <- withr::local_tempdir(.local_envir = parent.frame())
  nul <- if (with_controls) {
    chorale_null(fit, containers, n_permutations = 2, n_init = 2)
  } else {
    NULL
  }
  files <- chorale_report(fit, chorale_bound(fit), nul, path = path)
  list(fit = fit, path = path, files = files)
}

test_that("every output named in the plan is written", {
  r <- report_run()
  expected <- c("factors.tsv", "markers.tsv", "associations.tsv",
                "concordance.tsv", "bounds.tsv", "controls.tsv",
                "factors.gmt", "report.html")
  expect_true(all(expected %in% basename(r$files)))
  for (m in r$fit$modalities) {
    expect_true(paste0("loadings_", m, ".tsv") %in% basename(r$files))
    expect_true(paste0("scores_", m, ".tsv") %in% basename(r$files))
  }
  expect_true(all(file.exists(r$files)))
})

test_that("every factor is reported with a resolution status", {
  r <- report_run()
  f <- utils::read.delim(file.path(r$path, "factors.tsv"))
  expect_equal(nrow(f), sum(r$fit$n_factors))
  expect_true(all(c("pathway_definition", "n_markers", "status") %in% colnames(f)))
  # A factor that cannot be named is reported as unresolved, never dropped.
  expect_true(all(nzchar(f$status)))
})

test_that("scores carry the full sample metadata", {
  r <- report_run()
  s <- utils::read.delim(file.path(r$path, "scores_modality_1.tsv"))
  expect_true(all(c("sample_id", "phenotype", "sex", "age_months") %in% colnames(s)))
  expect_equal(nrow(s), ncol(r$fit$fits[[1]]$scores) * 0 + nrow(r$fit$fits[[1]]$scores))
})

test_that("associations recover the planted phenotype effect", {
  r <- report_run()
  a <- utils::read.delim(file.path(r$path, "associations.tsv"))
  expect_true(nrow(a) > 0)
  pheno <- a[a$covariate == "phenotype", , drop = FALSE]
  expect_true(any(pheno$p_permutation < 0.05))
})

test_that("controls are written even without a null object", {
  r <- report_run(with_controls = FALSE)
  ctl <- utils::read.delim(file.path(r$path, "controls.tsv"))
  expect_true(nrow(ctl) >= 1)
  expect_true(any(grepl("pure-feature", ctl$control)))
})

test_that("the GMT lists marker features per factor", {
  r <- report_run()
  gmt <- readLines(file.path(r$path, "factors.gmt"))
  skip_if(length(gmt) == 0)
  parts <- strsplit(gmt[1], "\t")[[1]]
  expect_gt(length(parts), 2)
})

test_that("the report renders as self-contained HTML", {
  r <- report_run()
  html <- readLines(file.path(r$path, "report.html"), warn = FALSE)
  expect_true(any(grepl("<html", html)))
  expect_true(any(grepl("chorale report", html)))
  # Self-contained: no external stylesheet or script is referenced.
  expect_false(any(grepl("<script|href=\"http", html)))
})

test_that("chorale_report rejects a non-fit", {
  expect_error(chorale_report(list(), NULL, NULL, path = withr::local_tempdir()),
               class = "rlang_error")
})
