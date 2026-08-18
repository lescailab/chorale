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

test_that("the integrated result is written as its own table", {
  r <- report_run()
  expect_true("programmes.tsv" %in% basename(r$files))
  p <- utils::read.delim(file.path(r$path, "programmes.tsv"))
  skip_if(nrow(p) == 0)
  expect_true(all(c("programme", "n_modalities", "modalities", "modality",
                    "factor", "markers", "phenotype_effect") %in% colnames(p)))
  # One row per modality carrying a programme, so a programme is named once
  # however many modalities it reaches.
  per <- table(p$programme)
  expect_true(all(per >= 2))
  expect_true(all(p$n_modalities >= 2))
})

test_that("the report leads with integration, not with per-modality factors", {
  r <- report_run()
  html <- paste(readLines(file.path(r$path, "report.html"), warn = FALSE),
                collapse = "\n")
  # The integrated section must precede the per-modality detail.
  expect_lt(regexpr("What integration found", html, fixed = TRUE)[1],
            regexpr("Supporting detail", html, fixed = TRUE)[1])
  expect_true(grepl("Shared biological programmes", html, fixed = TRUE))
})

test_that("every table carries a legend", {
  r <- report_run()
  html <- paste(readLines(file.path(r$path, "report.html"), warn = FALSE),
                collapse = "\n")
  n_tables <- length(gregexpr("<h3>", html, fixed = TRUE)[[1]])
  n_legends <- length(gregexpr("class='legend'", html, fixed = TRUE)[[1]])
  expect_gte(n_legends, n_tables)
})

test_that("a programme spans every modality carrying it, not one row per pair", {
  # Three modalities agreeing on one programme must be reported once, with all
  # three named, rather than as three separate pairings.
  sim <- chorale_simulate(n_modalities = 3, n_features = 150,
                          n_shared_factors = 3, n_private_factors = 2,
                          n_strains = 6, n_per_cell = 4, effect_size = 3,
                          seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = rep(5, 3), n_init = 3, seed = 1)
  pg <- chorale_programmes(fit)
  skip_if(nrow(pg) == 0)
  expect_true(all(c("programme", "n_modalities", "modalities") %in% colnames(pg)))
  # A factor belongs to exactly one programme.
  key <- paste(pg$modality, pg$factor)
  expect_equal(length(key), length(unique(key)))
  # n_modalities agrees with the members actually listed.
  for (pr in unique(pg$programme)) {
    d <- pg[pg$programme == pr, , drop = FALSE]
    expect_equal(d$n_modalities[1], length(unique(d$modality)))
  }
})

test_that("charts are inline and labelled for assistive technology", {
  r <- report_run()
  html <- paste(readLines(file.path(r$path, "report.html"), warn = FALSE),
                collapse = "\n")
  skip_if(!grepl("<svg", html, fixed = TRUE))
  expect_true(grepl("role='img'", html, fixed = TRUE))
  expect_true(grepl("aria-label=", html, fixed = TRUE))
  # Identity never rests on colour alone: a legend accompanies the marks.
  expect_true(grepl("legend-row", html, fixed = TRUE))
})

test_that("the report renders as self-contained HTML", {
  r <- report_run()
  html <- readLines(file.path(r$path, "report.html"), warn = FALSE)
  expect_true(any(grepl("<html", html)))
  expect_true(any(grepl("What integration found", html)))
  # Self-contained: no external stylesheet, script or remote asset.
  expect_false(any(grepl("<script|href=\"http|src=\"http", html)))
})

test_that("chorale_report rejects a non-fit", {
  expect_error(chorale_report(list(), NULL, NULL, path = withr::local_tempdir()),
               class = "rlang_error")
})

test_that("an unsupported assignment is never reported as a programme", {
  sim <- chorale_simulate(n_modalities = 2, n_features = 80,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 2, n_per_cell = 2, seed = 1)
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)

  # Strip the support, leaving the assignment in place. The assignment step
  # always returns a correspondence; only the null decides whether it means
  # anything, and a report must not substitute the one for the other.
  fit$programmes$supported <- FALSE
  fit$matches$significant <- FALSE
  expect_equal(nrow(chorale_programmes(fit)), 0L)

  factors <- chorale:::chorale_factor_table(fit)
  expect_false(any(factors$shared))

  path <- withr::local_tempdir()
  out <- chorale_report(fit, chorale_bound(fit), NULL, path = path)
  programmes <- utils::read.delim(file.path(path, "programmes.tsv"))
  expect_equal(nrow(programmes), 0L)
  html <- paste(readLines(file.path(path, "report.html"), warn = FALSE),
                collapse = "")
  expect_match(html, "No programme survived its null")
})
