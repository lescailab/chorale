#' Whether a collection meets the conditions the method rests on
#'
#' The estimator recovers a shared state from modalities measured on different
#' animals, and it can only do so where four conditions hold. This function
#' evaluates them on a collection before it is fitted, so a collection that
#' cannot support the estimand is refused rather than fitted and read.
#'
#' \describe{
#'   \item{non-Gaussianity}{Are a modality's recovered components further from
#'     normal than components recovered the same way from data Gaussian by
#'     construction? Independent component analysis maximises non-Gaussianity,
#'     so the question is only meaningful against that calibration.}
#'   \item{modality difference}{Do the component distributions of two
#'     modalities differ from one another? Identification rests on the
#'     modalities not being copies of one another in distribution, so agreement
#'     here is a failure.}
#'   \item{detectability}{How many components stand above what the same matrix
#'     produces with its covariance destroyed, and above the spiked-covariance
#'     threshold implied by its shape?}
#'   \item{anchor richness}{How many design strata are populated, and how many
#'     of those are shared across the modalities? The design is what stands in
#'     for matched individuals: no animal appears in two modalities, so the
#'     comparison is anchored on the strata the modalities have in common.}
#' }
#'
#' The two distributional conditions are evaluated in Python, which supplies
#' the Anderson-Darling and Kolmogorov-Smirnov machinery, but the components
#' they are evaluated on are recovered by this package's own estimator, passed
#' across as a callback. The estimator deciding whether the conditions hold is
#' therefore the estimator that will consume them.
#'
#' @param containers A named list of `SummarizedExperiment` objects, as
#'   [chorale_load()] returns, or of feature-by-sample matrices.
#' @param designs A named list of design tables. Required only where
#'   `containers` holds bare matrices; taken from the containers otherwise.
#' @param control A [chorale_control()] object. `alpha`, `n_init`,
#'   `n_factors_quantile` and `max_factors` are read from it.
#' @param n_factors Optional named integer vector fixing the component count
#'   per modality. Taken from the detectability condition where absent.
#' @param n_surrogate Integer number of Gaussian surrogates the non-Gaussianity
#'   condition is calibrated against.
#' @param n_perm Integer number of permutations for the detectability
#'   condition.
#' @param seed Integer seed.
#'
#' @returns An object of class `chorale_gates`: a list of one data frame per
#'   condition, plus `n_factors`, the component count each modality supports.
#'
#' @examples
#' \dontrun{
#' sim <- chorale_simulate(n_modalities = 2, n_samples = 60, n_features = 200)
#' containers <- lapply(sim$modalities, chorale_load, design = sim$design)
#' chorale_gates(containers)
#' }
#' @export
chorale_gates <- function(containers, designs = NULL,
                          control = chorale_control(),
                          n_factors = NULL, n_surrogate = 100L,
                          n_perm = 200L, seed = 1L) {
  if (!is.list(containers) || length(containers) < 1) {
    rlang::abort("`containers` must be a named list of modalities.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }

  assays <- lapply(containers, function(x) {
    if (inherits(x, "SummarizedExperiment")) SummarizedExperiment::assay(x) else as.matrix(x)
  })
  if (is.null(designs)) {
    designs <- lapply(containers, function(x) {
      if (!inherits(x, "SummarizedExperiment")) {
        rlang::abort("`designs` is required where `containers` holds matrices.")
      }
      as.data.frame(SummarizedExperiment::colData(x))
    })
  }

  # Samples by features, centred and scaled, which is what both the estimator
  # and the surrogate construction expect.
  xs <- lapply(assays, function(a) {
    x <- t(as.matrix(a))
    x[!is.finite(x)] <- 0
    x <- scale(x)
    x[, apply(x, 2, function(col) all(is.finite(col))), drop = FALSE]
  })

  detect <- chorale_gate_detectability(xs, n_perm = n_perm,
                                       quantile = control$n_factors_quantile,
                                       max_factors = control$max_factors,
                                       seed = seed)
  if (is.null(n_factors)) {
    n_factors <- stats::setNames(detect$n_factors, detect$modality)
  }

  ica_fn <- chorale_gate_ica(n_init = control$n_init, consensus = control$consensus)
  py <- chorale_gates_python()

  ng <- lapply(names(xs), function(m) {
    py$gate_nongaussianity(xs[[m]], as.integer(n_factors[[m]]), m, ica_fn,
                           seed = as.integer(seed),
                           n_surrogate = as.integer(n_surrogate),
                           alpha = control$alpha)
  })
  names(ng) <- names(xs)

  sources <- lapply(names(xs), function(m) {
    ica_fn(xs[[m]], as.integer(n_factors[[m]]), as.integer(seed))
  })
  names(sources) <- names(xs)
  difference <- py$gate_modality_difference(sources, alpha = control$alpha)

  out <- list(
    non_gaussianity = chorale_records_to_df(lapply(ng, function(r) r[[1]])),
    non_gaussianity_components = chorale_records_to_df(
      unlist(lapply(ng, function(r) r[[2]]), recursive = FALSE)),
    modality_difference = chorale_records_to_df(difference),
    detectability = detect,
    anchor_richness = chorale_gate_anchors(designs),
    n_factors = n_factors
  )
  class(out) <- "chorale_gates"
  out
}


#' The component estimator, wrapped for the gate module to call
#'
#' Returns a function of `(x, k, seed)` giving standardised sources, so the
#' components the conditions are judged on are recovered by the same estimator
#' that will consume them.
#' @keywords internal
#' @noRd
chorale_gate_ica <- function(n_init = 20L, consensus = TRUE) {
  function(x, k, seed) {
    x <- as.matrix(x)
    fit <- chorale_ica(x, n_factors = as.integer(k), n_init = n_init,
                       seed = as.integer(seed), consensus = consensus)
    unname(as.matrix(fit$scores))
  }
}


#' The detectability condition
#'
#' Two readings of the same question. Parallel analysis compares each
#' eigenvalue with what the matrix produces once each feature is permuted
#' independently, destroying covariance while preserving every marginal. The
#' spiked-covariance threshold asks the same of the shape of the matrix alone,
#' at `sigma^2 (1 + sqrt(gamma))` with `sigma^2` taken as the median bulk
#' eigenvalue.
#' @keywords internal
#' @noRd
chorale_gate_detectability <- function(xs, n_perm = 200L, quantile = 0.95,
                                       max_factors = 20L, seed = 1L) {
  rows <- lapply(names(xs), function(m) {
    x <- xs[[m]]
    n <- nrow(x)
    p <- ncol(x)
    ev <- svd(x, nu = 0, nv = 0)$d^2 / max(n - 1L, 1L)
    sigma2 <- stats::median(ev)
    gamma <- p / n
    threshold <- sigma2 * (1 + sqrt(gamma))
    k <- chorale_n_factors(x, n_perm = n_perm, quantile = quantile,
                           max_factors = max_factors, seed = seed)
    data.frame(
      modality = m,
      n_samples = n,
      n_features = p,
      gamma_p_over_n = round(gamma, 3),
      n_spikes_above_threshold = sum(ev > threshold),
      n_factors = k,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' The anchor-richness condition
#'
#' Reports, at each level of coarsening, how many design strata each modality
#' populates and how many of those every modality shares. Nothing beyond
#' `sample_id` and `phenotype` is required of a design, so the coarsening is
#' built from whatever covariates the designs actually carry.
#' @keywords internal
#' @noRd
chorale_gate_anchors <- function(designs) {
  designs <- lapply(designs, as.data.frame)
  shared_cols <- Reduce(intersect, lapply(designs, colnames))
  candidates <- setdiff(shared_cols, c("sample_id", "modality", "cohort"))
  varies <- vapply(candidates, function(cv) {
    all(vapply(designs, function(d) {
      length(unique(stats::na.omit(d[[cv]]))) >= 2
    }, logical(1)))
  }, logical(1))
  candidates <- candidates[varies]
  if (!"phenotype" %in% candidates) {
    rlang::abort("No shared covariate varies in every modality; nothing can be anchored on.")
  }
  # Coarsen from every shared covariate down to the phenotype alone, so a
  # collection that cannot support the finest anchoring is still reported at
  # the level it can support.
  others <- setdiff(candidates, "phenotype")
  levels_of <- lapply(seq(length(others), 0), function(k) {
    c("phenotype", utils::head(others, k))
  })

  rows <- list()
  for (keys in levels_of) {
    label <- paste(keys, collapse = " x ")
    per_mod <- lapply(designs, function(d) {
      sub <- d[stats::complete.cases(d[, keys, drop = FALSE]), keys, drop = FALSE]
      if (nrow(sub) == 0) return(character())
      unique(do.call(paste, c(sub, sep = "|")))
    })
    shared <- Reduce(intersect, per_mod)
    for (m in names(designs)) {
      rows[[length(rows) + 1L]] <- data.frame(
        coarsening = label, modality = m,
        n_strata = length(per_mod[[m]]),
        n_strata_shared = length(shared),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


#' Turn the gate module's records into a data frame
#'
#' The Python side returns lists of named scalars rather than a data frame, so
#' every value arrives as something R can compare. A record missing a field is
#' filled rather than dropped, so one modality cannot silently shorten a table.
#' @keywords internal
#' @noRd
chorale_records_to_df <- function(records) {
  if (is.null(records) || length(records) == 0) return(data.frame())
  if (!is.null(names(records)) && !is.list(records[[1]])) records <- list(records)
  records <- unname(records)
  cols <- unique(unlist(lapply(records, names)))
  out <- lapply(cols, function(cn) {
    vals <- lapply(records, function(r) {
      v <- r[[cn]]
      if (is.null(v) || length(v) == 0) NA else v[[1]]
    })
    if (all(vapply(vals, function(v) is.numeric(v) || is.na(v), logical(1)))) {
      as.numeric(unlist(vals))
    } else {
      as.character(unlist(lapply(vals, as.character)))
    }
  })
  names(out) <- cols
  as.data.frame(out, stringsAsFactors = FALSE)
}


#' Load the gate module into the bound Python session
#' @keywords internal
#' @noRd
chorale_gates_python <- function() {
  rlang::check_installed("reticulate")
  path <- system.file("python", package = "chorale")
  if (!nzchar(path)) {
    rlang::abort("Could not locate inst/python in the installed package.")
  }
  reticulate::import_from_path("gates", path = path, delay_load = FALSE)
}


#' @export
print.chorale_gates <- function(x, ...) {
  cat("<chorale_gates>\n")
  cat("\nnon-Gaussianity\n")
  print(x$non_gaussianity[, c("modality", "median_A2_observed",
                              "median_A2_surrogate", "p_value", "verdict")],
        row.names = FALSE)
  cat("\nmodality difference\n")
  print(x$modality_difference[, c("pair", "pooled_KS_D", "pooled_KS_p",
                                  "pct_pairs_indistinguishable", "verdict")],
        row.names = FALSE)
  cat("\ndetectability\n")
  print(x$detectability, row.names = FALSE)
  cat("\nanchor richness\n")
  print(utils::head(x$anchor_richness, 12), row.names = FALSE)
  invisible(x)
}
