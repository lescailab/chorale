#' Every decision the estimator makes, in one object
#'
#' A run is decided by more than its data: the threshold a programme must beat,
#' how many permutations calibrate it, what counts as a pure feature, how far a
#' lipid class must reach into a pathway before the link is used. Left inside
#' the code, those choices are invisible to the reader of a result and
#' unchangeable by anyone analysing different data. Gathered here they are
#' arguments, they travel with the fit that used them, and an analysis plan
#' fixed before a run is a call to this function rather than a description of
#' one.
#'
#' Every value has a default, so a caller who does not care sets nothing. A
#' caller who does care sets only what they mean to change, and the result
#' records the whole set.
#'
#' @param alpha Significance threshold a programme's joint evidence must beat,
#'   and the level at which false discovery is controlled.
#' @param n_perm Permutations calibrating the design channel. The smallest
#'   attainable p-value is `1 / (n_perm + 1)`, so a threshold the count cannot
#'   reach is not reportable.
#' @param n_pathway_perm Score-residual permutations calibrating the
#'   pathway channel. Zero skips the channel.
#' @param n_init Initialisations of the factorisation per modality.
#' @param consensus Deprecated compatibility switch. Factors are always taken
#'   from the medoid stable run; aligned-run averaging is not an ICA solution.
#' @param require_pure_features Optional diagnostic filter: report only
#'   programmes whose every member passes the requested loading-purity rule.
#'   This is not an identification condition and defaults to `FALSE`.
#' @param purity_ratio A feature is pure for a factor when its largest
#'   competing loading is at most this fraction of its own.
#' @param min_markers Features a factor needs to pass the optional purity
#'   diagnostic.
#' @param max_markers Markers retained per factor.
#' @param lambda Ridge penalty on the curated-set coefficients.
#' @param min_set_features Features of a modality a curated set must contain
#'   before it can describe that modality.
#' @param min_lipid_compounds Distinct compounds of a lipid class a pathway must
#'   be annotated with before the link is used.
#' @param min_lipid_specificity Share of a pathway's annotated compounds a lipid
#'   class must account for before the link is used, as an alternative to the
#'   count.
#' @param n_factors_quantile Quantile of the permuted eigenvalues a component
#'   must exceed when the factor count is chosen from the data.
#' @param max_factors Upper bound on a factor count chosen from the data.
#' @param n_grid Quantiles representing each marginal when bounding a coupling.
#' @param phenotype_column Name of the mandatory phenotype column shared by all
#'   modalities.
#' @param phenotype_reference Reference phenotype level. The default is
#'   `"control"`; an informative error is raised when that level is absent.
#' @param profile_covariates Optional covariates allowed to refine phenotype-led
#'   matching. `NULL` discovers every eligible covariate shared by all
#'   modalities. Phenotype is always included and always remains primary.
#' @param bound_strata Covariates used to condition coupling bounds. `NULL`
#'   uses the eligible shared design covariates. This is deliberately separate
#'   from the variables used for matching.
#' @param exchangeability_blocks Optional shared design columns within which
#'   permutation and bootstrap resampling must stay.
#' @param phenotype_alpha Family-wise significance threshold for phenotype
#'   support. `NULL` uses `alpha`.
#' @param ambiguity_level Confidence level for assignment-margin intervals.
#' @param n_ambiguity_boot Bootstrap replicates used to assess whether a
#'   phenotype-compatible assignment is separated from its alternatives. The
#'   default gives about 25 draws in each 2.5 per cent tail of a 95 per cent
#'   interval. Zero uses an analytic margin diagnostic without resampling.
#' @param n_cores Worker processes for bootstrap and permutation loops. Results
#'   are deterministic across worker counts. Windows uses a sequential fallback.
#'
#' @returns A named list of the decisions, of class `chorale_control`.
#' @export
#' @examples
#' # The defaults.
#' str(chorale_control()[c("alpha", "n_perm", "purity_ratio")])
#'
#' # An analysis plan fixed before a run is a call, not a paragraph.
#' plan <- chorale_control(alpha = 0.01, n_perm = 5000L)
#' plan$alpha
chorale_control <- function(alpha = 0.05,
                            n_perm = 200L,
                            n_pathway_perm = 200L,
                            n_init = 20L,
                            consensus = TRUE,
                            require_pure_features = FALSE,
                            purity_ratio = 0.25,
                            min_markers = 2L,
                            max_markers = 20L,
                            lambda = 1,
                            min_set_features = 5L,
                            min_lipid_compounds = 2L,
                            min_lipid_specificity = 0.05,
                            n_factors_quantile = 0.95,
                            max_factors = 20L,
                            n_grid = 200L,
                            phenotype_column = "phenotype",
                            phenotype_reference = "control",
                            profile_covariates = NULL,
                            bound_strata = NULL,
                            exchangeability_blocks = NULL,
                            phenotype_alpha = NULL,
                            ambiguity_level = 0.95,
                            n_ambiguity_boot = 999L,
                            n_cores = 1L) {
  phenotype_alpha <- phenotype_alpha %||% alpha
  out <- list(
    alpha = alpha, n_perm = as.integer(n_perm),
    n_pathway_perm = as.integer(n_pathway_perm), n_init = as.integer(n_init),
    consensus = isTRUE(consensus),
    require_pure_features = isTRUE(require_pure_features),
    purity_ratio = purity_ratio, min_markers = as.integer(min_markers),
    max_markers = as.integer(max_markers), lambda = lambda,
    min_set_features = as.integer(min_set_features),
    min_lipid_compounds = as.integer(min_lipid_compounds),
    min_lipid_specificity = min_lipid_specificity,
    n_factors_quantile = n_factors_quantile,
    max_factors = as.integer(max_factors), n_grid = as.integer(n_grid),
    phenotype_column = as.character(phenotype_column),
    phenotype_reference = as.character(phenotype_reference),
    profile_covariates = profile_covariates,
    bound_strata = bound_strata,
    exchangeability_blocks = exchangeability_blocks,
    phenotype_alpha = phenotype_alpha,
    ambiguity_level = ambiguity_level,
    n_ambiguity_boot = as.integer(n_ambiguity_boot),
    n_cores = as.integer(n_cores)
  )
  if (out$alpha <= 0 || out$alpha >= 1) {
    rlang::abort("`alpha` must lie strictly between 0 and 1.")
  }
  if (out$n_perm < 1) rlang::abort("`n_perm` must be at least 1.")
  if (length(out$phenotype_column) != 1L || !nzchar(out$phenotype_column)) {
    rlang::abort("`phenotype_column` must be one non-empty column name.")
  }
  if (length(out$phenotype_reference) != 1L ||
      !nzchar(out$phenotype_reference)) {
    rlang::abort("`phenotype_reference` must be one non-empty level.")
  }
  if (out$phenotype_alpha <= 0 || out$phenotype_alpha >= 1) {
    rlang::abort("`phenotype_alpha` must lie strictly between 0 and 1.")
  }
  if (out$ambiguity_level <= 0 || out$ambiguity_level >= 1) {
    rlang::abort("`ambiguity_level` must lie strictly between 0 and 1.")
  }
  if (out$n_ambiguity_boot < 0L) {
    rlang::abort("`n_ambiguity_boot` cannot be negative.")
  }
  if (length(out$n_cores) != 1L || is.na(out$n_cores) || out$n_cores < 1L) {
    rlang::abort("`n_cores` must be a positive integer.")
  }
  # A threshold the permutation count cannot reach can never be met, so a run
  # configured that way would report nothing for a reason invisible in its
  # output.
  if (1 / (out$n_perm + 1) > out$alpha) {
    rlang::abort(paste0(
      "With ", out$n_perm, " permutations the smallest attainable p-value is ",
      signif(1 / (out$n_perm + 1), 3), ", which cannot reach alpha = ",
      out$alpha, ". Raise `n_perm` to at least ",
      ceiling(1 / out$alpha) - 1L, "."
    ), class = "chorale_unreachable_alpha")
  }
  structure(out, class = c("chorale_control", "list"))
}

#' @export
print.chorale_control <- function(x, ...) {
  cat("<chorale_control>\n")
  for (nm in names(x)) {
    cat(sprintf("  %-24s %s\n", nm, format(x[[nm]])))
  }
  cat("  smallest attainable p-value:", signif(1 / (x$n_perm + 1), 3), "\n")
  invisible(x)
}

#' Apply named overrides to a control object
#'
#' A caller changing one decision should not have to restate the rest, and a
#' name that is not a decision should be reported rather than ignored.
#'
#' @keywords internal
#' @noRd
chorale_merge_control <- function(control, overrides) {
  if (!inherits(control, "chorale_control")) {
    if (is.list(control)) {
      control <- do.call(chorale_control, control)
    } else {
      rlang::abort("`control` must be a chorale_control object or a list.")
    }
  }
  if (length(overrides) == 0) return(control)
  unknown <- setdiff(names(overrides), names(control))
  if (length(unknown) > 0) {
    rlang::abort(paste0(
      "Unknown setting(s): ", paste(unknown, collapse = ", "),
      ". See chorale_control() for the decisions a run takes."
    ))
  }
  do.call(chorale_control, utils::modifyList(unclass(control), overrides))
}
