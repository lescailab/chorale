#' Permutation calibration and stability diagnostics
#'
#' Shuffles genotype labels within cohort and stratum, shuffles modality
#' labels, and refits across at least twenty random initialisations, since
#' ICA is non-convex and non-deterministic and a single fit is not a result
#' (`AGENT_PLAN.md` Section 8.5). Depends on a `chorale_fit` object, and so
#' is blocked on `chorale_fit()`, which is not yet implemented (see
#' `analysis/feasibility.md`).
#'
#' @param fit A `chorale_fit` object, as returned by [chorale_fit()].
#' @param n_permutations Integer, number of label permutations to run.
#' @param n_init Integer, number of random ICA initialisations per fit, at
#'   least 20 per `AGENT_PLAN.md` Section 8.5.
#'
#' @returns A `chorale_null` object (not yet defined).
#' @export
#' @examplesIf FALSE
#' # Not yet implemented; depends on chorale_fit().
#' chorale_null(fit, n_permutations = 1000, n_init = 20)
chorale_null <- function(fit, n_permutations = 1000, n_init = 20) {
  lifecycle::signal_stage("experimental", "chorale_null()")
  rlang::abort(
    "chorale_null() is not yet implemented; it depends on chorale_fit().",
    class = "chorale_not_implemented"
  )
}
