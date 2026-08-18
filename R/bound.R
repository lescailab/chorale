#' Bound the cross-modality coupling
#'
#' Reports Fréchet-style partial-identification bounds on the cross-modality
#' coupling alongside a point estimate, and how far the anchor covariates
#' (strain, phenotype, age, sex) narrow those bounds
#' (`MATHEMATICAL_FOUNDATION.md` Section 2; `AGENT_PLAN.md` Section 8.4).
#' Depends on a `chorale_fit` object, and so is blocked on `chorale_fit()`,
#' which is not yet implemented (see `analysis/feasibility.md`).
#'
#' @param fit A `chorale_fit` object, as returned by [chorale_fit()].
#'
#' @returns A `chorale_bound` object (not yet defined).
#' @export
#' @examplesIf FALSE
#' # Not yet implemented; depends on chorale_fit().
#' chorale_bound(fit)
chorale_bound <- function(fit) {
  lifecycle::signal_stage("experimental", "chorale_bound()")
  rlang::abort(
    "chorale_bound() is not yet implemented; it depends on chorale_fit().",
    class = "chorale_not_implemented"
  )
}
