#' Emit interpretable outputs from a chorale fit
#'
#' Writes the production outputs specified in `AGENT_PLAN.md` Section 9
#' (`factors.tsv`, `markers.tsv`, `loadings_<modality>.tsv`,
#' `scores_<modality>.tsv`, `associations.tsv`, `concordance.tsv`,
#' `bounds.tsv`, `controls.tsv`, `factors.gmt`, `chorale.h5ad`,
#' `chorale_mae.rds`, `report.html`) into `path`. No factor is reported
#' without a pathway definition and at least two marker features per
#' modality; a factor that cannot be named biologically is reported as
#' unresolved and excluded from interpretation. Depends on `chorale_fit()`,
#' `chorale_bound()` and `chorale_null()` output, none of which is yet
#' implemented (see `analysis/feasibility.md`).
#'
#' @param fit A `chorale_fit` object, as returned by [chorale_fit()].
#' @param bound A `chorale_bound` object, as returned by [chorale_bound()].
#' @param null A `chorale_null` object, as returned by [chorale_null()].
#' @param path Directory to write the output files into.
#'
#' @returns Invisibly, `path`.
#' @export
#' @examplesIf FALSE
#' # Not yet implemented; depends on chorale_fit(), chorale_bound() and
#' # chorale_null().
#' chorale_report(fit, bound, null, path = tempdir())
chorale_report <- function(fit, bound, null, path) {
  lifecycle::signal_stage("experimental", "chorale_report()")
  rlang::abort(
    paste(
      "chorale_report() is not yet implemented; it depends on chorale_fit(),",
      "chorale_bound() and chorale_null()."
    ),
    class = "chorale_not_implemented"
  )
}
