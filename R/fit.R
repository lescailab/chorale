#' Estimate the shared latent state across modalities
#'
#' The estimator specified in `AGENT_PLAN.md` Section 8.3: per-modality
#' non-Gaussian factor recovery by ICA, cross-modality matching of recovered
#' error distributions, and anchor weighting from the design covariates. Its
#' concrete form is gated by the Stage 3 feasibility checks (`chorale_null()`
#' preconditions; `AGENT_PLAN.md` Section 6): whether the linear
#' non-Gaussian identification route applies, or whether the
#' structural-sparsity route of `MATHEMATICAL_FOUNDATION.md` Section 4
#' becomes the primary argument instead. That decision has not yet been
#' made (see `analysis/feasibility.md`), so this function is not yet
#' implemented.
#'
#' @param containers A named list of [SummarizedExperiment::SummarizedExperiment]
#'   objects, one per modality, as returned by [chorale_load()].
#' @param n_factors Integer, or named integer vector per modality, giving the
#'   number of latent factors to recover per modality. Intended to default to
#'   the Stage 3 G3 detectability threshold once that gate has run.
#' @param gene_sets A list of curated gene sets used to constrain factor
#'   loadings towards a pathway definition at estimation (PLIER-style;
#'   `AGENT_PLAN.md` Section 8.3).
#'
#' @returns A `chorale_fit` object (not yet defined).
#' @export
#' @examplesIf FALSE
#' # Not yet implemented; see `analysis/feasibility.md` for the blocking
#' # Stage 3 gate results.
#' chorale_fit(containers, n_factors = 5, gene_sets = list())
chorale_fit <- function(containers, n_factors, gene_sets) {
  lifecycle::signal_stage("experimental", "chorale_fit()")
  rlang::abort(
    paste(
      "chorale_fit() is not yet implemented.",
      "Its design depends on the Stage 3 feasibility gate results",
      "(non-Gaussianity, pairwise difference, detectability; see",
      "AGENT_PLAN.md Section 6 and analysis/feasibility.md), which have",
      "not yet been produced."
    ),
    class = "chorale_not_implemented"
  )
}
