#' Required sample metadata columns for a chorale container
#'
#' @keywords internal
#' @noRd
chorale_required_col_data <- function() {
  c(
    "sample_id", "cohort", "modality", "strain", "genotype",
    "age_months", "sex", "region", "batch"
  )
}

#' Load a single-modality assay into a chorale container
#'
#' Builds a [SummarizedExperiment::SummarizedExperiment] from a
#' feature-by-sample assay matrix and its per-sample metadata, validating
#' that the metadata carries every column a chorale container requires:
#' `sample_id`, `cohort`, `modality`, `strain`, `genotype`, `age_months`,
#' `sex`, `region`, `batch` (`MATHEMATICAL_FOUNDATION.md`, `AGENT_PLAN.md`
#' Section 8.1). This is the common container every downstream chorale
#' function expects, one call per modality.
#'
#' @param assay A feature-by-sample numeric matrix. Column names must match
#'   `col_data$sample_id`.
#' @param col_data A data frame of per-sample metadata, one row per column of
#'   `assay`, carrying at least the required columns listed above.
#' @param assay_name Character scalar, the name to give the assay in the
#'   returned container.
#'
#' @returns A [SummarizedExperiment::SummarizedExperiment] with `assay_name`
#'   as its single assay and `col_data` as `colData`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 30, seed = 1)
#' se <- chorale_load(sim$modalities[[1]], sim$col_data[[1]])
#' se
chorale_load <- function(assay, col_data, assay_name = "counts") {
  if (!is.matrix(assay)) {
    rlang::abort("`assay` must be a matrix, features in rows and samples in columns.")
  }
  if (!is.data.frame(col_data)) {
    rlang::abort("`col_data` must be a data frame.")
  }

  missing_cols <- setdiff(chorale_required_col_data(), names(col_data))
  if (length(missing_cols) > 0) {
    rlang::abort(paste0(
      "`col_data` is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!identical(colnames(assay), col_data$sample_id)) {
    rlang::abort("`colnames(assay)` must match `col_data$sample_id`, in order.")
  }

  rownames(col_data) <- col_data$sample_id
  assay_list <- stats::setNames(list(assay), assay_name)
  SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = col_data
  )
}
