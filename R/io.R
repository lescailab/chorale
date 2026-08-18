#' Treat blank and placeholder strings as missing
#'
#' @keywords internal
#' @noRd
chorale_blank_to_na <- function(col_data) {
  placeholders <- c("", "NA", "na", "N/A", "n/a", "NaN", "null", "NULL",
                    "unknown", "Unknown", "not applicable", "not collected",
                    "missing", ".", "-")
  for (nm in setdiff(colnames(col_data), "sample_id")) {
    v <- col_data[[nm]]
    if (is.character(v) || is.factor(v)) {
      ch <- as.character(v)
      ch[trimws(ch) %in% placeholders] <- NA_character_
      col_data[[nm]] <- ch
    }
  }
  col_data
}

#' Required sample metadata columns for a chorale container
#'
#' @keywords internal
#' @noRd
chorale_required_col_data <- function() {
  c(
    "sample_id", "cohort", "modality", "strain", "phenotype",
    "age_months", "sex", "region", "batch"
  )
}

#' Load a single-modality assay into a chorale container
#'
#' Builds a [SummarizedExperiment::SummarizedExperiment] from a
#' feature-by-sample assay matrix and its per-sample metadata, validating
#' that the metadata carries every column a chorale container requires:
#' `sample_id`, `cohort`, `modality`, `strain`, `phenotype`, `age_months`,
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

  # Blank and placeholder strings are missing values, not levels. Reading a
  # design from a delimited file turns an absent value into "", and a covariate
  # carrying "" alongside one real level would otherwise look like a two-level
  # contrast and be analysed as one.
  col_data <- chorale_blank_to_na(col_data)

  rownames(col_data) <- col_data$sample_id
  assay_list <- stats::setNames(list(assay), assay_name)
  SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = col_data
  )
}

#' Read a committed test fixture
#'
#' Reads one of the fixture matrices under `tests/testthat/fixtures/`,
#' returning a feature-by-sample numeric matrix with its identifiers
#' restored. The fixtures are derived from the validated production layers by
#' `data-raw/fixtures.R`, so they carry real identifiers, real missingness and
#' real distributions at a size continuous integration can run.
#'
#' @param layer Character scalar, one of `"RNA"`, `"PROT"` or `"METAB"`.
#' @param path Directory holding the fixtures. Defaults to the installed
#'   fixture directory, and resolves under `testthat::test_path()` when called
#'   from a test.
#'
#' @returns A list with `assay`, a feature-by-sample numeric matrix, and
#'   `col_data`, the matching per-sample design table.
#' @export
#' @examplesIf dir.exists(system.file("fixtures", package = "chorale"))
#' fx <- chorale_fixture("RNA")
#' dim(fx$assay)
chorale_fixture <- function(layer = c("RNA", "PROT", "METAB"), path = NULL) {
  layer <- match.arg(layer)
  if (is.null(path)) {
    path <- system.file("fixtures", package = "chorale")
    if (!nzchar(path)) {
      path <- file.path("tests", "testthat", "fixtures")
    }
  }

  matrix_file <- file.path(path, paste0(layer, "_matrix.parquet"))
  design_file <- file.path(path, paste0(layer, "_design.tsv"))
  if (!file.exists(matrix_file)) {
    rlang::abort(paste0("Fixture not found: ", matrix_file))
  }

  tbl <- nanoparquet::read_parquet(matrix_file)
  feature_id <- tbl[["feature_id"]]
  tbl[["feature_id"]] <- NULL
  assay <- as.matrix(tbl)
  rownames(assay) <- feature_id

  col_data <- utils::read.delim(design_file, stringsAsFactors = FALSE)
  col_data <- col_data[match(colnames(assay), col_data$sample_id), , drop = FALSE]
  rownames(col_data) <- col_data$sample_id

  list(assay = assay, col_data = col_data)
}
