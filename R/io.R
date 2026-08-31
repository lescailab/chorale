#' The vocabulary the design covariates are read in
#'
#' Cohorts assembled by different groups label the same thing differently. One
#' records controls as `control`, another as `ctrl`; one records sex as `male`,
#' another as `M`. Nothing downstream can
#' see that these agree, so a comparison anchored on the design would find no
#' level in common and match on nothing.
#'
#' The registry holds only vocabulary that is general across biology: the
#' phenotype resolves to `case` and `control`, which is what the estimand is
#' about whatever the condition under study; sex resolves to `F` and `M`; a
#' genotype resolves to `carrier` and `non-carrier`, so a design carrying
#' genotype separately from phenotype status can use genotype as secondary
#' adjusted evidence.
#'
#' Names particular to one organism, model or study are deliberately absent. A
#' study whose labels the registry does not cover supplies its own entry to
#' [chorale_load()], which is where knowledge of that study belongs; a label the
#' registry does not recognise is left as it is rather than guessed at.
#'
#' @returns A named list, one entry per covariate, each a named character vector
#'   mapping a lower-cased label to its canonical value.
#' @export
#' @examples
#' names(chorale_label_registry())
#' # A study whose labels the registry does not cover extends it.
#' registry <- chorale_label_registry()
#' registry$phenotype <- c(registry$phenotype, "group_a" = "case")
#' registry$treatment <- c("intervention" = "treated", "baseline" = "untreated")
#' names(registry)
chorale_label_registry <- function() {
  list(
    phenotype = c(
      "case" = "case", "cases" = "case",
      "control" = "control", "controls" = "control", "ctrl" = "control",
      "control group" = "control"
    ),
    sex = c(
      "f" = "F", "female" = "F", "females" = "F", "woman" = "F", "women" = "F",
      "m" = "M", "male" = "M", "males" = "M", "man" = "M", "men" = "M"
    ),
    genotype = c(
      "carrier" = "carrier",
      "non-carrier" = "non-carrier", "noncarrier" = "non-carrier"
    )
  )
}

#' Put design labels in one vocabulary across cohorts
#'
#' Applies [chorale_label_registry()] to a design table. A label the registry
#' does not recognise is left as it is, since inventing a mapping for an
#' unrecognised label would be worse than reporting it.
#'
#' @param col_data A design table.
#' @param labels A registry, as returned by [chorale_label_registry()], or a
#'   list of entries to merge into it.
#'
#' @returns The design table with the covariates the registry covers put in one
#'   vocabulary.
#' @keywords internal
#' @noRd
chorale_canonical_labels <- function(col_data, labels = chorale_label_registry()) {
  for (cv in names(labels)) {
    if (!cv %in% colnames(col_data)) next
    map <- labels[[cv]]
    ch <- as.character(col_data[[cv]])
    hit <- match(tolower(trimws(ch)), names(map))
    ch[!is.na(hit)] <- unname(map[hit[!is.na(hit)]])
    col_data[[cv]] <- ch
  }
  col_data
}

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
  # Only what the estimand cannot do without: an identifier for each sample and
  # the contrast being estimated. Everything else is a covariate that sharpens
  # the comparison where it is present and is absent without consequence where
  # it is not, so requiring it would exclude designs the method can handle.
  c("sample_id", "phenotype")
}

#' Load a single-modality assay into a chorale container
#'
#' Builds a [SummarizedExperiment::SummarizedExperiment] from a
#' feature-by-sample assay matrix and its per-sample metadata. Two columns are
#' required, `sample_id` and `phenotype`, since an estimand defined as a
#' contrast cannot be formed without them. Every other column is a covariate
#' that sharpens the comparison where the modalities share it and is absent
#' without consequence where they do not, so none is required.
#'
#' Labels are put in one vocabulary as the data are read, through
#' [chorale_label_registry()], because cohorts assembled by different groups
#' record the same thing differently and nothing downstream can see that
#' `control` and `ctrl` agree.
#'
#' This is the common container every downstream chorale function expects, one
#' call per modality. [chorale_check_design()] reports what a collection of
#' designs can and cannot support before any of it is fitted.
#'
#' @param assay A feature-by-sample numeric matrix. Column names must match
#'   `col_data$sample_id`.
#' @param col_data A data frame of per-sample metadata, one row per column of
#'   `assay`, carrying at least `sample_id` and `phenotype`.
#' @param labels The vocabulary the design covariates are read in, as returned
#'   by [chorale_label_registry()]. Supply an extended registry where a study
#'   uses labels the default does not cover.
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
chorale_load <- function(assay, col_data, assay_name = "counts",
                         labels = chorale_label_registry()) {
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
  # Cohorts label the same thing differently, so the anchoring covariates are
  # put in one vocabulary before anything tries to compare them.
  col_data <- chorale_canonical_labels(col_data, labels)

  rownames(col_data) <- col_data$sample_id
  assay_list <- stats::setNames(list(assay), assay_name)
  SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = col_data
  )
}

#' Read a committed test fixture
#'
#' Reads one of the synthetic fixture matrices shipped with the package,
#' returning a feature-by-sample numeric matrix with its identifiers
#' restored. The fixtures are generated deterministically by
#' `data-raw/fixtures.R`; they contain generic identifiers and controlled
#' missingness at a size continuous integration can run.
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

  matrix_file <- file.path(path, paste0(layer, "_matrix.tsv"))
  design_file <- file.path(path, paste0(layer, "_design.tsv"))
  if (!file.exists(matrix_file)) {
    rlang::abort(paste0("Fixture not found: ", matrix_file))
  }

  tbl <- utils::read.delim(matrix_file, check.names = FALSE,
                           stringsAsFactors = FALSE)
  feature_id <- tbl[["feature_id"]]
  tbl[["feature_id"]] <- NULL
  assay <- as.matrix(tbl)
  rownames(assay) <- feature_id

  col_data <- utils::read.delim(design_file, stringsAsFactors = FALSE)
  col_data <- col_data[match(colnames(assay), col_data$sample_id), , drop = FALSE]
  rownames(col_data) <- col_data$sample_id

  list(assay = assay, col_data = col_data)
}

#' What a collection of designs can support, before anything is fitted
#'
#' The estimator anchors on the covariates the modalities share, so a design
#' that looks complete on its own can still leave nothing to compare: a
#' covariate present everywhere but recorded under different labels shares no
#' level, and a covariate constant within a cohort carries no contrast. Both
#' produce an empty comparison for reasons that are invisible in any one design.
#'
#' This reports, for every covariate, which modalities carry it, whether it
#' varies in each, and which of its levels are common to all. Reading it before
#' fitting is what turns a silent absence of results into a statement about the
#' data.
#'
#' @param designs A named list of design tables, or of containers from
#'   [chorale_load()].
#' @param labels The vocabulary to read the designs in.
#'
#' @returns A data frame with one row per covariate, carrying the modalities
#'   that hold it, whether it varies in all of them, its shared levels, and
#'   whether it can anchor a comparison. Attribute `usable` lists the covariates
#'   that can.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' chorale_check_design(containers)
chorale_check_design <- function(designs, labels = chorale_label_registry()) {
  if (!is.list(designs) || length(designs) < 1) {
    rlang::abort("`designs` must be a list of design tables or containers.")
  }
  if (is.null(names(designs))) {
    names(designs) <- paste0("modality_", seq_along(designs))
  }
  designs <- lapply(designs, function(d) {
    if (inherits(d, "SummarizedExperiment")) {
      d <- as.data.frame(SummarizedExperiment::colData(d))
    }
    chorale_canonical_labels(chorale_blank_to_na(as.data.frame(d)), labels)
  })

  covariates <- setdiff(unique(unlist(lapply(designs, colnames))),
                        c("sample_id", "modality"))
  rows <- lapply(covariates, function(cv) {
    present <- names(designs)[vapply(designs, function(d) cv %in% colnames(d),
                                     logical(1))]
    varies <- vapply(designs[present], function(d) {
      length(unique(stats::na.omit(d[[cv]]))) >= 2
    }, logical(1))
    # Reported as continuous only where a modality realises more than two
    # values, which is narrower than the rule the fit applies:
    # chorale_covariate_kind() calls any numeric covariate continuous and gives
    # it one standardised slope. The difference is in what is shown, not in what
    # is usable. A numeric covariate taking two values reads to an analyst as a
    # contrast, so its two values are listed rather than the word continuous,
    # and `can_anchor` agrees with the fit either way.
    numeric_all <- all(vapply(designs[present], function(d) {
      is.numeric(d[[cv]]) && length(unique(stats::na.omit(d[[cv]]))) > 2
    }, logical(1)))
    shared <- if (numeric_all) {
      "continuous"
    } else {
      lv <- lapply(designs[present], function(d) {
        unique(as.character(stats::na.omit(d[[cv]])))
      })
      paste(sort(Reduce(intersect, lv)), collapse = ", ")
    }
    n_shared <- if (numeric_all) Inf else {
      length(Reduce(intersect, lapply(designs[present], function(d) {
        unique(as.character(stats::na.omit(d[[cv]])))
      })))
    }
    data.frame(
      covariate = cv,
      in_all_modalities = length(present) == length(designs),
      varies_in_all = length(present) == length(designs) && all(varies),
      shared_levels = shared,
      can_anchor = length(present) == length(designs) && all(varies) &&
        n_shared >= 2,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$can_anchor, out$covariate), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "usable") <- out$covariate[out$can_anchor]
  out
}

#' Whether any individual appears in more than one modality
#'
#' The whole method rests on the modalities being measured on disjoint sets of
#' individuals: no sample is shared, so the comparison of phenotype effects is
#' anchored on design information rather than on matched individuals. Each
#' container is validated on its own by [chorale_load()], which cannot see that
#' assumption because it holds across containers rather than within one.
#'
#' This collects every modality's `sample_id` values and reports each identifier
#' carried by two or more of them. A collision is not necessarily a shared
#' individual: independent studies reuse plain identifiers, so `S01` in one
#' deposit and `S01` in another may be two animals with one name. It is
#' therefore reported and not refused, and the analyst decides which of the two
#' it is.
#'
#' @param containers A named list of [chorale_load()] containers, of design
#'   tables carrying `sample_id`, or of feature-by-sample matrices whose column
#'   names are the identifiers.
#'
#' @returns A data frame with one row per identifier carried by more than one
#'   modality, giving the identifier, the number of modalities carrying it and
#'   their names. Empty where the collection is disjoint. Attribute `unchecked`
#'   names the modalities that carry no identifier, which are reported on by
#'   neither an empty result nor a populated one.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 30, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' chorale_check_disjoint(containers)
chorale_check_disjoint <- function(containers) {
  if (!is.list(containers) || length(containers) < 1) {
    rlang::abort("`containers` must be a list of modalities.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  ids <- lapply(containers, function(x) {
    v <- if (inherits(x, "SummarizedExperiment")) {
      SummarizedExperiment::colData(x)$sample_id
    } else if (is.data.frame(x)) {
      x$sample_id
    } else {
      colnames(as.matrix(x))
    }
    unique(as.character(stats::na.omit(v)))
  })
  missing_ids <- names(ids)[vapply(ids, length, integer(1)) == 0L]
  if (length(missing_ids) == length(ids)) {
    rlang::abort("No modality carries a `sample_id`; disjointness cannot be checked.")
  }
  if (length(missing_ids) > 0) {
    # A modality with no identifier collides with nothing, which reads as a
    # disjoint collection and is instead a collection that was not checked.
    rlang::warn(paste0(
      "No `sample_id` in ", paste(missing_ids, collapse = ", "),
      ". Disjointness is reported for the remaining modalities only, so these ",
      "carry no verdict either way."),
      class = "chorale_unchecked_disjointness")
  }

  empty <- data.frame(sample_id = character(), n_modalities = integer(),
                      modalities = character(), stringsAsFactors = FALSE)
  attr(empty, "unchecked") <- missing_ids
  all_ids <- unlist(ids, use.names = FALSE)
  shared <- unique(all_ids[duplicated(all_ids)])
  if (length(shared) == 0) return(empty)

  shared <- sort(shared)
  carriers <- lapply(shared, function(id) names(ids)[vapply(ids, function(v) {
    id %in% v
  }, logical(1))])
  out <- data.frame(
    sample_id = shared,
    n_modalities = vapply(carriers, length, integer(1)),
    modalities = vapply(carriers, paste, character(1), collapse = ", "),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "unchecked") <- missing_ids
  out
}

#' Warn where a collection is not on disjoint individuals
#'
#' The check informs rather than decides: a collision may be an accidental reuse
#' of an identifier across independent studies, so the run continues and the
#' analyst reads the names.
#'
#' @keywords internal
#' @noRd
chorale_warn_shared_samples <- function(containers) {
  collisions <- chorale_check_disjoint(containers)
  if (nrow(collisions) == 0) return(invisible(collisions))
  shown <- utils::head(collisions, 5L)
  lines <- paste0("  ", shown$sample_id, ": ", shown$modalities)
  if (nrow(collisions) > nrow(shown)) {
    lines <- c(lines, paste0("  and ", nrow(collisions) - nrow(shown),
                             " further identifier(s)"))
  }
  rlang::warn(paste0(
    "The modalities are not on disjoint samples: ", nrow(collisions),
    " identifier(s) appear in more than one modality.\n",
    paste(lines, collapse = "\n"),
    "\nA shared identifier may be one individual measured twice, which the ",
    "method does not assume, or the same name given to two individuals by two ",
    "studies. Neither is decided here."),
    class = "chorale_shared_samples")
  invisible(collisions)
}
